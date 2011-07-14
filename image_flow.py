import vtk
import vtk.util.numpy_support as VN
import numpy as N
import math
import time
import os, sys
# sys.path.append("/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
# import libvtkvtgChartsPython as vtgCh
import vtkvtg as vtgCh
from constants import Direction

class ImageFlow(object):

	def __init__(self, data_source, input_link):
		"""Parallel coordinates view constructor needs a valid DataSource plus
		and external annotation link (from the icicle view).
		"""		
		
		self.ds = data_source
		self.input_link = input_link
		
		# Set up callback to listen for changes in IcicleView selections
		self.input_link.AddObserver("AnnotationChangedEvent", self.InputSelectionCallback)
		
		# Set up an annotation link for other views to monitor selected image
		self.output_link = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		self.output_link.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		# The chart seems to force INDEX selection, so I'm using vtkConvertSelection below to get
		# out PedigreeIds...
		self.output_link.GetCurrentSelection().GetNode(0).SetContentType(2)   # 2 = PedigreeIds, 4 = Indices
		# Going to manually create output_link selection list, so not setting up callback for it...
		
		if os.path.isfile(os.path.abspath('BlankImage.png')):
			blank_file = 'BlankImage.png'
		else:
			# For use in app bundles
			blank_file = os.path.join(sys.path[0],'BlankImage.png')
			
		self.blankImageReader = vtk.vtkPNGReader()
		self.blankImageReader.SetFileName(blank_file)
		self.blankImageReader.Update()
		
		tmp = self.blankImageReader.GetOutput().GetBounds()
		self.flowSpacing = float(tmp[1]-tmp[0])*1.1

		# Create a greyscale lookup table
		self.lut = vtk.vtkLookupTable()
		self.lut.SetRange(self.blankImageReader.GetOutput().GetPointData().GetArray('PNGImage').GetRange()) # image intensity range
		self.lut.SetValueRange(0.0, 1.0) # from black to white
		self.lut.SetSaturationRange(0.0, 0.0) # no color saturation
		self.lut.SetRampToLinear()
		self.lut.Build()
		
		# Map the image through the lookup table
		self.color = vtk.vtkImageMapToColors()
		self.color.SetLookupTable(self.lut)
		self.color.SetInput(self.blankImageReader.GetOutput())
		
		self.resliceList = []
		self.actorList = []
		self.numImages = 1

		# Map between indices of images and their Pedigree ID
		self.pedigree_id_dict = {}
		
		blankImageActor = vtk.vtkImageActor()
		blankImageActor.SetInput(self.color.GetOutput())
		blankImageActor.SetPickable(False)
		blankImageActor.SetOpacity(0.1)
		
		self.actorList.append(blankImageActor)

		self.expSpread = 0.5
		self.maxAngle = 80.0
		
		self.renderer = vtk.vtkRenderer()
		self.cam = self.renderer.GetActiveCamera()

		# renderer.SetBackground(0.15, 0.12, 0.1)
		# cc0,cc1,cc2 = [float(ccVal)/255.0 for ccVal in [68,57,53]]
		# self.renderer.SetBackground(cc0,cc1,cc2)
		# cc0,cc1,cc2 = [float(ccVal)/255.0 for ccVal in [60,60,60]]
		cc0,cc1,cc2 = [float(ccVal)/255.0 for ccVal in [160, 160, 160]]
		self.renderer.SetBackground(cc0,cc1,cc2)

		self.renderer.AddActor(self.actorList[0])
				
		self.highlightRect = vtk.vtkOutlineSource()
		self.highlightRect.SetBounds(self.actorList[0].GetBounds())
		self.highlightMapper = vtk.vtkPolyDataMapper()
		self.highlightMapper.SetInputConnection(self.highlightRect.GetOutputPort(0))
		self.highlightActor = vtk.vtkActor()
		self.highlightActor.SetMapper(self.highlightMapper)
		self.highlightActor.GetProperty().SetColor(0,0.5,1.0)
		self.highlightActor.GetProperty().SetLineWidth(6.0)
		self.highlightActor.GetProperty().SetOpacity(0.5)
		self.highlightActor.SetOrientation(self.actorList[0].GetOrientation())
		tmpPos = self.actorList[0].GetPosition()
		usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
		self.highlightActor.SetPosition(usePos)
		# Setting as if nothing picked even though initialized position & orientation to actor0
		self.highlightIndex = -1
		self.highlightActor.SetPickable(False)
		self.highlightActor.SetVisibility(False)
		self.renderer.AddActor(self.highlightActor)

		# Create the slider which will control the image positions
		self.sliderRep = vtk.vtkSliderRepresentation2D()
		self.sliderRep.SetMinimumValue(0)
		self.sliderRep.SetMaximumValue(self.numImages-1)
		self.sliderRep.SetValue(0)
		
		self.window = vtk.vtkRenderWindow()
		self.window.SetSize(600,300)

		self.window.AddRenderer(self.renderer)
				
		# Set up the interaction
		self.interactorStyle = vtk.vtkInteractorStyleImage()
		self.interactor = vtk.vtkRenderWindowInteractor()
		self.interactor.SetInteractorStyle(self.interactorStyle)
		self.window.SetInteractor(self.interactor)
				
		self.sliderWidget = vtk.vtkSliderWidget()
		self.sliderWidget.SetInteractor(self.interactor)
		self.sliderWidget.SetRepresentation(self.sliderRep)
		self.sliderWidget.SetAnimationModeToAnimate()
		self.sliderWidget.EnabledOn()
		
		self.sliderWidget.AddObserver("InteractionEvent", self.sliderCallback)
		self.sliderWidget.AddObserver("EndInteractionEvent", self.endSliderCallback)
		
		# Default flow direction Horizontal
		self.FlowDirection = Direction.Horizontal
		self.SetFlowDirection(self.FlowDirection)
				
		# Set up callback to toggle between inspect modes (manip axes & select data)
		# self.interactorStyle.AddObserver("SelectionChangedEvent", self.selectImage)

		# Create callbacks for mouse events
		self.mouseActions = {}
		self.mouseActions["LeftButtonDown"] = 0
		self.mouseActions["Picking"] = 0
		
		self.interactorStyle.AddObserver("MouseMoveEvent", self.MouseMoveCallback)
		self.interactorStyle.AddObserver("LeftButtonPressEvent", self.LeftButtonPressCallback)
		self.interactorStyle.AddObserver("LeftButtonReleaseEvent", self.LeftButtonReleaseCallback)

		self.cam.ParallelProjectionOn()
		self.renderer.ResetCamera(self.actorList[0].GetBounds())
		self.cam.Elevation(10)
		self.renderer.ResetCameraClippingRange()
		# self.window.Render()


	def LeftButtonPressCallback(self, obj, event):
		self.mouseActions["LeftButtonDown"] = 1
		self.mouseActions["Picking"] = 1

	def LeftButtonReleaseCallback(self, obj, event):
		self.mouseActions["LeftButtonDown"] = 0
		if self.mouseActions["Picking"] == 1:
			self.selectImage()
		self.mouseActions["Picking"] = 0
		self.sliderWidget.InvokeEvent("EndInteractionEvent")
		self.window.Render()
	
	def MouseMoveCallback(self, obj, event):
		(lastX, lastY) = self.interactor.GetLastEventPosition()
		(mouseX, mouseY) = self.interactor.GetEventPosition()
		if self.mouseActions["LeftButtonDown"] == 1:
			deltaX = mouseX - lastX
			deltaY = mouseY - lastY
			# If mouse movement greater than tolerance, then not picking
			if (abs(deltaX) + abs(deltaY)) > 0:
				self.mouseActions["Picking"] = 0
			slider_value = self.sliderWidget.GetRepresentation().GetValue()
			if self.FlowDirection == Direction.Horizontal:
				delta = deltaX
			else:
				delta = deltaY
			new_value = slider_value + float(delta)/7.0
			self.sliderWidget.GetRepresentation().SetValue(new_value)
			self.sliderWidget.InvokeEvent("InteractionEvent")
			self.window.Render()
			
		else:
			self.interactorStyle.OnMouseMove()
			
	def SetFlowDirection(self, dir):
	
		direction = int(dir)
		
		if (dir < 0) or (dir >= Direction.Max):
			print "not an allowed direction"
			return
			
		prev_dir = self.FlowDirection
		
		# Default flow direction horizontal
		self.FlowDirection = direction
		
		if self.FlowDirection == Direction.Horizontal:
			# Here we use normalized display coordinates (0,1) so that the
			#   slider will stay in the same proportionate location if the window
			#   is resized.
			self.sliderRep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
			self.sliderRep.GetPoint1Coordinate().SetValue(0.1 ,0.07)
			self.sliderRep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
			self.sliderRep.GetPoint2Coordinate().SetValue(0.9, 0.07)
		else:
			self.sliderRep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
			self.sliderRep.GetPoint1Coordinate().SetValue(0.93, 0.3)
			self.sliderRep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
			self.sliderRep.GetPoint2Coordinate().SetValue(0.93, 0.7)
		
		if direction != prev_dir:
			if direction == Direction.Horizontal:
				self.cam.Azimuth(-5)
				self.cam.Elevation(10)
			else:
				self.cam.Elevation(-10)
				self.cam.Azimuth(5)
		
		# Update the widget
		self.sliderRep.Modified()
		# self.window.Render()
		
	def InputSelectionCallback(self, caller, event):
		"""This is the callback that tracks changes in the parallel coordinates chart selection
		(pedigree ids) and sets the input images for the image flow view accordingly.
		Note: The PC selection should always contain a selection node, but it will
		have no tuples in the selection list if the PC selection has been cleared."""
		
		for prop in self.actorList:
			self.renderer.RemoveViewProp(prop)
		
		annSel = caller.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		if annSel.GetNumberOfNodes() > 0:
			idxVtk = annSel.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				
				# New selection list, so get new data
				idxArr = VN.vtk_to_numpy(idxVtk)
				print "PC input to image flow ", idxArr
				
				print "Getting image stack"
				self.imStack = self.ds.GetProjectedImages(idxArr.tolist())
				i_array_name = 'Intensity'
				
				print "Got image stack"
				
				# Need to manually generate a map of index to Pedigree ID since
				# vtkImageData returned doesn't contain Pedigree ID tag for each slice
				self.pedigree_id_dict = {}
				for index, ped_id in enumerate(idxArr.tolist()):
					self.pedigree_id_dict[index] = ped_id
				
			else:
				
				print "PC blank into Image Flow!"
				
				self.imStack = self.blankImageReader.GetOutputDataObject(0)
				i_array_name = 'PNGImage'
				
				# Clear out pedigree ID map
				self.pedigree_id_dict = {}
			
			# Calculate the center of the volume
			self.imStack.UpdateInformation()
			(xMin, xMax, yMin, yMax, zMin, zMax) = self.imStack.GetWholeExtent()
			(xSpacing, ySpacing, zSpacing) = self.imStack.GetSpacing()
			(x0, y0, z0) = self.imStack.GetOrigin()
			
			center = [x0 + xSpacing * 0.5 * (xMin + xMax),
					  y0 + ySpacing * 0.5 * (yMin + yMax),
					  z0 + zSpacing * 0.5 * (zMin + zMax)]
			
			
			# Adjust a greyscale lookup table
			self.lut.SetRange(self.imStack.GetPointData().GetArray(i_array_name).GetRange()) # image intensity range
			self.lut.Build()
			
			# Map the image through the lookup table
			self.color.SetInput(self.imStack)
			self.color.UpdateWholeExtent()		# since extent of images might be changing
			self.color.Update()
			
			self.resliceList = []
			self.actorList = []
			self.numImages = zMax-zMin+1
			
			for ii in range(self.numImages):
					
				zpos = z0 + zSpacing*(zMin+ii)
				axial = vtk.vtkMatrix4x4()
				axial.DeepCopy((1, 0, 0, center[0],
								0, 1, 0, center[1],
								0, 0, 1, zpos,
								0, 0, 0, 1))
								
				# Extract a slice in the desired orientation
				reslice = vtk.vtkImageReslice()
				reslice.SetInputConnection(self.color.GetOutputPort())
				reslice.SetOutputDimensionality(2)
				reslice.SetResliceAxes(axial)
				reslice.SetInterpolationModeToNearestNeighbor()
				
				if ii == 0:
					reslice.Update()
					tmp = reslice.GetOutput().GetBounds()
					if self.FlowDirection == Direction.Horizontal:
						self.flowSpacing = float(tmp[1]-tmp[0])*1.1
					else:
						self.flowSpacing = float(tmp[3]-tmp[2])*1.1
					sliceSpacing = reslice.GetOutput().GetSpacing()[2]
				
				self.resliceList.append(reslice)
				
				# Create a list of actors with each image already assigned
				actor = vtk.vtkImageActor()
				actor.SetInput(self.resliceList[ii].GetOutput())
				# NOTE: Fragile test for blank image
				if i_array_name == 'PNGImage':
					actor.SetPickable(False)
					actor.SetOpacity(0.1)
				else:
					actor.SetPickable(True)
					actor.SetOpacity(1.0)
				self.actorList.append(actor)
				
				self.renderer.AddActor(self.actorList[ii])
				
			self.highlightRect.SetBounds(self.actorList[0].GetBounds())
			self.highlightActor.SetOrientation(self.actorList[0].GetOrientation())
			tmpPos = self.actorList[0].GetPosition()
			usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
			self.highlightActor.SetPosition(usePos)
			# Setting as if nothing picked even though initialized position & orientation to actor0
			self.highlightIndex = -1
			self.highlightActor.SetPickable(False)
			self.highlightActor.SetVisibility(False)
			# self.renderer.AddActor(self.highlightActor)		
					
			# Now assign initial positions to the actors
			self.setImagesPosition(0.0)
			
			# Create the slider which will control the image positions
			self.sliderRep.SetMinimumValue(0)
			self.sliderRep.SetMaximumValue(self.numImages-1)
			self.sliderRep.SetValue(0)				
					
			self.renderer.ResetCamera(self.actorList[0].GetBounds())
			# self.cam.Elevation(10)
			self.renderer.ResetCameraClippingRange()
			self.window.Render()

			# Selection gets changed on image reload, so update output annotation link
			self.ImageFlowSelectionChanged()

		else:
			# If there is no SelectionNode in PC selection -- shouldn't reach here...
			print "PC no selection node to image flow called"


	def PrintCameraPosition(self,obj,event):
		print self.renderer.GetActiveCamera()
				
	def setImagesPosition(self,slider_value):
		xx = N.arange(self.numImages)-slider_value
		yy = 1.0/(1.0+N.exp(-xx/self.expSpread))
		zz = (-self.maxAngle/0.5)*(yy-0.5)
		p_inc = self.flowSpacing*N.cos(2*N.pi*(zz/360.0))/1.0
		pp = N.zeros(xx.shape)
		
		for ii in range(self.numImages):
			
			# Placement seems to be wrt center of image, so position needs to be incremented
			# by half of this one plus half of last one...
			if ii == 0:
				pp[ii] = 0
			else:
				pp[ii] = pp[ii-1] + (p_inc[ii] + p_inc[ii-1])/2.0
			
		# Find the xx=0 point for position offset
		pp0 = N.interp(0.0, xx, pp)
		
		for ii in range(self.numImages):
		
			if self.FlowDirection == Direction.Horizontal:
				# Set the position and rotation
				self.actorList[ii].SetOrientation(0, zz[ii], 0)
				self.actorList[ii].SetPosition(pp[ii]-pp0, 0, 0)
			else:
				# Set the position and rotation
				self.actorList[ii].SetOrientation(zz[ii], 0, 0)
				self.actorList[ii].SetPosition(0, pp[ii]-pp0, 0)
		
		# Also move highlight rectangle if something picked
		if self.highlightIndex >= 0:
			self.highlightActor.SetOrientation(self.actorList[self.highlightIndex].GetOrientation())
			tmpPos = self.actorList[self.highlightIndex].GetPosition()
			usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
			self.highlightActor.SetPosition(usePos)
			
	def sliderCallback(self,caller,event):
	
		slider_value = caller.GetRepresentation().GetValue()
		if (N.abs(slider_value - N.round(slider_value)) < 0.2):
			slider_value = N.round(slider_value)
			caller.GetRepresentation().SetValue(slider_value)
		self.setImagesPosition(slider_value)
		
	def endSliderCallback(self,caller,event):
	
		slider_value = caller.GetRepresentation().GetValue()
		slider_value = N.round(slider_value)
		caller.GetRepresentation().SetValue(slider_value)
		self.setImagesPosition(slider_value)

	def selectImage(self):
		(x0,y0) = self.interactor.GetLastEventPosition()
		(x,y) = self.interactor.GetEventPosition()
		picker = self.interactor.GetPicker()
		somethingPicked = picker.PickProp(x,y,self.renderer)
		if somethingPicked:
			pickedProp = picker.GetViewProp()
			self.highlightIndex = self.actorList.index(pickedProp)
			self.highlightActor.SetOrientation(pickedProp.GetOrientation())
			tmpPos = pickedProp.GetPosition()
			usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
			self.highlightActor.SetPosition(usePos)
			self.highlightActor.SetVisibility(True)
			
			slider_value = self.sliderWidget.GetRepresentation().GetValue()
			if self.highlightIndex != slider_value:
				# Animate 10 steps to new position if non-current image selected
				for vv in N.linspace(slider_value, self.highlightIndex, 10):
					time.sleep(0.02)
					self.sliderWidget.GetRepresentation().SetValue(vv)
					self.sliderWidget.InvokeEvent("InteractionEvent")
					self.window.Render()
			self.window.Render()
		else:
			self.highlightActor.SetVisibility(False)
			self.highlightIndex = -1
			self.window.Render()
		
		self.ImageFlowSelectionChanged()
		
	def ImageFlowSelectionChanged(self):
		"""Routine for adding pedigree ID to output annotation link selection list when 
		selection has been changed in image flow. Only allowing single selection for now."""
		if self.highlightIndex >= 0:
			ped_id_list = [self.pedigree_id_dict[self.highlightIndex]]
		else:
			ped_id_list = []
			
		id_array = N.array(ped_id_list, dtype='int64')
		id_vtk = VN.numpy_to_vtkIdTypeArray(id_array, deep=True)
		self.output_link.GetCurrentSelection().GetNode(0).SetSelectionList(id_vtk)
		self.output_link.InvokeEvent("AnnotationChangedEvent")
		
		
	def SetInteractorStyle(self, style):
		self.interactor = self.window.GetInteractor()
		self.interactorStyle = style
		self.window.GetInteractor().SetInteractorStyle(self.interactorStyle)
		self.interactorStyle.AddObserver("MouseMoveEvent", self.MouseMoveCallback)
		self.interactorStyle.AddObserver("LeftButtonPressEvent", self.LeftButtonPressCallback)
		self.interactorStyle.AddObserver("LeftButtonReleaseEvent", self.LeftButtonReleaseCallback)
		
		self.sliderWidget.SetInteractor(self.interactor)
		self.sliderWidget.SetRepresentation(self.sliderRep)
		self.sliderWidget.SetAnimationModeToAnimate()
		self.sliderWidget.EnabledOn()
		
		self.sliderWidget.AddObserver("InteractionEvent", self.sliderCallback)
		self.sliderWidget.AddObserver("EndInteractionEvent", self.endSliderCallback)
		

		
	def GetOutputAnnotationLink(self):
		return self.output_link
	
	def GetRenderWindow(self):
		return self.window
		
	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.view.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)

if __name__ == "__main__":

	# from tkFileDialog import askopenfilename
	from data_source import DataSource
	
	# data_file = askopenfilename()
	data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist1_5c_20100324.mat'

	# DataSource loads .mat file and can generate data from it for other views
	ds = DataSource(data_file)
		
	# Set up an annotation link as if selections were coming from another class
	dummy_link = vtk.vtkAnnotationLink()
	dummy_link.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
	dummy_link.GetCurrentSelection().GetNode(0).SetContentType(2)   # 2 = PedigreeIds, 4 = Indices
	
	if_class = ImageFlow(ds, dummy_link)
	if_class.GetRenderWindow().SetSize(600,300)
 	if_class.SetFlowDirection(Direction.Horizontal)
#	if_class.GetRenderWindow().SetSize(300,600)
#	if_class.SetFlowDirection(Direction.Vertical)
		
	# Fill selection link with dummy IDs
	id_array = N.array(range(50),dtype='int64')
	id_list = VN.numpy_to_vtkIdTypeArray(id_array)
	dummy_link.GetCurrentSelection().GetNode(0).SetSelectionList(id_list)
	dummy_link.InvokeEvent("AnnotationChangedEvent")
	
	# Only need to Start() interactor for one view
	if_class.GetRenderWindow().GetInteractor().Start()
	
	
	
	
	
