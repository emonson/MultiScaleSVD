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

class DetailImageFlow(object):

	def __init__(self, data_source, input_link):
		"""Parallel coordinates view constructor needs a valid DataSource plus
		and external annotation link (from the icicle view).
		"""		
		
		self.ds = data_source
		self.input_link = input_link
		
		# Set up callback to listen for changes in IcicleView selections
		self.input_link.AddObserver("AnnotationChangedEvent", self.InputSelectionCallback)
		
		# Set up an annotation link for other views to monitor selected image
		# This link will be created here, but used back and forth between this detail view
		# and the icicle view. It carries the current "scale" selected in both.
		self.output_link = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		self.output_link.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		# Here I'm outputting the "scale" of the selection, so I'm not sure it matters
		# whether output is index or other content type...
		self.output_link.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
		# Set up callback which will work either internally or triggered by change from icicle view.
		self.output_link.AddObserver("AnnotationChangedEvent", self.ScaleSelectionCallback)
				
		# Create a set of empty image stacks for use in empty selections
		# but use the dimensions of a "projected image" so scales match
		example_image = self.ds.GetProjectedImages([0])
		example_image.UpdateInformation()
		(xMin, xMax, yMin, yMax, zMin, zMax) = example_image.GetWholeExtent()
		(xSpacing, ySpacing, zSpacing) = example_image.GetSpacing()
		(x0, y0, z0) = example_image.GetOrigin()
		blankR = xMax - xMin + 1
		blankC = yMax - yMin + 1
		numScales = len(self.ds.ScaleMaxDim)		# Note: accessing member variable directly
		self.blank_image_list = []
		self.blank_image_weights = []
		self.numImagesList = []
		nDim = 3		# 3-dim for now...
		for dd in range(numScales):
			images_linear = N.zeros( blankR*blankC*nDim, dtype='float')	
			intensity = VN.numpy_to_vtk(images_linear, deep=True)
			intensity.SetName('PNGImage')
	
			imageData = vtk.vtkImageData()
			imageData.SetOrigin(x0, y0, z0)
			imageData.SetSpacing(xSpacing, ySpacing, zSpacing)
			imageData.SetExtent(xMin,xMax,yMin,yMax,0,nDim-1)
			imageData.GetPointData().AddArray(intensity)
			imageData.GetPointData().SetActiveScalars('PNGImage')
			
			self.blank_image_list.append(imageData)
			self.blank_image_weights.append(0.1*N.ones(nDim, dtype='float'))
			self.numImagesList.append(nDim)
		
		for dd in range(len(self.blank_image_list)):
			self.blank_image_list[dd].UpdateInformation()
				
		# Create a BrBg7 lookup table
		self.lut = self.ds.GetDivergingLUT('BrBg')
		
 		self.colorList = []
 		self.resliceList = []
 		self.assemblyList = []
 		self.numScales = len(self.ds.ScaleMaxDim)		# Note: accessing member variable

		self.expSpread = 0.5
		
		self.renderer = vtk.vtkRenderer()
		self.cam = self.renderer.GetActiveCamera()

		# renderer.SetBackground(0.15, 0.12, 0.1)
		# cc0,cc1,cc2 = [float(ccVal)/255.0 for ccVal in [40,40,40]]
		cc0,cc1,cc2 = [float(ccVal)/255.0 for ccVal in [60, 60, 60]]
		self.renderer.SetBackground(cc0,cc1,cc2)
				
		self.highlightRect = vtk.vtkOutlineSource()
		self.highlightMapper = vtk.vtkPolyDataMapper()
		self.highlightMapper.SetInputConnection(self.highlightRect.GetOutputPort(0))
		self.highlightActor = vtk.vtkActor()
		self.highlightActor.SetMapper(self.highlightMapper)
		self.highlightActor.GetProperty().SetColor(1,0.8,0.2)
		self.highlightActor.GetProperty().SetLineWidth(3.0)
		self.highlightActor.GetProperty().SetOpacity(0.6)
		# Setting as if nothing picked even though initialized position & orientation to actor0
		self.highlightIndex = -1
		self.highlightActor.SetPickable(False)
		self.highlightActor.SetVisibility(False)
		self.renderer.AddActor(self.highlightActor)

		# Create the slider which will control the image positions
		self.sliderRep = vtk.vtkSliderRepresentation2D()
		self.sliderRep.SetMinimumValue(0)
		self.sliderRep.SetMaximumValue(self.numScales-1)
		self.sliderRep.SetValue(0)
		
		# For remembering the previous slider setting when switching data sets
		self.prevSliderValue = int(self.numScales/2.0)
		
		# And need to keep track of whether to reset view
		self.needToResetCamera = True
		
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
		# Setting self.maxAngle in SetFlowDirection()
		self.FlowDirection = Direction.Vertical
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

		# Done setting up all of the image stuff, so call selection callback
		# to set up view with blank images
		# self.input_link.InvokeEvent('AnnotationChangedEvent')
		self.InputSelectionCallback(self.input_link,None)
		
		self.cam.ParallelProjectionOn()


	# --------------------------------------------------------
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
				delta = -deltaY		# Flowing top to bottom in detail view
			# Moving a little less with mouse for detail view
			new_value = slider_value + float(delta)/5.0
			self.sliderWidget.GetRepresentation().SetValue(new_value)
			self.sliderWidget.InvokeEvent("InteractionEvent")
			self.window.Render()
			
		else:
			self.interactorStyle.OnMouseMove()
			
	# --------------------------------------------------------
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
			self.maxAngle = 80.0

		else:
			self.sliderRep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
			self.sliderRep.GetPoint1Coordinate().SetValue(0.93, 0.6)
			self.sliderRep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
			self.sliderRep.GetPoint2Coordinate().SetValue(0.93, 0.4)
			self.maxAngle = 0.0
		
		if direction != prev_dir:
			if direction == Direction.Horizontal:
				self.cam.Azimuth(-0)
				self.cam.Elevation(10)
			else:
				self.cam.Elevation(-10)
				self.cam.Azimuth(0)
		
		# Update the widget
		self.sliderRep.Modified()
		
	# --------------------------------------------------------
	def InputSelectionCallback(self, caller, event):
		"""This is the callback that tracks changes in the parallel coordinates chart selection
		(pedigree ids) and sets the input images for the image flow view accordingly.
		Note: The image_flow selection should always contain a selection node, but it will
		have no tuples in the selection list if the PC selection has been cleared.
		This new version for variable dimensionality inverts the image stacks so there is
		as list of stacks with one stack per scale."""
		
		for prop in self.assemblyList:
			self.renderer.RemoveViewProp(prop)
		
		annSel = caller.GetCurrentSelection()
		# Note: When image_flow selection is cleared the current selection does NOT contain any tuples
		if annSel.GetNumberOfNodes() > 0:
			idxVtk = annSel.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				
				# New selection list, so get new data
				idxArr = VN.vtk_to_numpy(idxVtk)
				print "Image Flow input to Nup Flow ", idxArr
				
				# Only allowing single index [0]
				# Returning a list of image stacks
				self.imStackList = self.ds.GetDetailImages(idxArr.tolist()[0])
				i_array_name = 'DiffIntensity'
				
				self.imWeightArrList = self.ds.GetDetailWeights(idxArr.tolist()[0])

			else:
								
				self.imStackList = self.blank_image_list
				i_array_name = 'PNGImage'
				
				self.imWeightArrList = self.blank_image_weights

				# Comment this out if we don't want the view resetting on empty selections
				# self.needToResetCamera = True
										
			# Clear out scale ID map (I think this is just in case we need to reverse order...?)
			self.scale_dict = {}
			self.numScales = len(self.imStackList)		# Directly accessing member variable
			for index in range(self.numScales):
				self.scale_dict[index] = index
									
			self.colorList = []
			self.resliceList = []
			self.assemblyList = []
			self.numImagesList = []
			
			# Here is where we have to do things in a new order since there is one
			# stack per scale, so can form assemblies directly instead of adding
			# to each assembly in order as we go through the stacks.
			for nn, imStack in enumerate(self.imStackList):
				imStack.UpdateInformation()
				(xMin, xMax, yMin, yMax, zMin, zMax) = imStack.GetWholeExtent()
				(xSpacing, ySpacing, zSpacing) = imStack.GetSpacing()
				(x0, y0, z0) = imStack.GetOrigin()
				
				center = [x0 + xSpacing * 0.5 * (xMin + xMax),
						  y0 + ySpacing * 0.5 * (yMin + yMax),
						  z0 + zSpacing * 0.5 * (zMin + zMax)]
								
				# Adjust a blue-white-red lookup table
				# NOTE: Only basing on single image stack for now...
				i_range = N.array(imStack.GetPointData().GetArray(i_array_name).GetRange())
				i_ext = abs(i_range.min()) if (abs(i_range.min()) > abs(i_range.max())) else abs(i_range.max())
				if abs(i_range[1]-i_range[0]) < 1e-10: i_ext = 1024
				self.lut.SetRange(-i_ext,i_ext)
				# self.lut.Build()
				
				# Map the image stack through the lookup table
				color = vtk.vtkImageMapToColors()
				color.SetLookupTable(self.lut)
				color.SetInput(imStack)
				color.UpdateWholeExtent()		# since extent of images might be changing
				color.Update()
				self.colorList.append(color)
				
				assembly = vtk.vtkAssembly()
				self.assemblyList.append(assembly)

				numImages = zMax-zMin+1
				self.numImagesList.append(numImages)
				for ii in range(numImages):
						
					zpos = z0 + zSpacing*(zMin+ii)
					axial = vtk.vtkMatrix4x4()
					axial.DeepCopy((1, 0, 0, center[0],
									0, 1, 0, center[1],
									0, 0, 1, zpos,
									0, 0, 0, 1))
									
					# Extract a slice in the desired orientation
					reslice = vtk.vtkImageReslice()
					reslice.SetInputConnection(self.colorList[nn].GetOutputPort())
					reslice.SetOutputDimensionality(2)
					reslice.SetResliceAxes(axial)
					reslice.SetInterpolationModeToNearestNeighbor()
					
					# NOTE: Only basing spacing on first stack
					if ii == 0 and nn == 0:
						reslice.Update()
						tmp = reslice.GetOutput().GetBounds()
						if self.FlowDirection == Direction.Horizontal:
							self.flowSpacing = float(tmp[1]-tmp[0])*1.1
							self.nupSpacing = float(tmp[3]-tmp[2])*1.1
						else:
							self.flowSpacing = float(tmp[3]-tmp[2])*1.1
							self.nupSpacing = float(tmp[1]-tmp[0])*1.1
						sliceSpacing = reslice.GetOutput().GetSpacing()[2]
						
					# Not going to access these later, just keeping around so they
					# won't be garbage collected
					self.resliceList.append(reslice)
					
					# Create a list of actors with each image already assigned
					actor = vtk.vtkImageActor()
					actor.SetInput(self.resliceList[sum(self.numImagesList[:nn]) + ii].GetOutput())
					
					# NOTE: Fragile test for blank image
					if i_array_name == 'PNGImage':
						actor.SetPickable(False)
					else:
						actor.SetPickable(True)
					
					# Set opacity according to relative magnitude of abs(wavelet coeff)
					actor.SetOpacity(self.imWeightArrList[nn][ii])
					
					if self.FlowDirection == Direction.Horizontal:
						actor.SetPosition(0, ii*self.nupSpacing, 0)
					else:
						actor.SetPosition(ii*self.nupSpacing, 0, 0)
					
					self.assemblyList[nn].AddPart(actor)
			
			# Adding assemblies to renderer after all completed
			for ii in range(self.numScales):
				self.renderer.AddActor(self.assemblyList[ii])
				
			# Seems to work best if I set a temporary bounds here and don't do it again
			self.highlightRect.SetBounds(self.assemblyList[0].GetBounds())


			# Get the current scale, if there is one, from the output_link selection list
# 			scaleSel = self.output_link.GetCurrentSelection()
			
			# Note: Empty selection should still contain a node, but no tuples
# 			if (scaleSel.GetNumberOfNodes() > 0) and (scaleSel.GetNode(0).GetSelectionList().GetNumberOfTuples() > 0):
# 				# This should only contain a single value or none
# 				scaleVal = scaleSel.GetNode(0).GetSelectionList().GetValue(0)
# 			else:
# 				scaleVal = -1
# 			
# 			if (scaleVal > 0) and (scaleVal < len(self.assemblyList)):
# 				self.highlightRect.SetBounds(self.assemblyList[scaleVal].GetBounds())
# 				self.highlightActor.SetOrientation(self.assemblyList[scaleVal].GetOrientation())
# 				tmpPos = self.assemblyList[scaleVal].GetPosition()
# 				usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
# 				self.highlightActor.SetPosition(usePos)
# 				self.highlightIndex = scaleVal
# 				self.highlightActor.SetPickable(False)
# 				self.highlightActor.SetVisibility(True)
# 			else:
# 				self.highlightRect.SetBounds(self.assemblyList[0].GetBounds())
# 				self.highlightActor.SetOrientation(self.assemblyList[0].GetOrientation())
# 				tmpPos = self.assemblyList[0].GetPosition()
# 				usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
# 				self.highlightActor.SetPosition(usePos)
# 				self.highlightIndex = -1
# 				self.highlightActor.SetPickable(False)
# 				self.highlightActor.SetVisibility(False)

			# Create the slider which will control the image positions
			self.sliderRep.SetMinimumValue(0)
			self.sliderRep.SetMaximumValue(self.numScales-1)
			# self.sliderRep.SetValue(0)
			# self.sliderRep.SetValue(self.prevSliderValue)
			# In case prev value was greater than max
			# self.prevSliderValue = int(self.sliderRep.GetValue())
					
			# Now assign initial positions to the actors
			self.setImagesPosition(self.prevSliderValue)
			
			if self.needToResetCamera:
				# Camera reset based on height, not width...
				# (Xmin0,Xmax0,Ymin0,Ymax0,Zmin0,Zmax0) = self.assemblyList[0].GetBounds()
				(Xmin0,Xmax0,Ymin0,Ymax0,Zmin0,Zmax0) = self.assemblyList[self.prevSliderValue].GetBounds()
				(Xmin1,Xmax1,Ymin1,Ymax1,Zmin1,Zmax1) = self.assemblyList[-1].GetBounds()
				# eps = (Ymax0-Ymin1)*0.1
				eps = (Ymax0-Ymin0)*2.0
				# Scale for height of whole thing
				# self.renderer.ResetCamera(Xmin0,Xmax0,Ymin1-eps,Ymax0+eps,Zmin0,Zmax0)
				# Scale for width of one element
				# self.renderer.ResetCamera(self.assemblyList[self.prevSliderValue].GetBounds())
				self.renderer.ResetCamera(Xmin0,Xmax0,Ymin0-eps,Ymax0+eps,Zmin0,Zmax0)
				# self.cam.Elevation(10)
				self.renderer.ResetCameraClippingRange()
				# NOTE: Fragile test for blank image
				if i_array_name == 'PNGImage':
					self.needToResetCamera = False
				else:
					self.needToResetCamera = False
				
			# Use output_link callback to update highlight properly
			self.ScaleSelectionCallback()

			if event is not None:
				self.window.Render()

		else:
			# If there is no SelectionNode in PC selection -- shouldn't reach here...
			print "PC no selection node to image flow called"


	def ReloadBasisImages(self):
		"""Reload images when switching from wavelet to scaling function images
		"""
		annSel = self.input_link.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxVtk = annSel.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				
				# New selection list, so get new data
				idxArr = VN.vtk_to_numpy(idxVtk)
				print "Image Flow input to Nup Flow ", idxArr
				
				# Only allowing single index [0]
				# Returning a list of image stacks
				self.imStackList = self.ds.GetDetailImages(idxArr.tolist()[0])
				i_array_name = 'DiffIntensity'
				
				self.imWeightArrList = self.ds.GetDetailWeights(idxArr.tolist()[0])

				# Here is where we have to do things in a new order since there is one
				# stack per scale, so can form assemblies directly instead of adding
				# to each assembly in order as we go through the stacks.
				for nn, imStack in enumerate(self.imStackList):
					# Adjust a blue-white-red lookup table
					# NOTE: Only basing on single image stack for now...
					i_range = N.array(imStack.GetPointData().GetArray(i_array_name).GetRange())
					i_ext = abs(i_range.min()) if (abs(i_range.min()) > abs(i_range.max())) else abs(i_range.max())
					if abs(i_range[1]-i_range[0]) < 1e-10: i_ext = 1024
					self.lut.SetRange(-i_ext,i_ext)
					# self.lut.Build()
					
					self.colorList[nn].SetInput(imStack)
					
					prop3Dcollection = self.assemblyList[nn].GetParts()
					prop3Dcollection.InitTraversal()
					
					tmp_actor = vtk.vtkImageActor.SafeDownCast(prop3Dcollection.GetNextProp3D())

					for ii in range(prop3Dcollection.GetNumberOfItems()):							
						# Set opacity according to relative magnitude of abs(wavelet coeff)
						tmp_actor.SetOpacity(self.imWeightArrList[nn][ii])
						tmp_actor = vtk.vtkImageActor.SafeDownCast(prop3Dcollection.GetNextProp3D())
			
				self.window.Render()

	# --------------------------------------------------------
	def setImagesPosition(self,slider_value):
		xx = N.arange(self.numScales)-slider_value
		yy = 1.0/(1.0+N.exp(-xx/self.expSpread))
		zz = (-self.maxAngle/0.5)*(yy-0.5)
		p_inc = self.flowSpacing*N.cos(2*N.pi*(zz/360.0))/1.0
		pp = N.zeros(xx.shape)
		
		for ii in range(self.numScales):
			
			# Placement seems to be wrt center of image, so position needs to be incremented
			# by half of this one plus half of last one...
			if ii == 0:
				pp[ii] = 0
			else:
				pp[ii] = pp[ii-1] + (p_inc[ii] + p_inc[ii-1])/2.0
			
		# Find the xx=0 point for position offset
		pp0 = N.interp(0.0, xx, pp)
		
		for ii in range(self.numScales):
		
			if self.FlowDirection == Direction.Horizontal:
				# Set the position and rotation
				self.assemblyList[ii].SetOrientation(0, zz[ii], 0)
				self.assemblyList[ii].SetPosition(pp[ii]-pp0, 0, 0)
			else:
				# Set the position and rotation
				self.assemblyList[ii].SetOrientation(zz[ii], 0, 0)
				# NOTE: Negating Y value for Vertical so flows from top to bottom
				self.assemblyList[ii].SetPosition(0, -1.0*(pp[ii]-pp0), 0)
		
		# Also move highlight rectangle if something picked
		if self.highlightIndex >= 0:
			# Make sure not trying to highlight a scale that doesn't exist in this data point
			if self.highlightIndex > (len(self.assemblyList)-1):
				self.manuallySelectScale(len(self.assemblyList)-1)
			else:
				self.highlightActor.SetOrientation(self.assemblyList[self.highlightIndex].GetOrientation())
				tmpPos = self.assemblyList[self.highlightIndex].GetPosition()
				usePos = (tmpPos[0],tmpPos[1],tmpPos[2]+0.01)
				self.highlightActor.SetPosition(usePos)
			
	# --------------------------------------------------------
	def sliderCallback(self,caller,event):
	
		slider_value = caller.GetRepresentation().GetValue()
		if (N.abs(slider_value - N.round(slider_value)) < 0.2):
			slider_value = N.round(slider_value)
			caller.GetRepresentation().SetValue(slider_value)
		self.prevSliderValue = int(slider_value)
		self.setImagesPosition(slider_value)
		
	def endSliderCallback(self,caller,event):
	
		slider_value = caller.GetRepresentation().GetValue()
		slider_value = N.round(slider_value)
		caller.GetRepresentation().SetValue(slider_value)
		self.prevSliderValue = int(slider_value)
		self.setImagesPosition(slider_value)

	# --------------------------------------------------------
	def selectImage(self):
		"""Here just detect whether something was picked and pass the index (or empty list)
		to the output_link, then rely on the output_link callback to update the highlight
		actor and scroll the view to the correct index"""
		
		(x0,y0) = self.interactor.GetLastEventPosition()
		(x,y) = self.interactor.GetEventPosition()
		# DEBUG
		# print "start: ", x0,y0, " end: ", x,y
		picker = self.interactor.GetPicker()
		somethingPicked = picker.PickProp(x,y,self.renderer)
		if somethingPicked:
			pickedProp = picker.GetViewProp()
			index = self.assemblyList.index(pickedProp)
		else:
			index = -1
		
		if index >= 0:
			scale_list = [self.scale_dict[index]]
		else:
			scale_list = []
			
		print "scale picked in detail view: ", scale_list
		
		id_array = N.array(scale_list, dtype='int64')
		id_vtk = VN.numpy_to_vtkIdTypeArray(id_array, deep=True)
		self.output_link.GetCurrentSelection().GetNode(0).SetSelectionList(id_vtk)
		# This event should update selection highlight based on picked scale
		self.output_link.InvokeEvent("AnnotationChangedEvent")
		
	# --------------------------------------------------------
	def manuallySelectScale(self, index):
		"""Here just detect whether something was picked and pass the index (or empty list)
		to the output_link, then rely on the output_link callback to update the highlight
		actor and scroll the view to the correct index"""
		
		if index >= 0:
			scale_list = [self.scale_dict[index]]
		else:
			scale_list = []
					
		id_array = N.array(scale_list, dtype='int64')
		id_vtk = VN.numpy_to_vtkIdTypeArray(id_array, deep=True)
		self.output_link.GetCurrentSelection().GetNode(0).SetSelectionList(id_vtk)
		# This event should update selection highlight based on picked scale
		self.output_link.InvokeEvent("AnnotationChangedEvent")
		
	# --------------------------------------------------------
	def ScaleSelectionCallback(self,caller=None,event=None):
		"""Routine for adding scale value to output annotation link selection list when 
		selection has been changed in image flow. Only allowing single selection.
		Scale output_link will always contain a node, but may contain 0 tuples."""
		
		# Trying only doing highlight & positioning when have real data
		annSel = self.input_link.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		if annSel.GetNumberOfNodes() > 0:
			if annSel.GetNode(0).GetSelectionList().GetNumberOfTuples() == 0:
				return

		# Get current scale
		# The content type is Index for now...
		scaleSel = self.output_link.GetCurrentSelection()
		
		# Note: Empty selection should still contain a node, but no tuples
		if (scaleSel.GetNumberOfNodes() > 0) and (scaleSel.GetNode(0).GetSelectionList().GetNumberOfTuples() > 0):
			# This should only contain a single value or none
			scaleVal = scaleSel.GetNode(0).GetSelectionList().GetValue(0)
		else:
			scaleVal = -1
		
		print "Updating detail view with selected scale: ", scaleVal
		
		if (scaleVal >= 0) and (scaleVal < len(self.assemblyList)):
			self.highlightIndex = scaleVal
			self.highlightActor.SetVisibility(True)
			
			# Need to recenter Y-axis bounds for vertical flow
			(bx0,bx1,by0,by1,bz0,bz1) = self.assemblyList[self.highlightIndex].GetBounds()
			if self.FlowDirection == Direction.Horizontal:
				bxavg = (bx0+bx1)/2.0
				self.highlightRect.SetBounds((bx0-bxavg,bx1-bxavg,by0,by1,bz0,bz1))
			else:
				byavg = (by0+by1)/2.0
				self.highlightRect.SetBounds((bx0,bx1,by0-byavg,by1-byavg,bz0,bz1))
			
			slider_value = self.sliderWidget.GetRepresentation().GetValue()
			if self.highlightIndex != slider_value:
				# Animate 10 steps to new position if non-current image selected
				for vv in N.linspace(slider_value, self.highlightIndex, 10):
					time.sleep(0.02)
					self.sliderWidget.GetRepresentation().SetValue(vv)
					self.sliderWidget.InvokeEvent("InteractionEvent")
					self.window.Render()
		else:
			self.highlightIndex = -1
			self.highlightActor.SetVisibility(False)
			
			self.sliderRep.SetValue(self.prevSliderValue)
	
		# In case prev value was greater than max
		self.prevSliderValue = int(self.sliderRep.GetValue())
		
		# Don't want to call Render if this is the first setup call (internal) or will
		# have problems with extra render window popping up...
		if caller is not None:
			self.window.Render()
				
	# --------------------------------------------------------
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
		

	# --------------------------------------------------------
	def GetOutputAnnotationLink(self):
		return self.output_link
	
	def GetRenderWindow(self):
		return self.window
		
	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.view.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)


# --------------------------------------------------------
if __name__ == "__main__":

	# from tkFileDialog import askopenfilename
	from data_source import DataSource
	
	# data_file = askopenfilename()
	data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100825.mat'

	# DataSource loads .mat file and can generate data from it for other views
	ds = DataSource(data_file)
		
	# Set up an annotation link as if selections were coming from another class
	dummy_link = vtk.vtkAnnotationLink()
	dummy_link.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
	dummy_link.GetCurrentSelection().GetNode(0).SetContentType(2)   # 2 = PedigreeIds, 4 = Indices
	
	if_class = DetailImageFlow(ds, dummy_link)
# 	if_class.GetRenderWindow().SetSize(600,300)
#  	if_class.SetFlowDirection(Direction.Horizontal)
	if_class.GetRenderWindow().SetSize(600,300)
	if_class.SetFlowDirection(Direction.Vertical)
		
	# Fill selection link with dummy IDs
	id_array = N.array([1422],dtype='int64')
	id_list = VN.numpy_to_vtkIdTypeArray(id_array)
	dummy_link.GetCurrentSelection().GetNode(0).SetSelectionList(id_list)
	dummy_link.InvokeEvent("AnnotationChangedEvent")
	
	# Only need to Start() interactor for one view
	if_class.GetRenderWindow().GetInteractor().Start()
	
	
	
	
	
