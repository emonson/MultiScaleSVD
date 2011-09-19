import vtk
import vtk.util.numpy_support as VN
import numpy as N
import vtkvtg

class SelectionInspector(object):

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
		
		# self.renderer = vtk.vtkRenderer()

		# Setting as if nothing picked even though initialized position & orientation to actor0
		self.highlightIndex = -1
		
		# self.window = vtk.vtkRenderWindow()
		# self.window.SetSize(600,300)
		# self.window.AddRenderer(self.renderer)
				
		# Set up the interaction
		# self.interactorStyle = vtk.vtkInteractorStyleImage()
		# self.interactor = vtk.vtkRenderWindowInteractor()
		# self.interactor.SetInteractorStyle(self.interactorStyle)
		# self.window.SetInteractor(self.interactor)
				



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
				
				# Preprocess input selection
				
			else:
				
				# No input, blank out
				pass
			
			# Do what we need to do with the input selection




	def selectionMade(self):

		# self.highlightIndex = self.actorList.index(pickedProp)
		# self.ImageFlowSelectionChanged()
		
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
		
		
	def GetOutputAnnotationLink(self):
		return self.output_link
	
	def GetRenderWindow(self):
		return self.window
		
	def GetInteractor(self):
		return self.interactor
		
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
	
	if_class = SelectionInspector(ds, dummy_link)
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
	
	
	
	
	
