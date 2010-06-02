import vtk
import vtk.util.numpy_support as VN

class IcicleView(object):

	def __init__(self, input_file):
		reader1 = vtk.vtkXMLTreeReader()
		reader1.SetFileName("treetest.xml")
		
		self.view = vtk.vtkIcicleView()
		self.view.SetRepresentationFromInput(reader1.GetOutput())
		self.view.SetAreaSizeArrayName("size")
		self.view.SetAreaColorArrayName("level")
		self.view.SetAreaLabelArrayName("name")
		self.view.SetAreaLabelVisibility(True)
		self.view.SetAreaHoverArrayName("name")
		self.view.SetShrinkPercentage(0.05)
		self.view.SetLayerThickness(3.0)
		self.view.UseGradientColoringOff()
		self.view.Update()
		
		rep = self.view.GetRepresentation(0)
		
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		# enum		SelectionContent { 
		#	SELECTIONS, GLOBALIDS, PEDIGREEIDS, VALUES, 
		#	INDICES, FRUSTUM, LOCATIONS, THRESHOLDS, 
		#	BLOCKS 
		# }
		# enum		SelectionField { 
		#	CELL, POINT, FIELD, VERTEX, 
		#	EDGE, ROW 
		# }
		self.link = rep.GetAnnotationLink()
		self.link.GetCurrentSelection().GetNode(0).SetFieldType(3)		# Vertex
		self.link.GetCurrentSelection().GetNode(0).SetContentType(2)		# Pedigree Ids
		
		# Set up callback to update 3d render window when selections are changed in 
		#	parallel coordinates view
		self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)
		
		self.view.ResetCamera()
		self.view.Render()
		
		# self.view.GetInteractor().Start()

	def IcicleSelectionCallback(self, caller, event):
		# Handle updating of RenderWindow since it's not a "View"
		#	and so not covered by vtkViewUpdater
		# ren.ResetCamera()
		annSel = caller.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			print "ice ", VN.vtk_to_numpy(idxArr)
		
	def GetAnnotationLink(self):
		return self.link
	
	def GetView(self):
		return self.view
		
	def SetAnnotationLink(self, externalLink):
		self.link = externalLink
		self.view.GetRepresentation(0).SetAnnotationLink(self.link)
		self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)
