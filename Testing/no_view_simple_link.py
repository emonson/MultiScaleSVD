import vtk
import vtk.util.numpy_support as VN

class NoView(object):

	def __init__(self):
		self.link = vtk.vtkAnnotationLink()
		self.link.GetCurrentSelection().GetNode(0).SetFieldType(3)		# Vertex
		self.link.GetCurrentSelection().GetNode(0).SetContentType(2)		# Pedigree Ids
		
		# Set up callback to update 3d render window when selections are changed in 
		#	parallel coordinates view
		self.link.AddObserver("AnnotationChangedEvent", self.SelectionCallback)
		

	def SelectionCallback(self, caller, event):
		# Handle updating of RenderWindow since it's not a "View"
		#	and so not covered by vtkViewUpdater
		# ren.ResetCamera()
		annSel = caller.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			print "no ", VN.vtk_to_numpy(idxArr)
		
	def GetAnnotationLink(self):
		return self.link
	
	def SetAnnotationLink(self, externalLink):
		self.link = externalLink
		self.link.AddObserver("AnnotationChangedEvent", self.SelectionCallback)
