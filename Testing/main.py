from icicle_view_simple_link import IcicleView
from no_view_simple_link import NoView

import vtk
import vtk.util.numpy_support as VN


iclass = IcicleView('treetest.xml')
nclass = NoView()

iview = iclass.GetView()

ilink = vtk.vtkAnnotationLink()
ilink.GetCurrentSelection().GetNode(0).SetFieldType(3)		# Vertex
ilink.GetCurrentSelection().GetNode(0).SetContentType(2)	# Pedigree Ids

iclass.SetAnnotationLink(ilink)
nclass.SetAnnotationLink(ilink)

def IcicleSelectionCallback(caller, event):
	# Handle updating of RenderWindow since it's not a "View"
	#	and so not covered by vtkViewUpdater
	# ren.ResetCamera()
	annSel = caller.GetCurrentSelection()
	if annSel.GetNumberOfNodes() > 0:
		idxArr = annSel.GetNode(0).GetSelectionList()
		print "main ", VN.vtk_to_numpy(idxArr)

# Set up callback to update 3d render window when selections are changed in 
#	parallel coordinates view
ilink.AddObserver("AnnotationChangedEvent", IcicleSelectionCallback)

iview.GetInteractor().Start()