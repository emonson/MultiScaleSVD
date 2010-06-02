from data_source import DataSource
from icicle_noview_textured import IcicleNoView
from pcoords_chart import PCoordsChart
from image_flow import ImageFlow
from detail_image_flow import DetailImageFlow
from constants import Direction

from tkFileDialog import askopenfilename
import vtk
import vtk.util.numpy_support as VN
import os

print os.getcwd()

data_file = askopenfilename()
# data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100521.mat'

# DataSource loads .mat file and can generate data from it for other views
ds = DataSource(data_file)

# All view classes have access to an instance of that data source for internal queries
# Note that the only view which will pull and display data right away is the icicle view
#  the other views need to be able to initialize without any data and only pull and show
#  upon the first AnnotationChanged event...
ice_class = IcicleNoView(ds)
ice_class.GetRenderWindow().SetPosition(50,500)
ice_class.GetRenderWindow().SetSize(630,470)
ice_al_out = ice_class.GetOutputAnnotationLink()

# Note: With the way I've implemented the output annotation link in PCoords chart, 
#   it will always have a selection node, but the selection list may have no tuples if
#	it's an empty selection (and event gets fired on every empty selection

pc_class = PCoordsChart(ds)
pc_class.SetInputAnnotationLink(ice_al_out)
pc_class.GetView().GetRenderWindow().SetPosition(50,170)
pc_class.GetView().GetRenderWindow().SetSize(630,300)
pc_al_out = pc_class.GetOutputAnnotationLink()

if_class = ImageFlow(ds,pc_al_out)
if_class.GetRenderWindow().SetPosition(693,170)
if_class.GetRenderWindow().SetSize(600,300)
if_al_out = if_class.GetOutputAnnotationLink()

pc_class.SetHighlightAnnotationLink(if_al_out)

nf_class = DetailImageFlow(ds, if_al_out)
nf_class.GetRenderWindow().SetPosition(693,500)
nf_class.GetRenderWindow().SetSize(600,470)
nf_class.SetFlowDirection(Direction.Vertical)

def IcicleSelectionCallback(caller, event):

	annSel = caller.GetCurrentSelection()
	if annSel.GetNumberOfNodes() > 0:
		idxArr = annSel.GetNode(0).GetSelectionList()
		if idxArr.GetNumberOfTuples() > 0:
			print "if_out ", VN.vtk_to_numpy(idxArr)
		else:
			print "if back to main with selection node but no tuples"
	else:
		print "if back to main with no selection node"

# Set up callback to update 3d render window when selections are changed in 
#       parallel coordinates view
if_al_out.AddObserver("AnnotationChangedEvent", IcicleSelectionCallback)

# Only need to Start() interactor for one view
iview = pc_class.GetView()
iview.GetInteractor().Start()
