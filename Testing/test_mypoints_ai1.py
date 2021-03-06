# Translated to Python from [VTK]/Charts/Testing/Cxx/TestLinePlot.cxx

# This version is for testing reworked subclasses 8/13/2010

import vtk
from vtk.util import numpy_support as VN
import numpy as N
import math
import vtkvtg
from data_source import DataSource

data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist1_5c_20100324.mat'

# DataSource loads .mat file and can generate data from it for other views
ds = DataSource(data_file)


renwin = vtk.vtkRenderWindow()
renwin.SetMultiSamples(0)
renwin.SetSize(400, 600)

iren = vtk.vtkRenderWindowInteractor()
istyle = vtk.vtkInteractorStyleRubberBand2D()
iren.SetInteractorStyle(istyle)
iren.SetRenderWindow(renwin)

# setup the 2 view ports
viewports = [(0.0,0.0,1.0,0.5), (0.0,0.5,1.0,1.0)]

# CHART
chartRen = vtk.vtkRenderer()
chartRen.SetBackground(1.0,1.0,1.0)
chartRen.SetViewport(viewports[1])
renwin.AddRenderer(chartRen)

# Testing my custom chart class which has image hover tooltips
chart = vtkvtg.vtkMyChartXY()
chartScene = vtk.vtkContextScene()
chartActor = vtk.vtkContextActor()

chartScene.AddItem(chart)
chartActor.SetScene(chartScene)

# both needed
chartRen.AddActor(chartActor)
chartScene.SetRenderer(chartRen)
chartScene.SetInteractorStyle(istyle)

# AXIS IMAGES
aiRen = vtk.vtkRenderer()
aiRen.SetBackground(1.0,1.0,1.0)
aiRen.SetViewport(viewports[0])
renwin.AddRenderer(aiRen)

# Testing my custom chart class which has image hover tooltips
ai = vtkvtg.vtkAxisImageItem()
aiScene = vtk.vtkContextScene()
aiActor = vtk.vtkContextActor()

aiScene.AddItem(ai)
aiActor.SetScene(aiScene)

# both needed
aiRen.AddActor(aiActor)
aiScene.SetRenderer(aiRen)
aiScene.SetInteractorStyle(istyle)


# Create a annotation link to access selection in parallel coordinates view
annotationLink = vtk.vtkAnnotationLink()
# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
# See vtkSelectionNode doc for field and content type enum values
annotationLink.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
annotationLink.GetCurrentSelection().GetNode(0).SetContentType(4)   # Indices
# Connect the annotation link to the parallel coordinates representation
chart.SetAnnotationLink(annotationLink)

test_id = 68
table = ds.GetNodeOneScaleCoeffTable(test_id)

line1 = vtkvtg.vtkMyPlotPoints()
chart.AddPlot(line1)		# POINTS
line1.SetInput(table, 0, 1)
line1.SetMarkerStyle(2)
line1.SetColor(0, 0, 0, 255)


# Tooltip image stack will now be owned by the tooltip, so need to do that differently... 
id_list = ds.PIN[test_id]
image_stack = ds.GetProjectedImages(id_list)
chart.SetTooltipImageStack(image_stack)
chart.SetTooltipShowImage(True)
# chart.SetTooltipImageScalingFactor(2.0)
chart.SetTooltipImageTargetSize(40)

axis_images = ds.GetNodeBasisImages(test_id)
center_image = ds.GetNodeCenterImage(test_id)
ai.SetAxisImagesHorizontal()
ai.SetChartXY(chart)
ai.SetAxisImageStack(axis_images)
ai.SetCenterImage(center_image)

# Set up annotation link which will carry indices to parallel coordinates chart
# for highlighting outside selections (e.g. back from image_flow)
# This needs to carry indices, while image_flow link outputs pedigree ids
# so conversion happens in HighlightSelectionCallback
highlight_link_idxs = vtk.vtkAnnotationLink()
highlight_link_idxs.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
highlight_link_idxs.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
chart.SetHighlightLink(highlight_link_idxs)

# Finally render the scene and compare the image to a reference image
# view.GetRenderWindow().SetMultiSamples(0)

def selectionCallback(caller, event):
        annSel = annotationLink.GetCurrentSelection()
        if annSel.GetNumberOfNodes() > 0:
                idxArr = annSel.GetNode(0).GetSelectionList()
                if idxArr.GetNumberOfTuples() > 0:
                        print VN.vtk_to_numpy(idxArr)


annotationLink.AddObserver("AnnotationChangedEvent",selectionCallback)

# view.ResetCamera()
renwin.Render()

# Fill selection link with dummy IDs
id_array = N.array([0],dtype='int64')
id_list = VN.numpy_to_vtkIdTypeArray(id_array)
highlight_link_idxs.GetCurrentSelection().GetNode(0).SetSelectionList(id_list)
highlight_link_idxs.InvokeEvent("AnnotationChangedEvent")



# Start interaction event loop
renwin.GetInteractor().Start()
