# Translated to Python from [VTK]/Charts/Testing/Cxx/TestLinePlot.cxx

# This version is for testing reworked subclasses 8/13/2010

import vtk
from vtk.util import numpy_support as VN
import numpy as N
import math
import vtkvtg
from data_source import DataSource

data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20101119.mat'

# DataSource loads .mat file and can generate data from it for other views
ds = DataSource(data_file)

# Set up a 2D scene, add an XY chart to it
view = vtk.vtkContextView()
view.GetRenderWindow().SetSize(400, 300)

# Testing my custom chart class which has image hover tooltips
chart = vtkvtg.vtkMyChartXY()
chart.SetActionToButton(vtk.vtkChart.PAN, 2)
chart.SetActionToButton(vtk.vtkChart.ZOOM, 4)
chart.SetActionToButton(vtk.vtkChart.SELECT, 1)
view.GetScene().AddItem(chart)

# Create a annotation link to access selection in parallel coordinates view
annotationLink = vtk.vtkAnnotationLink()
# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
# See vtkSelectionNode doc for field and content type enum values
annotationLink.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
annotationLink.GetCurrentSelection().GetNode(0).SetContentType(4)   # Indices
# Connect the annotation link to the parallel coordinates representation
chart.SetAnnotationLink(annotationLink)

test_id = 3
table = ds.GetNodeOneScaleCoeffTable(test_id)

chart.ClearPlots()

line1 = vtkvtg.vtkMyPlotPoints()
chart.AddPlot(line1)		# POINTS
line1.SetInput(table, 0, 1)
line1.SetMarkerStyle(2)
line1.SetColor(0, 0, 0, 255)

# Tooltip image stack will now be owned by the tooltip, so need to do that differently... 
id_list = ds.PointsInNet[test_id]
image_stack = ds.GetProjectedImages(id_list)

# DEBUG
writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName('out.vti')
writer.SetInput(image_stack)
writer.Write()

chart.SetTooltipImageStack(image_stack)
chart.SetTooltipShowImage(True)
# chart.SetTooltipImageScalingFactor(2.0)
chart.SetTooltipImageTargetSize(40)

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
view.Render()

# Fill selection link with dummy IDs
id_array = N.array([0],dtype='int64')
id_list = VN.numpy_to_vtkIdTypeArray(id_array)
highlight_link_idxs.GetCurrentSelection().GetNode(0).SetSelectionList(id_list)
highlight_link_idxs.InvokeEvent("AnnotationChangedEvent")

# Set up annotation link which will carry indices to parallel coordinates chart
# for highlighting outside selections (e.g. back from image_flow)
# This needs to carry indices, while image_flow link outputs pedigree ids
# so conversion happens in HighlightSelectionCallback
data_col_idxs = vtk.vtkAnnotationLink()
data_col_idxs.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
data_col_idxs.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
chart.SetDataColumnsLink(data_col_idxs)
# Fill selection link with dummy IDs
col_array = N.array([1,2],dtype='int64')
col_list = VN.numpy_to_vtkIdTypeArray(col_array)
data_col_idxs.GetCurrentSelection().GetNode(0).SetSelectionList(col_list)
data_col_idxs.InvokeEvent("AnnotationChangedEvent")


# Start interaction event loop
view.GetInteractor().Start()
