# Translated to Python from [VTK]/Charts/Testing/Cxx/TestLinePlot.cxx

import vtk
from vtk.util import numpy_support as VN
import math
import vtkvtg
from data_source import DataSource

data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/olivetti_all_20100602.mat'

# DataSource loads .mat file and can generate data from it for other views
ds = DataSource(data_file)

# Set up a 2D scene, add an XY chart to it
view = vtk.vtkContextView()
view.GetRenderWindow().SetSize(400, 300)

# Testing my custom chart class which has image hover tooltips
chart = vtkvtg.vtkMyChartXY()
view.GetScene().AddItem(chart)

# Create a annotation link to access selection in parallel coordinates view
annotationLink = vtk.vtkAnnotationLink()
# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
# See vtkSelectionNode doc for field and content type enum values
annotationLink.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
annotationLink.GetCurrentSelection().GetNode(0).SetContentType(4)   # Indices
# Connect the annotation link to the parallel coordinates representation
chart.SetAnnotationLink(annotationLink)

test_id = 23
table = ds.GetNodeOneScaleCoeffTable(test_id)
id_list = ds.PIN[test_id]
image_stack = ds.GetProjectedImages(id_list)
axis_images = ds.GetNodeBasisImages(test_id)

line1 = chart.AddPlot(1)		# POINTS
line1.SetInput(table, 0, 1)
line1.SetMarkerStyle(2)
line1.SetColor(255, 0, 0, 255)

# Need to set the image stack for the plot which will get resliced 
chart.GetPlot(0).SetTooltipImageStack(image_stack)
chart.SetTooltipShowImage(True)
chart.SetTooltipImageScalingFactor(2.0)
chart.SetAxisImageStack(axis_images)

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

# Start interaction event loop
view.GetInteractor().Start()
