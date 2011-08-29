"""
simple_graph_view.py

Example using PyQt4 & VTK Infovis

10 Sept 2009 -- E Monson

"""

from PyQt4 import QtCore, QtGui
import vtk
from vtk import QVTKWidget
import sys
import math
import numpy as N
import vtk.util.numpy_support as VN
import vtkvtg
from data_source import DataSource

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(600, 400)
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.gridlayout = QtGui.QGridLayout(self.centralWidget)
		self.vtkWidget = QVTKWidget(self.centralWidget)
		self.gridlayout.addWidget(self.vtkWidget, 0, 0, 1, 1)
		MainWindow.setCentralWidget(self.centralWidget)

class SimpleView(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		   
		# Set up a 2D scene (later we'll add a pcoords chart to it)
		self.view = vtk.vtkContextView()
		# self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
		
		# // QVTKWidget *widget = new QVTKWidget;
		# // vtkContextView *view = vtkContextView::New();
		# // view->SetInteractor(widget->GetInteractor());
		# // widget->SetRenderWindow(view->GetRenderWindow());

		# Necessary for RenderView types
		self.view.SetInteractor(self.ui.vtkWidget.GetInteractor())
		self.ui.vtkWidget.SetRenderWindow(self.view.GetRenderWindow())

		data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20101119.mat'
		
		# DataSource loads .mat file and can generate data from it for other views
		ds = DataSource(data_file)

		# Testing my custom chart class which has image hover tooltips
		chart = vtkvtg.vtkMyChartXY()
		chart.SetActionToButton(vtk.vtkChart.PAN, 2)
		chart.SetActionToButton(vtk.vtkChart.ZOOM, 4)
		chart.SetActionToButton(vtk.vtkChart.SELECT, 1)
		self.view.GetScene().AddItem(chart)
		
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
		
		annotationLink.AddObserver("AnnotationChangedEvent", self.selectionCallback)
		
		# view.ResetCamera()
		# view.Render()
		
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
		self.ui.vtkWidget.show()
		
	def selectionCallback(self,caller, event):
		annSel = caller.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			if idxArr.GetNumberOfTuples() > 0:
				print VN.vtk_to_numpy(idxArr)
		
if __name__ == "__main__":

	app = QtGui.QApplication(sys.argv)
	window = SimpleView()
	window.show()
	sys.exit(app.exec_())
