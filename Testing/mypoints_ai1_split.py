from data_source import DataSource

# from tkFileDialog import askopenfilename
import vtk
import vtk.util.numpy_support as VN
import numpy as N
import os
import vtkvtg

# print os.getcwd()

from PyQt4 import QtCore, QtGui
from ui_mypoints_ai1 import Ui_MainWindow

class MultiScaleSVDViews(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		
		self.renWinList = []
		   
		# data_file = askopenfilename()
		data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist1_5c_20100324.mat'
# 		self.openFilesDefaultPath = QtCore.QDir.homePath()
# 		data_file = QtGui.QFileDialog.getOpenFileName(self,
# 				"Load Saved Matlab File", 
# 				self.openFilesDefaultPath,
# 				"All Files (*);;Matlab Files (*.mat)")
		
		# DataSource loads .mat file and can generate data from it for other views
		self.ds = DataSource(str(data_file))
		
		# All view classes have access to an instance of that data source for internal queries
		# Note that the only view which will pull and display data right away is the icicle view
		#  the other views need to be able to initialize without any data and only pull and show
		#  upon the first AnnotationChanged event...
		
		# View #0 -- Icicle View
		# Set up a 2D scene, add an XY chart to it
		self.chartView = vtk.vtkContextView()
		self.chartView.GetRenderer().SetBackground(1.0, 1.0, 1.0)
		self.renWinList.append(self.chartView.GetRenderWindow())

		self.chart = vtkvtg.vtkMyChartXY()
		self.chartView.GetScene().AddItem(self.chart)

		# View #1 -- AxisImageView
		self.axisView = vtk.vtkContextView()
		self.axisView.GetRenderer().SetBackground(1.0, 1.0, 1.0)
		self.renWinList.append(self.axisView.GetRenderWindow())

		self.ai = vtkvtg.vtkAxisImageItem()
		self.axisView.GetScene().AddItem(self.ai)
		
		# Set up all the render windows in the GUI
		self.ui.setupUi(self, self.renWinList)
		
		# Now need to get all the interactors working properly
		# XY
		style0 = vtk.vtkInteractorStyleRubberBand2D()
		self.chartView.GetInteractor().SetInteractorStyle(style0)
		self.chartView.GetScene().SetInteractorStyle(style0)

		# Axis images
		style1 = vtk.vtkInteractorStyleRubberBand2D()
		self.axisView.GetInteractor().SetInteractorStyle(style1)
		self.axisView.GetScene().SetInteractorStyle(style1)

		# Set sizes for veritcal splitters
		self.ui.splitter.setSizes([200,400])		
		
		# Connect signals and slots
		QtCore.QObject.connect(self.ui.actionExit, QtCore.SIGNAL("triggered()"), self.fileExit)

		
		
		# CORE setting up chart and axis images
		test_id = 68
		self.table = self.ds.GetNodeOneScaleCoeffTable(test_id)
		
		line1 = vtkvtg.vtkMyPlotPoints()
		self.chart.AddPlot(line1)		# POINTS
		line1.SetInput(self.table, 0, 1)
		line1.SetMarkerStyle(2)
		line1.SetColor(0, 0, 0, 255)
		
		
		# Tooltip image stack will now be owned by the tooltip, so need to do that differently... 
		id_list = self.ds.PIN[test_id]
		image_stack = self.ds.GetProjectedImages(id_list)
		self.chart.SetTooltipImageStack(image_stack)
		self.chart.SetTooltipShowImage(True)
		# self.chart.SetTooltipImageScalingFactor(2.0)
		self.chart.SetTooltipImageTargetSize(40)
		
		axis_images = self.ds.GetNodeBasisImages(test_id)
		center_image = self.ds.GetNodeCenterImage(test_id)
		self.ai.SetAxisImagesHorizontal()
		self.ai.SetChartXY(self.chart)
		self.ai.SetAxisImageStack(axis_images)
		self.ai.SetCenterImage(center_image)
		
		# Set up annotation link which will carry indices to parallel coordinates chart
		# for highlighting outside selections (e.g. back from image_flow)
		# This needs to carry indices, while image_flow link outputs pedigree ids
		# so conversion happens in HighlightSelectionCallback
		data_col_idxs = vtk.vtkAnnotationLink()
		data_col_idxs.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		data_col_idxs.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
		self.chart.SetDataColumnsLink(data_col_idxs)
		self.ai.SetDataColumnsLink(data_col_idxs)

		# Create a annotation link to access selection in parallel coordinates view
		annotationLink = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		annotationLink.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		annotationLink.GetCurrentSelection().GetNode(0).SetContentType(4)   # Indices
		# Connect the annotation link to the parallel coordinates representation
		self.chart.SetAnnotationLink(annotationLink)
		self.chart.GetAnnotationLink().AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)

		# Set up annotation link which will carry indices to parallel coordinates chart
		# for highlighting outside selections (e.g. back from image_flow)
		# This needs to carry indices, while image_flow link outputs pedigree ids
		# so conversion happens in HighlightSelectionCallback
		highlight_link_idxs = vtk.vtkAnnotationLink()
		highlight_link_idxs.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		highlight_link_idxs.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
		self.chart.SetHighlightLink(highlight_link_idxs)
		
		# Fill selection link with dummy IDs
		id_array = N.array([0],dtype='int64')
		id_list = VN.numpy_to_vtkIdTypeArray(id_array)
		highlight_link_idxs.GetCurrentSelection().GetNode(0).SetSelectionList(id_list)
		highlight_link_idxs.InvokeEvent("AnnotationChangedEvent")

		self.updater = vtk.vtkViewUpdater()
		self.updater.AddAnnotationLink(data_col_idxs)
		self.updater.AddView(self.axisView)
		self.updater.AddView(self.chartView)
		
		col_array = N.array([0,1],dtype='int64')
		col_vtk = VN.numpy_to_vtkIdTypeArray(col_array, deep=True)
		data_col_idxs.GetCurrentSelection().GetNode(0).SetSelectionList(col_vtk)
		data_col_idxs.InvokeEvent("AnnotationChangedEvent")
		
		
		# Only need to Start() interactor for one view
		# self.pc_class.GetView().GetInteractor().Start()
		# Shouldn't have to do this render...
		for rw in self.renWinList:
			rw.Render()

	def IcicleSelectionCallback(self, caller, event):
	
		annSel = caller.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			if idxArr.GetNumberOfTuples() > 0:
				print "if_out ", VN.vtk_to_numpy(idxArr)
			else:
				print "if back to main with selection node but no tuples"
		else:
			print "if back to main with no selection node"
	

	def fileExit(self):
		
		# Usually would use the qApp global variable qApp.quit(), but wasn't working...
		QtGui.QApplication.instance().quit()
		

from PyQt4.QtGui import QApplication
import sys

if __name__ == "__main__":

    app = QApplication(sys.argv)
    window = MultiScaleSVDViews()
    window.show()
    sys.exit(app.exec_())
