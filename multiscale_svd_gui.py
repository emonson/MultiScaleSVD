from data_source import DataSource
from icicle_noview_textured import IcicleNoView
from pcoords_chart import PCoordsChart
from image_flow import ImageFlow
from detail_image_flow import DetailImageFlow
from xy_chart import XYChart
from constants import Direction

# from tkFileDialog import askopenfilename
import vtk
import vtk.util.numpy_support as VN
import numpy as N
import os

# print os.getcwd()
from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

from ui_multiscale_svd import Ui_MainWindow

class MultiScaleSVDViews(QtGui.QMainWindow):

	def __init__(self, parent = None):

		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()

		if os.path.exists('/Users/emonson/Data/Fodava/EMoGWDataSets'):
			self.openFilesDefaultPath = '/Users/emonson/Data/Fodava/EMoGWDataSets'
		else:
			self.openFilesDefaultPath = QtCore.QDir.homePath()
		
		data_file = QtGui.QFileDialog.getOpenFileName(self,
				"Load Saved Matlab File",
				self.openFilesDefaultPath,
				"All Files (*);;Matlab Files (*.mat)")

		# DataSource loads .mat file and can generate data from it for other views
		print "Loading data from ", str(data_file)
		self.ds = DataSource(str(data_file))
		
		self.last_data_dir = os.path.dirname(str(data_file))
		
		# All view classes have access to an instance of that data source for internal queries
		# Note that the only view which will pull and display data right away is the icicle view
		#  the other views need to be able to initialize without any data and only pull and show
		#  upon the first AnnotationChanged event...

		# View #0 -- Icicle View
		self.ice_class = IcicleNoView(self.ds)
		self.ice_class.GetRenderWindow().SetPosition(50,500)
		self.ice_class.GetRenderWindow().SetSize(630,470)
		self.ice_al_out = self.ice_class.GetOutputAnnotationLink()

		# Note: With the way I've implemented the output annotation link in PCoords chart,
		#	it will always have a selection node, but the selection list may have no tuples if
		#	it's an empty selection (and event gets fired on every empty selection

		# View #1 -- PCoords View
		self.pc_class = PCoordsChart(self.ds)
		self.pc_class.SetInputAnnotationLink(self.ice_al_out)
		self.pc_al_out = self.pc_class.GetOutputAnnotationLink()
		self.pc_al = self.pc_class.GetAnnotationLink()

		# View #2 -- Image Flow View
		self.if_class = ImageFlow(self.ds, self.pc_al_out)
		self.if_al_out = self.if_class.GetOutputAnnotationLink()
		self.pc_class.SetHighlightAnnotationLink(self.if_al_out)

		# View #3 -- Detail View
		self.nf_class = DetailImageFlow(self.ds, self.if_al_out)
		self.nf_class.SetFlowDirection(Direction.Vertical)
		self.nf_al_out = self.nf_class.GetOutputAnnotationLink()

		# View #4 -- XY Chart View
		self.xy_class = XYChart(self.ds)
		self.xy_class.SetInputAnnotationLink(self.ice_al_out)
		self.xy_class.SetAnnotationLink(self.pc_al)
		self.xy_class.SetHighlightAnnotationLink(self.if_al_out)

		self.ice_class.SetGroupAnnotationLink(self.pc_al_out)
		self.ice_class.SetHighlightAnnotationLink(self.if_al_out)
		self.ice_class.SetScaleAnnotationLink(self.nf_al_out)

		# View #5 -- Axis Images (xy control) View

		# Set up all the render windows in the GUI
		self.ui.setupUi(self)

		self.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD :: Wavelets", None, QtGui.QApplication.UnicodeUTF8))

		# Now need to get all the interactors working properly
		# Icicle
		self.ui.qvtkWidget_0.SetRenderWindow(self.ice_class.GetRenderWindow())
		self.ice_class.GetInteractor().Initialize()
		# self.ui.qvtkWidget_0.show()
		# PCoords
		self.pc_class.GetView().SetInteractor(self.ui.qvtkWidget_1.GetInteractor())
		self.ui.qvtkWidget_1.SetRenderWindow(self.pc_class.GetView().GetRenderWindow())
		# self.ui.qvtkWidget_1.show()
		# Image Flow
		self.ui.qvtkWidget_2.SetRenderWindow(self.if_class.GetRenderWindow())
		self.if_class.GetInteractor().Initialize()
		# self.ui.qvtkWidget_2.show()
		# Detail Flow
		self.ui.qvtkWidget_3.SetRenderWindow(self.nf_class.GetRenderWindow())
		self.nf_class.GetInteractor().Initialize()
		# self.ui.qvtkWidget_3.show()
		# XYChart
		self.xy_class.GetChartView().SetInteractor(self.ui.qvtkWidget_4.GetInteractor())
		self.ui.qvtkWidget_4.SetRenderWindow(self.xy_class.GetChartView().GetRenderWindow())
		# self.ui.qvtkWidget_4.show()
		# AxisImages
		self.xy_class.GetAxisView().SetInteractor(self.ui.qvtkWidget_5.GetInteractor())
		self.ui.qvtkWidget_5.SetRenderWindow(self.xy_class.GetAxisView().GetRenderWindow())
		# self.ui.qvtkWidget_5.show()

		# Set sizes for veritcal splitters (left,right or top,bottom)
		self.ui.splitter_0.setSizes([280,220])
		self.ui.splitter_1.setSizes([260,240])
		self.ui.splitter_2.setSizes([130,370])
		self.ui.splitter_3.setSizes([490,280])
		self.ui.splitter_4.setSizes([330,770])




		# Deal with color_by_array menu items
		self.color_array_actions_list = []
		# Only want one array color to be set at a time
		self.colorActionGroup = QtGui.QActionGroup(self)
		# TODO: Switch this for different default
		self.ui.actionColorNone.setChecked(True)
		self.colorActionGroup.addAction(self.ui.actionColorNone)
		if self.ds.cat_labels_exist:
			for ii,label in enumerate(self.ds.label_names):
				actionTmp = QtGui.QAction(self)
				actionTmp.setCheckable(True)
				actionTmp.setObjectName(_fromUtf8(label))
				actionTmp.setText(QtGui.QApplication.translate("MainWindow", label, None, QtGui.QApplication.UnicodeUTF8))
				self.color_array_actions_list.append(actionTmp)
				self.ui.menuPlot_Colors.addAction(self.color_array_actions_list[ii])
				self.colorActionGroup.addAction(self.color_array_actions_list[ii])
				QtCore.QObject.connect(self.color_array_actions_list[ii], QtCore.SIGNAL("triggered()"), self.setColorByArray)

		# Explicitly set default here
		self.ds.SetCoeffSource('wavelets')
		# self.ds.SetCoeffSource('scaling')
		
		# TODO: Need to call menu actions to set color and coeff defaults rather than
		#   hard-coding them here...
		
		# Only want one type of coefficient to be set at a time
		basisActionGroup = QtGui.QActionGroup(self)
		self.ui.actionWavelet.setChecked(True)
		basisActionGroup.addAction(self.ui.actionWavelet)
		basisActionGroup.addAction(self.ui.actionScaling)

		# Only want one type of parallel coordinates scale range to be set at a time
		pcScaleRangeActionGroup = QtGui.QActionGroup(self)
		self.ui.actionPC_Coarsest_to_Current.setChecked(True)
		pcScaleRangeActionGroup.addAction(self.ui.actionPC_All_Scales)
		pcScaleRangeActionGroup.addAction(self.ui.actionPC_Current_Scale)
		pcScaleRangeActionGroup.addAction(self.ui.actionPC_Coarsest_to_Current)
		pcScaleRangeActionGroup.addAction(self.ui.actionPC_Current_to_Finest)

		# Connect signals and slots
		QtCore.QObject.connect(self.ui.actionExit, QtCore.SIGNAL("triggered()"), self.fileExit)
		QtCore.QObject.connect(self.ui.actionOpen, QtCore.SIGNAL("triggered()"), self.fileOpen)
		QtCore.QObject.connect(self.ui.actionWavelet, QtCore.SIGNAL("triggered()"), self.switchToWavelets)
		QtCore.QObject.connect(self.ui.actionScaling, QtCore.SIGNAL("triggered()"), self.switchToScaling)
		QtCore.QObject.connect(self.ui.actionPC_All_Scales, QtCore.SIGNAL("triggered()"), self.switchToPCAllScales)
		QtCore.QObject.connect(self.ui.actionPC_Current_Scale, QtCore.SIGNAL("triggered()"), self.switchToPCCurrentScale)
		QtCore.QObject.connect(self.ui.actionPC_Coarsest_to_Current, QtCore.SIGNAL("triggered()"), self.switchToPCCoarsestToCurrent)
		QtCore.QObject.connect(self.ui.actionPC_Current_to_Finest, QtCore.SIGNAL("triggered()"), self.switchToPCCurrentToFinest)
		QtCore.QObject.connect(self.ui.actionColorNone, QtCore.SIGNAL("triggered()"), self.setColorByArray)
		
		# NOTE: These are broken for now with QVTKWidget...
		# Trying to see whether I can pass selection bounds from xy chart to pcoords
		# self.xy_class.GetChartView().GetInteractor().AddObserver("LeftButtonReleaseEvent", self.XYSelectionReleaseCallback)
		# Trying to see whether I can pass selection bounds from xy chart to pcoords
		# self.xy_class.GetChartView().GetInteractor().AddObserver("LeftButtonPressEvent", self.XYSelectionPressCallback)

		# Testing out events for axis image modification so can have callback here
		self.xy_class.GetAxisImageItem().AddObserver("PropertyModifiedEvent", self.AIxyChangedCallback)

		# Only need to Start() interactor for one view
		# self.pc_class.GetView().GetInteractor().Start()
		# Shouldn't have to do this render...
		# for rw in self.renWinList:
		# 	rw.Render()

	def XYSelectionReleaseCallback(self, caller, event):
		x0,y0 = caller.GetEventPosition()
		print "Release (", x0, y0, ")"
		# Using this callback to get rid of parallel coordinates selection bars
		# if a selection has been made in the XY chart...
		self.pc_class.GetChart().ClearAxesSelections()
		self.pc_class.GetView().Render()

	def XYSelectionPressCallback(self, caller, event):
		x0,y0 = caller.GetEventPosition()
		print "Press (", x0, y0, ")"

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

	def AIxyChangedCallback(self, caller, event):
		xI = self.xy_class.GetAxisImageItem().GetXAxisIndex()
		yI = self.xy_class.GetAxisImageItem().GetYAxisIndex()
		# print "AI CALLBACK: (%d, %d)" % (xI,yI)
		self.xy_class.GetChartXY().SetPlotColumnIndices(xI,yI)
		self.xy_class.GetChartView().Render()
		self.pc_class.SetCurrentXY(xI,yI)
		self.pc_class.GetView().Render()

	# - - - - - - - - - - - - - - - - - - - - - -
	def switchToWavelets(self):
		if self.ds.GetCoeffSource().lower().startswith('sca'):
			self.ds.SetCoeffSource('wavelet')
			# Force other classes to reload their image and table data from new source
			# while keeping all selections...
			self.ice_class.ReloadTextureImages()
			# self.nf_class.ReloadBasisImages()
			# Instead of reloading basis images in detail view, try just firing event from image flow...
			self.if_al_out.InvokeEvent("AnnotationChangedEvent")
			# Axis images seems to reset properly on switchover, so follow that...
			xI = self.xy_class.GetAxisImageItem().GetXAxisIndex()
			yI = self.xy_class.GetAxisImageItem().GetYAxisIndex()
# 			self.xy_class.GetChartXY().SetPlotColumnIndices(xI,yI)
# 			self.xy_class.GetChartView().Render()
			self.pc_class.SetCurrentXY(xI,yI)
			self.pc_class.GetView().Render()
			self.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD :: Wavelets", None, QtGui.QApplication.UnicodeUTF8))

	def switchToScaling(self):
		if self.ds.GetCoeffSource().lower().startswith('wav'):
			self.ds.SetCoeffSource('scaling')
			# Force other classes to reload their image and table data from new source
			# while keeping all selections...
			self.ice_class.ReloadTextureImages()
			# self.nf_class.ReloadBasisImages()
			# Instead of reloading basis images in detail view, try just firing event from image flow...
			self.if_al_out.InvokeEvent("AnnotationChangedEvent")
			# Axis images seems to reset properly on switchover, so follow that...
			xI = self.xy_class.GetAxisImageItem().GetXAxisIndex()
			yI = self.xy_class.GetAxisImageItem().GetYAxisIndex()
# 			self.xy_class.GetChartXY().SetPlotColumnIndices(xI,yI)
# 			self.xy_class.GetChartView().Render()
			self.pc_class.SetCurrentXY(xI,yI)
			self.pc_class.GetView().Render()
			self.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD :: Scaling Functions", None, QtGui.QApplication.UnicodeUTF8))

	def switchToPCAllScales(self):
		self.pc_class.SetScaleRange('all')
		self.pc_class.UpdateChartWithCurrentData()

	def switchToPCCurrentScale(self):
		self.pc_class.SetScaleRange('current')
		self.pc_class.UpdateChartWithCurrentData()

	def switchToPCCoarsestToCurrent(self):
		self.pc_class.SetScaleRange('coarse')
		self.pc_class.UpdateChartWithCurrentData()

	def switchToPCCurrentToFinest(self):
		self.pc_class.SetScaleRange('fine')
		self.pc_class.UpdateChartWithCurrentData()

	def setColorByArray(self):
		sender = self.sender()
		
		if sender in self.color_array_actions_list:
			self.pc_class.SetColorByArray(str(sender.text()))
			self.pc_class.GetView().Render()
			self.xy_class.SetColorByArray(str(sender.text()))
			self.xy_class.GetChartView().Render()
		else:
			self.pc_class.SetColorByArrayOff()
			self.pc_class.GetView().Render()
			self.xy_class.SetColorByArrayOff()
			self.xy_class.GetChartView().Render()
		
	def generate_color_array_actions(self):
		self.color_array_actions_list = []
		# Only want one array color to be set at a time
		self.colorActionGroup = QtGui.QActionGroup(self)
		# TODO: Switch this for different default
		self.ui.actionColorNone.setChecked(True)
		self.colorActionGroup.addAction(self.ui.actionColorNone)
		if self.ds.cat_labels_exist:
			for ii,label in enumerate(self.ds.label_names):
				actionTmp = QtGui.QAction(self)
				actionTmp.setCheckable(True)
				actionTmp.setObjectName(_fromUtf8(label))
				actionTmp.setText(QtGui.QApplication.translate("MainWindow", label, None, QtGui.QApplication.UnicodeUTF8))
				self.color_array_actions_list.append(actionTmp)
				self.ui.menuPlot_Colors.addAction(self.color_array_actions_list[ii])
				self.colorActionGroup.addAction(self.color_array_actions_list[ii])
				QtCore.QObject.connect(self.color_array_actions_list[ii], QtCore.SIGNAL("triggered()"), self.setColorByArray)
		
	# - - - - - - - - - - - - - - - - - - - - - -
	def fileOpen(self):

		openFilesDefaultPath = ''

		# This dialog allows multiple files to be selected at once
		# If only want a single file, change to fileName = QtGui.QFileDialog.getOpenFileName(...)
		# Just change the string in the next to last line of the method for different file types
		file = QtGui.QFileDialog.getOpenFileName(self,
				"Load Saved Matlab File",
				self.last_data_dir,
				"All Files (*);;Matlab Files (*.mat)")

		if file:
			# Clear out selections
			empty_arr = N.array([],dtype='int64')
			empty_vtk = VN.numpy_to_vtkIdTypeArray(empty_arr, deep=True)
			self.pc_class.GetAnnotationLink().GetCurrentSelection().GetNode(0).SetSelectionList(empty_vtk)
			self.pc_class.GetAnnotationLink().InvokeEvent("AnnotationChangedEvent")
			self.if_class.GetOutputAnnotationLink().GetCurrentSelection().GetNode(0).SetSelectionList(empty_vtk)
			self.if_class.GetOutputAnnotationLink().InvokeEvent("AnnotationChangedEvent")
			self.nf_class.GetOutputAnnotationLink().GetCurrentSelection().GetNode(0).SetSelectionList(empty_vtk)
			self.nf_class.GetOutputAnnotationLink().InvokeEvent("AnnotationChangedEvent")

			self.ds.SetFileName(str(file))
			self.ds.LoadData()
			
			# Remove old color_by_array menu items
			for action in self.color_array_actions_list:
				self.ui.menuPlot_Colors.removeAction(action)
				self.colorActionGroup.removeAction(action)
				QtCore.QObject.connect(action, QtCore.SIGNAL("triggered()"), self.setColorByArray)
			
			# Add new color_by_array menu items
			self.generate_color_array_actions()
			
			self.ice_class.LoadData()
			self.ui.qvtkWidget_0.update()
			self.ui.qvtkWidget_1.update()

	def fileExit(self):

		# Usually would use the qApp global variable qApp.quit(), but wasn't working...
		QtGui.QApplication.instance().quit()

