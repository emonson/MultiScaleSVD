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
from ui_multiscale_svd import Ui_MainWindow

class MultiScaleSVDViews(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		
		self.renWinList = []
		   
		# data_file = askopenfilename()
		data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100909.mat'
# 		self.openFilesDefaultPath = QtCore.QDir.homePath()
# 		data_file = QtGui.QFileDialog.getOpenFileName(self,
# 				"Load Saved Matlab File", 
# 				self.openFilesDefaultPath,
# 				"All Files (*);;Matlab Files (*.mat)")
		
		# DataSource loads .mat file and can generate data from it for other views
		print "Loading data from ", str(data_file)
		self.ds = DataSource(str(data_file))
		
		# All view classes have access to an instance of that data source for internal queries
		# Note that the only view which will pull and display data right away is the icicle view
		#  the other views need to be able to initialize without any data and only pull and show
		#  upon the first AnnotationChanged event...
		
		# View #0 -- Icicle View
		self.ice_class = IcicleNoView(self.ds)
		self.ice_class.GetRenderWindow().SetPosition(50,500)
		self.ice_class.GetRenderWindow().SetSize(630,470)
		self.ice_al_out = self.ice_class.GetOutputAnnotationLink()
		
		self.renWinList.append(self.ice_class.GetRenderWindow())
		
		# Note: With the way I've implemented the output annotation link in PCoords chart, 
		#   it will always have a selection node, but the selection list may have no tuples if
		#	it's an empty selection (and event gets fired on every empty selection
		
		# View #1 -- PCoords View
		self.pc_class = PCoordsChart(self.ds)
		self.pc_class.SetInputAnnotationLink(self.ice_al_out)
		self.pc_al_out = self.pc_class.GetOutputAnnotationLink()
		self.pc_al = self.pc_class.GetAnnotationLink()
		
		self.renWinList.append(self.pc_class.GetView().GetRenderWindow())

		# View #2 -- Image Flow View
		self.if_class = ImageFlow(self.ds, self.pc_al_out)
		self.if_al_out = self.if_class.GetOutputAnnotationLink()
		self.pc_class.SetHighlightAnnotationLink(self.if_al_out)
		
		self.renWinList.append(self.if_class.GetRenderWindow())
		
		# View #3 -- Detail View
		self.nf_class = DetailImageFlow(self.ds, self.if_al_out)
		self.nf_class.SetFlowDirection(Direction.Vertical)
		self.nf_al_out = self.nf_class.GetOutputAnnotationLink()
		
		self.renWinList.append(self.nf_class.GetRenderWindow())

		# View #4 -- XY Chart View
		self.xy_class = XYChart(self.ds)
		self.xy_class.SetInputAnnotationLink(self.ice_al_out)
		self.xy_class.SetAnnotationLink(self.pc_al)
		self.xy_class.SetHighlightAnnotationLink(self.if_al_out)
		
		self.ice_class.SetGroupAnnotationLink(self.pc_al_out)
		self.ice_class.SetHighlightAnnotationLink(self.if_al_out)
		self.ice_class.SetScaleAnnotationLink(self.nf_al_out)
		
		self.renWinList.append(self.xy_class.GetChartView().GetRenderWindow())
		
		# View #5 -- Axis Images (xy control) View
		self.renWinList.append(self.xy_class.GetAxisView().GetRenderWindow())
				
		# Set up all the render windows in the GUI
		self.ui.setupUi(self, self.renWinList)
		
		# Explicitly set default here
		self.ds.SetCoeffSource('wavelets')
		# self.ds.SetCoeffSource('scaling')
		self.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD :: Wavelets", None, QtGui.QApplication.UnicodeUTF8))
		
		# Now need to get all the interactors working properly
		# Icicle
		self.style0 = vtk.vtkInteractorStyleImage()
		self.ice_class.SetInteractorStyle(self.style0)
		# PCoords
		style1 = vtk.vtkInteractorStyleRubberBand2D()
		self.pc_class.GetView().GetInteractor().SetInteractorStyle(style1)
		self.pc_class.GetView().GetScene().SetInteractorStyle(style1)
		# Image Flow
		self.style2 = vtk.vtkInteractorStyleImage()
		self.if_class.SetInteractorStyle(self.style2)		
		# Detail Flow
		self.style3 = vtk.vtkInteractorStyleImage()
		self.nf_class.SetInteractorStyle(self.style3)		
		# XYChart
		self.style4 = vtk.vtkInteractorStyleRubberBand2D()
		self.xy_class.GetChartView().GetInteractor().SetInteractorStyle(self.style4)
		self.xy_class.GetChartView().GetScene().SetInteractorStyle(self.style4)
		# AxisImages
		style5 = vtk.vtkInteractorStyleRubberBand2D()
		self.xy_class.GetAxisView().GetInteractor().SetInteractorStyle(style5)
		self.xy_class.GetAxisView().GetScene().SetInteractorStyle(style5)		

		# Set sizes for veritcal splitters (left,right or top,bottom)
		self.ui.splitter_0.setSizes([280,220])		
		self.ui.splitter_1.setSizes([260,240])
		self.ui.splitter_2.setSizes([130,370])
		self.ui.splitter_3.setSizes([490,280])
		self.ui.splitter_4.setSizes([330,770])
		
		# Only want one type of coefficient to be set at a time
		basisActionGroup = QtGui.QActionGroup(self)
		self.ui.actionWavelet.setChecked(True)
		basisActionGroup.addAction(self.ui.actionWavelet)
		basisActionGroup.addAction(self.ui.actionScaling)
		
		# Connect signals and slots
		QtCore.QObject.connect(self.ui.actionExit, QtCore.SIGNAL("triggered()"), self.fileExit)
		QtCore.QObject.connect(self.ui.actionOpen, QtCore.SIGNAL("triggered()"), self.fileOpen)
		QtCore.QObject.connect(self.ui.actionWavelet, QtCore.SIGNAL("triggered()"), self.switchToWavelets)
		QtCore.QObject.connect(self.ui.actionScaling, QtCore.SIGNAL("triggered()"), self.switchToScaling)
		
		# Trying to see whether I can pass selection bounds from xy chart to pcoords
		self.xy_class.GetChartView().GetInteractor().AddObserver("LeftButtonReleaseEvent", self.XYSelectionReleaseCallback)
		# Trying to see whether I can pass selection bounds from xy chart to pcoords
		self.xy_class.GetChartView().GetInteractor().AddObserver("LeftButtonPressEvent", self.XYSelectionPressCallback)
		
		# Testing out events for axis image modification so can have callback here
		self.xy_class.GetAxisImageItem().AddObserver("PropertyModifiedEvent", self.AIxyChangedCallback)

		# Only need to Start() interactor for one view
		# self.pc_class.GetView().GetInteractor().Start()
		# Shouldn't have to do this render...
		for rw in self.renWinList:
			rw.Render()

	def XYSelectionReleaseCallback(self, caller, event):
		x0,y0 = caller.GetEventPosition()
		print "Release (", x0, y0, ")"
	
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
		print "AI CALLBACK: (%d, %d)" % (xI,yI)
		self.xy_class.GetChartXY().SetPlotColumnIndices(xI,yI)
		self.xy_class.GetChartView().Render()
	
	# - - - - - - - - - - - - - - - - - - - - - -
	# TODO: Axis images by default switch back to 0,1 xy...
	def switchToWavelets(self):
		if self.ds.GetCoeffSource().lower().startswith('sca'):
			# Record which xy indices are being plotted to reset after switchover
			xI = self.xy_class.GetAxisImageItem().GetXAxisIndex()
			yI = self.xy_class.GetAxisImageItem().GetYAxisIndex()
			self.ds.SetCoeffSource('wavelet')
			# Force other classes to reload their image and table data from new source
			# while keeping all selections...
			self.ice_class.ReloadTextureImages()
			self.nf_class.ReloadBasisImages()
			# self.xy_class.GetAxisImageItem().SetAxisIndices(xI,yI)
			# self.xy_class.GetChartXY().SetPlotColumnIndices(xI,yI)
			# self.xy_class.GetAxisView().Render()
			# self.xy_class.GetChartView().Render()
			self.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD :: Wavelets", None, QtGui.QApplication.UnicodeUTF8))

	def switchToScaling(self):
		if self.ds.GetCoeffSource().lower().startswith('wav'):
			# Record which xy indices are being plotted to reset after switchover
			xI = self.xy_class.GetAxisImageItem().GetXAxisIndex()
			yI = self.xy_class.GetAxisImageItem().GetYAxisIndex()
			self.ds.SetCoeffSource('scaling')
			# Force other classes to reload their image and table data from new source
			# while keeping all selections...
			self.ice_class.ReloadTextureImages()
			self.nf_class.ReloadBasisImages()
			# self.xy_class.GetAxisImageItem().SetAxisIndices(xI,yI)
			# self.xy_class.GetChartXY().SetPlotColumnIndices(xI,yI)
			# self.xy_class.GetAxisView().Render()
			# self.xy_class.GetChartView().Render()
			self.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD :: Scaling Functions", None, QtGui.QApplication.UnicodeUTF8))

	# - - - - - - - - - - - - - - - - - - - - - -
	def fileOpen(self):
	
		openFilesDefaultPath = ''

		# This dialog allows multiple files to be selected at once
		# If only want a single file, change to fileName = QtGui.QFileDialog.getOpenFileName(...)
		# Just change the string in the next to last line of the method for different file types
		file = QtGui.QFileDialog.getOpenFileName(self,
				"Load Saved Matlab File", 
				openFilesDefaultPath,
				"All Files (*);;Matlab Files (*.mat)")
		
		if file:
			self.ds.SetFileName(str(file))
			self.ds.LoadData()
			self.ice_class.LoadData()

	def fileExit(self):
		
		# Usually would use the qApp global variable qApp.quit(), but wasn't working...
		QtGui.QApplication.instance().quit()
		
