from data_source import DataSource
from icicle_noview_textured import IcicleNoView
from pcoords_chart import PCoordsChart
from image_flow import ImageFlow
from detail_image_flow import DetailImageFlow
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
		data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100521.mat'
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
		
		# View #1 -- Icicle View
		self.ice_class = IcicleNoView(self.ds)
		self.ice_class.GetRenderWindow().SetPosition(50,500)
		self.ice_class.GetRenderWindow().SetSize(630,470)
		self.ice_al_out = self.ice_class.GetOutputAnnotationLink()
		
		self.renWinList.append(self.ice_class.GetRenderWindow())
		
		# Note: With the way I've implemented the output annotation link in PCoords chart, 
		#   it will always have a selection node, but the selection list may have no tuples if
		#	it's an empty selection (and event gets fired on every empty selection
		
		# View #2 -- PCoords View
		self.pc_class = PCoordsChart(self.ds)
		self.pc_class.SetInputAnnotationLink(self.ice_al_out)
		self.pc_class.GetView().GetRenderWindow().SetPosition(50,170)
		self.pc_class.GetView().GetRenderWindow().SetSize(630,300)
		self.pc_al_out = self.pc_class.GetOutputAnnotationLink()
		
		self.renWinList.append(self.pc_class.GetView().GetRenderWindow())

		# View #3 -- Image Flow View
		self.if_class = ImageFlow(self.ds, self.pc_al_out)
		self.if_class.GetRenderWindow().SetPosition(693,170)
		self.if_class.GetRenderWindow().SetSize(600,300)
		self.if_al_out = self.if_class.GetOutputAnnotationLink()
		
		self.renWinList.append(self.if_class.GetRenderWindow())
		
		# View #4 -- Detail View
		self.nf_class = DetailImageFlow(self.ds, self.if_al_out)
		self.nf_class.GetRenderWindow().SetPosition(693,500)
		self.nf_class.GetRenderWindow().SetSize(600,470)
		self.nf_class.SetFlowDirection(Direction.Vertical)
		
		self.renWinList.append(self.nf_class.GetRenderWindow())
		self.renWinList.append(None)
		
		self.pc_class.SetHighlightAnnotationLink(self.if_al_out)
		
		# Set up callback to update 3d render window when selections are changed in 
		#       parallel coordinates view
		self.if_al_out.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)
		
		# Set up all the render windows in the GUI
		self.ui.setupUi(self, self.renWinList)
		
		# Now need to get all the interactors working properly
		# Icicle
		self.style0 = vtk.vtkInteractorStyleImage()
		self.ice_class.SetInteractorStyle(self.style0)
		# PCoords
		style1 = vtk.vtkInteractorStyleRubberBand2D()
		self.pc_class.GetView().SetInteractorStyle(style1)
		self.pc_class.GetView().SetInteractionModeTo2D()
		# Image Flow
		self.style2 = vtk.vtkInteractorStyleImage()
		self.if_class.SetInteractorStyle(self.style2)		
		# Detail Flow
		self.style3 = vtk.vtkInteractorStyleImage()
		self.nf_class.SetInteractorStyle(self.style3)		

		# Set sizes for veritcal splitters
		self.ui.splitter_0.setSizes([360,240])		
		self.ui.splitter_1.setSizes([360,240])		
		
		# Connect signals and slots
		QtCore.QObject.connect(self.ui.actionExit, QtCore.SIGNAL("triggered()"), self.fileExit)
		QtCore.QObject.connect(self.ui.actionOpen, QtCore.SIGNAL("triggered()"), self.fileOpen)

		# Only need to Start() interactor for one view
		self.pc_class.GetView().GetInteractor().Start()

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
		
