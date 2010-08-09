from data_source import DataSource
from icicle_noview_textured import IcicleNoView

import vtk
from PyQt4 import QtCore, QtGui

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(400, 400)
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.gridlayout = QtGui.QGridLayout(self.centralWidget)
		self.vtkWidget = vtk.QVTKWidget(self.centralWidget)
		self.gridlayout.addWidget(self.vtkWidget, 0, 0, 1, 1)
		MainWindow.setCentralWidget(self.centralWidget)

class MultiScaleSVDViews(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
						   
		data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist1_5c_20100324.mat'
		
		# DataSource loads .mat file and can generate data from it for other views
		self.ds = DataSource(str(data_file))
		
		# View #0 -- Icicle View
		self.ice_class = IcicleNoView(self.ds)
		rw = self.ice_class.GetRenderWindow()
		
		# Now need to get all of the qvtkwidgets working properly
		# Icicle
		rw.SetInteractor(self.ui.vtkWidget.GetInteractor())
		self.ui.vtkWidget.SetRenderWindow(rw)
		istyle = vtk.vtkInteractorStyleImage()
		self.ice_class.SetInteractorStyle(istyle)

		
