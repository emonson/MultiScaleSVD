"""
simple_graph_view.py

Example using PyQt4 & VTK Infovis

10 Sept 2009 -- E Monson

"""

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import sys
import math
import vtk.util.numpy_support as VN

from pcoords_chart import PCoordsChart
from data_source import DataSource

class Ui_MainWindow(object):
	def setupUi(self, MainWindow,renWin):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(603, 553)
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.gridlayout = QtGui.QGridLayout(self.centralWidget)
		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget, rw=renWin)
		self.gridlayout.addWidget(self.vtkWidget, 0, 0, 1, 1)
		MainWindow.setCentralWidget(self.centralWidget)

class PCoordsView(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		   
		data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100521.mat'
		self.ds = DataSource(data_file)
		
		self.pc_class = PCoordsChart(self.ds)
		self.view = self.pc_class.GetView()
		
		self.ui.setupUi(self, self.view.GetRenderWindow())
		style = vtk.vtkInteractorStyleRubberBand2D()
		self.view.SetInteractorStyle(style)
		# Need this additional command (plus above set style) for this class
		self.view.SetInteractionModeTo2D()
		
		self.view.ResetCamera()
		self.view.GetInteractor().Start()
		
		
if __name__ == "__main__":

	app = QApplication(sys.argv)
	
	window = PCoordsView()
	
	window.show()
	sys.exit(app.exec_())
