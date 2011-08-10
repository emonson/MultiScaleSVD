#!/usr/local/bin/python
# -*- coding: utf-8 -*-

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import sys
import numpy as N
import vtk.util.numpy_support as VN
import scipy.io
import os

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
	def setupUi(self, MainWindow, renWin):
		MainWindow.setObjectName(_fromUtf8("MainWindow"))
		MainWindow.resize(400, 450)
		self.centralwidget = QtGui.QWidget(MainWindow)
		self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
		self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
		self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
		self.vtkWidget = QVTKRenderWindowInteractor(self.centralwidget, rw=renWin)
		self.vtkWidget.setObjectName(_fromUtf8("vtkWidget"))
		self.verticalLayout.addWidget(self.vtkWidget)
		self.horizontalSlider = QtGui.QSlider(self.centralwidget)
		self.horizontalSlider.setMinimum(0)
		self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
		self.horizontalSlider.setObjectName(_fromUtf8("horizontalSlider"))
		self.verticalLayout.addWidget(self.horizontalSlider)
		self.horizontalLayout = QtGui.QHBoxLayout()
		self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
		self.label = QtGui.QLabel(self.centralwidget)
		self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.label.setObjectName(_fromUtf8("label"))
		self.horizontalLayout.addWidget(self.label)
		self.comboBox = QtGui.QComboBox(self.centralwidget)
		self.comboBox.setObjectName(_fromUtf8("comboBox"))
		self.comboBox.addItem(_fromUtf8(""))
		self.comboBox.addItem(_fromUtf8(""))
		self.comboBox.addItem(_fromUtf8(""))
		self.horizontalLayout.addWidget(self.comboBox)
		self.verticalLayout.addLayout(self.horizontalLayout)

		MainWindow.setCentralWidget(self.centralwidget)
		self.menubar = QtGui.QMenuBar(MainWindow)
		self.menubar.setGeometry(QtCore.QRect(0, 0, 1100, 22))
		self.menubar.setObjectName(_fromUtf8("menubar"))
		self.menuFile = QtGui.QMenu(self.menubar)
		self.menuFile.setObjectName(_fromUtf8("menuFile"))
		MainWindow.setMenuBar(self.menubar)
		self.statusbar = QtGui.QStatusBar(MainWindow)
		self.statusbar.setObjectName(_fromUtf8("statusbar"))
		MainWindow.setStatusBar(self.statusbar)
		self.actionOpen = QtGui.QAction(MainWindow)
		self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
		self.actionExit = QtGui.QAction(MainWindow)
		self.actionExit.setObjectName(_fromUtf8("actionExit"))
		self.menuFile.addAction(self.actionOpen)
		self.menuFile.addAction(self.actionExit)
		self.menubar.addAction(self.menuFile.menuAction())

		MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Diffusion Embedding", None, QtGui.QApplication.UnicodeUTF8))
		self.label.setText(QtGui.QApplication.translate("MainWindow", "Transition style:", None, QtGui.QApplication.UnicodeUTF8))
		self.comboBox.setItemText(0, QtGui.QApplication.translate("MainWindow", "Alternating Axes", None, QtGui.QApplication.UnicodeUTF8))
		self.comboBox.setItemText(1, QtGui.QApplication.translate("MainWindow", "Diagonal", None, QtGui.QApplication.UnicodeUTF8))
		self.comboBox.setItemText(2, QtGui.QApplication.translate("MainWindow", "Simple", None, QtGui.QApplication.UnicodeUTF8))
		self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
		self.actionOpen.setText(QtGui.QApplication.translate("MainWindow", "Open", None, QtGui.QApplication.UnicodeUTF8))
		self.actionOpen.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))
		self.actionExit.setText(QtGui.QApplication.translate("MainWindow", "Exit", None, QtGui.QApplication.UnicodeUTF8))
		self.actionExit.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Q", None, QtGui.QApplication.UnicodeUTF8))

class SimpleView(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		# Number of slider divisions per integer value
		self.divs = 3

		# self.rot_method = 'alt_axis'
		self.rot_method = '111'
		# self.rot_method = 'simple'
		
		self.axis_rescale = False

		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		   
		# Set up a 2D scene, add an XY chart to it
		self.view = vtk.vtkContextView()
		self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
		
		self.ui.setupUi(self, self.view.GetRenderWindow())

		if self.rot_method == 'alt_axis':
			self.ui.comboBox.setCurrentIndex(0)
		elif self.rot_method == '111':
			self.ui.comboBox.setCurrentIndex(1)
		else:
			self.ui.comboBox.setCurrentIndex(2)
			
		style = vtk.vtkInteractorStyleRubberBand2D()
		self.view.GetInteractor().SetInteractorStyle(style)
		self.view.GetScene().SetInteractorStyle(style)
		# Need this additional command (plus above set style) for this class
		# self.view.SetInteractionModeTo2D()
		
		self.chart = vtk.vtkChartXY()
		self.chart.SetShowLegend(False)

		# Create a annotation link to access selection in parallel coordinates view
		annotationLink = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		annotationLink.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		annotationLink.GetCurrentSelection().GetNode(0).SetContentType(4)   # Indices
		# Connect the annotation link to the parallel coordinates representation
		self.chart.SetAnnotationLink(annotationLink)
		
		self.view.GetScene().AddItem(self.chart)		
		
		# Set up callback to update 3d render window when selections are changed in 
		#	parallel coordinates view
		annotationLink.AddObserver("AnnotationChangedEvent", self.selectionCallback)
						
		QtCore.QObject.connect(self.ui.actionExit, QtCore.SIGNAL("triggered()"), self.fileExit)
		QtCore.QObject.connect(self.ui.actionOpen, QtCore.SIGNAL("triggered()"), self.LoadData)
		QtCore.QObject.connect(self.ui.horizontalSlider, QtCore.SIGNAL("valueChanged(int)"), self.columnUpdate)
		QtCore.QObject.connect(self.ui.comboBox, QtCore.SIGNAL("currentIndexChanged(int)"), self.rotMethodChanged)
		
		self.LoadData()
		
		self.view.ResetCamera()
		self.ui.vtkWidget.GetRenderWindow().Render()
		self.view.GetInteractor().Start()
		
	def selectionCallback(self, caller, event):
		annSel = caller.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			if idxArr.GetNumberOfTuples() > 0:
				print VN.vtk_to_numpy(idxArr)

	# --------------------------------------
	def rotMethodChanged(self):
		
		cb_text = str(self.ui.comboBox.currentText())
		
		if 'Alternating' in cb_text:
			self.rot_method = 'alt_axis'
		elif 'Diagonal' in cb_text:
			self.rot_method = '111'
		else:
			self.rot_method = 'simple'
		
		self.columnUpdate()

	# --------------------------------------
# Simple non-orthogonal transition
	def F0_simple(self, t):
		return N.mat([[N.cos(t)],[N.sin(t)],[0.0]])
		
	def F1_simple(self, t):
		return N.mat([[0.0],[N.cos(t)],[N.sin(t)]])
			
# Rotate around (1/sqrt(3))(1,1,1)
	def F0_111(self, t):
		c = N.cos(t)
		s = N.sin(t)
		sqrt3 = 1.0/N.sqrt(3.0)
		inv3 = 1.0/3.0
		return N.mat([[c+inv3*(1.0-c)],[inv3*(1.0-c)+sqrt3*s],[inv3*(1.0-c)-sqrt3*s]])
		
	def F1_111(self, t):
		c = N.cos(t)
		s = N.sin(t)
		sqrt3 = 1.0/N.sqrt(3.0)
		inv3 = 1.0/3.0
		return N.mat([[inv3*(1.0-c)-sqrt3*s],[c+inv3*(1.0-c)],[inv3*(1.0-c)+sqrt3*s]])

# Rotate around (1,0,0)
	def F0_100(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[1.0],[0.0],[0.0]])
		
	def F1_100(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[0.0],[c],[s]])
			
# Rotate around (-1,0,0)
	def F0_n100(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[1.0],[0.0],[0.0]])
		
	def F1_n100(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[0.0],[c],[-s]])
			
# Rotate around (0,1,0)
	def F0_010(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[c],[0.0],[-s]])
		
	def F1_010(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[0.0],[1.0],[0.0]])
			
# Rotate around (0,-1,0)
	def F0_n010(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[c],[0.0],[s]])
		
	def F1_n010(self, t):
		c = N.cos(t)
		s = N.sin(t)
		return N.mat([[0.0],[1.0],[0.0]])
			
	def columnUpdate(self):
		slider_val = self.ui.horizontalSlider.value()
		x_axis = self.chart.GetAxis(1)
		y_axis = self.chart.GetAxis(0)
		
		val = float(slider_val)/float(self.divs)
		(frac_val, int_val) = N.modf(val)
		
		if (self.rot_method == '111'):
			t = (2.0*N.pi/3.0)*frac_val
		else:
			t = (N.pi/2.0)*frac_val
				
		# Special case for very last value
		if (int_val == self.data.shape[1]-2):
			this_data = N.concatenate((self.data[:, int_val:(int_val+3)],N.zeros_like(self.data[:,0])),axis=1)
		else:
			this_data = self.data[:, int_val:(int_val+3)]
		
		if (self.rot_method == 'alt_axis'):
			if N.mod(int_val,8) == 0:
				xt = self.F0_n010(t)
				yt = self.F1_n010(t)
			if N.mod(int_val,8) == 1:
				xt = self.F1_010(t)
				yt = self.F0_010(t)				
			if N.mod(int_val,8) == 2:
				xt = self.F0_010(t)
				yt = -self.F1_010(t)				
			if N.mod(int_val,8) == 3:
				xt = -self.F1_n010(t)
				yt = -self.F0_n010(t)				
			if N.mod(int_val,8) == 4:
				xt = -self.F0_n010(t)
				yt = -self.F1_n010(t)
			if N.mod(int_val,8) == 5:
				xt = -self.F1_010(t)
				yt = -self.F0_010(t)				
			if N.mod(int_val,8) == 6:
				xt = -self.F0_010(t)
				yt = self.F1_010(t)				
			if N.mod(int_val,8) == 7:
				xt = self.F1_n010(t)
				yt = self.F0_n010(t)
				
		elif (self.rot_method == '111'):
			xt = self.F0_111(t)
			yt = self.F1_111(t)
			
		else:
			xt = self.F0_simple(t)
			yt = self.F1_simple(t)		
			
		xo = int_val
		yo = int_val
		self.r0[:] = this_data*xt
		self.r1[:] = this_data*yt
		x_axis.SetTitle('%d  ~  %.2f   %.2f   %.2f' % (int_val, xt[0], xt[1], xt[2]))
		y_axis.SetTitle('%d  ~  %.2f   %.2f   %.2f' % (int_val, yt[0], yt[1], yt[2]))
		
		self.table.Modified()
		# Recalculate axis bounds
		if self.axis_rescale: 
			self.chart.RecalculateBounds()
		self.view.Render()

# -----------------
	def LoadData(self):
		
		# data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/yaleB_diffMap.mat'
		if os.path.exists('/Users/emonson/Data/Fodava/EMoGWDataSets'):
			self.openFilesDefaultPath = '/Users/emonson/Data/Fodava/EMoGWDataSets'
		else:
			self.openFilesDefaultPath = QtCore.QDir.homePath()
		
		picked_file = QtGui.QFileDialog.getOpenFileName(self,
				"Load Saved Matlab File",
				self.openFilesDefaultPath,
				"All Files (*);;Matlab Files (*.mat)")
		data_file = str(picked_file)
		
		# DataSource loads .mat file and can generate data from it for other views
		print "Loading data from ", str(data_file)
		
		self.last_data_dir = os.path.dirname(str(data_file))

		# Hack for labels
		if 'yaleB' in data_file:
			l_idx = 2
		else:
			l_idx = 0
			
		print 'Trying to really load now...'
		try:
			MatInput = scipy.io.loadmat(data_file, struct_as_record=True, chars_as_strings=True)
		except:
			print 'loadmat crapping out for some reason...'
			raise IOError, "Can't load supplied matlab file"
			# return

		# Create original multi-dimensional data matrix
		e_vecs = MatInput['EigenVecs']
		e_vals = MatInput['EigenVals']
		self.data = N.mat(e_vecs*e_vals.T[0])
		d = self.data.shape[1]
		
		# Try out fixed axis bounds for less distracting axis label switching
		# Still a problem that some points in in-between projections go out of these bounds...
		d_max = N.abs(self.data.max())
		d_min = N.abs(self.data.min())
		if d_max > d_min: 
			w_max = d_max
		else: 
			w_max = d_min
		
		# Rescale data so that max is 1.0
		self.data = self.data/w_max
		
		self.ui.horizontalSlider.setMaximum((d-2)*self.divs+1)
		self.ui.horizontalSlider.setPageStep(self.divs)
		self.ui.horizontalSlider.setValue(0)

		# Start by plotting the first two columns of the data
		# These column matrices will hold results of all later calculations
		# and share memory with the table columns plotted
		self.r0 = self.data[:,0].copy()
		self.r1 = self.data[:,1].copy()
		
		# Create a table with some points in it...
		self.table = vtk.vtkTable()
		self.rv0 = VN.numpy_to_vtk(self.r0)
		self.rv0.SetName('column0')
		self.rv1 = VN.numpy_to_vtk(self.r1)
		self.rv1.SetName('column1')
		self.table.AddColumn(self.rv0)
		self.table.AddColumn(self.rv1)
		
		# Color-by array
		labels_array = MatInput['Labels']
		self.cat_labels = N.zeros_like(labels_array)
		for ii in range(labels_array.shape[0]):
			cl_unique = set(labels_array[ii,:])
			cl_map = {}
			for jj,vv in enumerate(cl_unique):
				cl_map[vv] = jj
			self.cat_labels[ii,:] = N.array([cl_map[vv] for vv in labels_array[ii,:]])
		
		# Specifying a single label array here
		rvC = VN.numpy_to_vtk(self.cat_labels[l_idx,:].copy(),deep=True)
		rvC.SetName('color')		
		self.table.AddColumn(rvC)
		
		self.chart.ClearPlots()
		points1 = self.chart.AddPlot(vtk.vtkChart.POINTS)
		points1.SetInput(self.table, 0, 1)
		points1.SetColor(0, 0, 0, 255)
		points1.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)
		points1.SetScalarVisibility(1)
		points1.SetLookupTable(self.GetCategoryLUT(l_idx))
		points1.SelectColorArray('color')
		
		sc = 1.01
		x_axis = self.chart.GetAxis(1)
		y_axis = self.chart.GetAxis(0)

		if not self.axis_rescale: 
			y_axis.SetBehavior(1)
			y_axis.SetNotation(2)
		y_axis.SetRange(-sc, sc)
		y_axis.RecalculateTickSpacing()
		y_axis.SetPrecision(1)
		if not self.axis_rescale: 
			x_axis.SetBehavior(1)
			x_axis.SetNotation(2)
		x_axis.SetRange(-sc, sc)
		x_axis.RecalculateTickSpacing()
		x_axis.SetPrecision(1)

		x_axis.GetTitleProperties().BoldOff()
		y_axis.GetTitleProperties().BoldOff()
		x_axis.GetTitleProperties().SetFontSize(10)
		y_axis.GetTitleProperties().SetFontSize(10)
		
		self.columnUpdate()


# -----------------
	def GetCategoryLUT(self, idx = 0):
		"""Returns a LUT for category coloring. Result depends on number
		of categories. Pass an index to get a proper label map (otherwise
		it defaults to the first (zeroth) label."""

		cl = []
		lut = vtk.vtkLookupTable()
		
		num_labels = len(N.unique(self.cat_labels[idx,:]))
		
		if num_labels > 8 and num_labels <= 13:
			cl.append([float(cc)/255.0 for cc in [228, 26, 28]])  # Colorbrewer Set2 modY+4
			cl.append([float(cc)/255.0 for cc in [55, 126, 184]])
			cl.append([float(cc)/255.0 for cc in [77, 175, 74]])
			cl.append([float(cc)/255.0 for cc in [152, 78, 163]])
			cl.append([float(cc)/255.0 for cc in [255, 127, 0]])
			cl.append([float(cc)/255.0 for cc in [245, 193, 61]])
			cl.append([float(cc)/255.0 for cc in [166, 86, 40]])
			cl.append([float(cc)/255.0 for cc in [247, 129, 191]])
			cl.append([float(cc)/255.0 for cc in [153, 153, 153]])
			cl.append([float(cc)/255.0 for cc in [143, 0, 0]])
			cl.append([float(cc)/255.0 for cc in [22, 65, 110]])
			cl.append([float(cc)/255.0 for cc in [40, 115, 33]])
			cl.append([float(cc)/255.0 for cc in [61, 61, 61]])

			lutNum = len(cl)
			lut.SetNumberOfTableValues(lutNum)
			lut.Build()
			for ii,cc in enumerate(cl):
				lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
			lut.SetRange(0,len(cl)-1)
			lut.Build()
			return lut
							
		if num_labels <= 8:
			cl.append([float(cc)/255.0 for cc in [27, 158, 119]])	# Colorbrewer Dark2
			cl.append([float(cc)/255.0 for cc in [217, 95, 2]])
			cl.append([float(cc)/255.0 for cc in [117, 112, 179]])
			cl.append([float(cc)/255.0 for cc in [231, 41, 138]])
			cl.append([float(cc)/255.0 for cc in [102, 166, 30]])
			cl.append([float(cc)/255.0 for cc in [230, 171, 2]])
			cl.append([float(cc)/255.0 for cc in [166, 118, 29]])
			cl.append([float(cc)/255.0 for cc in [102, 102, 102]])
	
	# 		cl.append([float(cc)/255.0 for cc in [102, 102, 102]])	# Colorbrewer Dark2 (rev)
	# 		cl.append([float(cc)/255.0 for cc in [166, 118, 29]])
	# 		cl.append([float(cc)/255.0 for cc in [230, 171, 2]])
	# 		cl.append([float(cc)/255.0 for cc in [102, 166, 30]])
	# 		cl.append([float(cc)/255.0 for cc in [231, 41, 138]])
	# 		cl.append([float(cc)/255.0 for cc in [117, 112, 179]])
	# 		cl.append([float(cc)/255.0 for cc in [217, 95, 2]])
	# 		cl.append([float(cc)/255.0 for cc in [27, 158, 119]])
	
			lutNum = len(cl)
			lut.SetNumberOfTableValues(lutNum)
			lut.Build()
			for ii,cc in enumerate(cl):
				lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
			lut.SetRange(0,len(cl)-1)
			lut.Build()
			return lut

		if num_labels > 13:
			
			lut.SetNumberOfTableValues(num_labels)
			lut.SetValueRange(0.8, 0.8)
			lut.SetSaturationRange(0.8, 0.8)
			lut.SetHueRange(0.0,1.0)
			lut.SetRampToLinear()
			lut.SetRange(0,num_labels-1)
			lut.Build()
			return lut
	
		lut.SetNumberOfTableValues(256)
		lut.SetValueRange(1.0, 1.0)
		lut.SetSaturationRange(1.0, 1.0)
		lut.SetHueRange(0.0,1.0)
		lut.SetRampToLinear()
		lut.SetRange(0,10)
		lut.Build()
		return lut

# -------------------------------------------
	def fileOpen(self):

		openFilesDefaultPath = ''

		file = QtGui.QFileDialog.getOpenFileName(self,
				"Load Saved Matlab File",
				self.last_data_dir,
				"All Files (*);;Matlab Files (*.mat)")

		if file:
			pass

	def fileExit(self):

		# Usually would use the qApp global variable qApp.quit(), but wasn't working...
		QtGui.QApplication.instance().quit()

if __name__ == "__main__":

	app = QApplication(sys.argv)
	window = SimpleView()
	window.show()
	sys.exit(app.exec_())
