"""
simple_graph_view.py

Example using PyQt4 & VTK Infovis

10 Sept 2009 -- E Monson

"""

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication
import vtk
import sys
from vtk.util import numpy_support as VN

from QtSimpleView import Ui_MainWindow

class SimpleView(QtGui.QMainWindow):
	
	def __init__(self, parent = None):
	
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		
		self.ListView = vtk.vtkQtListView()
		
		# Create a annotation link to access selection in parallel coordinates view
		annotationLink = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		annotationLink.GetCurrentSelection().GetNode(0).SetFieldType(vtk.vtkSelectionNode.ROW)		 # Row
		annotationLink.GetCurrentSelection().GetNode(0).SetContentType(vtk.vtkSelectionNode.INDICES)	 # Indices
		# Connect the annotation link to the parallel coordinates representation

		# Set up callback to update 3d render window when selections are changed in 
		# parallel coordinates view
		annotationLink.AddObserver("AnnotationChangedEvent", self.selectionCallback)
								
		titles = vtk.vtkStringArray()
		titles.SetName('tick_labels')
		titles.SetNumberOfComponents(1)
		titles.InsertNextValue("War and Peace")
		titles.InsertNextValue("Barnaby Jones")
		titles.InsertNextValue("One Flew Over the Cookoo's Nest")
		titles.InsertNextValue("War and Peace")
		titles.InsertNextValue("Barnaby Jones")
		titles.InsertNextValue("One Flew Over the Cookoo's Nest")
		
		colors = vtk.vtkUnsignedIntArray()
		colors.SetName('color')
		colors.SetNumberOfComponents(1)
		colors.InsertNextValue(1)
		colors.InsertNextValue(2)
		colors.InsertNextValue(3)
		colors.InsertNextValue(4)
		colors.InsertNextValue(4)
		colors.InsertNextValue(6)
		
		# Create a table with some points in it...
		table = vtk.vtkTable()
		table.AddColumn(titles)
		table.AddColumn(colors)
		
		vt = vtk.vtkViewTheme()
		lut = vtk.vtkLookupTable()
		cl = []
		cl.append([float(cc)/255.0 for cc in [27, 158, 119]]) # Colorbrewer Dark2
		cl.append([float(cc)/255.0 for cc in [217, 95, 2]])
		cl.append([float(cc)/255.0 for cc in [117, 112, 179]])
		cl.append([float(cc)/255.0 for cc in [231, 41, 138]])
		cl.append([float(cc)/255.0 for cc in [102, 166, 30]])
		cl.append([float(cc)/255.0 for cc in [230, 171, 2]])
		cl.append([float(cc)/255.0 for cc in [166, 118, 29]])
		cl.append([float(cc)/255.0 for cc in [102, 102, 102]])
		
		lutNum = len(cl)
		lut.SetNumberOfTableValues(lutNum)
		lut.Build()
		for ii,cc in enumerate(cl):
			lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
		lut.SetRange(0,len(cl)-1)
		vt.SetPointLookupTable(lut)

		self.ListView.AddRepresentationFromInput(table)
		self.ListView.SetVisibleColumn(0)
		self.ListView.SetFieldType(vtk.vtkQtListView.ROW_DATA)
		
		self.ListView.SetColorByArray(True)
		self.ListView.SetColorArrayName('color')
		self.ListView.SetDecorationStrategy(0)
		# self.ListView.SetAlternatingRowColors(True)
		self.ListView.ApplyViewTheme(vt)
		
		self.ListView.GetRepresentation().SetAnnotationLink(annotationLink)
		self.ListView.Update()

		# Set widget for the list view
		self.ui.centralWidget.layout().addWidget(self.ListView.GetWidget())
		
	def selectionCallback(self,caller, event):
		annSel = caller.GetCurrentSelection()
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			if idxArr.GetNumberOfTuples() > 0:
				print VN.vtk_to_numpy(idxArr)
		
if __name__ == "__main__":

	app = QApplication(sys.argv)
	window = SimpleView()
	window.show()
	sys.exit(app.exec_())
