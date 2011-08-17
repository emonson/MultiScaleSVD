"""
simpleview.py

Example using PyQt4, Qt Designer & VTK Infovis
Plus M Maggioni multiscale graph representation

09 Jan 2009 -- E Monson

"""

from PyQt4 import QtCore, QtGui
from ui_simpleview import Ui_MainWindow
import vtk
import vtk.util.numpy_support as VN
import numpy as N
import scipy.io
import os

class SimpleView(QtGui.QMainWindow):
    
    def __init__(self, parent = None):
    
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.view = vtk.vtkGraphLayoutView()
        self.view.SetInteractor(self.ui.vtkWidget.GetRenderWindow().GetInteractor())
        self.ui.vtkWidget.SetRenderWindow(self.view.GetRenderWindow())

        self.level = 2
        self.basisNum = 0
        self.cutoff = 1.7e-5
        self.basisCutoff = 1.0e-6
        self.ebSumCutoff = 0.1
        self.loadGraph()
        
        self.cutoffValid = QtGui.QDoubleValidator(0.0, 1.0, 12, self.ui.lineEdit_cutoff)
        self.ui.lineEdit_cutoff.setValidator(self.cutoffValid)
        self.basisCutoffValid = QtGui.QDoubleValidator(0.0, 1.0, 12, self.ui.lineEdit_basisCutoff)
        self.ui.lineEdit_basisCutoff.setValidator(self.basisCutoffValid)

        # Connect signals and slots
        # QtCore.QObject.connect(self.ui.actionNew, QtCore.SIGNAL("triggered()"), self.fileOpen)
        QtCore.QObject.connect(self.ui.actionExit, QtCore.SIGNAL("triggered()"), self.fileExit)
        QtCore.QObject.connect(self.ui.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), self.basisUpdate)
        # QtCore.QObject.connect(self.ui.spinBox_basisIndex, QtCore.SIGNAL("valueChanged(int)"), self.basisUpdate)
        QtCore.QObject.connect(self.ui.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), self.levelUpdate)
        # QtCore.QObject.connect(self.ui.spinBox_level, QtCore.SIGNAL("valueChanged(int)"), self.levelUpdate)
        QtCore.QObject.connect(self.ui.listWidget, QtCore.SIGNAL("itemClicked(QListWidgetItem *)"), self.toggleVisibility)
        QtCore.QObject.connect(self.ui.lineEdit_cutoff, QtCore.SIGNAL("editingFinished()"), self.changeCutoff)
        QtCore.QObject.connect(self.ui.lineEdit_basisCutoff, QtCore.SIGNAL("editingFinished()"), self.changeBasisCutoff)

        self.ui.hSlider_basisIndex.setValue(self.basisNum)
        self.ui.spinBox_basisIndex.setValue(self.basisNum)
        self.ui.hSlider_level.setValue(self.level)
        self.ui.spinBox_level.setValue(self.level)
        
        # TODO: Problem with size of ExtBasis for levels below 2...
        self.ui.hSlider_level.setMinimum(2)
        self.ui.spinBox_level.setMinimum(2)
        
    def loadGraph(self):
        
        subgraph = vtk.vtkRandomGraphSource()

        self.view.SetRepresentationFromInputConnection(subgraph.GetOutputPort())
        self.view.ResetCamera()
        self.view.Render()
        self.view.GetInteractor().Start()


    def changeCutoff(self):
        
        tmpCutoffTxt = self.ui.lineEdit_cutoff.text()
        tmpCutoff = tmpCutoffTxt.toDouble()
        if tmpCutoff[1]:
            self.cutoff = tmpCutoff[0]
            # Update multi-resolution graph
            msGraph = self.buildSubGraph(self.level) 
            self.msMapper.SetInput(msGraph)
            self.msMapper.Modified()
            # self.view.Update()
            self.view.Render()
            
    def changeBasisCutoff(self):
        
        tmpCutoffTxt = self.ui.lineEdit_basisCutoff.text()
        tmpCutoff = tmpCutoffTxt.toDouble()
        if tmpCutoff[1]:
            self.basisCutoff = tmpCutoff[0]
            # Update multi-resolution graph
            self.selection.RemoveAllThresholds()
            self.selection.AddThreshold(self.basisCutoff, 10)
            self.selection.Modified()
            # self.view.Update()
            self.view.Render()
            
    
    def buildSubGraph(self, level):
        
        table = vtk.vtkTable()
        col0 = vtk.vtkIntArray()
        col0.SetName('index1')
        col1 = vtk.vtkIntArray()
        col1.SetName('index2')
        val = vtk.vtkDoubleArray()
        val.SetName('weight')
        
        Tmat = self.G.Tree[level,0].T
        minmax = "(%3.2e, %3.2e)" % (abs(Tmat.data).min(), abs(Tmat.data).max())
        self.ui.label_cutoff_minmax.setText(minmax)
        
        ebSum = self.ExtBasis.sum(0)
        idxs = []
        
        for ii in range(Tmat.nzmax):
            xx = Tmat.rowcol(ii)[0]
            yy = Tmat.rowcol(ii)[1]
            if ((ebSum[0,xx] > self.ebSumCutoff) and (ebSum[0,yy] > self.ebSumCutoff)):
                # print ebSum[0,xx], ebSum[0,yy]
                # Always include diagonal elements (no matter what cutoff) so all vertices are created
                if ((abs(Tmat.getdata(ii)) > self.cutoff) | (Tmat.rowcol(ii)[0]==Tmat.rowcol(ii)[1])):
                    col0.InsertNextValue(xx)
                    col1.InsertNextValue(yy)
                    # But, take out self-loops for edge color scale
                    if (xx == yy):
                        val.InsertNextValue(0.0)        
                        idxs.append(self.G.Tree[level,0].ExtIdxs[xx]-1)
                    else:
                        val.InsertNextValue(abs(Tmat.getdata(ii)))        
        
        table.AddColumn(col0)
        table.AddColumn(col1)
        table.AddColumn(val)
        
        tgraph = vtk.vtkTableToGraph()
        tgraph.SetInput(table)
        tgraph.AddLinkVertex('index2', 'stuff', False)
        tgraph.AddLinkVertex('index1', 'stuff', False)
        tgraph.AddLinkEdge('index2', 'index1')
        
        rawGraph = tgraph.GetOutput()
        rawGraph.Update()
        # print graph
        
        # Grab points from original set with ExtIdxs
        ptsArray = VN.vtk_to_numpy(self.pts.GetData())
        # Remember Matlab uses one-based indices (ExtIdxs-1 for Numpy)
        # print (self.G.Tree[level,0].ExtIdxs-1)
        # print ptsArray[self.G.Tree[level,0].ExtIdxs-1,:]
        ptsMatrix = N.matrix(ptsArray)
        ptsSubMatrix = ptsMatrix[idxs,:]
        ptsSub = VN.numpy_to_vtk(ptsSubMatrix)
        # print VN.vtk_to_numpy(ptsSub)
        pts = vtk.vtkPoints()
        pts.SetData(ptsSub)
        rawGraph.SetPoints(pts)
        # print pts
            
        strategy = vtk.vtkPassThroughLayoutStrategy()
        layout = vtk.vtkGraphLayout()
        layout.SetInput(rawGraph)
        layout.SetLayoutStrategy(strategy)
        
        edgeLayout = vtk.vtkEdgeLayout()
        edgeStrategy = vtk.vtkArcParallelEdgeStrategy()
        edgeStrategy.SetNumberOfSubdivisions(50)
        edgeLayout.SetInputConnection(layout.GetOutputPort())
        edgeLayout.SetLayoutStrategy(edgeStrategy)
        edgeLayout.Update()
        
        return edgeLayout.GetOutput()
         
    def basisUpdate(self, index):
        
        QtCore.QObject.disconnect(self.ui.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), self.basisUpdate)
        QtCore.QObject.disconnect(self.ui.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), self.levelUpdate)
        
        self.basisNum = index
        self.ui.spinBox_basisIndex.setValue(index)
        self.ui.hSlider_basisIndex.setValue(index)
        self.selection.RemoveAllIDs()
        
        # Try to directly swap ExtBasis vertex data
        # print 'ExtBasis shape: ' + str(self.ExtBasis[:,self.basisNum].data.shape) + ' (level,basis) (' + str(self.level) + ',' + str(self.basisNum) + ')'
        
        # Indexing of ExtBasis is backwards from graph, so need to reverse with [::-1]
        # and array not "contiguous" if don't do .copy()
        self.Esub[:] = self.ExtBasis[:,self.basisNum].data[::-1].copy()
        self.EsubSq[:] = self.Esub**2
        # self.Esub[:] = N.random.random(self.ExtBasis[:,self.basisNum].data.shape[0])
        # SubIdxs = (self.Esub > 0.001).nonzero()
        # SubIdxs = (self.Esub**2 > 0.8).nonzero()
        self.vertexData.Modified()
        self.degree.Update()
                
        # for ii in SubIdxs[0]:
        #     self.selection.AddID(0,ii)
        vertRange = self.vertMapper.GetInput().GetPointData().GetArray(self.vertMapper.GetArrayName()).GetRange()
        self.vertMapper.SetScalarRange(-1.0*vertRange[1], vertRange[1])
        self.selection.Modified()
        # tmpMin = self.EsubSq.min()
        # tmpMax = self.EsubSq.max()
        # self.map0.SetScalarRange(-100.0*tmpMax, 100.0*tmpMax)
        minmax = "(%3.2e, %3.2e)" % (self.EsubSq.min(), self.EsubSq.max())
        self.ui.label_basisCutoff_minmax.setText(minmax)
        # self.view.Update()
        self.view.Render()
        
        QtCore.QObject.connect(self.ui.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), self.basisUpdate)
        QtCore.QObject.connect(self.ui.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), self.levelUpdate)
        
    def levelUpdate(self, index):
        
        QtCore.QObject.disconnect(self.ui.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), self.basisUpdate)
        QtCore.QObject.disconnect(self.ui.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), self.levelUpdate)
        
        self.level = index
        self.ui.spinBox_level.setValue(index)
        self.ui.hSlider_level.setValue(index)
        # Start setting up basis function subgraph
        self.ExtBasis = self.G.Tree[self.level,0].ExtBasis
        self.basisMax = self.ExtBasis.shape[1]-1    # zero-based indices
        if self.basisNum > self.basisMax:
            self.basisNum = self.basisMax
            self.ui.hSlider_basisIndex.setValue(self.basisMax)
            self.ui.spinBox_basisIndex.setValue(self.basisMax)
        self.ui.hSlider_basisIndex.setMinimum(0)
        self.ui.hSlider_basisIndex.setMaximum(self.basisMax)
        self.ui.spinBox_basisIndex.setMinimum(0)
        self.ui.spinBox_basisIndex.setMaximum(self.basisMax)

        # self.selection.RemoveAllIDs()
        # print 'ExtBasis shape: ' + str(self.ExtBasis[:,self.basisNum].data.shape) + ' (level,basis) (' + str(self.level) + ',' + str(self.basisNum) + ')'
        
        # Indexing of ExtBasis is backwards from graph, so need to reverse with [::-1]
        # and array not "contiguous" if don't do .copy()
#         self.Esub[:] = self.ExtBasis[:,self.basisNum].data[::-1].copy()
#         self.EsubSq[:] = self.Esub**2
#         # self.Esub[:] = N.random.random(self.ExtBasis[:,self.basisNum].data.shape[0])
#         # SubIdxs = (self.Esub > 0.001).nonzero()
#         # SubIdxs = (self.Esub**2 > 0.8).nonzero()
#         self.vertexData.Modified()
#         self.degree.Update()
        
        # Update multi-resolution graph
        msGraph = self.buildSubGraph(index) 
        self.msMapper.SetInput(msGraph)
        self.msMapper.Modified()

        # for ii in SubIdxs[0]:
        #     self.selection.AddID(0,ii)
        self.selection.Modified()
        # tmpMin = self.Esub.min()
        # tmpMax = self.Esub.max()
        # self.map0.SetScalarRange(-100.0*tmpMax, 100.0*tmpMax)
        self.basisUpdate(self.ui.hSlider_basisIndex.value())
        
        QtCore.QObject.connect(self.ui.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), self.basisUpdate)
        QtCore.QObject.connect(self.ui.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), self.levelUpdate)
        
    def toggleVisibility(self, item):
        
        tmpQtActor = item.data(QtCore.Qt.UserRole)
        tmpActor = tmpQtActor.toPyObject()
        if tmpActor.GetClassName() == 'vtkActorCollection':
            # Loop through actor collection
            tmpActor.InitTraversal()
            singleActor = tmpActor.GetNextActor()
            while singleActor is not None:
                self.toggleActorVisibility(item,singleActor)
                singleActor = tmpActor.GetNextActor()
        else:
            # Single actor
            self.toggleActorVisibility(item,tmpActor)
            
    def toggleActorVisibility(self, item, tmpActor):
        actorVis = tmpActor.GetVisibility()
        itemCheck = item.checkState()   # QtCore.Qt.Checked = 2, QtCore.Qt.Unchecked = 0
        if (actorVis == 0) and (itemCheck == QtCore.Qt.Checked):
            tmpActor.SetVisibility(1)
            # self.view.Update()
            self.view.Render()
        elif (actorVis == 1) and (itemCheck == QtCore.Qt.Unchecked):
            tmpActor.SetVisibility(0)
            # self.view.Update()
            self.view.Render()
        
        
    def fileExit(self):
        
        # Usually would use the qApp global variable qApp.quit(), but wasn't working...
        QtGui.QApplication.instance().quit()
        
