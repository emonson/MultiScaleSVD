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
    
    def __init__(parent = None):
    
        QtGui.QMainWindow.__init__(parent)
        ui = Ui_MainWindow()
        ui_window.setupUi()
        
        view = vtk.vtkGraphLayoutView()
        view.SetInteractor(ui_window.vtkWidget.GetRenderWindow().GetInteractor())
        ui_window.vtkWidget.SetRenderWindow(view.GetRenderWindow())

        level = 2
        basisNum = 0
        cutoff = 1.7e-5
        basisCutoff = 1.0e-6
        ebSumCutoff = 0.1
        loadGraph()
        
        cutoffValid = QtGui.QDoubleValidator(0.0, 1.0, 12, ui_window.lineEdit_cutoff)
        ui_window.lineEdit_cutoff.setValidator(cutoffValid)
        basisCutoffValid = QtGui.QDoubleValidator(0.0, 1.0, 12, ui_window.lineEdit_basisCutoff)
        ui_window.lineEdit_basisCutoff.setValidator(basisCutoffValid)

        # Connect signals and slots
        # QtCore.QObject.connect(ui_window.actionNew, QtCore.SIGNAL("triggered()"), fileOpen)
        QtCore.QObject.connect(ui_window.actionExit, QtCore.SIGNAL("triggered()"), fileExit)
        QtCore.QObject.connect(ui_window.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
        # QtCore.QObject.connect(ui_window.spinBox_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
        QtCore.QObject.connect(ui_window.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
        # QtCore.QObject.connect(ui_window.spinBox_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
        QtCore.QObject.connect(ui_window.listWidget, QtCore.SIGNAL("itemClicked(QListWidgetItem *)"), toggleVisibility)
        QtCore.QObject.connect(ui_window.lineEdit_cutoff, QtCore.SIGNAL("editingFinished()"), changeCutoff)
        QtCore.QObject.connect(ui_window.lineEdit_basisCutoff, QtCore.SIGNAL("editingFinished()"), changeBasisCutoff)

        ui_window.hSlider_basisIndex.setValue(basisNum)
        ui_window.spinBox_basisIndex.setValue(basisNum)
        ui_window.hSlider_level.setValue(level)
        ui_window.spinBox_level.setValue(level)
        
        # TODO: Problem with size of ExtBasis for levels below 2...
        ui_window.hSlider_level.setMinimum(2)
        ui_window.spinBox_level.setMinimum(2)
        
    def loadGraph():
        
        # ----------
        # Load and construct whole graph and multi-resolution data from Matlab structure
        dataDir = '/Users/emonson/Programming/Matlab/EMonson/Fodava/DocumentAnalysis/Analysis/'
        filename = dataDir + 'X20_042709b.mat'
        # filename = '/Users/emonson/Programming/Python/VTK/X20_040609b.mat'
        X = scipy.io.loadmat(filename)
        # Get graph structure G out of matlab variables
        G = X['G']
        
        # ----------
        # Set multi-resolution level bounds in GUI sliders
        levelMax = G.Tree.shape[0]-1
        ui_window.hSlider_level.setMinimum(1)
        ui_window.hSlider_level.setMaximum(levelMax)
        ui_window.spinBox_level.setMinimum(1)
        ui_window.spinBox_level.setMaximum(levelMax)
        
        # Start setting up basis function for display as subgraph
        ExtBasis = GTree[level,0]['ExtBasis'][0][0][0]
        basisMax = ExtBasis.shape[1]-1    # zero-based indices
        
        # Set particular level basis function bounds in GUI sliders
        ui_window.hSlider_basisIndex.setMinimum(0)
        ui_window.hSlider_basisIndex.setMaximum(basisMax)
        ui_window.spinBox_basisIndex.setMinimum(0)
        ui_window.spinBox_basisIndex.setMaximum(basisMax)
        
        # Build table which will become graph
        table = vtk.vtkTable()
        col0 = vtk.vtkIntArray()
        col0.SetName('index1')
        col1 = vtk.vtkIntArray()
        col1.SetName('index2')
        val = vtk.vtkDoubleArray()
        val.SetName('weight')
        
        Tmat = G.T
        # Tmat = G.W
        
        for ii in range(Tmat.nzmax):
            col0.InsertNextValue(Tmat.rowcol(ii)[0])
            col1.InsertNextValue(Tmat.rowcol(ii)[1])
            val.InsertNextValue(abs(Tmat.getdata(ii)))
        
        table.AddColumn(col0)
        table.AddColumn(col1)
        table.AddColumn(val)
        
        # Vertex links need to be done with index2 first or indexing won't be right...
        # TODO: Make this foolproof so that graph always ends up with correct ordering of indices...
        tgraph = vtk.vtkTableToGraph()
        tgraph.SetInput(table)
        tgraph.AddLinkVertex('index2', 'stuff', False)
        tgraph.AddLinkVertex('index1', 'stuff', False)
        tgraph.AddLinkEdge('index2', 'index1')
        
        rawGraph = tgraph.GetOutput()
        rawGraph.Update()
        # print graph
        
        # Load and assign whole graph pre-layout coordinates
        ptsFile = os.path.splitext(filename)[0] + '_pts.vtp'
        if os.path.exists(ptsFile):
            polyreader = vtk.vtkXMLPolyDataReader()
            polyreader.SetFileName(ptsFile)
            polyreader.Update()
            pts = polyreader.GetOutput().GetPoints()
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
        
        graph = edgeLayout.GetOutput()
        graph.Update()
        
        # --------
        # Add ExtBasis to graph data & Select particular basis function
        
        # print 'ExtBasis shape: ' + str(ExtBasis[:,basisNum].data.shape) + ' (level,basis) (' + str(level) + ',' + str(basisNum) + ')'
        
        # Indexing of ExtBasis is backwards from graph, so need to reverse with [::-1]
        # and array not "contiguous" if don't do .copy()
        Esub = ExtBasis[:,basisNum].data[::-1].copy()
        EsubSq = Esub**2
        # Esub = N.random.random(ExtBasis[:,basisNum].data.shape[0])
        # SubIdxs = (Esub > 0.001).nonzero()
        # SubIdxs = (Esub**2 > 0.8).nonzero()
        
        # Set ExtBasis vertex data from numpy array
        basisFunc = VN.numpy_to_vtk(Esub)
        basisFunc.SetName('ExtBasis')
        basisFuncSq = VN.numpy_to_vtk(EsubSq)
        basisFuncSq.SetName('ExtBasisSq')
        
        vertexData = graph.GetVertexData()
        vertexData.AddArray(basisFunc)
        vertexData.AddArray(basisFuncSq)
        
        selection = vtk.vtkSelectionSource()
        selection.SetContentType(7) # vtkSelection::THRESHOLDS
        # selection.SetContentType(2) # vtkSelection::PEDIGREE_IDS
        selection.SetFieldType(3) # vtkSelection::VERTEX
        selection.SetArrayName("ExtBasisSq")
        selection.AddThreshold(basisCutoff, 10)
        # TODO: There was something wrong with the indexing in the PEDIGREE_IDS selection...
        # for ii in SubIdxs[0]:
        #     selection.AddID(0,ii)
        minmax = "(%3.2e, %3.2e)" % (EsubSq.min(), EsubSq.max())
        ui_window.label_basisCutoff_minmax.setText(minmax)
        selection.Update()
        
        # ----------
        # Back to pipeline
        degree = vtk.vtkVertexDegree()
        degree.SetInput(graph)
        
        subgraph = vtk.vtkExtractSelectedGraph()
        subgraph.SetRemoveIsolatedVertices(False)
        subgraph.SetInputConnection(degree.GetOutputPort())
        subgraph.SetSelectionConnection(selection.GetOutputPort())
        
        # +++++++++++++
        graphToPoly = vtk.vtkGraphToPolyData()
        graphToPoly.SetInputConnection(subgraph.GetOutputPort())
        edgeMapper = vtk.vtkPolyDataMapper()
        edgeMapper.SetInputConnection(graphToPoly.GetOutputPort())
        edgeMapper.SetScalarModeToUseCellData()
        edgeMapper.SetScalarVisibility(False)
        edgeMapper.SetImmediateModeRendering(True)
        edgeActor = vtk.vtkActor()
        edgeActor.SetMapper(edgeMapper)
        edgeActor.SetPosition(0, 0, -0.003);
        
        lut = vtk.vtkLookupTable()
        lutNum = 256
        lut.SetNumberOfTableValues(lutNum)
        ctf = vtk.vtkColorTransferFunction()
        ctf.SetColorSpaceToDiverging()
        ctf.AddRGBPoint(0.0, 0, 0, 1.0)
        ctf.AddRGBPoint(1.0, 1.0, 0, 0)
        for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
            cc = ctf.GetColor(ss)
            lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
        
        vertGlyph = vtk.vtkVertexGlyphFilter()
        vertGlyph.SetInputConnection(subgraph.GetOutputPort())
        vertMapper = vtk.vtkPolyDataMapper()
        vertMapper.SetInputConnection(vertGlyph.GetOutputPort())
        vertMapper.SetImmediateModeRendering(True)
        vertMapper.SetScalarModeToUsePointFieldData()
        vertMapper.SetLookupTable(lut)
        vertMapper.SelectColorArray('ExtBasis')
        vertMapper.Update()
        vertRange = vertMapper.GetInput().GetPointData().GetArray(vertMapper.GetArrayName()).GetRange()
        vertMapper.SetScalarRange(-1.0*vertRange[1], vertRange[1])
        vertActor = vtk.vtkActor()
        vertActor.SetMapper(vertMapper)
        
        outlineMapper = vtk.vtkPolyDataMapper()
        outlineMapper.SetScalarVisibility(False)
        outlineMapper.SetImmediateModeRendering(True)
        outlineMapper.SetInputConnection(vertGlyph.GetOutputPort())
        outlineActor = vtk.vtkActor()
        outlineActor.PickableOff()
        outlineActor.SetPosition(0, 0, -0.001)
        outlineActor.GetProperty().SetRepresentationToWireframe()
        outlineActor.SetMapper(outlineMapper)
        
        # Create an Actor Collection for applying visibility to group
        basisActorCollection = vtk.vtkActorCollection()
        basisActorCollection.AddItem(vertActor)
        # basisActorCollection.AddItem(edgeActor)
        basisActorCollection.AddItem(outlineActor)
        
        
        # Apply a theme to the views
        theme = vtk.vtkViewTheme.CreateMellowTheme()
        theme.SetLineWidth(3)
        theme.SetPointSize(5)
        theme.SetSelectedCellColor(1,1,1)
        theme.SetSelectedPointColor(1,1,1)
        theme.SetOutlineColor(0.8, 0.8, 0.8)
        # theme.SetPointColor(0.9, 0.7, 0.3)
        theme.SetCellColor(0.9, 0.7, 0.3)
        # theme.SetPointOpacity(0.5)
        # theme.SetPointHueRange(0.0, 0.15)
        # theme.SetPointSaturationRange(0.6, 0.8)
        # theme.SetPointValueRange(0.4,0.8)
        # theme.SetPointAlphaRange(0.2,0.8)
        # theme.SetPointAlphaRange(1.0,1.0)

        
        # Apply theme
        # vertActor.GetProperty().SetColor(theme.GetPointColor())
        # vertActor.GetProperty().SetOpacity(theme.GetPointOpacity())
        vertActor.GetProperty().SetPointSize(theme.GetPointSize())
        outlineActor.GetProperty().SetPointSize(vertActor.GetProperty().GetPointSize()+2)
        outlineActor.GetProperty().SetColor(theme.GetOutlineColor())
        outlineActor.GetProperty().SetOpacity(theme.GetPointOpacity())
        edgeActor.GetProperty().SetColor(theme.GetCellColor())
        edgeActor.GetProperty().SetOpacity(theme.GetCellOpacity())
        edgeActor.GetProperty().SetLineWidth(theme.GetLineWidth())
        
        
        # ----------
        # Background graph skeleton
        graphMapper = vtk.vtkGraphMapper()
        graphMapper.SetInputConnection(0, degree.GetOutputPort(0))
        
        # Apply a theme to the background graph
        gtheme = vtk.vtkViewTheme()
        gtheme.SetLineWidth(1)
        gtheme.SetPointSize(0)
        gtheme.SetCellColor(0.8, 0.8, 0.8)
        gtheme.SetCellOpacity(0.2)
        gtheme.SetOutlineColor(0.8, 0.8, 0.8)
        gtheme.SetPointColor(0.8, 0.8, 0.8)
        gtheme.SetPointOpacity(0.0)
        graphMapper.ApplyViewTheme(gtheme)

        graphActor = vtk.vtkActor()
        graphActor.SetMapper(graphMapper)
        graphActor.SetPosition(0,0,-0.005)
        
        # ----------
        # Background vertices
        graphPoly = vtk.vtkGraphToPolyData()
        graphPoly.SetInputConnection(0, tgraph.GetOutputPort(0))
        
        vertGlyph = vtk.vtkGlyph3D()
        vertGlyph.SetInputConnection(0, graphPoly.GetOutputPort())
        glyphSource = vtk.vtkGlyphSource2D()
        glyphSource.SetGlyphTypeToVertex()
        # glyphSource.SetGlyphTypeToCircle()
        # glyphSource.SetScale(0.025)
        vertGlyph.SetInputConnection(1, glyphSource.GetOutputPort())
        
        vertexMapper = vtk.vtkPolyDataMapper()
        vertexMapper.SetInputConnection(vertGlyph.GetOutputPort())
        vertexActor = vtk.vtkActor()
        vertexActor.SetMapper(vertexMapper)
        vertexActor.GetProperty().SetPointSize(4)
        vertexActor.GetProperty().SetOpacity(0.5)
        vertexActor.GetProperty().SetColor(0.6, 0.6, 0.6)
        vertexActor.SetPosition(0, 0, -0.004)
        
        # ----------
        # Vertex index labels
        labelMapper = vtk.vtkDynamic2DLabelMapper()
        labelMapper.SetInputConnection(0, graphPoly.GetOutputPort(0))
        labelMapper.SetLabelModeToLabelFieldData()
        labelMapper.SetFieldDataName("label")
        labelMapper.SetLabelFormat("%s")
        labelMapper.GetLabelTextProperty().SetColor(0.0, 0.0, 0.0)
        
        labelActor = vtk.vtkActor2D()
        labelActor.SetMapper(labelMapper)
        
        # ----------
        # MultiScale Graph
        msGraph = buildSubGraph(level) 
        
        msMapper = vtk.vtkGraphMapper()
        msMapper.SetInput(msGraph)
        msMapper.SetColorEdges(True)
        msMapper.SetEdgeColorArrayName('weight')
        
        # Apply a theme to the background graph
        mtheme = vtk.vtkViewTheme()
        mtheme.SetLineWidth(3)
        mtheme.SetPointSize(11)
        # mtheme.SetCellColor(0.5, 0.5, 0.7)
        # mtheme.SetCellOpacity(0.5)
        mtheme.SetOutlineColor(0.8, 0.8, 0.8)
        mtheme.SetPointColor(0.3, 0.3, 0.6)
        mtheme.SetPointOpacity(1.0)
        mtheme.SetCellHueRange(0.67, 0.67)
        mtheme.SetCellSaturationRange(0.6, 0.1)
        mtheme.SetCellValueRange(0.5,1.0)
        mtheme.SetCellAlphaRange(0.2,0.8)
        msMapper.ApplyViewTheme(mtheme)

        msActor = vtk.vtkActor()
        msActor.SetMapper(msMapper)
        msActor.SetPosition(0,0,-0.002)
        
        # ----------
        # Set up window and add actors        
        view.SetLayoutStrategyToPassThrough()
        # view.ApplyViewTheme(theme)                
        # view.SetupRenderWindow(win)
        view.GetRenderer().SetBackground(theme.GetBackgroundColor())
        view.GetRenderer().SetBackground2(theme.GetBackgroundColor2())
        view.GetRenderer().SetGradientBackground(True)
        view.GetRenderer().AddActor(vertActor)
        view.GetRenderer().AddActor(outlineActor)
        view.GetRenderer().AddActor(edgeActor)
        view.GetRenderer().AddActor(graphActor)
        view.GetRenderer().AddActor(vertexActor)
        view.GetRenderer().AddActor(labelActor)
        view.GetRenderer().AddActor(msActor)
        
        # ----------
        # General interactor
        isty = vtk.vtkInteractorStyleRubberBand2D()
        # RubberBand2D assumes/needs parallel projection ON
        view.GetRenderer().GetActiveCamera().ParallelProjectionOn()
        iren = view.GetRenderWindow().GetInteractor()
        iren.SetInteractorStyle(isty)        
        # Interactor style must be set before scalar bar can be shown
        # view.SetVertexScalarBarVisibility(True)

        sbActor = vtk.vtkScalarBarActor()
        sbActor.SetLookupTable(vertMapper.GetLookupTable())
        sbActor.SetTitle(vertMapper.GetArrayName())
        sbActor.SetNumberOfLabels(3)
        vertexScalarBar = vtk.vtkScalarBarWidget()
        vertexScalarBar.SetScalarBarActor(sbActor)
        vertexScalarBar.SetInteractor(iren)
        vertexScalarBar.SetDefaultRenderer(view.GetRenderer())
        vertexScalarBar.SetCurrentRenderer(view.GetRenderer())
        vertexScalarBar.SetEnabled(True)
        scalarBarRep = vertexScalarBar.GetRepresentation()
        scalarBarRep.SetOrientation(1)  # 0 = Horizontal, 1 = Vertical
        scalarBarRep.GetPositionCoordinate().SetValue(0.05,0.05)
        scalarBarRep.GetPosition2Coordinate().SetValue(0.15,0.25)
        
        # Adding it this way gets it to show up, but it's not interactive
        view.GetRenderer().AddActor(sbActor)
        view.ResetCamera()
        view.Render()
        
        # ----------
        # Add Actors to QListWidget to allow check and uncheck for visibility
        listItem0 = QtGui.QListWidgetItem()
        listItem0.setText('Index Labels')
        listItem0.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        labelActor.SetVisibility(0)
        listItem0.setCheckState(QtCore.Qt.Unchecked)
        # Put actor it in as data in the list widget item
        listItem0.setData(QtCore.Qt.UserRole, QtCore.QVariant(labelActor))
        ui_window.listWidget.insertItem(0,listItem0)
        
        # Test retrieval of actor from list widget item
        # tmpItem = ui_window.listWidget.item(0)
        # tmpQtActor = tmpItem.data(QtCore.Qt.UserRole)
        # tmpActor = tmpQtActor.toPyObject()
        # tmpActor.SetVisibility(0)
        
        # Shorter way to add item to list widget
        listItem1 = QtGui.QListWidgetItem('Vertices (background)', ui_window.listWidget)
        listItem1.setData(QtCore.Qt.UserRole, QtCore.QVariant(vertexActor))
        listItem1.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        listItem1.setCheckState(QtCore.Qt.Checked)

        listItem2 = QtGui.QListWidgetItem('Graph (background)', ui_window.listWidget)
        listItem2.setData(QtCore.Qt.UserRole, QtCore.QVariant(graphActor))
        listItem2.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        listItem2.setCheckState(QtCore.Qt.Checked)

        listItem3 = QtGui.QListWidgetItem()
        listItem3.setText('Basis Function Vertices')
        listItem3.setData(QtCore.Qt.UserRole, QtCore.QVariant(basisActorCollection))
        listItem3.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        listItem3.setCheckState(QtCore.Qt.Checked)
        ui_window.listWidget.insertItem(1,listItem3)
        
        listItem6 = QtGui.QListWidgetItem()
        listItem6.setText('Basis Function Edges')
        listItem6.setData(QtCore.Qt.UserRole, QtCore.QVariant(edgeActor))
        listItem6.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        listItem6.setCheckState(QtCore.Qt.Checked)
        ui_window.listWidget.insertItem(2,listItem6)
        
        listItem4 = QtGui.QListWidgetItem()
        listItem4.setText('MultiScale Graph')
        listItem4.setData(QtCore.Qt.UserRole, QtCore.QVariant(msActor))
        listItem4.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        listItem4.setCheckState(QtCore.Qt.Checked)
        ui_window.listWidget.insertItem(3,listItem4)
                
        listItem5 = QtGui.QListWidgetItem()
        listItem5.setText('Basis Function Scale Bar')
        listItem5.setData(QtCore.Qt.UserRole, QtCore.QVariant(sbActor))
        listItem5.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        listItem5.setCheckState(QtCore.Qt.Checked)
        ui_window.listWidget.insertItem(3,listItem5)
        
        iren.Initialize()
        iren.Start()

    def changeCutoff():
        
        tmpCutoffTxt = ui_window.lineEdit_cutoff.text()
        tmpCutoff = tmpCutoffTxt.toDouble()
        if tmpCutoff[1]:
            cutoff = tmpCutoff[0]
            # Update multi-resolution graph
            msGraph = buildSubGraph(level) 
            msMapper.SetInput(msGraph)
            msMapper.Modified()
            # view.Update()
            view.Render()
            
    def changeBasisCutoff():
        
        tmpCutoffTxt = ui_window.lineEdit_basisCutoff.text()
        tmpCutoff = tmpCutoffTxt.toDouble()
        if tmpCutoff[1]:
            basisCutoff = tmpCutoff[0]
            # Update multi-resolution graph
            selection.RemoveAllThresholds()
            selection.AddThreshold(basisCutoff, 10)
            selection.Modified()
            # view.Update()
            view.Render()
            
    
    def buildSubGraph(level):
        
        table = vtk.vtkTable()
        col0 = vtk.vtkIntArray()
        col0.SetName('index1')
        col1 = vtk.vtkIntArray()
        col1.SetName('index2')
        val = vtk.vtkDoubleArray()
        val.SetName('weight')
        
        Tmat = G.Tree[level,0].T
        minmax = "(%3.2e, %3.2e)" % (abs(Tmat.data).min(), abs(Tmat.data).max())
        ui_window.label_cutoff_minmax.setText(minmax)
        
        ebSum = ExtBasis.sum(0)
        idxs = []
        
        for ii in range(Tmat.nzmax):
            xx = Tmat.rowcol(ii)[0]
            yy = Tmat.rowcol(ii)[1]
            if ((ebSum[0,xx] > ebSumCutoff) and (ebSum[0,yy] > ebSumCutoff)):
                # print ebSum[0,xx], ebSum[0,yy]
                # Always include diagonal elements (no matter what cutoff) so all vertices are created
                if ((abs(Tmat.getdata(ii)) > cutoff) | (Tmat.rowcol(ii)[0]==Tmat.rowcol(ii)[1])):
                    col0.InsertNextValue(xx)
                    col1.InsertNextValue(yy)
                    # But, take out -loops for edge color scale
                    if (xx == yy):
                        val.InsertNextValue(0.0)        
                        idxs.append(G.Tree[level,0].ExtIdxs[xx]-1)
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
        ptsArray = VN.vtk_to_numpy(pts.GetData())
        # Remember Matlab uses one-based indices (ExtIdxs-1 for Numpy)
        # print (G.Tree[level,0].ExtIdxs-1)
        # print ptsArray[G.Tree[level,0].ExtIdxs-1,:]
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
         
    def basisUpdate(index):
        
        QtCore.QObject.disconnect(ui_window.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
        QtCore.QObject.disconnect(ui_window.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
        
        basisNum = index
        ui_window.spinBox_basisIndex.setValue(index)
        ui_window.hSlider_basisIndex.setValue(index)
        selection.RemoveAllIDs()
        
        # Try to directly swap ExtBasis vertex data
        # print 'ExtBasis shape: ' + str(ExtBasis[:,basisNum].data.shape) + ' (level,basis) (' + str(level) + ',' + str(basisNum) + ')'
        
        # Indexing of ExtBasis is backwards from graph, so need to reverse with [::-1]
        # and array not "contiguous" if don't do .copy()
        Esub[:] = ExtBasis[:,basisNum].data[::-1].copy()
        EsubSq[:] = Esub**2
        # Esub[:] = N.random.random(ExtBasis[:,basisNum].data.shape[0])
        # SubIdxs = (Esub > 0.001).nonzero()
        # SubIdxs = (Esub**2 > 0.8).nonzero()
        vertexData.Modified()
        degree.Update()
                
        # for ii in SubIdxs[0]:
        #     selection.AddID(0,ii)
        vertRange = vertMapper.GetInput().GetPointData().GetArray(vertMapper.GetArrayName()).GetRange()
        vertMapper.SetScalarRange(-1.0*vertRange[1], vertRange[1])
        selection.Modified()
        # tmpMin = EsubSq.min()
        # tmpMax = EsubSq.max()
        # map0.SetScalarRange(-100.0*tmpMax, 100.0*tmpMax)
        minmax = "(%3.2e, %3.2e)" % (EsubSq.min(), EsubSq.max())
        ui_window.label_basisCutoff_minmax.setText(minmax)
        # view.Update()
        view.Render()
        
        QtCore.QObject.connect(ui_window.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
        QtCore.QObject.connect(ui_window.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
        
    def levelUpdate(index):
        
        QtCore.QObject.disconnect(ui_window.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
        QtCore.QObject.disconnect(ui_window.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
        
        level = index
        ui_window.spinBox_level.setValue(index)
        ui_window.hSlider_level.setValue(index)
        # Start setting up basis function subgraph
        ExtBasis = G.Tree[level,0].ExtBasis
        basisMax = ExtBasis.shape[1]-1    # zero-based indices
        if basisNum > basisMax:
            basisNum = basisMax
            ui_window.hSlider_basisIndex.setValue(basisMax)
            ui_window.spinBox_basisIndex.setValue(basisMax)
        ui_window.hSlider_basisIndex.setMinimum(0)
        ui_window.hSlider_basisIndex.setMaximum(basisMax)
        ui_window.spinBox_basisIndex.setMinimum(0)
        ui_window.spinBox_basisIndex.setMaximum(basisMax)

        # selection.RemoveAllIDs()
        # print 'ExtBasis shape: ' + str(ExtBasis[:,basisNum].data.shape) + ' (level,basis) (' + str(level) + ',' + str(basisNum) + ')'
        
        # Indexing of ExtBasis is backwards from graph, so need to reverse with [::-1]
        # and array not "contiguous" if don't do .copy()
#         Esub[:] = ExtBasis[:,basisNum].data[::-1].copy()
#         EsubSq[:] = Esub**2
#         # Esub[:] = N.random.random(ExtBasis[:,basisNum].data.shape[0])
#         # SubIdxs = (Esub > 0.001).nonzero()
#         # SubIdxs = (Esub**2 > 0.8).nonzero()
#         vertexData.Modified()
#         degree.Update()
        
        # Update multi-resolution graph
        msGraph = buildSubGraph(index) 
        msMapper.SetInput(msGraph)
        msMapper.Modified()

        # for ii in SubIdxs[0]:
        #     selection.AddID(0,ii)
        selection.Modified()
        # tmpMin = Esub.min()
        # tmpMax = Esub.max()
        # map0.SetScalarRange(-100.0*tmpMax, 100.0*tmpMax)
        basisUpdate(ui_window.hSlider_basisIndex.value())
        
        QtCore.QObject.connect(ui_window.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
        QtCore.QObject.connect(ui_window.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
        
    def toggleVisibility(item):
        
        tmpQtActor = item.data(QtCore.Qt.UserRole)
        tmpActor = tmpQtActor.toPyObject()
        if tmpActor.GetClassName() == 'vtkActorCollection':
            # Loop through actor collection
            tmpActor.InitTraversal()
            singleActor = tmpActor.GetNextActor()
            while singleActor is not None:
                toggleActorVisibility(item,singleActor)
                singleActor = tmpActor.GetNextActor()
        else:
            # Single actor
            toggleActorVisibility(item,tmpActor)
            
    def toggleActorVisibility(item, tmpActor):
        actorVis = tmpActor.GetVisibility()
        itemCheck = item.checkState()   # QtCore.Qt.Checked = 2, QtCore.Qt.Unchecked = 0
        if (actorVis == 0) and (itemCheck == QtCore.Qt.Checked):
            tmpActor.SetVisibility(1)
            # view.Update()
            view.Render()
        elif (actorVis == 1) and (itemCheck == QtCore.Qt.Unchecked):
            tmpActor.SetVisibility(0)
            # view.Update()
            view.Render()
        
        
    def fileExit():
        
        # Usually would use the qApp global variable qApp.quit(), but wasn't working...
        QtGui.QApplication.instance().quit()
        
