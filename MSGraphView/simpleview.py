"""
simpleview.py

Example using PyQt4, Qt Designer & VTK Infovis
Plus M Maggioni multiscale graph representation

09 Jan 2009 -- E Monson

Modified to no use QVTKWidget, but open separate vtkGraphLayoutView & Qt control GUI
Also changed way selection is done on subgraph for thresholding ExtBasisSq values
04 Dec 2009 -- EMonson

"""

from PyQt4 import QtCore, QtGui
from PyQt4 import QtCore
from PyQt4 import QtGui
from PyQt4 import uic

import vtk
import vtk.util.numpy_support as VN
import numpy as N
import scipy.io
import os
import sys

# Create a PyQt window using a .ui file generated with Qt Designer ...
application = QtGui.QApplication(sys.argv)

ui_window = uic.loadUi("qsimpleview.ui")

view = vtk.vtkGraphLayoutView()
view.GetRenderWindow().SetSize(600, 600)

level = 3
minMSLevel = 0
basisNum = 0
cutoff = 1.7e-5
basisCutoff = 1.0e-4
ebSumCutoff = 0.02

# ----------
# Load and construct whole graph and multi-resolution data from Matlab structure
dataDir = '/Users/emonson/Data/Fodava/EMoDocMatlabData/'
# filename = dataDir + 'X20_042709b.mat'
filename = dataDir + 'SN_mauro_TFcorr_121509.mat'
# filename = dataDir + 'n20_set5train_TFcorr_121009.mat'

print 'Loading data set from .mat file...'
X = scipy.io.loadmat(filename,struct_as_record=True)
# Get graph structure G out of matlab variables
G = X['G']
classes = X['classes']
del X

GTree = G['Tree'][0,0]
T = G['T'][0][0]

# ----------
# Set multi-resolution level bounds in GUI sliders
levelMax = GTree.shape[0]-1
ui_window.hSlider_level.setMinimum(1)
ui_window.hSlider_level.setMaximum(levelMax)
ui_window.spinBox_level.setMinimum(1)
ui_window.spinBox_level.setMaximum(levelMax)
ui_window.lineEdit_basisCutoff.setText("%2.1e" % basisCutoff)

# Start setting up basis function for display as subgraph
# ExtBasis = GTree[level,0]['ExtBasis'][0][0]
ExtIdxs = GTree[level,0]['ExtIdxs'][0][0][0]
# TODO: Maybe limit the number of basis functions to have the same EbSumCutoff as level T
basisMax = GTree[level,0]['ExtBasis'][0][0].shape[1]-1    # zero-based indices

# Set particular level basis function bounds in GUI sliders
ui_window.hSlider_basisIndex.setMinimum(0)
ui_window.hSlider_basisIndex.setMaximum(basisMax)
ui_window.spinBox_basisIndex.setMinimum(0)
ui_window.spinBox_basisIndex.setMaximum(basisMax)

print 'Building full scale 0 graph...'

Tdata = abs(T.data)
Trows = T.nonzero()[0]
Tcols = T.nonzero()[1]
Tdata.dtype = 'float64'
Trows.dtype = 'int32'
Tcols.dtype = 'int32'

# Set all self-loops to zero weight
Tdata[N.nonzero(Trows==Tcols)[0]] = 0.0

col0 = VN.numpy_to_vtk(Trows)
col0.SetName('index1')
col1 = VN.numpy_to_vtk(Tcols)
col1.SetName('index2')
val = VN.numpy_to_vtk(Tdata)
val.SetName('weight')

# Make set of indices that made it through the T cutoff filter
# VN.vtk_to_numpy(table.GetColumnByName('index1')).tolist()
s1 = set(Trows)
s2 = set(Tcols)
ss = s1.union(s2)
# Then figure out which indices (vertices) from the thresholded matrix haven't made it through
s_skipped = set(range(T.shape[0])).difference(ss)

# And add these skipped (isolated) vertices back in
for index in s_skipped:
	col0.InsertNextValue(index)
	col1.InsertNextValue(index)
	val.InsertNextValue(0.0)

table = vtk.vtkTable()	
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
# if os.path.exists(ptsFile):
polyreader = vtk.vtkXMLPolyDataReader()
polyreader.SetFileName(ptsFile)
polyreader.Update()
pts = polyreader.GetOutput().GetPoints()
rawGraph.SetPoints(pts)
# print pts

# Grab points from original set with ExtIdxs
ptsArray = VN.vtk_to_numpy(pts.GetData())
ptsMatrix = N.matrix(ptsArray)


def changeTAlphaMinValue(value):
		
	newMin = float(value)/400.0
	msEdgeMapper.GetEdgeLookupTable().SetAlphaRange(newMin,1.0)
	view.Render()

def changeExtBasisAlphaMinValue(value):
		
	newMin = float(value)/400.0
	edgeMapper.GetLookupTable().SetAlphaRange(newMin,0.8)
	view.Render()

def changeGraphAlphaMinValue(value):
		
	newMin = float(value)/400.0
	graphMapper.GetEdgeLookupTable().SetAlphaRange(newMin,0.5)
	view.Render()

def changeBasisCutoff():

	global basisCutoff, ui_window, selArr
	
	tmpCutoffTxt = ui_window.lineEdit_basisCutoff.text()
	tmpCutoff = tmpCutoffTxt.toDouble()
	if tmpCutoff[1]:
		basisCutoff = tmpCutoff[0]
		# Reset the threshold minimum in original array directly
		selArr.SetValue(0,basisCutoff)
		selection.Modified()
		view.Render()
		

def buildSubGraph(level):
	
	Tmat_orig = GTree[level,0]['T'][0][0][0][0]
	Current_ExtIdxs_orig = GTree[level,0]['ExtIdxs'][0][0][0]
	
	# Only pass rows and columns of T which have large enough sum (L1 norm)
	ebSum = GTree[level,0]['ExtBasis'][0][0].sum(0)
	# ebSum is e.g. (1,1046) matrix 
	# ebSum.A[0] is e.g. (1046,) array
	# Convert ebSum to 1D array (.A[0]) from 2D matrix to get this to work right for indexing matrices
	over = N.nonzero(ebSum.A[0]>ebSumCutoff)[0]
	Tmat = Tmat_orig[N.ix_(over,over)]
	Current_ExtIdxs = Current_ExtIdxs_orig[over]
	# del Tmat_orig, Current_ExtIdxs_orig
	
	Tdata_orig = abs(Tmat.data)
	Trows_orig = Tmat.nonzero()[0]
	Tcols_orig = Tmat.nonzero()[1]
	
	tover = N.nonzero(Tdata_orig>cutoff)[0]
	Tdata = Tdata_orig[tover]
	Tdata.dtype = 'float64'
	Trows = Trows_orig[tover]
	Trows.dtype = 'int32'		# was getting TypeError on numpy_to_vtk if didn't set this explicitly...
	Tcols = Tcols_orig[tover]
	Tcols.dtype = 'int32'
	# del Tdata_orig, Trows_orig, Tcols_orig
	# Set all self-loops to zero weight
	Tdata[N.nonzero(Trows==Tcols)[0]] = 0.0
		
	col0a = VN.numpy_to_vtk(Trows)
	col0a.SetName('index1')
	col1a = VN.numpy_to_vtk(Tcols)
	col1a.SetName('index2')
	vala = VN.numpy_to_vtk(Tdata)
	vala.SetName('weight')
	
	# Make set of indices that made it through the T cutoff filter
	# VN.vtk_to_numpy(table.GetColumnByName('index1')).tolist()
	s1 = set(Trows)
	s2 = set(Tcols)
	ss = s1.union(s2)
	# Then figure out which indices (vertices) from the thresholded matrix haven't made it through
	s_skipped = set(range(Current_ExtIdxs.shape[0])).difference(ss)
	
	# And add these skipped (isolated) vertices back in
	for index in s_skipped:
		col0a.InsertNextValue(index)
		col1a.InsertNextValue(index)
		vala.InsertNextValue(0.0)
	
	table = vtk.vtkTable()	
	table.AddColumn(col0a)
	table.AddColumn(col1a)
	table.AddColumn(vala)
		
	tgraph = vtk.vtkTableToGraph()
	tgraph.SetInput(table)
	tgraph.AddLinkVertex('index2', 'stuff', False)
	tgraph.AddLinkVertex('index1', 'stuff', False)
	tgraph.AddLinkEdge('index2', 'index1')
	tgraph.Update()
	
	rawGraph = tgraph.GetOutputDataObject(0)
	# rawGraph.Update()
	
	# Remember Matlab uses one-based indices (ExtIdxs-1 for Numpy)
	ptsSubMatrix = ptsMatrix[Current_ExtIdxs-1,:]
	ptsSub = VN.numpy_to_vtk(ptsSubMatrix)
	pts = vtk.vtkPoints()
	pts.SetData(ptsSub)
	rawGraph.SetPoints(pts)
	rawGraph.Update()
			
	strategy = vtk.vtkPassThroughLayoutStrategy()
	layout = vtk.vtkGraphLayout()
	layout.SetInput(rawGraph)
	layout.SetLayoutStrategy(strategy)
	
	edgeLayout = vtk.vtkEdgeLayout()
	edgeStrategy = vtk.vtkArcParallelEdgeStrategy()
	edgeStrategy.SetNumberOfSubdivisions(10)
	edgeLayout.SetInputConnection(layout.GetOutputPort())
	edgeLayout.SetLayoutStrategy(edgeStrategy)
	edgeLayout.Update()
	
	return edgeLayout.GetOutputDataObject(0)
	 
def basisUpdate(index):

	global basisNum, Esub, EsubSq, minmax, vertRange, nodeGraph, level
	
	ui_window.hSlider_basisIndex.blockSignals(True)
	ui_window.hSlider_level.blockSignals(True)
	ui_window.spinBox_basisIndex.blockSignals(True)
	
	basisNum = index
	ui_window.spinBox_basisIndex.setValue(index)
	ui_window.hSlider_basisIndex.setValue(index)
	
	# Here need to deal with the fact that ExtBasis can be sparse, so using .toarray()[:,0] rather than .data
	Esub[:] = GTree[level,0]['ExtBasis'][0][0][:,basisNum].toarray()[:,0]
	EsubSq[:] = Esub**2
	vertexData.Modified()
			
# 	if (vertMapper.GetInput().GetPointData().GetNumberOfTuples() > 0):
# 		vertRange = vertMapper.GetInput().GetPointData().GetArray(vertMapper.GetArrayName()).GetRange()
# 	else:
# 		vertRange = (0, 0.1)
# 	maxabs = max(abs(N.array(vertRange)))
	maxabs = max(abs(Esub))
	vertMapper.SetScalarRange(-1.0*maxabs, maxabs)
	vertMapper.Modified()
	selection.Modified()
	
	# Update current single node T graph indicator
	# Remember Matlab uses one-based indices (ExtIdxs-1 for Numpy)
	currIdx = ExtIdxs[basisNum]-1
	currPt = VN.numpy_to_vtk(ptsMatrix[currIdx,:])
	nodeGraph.GetPoints().SetData(currPt)
	nodeGraph.Modified()

	minmax = "(%3.2e, %3.2e)" % (EsubSq.min(), EsubSq.max())
	ui_window.label_basisCutoff_minmax.setText(minmax)
	view.Render()
	
	ui_window.hSlider_basisIndex.blockSignals(False)
	ui_window.hSlider_level.blockSignals(False)
	ui_window.spinBox_basisIndex.blockSignals(False)
		
def levelUpdate(index):

	global basisNum, ExtIdxs, basisMax, level
	
	changeLevel = False
	
	ui_window.hSlider_basisIndex.blockSignals(True)
	ui_window.hSlider_level.blockSignals(True)
	ui_window.spinBox_basisIndex.blockSignals(True)
	ui_window.spinBox_level.blockSignals(True)
	
	PrevIdx = ExtIdxs[basisNum]
	TmpNewExtIdxs = GTree[index,0]['ExtIdxs'][0][0][0]
		
	# Try to find same index in next ExtIdxs set so can track ExtBasis fct across scales
	# For now, don't let increase MS level if basis function doesn't exist in new set...
	newBasisNumArr = (TmpNewExtIdxs==PrevIdx).nonzero()[0]
	if (newBasisNumArr.size > 0):
		changeLevel = True
		if (newBasisNumArr[0]!=basisNum):
			basisNum = newBasisNumArr[0]
			ui_window.hSlider_basisIndex.setValue(basisNum)
			ui_window.spinBox_basisIndex.setValue(basisNum)

	if changeLevel:
		ui_window.spinBox_level.setValue(index)
		ui_window.hSlider_level.setValue(index)
		# Start setting up basis function subgraph
		level = index
		# ExtBasis = GTree[level,0]['ExtBasis'][0][0]
		ExtIdxs = GTree[level,0]['ExtIdxs'][0][0][0]
	
		basisMax = GTree[level,0]['ExtBasis'][0][0].shape[1]-1    # zero-based indices
		if basisNum > basisMax:
			basisNum = basisMax
			ui_window.hSlider_basisIndex.setValue(basisMax)
			ui_window.spinBox_basisIndex.setValue(basisMax)
		ui_window.hSlider_basisIndex.setMinimum(0)
		ui_window.hSlider_basisIndex.setMaximum(basisMax)
		ui_window.spinBox_basisIndex.setMinimum(0)
		ui_window.spinBox_basisIndex.setMaximum(basisMax)
		
		# Update multi-resolution graph
		# msGraph = buildSubGraph(index) 
		msEdgeMapper.SetInput(msGraphList[level])
		msEdgeMapper.Modified()
		msVertGlyphs.SetInput(msGraphList[level])
		msVertGlyphs.Modified()
	
		basisUpdate(ui_window.hSlider_basisIndex.value())
	else:
		print 'Not allowing level increase because of basis function mismatch...'
		ui_window.spinBox_level.setValue(level)
		ui_window.hSlider_level.setValue(level)
	
	ui_window.hSlider_basisIndex.blockSignals(False)
	ui_window.hSlider_level.blockSignals(False)
	ui_window.spinBox_basisIndex.blockSignals(False)
	ui_window.spinBox_level.blockSignals(False)
	
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
		view.Render()
	elif (actorVis == 1) and (itemCheck == QtCore.Qt.Unchecked):
		tmpActor.SetVisibility(0)
		view.Render()
	
	
def fileExit():
	
	# Usually would use the qApp global variable qApp.quit(), but wasn't working...
	QtGui.QApplication.instance().quit()
	



	
strategy = vtk.vtkPassThroughLayoutStrategy()
layout = vtk.vtkGraphLayout()
layout.SetInput(rawGraph)
layout.SetLayoutStrategy(strategy)

edgeLayout = vtk.vtkEdgeLayout()
# edgeStrategy = vtk.vtkPassThroughEdgeStrategy()
edgeStrategy = vtk.vtkArcParallelEdgeStrategy()
edgeStrategy.SetNumberOfSubdivisions(10)
edgeLayout.SetInputConnection(layout.GetOutputPort())
edgeLayout.SetLayoutStrategy(edgeStrategy)

graph = edgeLayout.GetOutput()
graph.Update()

print 'Building support structures...'
# --------
# Add ExtBasis to graph data & Select particular basis function
# Here need to deal with the fact that ExtBasis can be sparse, so using .toarray()[:,0] rather than .data
Esub = GTree[level,0]['ExtBasis'][0][0][:,basisNum].toarray()[:,0]
EsubSq = Esub**2

# Set ExtBasis vertex data from numpy array
basisFunc = VN.numpy_to_vtk(Esub)
basisFunc.SetName('ExtBasis')
basisFuncSq = VN.numpy_to_vtk(EsubSq)
basisFuncSq.SetName('ExtBasisSq')

vertexData = graph.GetVertexData()
vertexData.AddArray(basisFunc)
vertexData.AddArray(basisFuncSq)

# Do selection the new way rather than using vtkSelectionSource...
selArr = vtk.vtkDoubleArray()
selArr.InsertNextValue(basisCutoff)	# min
selArr.InsertNextValue(100)			# max
selArr.SetName("ExtBasisSq")		# array name matches vertex array name

selNode = vtk.vtkSelectionNode()
selNode.SetContentType(7) 			# vtkSelectionNode::THRESHOLDS
selNode.SetFieldType(3) 			# vtkSelectionNode::VERTEX
selNode.SetSelectionList(selArr)

selection = vtk.vtkSelection()
selection.AddNode(selNode)
selection.Update()

minmax = "(%3.2e, %3.2e)" % (EsubSq.min(), EsubSq.max())
ui_window.label_basisCutoff_minmax.setText(minmax)

# ----------
# Back to pipeline
subgraph = vtk.vtkExtractSelectedGraph()
subgraph.SetRemoveIsolatedVertices(False)
subgraph.SetInput(0,graph)
subgraph.SetInput(1,selection)

# Apply a theme to the views
theme = vtk.vtkViewTheme.CreateMellowTheme()
theme.SetLineWidth(3)
theme.SetPointSize(5)
theme.SetSelectedCellColor(1,1,1)
theme.SetSelectedPointColor(1,1,1)
theme.SetOutlineColor(0.8, 0.8, 0.8)
# theme.SetPointColor(0.9, 0.7, 0.3)
# theme.SetCellColor(0.9, 0.7, 0.3)		# orig
theme.SetCellHueRange(0.11111, 0.11111)
theme.SetCellSaturationRange(0.66667, 0.66667)
theme.SetCellValueRange(0.8, 0.8)
theme.SetCellAlphaRange(0.075, 0.8)
ui_window.hSlider_ExtBasisAlpha.setValue(30)	# corresponds to min alpha 0.075
# theme.SetPointOpacity(0.5)
# theme.SetPointHueRange(0.0, 0.15)
# theme.SetPointSaturationRange(0.6, 0.8)
# theme.SetPointValueRange(0.4,0.8)
# theme.SetPointAlphaRange(0.2,0.8)
# theme.SetPointAlphaRange(1.0,1.0)

# +++++++++++++
graphToPoly = vtk.vtkGraphToPolyData()
graphToPoly.SetInputConnection(subgraph.GetOutputPort())
edgeMapper = vtk.vtkPolyDataMapper()
edgeMapper.SetInputConnection(graphToPoly.GetOutputPort())
edgeMapper.SetScalarVisibility(True)
edgeMapper.SetScalarModeToUseCellFieldData()
edgeMapper.SelectColorArray('weight')
edgeMapper.SetLookupTable(theme.GetCellLookupTable())
edgeMapper.SetImmediateModeRendering(True)
edgeActor = vtk.vtkActor()
edgeActor.SetMapper(edgeMapper)
# edgeActor.GetProperty().SetColor(theme.GetCellColor())
# edgeActor.GetProperty().SetOpacity(theme.GetCellOpacity())
edgeActor.GetProperty().SetLineWidth(theme.GetLineWidth())
edgeActor.SetPosition(0, 0, -0.003);
edgeActor.SetVisibility(0)

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
# vertActor.GetProperty().SetColor(theme.GetPointColor())
# vertActor.GetProperty().SetOpacity(theme.GetPointOpacity())
vertActor.GetProperty().SetPointSize(theme.GetPointSize())

outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetScalarVisibility(False)
outlineMapper.SetImmediateModeRendering(True)
outlineMapper.SetInputConnection(vertGlyph.GetOutputPort())
outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.PickableOff()
outlineActor.GetProperty().SetRepresentationToWireframe()
outlineActor.GetProperty().SetPointSize(vertActor.GetProperty().GetPointSize()+2)
outlineActor.GetProperty().SetColor(theme.GetOutlineColor())
outlineActor.GetProperty().SetOpacity(theme.GetPointOpacity())
outlineActor.SetPosition(0, 0, -0.001)

# ----------
# Current ExtIdxs MS Graph Node
nodeGraph = vtk.vtkMutableUndirectedGraph()
nodeGraph.AddVertex()
currIdx = ExtIdxs[basisNum]-1	# Matlab 1-based and Numpy 0-based indexing
currPt = VN.numpy_to_vtk(ptsMatrix[currIdx,:])
nodeGraph.GetPoints().SetData(currPt)

ntheme = vtk.vtkViewTheme()
ntheme.SetPointSize(13)
ntheme.SetOutlineColor(0.2, 0.2, 0.2)
ntheme.SetPointColor(0.1, 0.1, 0.1)
ntheme.SetPointOpacity(1.0)

nodeMapper = vtk.vtkGraphMapper()
nodeMapper.SetInput(nodeGraph)
nodeMapper.SetColorEdges(False)
nodeMapper.ApplyViewTheme(ntheme)

nodeActor = vtk.vtkActor()
nodeActor.SetMapper(nodeMapper)
nodeActor.SetPosition(0,0,-0.0025)

# ----------
# Create an Actor Collection for applying visibility to group
basisActorCollection = vtk.vtkActorCollection()
basisActorCollection.AddItem(vertActor)
basisActorCollection.AddItem(outlineActor)

# ----------
# Background graph skeleton
graphMapper = vtk.vtkGraphMapper()
graphMapper.SetInputConnection(0, edgeLayout.GetOutputPort(0))
graphMapper.SetColorEdges(True)
graphMapper.SetEdgeColorArrayName('weight')

# Apply a theme to the background graph
gtheme = vtk.vtkViewTheme()
gtheme.SetLineWidth(1)
gtheme.SetPointSize(0)
gtheme.SetOutlineColor(0.8, 0.8, 0.8)
gtheme.SetPointColor(0.8, 0.8, 0.8)
gtheme.SetPointOpacity(0.0)
gtheme.SetCellHueRange(0.0, 0.0)
gtheme.SetCellSaturationRange(0.0, 0.0)
gtheme.SetCellValueRange(1.0, 1.0)
gtheme.SetCellAlphaRange(0.1,0.5)
ui_window.hSlider_GraphAlpha.setValue(40)	# corresponds to min alpha 0.1
graphMapper.ApplyViewTheme(gtheme)

graphActor = vtk.vtkActor()
graphActor.SetMapper(graphMapper)
graphActor.SetPosition(0,0,-0.005)
graphActor.SetVisibility(0)

# ----------
# Background vertices
vertGlyph = vtk.vtkVertexGlyphFilter()
vertGlyph.SetInputConnection(0, edgeLayout.GetOutputPort(0))

vertexMapper = vtk.vtkPolyDataMapper()
vertexMapper.SetInputConnection(vertGlyph.GetOutputPort(0))
vertexActor = vtk.vtkActor()
vertexActor.SetMapper(vertexMapper)
vertexActor.GetProperty().SetPointSize(4)
vertexActor.GetProperty().SetOpacity(0.5)
vertexActor.GetProperty().SetColor(0.6, 0.6, 0.6)
vertexActor.SetPosition(0, 0, -0.004)

# ----------
# Vertex index labels
labelMapper = vtk.vtkDynamic2DLabelMapper()
# labelMapper.SetInputConnection(0, graphPoly.GetOutputPort(0))		# label background vertices
labelMapper.SetInputConnection(0, subgraph.GetOutputPort(0))	# label basis function vertices
labelMapper.SetLabelModeToLabelFieldData()
labelMapper.SetFieldDataName("label")
labelMapper.SetLabelFormat("%s")
labelMapper.GetLabelTextProperty().SetColor(0.0, 0.0, 0.0)

labelActor = vtk.vtkActor2D()
labelActor.SetMapper(labelMapper)

# ----------
# MultiScale Graph Edges
print 'Building all multi-scale graphs...'
msGraphList = range(minMSLevel)	# will not be using 0 or 1 element for now
for ll in range(minMSLevel,levelMax+1):
	print ll
	msGraphList.append(buildSubGraph(ll))

print 'Continuing with rest of render pipeline setup...'
msEdgeMapper = vtk.vtkGraphMapper()
msEdgeMapper.SetInput(msGraphList[level])
msEdgeMapper.SetColorEdges(True)
msEdgeMapper.SetEdgeColorArrayName('weight')

# Apply a theme to the background graph
mtheme = vtk.vtkViewTheme()
mtheme.SetLineWidth(3)
mtheme.SetPointSize(0)
# mtheme.SetCellColor(0.5, 0.5, 0.7)
# mtheme.SetCellOpacity(0.5)
mtheme.SetOutlineColor(0.8, 0.8, 0.8)
mtheme.SetPointColor(0.3, 0.3, 0.6)
mtheme.SetPointOpacity(1.0)
mtheme.SetCellHueRange(0.67, 0.67)
mtheme.SetCellSaturationRange(0.6, 0.6)	# original (0.6, 0.1)
mtheme.SetCellValueRange(0.5,0.5)		# original (0.5, 1.0)
mtheme.SetCellAlphaRange(0.2,1.0)
ui_window.hSlider_TAlpha.setValue(80)	# corresponds to min alpha 0.2
msEdgeMapper.ApplyViewTheme(mtheme)

msEdgeActor = vtk.vtkActor()
msEdgeActor.SetMapper(msEdgeMapper)
msEdgeActor.SetPosition(0,0,-0.002)
msEdgeActor.SetVisibility(0)

# ----------
# MultiScale Graph Vertices
msVertGlyphs = vtk.vtkVertexGlyphFilter()
msVertGlyphs.SetInput(msGraphList[level])

msVertMapper = vtk.vtkPolyDataMapper()
msVertMapper.SetInputConnection(msVertGlyphs.GetOutputPort(0))
msVertActor = vtk.vtkActor()
msVertActor.SetMapper(msVertMapper)
msVertActor.GetProperty().SetPointSize(11)
msVertActor.GetProperty().SetOpacity(1.0)
msVertActor.GetProperty().SetColor(0.3, 0.3, 0.6)
msVertActor.SetPosition(0, 0, -0.0015)

# ----------
# Set up window and add actors        
view.SetLayoutStrategyToPassThrough()
view.GetRenderer().SetBackground(theme.GetBackgroundColor())
view.GetRenderer().SetBackground2(theme.GetBackgroundColor2())
view.GetRenderer().SetGradientBackground(True)
view.GetRenderer().AddActor(vertActor)
view.GetRenderer().AddActor(nodeActor)
view.GetRenderer().AddActor(outlineActor)
view.GetRenderer().AddActor(edgeActor)
view.GetRenderer().AddActor(graphActor)
view.GetRenderer().AddActor(vertexActor)
view.GetRenderer().AddActor(labelActor)
view.GetRenderer().AddActor(msEdgeActor)
view.GetRenderer().AddActor(msVertActor)

# ----------
# General interactor
isty = vtk.vtkInteractorStyleRubberBand2D()
# RubberBand2D assumes/needs parallel projection ON
view.GetRenderer().GetActiveCamera().ParallelProjectionOn()
iren = view.GetRenderWindow().GetInteractor()
iren.SetInteractorStyle(isty)        

vertexScalarBar = vtk.vtkScalarBarWidget()
vertexScalarBar.SetInteractor(iren)
vertexScalarBar.SetDefaultRenderer(view.GetRenderer())
vertexScalarBar.SetCurrentRenderer(view.GetRenderer())
vertexScalarBar.SetEnabled(True)
sbActor = vertexScalarBar.GetScalarBarActor()
sbActor.SetLookupTable(vertMapper.GetLookupTable())
sbActor.SetTitle(vertMapper.GetArrayName())
sbActor.SetNumberOfLabels(3)
scalarBarRep = vertexScalarBar.GetRepresentation()
scalarBarRep.SetOrientation(1)  # 0 = Horizontal, 1 = Vertical
scalarBarRep.GetPositionCoordinate().SetValue(0.05,0.05)
scalarBarRep.GetPosition2Coordinate().SetValue(0.15,0.25)

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
listItem2.setCheckState(QtCore.Qt.Unchecked)

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
listItem6.setCheckState(QtCore.Qt.Unchecked)
ui_window.listWidget.insertItem(2,listItem6)

listItem4 = QtGui.QListWidgetItem()
listItem4.setText('MS Graph Edges')
listItem4.setData(QtCore.Qt.UserRole, QtCore.QVariant(msEdgeActor))
listItem4.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
listItem4.setCheckState(QtCore.Qt.Unchecked)
ui_window.listWidget.insertItem(3,listItem4)
		
listItem4a = QtGui.QListWidgetItem()
listItem4a.setText('MS Graph Vertices')
listItem4a.setData(QtCore.Qt.UserRole, QtCore.QVariant(msVertActor))
listItem4a.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
listItem4a.setCheckState(QtCore.Qt.Checked)
ui_window.listWidget.insertItem(3,listItem4a)
		
listItem4a = QtGui.QListWidgetItem()
listItem4a.setText('MS Basis Node')
listItem4a.setData(QtCore.Qt.UserRole, QtCore.QVariant(nodeActor))
listItem4a.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
listItem4a.setCheckState(QtCore.Qt.Checked)
ui_window.listWidget.insertItem(3,listItem4a)
		
listItem5 = QtGui.QListWidgetItem()
listItem5.setText('Scale Bar (basis function)')
listItem5.setData(QtCore.Qt.UserRole, QtCore.QVariant(scalarBarRep))
listItem5.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
listItem5.setCheckState(QtCore.Qt.Checked)
ui_window.listWidget.insertItem(3,listItem5)


basisCutoffValid = QtGui.QDoubleValidator(0.0, 1.0, 12, ui_window.lineEdit_basisCutoff)
ui_window.lineEdit_basisCutoff.setValidator(basisCutoffValid)

# Connect signals and slots
# QtCore.QObject.connect(ui_window.actionNew, QtCore.SIGNAL("triggered()"), fileOpen)
QtCore.QObject.connect(ui_window.actionExit, QtCore.SIGNAL("triggered()"), fileExit)
QtCore.QObject.connect(ui_window.hSlider_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
QtCore.QObject.connect(ui_window.spinBox_basisIndex, QtCore.SIGNAL("valueChanged(int)"), basisUpdate)
QtCore.QObject.connect(ui_window.hSlider_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
QtCore.QObject.connect(ui_window.spinBox_level, QtCore.SIGNAL("valueChanged(int)"), levelUpdate)
QtCore.QObject.connect(ui_window.listWidget, QtCore.SIGNAL("itemClicked(QListWidgetItem *)"), toggleVisibility)
QtCore.QObject.connect(ui_window.lineEdit_basisCutoff, QtCore.SIGNAL("editingFinished()"), changeBasisCutoff)

QtCore.QObject.connect(ui_window.hSlider_TAlpha, QtCore.SIGNAL("valueChanged(int)"), changeTAlphaMinValue)
QtCore.QObject.connect(ui_window.hSlider_ExtBasisAlpha, QtCore.SIGNAL("valueChanged(int)"), changeExtBasisAlphaMinValue)
QtCore.QObject.connect(ui_window.hSlider_GraphAlpha, QtCore.SIGNAL("valueChanged(int)"), changeGraphAlphaMinValue)

ui_window.hSlider_basisIndex.setValue(basisNum)
ui_window.spinBox_basisIndex.setValue(basisNum)
ui_window.hSlider_level.setValue(level)
ui_window.spinBox_level.setValue(level)

# TODO: Problem with size of ExtBasis for levels below 2...
ui_window.hSlider_level.setMinimum(minMSLevel)
ui_window.spinBox_level.setMinimum(minMSLevel)

ui_window.show()


# This initializes the VTK window for interaction, but doesn't start an event-loop ...
view.GetRenderWindow().Start()

# Start the Qt event-loop ...
application.exec_()

	

