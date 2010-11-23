import vtk
import vtk.util.numpy_support as VN
import numpy as N
import sys
# sys.path.append("/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
# import libvtkvtgChartsPython as vtgCh
import vtkvtg

class XYChart(object):

	def __init__(self, data_source, input_link=None, highlight_link=None):
		"""Parallel coordinates view constructor needs a valid DataSource plus
		and external annotation link (from the icicle view).
		"""

		self.ds = data_source

		self.input_link = None
		if input_link is not None:
			self.SetInputAnnotationLink(input_link)

		# Set up a 2D scene, add an XY chart to it
		# Relies on changes to VTK (drawimage branch) that allow DrawImage() scaling factor
		self.chartView = vtk.vtkContextView()
		self.chartView.GetRenderer().SetBackground(1.0, 1.0, 1.0)

		self.chart = vtkvtg.vtkMyChartXY()
		self.chartView.GetScene().AddItem(self.chart)

		self.axisView = vtk.vtkContextView()
		self.axisView.GetRenderer().SetBackground(1.0, 1.0, 1.0)

		self.ai = vtkvtg.vtkAxisImageItem()
		self.axisView.GetScene().AddItem(self.ai)

		self.highlight_link = None
		if highlight_link is not None:
			self.SetHighlightAnnotationLink(highlight_link)
			# Set up callback to listen for changes in IcicleView selections
			self.highlight_link.AddObserver("AnnotationChangedEvent", self.HighlightSelectionCallback)

			# Set up annotation link which will carry indices to parallel coordinates chart
			# for highlighting outside selections (e.g. back from image_flow)
			# This needs to carry indices, while image_flow link outputs pedigree ids
			# so conversion happens in HighlightSelectionCallback
			self.highlight_link_idxs = vtk.vtkAnnotationLink()
			self.highlight_link_idxs.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
			self.highlight_link_idxs.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
			self.chart.SetHighlightLink(self.highlight_link_idxs)

		# Create a annotation link to access selection in parallel coordinates view
		self.link = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		self.link.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		# The chart seems to force INDEX selection, so I'm using vtkConvertSelection below to get
		# out PedigreeIds...
		self.link.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
		# Connect the annotation link to the parallel coordinates representation
		self.chart.SetAnnotationLink(self.link)

		# Set up callback for ID -> Pedigree ID conversion & copy to output link
		# self.link.AddObserver("AnnotationChangedEvent", self.PCoordsSelectionCallback)

		# Set up output annotation link which will carry pedigree ids to image flow view
		# Type and field will be set during conversion to pedigree ids in PCoordsSelectionCallback
		# self.output_link = vtk.vtkAnnotationLink()

		# self.chartView.ResetCamera()
		# self.chartView.Render()

		# Want to keep track of whether the node coming in on the input_link
		# is new or not
		self.input_link_idx = 0

		# Flag for whether to color by 'category_ids'
		self.color_by_array = True
		self.color_array_name = 'category_ids'

		cl = []
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

# 		cl.append([float(cc)/255.0 for cc in [228, 26, 28]])  # Colorbrewer Set2 modY+4
# 		cl.append([float(cc)/255.0 for cc in [55, 126, 184]])
# 		cl.append([float(cc)/255.0 for cc in [77, 175, 74]])
# 		cl.append([float(cc)/255.0 for cc in [152, 78, 163]])
# 		cl.append([float(cc)/255.0 for cc in [255, 127, 0]])
# 		cl.append([float(cc)/255.0 for cc in [245, 193, 61]])
# 		cl.append([float(cc)/255.0 for cc in [166, 86, 40]])
# 		cl.append([float(cc)/255.0 for cc in [247, 129, 191]])
# 		cl.append([float(cc)/255.0 for cc in [153, 153, 153]])
# 		cl.append([float(cc)/255.0 for cc in [143, 0, 0]])
# 		cl.append([float(cc)/255.0 for cc in [22, 65, 110]])
# 		cl.append([float(cc)/255.0 for cc in [40, 115, 33]])
# 		cl.append([float(cc)/255.0 for cc in [61, 61, 61]])

		self.lut = vtk.vtkLookupTable()
		lutNum = len(cl)
		self.lut.SetNumberOfTableValues(lutNum)
		self.lut.Build()
		for ii,cc in enumerate(cl):
			self.lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
		self.lut.SetRange(0,len(cl)-1)
		self.lut.SetAlpha(0.8)

	def SetInputAnnotationLink(self, link):

		self.input_link = link

		# Set up callback to listen for changes in IcicleView selections
		self.input_link.AddObserver("AnnotationChangedEvent", self.InputSelectionCallback)

	def SetAnnotationLink(self, link):

		self.link = link
		self.link.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		self.link.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
		self.chart.SetAnnotationLink(self.link)

		# Set up callback for ID -> Pedigree ID conversion & copy to output link
		self.link.AddObserver("AnnotationChangedEvent", self.XYSelectionCallback)

	def SetHighlightAnnotationLink(self, link):

		self.highlight_link = link

		# Set up callback to listen for changes in IcicleView selections
		self.highlight_link.AddObserver("AnnotationChangedEvent", self.HighlightSelectionCallback)

		# Set up annotation link which will carry indices to parallel coordinates chart
		# for highlighting outside selections (e.g. back from image_flow)
		# This needs to carry indices, while image_flow link outputs pedigree ids
		# so conversion happens in HighlightSelectionCallback
		self.highlight_link_idxs = vtk.vtkAnnotationLink()
		self.highlight_link_idxs.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
		self.highlight_link_idxs.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
		self.chart.SetHighlightLink(self.highlight_link_idxs)

	def SetColorByArrayOn(self):

		self.color_by_array = True
		self.color_array_name = 'category_ids'

	def SetColorByArray(self, array_name):

		self.color_by_array = True
		if type(array_name).__name__ == 'str':
			self.color_array_name = array_name
		else:
			self.color_array_name = ''

	def SetColorByCategoryOff(self):

		self.color_by_array = False
		self.color_array_name = ''

	def XYSelectionCallback(self, caller, event):

		# Right now just taking in selections, so not converting output_link to pedigree_ids
		self.chartView.Render()

	def InputSelectionCallback(self, caller, event):
		"""This is the callback that tracks changes in the icicle view and
		sets the input table for the parallel coordinates chart accordingly.
		Note: This can only handle a single input nodeID for now..."""

		annSel = caller.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		if annSel.GetNumberOfNodes() > 0:
			idxVtk = annSel.GetNode(0).GetSelectionList()
			idxArr = VN.vtk_to_numpy(idxVtk)
			print "ice input to XY ", idxArr

			# Here is where I'm limiting the selection to _one_ nodeID for now...
			node_id = idxArr[0]
			self.table = self.ds.GetNodeOneScaleCoeffTable(node_id)
			id_list = self.ds.PointsInNet[node_id]	# Directly accessing member variable
			self.image_stack = self.ds.GetProjectedImages(id_list)
			self.axis_images = self.ds.GetNodeBasisImages(node_id)
			self.center_image = self.ds.GetNodeCenterImage(node_id)

			# Get the axis image XY indices in case resetting to those values
			# and the number of dimensions has changed, and xI or yI are over the limit
			xI = self.ai.GetXAxisIndex()
			yI = self.ai.GetYAxisIndex()
			if yI > xI:
				y_bigger = 1
			else:
				y_bigger = 0

			# Check whether they're out of range for the new data
			(dX,dY,dZ) = self.axis_images.GetDimensions()
			max_dim = dZ - 1
			if yI > max_dim:
				yI = max_dim
			if xI > max_dim:
				xI = max_dim
			if xI == yI:
				if y_bigger:
					xI = yI-1
				else:
					yI = xI-1
			if xI < 0 or yI < 0:
				xI = 0
				yI = 1

			self.chart.ClearPlots()
			self.ai.ClearAxisImages()

			line1 = vtkvtg.vtkMyPlotPoints()
			self.chart.AddPlot(line1)		# POINTS
			if (self.table.GetNumberOfColumns() > 2):
				line1.SetInput(self.table, 0, 1)
			else:
				line1.SetInput(self.table, 0, 0)
			line1.SetMarkerStyle(vtkvtg.vtkMyPlotPoints.CIRCLE)

			if self.color_by_array:
				line1.SetScalarVisibility(1)
				line1.SetLookupTable(self.lut)
				line1.SelectColorArray(self.color_array_name)
			else:
				line1.SetScalarVisibility(0)
				line1.SetColor(0, 0, 0, 255)

			# If this is the same icicle node as before, then reset to original XY indices
			# before view is updated
			if node_id == self.input_link_idx:
				self.chart.SetPlotColumnIndices(xI,yI)

			# Need to set the image stack for the plot which will get resliced
			self.chart.SetTooltipImageStack(self.image_stack)
			self.chart.SetTooltipShowImage(True)
			# self.chart.SetTooltipImageScalingFactor(2.0)
			self.chart.SetTooltipImageTargetSize(50)
			self.chart.Update()

			self.ai.SetAxisImagesHorizontal()
			self.ai.SetAxisImageStack(self.axis_images)
			self.ai.SetCenterImage(self.center_image)

			# If this is the same icicle node as before, then reset to original XY indices
			# before view is updated
			if node_id == self.input_link_idx:
				self.ai.SetAxisIndices(xI,yI)

			self.ai.Update()

			self.PedIdToIndexSelection()

			# self.chartView.ResetCamera()
			self.input_link_idx = node_id
			self.chartView.Render()
			self.axisView.Render()

		else:
			self.chart.ClearPlots()
			self.table = None
			# self.chart.Update()
			# self.chartView.ResetCamera()
			self.chartView.Render()

			self.ai.ClearAxisImages()
			self.axisView.Render()


	def PedIdToIndexSelection(self):

		# Convert pedigree id highlight_link to highlight_link_idxs selection
		# (Don't have to do this with output_link & link right now because this chart
		#  is sharing the internal index link with the pcoords_chart which is doing the calc.
		#  In principle could just share the same highlight_link_idxs, too and not have to do
		#  this...)
		pedIdSel = self.highlight_link.GetCurrentSelection()
		cs = vtk.vtkConvertSelection()
		idxSel = cs.ToIndexSelection(pedIdSel, self.table)
		if idxSel.GetNumberOfNodes() > 0:
			idxVtk = idxSel.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				print "XY highlight_link_idxs indices: ", VN.vtk_to_numpy(idxVtk)
			else:
				print "XY highlight_link_idxs NO TUPLES"
		else:
			print "XY highlight_link_idxs NO NODES"
		self.highlight_link_idxs.SetCurrentSelection(idxSel)

	def HighlightSelectionCallback(self, caller, event):
		# Need to convert pedigree IDs that we get back from image_flow into indices
		annSel = caller.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		if annSel.GetNumberOfNodes() > 0:
			idxVtk = annSel.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				cs = vtk.vtkConvertSelection()
				idxSelection = cs.ToIndexSelection(annSel, self.table)
				# Copy converted selection to output annotation link (fires AnnotationChangedEvent)
				self.highlight_link_idxs.SetCurrentSelection(idxSelection)
			else:
				empty_arr = N.array([],dtype='int64')
				empty_vtk = VN.numpy_to_vtkIdTypeArray(empty_arr, deep=True)
				self.highlight_link_idxs.GetCurrentSelection().GetNode(0).SetSelectionList(empty_vtk)

		self.chartView.Render()

# 	def GetOutputAnnotationLink(self):
# 		# The point here would be, like with the pcoords chart, to output pedigree_ids
# 		# so that other views could collect the correct original data based on selections here
# 		return self.output_link

	def GetChartView(self):
		return self.chartView

	def GetAxisView(self):
		return self.axisView

	def GetAxisImageItem(self):
		return self.ai

	def GetChartXY(self):
		return self.chart

	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.chartView.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)




