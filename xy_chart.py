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
		self.view = vtk.vtkContextView()
		self.view.GetRenderWindow().SetSize(600,300)
		
		self.chart = vtkvtg.vtkMyChartXY()
		
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
		
		self.view.GetScene().AddItem(self.chart)
		# self.view.ResetCamera()
		# self.view.Render()


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

	def XYSelectionCallback(self, caller, event):
	
		# Right now just taking in selections, so not converting output_link to pedigree_ids
		print "XY ANNOTATION LINK"
		print caller
		self.view.Render()
	
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
			id_list = self.ds.PIN[node_id]
			self.image_stack = self.ds.GetProjectedImages(id_list)

			self.chart.ClearPlots()
			line1 = self.chart.AddPlot(1)		# POINTS
			line1.SetInput(self.table, 0, 1)
			line1.SetMarkerStyle(2)
			line1.SetColor(0, 0, 0, 255)

			# Need to set the image stack for the plot which will get resliced 
			self.chart.GetPlot(0).SetImageStack(self.image_stack)
			self.chart.SetTooltipShowImage(True)
			self.chart.SetTooltipImageScalingFactor(2.0)

			# self.view.ResetCamera()
			self.view.Render()

		else:
			self.chart.ClearPlots()
			self.table = None
			# self.chart.Update()
			# self.view.ResetCamera()
			self.view.Render()
			
			
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
			
		self.view.Render()
		
# 	def GetOutputAnnotationLink(self):
# 		# The point here would be, like with the pcoords chart, to output pedigree_ids
# 		# so that other views could collect the correct original data based on selections here
# 		return self.output_link
	
	def GetView(self):
		return self.view
		
	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.view.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)



