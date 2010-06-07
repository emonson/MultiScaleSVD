import vtk
import vtk.util.numpy_support as VN
import numpy as N
import sys
# sys.path.append("/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
# import libvtkvtgChartsPython as vtgCh
import vtkvtg

class XYChart(object):

	def __init__(self, data_source, input_link=None):
		"""Parallel coordinates view constructor needs a valid DataSource plus
		and external annotation link (from the icicle view).
		"""
		
		self.ds = data_source
		
		self.input_link = None
		if input_link is not None:
			self.SetInputAnnotationLink(input_link)
			
		# Set up a 2D scene, add an XY chart to it
		self.view = vtkvtg.vtkMyContextView()
		self.view.GetRenderWindow().SetSize(600,300)
		
		self.chart = vtkvtg.vtkMyChartXY()
		
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
		
	def InputSelectionCallback(self, caller, event):
		"""This is the callback that tracks changes in the icicle view and
		sets the input table for the parallel coordinates chart accordingly.
		Note: This can only handle a single input nodeID for now..."""
		
		annSel = caller.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		if annSel.GetNumberOfNodes() > 0:
			idxVtk = annSel.GetNode(0).GetSelectionList()
			idxArr = VN.vtk_to_numpy(idxVtk)
			print "ice input to PC ", idxArr
			
			# Here is where I'm limiting the selection to _one_ nodeID for now...
			node_id = idxArr[0]
			self.table = self.ds.GetNodeOneScaleCoeffTable(node_id)
			id_list = self.ds.PIN[node_id]
			image_stack = self.ds.GetProjectedImages(id_list)

			self.chart.ClearPlots()
			line1 = self.chart.AddPlot(1)		# POINTS
			line1.SetInput(self.table, 0, 1)
			line1.SetMarkerStyle(2)
			line1.SetColor(255, 0, 0, 255)

			# Need to set the image stack for the plot which will get resliced 
			self.chart.GetPlot(0).SetImageStack(image_stack)
			self.chart.SetTooltipShowImage(True)
			
			# self.view.ResetCamera()
			# self.view.Render()

		else:
			self.table = None
			self.chart.ClearPlots()
			# self.chart.Update()
			# self.view.ResetCamera()
			# self.view.Render()
			
			
# 	def GetOutputAnnotationLink(self):
# 		return self.output_link
	
	def GetView(self):
		return self.view
		
	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.view.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)




