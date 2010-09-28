import vtk
import vtk.util.numpy_support as VN
import numpy as N
import sys
# sys.path.append("/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
# import libvtkvtgChartsPython as vtgCh
import vtkvtg

class PCoordsChart(object):

	def __init__(self, data_source, input_link=None, highlight_link=None):
		"""Parallel coordinates view constructor needs a valid DataSource plus
		and external annotation link (from the icicle view).
		"""
		
		self.ds = data_source
		
		self.input_link = None
		if input_link is not None:
			self.SetInputAnnotationLink(input_link)
			
		# Set up a 2D scene, add an XY chart to it
		self.view = vtk.vtkContextView()
		self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
		self.view.GetRenderWindow().SetSize(600,300)
		
		self.chart = vtkvtg.vtkMyChartParallelCoordinates()
		
		self.highlight_link = None
		if highlight_link is not None:
			self.SetHighlightAnnotationLink(highlight_link)
		
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
		self.link.AddObserver("AnnotationChangedEvent", self.PCoordsSelectionCallback)
		
		# Set up output annotation link which will carry pedigree ids to image flow view
		# Type and field will be set during conversion to pedigree ids in PCoordsSelectionCallback
		self.output_link = vtk.vtkAnnotationLink()
		
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
		self.link.AddObserver("AnnotationChangedEvent", self.PCoordsSelectionCallback)
	
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

	def SetScaleDimLimit(self, dim_limit):
		self.dim_limit = dim_limit
		self.ReloadData()

	def ReloadData(self):
		self.InputSelectionCallback(self.input_link, None)

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
			
			# Manually clear out selections if changing data
			# self.chart.GetPlot(0).ResetSelectionRange()

			# Here is where I'm limiting the selection to _one_ nodeID for now...
			node_id = idxArr[0]
			
			# NOTE: Hard coding for testing!!!
			# numPerSet = self.ds.ManifoldDim			# Directly accessing member variable...	
			# numPerSet = 3
			
			# NOTE: Using only first few dimensions for testing!!!
			self.table, self.scale_dims = self.ds.GetNodeAllScaleCoeffTable(node_id)
			# self.table = self.ds.GetNodeAllScaleDDimCoeffTable(node_id,numPerSet)
			begin_idxs = N.array(self.scale_dims).cumsum()

			currentScale = self.ds.Scales[node_id]	# Directly accessing member variable...
			
			# Since can't really clear out and reload new plot on PCoords chart, 
			# and much is based on chart visible columns, trying to clear out visibility
			# and reload new before change data. Column visibility is just based on names
			# in a vtkStringArray
			self.chart.SetAllColumnsInvisible()
			
			# By default only 10 axes (columns). Set these before try to change data
			# since plot update depends on column visibility numbers
			for ii in range(self.table.GetNumberOfColumns()):
				col_name = self.table.GetColumnName(ii)
				if not col_name.lower().endswith('_ids'):
					self.chart.SetColumnVisibility(col_name, True)
			
			self.chart.GetPlot(0).SetInput(self.table)
			self.chart.GetPlot(0).Modified()
			self.chart.SetDrawSets(True)
			# self.chart.SetNumPerSet(numPerSet)
			self.chart.SetNumberOfScales(len(self.scale_dims))
			for ii in range(len(self.scale_dims)):
				self.chart.SetScaleDim(ii, self.scale_dims[ii])
			self.chart.SetCurrentScale(currentScale)
			self.chart.Update()
			
			# Get extrema from _all_ wavelet coefficients instead of just used columns
			exRange = self.ds.GetCoeffRange()
			# Move max/min in a bit on plot
			chMax = exRange[1]*1.05	
			chMin = exRange[0]*1.05
			tickPos = VN.numpy_to_vtk(N.array([chMin, 0.0, chMax]), deep=True)
			tickNoPos = VN.numpy_to_vtk(N.array([]), deep=True)
			
			tickLabels = vtk.vtkStringArray()
			tickLabels.SetName('tick_labels')
			tickLabels.SetNumberOfComponents(1)
			tickLabels.InsertNextValue(str(int(N.round(chMin))))
			tickLabels.InsertNextValue(str(int(0)))
			tickLabels.InsertNextValue(str(int(N.round(chMax))))
			
			noLabels = vtk.vtkStringArray()
			noLabels.SetName('no_labels')
			noLabels.SetNumberOfComponents(1)
			
			for ii in range(self.table.GetNumberOfColumns()):
				col_name = self.table.GetColumnName(ii)
				if not col_name.lower().endswith('_ids'):
					self.chart.GetAxis(ii).SetBehavior(1)	# fixed
					# Set same range for all axes
					self.chart.GetAxis(ii).SetMaximum(chMax)
					self.chart.GetAxis(ii).SetMinimum(chMin)
					# Only put ticks on first axis, plus last in any set
					if (ii==0): # or (ii not in begin_idxs):
						self.chart.GetAxis(ii).SetTickPositions(tickPos)
						self.chart.GetAxis(ii).SetTickLabels(tickLabels)
					else:
						self.chart.GetAxis(ii).SetTickPositions(tickNoPos)
						self.chart.GetAxis(ii).SetTickLabels(noLabels)
					# Only put range labels on first set of axes
					if ii != 0:
						self.chart.GetAxis(ii).SetLabelsVisible(False)
					# Only put title at bottom for first in set, and remove the .0 from that
					if not col_name.endswith('.0'):
						self.chart.GetAxis(ii).SetTitle('')
					else:
						self.chart.GetAxis(ii).SetTitle(col_name[:-2])
			
			self.PedIdToIndexSelection()

			self.chart.Modified()
			self.view.Render()

		else:
			# Need to clear out annotation link so downstream views will be cleared
			# before table is destroyed
			empty_vtk = VN.numpy_to_vtkIdTypeArray(N.array([],dtype='int64'), deep=True)
			self.link.GetCurrentSelection().GetNode(0).SetSelectionList(empty_vtk)
 			empty_vtk2 = VN.numpy_to_vtkIdTypeArray(N.array([],dtype='int64'), deep=True)
 			self.highlight_link.GetCurrentSelection().GetNode(0).SetSelectionList(empty_vtk2)
			
			# And make sure the output_link knows the selection has changed
			# NOTE: Commented out for now so we can preserve selections across icicle_view node changes
			# self.link.InvokeEvent("AnnotationChangedEvent")
 			# self.highlight_link.InvokeEvent("AnnotationChangedEvent")

			# If there is no SelectionNode in IcicleView selection
			print "Ice cleared PC called"
			# For right now the only way I can get it to clear out the data is
			# to set the table=None. Should try creating a table with all the right
			# columns but no data in each column...
			self.table = None
			self.chart.GetPlot(0).SetInput(self.table)
			# self.chart.Update()
			# self.view.ResetCamera()
			# self.view.Render()
			
	def PedIdToIndexSelection(self):
		
		# Get pedigree ids from current output_link selection
		pedIdSel = self.output_link.GetCurrentSelection()
		cs = vtk.vtkConvertSelection()
		idxSel = cs.ToIndexSelection(pedIdSel, self.table)
		
		# Trying to disable event invocation while chaging this selection
		self.link.RemoveObservers("AnnotationChangedEvent")
		self.link.SetCurrentSelection(idxSel)
		self.link.AddObserver("AnnotationChangedEvent", self.PCoordsSelectionCallback)
		
		# Now do the same for the highlight_link_idxs selection
		pedIdSel = self.highlight_link.GetCurrentSelection()		
		cs = vtk.vtkConvertSelection()
		idxSel = cs.ToIndexSelection(pedIdSel, self.table)
		if idxSel.GetNumberOfNodes() > 0:
			idxVtk = idxSel.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				print "PC highlight_link_idxs indices: ", VN.vtk_to_numpy(idxVtk)
			else:
				print "PC highlight_link_idxs NO TUPLES"
		else:
			print "PC highlight_link_idxs NO NODES"
		self.highlight_link_idxs.SetCurrentSelection(idxSel)
		print "_>_>_>_> set PC current selection"

	def PCoordsSelectionCallback(self, caller, event):
		# Defined for testing ID picking
		annSel = caller.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		cs = vtk.vtkConvertSelection()
		pedIdSelection = cs.ToPedigreeIdSelection(annSel, self.table)
		
		# Copy converted selection to output annotation link (fires AnnotationChangedEvent)
		self.output_link.SetCurrentSelection(pedIdSelection)
		
		# Test conversion by printing
		if pedIdSelection.GetNumberOfNodes() > 0:
			idxVtk = pedIdSelection.GetNode(0).GetSelectionList()
			if idxVtk.GetNumberOfTuples() > 0:
				print "PC ped ids: ", VN.vtk_to_numpy(idxVtk)
		else:
			print "PC empty selection"
		
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
		
	def GetOutputAnnotationLink(self):
		# This one contains pedigree_ids so that things like the image_flow can collect
		# the correct original data based on pcoords selections
		return self.output_link
	
	def GetAnnotationLink(self):
		# This one contain indices, which internally is how the pcoords chart
		# handles selections
		return self.link

	def GetView(self):
		return self.view
		
	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.view.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)

if __name__ == "__main__":

	# from tkFileDialog import askopenfilename
	from data_source import DataSource
	
	# data_file = askopenfilename()
	data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100909.mat'

	# DataSource loads .mat file and can generate data from it for other views
	ds = DataSource(data_file)
		
	# Set up an annotation link as if selections were coming from another class
	ice_output_link = vtk.vtkAnnotationLink()
	# ice_output_link.GetCurrentSelection().GetNode(0).SetFieldType(3)		# Vertex
	# ice_output_link.GetCurrentSelection().GetNode(0).SetContentType(2)	# Pedigree Ids
	# Above are the correct types, but we need a blank selection to begin
	ice_output_link.GetCurrentSelection().RemoveAllNodes()

	pc_class = PCoordsChart(ds)
	pc_class.SetInputAnnotationLink(ice_output_link)
	pc_class.GetView().GetRenderWindow().SetSize(600,300)
		
	# Set up an annotation link as if selections were coming from another class
	dummy_link2 = vtk.vtkAnnotationLink()
	dummy_link2.GetCurrentSelection().GetNode(0).SetFieldType(1)     # Point
	dummy_link2.GetCurrentSelection().GetNode(0).SetContentType(4)   # 2 = PedigreeIds, 4 = Indices
	pc_class.SetHighlightAnnotationLink(dummy_link2)
	
	# Fill selection link with dummy IDs
	id_array = N.array([102],dtype='int64')
	id_list = VN.numpy_to_vtkIdTypeArray(id_array)
	node = vtk.vtkSelectionNode()
	node.SetFieldType(3)		# Vertex
	node.SetContentType(2)	# Pedigree Ids
	node.SetSelectionList(id_list)
	ice_output_link.GetCurrentSelection().AddNode(node)
	ice_output_link.InvokeEvent("AnnotationChangedEvent")
	
	# Only need to Start() interactor for one view
	pc_class.GetView().GetRenderWindow().GetInteractor().Start()



