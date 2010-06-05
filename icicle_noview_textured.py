import vtk
import vtk.util.numpy_support as VN
import numpy as N

class NavMenu(object):
	def __init__(self):
		self.poly_list = []
		self.actor_list = []
		
		# Create shapes and actors
		boxUp = vtk.vtkRegularPolygonSource()
		boxUp.GeneratePolygonOn()
		boxUp.GeneratePolylineOff()
		boxUp.SetNumberOfSides(3)
		boxUp.SetCenter(0,14,0)
		boxUp.SetRadius(11)
		self.poly_list.append(boxUp)
		
		boxL = vtk.vtkRegularPolygonSource()
		boxL.GeneratePolygonOn()
		boxL.GeneratePolylineOff()
		boxL.SetNumberOfSides(4)
		boxL.SetCenter(-15,-2,0)
		boxL.SetRadius(9)
		self.poly_list.append(boxL)
		
		boxR = vtk.vtkRegularPolygonSource()
		boxR.GeneratePolygonOn()
		boxR.GeneratePolylineOff()
		boxR.SetNumberOfSides(4)
		boxR.SetCenter(15,-2,0)
		boxR.SetRadius(9)
		self.poly_list.append(boxR)
		
		boxLDown = vtk.vtkRegularPolygonSource()
		boxLDown.GeneratePolygonOn()
		boxLDown.GeneratePolylineOff()
		boxLDown.SetNumberOfSides(6)
		boxLDown.SetCenter(-12,-22,0)
		boxLDown.SetRadius(9)
		self.poly_list.append(boxLDown)
		
		boxRDown = vtk.vtkRegularPolygonSource()
		boxRDown.GeneratePolygonOn()
		boxRDown.GeneratePolylineOff()
		boxRDown.SetNumberOfSides(6)
		boxRDown.SetCenter(12,-22,0)
		boxRDown.SetRadius(9)
		self.poly_list.append(boxRDown)
		
		for ii, poly in enumerate(self.poly_list):
			map = vtk.vtkPolyDataMapper2D()
			map.SetInputConnection(poly.GetOutputPort(0))
			act = vtk.vtkActor2D()
			act.SetMapper(map)
			act.SetPickable(True)
			act.GetProperty().SetColor(0.5, 0.45, 0.35)
			act.GetProperty().SetLineWidth(0.0)
			act.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
			act.GetPositionCoordinate().SetValue(0.075 ,0.1)
			act.GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
			act.GetPosition2Coordinate().SetValue(0.5, 0.5)
			self.actor_list.append(act)
			
	def GetPolySourceList(self):
		return self.poly_list
	
	def GetActorList(self):
		return self.actor_list
		

class IcicleNoView(object):

	def __init__(self, data_source):
		"""Icicle view constructor needs a valid DataSource and will pull data
		from it immediately."""
		
		self.ds = data_source
		
		self.SHRINK = 0.1
		self.THICK = 1.0
		
		# Create the RenderWindow, Renderer and both Actors
		self.renderer = vtk.vtkRenderer()
		self.renWin = vtk.vtkRenderWindow()
		self.renWin.AddRenderer(self.renderer)
		self.interactor = vtk.vtkRenderWindowInteractor()
		self.interactor.SetRenderWindow(self.renWin)
		self.istyle = vtk.vtkInteractorStyleImage()
		self.interactor.SetInteractorStyle(self.istyle)
		self.renderer.GetActiveCamera().ParallelProjectionOn()

		# Create callbacks for mouse events
		self.mouseActions = {}
		self.mouseActions["LeftButtonDown"] = 0
		self.mouseActions["Picking"] = 0
		
		# self.istyle.AddObserver("MouseMoveEvent", self.MouseMoveCallback)
		self.istyle.AddObserver("LeftButtonPressEvent", self.LeftButtonPressCallback)
		self.istyle.AddObserver("LeftButtonReleaseEvent", self.LeftButtonReleaseCallback)

		
		# Apply a theme to the views
		self.theme = vtk.vtkViewTheme.CreateMellowTheme()
		self.theme.SetPointColor(0,0,0)
		self.theme.SetPointOpacity(0.1)
		c = N.array([255,204,0])/255.0		# nice yellow
		self.theme.SetSelectedPointColor(c[0],c[1],c[2])
		self.theme.SetBackgroundColor(0.1, 0.1, 0.06)
		self.theme.SetBackgroundColor2(0.25, 0.25, 0.2)
		
		self.renderer.SetBackground(self.theme.GetBackgroundColor())
		self.renderer.SetBackground2(self.theme.GetBackgroundColor2())
		self.renderer.SetGradientBackground(True)
		
		# Grab annotation link to monitor selection changes
		# self.output_link = rep.GetAnnotationLink()
		
		self.output_link = vtk.vtkAnnotationLink()
		# If you don't set the FieldType explicitly it ends up as UNKNOWN (as of 21 Feb 2010)
		# See vtkSelectionNode doc for field and content type enum values
		# enum		SelectionContent { 
		#	SELECTIONS, GLOBALIDS, PEDIGREEIDS, VALUES, 
		#	INDICES, FRUSTUM, LOCATIONS, THRESHOLDS, 
		#	BLOCKS 
		# }
		# enum		SelectionField { 
		#	CELL, POINT, FIELD, VERTEX, 
		#	EDGE, ROW 
		# }
		
		# Note: This may be the defaults, anyway, for the IcicleView annotation link
		# self.output_link.GetCurrentSelection().GetNode(0).SetFieldType(3)		# Vertex
		# self.output_link.GetCurrentSelection().GetNode(0).SetContentType(2)	# Pedigree Ids
		
		# Above are the correct types, but we need a blank selection to begin
		self.output_link.GetCurrentSelection().RemoveAllNodes()
		
		# TEST: Set up callback to test icicle view selection IDs
		self.output_link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)
		
		# Connect the annotation link to the icicle representation
		# rep = self.view.GetRepresentation(0)
		
		# Create Navigation Menu
		self.menu = NavMenu()
		menu_actor_list = self.menu.GetActorList()
		for m_actor in menu_actor_list:
			self.renderer.AddActor(m_actor)
		
		# Create a list of functions to be called depending on which nav button is clicked
		# Note: this is a bit fragile since it needs to correspond in number and routine
		# to the same order as created in nav menu class...
		self.menu_functions = []
		self.menu_functions.append(self.OnNavUp)
		self.menu_functions.append(self.OnNavL)
		self.menu_functions.append(self.OnNavR)
		self.menu_functions.append(self.OnNavLDown)
		self.menu_functions.append(self.OnNavRDown)
		
		# Load Data
		self.LoadData()

		self.renderer.ResetCamera()
		# self.renWin.Render()

	def LoadData(self):
		
		# Remove all old actors from renderer
		self.renderer.RemoveAllViewProps()
		# Put back nav menu
		menu_actor_list = self.menu.GetActorList()
		for m_actor in menu_actor_list:
			self.renderer.AddActor(m_actor)
		
		tree = self.ds.GetTree()
		
		# Parallel pipeline with no shrinkage to get TCoords
		self.TreeLevels = vtk.vtkTreeLevelsFilter()
		self.TreeLevels.SetInput(tree)
		
		VertexDegree = vtk.vtkVertexDegree()
		VertexDegree.SetInputConnection(self.TreeLevels.GetOutputPort(0))
		
		TreeAggregation = vtk.vtkTreeFieldAggregator()
		TreeAggregation.LeafVertexUnitSizeOff()
		TreeAggregation.SetField('size')
		TreeAggregation.SetInputConnection(VertexDegree.GetOutputPort(0))
		
		# Layout without shrinkage for generating texture coordinates
		strategy = vtk.vtkStackedTreeLayoutStrategy()
		strategy.UseRectangularCoordinatesOn()
		strategy.SetRootStartAngle(0.0)
		strategy.SetRootEndAngle(15.0)
		strategy.SetRingThickness(self.THICK)	# layer thickness
		strategy.ReverseOn()
		strategy.SetShrinkPercentage(0.0)
		layout = vtk.vtkAreaLayout()
		layout.SetLayoutStrategy(strategy)
		layout.SetInputConnection(TreeAggregation.GetOutputPort(0))
		layout.SetAreaArrayName("area")
		layout.SetSizeArrayName("num_in_vertex")
		areapoly = vtk.vtkTreeMapToPolyData()
		areapoly.SetInputConnection(layout.GetOutputPort(0))
		areapoly.SetAddNormals(0)
		areapoly.SetInputArrayToProcess( 0, 0, 0, 4, "area")  # 4 = vtkDataObject::FIELD_ASSOCIATION_VERTICES
		texPlane = vtk.vtkTextureMapToPlane()
		texPlane.SetInputConnection(areapoly.GetOutputPort(0))
		texPlane.AutomaticPlaneGenerationOn()
		texPlane.Update()
		
		# Layout with shrinkage for generating geometry
		strategy0 = vtk.vtkStackedTreeLayoutStrategy()
		strategy0.UseRectangularCoordinatesOn()
		strategy0.SetRootStartAngle(0.0)
		strategy0.SetRootEndAngle(15.0)
		strategy0.SetRingThickness(self.THICK)	# layer thickness
		strategy0.ReverseOn()
		strategy0.SetShrinkPercentage(self.SHRINK)
		layout0 = vtk.vtkAreaLayout()
		layout0.SetLayoutStrategy(strategy0)
		layout0.SetInputConnection(TreeAggregation.GetOutputPort(0))
		layout0.SetAreaArrayName("area")
		layout0.SetSizeArrayName("num_in_vertex")
		areapoly0 = vtk.vtkTreeMapToPolyData()
		areapoly0.SetAddNormals(0)
		areapoly0.SetInputConnection(layout0.GetOutputPort(0))
		areapoly0.SetInputArrayToProcess( 0, 0, 0, 4, "area")  # 4 = vtkDataObject::FIELD_ASSOCIATION_VERTICES
		areapoly0.Update()
		
		# Copy over texture coordinates
		def transferTCoords():
			input = paf.GetInputDataObject(0,0)
			refin = paf.GetInputList().GetItem(0)
			output = paf.GetPolyDataOutput()
			
			TCorig = refin.GetPointData().GetTCoords()
			
			TC = vtk.vtkFloatArray()
			TC.SetNumberOfComponents(TCorig.GetNumberOfComponents())
			TC.SetNumberOfTuples(TCorig.GetNumberOfTuples())
			TC.SetName('Texture Coordinates')
			for ii in range(TCorig.GetNumberOfTuples()):
				ff = TCorig.GetTuple2(ii)
				TC.SetTuple2(ii,ff[0],ff[1])
			
			output.GetPointData().AddArray(TC)
			output.GetPointData().SetActiveTCoords('Texture Coordinates')
			
		paf = vtk.vtkProgrammableAttributeDataFilter()
		paf.SetInput(areapoly0.GetOutput())
		paf.AddInput(texPlane.GetOutput())
		paf.SetExecuteMethod(transferTCoords)
		
		# Need to find proper ordering of wavelet coeffs based on icicle layout
		# tree.GetVertexData().GetArray('area') is 4-component (Xmin,Xmax,Ymin,Ymax)
		print 'Reordering wavelet coeffs'
		out_polys = areapoly.GetOutputDataObject(0)
		isleaf = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('leaf'))
		poly_bounds = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('area'))
		vertex_ids = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('vertex_ids'))
		
		LeafIds = vertex_ids[isleaf>0]
		LeafXmins = poly_bounds[isleaf>0,0]
		XOrderedLeafIds = LeafIds[LeafXmins.argsort()]
					
		# And then grab the Wavelet Coefficients matrix sorted according to this
		self.WCimageData = self.ds.GetWaveletCoeffImage(XOrderedLeafIds)
		
		WCrange = N.array(self.WCimageData.GetPointData().GetScalars().GetRange())
		WCext = abs(WCrange.min()) if (abs(WCrange.min()) > abs(WCrange.max())) else abs(WCrange.max())
		# print WCext
		
		# Create blue to white to red lookup table
		self.lut = vtk.vtkLookupTable()
		lutNum = 256
		self.lut.SetNumberOfTableValues(lutNum)
		self.lut.Build()
		ctf = vtk.vtkColorTransferFunction()
		ctf.SetColorSpaceToDiverging()
		ctf.AddRGBPoint(0.0, 0, 0, 1.0)
		ctf.AddRGBPoint(1.0, 1.0, 0, 0)
		for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
			cc = ctf.GetColor(ss)
			self.lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
		self.lut.SetRange(-WCext,WCext)
		
		# Set up texture with lookup table for matrix polys
		tex = vtk.vtkTexture()
		tex.SetInput(self.WCimageData)
		tex.SetLookupTable(self.lut)
		
		
		# Separate mapper and actor for textured polys
		map2 = vtk.vtkPolyDataMapper()
		map2.SetInputConnection(paf.GetOutputPort(0))
		map2.ScalarVisibilityOff()
		act2 = vtk.vtkActor()
		act2.SetMapper(map2)
		act2.SetTexture(tex)
		act2.GetProperty().SetColor(1,1,1)
		act2.SetPickable(0)
		act2.SetPosition(0,0,0.1)
		
		# Add textured polys to the view
		self.renderer.AddActor(act2)

				
		# Layout with shrinkage for generating outline geometry for showing selections
		self.applycolors1 = vtk.vtkApplyColors()
		self.applycolors1.SetInputConnection(0,layout0.GetOutputPort(0))
		self.applycolors1.AddInputConnection(1,self.output_link.GetOutputPort(0))
		self.applycolors1.SetDefaultPointColor(self.theme.GetPointColor())
		self.applycolors1.SetDefaultPointOpacity(self.theme.GetPointOpacity())
		self.applycolors1.SetSelectedPointColor(self.theme.GetSelectedPointColor())
		self.applycolors1.SetSelectedPointOpacity(self.theme.GetSelectedPointOpacity())

		self.areapoly1 = vtk.vtkTreeMapToPolyData()
		self.areapoly1.SetInputConnection(self.applycolors1.GetOutputPort(0))
		self.areapoly1.SetAddNormals(0)
		self.areapoly1.SetInputArrayToProcess( 0, 0, 0, 4, "area")  # 4 = vtkDataObject::FIELD_ASSOCIATION_VERTICES
		
		# Separate mapper and actor for icicle polys outlines (pickable)
		map = vtk.vtkPolyDataMapper()
		map.SetInputConnection(self.areapoly1.GetOutputPort(0))
		map.SetScalarModeToUseCellFieldData()
		map.SelectColorArray("vtkApplyColors color")
		map.SetScalarVisibility(True)
		act = vtk.vtkActor()
		act.SetMapper(map)
		act.GetProperty().SetColor(1,1,1)
		act.SetPickable(True)
		act.SetPosition(0,0,0)
		act.GetProperty().SetRepresentationToWireframe()
		act.GetProperty().SetLineWidth(3.0)
		
		self.icicle_actor = act
		
		# Add actor for selection highlight outlines
		self.renderer.AddActor(act)
		
		# Clear out selections on data change
		self.output_link.GetCurrentSelection().RemoveAllNodes()
		self.output_link.InvokeEvent("AnnotationChangedEvent")


	# Navigation "buttons" callbacks
	# NOTE: Directly accessing data source member variables here...
	def OnNavUp(self, prev_ped_id):
		new_ped_id = self.ds.CP[prev_ped_id]
		if new_ped_id >= 0:
			return new_ped_id
		else:
			return prev_ped_id
		
	def OnNavL(self, prev_ped_id):
		out_polys = self.areapoly1.GetOutputDataObject(0)
		poly_bounds = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('area'))
		vertex_ids = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('vertex_ids'))
		
		same_scale_ids = N.nonzero(self.ds.Scales==self.ds.Scales[prev_ped_id])[0]
		same_scale_xmins = poly_bounds[self.ds.Scales==self.ds.Scales[prev_ped_id], 0]
		ordered_scale_ids = same_scale_ids[same_scale_xmins.argsort()]
		prev_idx = N.nonzero(ordered_scale_ids==prev_ped_id)[0]
		
		if prev_idx > 0:
			return ordered_scale_ids[prev_idx - 1]
		else:
			return ordered_scale_ids[prev_idx]
		
	def OnNavR(self, prev_ped_id):
		out_polys = self.areapoly1.GetOutputDataObject(0)
		poly_bounds = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('area'))
		vertex_ids = VN.vtk_to_numpy(out_polys.GetCellData().GetArray('vertex_ids'))
		
		same_scale_ids = N.nonzero(self.ds.Scales==self.ds.Scales[prev_ped_id])[0]
		same_scale_xmins = poly_bounds[self.ds.Scales==self.ds.Scales[prev_ped_id], 0]
		ordered_scale_ids = same_scale_ids[same_scale_xmins.argsort()]
		prev_idx = N.nonzero(ordered_scale_ids==prev_ped_id)[0]
		
		if prev_idx < (ordered_scale_ids.size - 1):
			return ordered_scale_ids[prev_idx + 1]
		else:
			return ordered_scale_ids[prev_idx]
		
		
	def OnNavLDown(self, prev_ped_id):
		new_ped_id_arr = N.nonzero(self.ds.CP==prev_ped_id)[0]
		if new_ped_id_arr.size == 2:
			return new_ped_id_arr[0]
		elif new_ped_id_arr.size == 1:
			print "!!Only one child!!"
			return new_ped_id_arr[0]
		else:
			return prev_ped_id
		
	def OnNavRDown(self, prev_ped_id):
		new_ped_id_arr = N.nonzero(self.ds.CP==prev_ped_id)[0]
		if new_ped_id_arr.size == 2:
			return new_ped_id_arr[1]
		elif new_ped_id_arr.size == 1:
			print "!!Only one child!!"
			return new_ped_id_arr[0]
		else:
			return prev_ped_id
		
	# Mouse button callbacks
	def LeftButtonPressCallback(self, obj, event):
		self.mouseActions["LeftButtonDown"] = 1
		self.mouseActions["Picking"] = 1

	def LeftButtonReleaseCallback(self, obj, event):
		self.mouseActions["LeftButtonDown"] = 0
		if self.mouseActions["Picking"] == 1:
			self.selectNode()
		self.mouseActions["Picking"] = 0
	
	def MouseMoveCallback(self, obj, event):
		# Dummy routine for now...
		(lastX, lastY) = self.interactor.GetLastEventPosition()
		(mouseX, mouseY) = self.interactor.GetEventPosition()
		if self.mouseActions["LeftButtonDown"] == 1:
			self.interactorStyle.OnMouseMove()			
		else:
			self.interactorStyle.OnMouseMove()
			
	# Callback routine for any "Pick" event
	def selectNode(self):
		# (x0,y0) = self.interactor.GetLastEventPosition()
		(x,y) = self.interactor.GetEventPosition()
		cellPicker = vtk.vtkCellPicker()
		someCellPicked = cellPicker.Pick(x,y,0,self.renderer)
		pickedCellId = cellPicker.GetCellId()
		propPicker = vtk.vtkPropPicker()
		somePropPicked = propPicker.PickProp(x,y,self.renderer)
		pickedProp = propPicker.GetViewProp()
		
		if somePropPicked and (pickedProp != self.icicle_actor):
			# First, test whether there is a current selection
			if self.output_link.GetCurrentSelection().GetNumberOfNodes() > 0:
				# If so, get the current pedigree_id
				prev_ped_vtk = self.output_link.GetCurrentSelection().GetNode(0).GetSelectionList()
				# NOTE: Counting on a single ID
				prev_ped_id = prev_ped_vtk.GetValue(0)
				# Now figure out which menu item was picked
				for ii, m_actor in enumerate(self.menu.GetActorList()):
					if pickedProp == m_actor:
						# print "Menu Actor picked, index = ", ii
						# Call the proper nav routine and get back the new ped_id
						new_ped_id = self.menu_functions[ii](prev_ped_id)
						new_ped_n = N.array([new_ped_id], dtype='int64')
						new_ped_vtk = VN.numpy_to_vtkIdTypeArray(new_ped_n, deep=True)
						self.output_link.GetCurrentSelection().GetNode(0).SetSelectionList(new_ped_vtk)
						self.output_link.InvokeEvent("AnnotationChangedEvent")
						self.applycolors1.Update()
						self.renWin.Render()					
					
		# Cell picker doesn't work with Actor2D, so nav menu won't report any cells
		if someCellPicked:
			print "Icicle picked cell index: ", pickedCellId
			# Assuming for now that this is a cell "Index", so getting pedigree ID
			sel = vtk.vtkSelection()
			node = vtk.vtkSelectionNode()
			node.SetFieldType(0)		# Cell
			node.SetContentType(4)		# Indices
			id_array = N.array([pickedCellId], dtype='int64')
			id_vtk = VN.numpy_to_vtkIdTypeArray(id_array, deep=True)
			node.SetSelectionList(id_vtk)
			sel.AddNode(node)
			# Note: When selection is cleared, the current selection does NOT contain any nodes
			self.areapoly1.Update()
			cs = vtk.vtkConvertSelection()
			pedIdSelection = cs.ToPedigreeIdSelection(sel, self.areapoly1.GetOutput())
			pedIdSelection.GetNode(0).SetFieldType(3)	# convert to vertext selection
			
			# Copy converted selection to output annotation link (fires AnnotationChangedEvent)
			self.output_link.SetCurrentSelection(pedIdSelection)
			self.applycolors1.Update()
			self.renWin.Render()
		
		if not someCellPicked and not somePropPicked:
			# reset selection to blank
			print "Blank selection"
			self.output_link.GetCurrentSelection().RemoveAllNodes()
			self.output_link.InvokeEvent("AnnotationChangedEvent")
			self.applycolors1.Update()
			self.renWin.Render()
			
	def IcicleSelectionCallback(self, caller, event):
		# Defined for testing ID picking
		annSel = caller.GetCurrentSelection()
		# Note: When selection is cleared, the current selection does NOT contain any nodes
		if annSel.GetNumberOfNodes() > 0:
			idxArr = annSel.GetNode(0).GetSelectionList()
			print "ice ", VN.vtk_to_numpy(idxArr)
		else:
			print "ice Empty selection"
		
	def GetOutputAnnotationLink(self):
		return self.output_link
	
	def GetRenderWindow(self):
		return self.renWin
		
	def SetInteractorStyle(self, style):
		self.interactor = self.renWin.GetInteractor()
		self.istyle = style
		self.renWin.GetInteractor().SetInteractorStyle(self.istyle)
		# self.istyle.AddObserver("MouseMoveEvent", self.MouseMoveCallback)
		self.istyle.AddObserver("LeftButtonPressEvent", self.LeftButtonPressCallback)
		self.istyle.AddObserver("LeftButtonReleaseEvent", self.LeftButtonReleaseCallback)
		
	# def SetAnnotationLink(self, externalLink):
	#	self.output_link = externalLink
	#	self.view.GetRepresentation(0).SetAnnotationLink(self.output_link)
	#	self.output_link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)

if __name__ == "__main__":

	from data_source import DataSource

	# from tkFileDialog import askopenfilename
	# data_file = askopenfilename()
	data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/mnist12_1k_20100521.mat'
	
	# DataSource loads .mat file and can generate data from it for other views
	ds = DataSource(data_file)
	
	# All view classes have access to an instance of that data source for internal queries
	# Note that the only view which will pull and display data right away is the icicle view
	#  the other views need to be able to initialize without any data and only pull and show
	#  upon the first AnnotationChanged event...
	ice_class = IcicleNoView(ds)
	ice_class.GetRenderWindow().SetPosition(50,500)
	ice_class.GetRenderWindow().SetSize(630,470)
	ice_al_out = ice_class.GetOutputAnnotationLink()
	
	# Only need to Start() interactor for one view
	ice_class.GetRenderWindow().GetInteractor().Start()