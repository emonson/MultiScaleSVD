import vtk
import vtk.util.numpy_support as VN
import numpy as N

class IcicleView(object):

	def __init__(self, data_source):
		"""Icicle view constructor needs a valid DataSource and will pull data
		from it immediately."""
		
		self.ds = data_source
		
		SHRINK = 0.1
		THICK = 1.0
		
		tree = self.ds.GetTree()
		
		# Build view
		self.view = vtk.vtkIcicleView()
		self.view.SetRepresentationFromInput(tree)
		self.view.SetAreaSizeArrayName("num_in_vertex")
		self.view.SetAreaColorArrayName("scale")
		self.view.SetAreaLabelArrayName("blank")
		self.view.SetLabelPriorityArrayName("VertexDegree")
		self.view.SetAreaLabelVisibility(True)
		self.view.SetAreaHoverArrayName("vertex_ids")
		self.view.SetDisplayHoverText(True)
		self.view.SetShrinkPercentage(SHRINK)
		self.view.SetLayerThickness(THICK)
		self.view.UseGradientColoringOff()
		
		self.style = vtk.vtkInteractorStyleImage()
		self.view.SetInteractorStyle(self.style)
		
		# Parallel pipeline with no shrinkage to get TCoords
		TreeLevels = vtk.vtkTreeLevelsFilter()
		TreeLevels.SetInput(tree)
		
		VertexDegree = vtk.vtkVertexDegree()
		VertexDegree.SetInputConnection(TreeLevels.GetOutputPort(0))
		
		TreeAggregation = vtk.vtkTreeFieldAggregator()
		TreeAggregation.LeafVertexUnitSizeOff()
		TreeAggregation.SetField('size')
		TreeAggregation.SetInputConnection(VertexDegree.GetOutputPort(0))
		
		# Layout with shrinkage for generating geometry
		strategy = vtk.vtkStackedTreeLayoutStrategy()
		strategy.UseRectangularCoordinatesOn()
		strategy.SetRootStartAngle(0.0)
		strategy.SetRootEndAngle(15.0)
		strategy.SetRingThickness(THICK)	# layer thickness
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
		
		# Layout without shrinkage for generating texture coordinates
		strategy0 = vtk.vtkStackedTreeLayoutStrategy()
		strategy0.UseRectangularCoordinatesOn()
		strategy0.SetRootStartAngle(0.0)
		strategy0.SetRootEndAngle(15.0)
		strategy0.SetRingThickness(THICK)	# layer thickness
		strategy0.ReverseOn()
		strategy0.SetShrinkPercentage(SHRINK)
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
		WCimageData = self.ds.GetWaveletCoeffImage(XOrderedLeafIds)
		
		WCrange = N.array(WCimageData.GetPointData().GetScalars().GetRange())
		WCext = abs(WCrange.min()) if (abs(WCrange.min()) > abs(WCrange.max())) else abs(WCrange.max())
		# print WCext
		
		# Create blue to white to red lookup table
		lut = vtk.vtkLookupTable()
		lutNum = 256
		lut.SetNumberOfTableValues(lutNum)
		lut.Build()
		ctf = vtk.vtkColorTransferFunction()
		ctf.SetColorSpaceToDiverging()
		ctf.AddRGBPoint(0.0, 0, 0, 1.0)
		ctf.AddRGBPoint(1.0, 1.0, 0, 0)
		for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
			cc = ctf.GetColor(ss)
			lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
		lut.SetRange(-WCext,WCext)
		
		# Set up texture with lookup table for matrix polys
		tex = vtk.vtkTexture()
		tex.SetInput(WCimageData)
		tex.SetLookupTable(lut)
		
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
		ren = self.view.GetRenderer()
		ren.AddActor(act2)
		
		# NOTE: This is the one place I'm still hacking into the view -- to set the
		#   normal icicle view to wireframe and set the color...
		
		# Ren has an actor2d which is a scalar bar (edge)
		acts = ren.GetActors()
		# Acts has two actors -- graph blocks (0) and labels (1)
		act0 = acts.GetItemAsObject(0)
		
		act0.GetProperty().SetRepresentationToWireframe()
		act0.GetProperty().SetLineWidth(3.0)
		
		# Apply a theme to the views
		theme = vtk.vtkViewTheme.CreateMellowTheme()
		theme.SetPointHueRange(0,0)
		theme.SetPointSaturationRange(0.2,0.5)
		theme.SetPointValueRange(0.0,0.0)
		theme.SetPointAlphaRange(0.0,0.0)
		c = N.array([255,204,0])/255.0
		theme.SetSelectedPointColor(c[0],c[1],c[2])
		theme.SetBackgroundColor(0.1, 0.1, 0.06)
		theme.SetBackgroundColor2(0.25, 0.25, 0.2)
		self.view.ApplyViewTheme(theme)
		theme.FastDelete()
		
		# Connect the annotation link to the icicle representation
		rep = self.view.GetRepresentation(0)
		
		# Grab annotation link to monitor selection changes
		self.link = rep.GetAnnotationLink()
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
		self.link.GetCurrentSelection().GetNode(0).SetFieldType(3)		# Vertex
		self.link.GetCurrentSelection().GetNode(0).SetContentType(2)	# Pedigree Ids
		
		# TEST: Set up callback to test icicle view selection IDs
		self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)
		
		self.view.GetRenderer().GetActiveCamera().ParallelProjectionOn()		
		self.view.ResetCamera()
		self.view.Render()

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
		return self.link
	
	def GetView(self):
		return self.view
		
	# def SetAnnotationLink(self, externalLink):
	# 	self.link = externalLink
	# 	self.view.GetRepresentation(0).SetAnnotationLink(self.link)
	# 	self.link.AddObserver("AnnotationChangedEvent", self.IcicleSelectionCallback)

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
	ice_class = IcicleView(ds)
	ice_class.GetView().GetRenderWindow().SetPosition(50,500)
	ice_class.GetView().GetRenderWindow().SetSize(630,470)
	ice_al_out = ice_class.GetOutputAnnotationLink()
	
	# Only need to Start() interactor for one view
	ice_class.GetView().GetInteractor().Start()
