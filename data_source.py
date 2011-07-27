# With this version trying to set a separate actor in front of the icicle view which
# contains a "non-pickable" version of the textured icicle polys

import vtk
import scipy.io
import vtk.util.numpy_support as VN
import numpy as N
import os
import vtkvtg


class DataSource(object):
	"""Class that loads MultiScale SVD data from Matlab file and returns various pieces
	for visualization"""

	def __init__(self, filename=''):

		self.data_loaded = False

		# Built so it will automatically load a valid matlab file if given in constructor
		# Otherwise, call SetFileName('file.mat') and LoadData() separately

		if len(filename) > 0:
			self.data_file = os.path.abspath(filename)
		else:
			self.data_file = ''

		if os.path.isfile(self.data_file):
			try:
				self.LoadData()
			except:
				print "Came back from failed LoadData"

	def SetFileName(self, filename):
		"""Set file name manually for Matlab file. Can also do this in constructor."""

		if len(filename) > 0:
			self.data_file = os.path.abspath(filename)
		else:
			self.data_file = ''

		if not os.path.isfile(self.data_file):
			self.data_file = ''
			print "Supplied file does not exist"

	def LoadData(self):
		"""Routine that does the actual data loading and some format conversion.
		If a valid file name is given in the constructor, then this routine is called
		automatically. If you haven't given a file name in the constructor then you'll
		have to call SetFileName() before calling this."""

		# ----------
		# Load and construct whole graph and multi-resolution data from Matlab structure
		print 'Trying to load data set from .mat file... ', self.data_file

		if len(self.data_file) == 0:
			print "No data file name error!!!"
			raise IOError, "No data file name: Use SetFileName('file.mat') before LoadData()"

		print 'Trying to really load now...'
		try:
			MatInput = scipy.io.loadmat(self.data_file, struct_as_record=True, chars_as_strings=True)
		except:
			print 'loadmat crapping out for some reason...'
			raise IOError, "Can't load supplied matlab file"
			# return

		# Get variables out of Matlab structure
		print 'Transferring variables from Matlab structures'

		# GWTopts structure giving problems on loadmat when running in standalone app...
		# self.GWTopts = MatInput['GWTopts']

		# Test if original images are downsampled in this data
		if MatInput.has_key('X_down'):
			self.X_down = N.mat(MatInput['X_down'])
			self.downsampled = True
		else:
			# X has already been projected to D dim by PCA
			self.X = N.mat(MatInput['X'])
			self.downsampled = False
			# V = SVD result (X = (X0-mean)*V[:,:D])
			# This V already is only V(:,1:D) so we can use it whole
			self.V = N.mat(MatInput['V'])
			self.cm = N.mat(MatInput['cm'])	# not sure if should be matrix or array...

		# Various plain matrices
		# NOTE: Have to be careful of anything that can have a 0 value in Matlab
		# because it might be naturally imported as an unsigned int, and then
		# when you subtract 1 from it you don't get a negative number as you'd expect
		self.cp = (MatInput['cp'][0].astype('int16') - 1)	# change here to zero-based indexing
		self.IniLabels = (MatInput['IniLabels'][0] - 1)		# change here to zero-based indexing
		self.NumberInNet = MatInput['NumberInNet'][0]
		self.Scales = (MatInput['Scales'][0] - 1)					# zero-based
		self.IsALeaf = MatInput['IsALeaf'][0].astype('bool')
		self.LeafNodes = (MatInput['LeafNodes'][0] - 1)		# zero-based
		self.LeafNodesImap = (MatInput['LeafNodesImap'][0].astype('int16') - 1)		# zero-based
		self.EigenVecs = MatInput['EigenVecs'][0]
		self.CelWavCoeffs = MatInput['CelWavCoeffs']
		self.CelScalCoeffs = MatInput['CelScalCoeffs']

		# Load in category labels, but map them to sequential integers starting at 0
		if 'Labels' in MatInput:
			labels_array = MatInput['Labels'] # ncats x npoints 2d array
			self.cat_labels = N.zeros_like(labels_array)
			for ii in range(labels_array.shape[0]):
				cl_unique = set(labels_array[ii,:])
				cl_map = {}
				for jj,vv in enumerate(cl_unique):
					cl_map[vv] = jj
				self.cat_labels[ii,:] = N.array([cl_map[vv] for vv in labels_array[ii,:]])
			self.cat_labels_exist = True
		else:
			self.cat_labels_exist = False
		
		if self.cat_labels_exist:
			self.label_names = []
			# Check whether there are labels names and if there are the right number
			if ('LabelNames' in MatInput) and (MatInput['LabelNames'][0].size == self.cat_labels.shape[0]):
				names_array = MatInput['LabelNames'][0]
				for name_ar in names_array:
					self.label_names.append(name_ar[0] + '_ids')
			# Else generate fake names
			else:
				for ii in range(self.cat_labels.shape[0]):
					self.label_names.append('label_' + str(ii) + '_ids')

		# Need the max number of dims at each scale to fill in zeros for pcoords plot
		# Create copies for both wavelet coeffs and scaling functions since these often
		#  have different dimensionalities
		self.WavMaxDim = N.zeros(self.CelWavCoeffs.shape[1],dtype='int')
		for row in range(self.CelWavCoeffs.shape[0]):
			for col in range(self.CelWavCoeffs.shape[1]):
				if self.WavMaxDim[col] < self.CelWavCoeffs[row,col].shape[1]:
					self.WavMaxDim[col] = self.CelWavCoeffs[row,col].shape[1]
		self.ScalMaxDim = N.zeros(self.CelScalCoeffs.shape[1],dtype='int')
		for row in range(self.CelScalCoeffs.shape[0]):
			for col in range(self.CelScalCoeffs.shape[1]):
				if self.ScalMaxDim[col] < self.CelScalCoeffs[row,col].shape[1]:
					self.ScalMaxDim[col] = self.CelScalCoeffs[row,col].shape[1]

		# Gather helpful statistics to be used by other classes
		print 'Calulating extrema of coefficients'
		self.WavCoeffMax = -1e200
		self.WavCoeffMin = 1e200
		self.ScalCoeffMax = -1e200
		self.ScalCoeffMin = 1e200
		for ii in range(self.CelWavCoeffs.shape[0]):
			for jj in range(self.CelWavCoeffs.shape[1]):
				if (self.CelWavCoeffs[ii,jj].size > 0):
					wmax = N.amax(self.CelWavCoeffs[ii,jj])
					wmin = N.amin(self.CelWavCoeffs[ii,jj])
					if (wmax > self.WavCoeffMax): self.WavCoeffMax = wmax
					if (wmin < self.WavCoeffMin): self.WavCoeffMin = wmin
				if (self.CelScalCoeffs[ii,jj].size > 0):
					smax = N.amax(self.CelScalCoeffs[ii,jj])
					smin = N.amin(self.CelScalCoeffs[ii,jj])
					if (smax > self.ScalCoeffMax): self.ScalCoeffMax = smax
					if (smin < self.ScalCoeffMin): self.ScalCoeffMin = smin

		# NOTE: gW and Data are class numpy.ndarray
		#		MatInput is just a dict, so can directly look for variables there

		if MatInput.has_key('imR') and MatInput.has_key('imC'):
			print 'Grabbing image dimensions from matlab file'
			self.imR = MatInput['imR'][0,0]
			self.imC = MatInput['imC'][0,0]
		else:
			# Right now Matlab data doesn't have any record of original image dimensions
			# NOTE: Hard coding shape for now!
			print 'Hacking image dimensions from file name'
			if (self.data_file.find('mnist') >= 0):
				self.imR = 28	# rows
				self.imC = 28	# cols
			elif (self.data_file.find('frey') >= 0):
				self.imR = 20	# rows
				self.imC = 28	# cols
			elif (self.data_file.find('olivetti') >= 0):
				self.imR = 64	# rows
				self.imC = 64	# cols
			else:
				self.imR = 20
				self.imC = 20
				print 'Could not find matching file name -- probably wrong image dimensions!!!!'
		if self.downsampled:
			self.imR_down = MatInput['imR_down'][0,0]
			self.imC_down = MatInput['imC_down'][0,0]

		# NumPts = Number of data points (here number of individual images)
		self.NumPts = self.IniLabels.shape[0]
		self.AmbientDimension = MatInput['AmbientDimension'][0,0]	# used to call this D

		# Manifold dimensionality variable now, not fixed...
		# self.ManifoldDim = self.gW['ManifoldDimension'][0,0][0,0]

		# Converting cell arrays to lists of numpy arrays
		self.PointsInNet = []	# Points In Net
		for ii in range(MatInput['PointsInNet'].shape[1]):
			self.PointsInNet.append(MatInput['PointsInNet'][0,ii][0]-1)	# 0-based indices

		if self.downsampled:
			self.Centers_down = []
			self.WavBases_down = []
			self.ScalFuns_down = []
			for ii in range(MatInput['PointsInNet'].shape[1]):
				self.Centers_down.append(N.mat(MatInput['Centers_down'][0,ii]))			# matrix
				self.WavBases_down.append(N.mat(MatInput['WavBases_down'][0,ii]))			# matrix
				self.ScalFuns_down.append(N.mat(MatInput['ScalFuns_down'][0,ii]))			# matrix
		else:
			self.ScalFuns = []	# Scaling functions
			self.WavBases = []	# Wavelet bases
			self.Centers = []	# Center of each node
			for ii in range(MatInput['PointsInNet'].shape[1]):
				self.ScalFuns.append(N.mat(MatInput['ScalFuns'][0,ii]))			# matrix
				self.WavBases.append(N.mat(MatInput['WavBases'][0,ii]))			# matrix
				self.Centers.append(N.mat(MatInput['Centers'][0,ii][0])) 		# matrix


		# J = Total number of scales
		# self.J = self.Scales.max()

		# Creating a storage space for ordering of leaf nodes in icicle view with
		# default values of ordering according to Matlab-saved LeafNodes
		ice_leaf_ids = self.LeafNodes
		ice_leaf_xmins = N.arange(ice_leaf_ids.size)
		ice_ids_mapped = self.LeafNodesImap[ice_leaf_ids]
		self.mapped_leaf_pos = N.zeros_like(ice_leaf_xmins)
		self.mapped_leaf_pos[ice_ids_mapped] = ice_leaf_xmins

		# Flag to set whether generic routines should return Wavelet or Scaling Function
		# coefficients / images -- "wav" or "scal"
		self.SetCoeffSource("wavelet")
		
		# -- Wordle --
		# Flag to indicate whether images should be generated directly
		# or by passing terms through QtWordleView
		self.WordleImages = False
		self.qinit = None
		self.WordleView = None
		self.WordleTable = None
		self.Terms = None
		
		if 'terms' in MatInput.keys():
			self.WordleImages = True
			mat_terms = MatInput['terms'].T[0]
			self.Terms = vtk.vtkStringArray()
			self.Terms.SetName('dictionary')
			self.Terms.SetNumberOfComponents(1)
			for term in mat_terms:
				self.Terms.InsertNextValue(term[0])

			# Init Table and put in some sample data that will be replaced later
			basis_idx = 0
			
			coeffs = VN.numpy_to_vtk(self.WavBases[basis_idx][:,0]*100, deep=True)
			coeffs.SetName('coefficient')
			c_sign = VN.numpy_to_vtk(N.sign(self.WavBases[basis_idx][:,0]), deep=True)
			c_sign.SetName('sign')
			
			# Create a table with some points in it...
			self.WordleTable = vtk.vtkTable()
			self.WordleTable.AddColumn(self.Terms)
			self.WordleTable.AddColumn(coeffs)
			self.WordleTable.AddColumn(c_sign)
			
			# self.qinit = vtk.vtkQtInitialization()
			self.WordleView = vtkvtg.vtkQtWordleView()

			vt = vtk.vtkViewTheme()
			lut = vtk.vtkLookupTable()
			lut.SetHueRange(0, 0.66)
			lut.SetValueRange(0.7, 0.7)
			lut.SetSaturationRange(1, 1)
			lut.Build()
			# Set value for no color by array
			vt.SetPointColor(0,0,0)
			# Set LUT for color by array
			vt.SetPointLookupTable(lut)
			# ViewTheme Background color is black by default
			vt.SetBackgroundColor(1,1,1)
			
			self.WordleView.SetFieldType(vtkvtg.vtkQtWordleView.ROW_DATA)
			self.WordleView.AddRepresentationFromInput(self.WordleTable)
			self.WordleView.SetColorByArray(True)
			self.WordleView.ApplyViewTheme(vt)
			self.WordleView.SetColorArrayName('sign')
			self.WordleView.SetTermsArrayName('dictionary')
			self.WordleView.SetSizeArrayName('coefficient')
			self.WordleView.SetOutputImageDataDimensions(200, 200)
			self.WordleView.SetMaxNumberOfWords(100);
			self.WordleView.SetFontFamily("Rockwell")
			self.WordleView.SetFontStyle(vtkvtg.vtkQtWordleView.StyleNormal)
			self.WordleView.SetFontWeight(99)
			
			# self.WordleView.SetOrientation(vtkvtg.vtkQtWordleView.HORIZONTAL)
			self.WordleView.SetOrientation(vtkvtg.vtkQtWordleView.MOSTLY_HORIZONTAL)
			# self.WordleView.SetOrientation(vtkvtg.vtkQtWordleView.HALF_AND_HALF)
			# self.WordleView.SetOrientation(vtkvtg.vtkQtWordleView.MOSTLY_VERTICAL)
			# self.WordleView.SetOrientation(vtkvtg.vtkQtWordleView.VERTICAL)
	
			# self.WordleView.SetLayoutPathShape(vtkvtg.vtkQtWordleView.CIRCULAR_PATH)
			self.WordleView.SetLayoutPathShape(vtkvtg.vtkQtWordleView.SQUARE_PATH)
			
		self.data_loaded = True

	# ---------------------------------------
	# ---------------------------------------
	def GetTree(self):
		"""Returns a full vtkTree based on data loaded in LoadData()."""

		if self.data_loaded:

			vertex_id = vtk.vtkIdTypeArray()
			vertex_id.SetName('vertex_ids')
			for ii in range(len(self.cp)):
				vertex_id.InsertNextValue(ii)

			NINvtk = VN.numpy_to_vtk(self.NumberInNet, deep=True)
			NINvtk.SetName('num_in_vertex')
			SCALESvtk = VN.numpy_to_vtk(self.Scales, deep=True)
			SCALESvtk.SetName('scale')

			# This array will default to empty strings
			BLANKvtk = vtk.vtkStringArray()
			BLANKvtk.SetNumberOfComponents(1)
			BLANKvtk.SetNumberOfTuples(self.NumberInNet.shape[0])
			BLANKvtk.SetName('blank')

			# Build tree out of CP list of "is a child of"
			#	remembering that Matlab indices are 1-based and numpy/VTK 0-based
			print 'Building graph'
			dg = vtk.vtkMutableDirectedGraph()
			edge_id = vtk.vtkIdTypeArray()
			edge_id.SetName('edge_ids')
			for ii in range(self.cp.size):
				dg.AddVertex()
			for ii in range(self.cp.size):
				if self.cp[ii] > 0:		# CP already zero-based
					dg.AddGraphEdge(self.cp[ii],ii)		# Method for use with wrappers -- AddEdge() in C++
					edge_id.InsertNextValue(ii)

			dg.GetVertexData().AddArray(NINvtk)
			dg.GetVertexData().AddArray(SCALESvtk)
			dg.GetVertexData().AddArray(vertex_id)
			dg.GetVertexData().SetActiveScalars('scale')
			dg.GetVertexData().SetActivePedigreeIds('vertex_ids')
			dg.GetEdgeData().AddArray(edge_id)
			dg.GetEdgeData().SetActivePedigreeIds('edge_ids')

			tree = vtk.vtkTree()
			tree.CheckedShallowCopy(dg)

			return tree

		else:
			raise IOError, "Can't get tree until data is loaded successfully"

	# ---------------------------------------
	def SetCoeffSource(self, source_name):
		"""Set whether generic routines should return Wavelet or Scaling Function
		coefficents. Use "wav" or "scal", but it will work with longer, capitalized versions.
		"""
		if source_name.lower().startswith('wav'):
			self.coeff_source = 'wav'
			self.ScaleMaxDim = self.WavMaxDim
			self.CelCoeffs = self.CelWavCoeffs
			if self.downsampled:
				self.Bases_down = self.WavBases_down
			else:
				self.Bases = self.WavBases
		elif source_name.lower().startswith('sca'):
			self.coeff_source = 'scal'
			self.ScaleMaxDim = self.ScalMaxDim
			self.CelCoeffs = self.CelScalCoeffs
			if self.downsampled:
				self.Bases_down = self.ScalFuns_down
			else:
				self.Bases = self.ScalFuns
		else:
			print "Error: Unknown coefficient source. Use 'wavelet' or 'scaling'."

	# ---------------------------------------
	def GetCoeffSource(self):
		"""Get whether generic routines will return Wavelet or Scaling Function
		coefficents: 'wavelet' or 'scaling'
		"""
		if self.coeff_source == 'wav':
			return 'wavelet'
		else:
			return 'scaling'

	# ---------------------------------------
	def GetCoeffRange(self):
		"""Returns a tuple containing the range of values (min,max) of all the wavelet coefficients
		or scaling coefficients depending on which has been set in SetCoeffSource. Default is Wavelet.
		"""
		if self.coeff_source == 'wav':
			return (self.WavCoeffMin,self.WavCoeffMax)
		else:
			return (self.ScalCoeffMin,self.ScalCoeffMax)

	# ---------------------------------------
	def GetCategoryLabelRange(self, idx=0):
		"""Returns a tuple containing the range of values (min,max) of
		the category labels (which have been mapped above to sequential integers).
		For multiple category labels you can pass an index value (0-based), which
		defaults to zero.
		"""
		if self.cat_labels_exist and (idx < self.cat_labels.shape[0]):
			return (N.min(self.cat_labels[idx,:]), N.max(self.cat_labels[idx,:]))
		else:
			# TODO: Should think about what to return here...
			# TODO: Should also probably put up an error...
			return (0,0)

	# ---------------------------------------
	def GetCoeffImages(self, ice_leaf_ids=None, ice_leaf_xmins=None ):
		"""Returns a list of vtkImageData 2D image with the wavelet or scaling function
		coefficients at all dimensions for all nodes. 
		If you give the positions and IDs of the leaf nodes, as laid out by
		the icicle view, then the matrix will be sorted accordingly
		so the image should be correct for the icicle view."""

		if self.data_loaded:

			if (ice_leaf_ids is not None) and (len(ice_leaf_ids) != len(self.LeafNodes)):

				raise ValueError, "Number of leaves in icicle view doesn't match leaf nodes in tree"

			else:

				# Matlab method for appending wavelet coeffients together
				# Need to use this if the icicle view layout may have reordered the nodes
				#
				# num_nodes = length(gW.cp);
				# LeafNodesImap(gW.LeafNodes) = 1:length(gW.LeafNodes);
				# NodeWavCoeffs = cell(1,num_nodes);
				#
				# for node_idx = 1:num_nodes,
				# 	offspring = [node_idx get_offspring(gW.cp, node_idx)];
				# 	relevantLeafNodes = offspring(logical(gW.isaleaf(offspring)));
				# 	NodeWavCoeffs{node_idx} = cat(1, Data.CelWavCoeffs{LeafNodesImap(relevantLeafNodes), gW.Scales(node_idx)});
				# end

				# NOTE: The complication with this version is that we need to be able to handle
				# cases where the icicle view has arranged the tree & leaf nodes in a different
				# order than Matlab stored the leaf nodes (and indexed the wavelet coeffs)

				# If the caller wants the wavelet coeffs returned in a certain order, then
				# they need to supply positions with LeafIds, otherwise just use the
				# original order of the leaf nodes
				if ice_leaf_ids is None:
					ice_leaf_ids = self.LeafNodes
					ice_leaf_xmins = N.arange(ice_leaf_ids.size)

				# Create an array with correct assignments of ice_pos for LeafNodes index
				# and store for use by other routines
				ice_ids_mapped = self.LeafNodesImap[ice_leaf_ids]
				self.mapped_leaf_pos = N.zeros_like(ice_leaf_xmins)
				self.mapped_leaf_pos[ice_ids_mapped] = ice_leaf_xmins

				WC_imagedata_list = []

				# Looping through by node_id in tree
				for node_id in range(self.cp.size):
					leaf_offspring = N.array(list(self.get_leaf_children(self.cp, node_id)))
					offspring_idxs = self.LeafNodesImap[leaf_offspring]
					offspring_pos = self.mapped_leaf_pos[offspring_idxs]
					sorted_offspring_pos_idxs = N.argsort(offspring_pos)
					sorted_offspring_idxs = offspring_idxs[sorted_offspring_pos_idxs]
					# Need to reverse the order up-down of wav coeffs
					img_list = list(pp[::-1,:] for pp in self.CelCoeffs[sorted_offspring_idxs, self.Scales[node_id]])
					# Need to reverse the list to get in the right order
					img_list.reverse()
					# Need to transpose the concatenated matrices
					img = N.concatenate(tuple(img_list), axis=0).T

					# Create vtkImageData out of WavCoeffs for texturing icicle view tree
					# .copy() is to force the array to be contiguous for numpy_to_vtk
					# deep=True should keep reference around even after numpy array is destroyed

					# Need to reverse the order again
					WCvtk = VN.numpy_to_vtk(img.ravel()[::-1].copy(), deep=True)
					WCvtk.SetName('Coeffs')

					WCimageData = vtk.vtkImageData()
					WCimageData.SetOrigin(0,0,0)
					WCimageData.SetSpacing(1,1,1)
					WCimageData.SetDimensions(img.shape[1],img.shape[0],1)
					WCimageData.GetPointData().AddArray(WCvtk)
					WCimageData.GetPointData().SetActiveScalars('Coeffs')

					WC_imagedata_list.append(WCimageData)

				return WC_imagedata_list

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetIdsFractionalPosition(self, XOrderedLeafIds=None ):
		"""Returns a vtkImageData 2D image with the wavelet coefficients at all dimensions
		and scales. If you give an ordered list of the leaf node IDs, then the matrix will
		be sorted accordingly so the image should be correct for the icicle view."""

		if self.data_loaded:

			SortedLeafIdxArray = N.array([],dtype='uint16')
			# If the caller wants the wavelet coeffs returned in a certain order, then
			# they need to supply a sorted list of the LeafIds
			if XOrderedLeafIds is not None:
				# Create an array holding the indices of the leaf vertices in the proper order
				for ii in range(XOrderedLeafIds.size):
					SortedLeafIdxArray = N.concatenate((SortedLeafIdxArray,self.PointsInNet[XOrderedLeafIds[ii]]))

			else:
				# Assume that self.leafNodes is in the proper order
				for ii in range(self.LeafNodes.size):
					SortedLeafIdxArray = N.concatenate((SortedLeafIdxArray,self.PointsInNet[self.LeafNodes[ii]]))

			return (SortedLeafIdxArray.argsort()+0.5)/float(SortedLeafIdxArray.size)

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetNodeAllScaleCoeffTable(self, node_id, dim_limit=0):
		"""Returns a table of all the wavelet coefficients for a single tree
		icicle view) node at all scales for plotting on a parallel coodinates plot.
		For right now, padding with zeros for unused dimensions in some nodes...
		This version supports variable dimensionality and category labels."""

		# NOTE: Right now filling with zeros to max scale even if chosen IDs don't
		# include leaf nodes which run to the max scale...

		if self.data_loaded:

			# For a given node_id, get PIN and then extract all coeffs at every scale
			# Columns of table will be rows of the WavCoeffsOrig matrix

			IDarray = self.PointsInNet[node_id]

			# Really want zeros padded only to max dim at each scale for _this_ particular
			# node and its children, not the max that exists in the entire data set.
			leaf_children = list(self.get_leaf_children(self.cp, node_id))
			mapped_leaf_children = [self.LeafNodesImap[nod] for nod in leaf_children]
			# Even empty arrays in self.CelCoeffs give back okay shape[1] (0)
			dims_arr_lin = N.array(list(arr.shape[1] for id in mapped_leaf_children for arr in self.CelCoeffs[id,:]))
			dims_arr = dims_arr_lin.reshape((-1,self.ScaleMaxDim.size))
			scale_max_dim = dims_arr.max(axis=0)

			# Need to trim off zeros at the end or they cause problems...
			scale_max_dim = scale_max_dim[N.nonzero(scale_max_dim)]

			# TEST: Try limiting to a max number of dims...
			if (dim_limit > 0):
				scale_max_dim[scale_max_dim > dim_limit] = dim_limit

			wav_coeffs = N.zeros((len(IDarray),sum(scale_max_dim)))
			# Append all zero-padded rows gathered from CelWavCoeffs
			for ii,data_id in enumerate(IDarray):
				leaf_node = self.IniLabels[data_id]
				row = N.nonzero(self.PointsInNet[leaf_node]==data_id)[0][0]	# final zero turns array->scalar
				mapped_node_idx = self.LeafNodesImap[leaf_node]

				# Skip any empty arrays
				wav_row_tuple = tuple(arr[row,:] for arr in self.CelCoeffs[mapped_node_idx,:] if arr.size != 0)
				# Create zero-padded arrays
				zero_row_tuple = tuple(N.zeros(sc) for sc in scale_max_dim)
				# And transfer over values
				for zz in range(len(wav_row_tuple)):
					sz_limit = min([scale_max_dim[zz], len(wav_row_tuple[zz])])
					zero_row_tuple[zz][:sz_limit] = wav_row_tuple[zz][:sz_limit]

				wav_coeffs[ii,:] = N.concatenate(zero_row_tuple, axis=1)

			table = vtk.vtkTable()
			col_idx = 0
			for scale,maxdim in enumerate(scale_max_dim):
				for ii in range(maxdim):
					column = VN.numpy_to_vtk(wav_coeffs[:,col_idx].copy(), deep=True)
					column.SetName(str(scale) + '.' + str(ii))
					table.AddColumn(column)
					col_idx += 1

			# Trying to set PedigreeIds to that parallel coords selections have correct IDs
			IDvtk = VN.numpy_to_vtk(IDarray, deep=True)
			IDvtk.SetName('pedigree_ids')
			table.AddColumn(IDvtk)
			table.GetRowData().SetActivePedigreeIds('pedigree_ids')

			# Adding in category labels
			# Name needs to end in _ids so plots will ignore it
			if self.cat_labels_exist:
				for ii in range(self.cat_labels.shape[0]):
					CATvtk = VN.numpy_to_vtk(self.cat_labels[ii,IDarray], deep=True)
					CATvtk.SetName(self.label_names[ii])
					table.AddColumn(CATvtk)

			return table, scale_max_dim.tolist()

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetNodeAllScaleDDimCoeffTable(self, node_id, D):
		"""Returns a table of all the wavelet coefficients for a single tree
		icicle view) node at all scales for plotting on a parallel coodinates plot.
		* * This version used for testing PCoords with each scale abbreviated to first
		D dimensions... * *"""

		# NOTE: Right now filling with zeros to max scale even if chosen IDs don't
		# include leaf nodes which run to the max scale...

		if self.data_loaded:

			# For a given node_id, get PIN and then extract all coeffs at every scale
			# Columns of table will be rows of the WavCoeffsOrig matrix

			IDarray = self.PointsInNet[node_id]

			wav_coeffs = N.zeros((len(IDarray),D*len(self.ScaleMaxDim)))
			# Append all zero-padded rows gathered from CelWavCoeffs
			for ii,data_id in enumerate(IDarray):
				leaf_node = self.IniLabels[data_id]
				row = N.nonzero(self.PointsInNet[leaf_node]==data_id)[0][0]	# final zero turns array->scalar
				mapped_node_idx = self.LeafNodesImap[leaf_node]

				# Skip any empty arrays
				wav_row_tuple = tuple(arr[row,:] for arr in self.CelCoeffs[mapped_node_idx,:] if arr.size != 0)
				# Create zero-padded arrays
				zero_row_tuple = tuple(N.zeros(D) for sc in self.ScaleMaxDim)
				# And transfer over values
				for zz in range(len(wav_row_tuple)):
					if len(wav_row_tuple[zz]) >= D:
						zero_row_tuple[zz][:] = wav_row_tuple[zz][:D]
					else:
						d = len(wav_row_tuple[zz])
						zero_row_tuple[zz][:d] = wav_row_tuple[zz][:]

				wav_coeffs[ii,:] = N.concatenate(zero_row_tuple, axis=1)

			table = vtk.vtkTable()
			col_idx = 0
			for scale,maxdim in enumerate(self.ScaleMaxDim):
				for ii in range(D):
					column = VN.numpy_to_vtk(wav_coeffs[:,col_idx].copy(), deep=True)
					column.SetName(str(scale) + '.' + str(ii))
					table.AddColumn(column)
					col_idx += 1

			# Trying to set PedigreeIds to that parallel coords selections have correct IDs
			IDvtk = VN.numpy_to_vtk(IDarray, deep=True)
			IDvtk.SetName('pedigree_ids')
			table.AddColumn(IDvtk)
			table.GetRowData().SetActivePedigreeIds('pedigree_ids')

			return table

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetNodeOneScaleCoeffTable(self, node_id):
		"""Returns a table of the wavelet coefficients at a single node at a single
		scale for plotting on a scatter plot. Relying on icicle_view already having
		called GetWaveletCoeffImages() with correct positions of leaf nodes in view,
		otherwise just using original Matlab-saved LeafNodes ordering.
		This version supports category labels."""

		if self.data_loaded:

			# For a given node_id, concatenate wavelet coeffs in proper order
			# (according to leaf node positions in icicle view if it was set already)
			# Columns of table will be rows of the wavelet coeffs image

			scale = self.Scales[node_id]

			leaf_offspring = N.array(list(self.get_leaf_children(self.cp, node_id)))
			offspring_idxs = self.LeafNodesImap[leaf_offspring]
			offspring_pos = self.mapped_leaf_pos[offspring_idxs]
			sorted_offspring_pos_idxs = N.argsort(offspring_pos)
			sorted_offspring_idxs = offspring_idxs[sorted_offspring_pos_idxs]
			img_tuple = tuple(pp for pp in self.CelCoeffs[sorted_offspring_idxs, scale])
			# The image comes out with shape (npts, ndims)
			# May need to reorder (reverse) this...?
			img = N.concatenate(img_tuple, axis=0)

			table = vtk.vtkTable()
			for ii in range(img.shape[1]):
				column = VN.numpy_to_vtk(img[:,ii].copy(), deep=True)
				column.SetName(str(scale) + '.' + str(ii))
				table.AddColumn(column)

			IDtuple = tuple(self.PointsInNet[xx] for xx in self.LeafNodes[sorted_offspring_idxs])
			IDarray = N.concatenate(IDtuple)
			# Trying to set PedigreeIds to that parallel coords selections have correct IDs
			IDvtk = VN.numpy_to_vtk(IDarray, deep=True)
			IDvtk.SetName('pedigree_ids')
			table.AddColumn(IDvtk)
			table.GetRowData().SetActivePedigreeIds('pedigree_ids')

			# Adding in category labels
			# Name needs to end in _ids so plots will ignore it
			if self.cat_labels_exist:
				for ii in range(self.cat_labels.shape[0]):
					CATvtk = VN.numpy_to_vtk(self.cat_labels[ii,IDarray], deep=True)
					CATvtk.SetName(self.label_names[ii])
					table.AddColumn(CATvtk)

			return table


		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetNodeBasisImages(self, node_id):
		"""Returns a vtkImageData of all wavelet or scaling function
		basis images for a given node."""

		if self.data_loaded:
			if self.WordleImages:
				# Scaling functions coeffs are defined wrt parent node scaling functions...
				# TODO: Switch this back when back to computing CelScalCoeffs rather than
				#   CelTangCoeffs...
				if self.coeff_source == 'scal' and self.cp[node_id] >= 0:
					node_id = self.cp[node_id]
	
				# Need to create separate images (Z) for each column of matrix result
				# Bases is D x N matrix
				image_cols = self.Bases[node_id]
				
				imgAppend = vtk.vtkImageAppend()
				imgAppend.SetAppendAxis(2)	# Z

				for ii in range(self.Bases[node_id].shape[1]):
					
					coeffs = VN.numpy_to_vtk(self.Bases[node_id][:,ii]*100, deep=True)
					coeffs.SetName('coefficient')
					c_sign = VN.numpy_to_vtk(N.sign(self.Bases[node_id][:,ii]), deep=True)
					c_sign.SetName('sign')
					
					self.WordleTable.RemoveColumn(2)
					self.WordleTable.RemoveColumn(1)
					self.WordleTable.AddColumn(coeffs)
					self.WordleTable.AddColumn(c_sign)
					self.WordleView.RemoveAllRepresentations()
					self.WordleView.AddRepresentationFromInput(self.WordleTable)
					
					self.WordleTable.Modified()
					
					img = vtk.vtkImageData()
					img.DeepCopy(self.WordleView.GetImageData())
					img.GetPointData().GetScalars().SetName('DiffIntensity')
					imgAppend.AddInput(img)
				
				imgAppend.Update()
				out_img = vtk.vtkImageData()
				out_img.DeepCopy(imgAppend.GetOutput())
				writer = vtk.vtkXMLImageDataWriter()
				writer.SetFileName("out.vti")
				writer.SetInput(out_img)
				writer.Write()
				return out_img
				
			else:
				
				# Scaling functions coeffs are defined wrt parent node scaling functions...
				# TODO: Switch this back when back to computing CelScalCoeffs rather than
				#   CelTangCoeffs...
				if self.coeff_source == 'scal' and self.cp[node_id] >= 0:
					node_id = self.cp[node_id]
	
				# %% Display all detail coordinates for a given leaf node
				#
				# imagesc(reshape(V(:,1:D)*gW.ScalFuns{i}, self.imR,[]));
				#
				# Need to create separate images (Z) for each column of matrix result
	
				if self.downsampled:
					image_cols = self.WavBases_down[node_id]
					imR = self.imR_down
					imC = self.imC_down
				else:
					# V now already chopped to AmbientDimension
					# Compute all detail images for that dimension
					# print "DS Calculating center image"
					# print node_id, self.Centers[node_id].shape, self.V.T.shape, self.cm.shape
					image_cols = self.V*self.Bases[node_id]
					imR = self.imR
					imC = self.imC
	
				# To make it linear, it is the correct order (one image after another) to .ravel()
				images_linear = N.asarray(image_cols.T).ravel()
	
				intensity = VN.numpy_to_vtk(images_linear, deep=True)
				intensity.SetName('DiffIntensity')
	
				imageData = vtk.vtkImageData()
				imageData.SetOrigin(0,0,0)
				imageData.SetSpacing(1,1,1)
				imageData.SetDimensions(imR, imC, image_cols.shape[1])
				imageData.GetPointData().AddArray(intensity)
				imageData.GetPointData().SetActiveScalars('DiffIntensity')
	
				return imageData
	
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetNodeCenterImage(self, node_id):
		"""Returns a vtkImageData of the center image for a given node."""

		if self.data_loaded:
			# if self.WordleImages:
			if False:
				pass
				
			else:
				
				# imagesc(reshape(gW.Centers{1}*V(:,1:D)'+cm,28,[]))
	
				if self.downsampled:
					image_col = self.Centers_down[node_id]
					imR = self.imR_down
					imC = self.imC_down
				else:
					# V now already chopped to AmbientDimension
					# Compute all detail images for that dimension
					# print "DS Calculating center image"
					# print node_id, self.Centers[node_id].shape, self.V.T.shape, self.cm.shape
					image_col = self.Centers[node_id]*self.V.T + self.cm
					imR = self.imR
					imC = self.imC
	
				# print "DS done calculating center image"
				# To make it linear, it is the correct order (one image after another) to .ravel()
				image_linear = N.asarray(image_col.T).ravel()
	
				intensity = VN.numpy_to_vtk(image_linear, deep=True)
				intensity.SetName('Intensity')
	
				imageData = vtk.vtkImageData()
				imageData.SetOrigin(0,0,0)
				imageData.SetSpacing(1,1,1)
				imageData.SetDimensions(imR, imC, 1)
				imageData.GetPointData().AddArray(intensity)
				imageData.GetPointData().SetActiveScalars('Intensity')
	
				return imageData
	
		else:
			raise IOError, "Can't get image until data is loaded successfully"
	
	# ---------------------------------------
	def GetProjectedImages(self, IDlist):
		"""Given a list of IDs selected from a parallel coordinates plot, returns
		a vtkImageData with all of the projected (reduced dimensionality by SVD) images
		for those IDs. (e.g. typically 120 dim rather than original 768 dim for MNIST digits)"""

		if self.data_loaded:
			# if self.WordleImages:
			if False:
				pass
				
			else:
				
				# X_orig = X*V(:,1:GWTopts.AmbientDimension)'+repmat(cm, size(X,1),1);
	
				if self.downsampled:
					X_orig = self.X_down[IDlist,:]
					imR = self.imR_down
					imC = self.imC_down
				else:
					# V now already chopped to AmbientDimension
					Xtmp = self.X[IDlist,:]*self.V.T
		
					# numpy should automatically do tiling!!
					X_orig = Xtmp + self.cm
					# X_orig = Xtmp + N.tile(self.cm,(Xtmp.shape[0],1))	# tile ~ repmat
					imR = self.imR
					imC = self.imC
	
				# To make it linear, it is the correct order (one image after another) to .ravel()
				X_linear = N.asarray(X_orig).ravel()
	
				# If we want to rearrange it into a stack of images in numpy
				# X_im = N.asarray(X_orig).reshape(Xtmp.shape[0],self.imR,-1)
	
				# Going ahead and using numpy_support here...  Much faster!!!
				Xvtk = VN.numpy_to_vtk(X_linear, deep=True)	# even with the (necessary) deep copy
				Xvtk.SetName('Intensity')
	
				imageData = vtk.vtkImageData()
				imageData.SetOrigin(0,0,0)
				imageData.SetSpacing(1,1,1)
				imageData.SetDimensions(imR, imC, X_orig.shape[0])
				imageData.GetPointData().AddArray(Xvtk)
				imageData.GetPointData().SetActiveScalars('Intensity')
	
				return imageData
	
		else:
			raise IOError, "Can't get image until data is loaded successfully"
	
	# ---------------------------------------
	def GetOriginalImages(self, IDlist):
		"""Given a list of IDs selected from a parallel coordinates plot, returns
		a vtkImageData with all of the original images for those IDs."""

		# * * * OLD FIXED-DIM VERSION * * *

		if self.data_loaded:

			for ii in IDlist:
				pass
				# i = sample; % 39 when digit =1
				#
				# %% Add original image
				#
				# imagesc(reshape(X0(i,:), self.imR,[]));

			return

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetMultiResolutionImages(self, ID):
		"""Given a single ID selected from an image flow view, returns
		a vtkImageData with all of the image approximations at all scales for
		that ID."""

		# * * * OLD FIXED-DIM VERSION * * *

		if self.data_loaded:


			# NEW Matlab version of how projections are calculated
			#
			# for i = 1:nLeafNodes
			#
			# 	net = leafNodes(i);
			# 	netPts = find(gW.IniLabels == net);
			# 	nPts = length(netPts);
			# 	j_max = gW.Scales(net);
			# 	gWCentersnet = repmat(gW.Centers{net},nPts,1);
			#
			# 	Data.Projections(netPts,:,j_max) = Data.CelScalCoeffs{i,j_max}*gW.ScalFuns{net}' + gWCentersnet;

			# OLD % Calculate and display approximations at all scales
			#
			# i = sample; % 39 when digit =1
			#
			# for j = J:-1:1
			# 	X_approx = Data.Projections(:,:,j);
			# 	X_img = X_approx*V(:,1:D)'+repmat(cm, size(X_approx,1),1);
			# 	imagesc(reshape(X_img(i,:), self.imR,[]));
			# end
			#

			# TODO: Not saving Projections now (for space and loading) so need to calculate
			# 	them from scaling functions and coefficients as in NEW calc above...

			# In principle, can calculate projection at scale which doesn't exist for a particular
			# image, so I'm testing which scales are in the chain of nodes for that image
			# and only calulating multires for those scales so the multi-res and detail images
			# will match up...

			leafNode = self.IniLabels[ID]
			chain = self.find_path_down_the_tree(leafNode)
			scales = self.Scales[chain]

			# PROJ[ID][:,scale]

			Xtmp = N.mat(self.PROJ[ID][:,scales]).T*self.V[:,:self.D].T		# Xtmp.shape = (scales,imR*imC)
			X_img = Xtmp + N.tile(self.cm,(Xtmp.shape[0],1))	# tile ~ repmat

			# To make it linear, it is the correct order (one image after another) to .ravel()
			X_linear = N.asarray(X_img).ravel()

			# Going ahead and using numpy_support here...  Much faster!!!
			Xvtk = VN.numpy_to_vtk(X_linear, deep=True)	# even with the (necessary) deep copy
			Xvtk.SetName('Intensity')

			imageData = vtk.vtkImageData()
			imageData.SetOrigin(0,0,0)
			imageData.SetSpacing(1,1,1)
			imageData.SetDimensions(self.imR, self.imC, Xtmp.shape[0])
			imageData.GetPointData().AddArray(Xvtk)
			imageData.GetPointData().SetActiveScalars('Intensity')

			return imageData

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetDetailImages(self, data_id):
		"""Returns a vtkImageData of all detail images up the tree for a given data ID.
		Right now this is the same as the Wavelet Basis Images for each node up the tree.
		OLD fixed-dim version returned a list of one imagedata for each dimension.
		NEW variable-dim version returns a list of one imagedata for each _scale_"""

		if self.data_loaded:

			leafNode = self.IniLabels[data_id]
			chain = self.find_path_down_the_tree(leafNode)
			chain.reverse()		# change order so index will correspond to scale

			images_list = []
			# Need to separate out images for each dimension
			for node_id in chain:
				# Already a switch for wavelet or scaling functions in GetNodeBasisImages()
				images_list.append(self.GetNodeBasisImages(node_id))

			return images_list

		else:
			raise IOError, "Can't get image until data is loaded successfully"

	# ---------------------------------------
	def GetDetailWeights(self, data_id):
		"""Returns a list of arrays corresponding to the weights associated with
		the detail images for this particular data ID. (This is equivalent to the
		magnitudes of the wavelet or scaling function coefficients,
		reordered like the detail image stacks (one list item array for each scale)."""

		if self.data_loaded:

			leaf_node = self.IniLabels[data_id]
			row = N.nonzero(self.PointsInNet[leaf_node]==data_id)[0][0]	# final zero turns array->scalar
			mapped_node_idx = self.LeafNodesImap[leaf_node]

			# Skip any empty arrays
			wav_row_tuple = tuple(arr[row,:] for arr in self.CelCoeffs[mapped_node_idx,:] if arr.size != 0)
			wav_row = N.concatenate(wav_row_tuple, axis=1)

			# Right now doing fractional magnitudes only relative to this row's (data point's) values
			rel_mag_row_list = [N.abs(arr)/N.max(N.abs(wav_row)) for arr in wav_row_tuple]

			return rel_mag_row_list

		else:
			raise IOError, "Can't get image until data is loaded successfully"


	# ---------------------------------------
	def GetCategoryLUT(self, idx = 0):
		"""Returns a LUT for category coloring. Result depends on number
		of categories. Pass an index to get a proper label map (otherwise
		it defaults to the first (zeroth) label."""

		if self.data_loaded:

			cl = []
			lut = vtk.vtkLookupTable()
			
			if self.cat_labels_exist:
				num_labels = len(N.unique(self.cat_labels[idx,:]))
				
				if num_labels > 8 and num_labels <= 13:
					cl.append([float(cc)/255.0 for cc in [228, 26, 28]])  # Colorbrewer Set2 modY+4
					cl.append([float(cc)/255.0 for cc in [55, 126, 184]])
					cl.append([float(cc)/255.0 for cc in [77, 175, 74]])
					cl.append([float(cc)/255.0 for cc in [152, 78, 163]])
					cl.append([float(cc)/255.0 for cc in [255, 127, 0]])
					cl.append([float(cc)/255.0 for cc in [245, 193, 61]])
					cl.append([float(cc)/255.0 for cc in [166, 86, 40]])
					cl.append([float(cc)/255.0 for cc in [247, 129, 191]])
					cl.append([float(cc)/255.0 for cc in [153, 153, 153]])
					cl.append([float(cc)/255.0 for cc in [143, 0, 0]])
					cl.append([float(cc)/255.0 for cc in [22, 65, 110]])
					cl.append([float(cc)/255.0 for cc in [40, 115, 33]])
					cl.append([float(cc)/255.0 for cc in [61, 61, 61]])

					lutNum = len(cl)
					lut.SetNumberOfTableValues(lutNum)
					lut.Build()
					for ii,cc in enumerate(cl):
						lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
					lut.SetRange(0,len(cl)-1)
					lut.Build()
					return lut
									
				if num_labels <= 8:
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
			
					lutNum = len(cl)
					lut.SetNumberOfTableValues(lutNum)
					lut.Build()
					for ii,cc in enumerate(cl):
						lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
					lut.SetRange(0,len(cl)-1)
					lut.Build()
					return lut

				if num_labels > 13:
					
					lut.SetNumberOfTableValues(num_labels)
					lut.SetValueRange(0.8, 0.8)
					lut.SetSaturationRange(0.8, 0.8)
					lut.SetHueRange(0.0,1.0)
					lut.SetRampToLinear()
					lut.SetRange(0,num_labels-1)
					lut.Build()
					return lut
			
			lut.SetNumberOfTableValues(256)
			lut.SetValueRange(1.0, 1.0)
			lut.SetSaturationRange(1.0, 1.0)
			lut.SetHueRange(0.0,1.0)
			lut.SetRampToLinear()
			lut.SetRange(0,10)
			lut.Build()
			return lut
	
	# ---------------------------------------
	def GetDivergingLUT(self, map_colors = 'BrBg'):
		"""Returns a LUT for detail/axis images Default is ColorBrewer BrBg7
		(map_colors = 'BrBg'). Other choice is 'BR', which is a cheat on a blue-red
		binary map which isn't really diverging yet..."""

		if self.data_loaded:

			lut = vtk.vtkLookupTable()
		
			if map_colors == 'BR':
				
				lut.SetHueRange(0, 0.66)
				lut.SetValueRange(0.7, 0.7)
				lut.SetSaturationRange(1, 1)
				lut.Build()
				
			else:
			
				lutNum = 256
				lut.SetNumberOfTableValues(lutNum)
				ctf = vtk.vtkColorTransferFunction()
				ctf.SetColorSpaceToDiverging()
	
				cl = []
				cl.append([float(cc)/255.0 for cc in [140, 81, 10]])	# Colorbrewer BrBG 7
				cl.append([float(cc)/255.0 for cc in [216, 179, 101]])
				cl.append([float(cc)/255.0 for cc in [246, 232, 195]])
				cl.append([float(cc)/255.0 for cc in [245, 245, 245]])
				cl.append([float(cc)/255.0 for cc in [199, 234, 229]])
				cl.append([float(cc)/255.0 for cc in [90, 180, 172]])
				cl.append([float(cc)/255.0 for cc in [1, 102, 94]])
				vv = [float(xx)/float(len(cl)-1) for xx in range(len(cl))]
				vv.reverse()
				for pt,color in zip(vv,cl):
					ctf.AddRGBPoint(pt, color[0], color[1], color[2])
				
				for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
					cc = ctf.GetColor(ss)
					lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
				lut.SetRange(-1024,1024)
			
			return lut

	# ---------------------------------------
	def GetGrayscaleLUT(self, map_colors = 'gray'):
		"""Returns a linear LUT center and projected images. Grayscale default."""

		if self.data_loaded:

			lut = vtk.vtkLookupTable()
		
			if map_colors == 'R':
				
				lut.SetHueRange(0,0)
				lut.SetValueRange(0,1)
				lut.SetSaturationRange(1,1)
				lut.SetRampToLinear()
				lut.Build()
						
			else:
			
				lut.SetHueRange(0,0)
				lut.SetValueRange(0,1)
				lut.SetSaturationRange(0,0)
				lut.SetRampToLinear()
				lut.Build()
			
			return lut

# ==================
# Internal utility methods

	def find_path_down_the_tree(self, leafNodeID):
		""" Internal method returning chain: vector of the path down the tree,
		the first element of chain is the root and the last element of chain
		is the current node n"""

		# chain = [];
		# while n>0
		#    chain = [chain n];
		#    n = cp(n);
		# end

		n = leafNodeID
		chain = []

		while n >= 0:
			chain.append(n)
			n = self.cp[n]

		return chain

	def get_offspring(self, cp, node_id):
		"""Internal method finds all the offspring of (but not including)
		the given node in the tree cp."""

		# function offspring = get_offspring(cp, node)
		# offspring = [];
		# currentNodes = node;
		#
		# while ~isempty(currentNodes)
		#     newNodes = []; % collects all the children of currentNodes
		#     for i = 1:length(currentNodes)
		#         children = find(cp == currentNodes(i));
		#         newNodes = [newNodes children];
		#     offspring = [offspring newNodes];
		#     currentNodes = newNodes;

		offspring = N.array([],dtype='int32')
		current_nodes = N.array([node_id],dtype='int32')

		while (current_nodes.size > 0):
			new_nodes = N.array([],dtype='int32')		# collects children of current_nodes
			for node in current_nodes:
				children = N.nonzero(cp == node)[0]
				new_nodes = N.concatenate((new_nodes,children))
			offspring = N.concatenate((offspring,new_nodes))
			current_nodes = new_nodes

		return offspring

	def get_offspring2(self, cp, node_id):
		"""Returns a nested list of all offspring. Can be "flattened" by list(xflatten(result))"""
		children = N.nonzero(cp == node_id)[0]

		if children.size > 0:
			return [children.tolist(),[get_offspring2(cp,child) for child in children]]
		else:
			return []

	def has_children(self, cp, node_id):
		if N.nonzero(cp == node_id)[0].size > 0:
			return True
		else:
			return False

	def get_leaf_children(self, cp, node_id):
		"""This reaturns all leaf nodes that are descendants of a given node.
		It is a generator expression, so call: list(get_leaf_children(cp, node_id))"""
		children = N.nonzero(cp == node_id)[0]

		# If originally supplied with a leaf, then return it
		if len(children)==0:
			yield node_id
		else:
			for child_id in children:
				for sub_child in self.get_leaf_children(cp, child_id):
					yield sub_child

	def xflatten(self, seq):
		"""a generator to flatten a nested list
		from vegaseat on http://www.daniweb.com/forums/thread66694.html
		usage: flat_list = list(xflatten(nested_list))
		"""
		for x in seq:
			if type(x) is list:
				for y in xflatten(x):
					yield y
			else:
					yield x
