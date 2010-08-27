# With this version trying to set a separate actor in front of the icicle view which
# contains a "non-pickable" version of the textured icicle polys

import vtk
import scipy.io
import vtk.util.numpy_support as VN
import numpy as N
import os


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
		print 'Trying to load data set from .mat file...'
		
		if len(self.data_file) == 0:
			raise IOError, "No data file name: Use SetFileName('file.mat') before LoadData()"
		
		try:
			MatInput = scipy.io.loadmat(self.data_file, struct_as_record=True)
		except:
			raise IOError, "Can't load supplied matlab file"
			# return
			
		# Get variables out of Matlab structure
		print 'Transferring variables from Matlab structures'
		
		# gW = 
		#                   knn: 20
		#           knnAutotune: 10
		#      SmallestMetisNet: 10
		#            threshold0: 0
		#            threshold1: 0.0100
		#            threshold2: 0.0100
		#     ManifoldDimension: 4
		#   gW: a structure with the following fields:
		#       .cp: vector encoding the metis tree structure with .cp(x) = the parent of x.
		#             Each node is a subset of points and a parent is the union of its children
		#       .IniLabels: the labels of the data points w.r.t. the leaf nodes in the
		#                  METIS tree (all the leaves are a full partition of the data)
		#       .PointsInNet: cell array, PointsInNet{i} contains all the points that
		#                           is in the node i of the metis tree
		#       .Sizes: vector of number of points in each node
		#       .Radii: vector of the radius of each node
		#       .Scales: vector of scales of all the nodes in the metis tree;
		#                    the root has scale 1 and a larger scale implies a finer
		#                    approximation
		#       .Centers: cell array of the centers of the nodes
		#       .ScalFuns: cell array of the local bases of the nodes
		#       .Sigmas: cell array of the local singular values 
		#       .WavBasis: cell array of wavelet bases; each cell encodes the basis vectors that are
		#                        present in the current node but orthogonal to the
		#                        parent. For simplicity, the wavelet basis at the root
		#                        is identified with the scaling functions at the root
		#       .WavConsts: cell array of translations that are needed to move to a
		#                           node from its parent                       
		#       .WavSingVals: cell array of corresponding singular values associated
		#                             with the wavelet bases. At the root, it coincides
		#                             with .Sigmas
		
		# NEW gW version in pruning code (8/25/2010)
		# 		               X: [2000x120 double]
		#              X_clean: [2000x120 double]
		#                   X0: [2000x784 double]
		#                 opts: [1x1 struct]
		#                   cp: [1x103 double]
		#               Scales: [1x103 double]
		#              isaleaf: [1x103 double]
		#            LeafNodes: [1x52 double]
		#            IniLabels: [1x2000 double]
		#          PointsInNet: {1x103 cell}
		#                Radii: [1x103 double]
		#                Sizes: [1x103 double]
		#              Centers: {1x103 cell}
		#             ScalFuns: {1x103 cell}
		#               Sigmas: {1x103 cell}
		#             WavBases: {1x103 cell}
		#            WavConsts: {1x103 cell}
		#          WavSingVals: {1x103 cell}
		#     epsEncodingCosts: [1x103 double]
		#              cp_orig: [1x103 double]
		#            DictCosts: 420240
		
		# self.gW = MatInput['gW']
		
		#   Data: a structure of the following fields:
		#       .ScalCoeffs: N by k*J matrix of coefficients relative to the
		#                    scaling functions at the local nets
		#       .Projections: N-by-D-by-J array of projections of X onto the 
		#                     local tangent planes at all J scales
		#       .MatWavCoeffs: N by k*J matrix of wavelet coefficients relative to the
		#                      wavelet bases at the local nets
		#       .CelWavCoeffs: N by J cell array of wavelet coefficients relative to the
		#                      wavelet bases at the local nets
		#       .Wavelets: N-by-D-by-J array of differences between the projections 
		#                  at consecutive scales
		#        (In the above, J is the number of scales, and k is the manifold
		#         dimension)
		
		# NEW Data version in pruning code (8/25/2010)
		# 
		#              CelScalCoeffs: {52x8 cell}
		#                Projections: [2000x120x8 double]
		#               CelWavCoeffs: {52x8 cell}
		#                   Wavelets: [2000x120x8 double]
		#      TangentialCorrections: [2000x120x8 double]
		#     ProjectedInitialErrors: [2000x120x8 double]
		#               MatWavCoeffs: [2000x187 double]
		#                 maxWavDims: [6 6 7 8 7 42 64 47]
		#                 MatWavDims: [2000x8 double]
		#                CoeffsCosts: 60208
		#                      Xmean: [1x784 double]
		#                          V: [784x784 double]

		# self.Data = MatInput['Data']
		
		# 		GWTopts = 
		# 
		#             AmbientDimension: 120
		#                          knn: 50
		#                  knnAutotune: 30
		#             smallestMetisNet: 20
		#            ManifoldDimension: 0
		#                    errorType: 'relative'
		#                   threshold0: 0.5000
		#                    precision: 0.0100
		#                   threshold1: 0.1000
		#                   threshold2: 0.0100
		#                      pruning: 1
		#     addTangentialCorrections: 1
		#                  sparsifying: 0
		#                    splitting: 0
		
		# NEW: Poofed out variables
		#
		#     % "poofing" out variables for easier loading in Python code 
		#     % and for file size since don't need nearly all of this stuff...
		#     % I know it makes it less flexible later...
		# 
		#     AmbientDimension = GWTopts.AmbientDimension;
		#     X = gW.X;
		#     % cm
		#     % imR
		#     % imC
		#     
		#     % Redundant for now...
		#     CelWavCoeffs = Data.CelWavCoeffs;
		#     
		#     num_nodes = length(gW.cp);
		#     LeafNodesImap(gW.LeafNodes) = 1:length(gW.LeafNodes);
		#     NodeWavCoeffs = cell(1,num_nodes);
		# 
		#     for node_idx = 1:num_nodes,
		#         offspring = [node_idx get_offspring(gW.cp, node_idx)];
		#         relevantLeafNodes = offspring(logical(gW.isaleaf(offspring)));
		#         NodeWavCoeffs{node_idx} = cat(1, Data.CelWavCoeffs{LeafNodesImap(relevantLeafNodes), gW.Scales(node_idx)});
		#     end
		#     
		#     CelScalCoeffs = Data.CelScalCoeffs;
		#     NodeScalCoeffs = cell(1,num_nodes);
		# 
		#     for node_idx = 1:num_nodes,
		#         offspring = [node_idx get_offspring(gW.cp, node_idx)];
		#         relevantLeafNodes = offspring(logical(gW.isaleaf(offspring)));
		#         NodeScalCoeffs{node_idx} = cat(1, Data.CelWavCoeffs{LeafNodesImap(relevantLeafNodes), gW.Scales(node_idx)});
		#     end
		#     
		#     % Should calculate Projections rather than storing -- it's big...
		#     
		#     % node_idx = leafNodes(leaf_node_idx);
		#     % data_idxs = find(gW.IniLabels == node_idx); % same as PointsInNet{net}
		#     % nPts = length(data_idxs);
		#     % j_max = gW.Scales(node_idx);
		#     % gWCentersnet = repmat(gW.Centers{node_idx},nPts,1);
		#     % Data.Projections(data_idxs,:,j_max) = Data.CelScalCoeffs{i,j_max}*gW.ScalFuns{node_idx}' + gWCentersnet;
		#     % X_approx = Data.Projections(:,:,scale);
		#     % X_img = X_approx*V(:,1:GWTopts.AmbientDimension)'+repmat(cm, size(X_approx,1),1);
		# 
		#     % Projections = Data.Projections;
		#     
		#     % May be able to get away with only saving
		#     % V(:,1:GWTopts.AmbientDimension)
		#     V = Data.V(:,1:GWTopts.AmbientDimension);
		#     
		#     cp = gW.cp;
		#     IniLabels = gW.IniLabels;
		#     PointsInNet = gW.PointsInNet;
		#     NumberInNet = gW.Sizes;
		#     ScalFuns = gW.ScalFuns;
		#     WavBases = gW.WavBases;
		#     Centers = gW.Centers;
		#     Scales = gW.Scales;
		#     IsALeaf = gW.isaleaf;
		#     LeafNodes = gW.LeafNodes;
		
		self.GWTopts = MatInput['GWTopts']

		# X has already been projected to D dim by PCA
		self.X = N.mat(MatInput['X'])
		self.cm = N.mat(MatInput['cm'])	# not sure if should be matrix or array...
		
		# Various plain matrices
		# NOTE: Have to be careful of anything that can have a 0 value in Matlab
		# because it might be naturally imported as an unsigned int, and then
		# when you subtract 1 from it you don't get a negative number as you'd expect
		self.V = N.mat(MatInput['V'])
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
		
		# NumPts = Number of data points (here number of individual images)
		self.NumPts = self.IniLabels.shape[0]
		self.AmbientDimension = MatInput['AmbientDimension'][0,0]	# used to call this D
				
		# Manifold dimensionality variable now, not fixed...
		# self.ManifoldDim = self.gW['ManifoldDimension'][0,0][0,0]
		
		# V = SVD result (X = (X0-mean)*V[:,:D])
		# This V already is only V(:,1:D) so we can use it whole
		self.V = N.mat(MatInput['V'])
				
		#     PointsInNet = gW.PointsInNet;
		#     ScalFuns = gW.ScalFuns;
		#     WavBases = gW.WavBases;
		#     Centers = gW.Centers;

		# Converting cell arrays to lists of numpy arrays
		self.PointsInNet = []	# Points In Net
		self.ScalFuns = []	# Scaling functions
		self.WavBases = []	# Wavelet bases
		self.Centers = []	# Center of each node
		# self.NodeWavCoeffs = []
		# self.NodeScalCoeffs = []
		for ii in range(MatInput['PointsInNet'].shape[1]):
			self.PointsInNet.append(MatInput['PointsInNet'][0,ii][0]-1)	# 0-based indices
			self.ScalFuns.append(N.mat(MatInput['ScalFuns'][0,ii]))			# matrix
			self.WavBases.append(N.mat(MatInput['WavBases'][0,ii]))			# matrix
			self.Centers.append(N.mat(MatInput['Centers'][0,ii][0])) 		# matrix
			# self.NodeWavCoeffs.append(N.mat(MatInput['NodeWavCoeffs'][0,ii])) 		# matrix
			# self.NodeScalCoeffs.append(N.mat(MatInput['NodeScalCoeffs'][0,ii])) 		# matrix
				
		# J = Total number of scales
		# self.J = self.Scales.max()
		
		# Creating a storage space for ordering of leaf nodes in icicle view with
		# default values of ordering according to Matlab-saved LeafNodes
		ice_leaf_ids = self.LeafNodes
		ice_leaf_xmins = N.arange(ice_leaf_ids.size)
		ice_ids_mapped = self.LeafNodesImap[ice_leaf_ids]
		self.mapped_leaf_pos = N.zeros_like(ice_leaf_xmins)
		self.mapped_leaf_pos[ice_ids_mapped] = ice_leaf_xmins

		self.data_loaded = True

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
		
	def GetWaveletCoeffImages(self, ice_leaf_ids=None, ice_leaf_xmins=None ):
		"""Returns a list of vtkImageData 2D image with the wavelet coefficients at all dimensions
		for all nodes. If you give the positions and IDs of the leaf nodes, as laid out by
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
					img_tuple = tuple(pp for pp in self.CelWavCoeffs[sorted_offspring_idxs, self.Scales[node_id]])
					# May need to transpose this...
					img = N.concatenate(img_tuple, axis=0).T
					
					# Create vtkImageData out of WavCoeffs for texturing icicle view tree
					# .copy() is to force the array to be contiguous for numpy_to_vtk
					# deep=True should keep reference around even after numpy array is destroyed
					
					WCvtk = VN.numpy_to_vtk(img.ravel()[::-1].copy(), deep=True)
					WCvtk.SetName('Coeffs')
					
					WCimageData = vtk.vtkImageData()
					WCimageData.SetOrigin(0,0,0)
					WCimageData.SetSpacing(1,1,1)
					WCimageData.SetDimensions(img.shape[0],img.shape[1],1)
					WCimageData.GetPointData().AddArray(WCvtk)
					WCimageData.GetPointData().SetActiveScalars('Coeffs')
					
					WC_imagedata_list.append(WCimageData)
					
				return WC_imagedata_list
		
		else:
			raise IOError, "Can't get image until data is loaded successfully"

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

	def GetNodeAllScaleCoeffTable(self, node_id):
		"""Returns a table of all the wavelet coefficients for a single tree
		icicle view) node at all scales for plotting on a parallel coodinates plot."""
				
		# * * * OLD FIXED-DIM VERSION * * *
		
		if self.data_loaded:
			
			# For a given node_id, get PIN and then extract all coeffs at every scale
			# Columns of table will be rows of the WavCoeffsOrig matrix
			
			IDarray = self.PIN[node_id]
			rows = N.arange(self.WavCoeffsOrig.shape[0])
			scales = (rows/self.ManifoldDim)
			rems = N.mod(rows, self.ManifoldDim)
			
			table = vtk.vtkTable()
			for ii in rows:
				column = VN.numpy_to_vtk(self.WavCoeffsOrig[ii,IDarray].T, deep=True)
				column.SetName(str(scales[ii]) + '.' + str(rems[ii])) 
				table.AddColumn(column)
			
			# Trying to set PedigreeIds to that parallel coords selections have correct IDs
			IDvtk = VN.numpy_to_vtk(IDarray, deep=True)
			IDvtk.SetName('pedigree_ids')
			table.AddColumn(IDvtk)
			table.GetRowData().SetActivePedigreeIds('pedigree_ids')
			
			return table
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetNodeOneScaleCoeffTable(self, node_id):
		"""Returns a table of the wavelet coefficients at a single node at a single
		scale for plotting on a scatter plot. Relying on icicle_view already having
		called GetWaveletCoeffImages() with correct positions of leaf nodes in view,
		otherwise just using original Matlab-saved LeafNodes ordering."""
				
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
			img_tuple = tuple(pp for pp in self.CelWavCoeffs[sorted_offspring_idxs, scale])
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
			
			return table
			
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetNodeWaveletImages(self, node_id):
		"""Returns a vtkImageData of all wavelet basis images for a given node."""
				
		if self.data_loaded:
			
			# %% Display all detail coordinates for a given leaf node
			# 
			# imagesc(reshape(V(:,1:D)*gW.ScalFuns{i}, self.imR,[])); 
			#
			# Need to create separate images (Z) for each column of matrix result

			# Compute all detail images for that node
			# Now V already chopped to AmbientDimension
			image_cols = self.V*self.WavBases[node_id]
			# To make it linear, it is the correct order (one image after another) to .ravel()
			images_linear = N.asarray(image_cols.T).ravel()
			
			intensity = VN.numpy_to_vtk(images_linear, deep=True)
			intensity.SetName('DiffIntensity')
	
			imageData = vtk.vtkImageData()
			imageData.SetOrigin(0,0,0)
			imageData.SetSpacing(1,1,1)
			imageData.SetDimensions(self.imR, self.imC, image_cols.shape[1])
			imageData.GetPointData().AddArray(intensity)
			imageData.GetPointData().SetActiveScalars('DiffIntensity')
			
			return imageData
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetNodeScalingFunctionImages(self, node_id):
		"""Returns a vtkImageData of all scaling function basis images for a given node."""
				
		if self.data_loaded:
			
			image_cols = self.V*self.ScalFuns[node_id]
			images_linear = N.asarray(image_cols.T).ravel()
			
			intensity = VN.numpy_to_vtk(images_linear, deep=True)
			intensity.SetName('DiffIntensity')
	
			imageData = vtk.vtkImageData()
			imageData.SetOrigin(0,0,0)
			imageData.SetSpacing(1,1,1)
			imageData.SetDimensions(self.imR, self.imC, image_cols.shape[1])
			imageData.GetPointData().AddArray(intensity)
			imageData.GetPointData().SetActiveScalars('DiffIntensity')
			
			return imageData
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetNodeCenterImage(self, node_id):
		"""Returns a vtkImageData of the center image for a given node."""
				
		if self.data_loaded:
			
			# imagesc(reshape(gW.Centers{1}*V(:,1:D)'+cm,28,[]))
			
			# Compute all detail images for that dimension
			image_col = self.Centers[node_id]*self.V.T + self.cm
			# To make it linear, it is the correct order (one image after another) to .ravel()
			image_linear = N.asarray(image_col.T).ravel()
			
			intensity = VN.numpy_to_vtk(image_linear, deep=True)
			intensity.SetName('Intensity')
	
			imageData = vtk.vtkImageData()
			imageData.SetOrigin(0,0,0)
			imageData.SetSpacing(1,1,1)
			imageData.SetDimensions(self.imR, self.imC, 1)
			imageData.GetPointData().AddArray(intensity)
			imageData.GetPointData().SetActiveScalars('Intensity')
			
			return imageData
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetProjectedImages(self, IDlist):
		"""Given a list of IDs selected from a parallel coordinates plot, returns
		a vtkImageData with all of the projected (reduced dimensionality by SVD) images
		for those IDs. (e.g. typically 120 dim rather than original 768 dim for MNIST digits)"""
				
		if self.data_loaded:
			
			# X_orig = X*V(:,1:GWTopts.AmbientDimension)'+repmat(cm, size(X,1),1);

			# V now already chopped to AmbientDimension
			Xtmp = self.X[IDlist,:]*self.V.T
						
			# numpy should automatically do tiling!!
			X_orig = Xtmp + self.cm
			# X_orig = Xtmp + N.tile(self.cm,(Xtmp.shape[0],1))	# tile ~ repmat
			
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
			imageData.SetDimensions(self.imR, self.imC, Xtmp.shape[0])
			imageData.GetPointData().AddArray(Xvtk)
			imageData.GetPointData().SetActiveScalars('Intensity')
			
			return imageData
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

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
				images_list.append(self.GetNodeWaveletImages(node_id))

			return images_list
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetDetailWeights(self, data_id):
		"""Returns a list of arrays corresponding to the weights associated with
		the detail images for this particular data ID. (This is equivalent to the
		magnitudes of the wavelet coefficients, reordered like the detail image stacks
		(one list item array for each scale)."""
				
		if self.data_loaded:
			
			leaf_node = self.IniLabels[data_id]
			row = N.nonzero(self.PointsInNet[leaf_node]==data_id)[0][0]	# final zero turns array->scalar
			mapped_node_idx = self.LeafNodesImap[leaf_node]
			
			# Skip any empty arrays
			wav_row_tuple = tuple(arr[row,:] for arr in self.CelWavCoeffs[mapped_node_idx,:] if arr.size != 0)
			wav_row = N.concatenate(wav_row_tuple, axis=1)
						
			# Right now doing fractional magnitudes only relative to this row's (data point's) values
			rel_mag_row_list = [N.abs(arr)/N.max(N.abs(wav_row)) for arr in wav_row_tuple]
			
			return rel_mag_row_list
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

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
					