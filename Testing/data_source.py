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
		self.gW = MatInput['gW']
		
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
		self.Data = MatInput['Data']
		
		# X has already been projected to D dim by PCA
		self.X = N.mat(MatInput['X'])
		
		# X0 is original data
		self.X0 = N.mat(MatInput['X0'])
		self.cm = self.X0.mean(0)
		
		# Right now Matlab data doesn't have any record of original image dimensions
		# NOTE: Hard coding shape for now!
		self.imR = 28	# rows
		self.imC = 28	# cols
		
		# WC = MatInput['WavCoeffs']
		# Instead of using WC, which has already been ordered within Matlab
		self.WavCoeffsOrig = self.Data['MatWavCoeffs'][0,0].T
		
		# Structures are accessed like this:
		
		# PROJ.shape = (1000,120,8) for 1000 pts, 120 PCs, 8 scales -- only plot first few PCs
		self.PROJ = self.Data['Projections'][0,0]
		
		# NumPts = Number of data points (here number of individual images)
		self.NumPts = self.PROJ.shape[0]
		
		# D = dimensionality of original space before transformation
		#   (The original dimensionality of the data before SVD is higher than this)
		self.D = self.PROJ.shape[1]
		
		# Manifold dimensionality
		self.ManifoldDim = self.gW['ManifoldDimension'][0,0][0,0]
		
		# V = SVD result (X = (X0-mean)*V[:,:D])
		self.V = N.mat(MatInput['V'])
		
		# Hopefully saving memory by deleting this, but probably doesn't matter since it's
		# a local variable to the constructor...
		# del MatInput
		
		# CP is the child-parent array that defines the tree: CP[i]=(i's parent node)
		# Without cast CP ends up as 'uint16' or 'uint8' depending on max value
		self.CP = (self.gW['cp'][0,0][0].astype('int16') - 1)		# change here to zero-based indexing
		
		# IniLabels holds the node ID for each data point (change to zero-based)
		self.IniLabels = self.gW['IniLabels'][0,0][0] - 1
		
		# Making a list of numpy arrays because cell array structure is a pain...
		self.PIN = []	# Points In Net
		self.NIN = []	# Number In Net
		self.SCFUNS = []	# Scaling functions
		for ii in range(self.gW['PointsInNet'][0,0].shape[1]):
			self.PIN.append(self.gW['PointsInNet'][0,0][0,ii][0]-1)	# 0-based indices
			self.NIN.append(self.gW['PointsInNet'][0,0][0,ii][0].size)
			self.SCFUNS.append(N.mat(self.gW['ScalFuns'][0,0][0,ii]))	# matrix
		
		# Scale of each node
		self.Scales = self.gW['Scales'][0,0][0] - 1	# zero-based
		
		# J = Total number of scales
		self.J = self.Scales.max() + 1				# compensate for zero-based
		
		# nAllNets = Number of nodes in the METIS tree
		# nAllNets = length(gW.cp); 
		
		# Note that not all leaf nodes have scale J due to the possible 
		#   incompleteness of the metis tree
		# isaleaf = ones(1,nAllNets); 
		# isaleaf(gW.cp(gW.cp>0)) = 0; 
		# leafNodes = find(isaleaf>0); 

		self.data_loaded = True

	def GetTree(self):
		"""Returns a full vtkTree based on data loaded in LoadData()."""
		
		if self.data_loaded:
			
			vertex_id = vtk.vtkIdTypeArray()
			vertex_id.SetName('vertex_ids')
			for ii in range(self.gW['PointsInNet'][0,0].shape[1]):
				vertex_id.InsertNextValue(ii)
	
			# NINarray = N.array(NIN)
			# NINvtk = VN.numpy_to_vtk(NINarray)
			# SCALESvtk = VN.numpy_to_vtk(gW['Scales'][0,0][0])
			
			# Trying to avoid numpy_to_vtk for now... had a few troubles once...
			NINvtk = vtk.vtkIntArray()
			NINvtk.SetNumberOfComponents(1)
			NINvtk.SetNumberOfTuples(len(self.NIN))
			NINvtk.SetName('num_in_vertex')
			SCALESvtk = vtk.vtkIntArray()
			SCALESvtk.SetNumberOfComponents(1)
			SCALESvtk.SetNumberOfTuples(len(self.NIN))
			SCALESvtk.SetName('scale')
			BLANKvtk = vtk.vtkStringArray()
			BLANKvtk.SetNumberOfComponents(1)
			BLANKvtk.SetNumberOfTuples(len(self.NIN))
			BLANKvtk.SetName('blank')
			for ii in range(len(self.NIN)):
				NINvtk.SetTuple1(ii,self.NIN[ii])
				SCALESvtk.SetTuple1(ii,self.gW['Scales'][0,0][0][ii])
				BLANKvtk.SetValue(ii,"")
			
			# Build tree out of CP list of "is a child of"
			#	remembering that Matlab indices are 1-based and numpy/VTK 0-based
			print 'Building graph'
			dg = vtk.vtkMutableDirectedGraph()
			edge_id = vtk.vtkIdTypeArray()
			edge_id.SetName('edge_ids')
			for ii in range(self.CP.size):
				dg.AddVertex()
			for ii in range(self.CP.size):
				if self.CP[ii] > 0:		# CP already zero-based
					dg.AddGraphEdge(self.CP[ii],ii)		# Method for use with wrappers -- AddEdge() in C++
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
		
	def GetWaveletCoeffImage(self, XOrderedLeafIds=None ):
		"""Returns a vtkImageData 2D image with the wavelet coefficients at all dimensions
		and scales. If you give an ordered list of the leaf node IDs, then the matrix will
		be sorted accordingly so the image should be correct for the icicle view."""
		
		if self.data_loaded:
			
			# If the caller wants the wavelet coeffs returned in a certain order, then
			# they need to supply a sorted list of the LeafIds
			if XOrderedLeafIds is not None:
				# Create an array holding the indices of the leaf vertices in the proper order
				SortedLeafIdxArray = N.array([],dtype='uint16')
				for ii in range(XOrderedLeafIds.size):
					SortedLeafIdxArray = N.concatenate((SortedLeafIdxArray,self.PIN[XOrderedLeafIds[ii]]))
					
				# And then reorder the Wavelet Coefficients matrix according to this
				WCsorted = self.WavCoeffsOrig[:,SortedLeafIdxArray]
			else:
				WCsorted = self.WavCoeffsOrig
				
			# Create vtkImageData out of WavCoeffs for texturing icicle view tree
			# WCvtk = VN.numpy_to_vtk(WC.ravel()[::-1].copy())
			WCvtk = vtk.vtkDoubleArray()
			WCvtk.SetNumberOfComponents(1)
			WCvtk.SetNumberOfTuples(WCsorted.size)
			ind = 0
			WCtmp = WCsorted[::-1,:].copy()		# not really sure if need copy() here...
			for ii in range(WCtmp.shape[0]):
				for jj in range(WCtmp.shape[1]):
					zz = WCtmp[ii,jj].copy()	# or here...
					WCvtk.SetTuple1(ind, zz)
					ind += 1
			WCvtk.SetName('Coeffs')
			
			WCimageData = vtk.vtkImageData()
			WCimageData.SetOrigin(0,0,0)
			WCimageData.SetSpacing(1,1,1)
			WCimageData.SetDimensions(self.WavCoeffsOrig.shape[1],self.WavCoeffsOrig.shape[0],1)
			WCimageData.GetPointData().AddArray(WCvtk)
			WCimageData.GetPointData().SetActiveScalars('Coeffs')
			
			return WCimageData
		
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetNodeAllScaleCoeffTable(self, nodeID):
		"""Returns a table of all the wavelet coefficients for a single tree
		icicle view) node at all scales for plotting on a parallel coodinates plot."""
				
		if self.data_loaded:
			
			# For a given nodeID, get PIN and then extract all coeffs at every scale
			# Columns of table will be rows of the WavCoeffsOrig matrix
			
			IDarray = self.PIN[nodeID]
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

	def GetNodeOneScaleCoeffTable(self, nodeID):
		"""Returns a table of the wavelet coefficients at a single node at a single
		scale for plotting on a scatter plot."""
		
		if self.data_loaded:
			
			# For a given nodeID, get PIN and then extract only coeffs at correct scale
			# Columns of table will be rows of the WavCoeffsOrig matrix
			
			IDarray = self.PIN[nodeID]
			scale = self.Scales[nodeID]
			
			rows = N.arange(self.WavCoeffsOrig.shape[0])
			scales = (rows/self.ManifoldDim)
			rems = N.mod(rows, self.ManifoldDim)
			
			table = vtk.vtkTable()
			for ii in rows:
				# Only pull out correct scale
				if scales[ii] == scale:
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

	def GetNodeBasisImages(self, nodeID):
		"""Returns a vtkImageData of all basis images for a given node."""
		
		if self.data_loaded:
			
			# %% Display all detail coordinates for a given leaf node
			# 
			# imagesc(reshape(V(:,1:D)*gW.ScalFuns{i}, self.imR,[])); 
			#
			# Need to create separate images (Z) for each column of matrix result

			return
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetProjectedImages(self, IDlist):
		"""Given a list of IDs selected from a parallel coordinates plot, returns
		a vtkImageData with all of the projected (reduced dimensionality) images
		for those IDs."""
		
		if self.data_loaded:
			
			# %% Add original but projected
			# i = sample; % 39 when digit =1
			# 
			# X_orig = X*V(:,1:D)'+repmat(cm, size(X,1),1);
			# imagesc(reshape(X_orig(i,:), self.imR,[]));

			Xtmp = self.X[IDlist,:]*self.V[:,:self.D].T
			X_orig = Xtmp + N.tile(self.cm,(Xtmp.shape[0],1))	# tile ~ repmat
			
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
		that ID. (Not sure whether will also return original...)"""
		
		if self.data_loaded:
			
			
			# %% Calculate and display approximations at all scales
			# 
			# i = sample; % 39 when digit =1
			# 
			# for j = J:-1:1
			# 	X_approx = Data.Projections(:,:,j);
			# 	X_img = X_approx*V(:,1:D)'+repmat(cm, size(X_approx,1),1);
			# 	imagesc(reshape(X_img(i,:), self.imR,[]));
			# end
			# 

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

	def GetDetailImages(self, ID):
		"""Returns a vtkImageData of all detail images up the tree for a given data ID.
		Right now this is the same as the Basis Images for each node up the tree.
		Returning list of vtkImageData, one for each dimension."""
		
		if self.data_loaded:
			
			# %% Display all detail coordinates for a given leaf node
			# 
			# leafNode = gW.IniLabels(i);
			# chain = find_path_down_the_tree(gW.cp, leafNode);
			# 
			# for i = 1:length(chain)
			#     imagesc(reshape(V(:,1:D)*gW.ScalFuns{chain(i)}, self.imR,[])); 
			# end
			
			leafNode = self.IniLabels[ID]
			chain = self.find_path_down_the_tree(leafNode)
			chain.reverse()		# change order so index will correspond to scale
			
			images_list = []
			# Need to separate out images for each dimension
			for dd in range(self.ManifoldDim):
				# Concatenate scaling functions for all nodes in chain
				ScFuncts = N.concatenate([self.SCFUNS[node_id][:,dd] for node_id in chain],1)
				# Compute all detail images for that dimension
				image_cols = self.V[:,:self.D]*ScFuncts
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
				
				images_list.append(imageData)

			return images_list
			
		else:
			raise IOError, "Can't get image until data is loaded successfully"

	def GetDetailWeights(self, ID):
		"""Returns an array corresponding to the weights associated with
		the detail images for this particular data ID. (This is equivalent to the
		magnitudes of the wavelet coefficients, reordered like the detail image stacks."""
		
		if self.data_loaded:
			
			wav_col = self.WavCoeffsOrig[:,ID]
						
			# Right now doing fractional magnitudes only relative to this column's values
			rel_mag_col = N.abs(wav_col)/N.max(N.abs(wav_col))
			wav_tmp = rel_mag_col.reshape(-1,self.ManifoldDim)
			
			return wav_tmp
			
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
			n = self.CP[n]
			
		return chain