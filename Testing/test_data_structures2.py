# With this version trying to set a separate actor in front of the icicle view which
# contains a "non-pickable" version of the textured icicle polys

import vtk
import scipy.io
import vtk.util.numpy_support as VN
import numpy as N

data_dir = '/Users/emonson/Data/Fodava/EMoGWDataSets/'
data_file = data_dir + 'mnist12_1k_20100825.mat'


print 'Trying to load data set from .mat file...'

if len(data_file) == 0:
	raise IOError, "No data file name: Use SetFileName('file.mat') before LoadData()"

try:
	MatInput = scipy.io.loadmat(data_file, struct_as_record=True)
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

# gW = MatInput['gW']

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

# Data = MatInput['Data']

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

GWTopts = MatInput['GWTopts']

# X has already been projected to D dim by PCA
X = N.mat(MatInput['X'])
cm = N.mat(MatInput['cm'])	# not sure if should be matrix or array...

# Various plain matrices
# NOTE: Have to be careful of anything that can have a 0 value in Matlab
# because it might be naturally imported as an unsigned int, and then
# when you subtract 1 from it you don't get a negative number as you'd expect
V = N.mat(MatInput['V'])
cp = (MatInput['cp'][0].astype('int16') - 1)	# change here to zero-based indexing
IniLabels = (MatInput['IniLabels'][0] - 1)		# change here to zero-based indexing
NumberInNet = MatInput['NumberInNet'][0]
Scales = (MatInput['Scales'][0] - 1)					# zero-based
IsALeaf = MatInput['IsALeaf'][0].astype('bool')
LeafNodes = (MatInput['LeafNodes'][0] - 1)		# zero-based
LeafNodesImap = (MatInput['LeafNodesImap'][0].astype('int16') - 1)		# zero-based
CelWavCoeffs = MatInput['CelWavCoeffs']
CelScalCoeffs = MatInput['CelScalCoeffs']

# NOTE: gW and Data are class numpy.ndarray
#		MatInput is just a dict, so can directly look for variables there

if MatInput.has_key('imR') and MatInput.has_key('imC'):
	print 'Grabbing image dimensions from matlab file'
	imR = MatInput['imR']
	imC = MatInput['imC']
else:
	# Right now Matlab data doesn't have any record of original image dimensions
	# NOTE: Hard coding shape for now!
	print 'Hacking image dimensions from file name'
	if (data_file.find('mnist') >= 0):
		imR = 28	# rows
		imC = 28	# cols
	elif (data_file.find('frey') >= 0):
		imR = 20	# rows
		imC = 28	# cols
	elif (data_file.find('olivetti') >= 0):
		imR = 64	# rows
		imC = 64	# cols
	else:
		imR = 20
		imC = 20
		print 'Could not find matching file name -- probably wrong image dimensions!!!!'

# NumPts = Number of data points (here number of individual images)
NumPts = IniLabels.shape[0]
AmbientDimension = MatInput['AmbientDimension'][0,0]	# used to call this D
		
# Manifold dimensionality variable now, not fixed...
# ManifoldDim = gW['ManifoldDimension'][0,0][0,0]

# V = SVD result (X = (X0-mean)*V[:,:D])
# This V already is only V(:,1:D) so we can use it whole
V = N.mat(MatInput['V'])

#     PointsInNet = gW.PointsInNet;
#     ScalFuns = gW.ScalFuns;
#     WavBases = gW.WavBases;
#     Centers = gW.Centers;

# Converting cell arrays to lists of numpy arrays
PointsInNet = []	# Points In Net
ScalFuns = []	# Scaling functions
WavBases = []	# Wavelet bases
Centers = []	# Center of each node
NodeWavCoeffs = []
NodeScalCoeffs = []
for ii in range(MatInput['PointsInNet'].shape[1]):
	PointsInNet.append(MatInput['PointsInNet'][0,ii][0]-1)	# 0-based indices
	ScalFuns.append(N.mat(MatInput['ScalFuns'][0,ii]))			# matrix
	WavBases.append(N.mat(MatInput['WavBases'][0,ii]))			# matrix
	Centers.append(N.mat(MatInput['Centers'][0,ii][0])) 		# matrix
	NodeWavCoeffs.append(N.mat(MatInput['NodeWavCoeffs'][0,ii])) 		# matrix
	NodeScalCoeffs.append(N.mat(MatInput['NodeScalCoeffs'][0,ii])) 		# matrix

# J = Total number of scales
J = Scales.max()

def get_offspring(cp, node_id):
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
			
def get_offspring2(cp, node_id):
	"""Returns a nested list of all offspring. Can be "flattened" by xflatten(result)"""
	children = N.nonzero(cp == node_id)[0]
	
	if children.size > 0:
		return [children.tolist(),[get_offspring2(cp,child) for child in children]]
	else:
		return []
		
def has_children(cp, node_id):
	if N.nonzero(cp == node_id)[0].size > 0:
		return True
	else:
		return False

def get_leaf_children(cp, node_id):
	"""This reaturns all leaf nodes that are descendants of a given node"""
	children = N.nonzero(cp == node_id)[0]
	
	if len(children)==0:
		yield node_id
	else:
		for x in children:
			for y in get_leaf_children(cp, x):
				yield y
		
def xflatten(seq):
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
				