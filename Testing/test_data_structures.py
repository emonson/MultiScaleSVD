# With this version trying to set a separate actor in front of the icicle view which
# contains a "non-pickable" version of the textured icicle polys

import vtk
import scipy.io
import vtk.util.numpy_support as VN
import numpy as N

SHRINK = 0.1
THICK = 1.0

# ----------
# Load and construct whole graph and multi-resolution data from Matlab structure
dataDir = '/Users/emonson/Data/Fodava/EMoGWDataSets/'
# filename = dataDir + 'mnist1_1k_20100324.mat'
filename = dataDir + 'mnist1_5c_20100324.mat'
# filename = dataDir + 'mnist1_5k_20100320.mat'

print 'Loading data set from .mat file...'
MatInput = scipy.io.loadmat(filename,struct_as_record=True)

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
gW = MatInput['gW']

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
Data = MatInput['Data']

# X has already been projected to D dim by PCA
X = N.mat(MatInput['X'])

# X0 is original data
X0 = N.mat(MatInput['X0'])
cm = X0.mean(0)

# WC = MatInput['WavCoeffs']
# Instead of using WC, which has already been ordered within Matlab
WavCoeffsOrig = Data['MatWavCoeffs'][0,0].T

# Structures are accessed like this:

# PROJ.shape = (1000,120,8) for 1000 pts, 120 PCs, 8 scales -- only plot first few PCs
PROJ = Data['Projections'][0,0]

# NumPts = Number of data points (here number of individual images)
NumPts = PROJ.shape[0]

# D = dimensionality of original space before transformation
D = PROJ.shape[1]

# V = SVD result (X = (X0-mean)*V[:,:D])
V = N.mat(MatInput['V'])

del MatInput

# CP is the child-parent array that defines the tree: CP[i]=(i's parent node)
# Without cast CP ends up as 'uint16' or 'uint8' depending on max value
CP = (gW['cp'][0,0][0].astype('int16') - 1)		# change here to zero-based indexing

# Making a list of numpy arrays because cell array structure is a pain...
PIN = []	# Points In Net
NIN = []	# Number In Net
SCFUNS = []	# Scaling functions
vertex_id = vtk.vtkIdTypeArray()
vertex_id.SetName('vertex_ids')
for ii in range(gW['PointsInNet'][0,0].shape[1]):
	PIN.append(gW['PointsInNet'][0,0][0,ii][0]-1)	# 0-based indices
	NIN.append(gW['PointsInNet'][0,0][0,ii][0].size)
	SCFUNS.append(gW['ScalFuns'][0,0][0,ii])
	vertex_id.InsertNextValue(ii)
	
# NINarray = N.array(NIN)
# NINvtk = VN.numpy_to_vtk(NINarray)
# SCALESvtk = VN.numpy_to_vtk(gW['Scales'][0,0][0])

# Trying to avoid numpy_to_vtk for now... had a few troubles once...
NINvtk = vtk.vtkIntArray()
NINvtk.SetNumberOfComponents(1)
NINvtk.SetNumberOfTuples(len(NIN))
NINvtk.SetName('num_in_vertex')
SCALESvtk = vtk.vtkIntArray()
SCALESvtk.SetNumberOfComponents(1)
SCALESvtk.SetNumberOfTuples(len(NIN))
SCALESvtk.SetName('scale')
for ii in range(len(NIN)):
	NINvtk.SetTuple1(ii,NIN[ii])
	SCALESvtk.SetTuple1(ii,gW['Scales'][0,0][0][ii])

# Build tree out of CP list of "is a child of"
#	remembering that Matlab indices are 1-based and numpy/VTK 0-based
print 'Building graph'
dg = vtk.vtkMutableDirectedGraph()
edge_id = vtk.vtkIdTypeArray()
edge_id.SetName('edge_ids')
for ii in range(CP.size):
	dg.AddVertex()
for ii in range(CP.size):
	if CP[ii] > 0:		# CP already zero-based
		dg.AddGraphEdge(CP[ii],ii)		# Method for use with wrappers -- AddEdge() in C++
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

