from data_source import DataSource
ds = DataSource('/Users/emonson/Data/Fodava/EMoGWDataSets/mnist1_1k_20100324.mat')
nodeID = 60
tableAll = ds.GetNodeAllScaleCoeffTable(nodeID)
tableOne = ds.GetNodeOneScaleCoeffTable(nodeID)
print tableAll
print tableOne

