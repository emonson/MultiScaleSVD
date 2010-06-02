from data_source import DataSource
ds = DataSource('/Users/emonson/Data/Fodava/EMoGWDataSets/mnist1_5c_20100324.mat')
nodeID = 60
image_list = ds.GetDetailImages(nodeID)
multires_images = ds.GetMultiResolutionImages(nodeID)
print image_list
print image_list[0]
print multires_images

