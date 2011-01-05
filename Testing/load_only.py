import scipy.io

data_file = '/Users/emonson/Data/Fodava/EMoGWDataSets/cell_array_of_strings.mat'
MatInput = scipy.io.loadmat(data_file, struct_as_record=True, chars_as_strings=True)

xx = MatInput['xx']
print xx
