### makeapplication.py
from bundlebuilder import buildapp
import glob, os

package_root = '/Users/emonson/Programming/Python/VTK/MultiScaleSVD'
qt_root = '/usr/local/Trolltech/Qt-4.7.0/lib'

lib_list = []
lib_list.append(os.path.join(qt_root,'QtCore.framework'))
lib_list.append(os.path.join(qt_root,'QtGui.framework'))

resource_list = []
resource_list.append(os.path.join(package_root,'BlankImage.png'))

# NOTE: Don't need to copy all of these libraries here, because it'll be taken
#       care of during dependency checks in MSSVD_OSX_MakeStandAloneBundle.cmake

# libpath = '/usr/local/lib/vtk-5.7/'
# for libfile in glob.glob( os.path.join(libpath, 'libvtk*.dylib') ):
#     liblist.append(libfile)
# 
# libpath = '/Users/emonson/Programming/VTK_git/vtkVTG/build/bin'
# for libfile in glob.glob( os.path.join(libpath, 'libvtk*.dylib') ):
#     lib_list.append(libfile)

# liblist.append('/Users/emonson/Programming/VTK_git/vtkVTG/build/bin/libvtkvtgCharts.dylib')
# liblist.append('/Users/emonson/Programming/VTK_git/vtkVTG/build/bin/libvtkvtgChartsPython.so')
# liblist.append('/Users/emonson/Programming/VTK_git/vtkVTG/build/bin/libvtkvtgChartsPythonD.dylib')

buildapp(
    name='MS_SVD_vis_color_0725.app', # what to build
    mainprogram='main.py', # your app's main()
    # argv_emulation=1, # drag&dropped filenames show up in sys.argv
    # iconfile='myapp.icns', # file containing your app's icons
    standalone=1, # make this app self contained.
    # encodings.utf_8 is necessary for loading cell arrays of strings
    includeModules=['encodings.ascii','encodings.utf_8','sip'], # list of additional Modules to force in
    includePackages=['vtk','vtkvtg','PyQt4'], # list of additional Packages to force in
    resources=resource_list,
    libs=lib_list, # list of shared libs or Frameworks to include
)

### end of makeapplication.py
