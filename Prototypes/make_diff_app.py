### makeapplication.py
from bundlebuilder import buildapp
import glob, os

package_root = '/Users/emonson/Programming/Python/VTK/QtTests'
qt_root = '/usr/local/Trolltech/Qt-4.7.0/lib'

lib_list = []
lib_list.append(os.path.join(qt_root,'QtCore.framework'))
lib_list.append(os.path.join(qt_root,'QtGui.framework'))

resource_list = []
# resource_list.append(os.path.join(package_root,'BlankImage.png'))

buildapp(
    name='DiffEmbedNav_0117.app',
    mainprogram='diff_w_slider.py', # your app's main()
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
