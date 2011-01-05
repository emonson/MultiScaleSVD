### makeapplication.py
from bundlebuilder import buildapp
import glob, os

package_root = '/Users/emonson/Programming/Python/VTK/MultiScaleSVD'
qt_root = '/usr/local/Trolltech/Qt-4.7.0/lib'

lib_list = []
resource_list = []

buildapp(
    name='load_only.app', # what to build
    mainprogram='load_only.py', # your app's main()
    # argv_emulation=1, # drag&dropped filenames show up in sys.argv
    # iconfile='myapp.icns', # file containing your app's icons
    standalone=1, # make this app self contained.
    includeModules=['encodings.ascii','encodings.utf_8'], # list of additional Modules to force in
    includePackages=[], # list of additional Packages to force in
    resources=resource_list,
    libs=lib_list, # list of shared libs or Frameworks to include
)

### end of makeapplication.py
