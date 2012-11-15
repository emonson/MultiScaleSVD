set(InstallPrefix "/Users/emonson/Programming/Python/VTK/MultiScaleSVD/dist")
set(ENV{InstallPrefix} "${InstallPrefix}")
 
set(bundle "${InstallPrefix}/main.app")
 
if(NOT EXISTS "${bundle}")
  message(FATAL_ERROR "error: have to generate bundle with bundlebuilder first: ${bundle}")
endif()
 
# Fixup the .app bundle in the install tree:
#
include(BundleUtilities)

# GLOB the list of Python.so files (treat them like plugins, too, for
# fixup_bundle purposes since they will not be pulled in automatically
# by dependency analysis)
#
file(GLOB VTK_Python_Libs "${bundle}/Contents/Resources/ExtensionModules/*Python.so")
file(GLOB VTKVTG_Python_Libs "${bundle}/Contents/Resources/ExtensionModules/vtkvtg/*Python.so")
file(GLOB PyQt_Python_Libs "${bundle}/Contents/Resources/ExtensionModules/PyQt/*Python.so")

# Additional libs may be found in:
set(libs_path "/Users/emonson/Programming/VTK_git/VTK/build/bin")
list(APPEND libs_path "/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
list(APPEND libs_path "/Library/Python/2.7/site-packages/PyQt4")

# list(REMOVE_DUPLICATES libs_path)
 
# Fix it!
#
fixup_bundle( "${bundle}" "" "" )
