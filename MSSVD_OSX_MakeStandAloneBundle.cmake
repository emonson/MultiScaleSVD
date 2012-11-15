# Make sure this CMake has BundleUtilities.cmake:
#
if(NOT EXISTS "${CMAKE_ROOT}/Modules/BundleUtilities.cmake")
  message(FATAL_ERROR "error: BundleUtilities.cmake not found. Use CMake 2.6.4 or later.")
endif(NOT EXISTS "${CMAKE_ROOT}/Modules/BundleUtilities.cmake")
 
 
# Avoid following symlinks encountered during FILE GLOB_RECURSE calls:
#
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0009 NEW)
endif(COMMAND CMAKE_POLICY)
 
# Allow include to do cmake_policy push/pops:
#
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0011 NEW)
endif(COMMAND CMAKE_POLICY)

# set (PluginList)
# foreach (pluginname ${packaged_plugin_names})
#   list (APPEND PluginList "/Users/emonson/Programming/ParaView_git/ParaView/serial/bin/lib${pluginname}.dylib")
# endforeach()

# list(APPEND PluginList "/usr/local/Trolltech/Qt-4.7.4/plugins/sqldrivers/libqsqlite.dylib")
 
 
# gp_item_default_embedded_path_override item default_embedded_path_var
#
# Return the path that others should refer to the item by when the item
# is embedded inside a bundle.
#
# This is a project-specific override of BundleUtilities.cmake's
# gp_item_default_embedded_path
#
function(gp_item_default_embedded_path_override item default_embedded_path_var)
  # By default, embed items as set by gp_item_default_embedded_path:
  #
  set(path "${${default_embedded_path_var}}")
 
  # But for ParaView...
  #
  # ...embed *.dylib in the Frameworks folder:
  #
  if(item MATCHES "\\.dylib$")
    message(" * * dylib match " ${item} )
    set(path "@executable_path/../Frameworks")
  endif(item MATCHES "\\.dylib$")
 
  # ...default embed *.so in the lib-dynload folder:
  #
  if(item MATCHES "\\.so$")
    message(" o o so lib " ${item} )
    set(path "@executable_path/../Resources/lib/python2.7/lib-dynload")
  endif(item MATCHES "\\.so$")
 
  # ...default embed *PythonSIP.so in the Resources/lib/python2.7/lib-dynload/vtk folder:
  #
  if(item MATCHES "PythonSIP\\.so$")
    message(" s s sip lib " ${item} )
    set(path "@executable_path/../Resources/lib/python2.7/lib-dynload/vtk")
  endif(item MATCHES "PythonSIP\\.so$")
 
  # ...default embed QVTKPython.so in the Resources/lib/python2.7/lib-dynload/vtk folder:
  #
  if(item MATCHES "/QVTKPython\\.so$")
    message(" - - qvtk lib " ${item} )
    set(path "@executable_path/../Resources/lib/python2.7/lib-dynload/vtk")
  endif(item MATCHES "/QVTKPython\\.so$")
 
  # ...default embed QVTKPython.so in the Resources/lib/python2.7/lib-dynload/vtk folder:
  #
  if(item MATCHES "/Qt[A-Za-z]*\\.so$")
    message(" + + pyqt lib " ${item} )
    set(path "@executable_path/../Resources/lib/python2.7/lib-dynload/PyQt4")
  endif(item MATCHES "/Qt[A-Za-z]*\\.so$")
 
  # ...default embed QVTKPython.so in the Resources/lib/python2.7/lib-dynload/vtk folder:
  #
  if(item MATCHES "Python$")
    message("* * Python lib " ${item} )
    set(path "@executable_path/../Frameworks/Python.framework/Versions/2.7")
  endif(item MATCHES "Python$")
 
 
  set(${default_embedded_path_var} "${path}" PARENT_SCOPE)
endfunction(gp_item_default_embedded_path_override)
 
 
# Copy the .app bundle from the build tree to the install tree.
# Set up the InstallPrefix ENV var and execute the shell script:
#
set(InstallPrefix "/Users/emonson/Programming/Python/VTK/MultiScaleSVD/dist")
set(ENV{InstallPrefix} "${InstallPrefix}")
 
# If using a shell script to generate initial bundle 
# -- should put bundlebuilder command here
# execute_process(WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                 COMMAND python2.7 makeapplication.py build
#                 )

# Relying on python setup.py py2app already being run...

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
file(GLOB VTK_Python_Libs "${bundle}/Contents/Resources/lib/python2.7/lib-dynload/*.so")
file(GLOB VTK_Qt_Libs "${bundle}/Contents/Resources/lib/python2.7/lib-dynload/vtk/*.so")
file(GLOB PyQt_Python_Libs "${bundle}/Contents/Resources/lib/python2.7/lib-dynload/PyQt4/Qt*.so")
file(GLOB QtCore_Framework_Lib "${bundle}/Contents/Frameworks/QtCore.framework/Versions/Current/QtCore")
file(GLOB QtGui_Framework_Lib "${bundle}/Contents/Frameworks/QtGui.framework/Versions/Current/QtGui")
file(GLOB Python_Lib "${bundle}/Contents/Frameworks/Python.framework/Versions/2.7/Python")
file(GLOB Python_AppLib "${bundle}/Contents/Frameworks/Python.framework/Versions/2.7/Python.framework/Resources/Python.app/Contents/MacOS/Python")
file(GLOB Python_OrigLib "/usr/local/Cellar/python/2.7.3/Frameworks/Python.framework/Versions/2.7/Python")
message(  "${VTK_Python_Libs}\n\n
           ${VTK_Qt_Libs}\n\n
           ${PyQt_Python_Libs}\n\n
           ${QtCore_Framework_Lib}\n\n
           ${QtGui_Framework_Lib}\n\n
           ${Python_Lib}\n\n
           ${Python_AppLib}\n\n
           ${Python_OrigLib}\n\n"
)

# Additional libs may be found in:
set(libs_path "/Users/emonson/Programming/VTK_git/VTK/build/bin")
list(APPEND libs_path "/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
list(APPEND libs_path "/usr/local/lib/python2.7/site-packages/PyQt4")
# list(APPEND libs_path "/usr/local/opt/python/Frameworks/Python.framework/Versions/2.7")

list(REMOVE_DUPLICATES libs_path)

# Fix it!
#
fixup_bundle(
  "${bundle}"
  "${VTK_Python_Libs};${VTK_Qt_Libs};${PyQt_Python_Libs};${QtCore_Framework_Lib};${QtGui_Framework_Lib};${Python_Lib};${Python_AppLib};${Python_OrigLib}"
  "${libs_path}"
  )

# install_name_tool -change /usr/local/Cellar/python/2.7.3/Frameworks/Python.framework/Versions/2.7/Python @executable_path/../Frameworks/Python.framework/Versions/2.7/Python dist/main.app/Contents/Frameworks/Python.framework/Versions/2.7/Python.framework/Resources/Python.app/Contents/MacOS/Python