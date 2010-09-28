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


# GLOB the list of Python.so files (treat them like plugins, too, for
# fixup_bundle purposes since they will not be pulled in automatically
# by dependency analysis)
#
file(GLOB VTK_Python_Libs "/Users/emonson/Programming/VTK_git/VTK/build/bin/*Python.so")
file(GLOB VTKVTG_Python_Libs "/Users/emonson/Programming/VTK_git/vtkVTG/build/bin/*Python.so")
file(GLOB PyQt_Python_Libs "/Library/Python/2.6/site-packages/PyQt4/Qt*.so")

# file(GLOB VTK_Libs "/Users/emonson/Programming/VTK_git/VTK/build/bin/*.dylib")
# file(GLOB VTKVTG_Libs "/Users/emonson/Programming/VTK_git/vtkVTG/build/bin/*.dylib")

set (PluginList)
foreach (pluginname ${packaged_plugin_names})
  list (APPEND PluginList "/Users/emonson/Programming/ParaView_git/ParaView/serial/bin/lib${pluginname}.dylib")
endforeach()
 
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
    set(path "@executable_path/../Frameworks")
  endif(item MATCHES "\\.dylib$")
 
  # ...default embed *.so in the Frameworks folder:
  #
  if(item MATCHES "\\.so$")
    set(path "@executable_path/../Frameworks")
  endif(item MATCHES "\\.so$")
 
  # ...default embed *Python.so in the Resources/ExtensionModules folder:
  #
  if(item MATCHES "Python\\.so$")
    set(path "@executable_path/../Resources/ExtensionModules")
  endif(item MATCHES "Python\\.so$")
 
  # ...embed VTK Python libs from ${VTK_Python_Libs} 
  #    in the Resources/ExtensionModules/vtk folder:
  #
  list(FIND VTK_Python_Libs ${item} libFound)
  if(libFound GREATER -1)
    set(path "@executable_path/../Resources/ExtensionModules/vtk")
  endif(libFound GREATER -1)

  # ...embed VTKVTG Python libs from ${VTKVTG_Python_Libs} 
  #    in the Resources/ExtensionModules/vtkvtg folder:
  #
  list(FIND VTKVTG_Python_Libs ${item} libFound)
  if(libFound GREATER -1)
    set(path "@executable_path/../Resources/ExtensionModules/vtkvtg")
  endif(libFound GREATER -1)

  # ...embed PyQt4 Python libs from ${PyQt_Python_Libs} 
  #    in the Resources/ExtensionModules/PyQt4 folder:
  #
  list(FIND PyQt_Python_Libs ${item} libFound)
  if(libFound GREATER -1)
    set(path "@executable_path/../Resources/ExtensionModules/PyQt4")
  endif(libFound GREATER -1)

  # ...embed libqsqlite.dylib in the Plugins/sqldrivers folder:
  #
  # if(item MATCHES "libqsqlite\\.dylib$")
  #   set(path "@executable_path/../Plugins/sqldrivers")
  # endif(item MATCHES "libqsqlite\\.dylib$")
 
  set(${default_embedded_path_var} "${path}" PARENT_SCOPE)
endfunction(gp_item_default_embedded_path_override)
 
 
# Copy the .app bundle from the build tree to the install tree.
# Set up the InstallPrefix ENV var and execute the shell script:
#
set(InstallPrefix "/Users/emonson/Programming/Python/VTK/MultiScaleSVD/build")
set(ENV{InstallPrefix} "${InstallPrefix}")
 
# If using a shell script to generate initial bundle 
# -- should put bundlebuilder command here
execute_process(WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                COMMAND python2.6 makeapplication.py build
                )

set(bundle "${InstallPrefix}/MS_SVD_vis.app")
 
if(NOT EXISTS "${bundle}")
  message(FATAL_ERROR "error: have to generate bundle with bundlebuilder first: ${bundle}")
endif()
 
 
# Fixup the .app bundle in the install tree:
#
include(BundleUtilities)

# list(APPEND PluginList "/usr/local/Trolltech/Qt-4.6.2/plugins/sqldrivers/libqsqlite.dylib")
 
# Additional libs may be found in:
#   (not sure if I can get rid of above GLOB for these directories, then...
#
set(libs_path "/Users/emonson/Programming/VTK_git/VTK/build/bin")
list(APPEND libs_path "/Users/emonson/Programming/VTK_git/vtkVTG/build/bin")
list(APPEND libs_path "/Library/Python/2.6/site-packages/PyQt4")

list(REMOVE_DUPLICATES libs_path)
 
# Fix it!
#
fixup_bundle(
  "${bundle}"
  "${VTK_Python_Libs};${VTKVTG_Python_Libs};${PyQt_Python_Libs};"
  "${libs_path}"
  )
