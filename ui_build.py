#!/usr/bin/env python

import os, re

command = "pyuic4 -o temp_multiscale_svd.py multiscale_svd.ui"
print 'Generating ui_multiscale_svd.py from multiscale_svd.ui'
os.system(command)

# Need to do some substitution in generated file to make the UI work in python
# when the QVTKWidget was created to be used in C++

f = open('temp_multiscale_svd.py', 'r')
s = f.read()
f.close()

s = s.replace('QVTKWidget','QVTKRenderWindowInteractor')
s = s.replace('from QVTK', 'from vtk.qt4.QVTK')

# For vtkRenderView subclasses need to pass the render window to constructor on setup
s = re.sub(r'def setupUi\((.+?)\):', 'def setupUi(\g<1>, renWinList):', s)
s = re.sub(r'qvtkWidget_(\d) = QVTKRenderWindowInteractor\((.+?)\)', 'qvtkWidget_\g<1> = QVTKRenderWindowInteractor(\g<2>, rw=renWinList[\g<1>])', s)

w = open('ui_multiscale_svd.py', 'w')
w.write(s)
w.close()

os.remove('temp_multiscale_svd.py')
