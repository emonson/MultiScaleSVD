#!/usr/bin/env python

import os

command = "pyuic4 -o temp_simpleview.py qsimpleview.ui"
print 'Generating ui_simpleview.py from qsimpleview.ui'
os.system(command)

# Need to do some substitution in generated file to make the UI work in python
# when the QVTKWidget was created to be used in C++

f = open('temp_simpleview.py', 'r')
s = f.read()
f.close()

s = s.replace('QVTKWidget','QVTKRenderWindowInteractor')
s = s.replace('from QVTK', 'from vtk.qt4.QVTK')

w = open('ui_simpleview.py', 'w')
w.write(s)
w.close()

os.remove('temp_simpleview.py')
