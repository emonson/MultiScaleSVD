#!/usr/bin/python

"""
main.pyw

Example application using PyQt4 (with Qt Designer) and VTK

14 Oct 2008 -- E Monson
"""

import sys
from PyQt4 import QtGui
from multiscale_svd_gui import MultiScaleSVDViews

def dependencies_for_myprogram():
    from scipy.sparse.csgraph import _validation

if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    window = MultiScaleSVDViews()
    window.show()
    sys.exit(app.exec_())
