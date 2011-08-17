#!/usr/bin/env python

"""
main.pyw

Example application using PyQt4 (with Qt Designer) and VTK

14 Oct 2008 -- E Monson
"""

import sys
from PyQt4.QtGui import QApplication
from simpleview import SimpleView

if __name__ == "__main__":

    app = QApplication(sys.argv)
    window = SimpleView()
    window.show()
    sys.exit(app.exec_())
