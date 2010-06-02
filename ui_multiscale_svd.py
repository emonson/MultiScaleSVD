# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'multiscale_svd.ui'
#
# Created: Tue Jun  1 15:13:37 2010
#      by: PyQt4 UI code generator 4.7
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow, renWinList):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(900, 600)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.splitter_2 = QtGui.QSplitter(self.centralwidget)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.splitter_0 = QtGui.QSplitter(self.splitter_2)
        self.splitter_0.setOrientation(QtCore.Qt.Vertical)
        self.splitter_0.setObjectName("splitter_0")
        self.qvtkWidget_0 = QVTKRenderWindowInteractor(self.splitter_0, rw=renWinList[0])
        self.qvtkWidget_0.setObjectName("qvtkWidget_0")
        self.qvtkWidget_1 = QVTKRenderWindowInteractor(self.splitter_0, rw=renWinList[1])
        self.qvtkWidget_1.setObjectName("qvtkWidget_1")
        self.splitter_1 = QtGui.QSplitter(self.splitter_2)
        self.splitter_1.setOrientation(QtCore.Qt.Vertical)
        self.splitter_1.setObjectName("splitter_1")
        self.qvtkWidget_3 = QVTKRenderWindowInteractor(self.splitter_1, rw=renWinList[3])
        self.qvtkWidget_3.setObjectName("qvtkWidget_3")
        self.qvtkWidget_2 = QVTKRenderWindowInteractor(self.splitter_1, rw=renWinList[2])
        self.qvtkWidget_2.setObjectName("qvtkWidget_2")
        self.gridLayout.addWidget(self.splitter_2, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 900, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setText(QtGui.QApplication.translate("MainWindow", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setText(QtGui.QApplication.translate("MainWindow", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Q", None, QtGui.QApplication.UnicodeUTF8))

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
