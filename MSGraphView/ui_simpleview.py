# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'qsimpleview.ui'
#
# Created: Tue Jan 13 10:21:29 2009
#      by: PyQt4 UI code generator 4.4.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(790, 571)
        self.centralWidget = QtGui.QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")
        self.gridlayout = QtGui.QGridLayout(self.centralWidget)
        self.gridlayout.setObjectName("gridlayout")
        self.splitter_2 = QtGui.QSplitter(self.centralWidget)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.splitter = QtGui.QSplitter(self.splitter_2)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName("splitter")
        self.listWidget = QtGui.QListWidget(self.splitter)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.listWidget.sizePolicy().hasHeightForWidth())
        self.listWidget.setSizePolicy(sizePolicy)
        self.listWidget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.listWidget.setSelectionRectVisible(False)
        self.listWidget.setObjectName("listWidget")
        self.tabWidget = QtGui.QTabWidget(self.splitter)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(5)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setTabShape(QtGui.QTabWidget.Rounded)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtGui.QWidget()
        self.tab.setObjectName("tab")
        self.hSlider_basisIndex = QtGui.QSlider(self.tab)
        self.hSlider_basisIndex.setGeometry(QtCore.QRect(10, 160, 160, 22))
        self.hSlider_basisIndex.setMaximum(10)
        self.hSlider_basisIndex.setOrientation(QtCore.Qt.Horizontal)
        self.hSlider_basisIndex.setObjectName("hSlider_basisIndex")
        self.label_basisIndex = QtGui.QLabel(self.tab)
        self.label_basisIndex.setGeometry(QtCore.QRect(10, 140, 161, 17))
        self.label_basisIndex.setObjectName("label_basisIndex")
        self.spinBox_basisIndex = QtGui.QSpinBox(self.tab)
        self.spinBox_basisIndex.setGeometry(QtCore.QRect(180, 150, 63, 27))
        self.spinBox_basisIndex.setReadOnly(True)
        self.spinBox_basisIndex.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.spinBox_basisIndex.setObjectName("spinBox_basisIndex")
        self.label_level = QtGui.QLabel(self.tab)
        self.label_level.setGeometry(QtCore.QRect(10, 10, 171, 17))
        self.label_level.setObjectName("label_level")
        self.spinBox_level = QtGui.QSpinBox(self.tab)
        self.spinBox_level.setGeometry(QtCore.QRect(180, 20, 63, 27))
        self.spinBox_level.setReadOnly(True)
        self.spinBox_level.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.spinBox_level.setObjectName("spinBox_level")
        self.hSlider_level = QtGui.QSlider(self.tab)
        self.hSlider_level.setGeometry(QtCore.QRect(10, 30, 160, 22))
        self.hSlider_level.setOrientation(QtCore.Qt.Horizontal)
        self.hSlider_level.setObjectName("hSlider_level")
        self.label_cutoff = QtGui.QLabel(self.tab)
        self.label_cutoff.setGeometry(QtCore.QRect(10, 64, 221, 17))
        self.label_cutoff.setObjectName("label_cutoff")
        self.lineEdit_cutoff = QtGui.QLineEdit(self.tab)
        self.lineEdit_cutoff.setGeometry(QtCore.QRect(10, 84, 91, 22))
        self.lineEdit_cutoff.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_cutoff.setObjectName("lineEdit_cutoff")
        self.label_cutoff_minmax = QtGui.QLabel(self.tab)
        self.label_cutoff_minmax.setGeometry(QtCore.QRect(110, 90, 121, 17))
        font = QtGui.QFont()
        font.setPointSize(11)
        self.label_cutoff_minmax.setFont(font)
        self.label_cutoff_minmax.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.label_cutoff_minmax.setObjectName("label_cutoff_minmax")
        self.line = QtGui.QFrame(self.tab)
        self.line.setGeometry(QtCore.QRect(10, 110, 231, 16))
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName("line")
        self.label_basisCutoff_minmax = QtGui.QLabel(self.tab)
        self.label_basisCutoff_minmax.setGeometry(QtCore.QRect(110, 221, 121, 17))
        font = QtGui.QFont()
        font.setPointSize(11)
        self.label_basisCutoff_minmax.setFont(font)
        self.label_basisCutoff_minmax.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.label_basisCutoff_minmax.setObjectName("label_basisCutoff_minmax")
        self.label_basisCutoff = QtGui.QLabel(self.tab)
        self.label_basisCutoff.setGeometry(QtCore.QRect(10, 195, 221, 17))
        self.label_basisCutoff.setObjectName("label_basisCutoff")
        self.lineEdit_basisCutoff = QtGui.QLineEdit(self.tab)
        self.lineEdit_basisCutoff.setGeometry(QtCore.QRect(10, 215, 91, 22))
        self.lineEdit_basisCutoff.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_basisCutoff.setObjectName("lineEdit_basisCutoff")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
        self.vtkWidget = QVTKRenderWindowInteractor(self.splitter_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(10)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.vtkWidget.sizePolicy().hasHeightForWidth())
        self.vtkWidget.setSizePolicy(sizePolicy)
        self.vtkWidget.setObjectName("vtkWidget")
        self.gridlayout.addWidget(self.splitter_2, 0, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 790, 22))
        self.menuBar.setObjectName("menuBar")
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menuBar)
        self.actionNew = QtGui.QAction(MainWindow)
        self.actionNew.setObjectName("actionNew")
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.menuFile.addAction(self.actionNew)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuBar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.label_basisIndex.setText(QtGui.QApplication.translate("MainWindow", "Basis function index", None, QtGui.QApplication.UnicodeUTF8))
        self.label_level.setText(QtGui.QApplication.translate("MainWindow", "Multi-resolution level", None, QtGui.QApplication.UnicodeUTF8))
        self.label_cutoff.setText(QtGui.QApplication.translate("MainWindow", "Multi-resolution graph T cutoff", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_cutoff.setText(QtGui.QApplication.translate("MainWindow", "1.7e-5", None, QtGui.QApplication.UnicodeUTF8))
        self.label_cutoff_minmax.setText(QtGui.QApplication.translate("MainWindow", "(min, max)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_basisCutoff_minmax.setText(QtGui.QApplication.translate("MainWindow", "(min, max)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_basisCutoff.setText(QtGui.QApplication.translate("MainWindow", "Basis function value cutoff", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_basisCutoff.setText(QtGui.QApplication.translate("MainWindow", "1.0e-6", None, QtGui.QApplication.UnicodeUTF8))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), QtGui.QApplication.translate("MainWindow", "Basis", None, QtGui.QApplication.UnicodeUTF8))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), QtGui.QApplication.translate("MainWindow", "Display", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionNew.setText(QtGui.QApplication.translate("MainWindow", "New", None, QtGui.QApplication.UnicodeUTF8))
        self.actionNew.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+N", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setText(QtGui.QApplication.translate("MainWindow", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Q", None, QtGui.QApplication.UnicodeUTF8))

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor