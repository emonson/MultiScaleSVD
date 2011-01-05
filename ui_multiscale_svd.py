# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'multiscale_svd.ui'
#
# Created: Thu Dec 30 12:16:37 2010
#      by: PyQt4 UI code generator 4.8.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow, renWinList):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1100, 500)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.splitter_4 = QtGui.QSplitter(self.centralwidget)
        self.splitter_4.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_4.setObjectName(_fromUtf8("splitter_4"))
        self.splitter_2 = QtGui.QSplitter(self.splitter_4)
        self.splitter_2.setOrientation(QtCore.Qt.Vertical)
        self.splitter_2.setObjectName(_fromUtf8("splitter_2"))
        self.qvtkWidget_5 = QVTKRenderWindowInteractor(self.splitter_2, rw=renWinList[5])
        self.qvtkWidget_5.setObjectName(_fromUtf8("qvtkWidget_5"))
        self.qvtkWidget_4 = QVTKRenderWindowInteractor(self.splitter_2, rw=renWinList[4])
        self.qvtkWidget_4.setObjectName(_fromUtf8("qvtkWidget_4"))
        self.splitter_3 = QtGui.QSplitter(self.splitter_4)
        self.splitter_3.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_3.setObjectName(_fromUtf8("splitter_3"))
        self.splitter_1 = QtGui.QSplitter(self.splitter_3)
        self.splitter_1.setOrientation(QtCore.Qt.Vertical)
        self.splitter_1.setObjectName(_fromUtf8("splitter_1"))
        self.qvtkWidget_0 = QVTKRenderWindowInteractor(self.splitter_1, rw=renWinList[0])
        self.qvtkWidget_0.setObjectName(_fromUtf8("qvtkWidget_0"))
        self.qvtkWidget_1 = QVTKRenderWindowInteractor(self.splitter_1, rw=renWinList[1])
        self.qvtkWidget_1.setObjectName(_fromUtf8("qvtkWidget_1"))
        self.splitter_0 = QtGui.QSplitter(self.splitter_3)
        self.splitter_0.setOrientation(QtCore.Qt.Vertical)
        self.splitter_0.setObjectName(_fromUtf8("splitter_0"))
        self.qvtkWidget_3 = QVTKRenderWindowInteractor(self.splitter_0, rw=renWinList[3])
        self.qvtkWidget_3.setObjectName(_fromUtf8("qvtkWidget_3"))
        self.qvtkWidget_2 = QVTKRenderWindowInteractor(self.splitter_0, rw=renWinList[2])
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.qvtkWidget_2.sizePolicy().hasHeightForWidth())
        self.qvtkWidget_2.setSizePolicy(sizePolicy)
        self.qvtkWidget_2.setObjectName(_fromUtf8("qvtkWidget_2"))
        self.gridLayout.addWidget(self.splitter_4, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1100, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuView = QtGui.QMenu(self.menubar)
        self.menuView.setObjectName(_fromUtf8("menuView"))
        self.menuPC_Scales = QtGui.QMenu(self.menuView)
        self.menuPC_Scales.setObjectName(_fromUtf8("menuPC_Scales"))
        self.menuPlot_Colors = QtGui.QMenu(self.menuView)
        self.menuPlot_Colors.setObjectName(_fromUtf8("menuPlot_Colors"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionWavelet = QtGui.QAction(MainWindow)
        self.actionWavelet.setCheckable(True)
        self.actionWavelet.setObjectName(_fromUtf8("actionWavelet"))
        self.actionScaling = QtGui.QAction(MainWindow)
        self.actionScaling.setCheckable(True)
        self.actionScaling.setObjectName(_fromUtf8("actionScaling"))
        self.actionPC_All_Scales = QtGui.QAction(MainWindow)
        self.actionPC_All_Scales.setCheckable(True)
        self.actionPC_All_Scales.setObjectName(_fromUtf8("actionPC_All_Scales"))
        self.actionPC_Current_Scale = QtGui.QAction(MainWindow)
        self.actionPC_Current_Scale.setCheckable(True)
        self.actionPC_Current_Scale.setObjectName(_fromUtf8("actionPC_Current_Scale"))
        self.actionPC_Coarsest_to_Current = QtGui.QAction(MainWindow)
        self.actionPC_Coarsest_to_Current.setCheckable(True)
        self.actionPC_Coarsest_to_Current.setObjectName(_fromUtf8("actionPC_Coarsest_to_Current"))
        self.actionPC_Current_to_Finest = QtGui.QAction(MainWindow)
        self.actionPC_Current_to_Finest.setCheckable(True)
        self.actionPC_Current_to_Finest.setObjectName(_fromUtf8("actionPC_Current_to_Finest"))
        self.actionColorNone = QtGui.QAction(MainWindow)
        self.actionColorNone.setCheckable(True)
        self.actionColorNone.setObjectName(_fromUtf8("actionColorNone"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionExit)
        self.menuPC_Scales.addAction(self.actionPC_All_Scales)
        self.menuPC_Scales.addAction(self.actionPC_Coarsest_to_Current)
        self.menuPC_Scales.addAction(self.actionPC_Current_Scale)
        self.menuPC_Scales.addAction(self.actionPC_Current_to_Finest)
        self.menuPlot_Colors.addAction(self.actionColorNone)
        self.menuView.addAction(self.actionWavelet)
        self.menuView.addAction(self.actionScaling)
        self.menuView.addSeparator()
        self.menuView.addAction(self.menuPC_Scales.menuAction())
        self.menuView.addSeparator()
        self.menuView.addAction(self.menuPlot_Colors.menuAction())
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuView.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Multi-scale SVD", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.menuView.setTitle(QtGui.QApplication.translate("MainWindow", "View", None, QtGui.QApplication.UnicodeUTF8))
        self.menuPC_Scales.setTitle(QtGui.QApplication.translate("MainWindow", "PC Scales", None, QtGui.QApplication.UnicodeUTF8))
        self.menuPlot_Colors.setTitle(QtGui.QApplication.translate("MainWindow", "Plot Color by Array", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setText(QtGui.QApplication.translate("MainWindow", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setText(QtGui.QApplication.translate("MainWindow", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Q", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWavelet.setText(QtGui.QApplication.translate("MainWindow", "Wavelet Bases", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWavelet.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Alt+W", None, QtGui.QApplication.UnicodeUTF8))
        self.actionScaling.setText(QtGui.QApplication.translate("MainWindow", "Scaling Function Bases", None, QtGui.QApplication.UnicodeUTF8))
        self.actionScaling.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Alt+S", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_All_Scales.setText(QtGui.QApplication.translate("MainWindow", "PC All Scales", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_All_Scales.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Alt+1", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_Current_Scale.setText(QtGui.QApplication.translate("MainWindow", "PC Current Scale", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_Current_Scale.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Alt+3", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_Coarsest_to_Current.setText(QtGui.QApplication.translate("MainWindow", "PC Coarsest to Current", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_Coarsest_to_Current.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Alt+2", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_Current_to_Finest.setText(QtGui.QApplication.translate("MainWindow", "PC Current to Finest", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPC_Current_to_Finest.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Alt+4", None, QtGui.QApplication.UnicodeUTF8))
        self.actionColorNone.setText(QtGui.QApplication.translate("MainWindow", "None", None, QtGui.QApplication.UnicodeUTF8))

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
