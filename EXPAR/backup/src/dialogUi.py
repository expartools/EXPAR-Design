# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dialog2.ui'
#
# Created: Tue Mar 13 07:45:03 2012
#      by: PyQt4 UI code generator 4.8.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(966, 545)
        self.list = QtGui.QTreeWidget(Dialog)
        self.list.setGeometry(QtCore.QRect(0, 0, 956, 449))
        self.list.setMinimumSize(QtCore.QSize(497, 0))
        self.list.setObjectName(_fromUtf8("list"))
        self.pushButton = QtGui.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(0, 450, 956, 23))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton_2 = QtGui.QPushButton(Dialog)
        self.pushButton_2.setGeometry(QtCore.QRect(0, 480, 956, 23))
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.pushButton_3 = QtGui.QPushButton(Dialog)
        self.pushButton_3.setGeometry(QtCore.QRect(0, 510, 956, 23))
        self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(0, QtGui.QApplication.translate("Dialog", "ID", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(1, QtGui.QApplication.translate("Dialog", "Type", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(2, QtGui.QApplication.translate("Dialog", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(3, QtGui.QApplication.translate("Dialog", "Length", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(4, QtGui.QApplication.translate("Dialog", "Finger", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(5, QtGui.QApplication.translate("Dialog", "Exclude", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(6, QtGui.QApplication.translate("Dialog", "Include", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Dialog", "Submit", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_2.setText(QtGui.QApplication.translate("Dialog", "Reset", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_3.setText(QtGui.QApplication.translate("Dialog", "Save this result", None, QtGui.QApplication.UnicodeUTF8))

