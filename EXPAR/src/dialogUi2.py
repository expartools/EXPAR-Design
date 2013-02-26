# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dialog2.ui'
#
# Created: Tue Mar 20 11:38:35 2012
#      by: PyQt4 UI code generator 4.8.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Dialog2(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(966, 518)
        self.list = QtGui.QTreeWidget(Dialog)
        self.list.setGeometry(QtCore.QRect(0, 0, 956, 449))
        self.list.setMinimumSize(QtCore.QSize(497, 0))
        self.list.setObjectName(_fromUtf8("list"))
        self.pushButton = QtGui.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(0, 470, 956, 23))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(0, QtGui.QApplication.translate("Dialog", "Finger_ID", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(1, QtGui.QApplication.translate("Dialog", "Type", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(2, QtGui.QApplication.translate("Dialog", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(3, QtGui.QApplication.translate("Dialog", "Tigger_generated", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(4, QtGui.QApplication.translate("Dialog", "Trigger", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(5, QtGui.QApplication.translate("Dialog", "Template", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(6, QtGui.QApplication.translate("Dialog", "Bayes_class", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(7, QtGui.QApplication.translate("Dialog", "Pwm_class", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(8, QtGui.QApplication.translate("Dialog", "P90_score", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(9, QtGui.QApplication.translate("Dialog", "Diff_score", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(10, QtGui.QApplication.translate("Dialog", "Tri-temp Tm", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(11, QtGui.QApplication.translate("Dialog", "Temp-temp Tm", None, QtGui.QApplication.UnicodeUTF8))
        self.list.headerItem().setText(12, QtGui.QApplication.translate("Dialog", "Bonds", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Dialog", "Save as", None, QtGui.QApplication.UnicodeUTF8))

