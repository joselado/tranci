# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'interface.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1187, 547)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_21 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_21.setObjectName("gridLayout_21")
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName("groupBox_3")
        self.gridLayout_20 = QtWidgets.QGridLayout(self.groupBox_3)
        self.gridLayout_20.setObjectName("gridLayout_20")
        self.logo = QtWidgets.QLabel(self.groupBox_3)
        self.logo.setMinimumSize(QtCore.QSize(200, 370))
        self.logo.setText("")
        self.logo.setObjectName("logo")
        self.gridLayout_20.addWidget(self.logo, 0, 0, 1, 1)
        self.gridLayout_21.addWidget(self.groupBox_3, 0, 0, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_17 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_17.setObjectName("gridLayout_17")
        self.tabWidget_3 = QtWidgets.QTabWidget(self.groupBox_2)
        self.tabWidget_3.setObjectName("tabWidget_3")
        self.tab_7 = QtWidgets.QWidget()
        self.tab_7.setObjectName("tab_7")
        self.gridLayout_16 = QtWidgets.QGridLayout(self.tab_7)
        self.gridLayout_16.setObjectName("gridLayout_16")
        self.initialize_one_shot = QtWidgets.QPushButton(self.tab_7)
        self.initialize_one_shot.setObjectName("initialize_one_shot")
        self.gridLayout_16.addWidget(self.initialize_one_shot, 0, 0, 1, 1)
        self.show_pdf = QtWidgets.QPushButton(self.tab_7)
        self.show_pdf.setObjectName("show_pdf")
        self.gridLayout_16.addWidget(self.show_pdf, 1, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_7, "")
        self.tab_8 = QtWidgets.QWidget()
        self.tab_8.setObjectName("tab_8")
        self.gridLayout_15 = QtWidgets.QGridLayout(self.tab_8)
        self.gridLayout_15.setObjectName("gridLayout_15")
        self.gridLayout_13 = QtWidgets.QGridLayout()
        self.gridLayout_13.setObjectName("gridLayout_13")
        self.label_24 = QtWidgets.QLabel(self.tab_8)
        self.label_24.setAlignment(QtCore.Qt.AlignCenter)
        self.label_24.setObjectName("label_24")
        self.gridLayout_13.addWidget(self.label_24, 0, 0, 1, 1)
        self.operator_name = QtWidgets.QComboBox(self.tab_8)
        self.operator_name.setObjectName("operator_name")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.operator_name.addItem("")
        self.gridLayout_13.addWidget(self.operator_name, 0, 1, 1, 1)
        self.gridLayout_15.addLayout(self.gridLayout_13, 1, 0, 1, 2)
        self.gridLayout_11 = QtWidgets.QGridLayout()
        self.gridLayout_11.setObjectName("gridLayout_11")
        self.plot_spectrum = QtWidgets.QPushButton(self.tab_8)
        self.plot_spectrum.setObjectName("plot_spectrum")
        self.gridLayout_11.addWidget(self.plot_spectrum, 1, 0, 1, 1)
        self.plot_excitations = QtWidgets.QPushButton(self.tab_8)
        self.plot_excitations.setObjectName("plot_excitations")
        self.gridLayout_11.addWidget(self.plot_excitations, 2, 0, 1, 1)
        self.plot_degeneracy = QtWidgets.QPushButton(self.tab_8)
        self.plot_degeneracy.setObjectName("plot_degeneracy")
        self.gridLayout_11.addWidget(self.plot_degeneracy, 3, 0, 1, 1)
        self.plot_operator = QtWidgets.QPushButton(self.tab_8)
        self.plot_operator.setObjectName("plot_operator")
        self.gridLayout_11.addWidget(self.plot_operator, 4, 0, 1, 1)
        self.label_20 = QtWidgets.QLabel(self.tab_8)
        self.label_20.setAlignment(QtCore.Qt.AlignBottom|QtCore.Qt.AlignHCenter)
        self.label_20.setObjectName("label_20")
        self.gridLayout_11.addWidget(self.label_20, 0, 0, 1, 1)
        self.gridLayout_15.addLayout(self.gridLayout_11, 2, 0, 1, 1)
        self.gridLayout_10 = QtWidgets.QGridLayout()
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.label_17 = QtWidgets.QLabel(self.tab_8)
        self.label_17.setObjectName("label_17")
        self.gridLayout_10.addWidget(self.label_17, 2, 0, 1, 1)
        self.initial_value = QtWidgets.QLineEdit(self.tab_8)
        self.initial_value.setObjectName("initial_value")
        self.gridLayout_10.addWidget(self.initial_value, 2, 1, 1, 1)
        self.final_value = QtWidgets.QLineEdit(self.tab_8)
        self.final_value.setObjectName("final_value")
        self.gridLayout_10.addWidget(self.final_value, 3, 1, 1, 1)
        self.steps = QtWidgets.QLineEdit(self.tab_8)
        self.steps.setObjectName("steps")
        self.gridLayout_10.addWidget(self.steps, 4, 1, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.tab_8)
        self.label_18.setObjectName("label_18")
        self.gridLayout_10.addWidget(self.label_18, 3, 0, 1, 1)
        self.label_19 = QtWidgets.QLabel(self.tab_8)
        self.label_19.setObjectName("label_19")
        self.gridLayout_10.addWidget(self.label_19, 4, 0, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.tab_8)
        self.label_21.setAlignment(QtCore.Qt.AlignBottom|QtCore.Qt.AlignHCenter)
        self.label_21.setObjectName("label_21")
        self.gridLayout_10.addWidget(self.label_21, 1, 0, 1, 2)
        self.gridLayout_15.addLayout(self.gridLayout_10, 2, 1, 1, 1)
        self.gridLayout_9 = QtWidgets.QGridLayout()
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.label_16 = QtWidgets.QLabel(self.tab_8)
        self.label_16.setAlignment(QtCore.Qt.AlignCenter)
        self.label_16.setObjectName("label_16")
        self.gridLayout_9.addWidget(self.label_16, 0, 0, 2, 1)
        self.sweep_variable = QtWidgets.QComboBox(self.tab_8)
        self.sweep_variable.setObjectName("sweep_variable")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.sweep_variable.addItem("")
        self.gridLayout_9.addWidget(self.sweep_variable, 0, 1, 2, 1)
        self.gridLayout_15.addLayout(self.gridLayout_9, 0, 0, 1, 2)
        self.tabWidget_3.addTab(self.tab_8, "")
        self.tab_6 = QtWidgets.QWidget()
        self.tab_6.setObjectName("tab_6")
        self.gridLayout_14 = QtWidgets.QGridLayout(self.tab_6)
        self.gridLayout_14.setObjectName("gridLayout_14")
        self.gridLayout_12 = QtWidgets.QGridLayout()
        self.gridLayout_12.setObjectName("gridLayout_12")
        self.lineEdit = QtWidgets.QLineEdit(self.tab_6)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout_12.addWidget(self.lineEdit, 0, 1, 1, 1)
        self.tol_ene = QtWidgets.QLineEdit(self.tab_6)
        self.tol_ene.setObjectName("tol_ene")
        self.gridLayout_12.addWidget(self.tol_ene, 1, 1, 1, 1)
        self.label_22 = QtWidgets.QLabel(self.tab_6)
        self.label_22.setObjectName("label_22")
        self.gridLayout_12.addWidget(self.label_22, 0, 0, 1, 1)
        self.label_23 = QtWidgets.QLabel(self.tab_6)
        self.label_23.setObjectName("label_23")
        self.gridLayout_12.addWidget(self.label_23, 1, 0, 1, 1)
        self.save_tranci = QtWidgets.QPushButton(self.tab_6)
        self.save_tranci.setObjectName("save_tranci")
        self.gridLayout_12.addWidget(self.save_tranci, 3, 0, 1, 2)
        self.label_26 = QtWidgets.QLabel(self.tab_6)
        self.label_26.setObjectName("label_26")
        self.gridLayout_12.addWidget(self.label_26, 2, 0, 1, 1)
        self.nwf_heff = QtWidgets.QLineEdit(self.tab_6)
        self.nwf_heff.setObjectName("nwf_heff")
        self.gridLayout_12.addWidget(self.nwf_heff, 2, 1, 1, 1)
        self.gridLayout_14.addLayout(self.gridLayout_12, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_6, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.gridLayout_22 = QtWidgets.QGridLayout(self.tab_5)
        self.gridLayout_22.setObjectName("gridLayout_22")
        self.gridLayout_19 = QtWidgets.QGridLayout()
        self.gridLayout_19.setObjectName("gridLayout_19")
        self.label_25 = QtWidgets.QLabel(self.tab_5)
        self.label_25.setAlignment(QtCore.Qt.AlignCenter)
        self.label_25.setObjectName("label_25")
        self.gridLayout_19.addWidget(self.label_25, 0, 0, 1, 3)
        self.zaxis_y = QtWidgets.QLineEdit(self.tab_5)
        self.zaxis_y.setObjectName("zaxis_y")
        self.gridLayout_19.addWidget(self.zaxis_y, 1, 1, 1, 1)
        self.zaxis_z = QtWidgets.QLineEdit(self.tab_5)
        self.zaxis_z.setObjectName("zaxis_z")
        self.gridLayout_19.addWidget(self.zaxis_z, 1, 2, 1, 1)
        self.zaxis_x = QtWidgets.QLineEdit(self.tab_5)
        self.zaxis_x.setObjectName("zaxis_x")
        self.gridLayout_19.addWidget(self.zaxis_x, 1, 0, 1, 1)
        self.gridLayout_22.addLayout(self.gridLayout_19, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_5, "")
        self.gridLayout_17.addWidget(self.tabWidget_3, 0, 0, 1, 1)
        self.gridLayout_21.addWidget(self.groupBox_2, 0, 1, 1, 1)
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_18 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_18.setObjectName("gridLayout_18")
        self.tabWidget = QtWidgets.QTabWidget(self.groupBox)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label = QtWidgets.QLabel(self.tab)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 1, 0, 1, 1)
        self.U = QtWidgets.QLineEdit(self.tab)
        self.U.setObjectName("U")
        self.gridLayout_2.addWidget(self.U, 1, 1, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.tab)
        self.label_15.setObjectName("label_15")
        self.gridLayout_2.addWidget(self.label_15, 0, 0, 1, 1)
        self.n = QtWidgets.QLineEdit(self.tab)
        self.n.setObjectName("n")
        self.gridLayout_2.addWidget(self.n, 0, 1, 1, 1)
        self.gridLayout_8.addLayout(self.gridLayout_2, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.tab_2)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_2 = QtWidgets.QLabel(self.tab_2)
        self.label_2.setObjectName("label_2")
        self.gridLayout_3.addWidget(self.label_2, 0, 0, 1, 1)
        self.soc = QtWidgets.QLineEdit(self.tab_2)
        self.soc.setObjectName("soc")
        self.gridLayout_3.addWidget(self.soc, 0, 1, 1, 1)
        self.gridLayout_7.addLayout(self.gridLayout_3, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.tab_3)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.gridLayout_5 = QtWidgets.QGridLayout()
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label_3 = QtWidgets.QLabel(self.tab_3)
        self.label_3.setObjectName("label_3")
        self.gridLayout_5.addWidget(self.label_3, 0, 0, 1, 1)
        self.D = QtWidgets.QLineEdit(self.tab_3)
        self.D.setObjectName("D")
        self.gridLayout_5.addWidget(self.D, 0, 1, 1, 1)
        self.E = QtWidgets.QLineEdit(self.tab_3)
        self.E.setObjectName("E")
        self.gridLayout_5.addWidget(self.E, 1, 1, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.tab_3)
        self.label_4.setObjectName("label_4")
        self.gridLayout_5.addWidget(self.label_4, 1, 0, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.tab_3)
        self.label_5.setObjectName("label_5")
        self.gridLayout_5.addWidget(self.label_5, 2, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.tab_3)
        self.label_7.setObjectName("label_7")
        self.gridLayout_5.addWidget(self.label_7, 4, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.tab_3)
        self.label_6.setObjectName("label_6")
        self.gridLayout_5.addWidget(self.label_6, 3, 0, 1, 1)
        self.trigonal = QtWidgets.QLineEdit(self.tab_3)
        self.trigonal.setObjectName("trigonal")
        self.gridLayout_5.addWidget(self.trigonal, 3, 1, 1, 1)
        self.O = QtWidgets.QLineEdit(self.tab_3)
        self.O.setObjectName("O")
        self.gridLayout_5.addWidget(self.O, 2, 1, 1, 1)
        self.z4 = QtWidgets.QLineEdit(self.tab_3)
        self.z4.setObjectName("z4")
        self.gridLayout_5.addWidget(self.z4, 4, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.tab_3)
        self.label_11.setObjectName("label_11")
        self.gridLayout_5.addWidget(self.label_11, 5, 0, 1, 1)
        self.x2y2 = QtWidgets.QLineEdit(self.tab_3)
        self.x2y2.setObjectName("x2y2")
        self.gridLayout_5.addWidget(self.x2y2, 5, 1, 1, 1)
        self.gridLayout_6.addLayout(self.gridLayout_5, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_3, "")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.gridLayout = QtWidgets.QGridLayout(self.tab_4)
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_4 = QtWidgets.QGridLayout()
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_12 = QtWidgets.QLabel(self.tab_4)
        self.label_12.setObjectName("label_12")
        self.gridLayout_4.addWidget(self.label_12, 0, 2, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.tab_4)
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 0, 0, 1, 1)
        self.B = QtWidgets.QLineEdit(self.tab_4)
        self.B.setObjectName("B")
        self.gridLayout_4.addWidget(self.B, 0, 1, 1, 1)
        self.theta_b = QtWidgets.QLineEdit(self.tab_4)
        self.theta_b.setObjectName("theta_b")
        self.gridLayout_4.addWidget(self.theta_b, 1, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.tab_4)
        self.label_10.setObjectName("label_10")
        self.gridLayout_4.addWidget(self.label_10, 2, 0, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.tab_4)
        self.label_9.setObjectName("label_9")
        self.gridLayout_4.addWidget(self.label_9, 1, 0, 1, 1)
        self.phi_b = QtWidgets.QLineEdit(self.tab_4)
        self.phi_b.setObjectName("phi_b")
        self.gridLayout_4.addWidget(self.phi_b, 2, 1, 1, 1)
        self.j = QtWidgets.QLineEdit(self.tab_4)
        self.j.setObjectName("j")
        self.gridLayout_4.addWidget(self.j, 0, 3, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.tab_4)
        self.label_13.setObjectName("label_13")
        self.gridLayout_4.addWidget(self.label_13, 1, 2, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.tab_4)
        self.label_14.setObjectName("label_14")
        self.gridLayout_4.addWidget(self.label_14, 2, 2, 1, 1)
        self.theta_j = QtWidgets.QLineEdit(self.tab_4)
        self.theta_j.setObjectName("theta_j")
        self.gridLayout_4.addWidget(self.theta_j, 1, 3, 1, 1)
        self.phi_j = QtWidgets.QLineEdit(self.tab_4)
        self.phi_j.setObjectName("phi_j")
        self.gridLayout_4.addWidget(self.phi_j, 2, 3, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_4, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_4, "")
        self.gridLayout_18.addWidget(self.tabWidget, 0, 0, 1, 1)
        self.gridLayout_21.addWidget(self.groupBox, 0, 2, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1187, 20))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget_3.setCurrentIndex(0)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Tranci"))
        self.groupBox_3.setTitle(_translate("MainWindow", "TranCI"))
        self.groupBox_2.setTitle(_translate("MainWindow", "Calculation"))
        self.initialize_one_shot.setToolTip(_translate("MainWindow", "Prform a single shot calculation"))
        self.initialize_one_shot.setText(_translate("MainWindow", "Initialize and run"))
        self.show_pdf.setToolTip(_translate("MainWindow", "Show a pdf with the summary of the results"))
        self.show_pdf.setText(_translate("MainWindow", "Show pdf"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_7), _translate("MainWindow", "Single shot"))
        self.label_24.setText(_translate("MainWindow", "Operator"))
        self.operator_name.setToolTip(_translate("MainWindow", "Operator used to compute expectation values"))
        self.operator_name.setItemText(0, _translate("MainWindow", "Sx"))
        self.operator_name.setItemText(1, _translate("MainWindow", "Sy"))
        self.operator_name.setItemText(2, _translate("MainWindow", "Sz"))
        self.operator_name.setItemText(3, _translate("MainWindow", "Lx"))
        self.operator_name.setItemText(4, _translate("MainWindow", "Ly"))
        self.operator_name.setItemText(5, _translate("MainWindow", "Lz"))
        self.operator_name.setItemText(6, _translate("MainWindow", "Jx"))
        self.operator_name.setItemText(7, _translate("MainWindow", "Jy"))
        self.operator_name.setItemText(8, _translate("MainWindow", "Jz"))
        self.operator_name.setItemText(9, _translate("MainWindow", "LS"))
        self.operator_name.setItemText(10, _translate("MainWindow", "L2"))
        self.operator_name.setItemText(11, _translate("MainWindow", "S2"))
        self.operator_name.setItemText(12, _translate("MainWindow", "J2"))
        self.operator_name.setItemText(13, _translate("MainWindow", "x2"))
        self.operator_name.setItemText(14, _translate("MainWindow", "y2"))
        self.operator_name.setItemText(15, _translate("MainWindow", "z2"))
        self.operator_name.setItemText(16, _translate("MainWindow", "S_(111)"))
        self.operator_name.setItemText(17, _translate("MainWindow", "L_(111)"))
        self.plot_spectrum.setToolTip(_translate("MainWindow", "Plot the spectrum of the many body system"))
        self.plot_spectrum.setText(_translate("MainWindow", "Plot spectrum"))
        self.plot_excitations.setToolTip(_translate("MainWindow", "Plot the excitations with respect to the ground state"))
        self.plot_excitations.setText(_translate("MainWindow", "Plot excitations"))
        self.plot_degeneracy.setToolTip(_translate("MainWindow", "Plot the degeneracy of the ground state"))
        self.plot_degeneracy.setText(_translate("MainWindow", "Plot degeneracy"))
        self.plot_operator.setToolTip(_translate("MainWindow", "Plot the expectation value of a certain operator on the ground state manifold"))
        self.plot_operator.setText(_translate("MainWindow", "Plot operator"))
        self.label_20.setText(_translate("MainWindow", "Tasks"))
        self.label_17.setText(_translate("MainWindow", "Initial"))
        self.initial_value.setToolTip(_translate("MainWindow", "Initial value of the sweep"))
        self.initial_value.setText(_translate("MainWindow", "0.0"))
        self.final_value.setToolTip(_translate("MainWindow", "Final value of the sweep"))
        self.final_value.setText(_translate("MainWindow", "1.0"))
        self.steps.setToolTip(_translate("MainWindow", "Number of steps in the sweep"))
        self.steps.setText(_translate("MainWindow", "40"))
        self.label_18.setText(_translate("MainWindow", "Final"))
        self.label_19.setText(_translate("MainWindow", "Steps"))
        self.label_21.setText(_translate("MainWindow", "Values to sweep"))
        self.label_16.setText(_translate("MainWindow", "Parameter to sweep"))
        self.sweep_variable.setItemText(0, _translate("MainWindow", "soc"))
        self.sweep_variable.setItemText(1, _translate("MainWindow", "z^2"))
        self.sweep_variable.setItemText(2, _translate("MainWindow", "z^4"))
        self.sweep_variable.setItemText(3, _translate("MainWindow", "x^2-y^2"))
        self.sweep_variable.setItemText(4, _translate("MainWindow", "x^4+y^4+z^4"))
        self.sweep_variable.setItemText(5, _translate("MainWindow", "(x+y+z)^2"))
        self.sweep_variable.setItemText(6, _translate("MainWindow", "B"))
        self.sweep_variable.setItemText(7, _translate("MainWindow", "theta_B"))
        self.sweep_variable.setItemText(8, _translate("MainWindow", "phi_B"))
        self.sweep_variable.setItemText(9, _translate("MainWindow", "J"))
        self.sweep_variable.setItemText(10, _translate("MainWindow", "theta_J"))
        self.sweep_variable.setItemText(11, _translate("MainWindow", "phi_J"))
        self.sweep_variable.setItemText(12, _translate("MainWindow", "U"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_8), _translate("MainWindow", "Sweep"))
        self.lineEdit.setText(_translate("MainWindow", "40"))
        self.tol_ene.setText(_translate("MainWindow", "0.001"))
        self.label_22.setText(_translate("MainWindow", "Maximum number of eigenfunctions"))
        self.label_23.setText(_translate("MainWindow", "Energy tolerancy [eV]"))
        self.save_tranci.setToolTip(_translate("MainWindow", "Save data in the current folder"))
        self.save_tranci.setText(_translate("MainWindow", "Save data"))
        self.label_26.setText(_translate("MainWindow", "# of WF for effective Hamiltonian"))
        self.nwf_heff.setToolTip(_translate("MainWindow", "# of low energy wavefunctions in the effective Hamiltonian. This would be the dimensional of the local effective model (for example, 2S+1)"))
        self.nwf_heff.setText(_translate("MainWindow", "9"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_6), _translate("MainWindow", "Control"))
        self.label_25.setToolTip(_translate("MainWindow", "axis for the wavefunctions, controls what is the direction with respect to the initial axis to represent the harmonic and spin-part of the wavefunction"))
        self.label_25.setText(_translate("MainWindow", "z axis for the wavefunctions"))
        self.zaxis_y.setText(_translate("MainWindow", "0"))
        self.zaxis_z.setText(_translate("MainWindow", "1"))
        self.zaxis_x.setText(_translate("MainWindow", "0"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_5), _translate("MainWindow", "Advanced"))
        self.groupBox.setTitle(_translate("MainWindow", "Parameters"))
        self.label.setText(_translate("MainWindow", "U [eV]"))
        self.U.setToolTip(_translate("MainWindow", "Electron-electron interaction"))
        self.U.setText(_translate("MainWindow", "4.0"))
        self.label_15.setText(_translate("MainWindow", "# of electrons"))
        self.n.setToolTip(_translate("MainWindow", "Number of electrons in the d shell"))
        self.n.setText(_translate("MainWindow", "2"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "General"))
        self.label_2.setText(_translate("MainWindow", "SOC [eV]"))
        self.soc.setToolTip(_translate("MainWindow", "Value of the spin-orbit coupling"))
        self.soc.setText(_translate("MainWindow", "0.05"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Spin-orbit"))
        self.label_3.setText(_translate("MainWindow", "Uniaxial z^2 [eV]"))
        self.D.setToolTip(_translate("MainWindow", "Uniaxial crystal field"))
        self.D.setText(_translate("MainWindow", "0.0"))
        self.E.setText(_translate("MainWindow", "0.0"))
        self.label_4.setText(_translate("MainWindow", "Shear x^2-y^2 [eV]"))
        self.label_5.setText(_translate("MainWindow", "Octahedral x^4+y^4+z^4"))
        self.label_7.setText(_translate("MainWindow", "Uniaxial\' z^4 [eV]"))
        self.label_6.setText(_translate("MainWindow", "Trigonal (x+y+z)^2"))
        self.trigonal.setToolTip(_translate("MainWindow", "Trigonal crystal field"))
        self.trigonal.setText(_translate("MainWindow", "0.05"))
        self.O.setToolTip(_translate("MainWindow", "Octahedral crystal field"))
        self.O.setText(_translate("MainWindow", "0.3"))
        self.z4.setToolTip(_translate("MainWindow", "Uniaxial crystal field"))
        self.z4.setText(_translate("MainWindow", "0.0"))
        self.label_11.setText(_translate("MainWindow", "Shear\' x^2y^2"))
        self.x2y2.setText(_translate("MainWindow", "0.0"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Crystal field"))
        self.label_12.setText(_translate("MainWindow", "J [eV] "))
        self.label_8.setText(_translate("MainWindow", "B [eV]"))
        self.B.setToolTip(_translate("MainWindow", "Absolute value of the magnetic field"))
        self.B.setText(_translate("MainWindow", "0.0"))
        self.theta_b.setText(_translate("MainWindow", "0.0"))
        self.label_10.setText(_translate("MainWindow", "Phi_B [pi]"))
        self.label_9.setText(_translate("MainWindow", "Theta_B [pi]"))
        self.phi_b.setText(_translate("MainWindow", "0.0"))
        self.j.setToolTip(_translate("MainWindow", "Absolute value of an exchange field that couples only to the spin sector"))
        self.j.setText(_translate("MainWindow", "0.0"))
        self.label_13.setText(_translate("MainWindow", "Theta_J [pi]"))
        self.label_14.setText(_translate("MainWindow", "Phi_J [pi]"))
        self.theta_j.setText(_translate("MainWindow", "0.0"))
        self.phi_j.setText(_translate("MainWindow", "0.0"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("MainWindow", "Magnetic field"))

