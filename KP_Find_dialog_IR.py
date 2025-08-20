# -*- coding: utf-8 -*-

import os

from qgis.PyQt import uic
from qgis.PyQt import QtWidgets

# This loads your .ui file so that PyQt can populate your plugin with the elements from Qt Designer
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'KP_Find_dialog_IR.ui'))


class KpFindDialogInteractive(QtWidgets.QDialog, FORM_CLASS):
    KP_out = 0
    KP_DCC_out = 1
    DMM_KP_out = 2 
    KP_DMM_out = 3 #we will use these to compare later on in the program

    def __init__(self, parent=None):
        """Constructor."""
        super(KpFindDialogInteractive, self).__init__(parent)
        
        # Set up the user interface from Designer through FORM_CLASS.
        # After self.setupUi() you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.output_format_box.addItems([
            'KP 0.123', 
            'KP 0.123, DCC 12 m', 
            '42\xB0 42.7521\' N, 023\xB0 17.0568\' E (KP 0.123)', 
            'KP 0.123 (42\xB0 42.7521\' N, 023\xB0 17.0568\' E)'
        ]) #populate dropdown
        self.KP_prec.setValue(3) #set default values
        self.DCC_prec.setValue(0) #set default values
        self.Reverse_KP_box.setCheckState(False)
        self.Geodetic_box.setCheckState(2)
        self.Extendline_box.setCheckState(2)
        self.offset_m.setValue(0)
        