# -*- coding: utf-8 -*-

import os

from qgis.PyQt import uic
from qgis.PyQt import QtWidgets

# This loads your .ui file so that PyQt can populate your plugin with the elements from Qt Designer
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'KP_Find_dialog_IR.ui'))


class KpFindDialogIR(QtWidgets.QDialog, FORM_CLASS):
    KP_out = 0
    KP_DCC_Out = 1
    DMS_out = 2 #we will use these to compare later on in the program

    def __init__(self, parent=None):
        """Constructor."""
        super(KpFindDialogIR, self).__init__(parent)
        
        # Set up the user interface from Designer through FORM_CLASS.
        # After self.setupUi() you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.output_format_box.addItems(['KP 0.123', 'KP 0.123, DCC 12 m', '05\xB0 13.2680\' S, 009\xB0 33.0401\' E (KP 0.123)']) #populate dropdown
        self.KP_prec.setValue(3) #set default values
        self.DCC_prec.setValue(0) #set default values
        self.Reverse_KP.setCheckState(False)
        self.offset_m.setValue(0)
        