# -*- coding: utf-8 -*-

import os

from qgis.PyQt import uic
from qgis.PyQt import QtWidgets

# This loads your .ui file so that PyQt can populate your plugin with the elements from Qt Designer
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'KP_Find_dialog_kp4p.ui'))


class KpFindDialogKp4p(QtWidgets.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(KpFindDialogKp4p, self).__init__(parent)
        
        # Set up the user interface from Designer through FORM_CLASS.
        # After self.setupUi() you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.KP_prec_kp4p.setValue(3) #set default values
        self.DCC_prec_kp4p.setValue(0) #set default values
        self.Reverse_KP_kp4p.setCheckState(False)
        self.Geodetic_kp4p.setCheckState(2)
        self.offset_m_kp4p.setValue(0)
        