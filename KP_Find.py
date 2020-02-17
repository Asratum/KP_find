# -*- coding: utf-8 -*-
import site
import os
site.addsitedir(os.path.abspath(os.path.dirname(__file__) + '/ext-libs'))



from qgis.PyQt.QtCore import Qt, QSettings, QTranslator, QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker

from PyQt5.QtWidgets import QAction, QMessageBox #give messages

from qgis.core import (QgsCoordinateReferenceSystem,
                       QgsProject,
                       QgsMapLayerProxyModel,
                       QgsVectorLayer,
                       Qgis)

from .kpTool import KPTool


# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .KP_Find_dialog_IR import KpFindDialogIR
from .KP_Find_dialog_kp4p import KpFindDialogKp4p
import os.path


class KpFind:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        self.canvas = iface.mapCanvas()

        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'KpFind_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        
        lineLayer= None
        self.actions = []
        self.menu = self.tr(u'&KP Find')
        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_startIR = None
        self.first_startKp4p = None


    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('KpFind', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):


        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)
        return action
        
    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""
    
        icon_path_KP_Int = ':/plugins/KP_Find/KP_Interactive.png'
        self.add_action(
            icon_path_KP_Int,
            text=self.tr(u'Interactive KP Find'),
            callback=self.Interactive_run,
            parent=self.iface.mainWindow())
        icon_path_KP_4p = ':/plugins/KP_Find/KP_4p.png'    
        self.add_action(
            icon_path_KP_4p,
            text=self.tr(u'KP and DOL for point layer'),
            callback=self.Kp4Points,
            parent=self.iface.mainWindow())
        
        # will be set False in run()
        self.first_startIR = True
        self.first_startKp4p = True
        self.mapTool = KPTool(self.iface)

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&KP Find'),
                action)
            self.iface.removePluginMenu(
                self.tr(u'&KP and DCC for point layer'),
                action)
            self.iface.removeToolBarIcon(action)

    def Interactive_run(self):
        """Run method that performs all the real work"""
        # Create the dialog with elements (after translation) and keep reference
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        if self.first_startIR == True:
            self.first_startIR = False
            self.dlgIR = KpFindDialogIR()
            self.dlgIR.inputLineLayer.setFilters(QgsMapLayerProxyModel.LineLayer) #allows only line layers in the box
            self.dlgIR.setWindowFlags(Qt.WindowStaysOnTopHint) #dialog window should always stay on top


        # show the dialog
        
        self.dlgIR.show()
        # Run the dialog event loop
        result = self.dlgIR.exec_()

        # See if OK was pressed
        if result:
            
            lineLayer=self.dlgIR.inputLineLayer.currentLayer() #gets the layer from the input QT box
            
            self.proj_crs = self.canvas.mapSettings().destinationCrs()
            check = lineLayer.crs() == self.proj_crs
            
            if check is True:
                if lineLayer.featureCount() > 1:
                    self.iface.messageBar().pushMessage("More than 1 feature in line layer, using the first feature!",level=Qgis.Warning)
                
                self.mapTool = KPTool(self.iface, lineLayer) #initializes the copytool class
                self.mapTool.removeVertexMarker()
                self.mapTool.kpdec = self.dlgIR.KP_prec.value() #pass on variable for decimals for KP
                self.mapTool.dccdec = self.dlgIR.DCC_prec.value() #pass on variable for decimals for DCC
                self.mapTool.out_format = self.dlgIR.output_format_box.currentIndex() #pass on variable for decimals for DCC
                self.mapTool.lineLayer=lineLayer #push the line layer to mapTool
                self.canvas.setMapTool(self.mapTool) #tell QGIS that we use our tool
                
                #self.iface.messageBar().pushMessage("Processing " + lineLayer.name() ,duration=1)
            else:
                self.iface.messageBar().pushMessage("Line and project CRS do not match (and they have to)!",level=Qgis.Warning)
                self.iface.actionPan().trigger() #trigger pan tool to prevent user working with old layer if already loaded
            
      
    def Kp4Points(self):
        if self.first_startKp4p == True:
            self.first_startKp4p = False
            self.dlgKp4p = KpFindDialogKp4p()
            self.dlgKp4p.inputLineLayerKp4p.setFilters(QgsMapLayerProxyModel.LineLayer) #allows only line layers in the box
            self.dlgKp4p.inputPtLayerKp4p.setFilters(QgsMapLayerProxyModel.PointLayer) #allows only point layers in the box
            self.dlgKp4p.setWindowFlags(Qt.WindowStaysOnTopHint) #make window always stay on top
 
        # show the dialog

        self.dlgKp4p.show()
        # Run the dialog event loop
        result = self.dlgKp4p.exec_()
        # see if OK was pressed
        if result:
            ptLayer = self.dlgKp4p.inputPtLayerKp4p.currentLayer()
            lnLayer = self.dlgKp4p.inputLineLayerKp4p.currentLayer()
            check = ptLayer.crs() == lnLayer.crs()
            if check is True:
                if lnLayer.featureCount() > 1:
                    self.iface.messageBar().pushMessage("More than 1 feature in line layer, using the first feature!",level=Qgis.Warning)
                self.mapTool = KPTool(self.iface, lnLayer)
                self.mapTool.kpdeckp4p = self.dlgKp4p.KP_prec_kp4p.value() #pass on precision values
                self.mapTool.dccdeckp4p = self.dlgKp4p.DCC_prec_kp4p.value() #pass on precision values
                new_ptLayer =self.mapTool.KP_Iterate_pts(lnLayer, ptLayer) #function to iterate over points, calls other functions in it
                QgsProject.instance().addMapLayer(new_ptLayer) #add new layer to map
                self.iface.messageBar().pushMessage("Processed points along layer " + str(lnLayer.name()),level=Qgis.Info, duration = 2)
            else:
                self.iface.messageBar().pushMessage("Line and points CRS do not match (and they have to)!",level=Qgis.Warning)
                self.iface.actionPan().trigger() #trigger pan tool to prevent user working with old layer if already loaded
                

