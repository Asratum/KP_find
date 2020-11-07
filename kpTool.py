from geographiclib.geodesic import Geodesic #extra library included with plugin
from qgis.PyQt.QtCore import Qt, pyqtSignal, QVariant
from qgis.PyQt.QtGui import QColor
from qgis.PyQt.QtWidgets import QApplication
from qgis.core import (Qgis, 
                        QgsCoordinateTransform, 
                        QgsPointXY, 
                        QgsProject, 
                        QgsSettings, 
                        QgsDistanceArea, 
                        QgsCoordinateReferenceSystem, 
                        QgsFeature,
                        QgsGeometry,
                        QgsPoint,
                        QgsVectorLayer,
                        QgsField)
from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker

import math

from .KP_Find_dialog_IR import KpFindDialogIR

class KPTool(QgsMapToolEmitPoint):
    '''Class to interact with the map canvas to capture the coordinate
    when the mouse button is pressed and to display the coordinate in
    in the status bar.
    It will take all the things from Class QgsMapToolEmitPoint and
    overwrite some of the functions, activate, etc..
    '''

    captureStopped = pyqtSignal()
    epsg4326 = QgsCoordinateReferenceSystem("EPSG:4326")
    geod = Geodesic.WGS84

    
    def __init__(self, iface, linelayer=None):
        QgsMapToolEmitPoint.__init__(self, iface.mapCanvas())
        self.iface = iface
        self.canvas = iface.mapCanvas()
        self.marker = None
        self.vertex = None
        self.linelayer=linelayer #import a line layer from the main function
        self.proj_crs = iface.mapCanvas().mapSettings().destinationCrs()
    
    def __str__(self, linelayer):
        return linelayer.name()
    
    def activate(self):
        '''When activated set the cursor to a crosshair.'''
        self.canvas.setCursor(Qt.CrossCursor)
        self.snapcolor = QgsSettings().value( "/qgis/digitizing/snap_color" , QColor( Qt.magenta ) )

    def deactivate(self):
        self.removeMarker()
        self.removeVertexMarker()
        self.captureStopped.emit()

    def canvasMoveEvent(self, event):
        '''Capture the coordinate as the user moves the mouse over
        the canvas. Show it in the status bar.'''
        pt = self.snappoint(event.originalPixelPoint()) # input is QPoint, cursor snaps


    def snappoint(self, qpoint):
        match = self.canvas.snappingUtils().snapToMap(qpoint)
        if match.isValid():
            if self.vertex is None:
                self.vertex = QgsVertexMarker(self.canvas)
                self.vertex.setIconSize(12)
                self.vertex.setPenWidth(2)
                self.vertex.setColor(self.snapcolor)
                self.vertex.setIconType(QgsVertexMarker.ICON_BOX)
            self.vertex.setCenter(match.point())
            return (match.point()) # Returns QgsPointXY
        else:
            self.removeVertexMarker()
            return self.toMapCoordinates(qpoint) # QPoint input, returns QgsPointXY


    def canvasReleaseEvent(self, event):
        '''Capture the coordinate when the mouse button has been released,
        format it, and copy it to the clipboard. pt is QgsPointXY'''
        pt = self.snappoint(event.originalPixelPoint())
        self.removeVertexMarker()
        #if we enable this we will get some vertex markers:
#        if self.marker is None:
#            self.marker = QgsVertexMarker(self.canvas)
#            self.marker.setIconSize(18)
#            self.marker.setPenWidth(2)
#            self.marker.setIconType(QgsVertexMarker.ICON_CROSS)
#        self.marker.setCenter(pt)
        dist_off_L = self.closestPt(self.linelayer, pt, True)[1]
        line_length = self.measureLineGeod(self.linelayer, pt)
        whole_line_length = self.measureWholeLineGeod(self.linelayer)
        if self.reversekp_status == 0:            
            LL_message = '{:.{prec}f}'.format(line_length+self.offset, prec=self.kpdec) #round to specified decimal
        elif self.reversekp_status == 2:
            LL_message = '{:.{prec}f}'.format(whole_line_length-line_length+self.offset, prec=self.kpdec) #round to specified decimal
        DT_message = '{:.{prec}f}'.format(dist_off_L, prec=self.dccdec) #round to specified decimal
        #check for output format below:
        if self.out_format == KpFindDialogIR.KP_out:  # KP
            msg = "KP " + str(LL_message)
        elif self.out_format == KpFindDialogIR.KP_DCC_Out:  # KP and DCC
            msg = "KP " + str(LL_message) + " , DOL : " + str(DT_message) +" m"
        elif self.out_format == KpFindDialogIR.DMS_out:  # Lat Lon and KP
            msg = self.formatCoord(pt)[0] + ", "+ self.formatCoord(pt)[1] + " (KP " + str(LL_message) + ")"
        if msg is not None:
            clipboard = QApplication.clipboard()
            clipboard.setText(msg)
            self.iface.messageBar().pushMessage(msg +" on layer " + str(self.linelayer.name())+ " copied to clipboard",level=Qgis.Info, duration=2)            
        else:
            self.iface.messageBar().pushMessage("Something went wrong with the coordinate composition",level=Qgis.Info, duration=2)


    def removeMarker(self):
        if self.marker is not None:
            self.canvas.scene().removeItem(self.marker)
            self.marker = None

    def removeVertexMarker(self):
        if self.vertex is not None:
            self.canvas.scene().removeItem(self.vertex)
            self.vertex = None
            
    def magnitude(self, p1, p2): #not used as of current implementation
        if p1==p2: return 0
        else:
            vect_x = p2.x() - p1.x()
            vect_y = p2.y() - p1.y()
            return math.sqrt(vect_x**2 + vect_y**2)
            
    def name_collision_checker(self, new_names, old_layer):
        #check for name collisions in attribute fields
        fields_ptLayer = old_layer.dataProvider().fields()
        old_fields_names = [] #we will use this later to compare for collisions
        for f in fields_ptLayer: #iterate over all field names and add to new provider
            old_fields_names.append(f.name()) #create the collision test list
        trynumber = 1
        num_suffix = ''
        max_num = 1      
        for name in new_names:
            collision = True
            while collision:   # Iterate until there are no collissions
                collision = False
                if name+num_suffix in old_fields_names:
                    collision = True
                    num_suffix = '_'+str(trynumber)
                    if trynumber > max_num+1:
                        max_num = trynumber 
                        num_suffix = '_'+str(max_num)
                    trynumber = trynumber + 1
        new_list = [item + num_suffix for item in new_names] #add same suffix to all fields
        return new_list
        
    def closestPt(self, linelayer, clicked_point, lineext=False): #get the linelayer and QgsPointXY
        '''This will take a point and a line layer and find the closest point
        along a perpendicular line from the point to the line layer. It takes
        the first feature in that line layer. It will also extend the first and
        last line segments.'''
        line_feat = linelayer.getFeature(0) #get the only feature in the layer
        geomL = line_feat.geometry() #get the linear geometry from the feature
        if lineext == True:
            geomL = geomL.extendLine(1000000,10000000) #extend fist and last seg by 1000km
        distinit,mindistpt,afterveretexinit,leftoff=geomL.closestSegmentWithContext(clicked_point) #get min distance point to line
        ProjPoint = QgsPointXY(mindistpt[0],mindistpt[1]) #create projected point on line
        
        Distance = self.magnitude(clicked_point, ProjPoint) #measure planar distance, currently not used
        distancef = QgsDistanceArea() #define geodetic measurement function
        distancef.setEllipsoid('WGS84') #set ellipsoid
        distancef.setSourceCrs(self.proj_crs,QgsProject.instance().transformContext()) #set the CRS for measurement
        Distance_geod=distancef.measureLine(clicked_point, ProjPoint) #give geodetic distance
        #we'll try to implement geographiclib here
        srcCRS = linelayer.sourceCrs() #get CRS from line
        wgs84=KPTool.epsg4326 #define EPSG 43226
        if srcCRS != wgs84:
            geomTo4326 = QgsCoordinateTransform(srcCRS, wgs84, QgsProject.instance()) #convert if needed
            ptCP84 = geomTo4326.transform(clicked_point)
            ptPP84 = geomTo4326.transform(ProjPoint)
        else:
            ptCP84 = clicked_point
            ptPP84 = ProjPoint
        geod_ds = KPTool.geod.Inverse(ptCP84.y(), ptCP84.x(), ptPP84.y(), ptPP84.x())
        Distance_geod2 = geod_ds['s12']
        if leftoff < 0:
            Distance_geod = -Distance_geod
            Distance = -Distance
            Distance_geod2 = -Distance_geod2
            
        return ProjPoint, Distance_geod2
        
    def measureLineGeod(self, linetoMeasure, clicked_pt):
        '''This will take a line segment and the clicked point and return the
        length up till that point for the segment. Inspired by shape tools and closest point plugins'''
        srcCRS = linetoMeasure.sourceCrs() #get CRS from line
        feature = linetoMeasure.getFeature(0) #get first feature from line
        wgs84=KPTool.epsg4326 #define EPSG 4326
        proj_pt, dist_ext = self.closestPt(linetoMeasure, clicked_pt, True) #we'll need that point later to compare with projected points of subsegments
        proj_pt_on_line, dist_nonext = self.closestPt(linetoMeasure, clicked_pt, False) #get proj point on non-extended line
        geomF = feature.geometry()
        if srcCRS != wgs84:
            geomTo4326 = QgsCoordinateTransform(srcCRS, wgs84, QgsProject.instance()) #convert if needed
        if geomF.isMultipart(): #get nodes out of data
            ptdata = [geomF.asMultiPolyline()]
        else:
            ptdata = [[geomF.asPolyline()]]
        for seg in ptdata:
            if len(seg) < 1: #should never happen
                self.iface.messageBar().pushMessage("Something is strange with your line file",level=Qgis.Critical, duration=2)
                continue 
            for pts in seg: #now we get all nodes from start to end of line
                numpoints = len(pts)
                
                if numpoints < 2: #should never happen
                    self.iface.messageBar().pushMessage("Something is strange with your line file",level=Qgis.Critical, duration=2)
                    continue
                    
                ptStart = QgsPointXY(pts[0].x(), pts[0].y()) #get initial point coords
                # Calculate the total distance of this line segment
                distance = 0.0
                for x in range(1,numpoints): #from point one (since we preextend segment), cumulative
                    ptEnd = QgsPointXY(pts[x].x(), pts[x].y()) #coordinates of next point
                    PP = QgsPoint(ptStart.x(), ptStart.y())
                    PPP = QgsPoint(ptEnd.x(), ptEnd.y())
                    seg_geom = QgsGeometry.fromPolyline((PP, PPP)) #generate geometry for current subsegment
                    distinit,mindistpt,aftervertexinit,leftoff = seg_geom.closestSegmentWithContext(proj_pt)
                    TestPoint=QgsPointXY(mindistpt[0], mindistpt[1]) #create projected point on subsegment
                    
                    if srcCRS != wgs84: # Convert to 4326 - just to be safe when using geodetic measure
                        ptStart84 = geomTo4326.transform(ptStart)
                        ptEnd84 = geomTo4326.transform(ptEnd)
                        proj_pt84 = geomTo4326.transform(proj_pt)
                        proj_pt_on_line84 = geomTo4326.transform(proj_pt_on_line)
                        TestPoint84 = geomTo4326.transform(TestPoint)
                    else:
                        ptStart84 = ptStart
                        ptEnd84 = ptEnd
                        proj_pt84 = proj_pt
                        proj_pt_on_line84 = proj_pt_on_line
                        TestPoint84 = TestPoint
                    test_dist = KPTool.geod.Inverse(proj_pt_on_line84.y(), proj_pt_on_line84.x(), TestPoint84.y(), TestPoint84.x()) #check distance between two on-subsegment-points
                    l = KPTool.geod.Inverse(ptStart84.y(), ptStart84.x(), ptEnd84.y(), ptEnd84.x())  #geodetic distance between begin and end of subsegment
                    roundTD = round(test_dist['s12'], 6) #round so we can get zero
                    
                    if roundTD != 0: #not yet the last segment to be measured
                        distance += l['s12'] #add to comulative
                    else: #this is the segment, where we have to stop measuring
                        d_to_click = KPTool.geod.Inverse(ptStart84.y(), ptStart84.x(), proj_pt_on_line84.y(), proj_pt_on_line84.x())
                        distance += d_to_click['s12']
                        break #break loop and give distance so far
                    ptStart = ptEnd #the last shall be the first and the loop goes on
                if QgsPoint(proj_pt) != QgsPoint(proj_pt_on_line):
                    extra_dist = KPTool.geod.Inverse(proj_pt84.y(), proj_pt84.x(), proj_pt_on_line84.y(), proj_pt_on_line84.x())
                    if round(distance, 6) == 0:
                        distance = -extra_dist['s12']
                    else:
                        distance = distance + extra_dist['s12']
                out = distance/1000 # Distance converted KM
        
        return out
        
    def KP_Iterate_pts(self, linetoMeasurekp4p, ptLayer):
        """will iterate over all points in a layer, finding distance to and along line.
        Outputs a new layer. Inspired by shape tools and closest point plugins"""
        new_pt_layer = QgsVectorLayer("Point", "KPed_"+ str(ptLayer.name()), "memory") #define new layer we will return
        new_pt_layer.setCrs(ptLayer.crs()) #get crs from input layer
        prov_old = ptLayer.dataProvider() #provider to get the attribute field names and type
        provider_ptLayer = new_pt_layer.dataProvider() #provider for new layer to add the features to
        fields_ptLayer = prov_old.fields()
        whole_line_length = self.measureWholeLineGeod(linetoMeasurekp4p)
        for f in fields_ptLayer: #iterate over all field names and add to new provider
            znameField= f.name()
            Type= str(f.typeName())
            if Type == 'Integer': provider_ptLayer.addAttributes([ QgsField( znameField, QVariant.Int)])
            if Type == 'Real': provider_ptLayer.addAttributes([ QgsField( znameField, QVariant.Double)])
            if Type == 'String': provider_ptLayer.addAttributes([ QgsField( znameField, QVariant.String)])
            else : provider_ptLayer.addAttributes([ QgsField( znameField, QVariant.String)])
        iterate_names_list = ["KP","DOL","Latitude","Longitude"]
        new_attr_names = self.name_collision_checker(iterate_names_list, ptLayer)
        provider_ptLayer.addAttributes([QgsField(new_attr_names[0], QVariant.Double),
                                    QgsField(new_attr_names[1], QVariant.Double), 
                                    QgsField(new_attr_names[2], QVariant.String),
                                    QgsField(new_attr_names[3], QVariant.String)])  #four new fields we are calculating
        new_pt_layer.startEditing()
        
        for old_feat in ptLayer.getFeatures(): #iterate over all point features
            geom_of=old_feat.geometry()
            point_of=geom_of.asPoint() #get point object
            #geom_nf= QgsGeometry().fromPointXY(point_of) 
            new_feat = QgsFeature()
            new_feat.setGeometry(geom_of) #create and set geometry from old feature
            attributes_nf=old_feat.attributes()  #copy of old feat attributes
            if self.reversekp_statuskp4p == 0:
                kp_dist = self.measureLineGeod(linetoMeasurekp4p, point_of)+self.offsetkp4p #run measurements
            elif self.reversekp_statuskp4p == 2:
                kp_dist = whole_line_length - self.measureLineGeod(linetoMeasurekp4p, point_of)+self.offsetkp4p
            prjp, dcc_dist = self.closestPt(linetoMeasurekp4p, point_of, True) #run measurements. we won't need prjp

            attributes_nf.append(round(kp_dist, self.kpdeckp4p)) #round with precision values
            attributes_nf.append(round(dcc_dist,self.dccdeckp4p)) #round with precision values
            attributes_nf.append(self.formatCoord(point_of)[0]) #insert latitude
            attributes_nf.append(self.formatCoord(point_of)[1]) #insert longitude
            new_feat.setAttributes(attributes_nf) #set new attributes to new feat
            provider_ptLayer.addFeatures([ new_feat ]) #add new feature to provider
            
        new_pt_layer.commitChanges() #finish editing layer
        return new_pt_layer #return the new layer with new features with new attributes
        
    def formatCoord(self, pt):
        canvasCRS = self.canvas.mapSettings().destinationCrs()
        epsg4326 = QgsCoordinateReferenceSystem('EPSG:4326')
        #convert point to wgs84 for conversion
        if canvasCRS == epsg4326:
            pt4326 = pt
        else:
            transform = QgsCoordinateTransform(canvasCRS, epsg4326, QgsProject.instance())
            pt4326 = transform.transform(pt.x(), pt.y())

        lat = self.convertDD2DM(pt4326.y(), True, 4)
        lon = self.convertDD2DM(pt4326.x(), False, 4)
        return [lat,lon]
        
        
    def convertDD2DM(self, coord, islat, prec):
        '''Convert decimal degrees to DM - taken from latlontools plugin'''
        if islat:
            if coord < 0:
                unit = 'S'
            else:
                unit = 'N'
        else:
            if coord > 0:
                unit = 'E'
            else:
                unit = 'W'
        dmsSpace = " " #put some spaces in there
        zeroes = 1 #this will be used for padding
        coord = math.fabs(coord)
        deg = math.floor(coord)
        dmin = (coord - deg) * 60.0
        min = math.floor(dmin)
        sec = (dmin - min) * 60.0
        s = ""

        # Properly handle rounding based on the digit precision
        d = "{:.{prec}f}".format(dmin, prec=prec)
        if float(d) == 60:
            deg += 1
            dmin = 0
        if islat:
            s = "{:0{}.0f}\xB0{}{:0{}.0{prec}f}\'{}{}".format(deg, zeroes*2, dmsSpace, dmin, prec+zeroes*3, dmsSpace, unit, prec=prec)
        else:
            s = "{:0{}.0f}\xB0{}{:0{}.0{prec}f}\'{}{}".format(deg, zeroes*3, dmsSpace, dmin, prec+zeroes*3, dmsSpace, unit, prec=prec)
        return(s)
        
    def measureWholeLineGeod(self, linetoMeasure):
        '''This will take a line segment and return its geodetic length'''
        srcCRS = linetoMeasure.sourceCrs() #get CRS from line
        feature = linetoMeasure.getFeature(0) #get first feature from line
        wgs84=KPTool.epsg4326 #define EPSG 43226
        geomF = feature.geometry()
        if srcCRS != wgs84:
            geomTo4326 = QgsCoordinateTransform(srcCRS, wgs84, QgsProject.instance()) #convert if needed
        if geomF.isMultipart(): #get nodes out of data
            ptdata = [geomF.asMultiPolyline()]
        else:
            ptdata = [[geomF.asPolyline()]]
        for seg in ptdata:
            if len(seg) < 1: #should never happen
                self.iface.messageBar().pushMessage("Something is strange with your line file",level=Qgis.Critical, duration=2)
                continue 
            for pts in seg: #now we get all nodes from start to end of line
                numpoints = len(pts)
                
                if numpoints < 2: #should never happen
                    self.iface.messageBar().pushMessage("Something is strange with your line file",level=Qgis.Critical, duration=2)
                    continue
                    
                ptStart = QgsPointXY(pts[0].x(), pts[0].y()) #get initial point coords
                # Calculate the total distance of this line segment
                distance = 0.0
                for x in range(1,numpoints): #from point one (since we preextend segment), cumulative
                    ptEnd = QgsPointXY(pts[x].x(), pts[x].y()) #coordinates of next point                 
                    if srcCRS != wgs84: # Convert to 4326 - just to be safe when using geodetic measure
                        ptStart84 = geomTo4326.transform(ptStart)
                        ptEnd84 = geomTo4326.transform(ptEnd)
                    else:
                        ptStart84 = ptStart
                        ptEnd84 = ptEnd
                    l = KPTool.geod.Inverse(ptStart84.y(), ptStart84.x(), ptEnd84.y(), ptEnd84.x())  #geodetic distance between begin and end of subsegment                    
                    distance += l['s12'] #add to comulative

                    ptStart = ptEnd #the last shall be the first and the loop goes on
                out = distance/1000 # Distance converted KM
        
        return out