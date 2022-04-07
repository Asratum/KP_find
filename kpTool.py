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
                        QgsField,
                        QgsMultiLineString,
                        QgsLineString)
from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker
import math
from .KP_Find_dialog_IR import KpFindDialogInteractive

class KPTool(QgsMapToolEmitPoint):
    """Class to interact with the map canvas to capture the coordinate
    when the mouse button is pressed and to display the coordinate in
    in the status bar.
    It will take all the things from Class QgsMapToolEmitPoint and
    overwrite some of the functions, activate, etc..
    Other functions here do the measurements
    """
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
        the canvas. Show it in the status bar. Currently not used'''
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
        """Capture the coordinate when the mouse button has been released,
        format it, and copy it to the clipboard. pt is QgsPointXY"""
        pt = self.snappoint(event.originalPixelPoint())
        self.removeVertexMarker()
        # if we enable this we will get some vertex markers:
        # if self.marker is None:
        #     self.marker = QgsVertexMarker(self.canvas)
        #     self.marker.setIconSize(18)
        #     self.marker.setPenWidth(2)
        #     self.marker.setIconType(QgsVertexMarker.ICON_CROSS)
        # self.marker.setCenter(pt)        
        if self.geodetic_use == 0:
            # use cartesian
            line_length = self.measureLine(self.linelayer, pt, False)
            whole_line_length = self.measureWholeLine(self.linelayer, False)
            dist_off_Line = self.closestPt(self.linelayer, pt, True, False)[1]  #second False means we use Cartesian distance
        elif self.geodetic_use == 2:
            # use geodetic
            line_length = self.measureLine(self.linelayer, pt, True)
            whole_line_length = self.measureWholeLine(self.linelayer, True)
            dist_off_Line = self.closestPt(self.linelayer, pt, True, True)[1]
        if self.reversekp_status == 0:            
            LL_message = '{:.{prec}f}'.format(line_length+self.offset, prec=self.kpdec) #round to specified decimal
        elif self.reversekp_status == 2:
            LL_message = '{:.{prec}f}'.format(whole_line_length-line_length+self.offset, prec=self.kpdec) #round to specified decimal
        DT_message = '{:.{prec}f}'.format(dist_off_Line, prec=self.dccdec) #round to specified decimal
        # check for output format below:
        if self.out_format == KpFindDialogInteractive.KP_out:  # KP
            msg = 'KP ' + str(LL_message)
        elif self.out_format == KpFindDialogInteractive.KP_DCC_Out:  # KP and DCC
            msg = 'KP ' + str(LL_message) + ' , DOL : ' + str(DT_message) +' m'
        elif self.out_format == KpFindDialogInteractive.DMS_out:  # Lat Lon and KP
            msg = self.formatCoord(pt)[0] + ', '+ self.formatCoord(pt)[1] + ' (KP ' + str(LL_message) + ')'
        if msg is not None:
            clipboard = QApplication.clipboard()
            clipboard.setText(msg)
            self.iface.messageBar().pushMessage(msg + ' on layer ' + str(self.linelayer.name()) + ' copied to clipboard', level=Qgis.Info, duration=2)            
        else:
            self.iface.messageBar().pushMessage('Something went wrong with the coordinate composition', level=Qgis.Info, duration=2)

    def removeMarker(self):
        if self.marker is not None:
            self.canvas.scene().removeItem(self.marker)
            self.marker = None

    def removeVertexMarker(self):
        if self.vertex is not None:
            self.canvas.scene().removeItem(self.vertex)
            self.vertex = None
                        
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
        
    def closestPt(self, linelayer, clicked_point, extend_line=False, give_geod=True): # get the linelayer and QgsPointXY
        """This will take a point and a line layer and find the closest point
        along a perpendicular line from the point to the line layer. It takes
        the first feature in that line layer. It will also extend the first and
        last line segments. Cartesian and geodetic distance possible"""
        if linelayer.getFeature(0).geometry().isMultipart():
            line_feat = linelayer.getFeature(0)
            linear_geom = line_feat.geometry()
        else:
            line_feat = linelayer.getFeature(1) #for some reason linestring feature 0 is NULL
            linear_geom = line_feat.geometry()
            linear_geom.convertToMultiType()
    
        if extend_line == True:
            linear_geom = linear_geom.extendLine(1000000, 10000000) # extend fist and last seg by 1000km, bit of a fudge
        _, mindistpt, _, leftoff = linear_geom.closestSegmentWithContext(clicked_point) # get min distance point to line
        projected_point = QgsPointXY(mindistpt[0], mindistpt[1]) # create projected point on line        
        distance_cart = QgsDistanceArea().measureLine(clicked_point, projected_point) #cartesian distance
        srcCRS = linelayer.sourceCrs() # get CRS from line
        wgs84 = KPTool.epsg4326 # define EPSG 4326
        if srcCRS != wgs84:
            geomTo4326 = QgsCoordinateTransform(srcCRS, wgs84, QgsProject.instance()) #convert if needed
            ptCP84 = geomTo4326.transform(clicked_point)
            ptPP84 = geomTo4326.transform(projected_point)
        else:
            ptCP84 = clicked_point
            ptPP84 = projected_point
        geod_ds = KPTool.geod.Inverse(ptCP84.y(), ptCP84.x(), ptPP84.y(), ptPP84.x()) # use geographiclib
        distance_geod = geod_ds['s12']
        if leftoff < 0:
            distance_cart = -distance_cart #planar distance
            distance_geod = -distance_geod #return the geographiclib distance
        #if abs(distance_cart)>1000000 or abs(distance_geod)>1000000:
            #self.iface.messageBar().pushMessage('Line extended more than 1000 km, possible errors',level=Qgis.Warning) #warn user if they click more than 1000 km away from end of line
        if give_geod is True:    
            return projected_point, distance_geod
        else:
            return projected_point, distance_cart
        
    def measureLine(self, linetoMeasure, clicked_pt, geodetic_measure=True):
        """This will take a line segment and the clicked point and return the
        length up till that point for the segment in km. Inspired by shape tools and closest point plugins"""
        srcCRS = linetoMeasure.sourceCrs() #get CRS from line
        feature = linetoMeasure.getFeature(0) #get first feature from line
        wgs84 = KPTool.epsg4326 #define EPSG 4326
        projected_point_raw, dist_ext = self.closestPt(linetoMeasure, clicked_pt, True) #we'll need that point later to compare with projected points of subsegments
        projected_point_on_line_raw, dist_nonext = self.closestPt(linetoMeasure, clicked_pt, False) #get proj point on non-extended line
        if srcCRS!=wgs84:
            geomTo4326 = QgsCoordinateTransform(srcCRS, wgs84, QgsProject.instance()) #convert if needed
        if feature.geometry().isMultipart():
            # get nodes out of data
            ptdata = [feature.geometry().asMultiPolyline()]
        else:
            feature = linetoMeasure.getFeature(1) #for some reason singleline first feature is NULL
            ptdata = [[feature.geometry().asPolyline()]]
        for seg in ptdata:
            if len(seg) < 1: 
                # should never happen
                self.iface.messageBar().pushMessage('Something is strange with your line file',level=Qgis.Critical, duration=2)
                continue 
            for pts in seg:
                # now we get all nodes from start to end of segment line
                numpoints = len(pts)                
                if numpoints < 2:
                    # should never happen
                    self.iface.messageBar().pushMessage('Something is strange with your line file',level=Qgis.Critical, duration=2)
                    continue                    
                ptStart_raw = QgsPointXY(pts[0].x(), pts[0].y()) #get initial point coords
                # Calculate the total distance of this line segment
                distance = 0.0
                for x in range(1,numpoints):
                    # from point one (since we preextend segment) start measuring distance cumulative
                    ptEnd_raw = QgsPointXY(pts[x].x(), pts[x].y()) # coordinates of next point
                    P1 = QgsPoint(ptStart_raw.x(), ptStart_raw.y())
                    P2 = QgsPoint(ptEnd_raw.x(), ptEnd_raw.y())
                    seg_geom = QgsGeometry.fromPolyline((P1, P2)) # generate geometry for current subsegment
                    _, mindistpt, _, _ = seg_geom.closestSegmentWithContext(projected_point_raw)
                    TestPoint = QgsPointXY(mindistpt[0], mindistpt[1]) #create projected point on subsegment 
                    ptStart = ptStart_raw
                    ptEnd = ptEnd_raw
                    projected_point = projected_point_raw
                    projected_point_on_line = projected_point_on_line_raw
                    if geodetic_measure==True:
                        if srcCRS != wgs84:
                            # Convert to 4326 - just to be safe when using geodetic measure
                            ptStart = geomTo4326.transform(ptStart_raw)
                            ptEnd = geomTo4326.transform(ptEnd_raw)
                            projected_point = geomTo4326.transform(projected_point_raw)
                            projected_point_on_line = geomTo4326.transform(projected_point_on_line_raw)
                            TestPoint = geomTo4326.transform(TestPoint)
                        d_to_click = KPTool.geod.Inverse(ptStart.y(), ptStart.x(), projected_point_on_line.y(), projected_point_on_line.x())['s12']
                        len_subsegment = KPTool.geod.Inverse(ptStart.y(), ptStart.x(), ptEnd.y(), ptEnd.x())['s12']  # geodetic distance between begin and end of subsegment
                        test_distance = KPTool.geod.Inverse(projected_point_on_line.y(), projected_point_on_line.x(), TestPoint.y(), TestPoint.x())['s12'] # check distance between two on-subsegment-points
                        
                    else:
                        d_to_click = QgsDistanceArea().measureLine(ptStart, projected_point_on_line)
                        len_subsegment = QgsDistanceArea().measureLine(ptStart, ptEnd) # cartesian distance between begin and end point
                        test_distance = QgsDistanceArea().measureLine(projected_point_on_line, TestPoint) # check distance between two on-subsegment-points                                                
                    round_test_distance = round(test_distance, 6) # round so we can get zero, otherwise it will be a small fraction    
                    if round_test_distance != 0:
                        # not yet the last segment to be measured
                        distance += len_subsegment #add to comulative
                    else:
                        # this is the segment, where we have to stop measuring
                        distance += d_to_click
                        break # break loop and give distance so far
                    ptStart_raw = ptEnd_raw # make ready for next pair of points                    
                if QgsPoint(projected_point) != QgsPoint(projected_point_on_line):
                    # if our point on the extended line is different, the projected point is before start or after end of line layer
                    if geodetic_measure==True:
                        extra_dist = KPTool.geod.Inverse(projected_point.y(), projected_point.x(), projected_point_on_line.y(), projected_point_on_line.x())['s12']
                    else:
                        extra_dist = QgsDistanceArea().measureLine(projected_point, projected_point_on_line)
                    if round(distance, 6) == 0:
                        # we are at the start of the line layer, the for loops hasn't done anything, so negative distance to starting node
                        distance = -extra_dist
                    else:
                        # we are at the end of the line segment, so need to add more distance
                        distance = distance + extra_dist
                out = distance/1000 # Distance converted to KM        
        return out    
    
    def kpIteratePts(self, linetoMeasurekp4p, ptLayer):
        """will iterate over all points in a layer, finding distance to and along line.
        Outputs a new layer."""
        new_pt_layer = QgsVectorLayer('Point', 'KPed_'+ str(ptLayer.name()), 'memory') # define new layer we will return
        new_pt_layer.setCrs(ptLayer.crs()) # get crs from input layer
        prov_old = ptLayer.dataProvider() # provider to get the attribute field names and type
        provider_ptLayer = new_pt_layer.dataProvider() # provider for new layer to add the features to
        fields_ptLayer = prov_old.fields()
        for f in fields_ptLayer:
            # iterate over all field names and add to new provider
            znameField = f.name()
            type_field = str(f.typeName())
            if type_field == 'Integer': provider_ptLayer.addAttributes([ QgsField(znameField, QVariant.Int)])
            if type_field == 'Real': provider_ptLayer.addAttributes([ QgsField(znameField, QVariant.Double)])
            if type_field == 'String': provider_ptLayer.addAttributes([ QgsField(znameField, QVariant.String)])
            else: provider_ptLayer.addAttributes([ QgsField( znameField, QVariant.String)])

        iterate_names_list = ["KP","DOL","Latitude","Longitude"]
        new_attr_names = self.name_collision_checker(iterate_names_list, ptLayer) # check if the new names we add already exist
        provider_ptLayer.addAttributes([QgsField(new_attr_names[0], QVariant.Double),
                                    QgsField(new_attr_names[1], QVariant.Double), 
                                    QgsField(new_attr_names[2], QVariant.String),
                                    QgsField(new_attr_names[3], QVariant.String)])  # four new fields we are calculating
        new_pt_layer.startEditing()        
        for old_feat in ptLayer.getFeatures():
            # iterate over all point features
            geom_of = old_feat.geometry()
            point_of = geom_of.asPoint() # get point object
            new_feat = QgsFeature()
            new_feat.setGeometry(geom_of) # create and set geometry from old feature
            attributes_nf = old_feat.attributes()  # copy of old feat attributes
            
            if self.geodetic_usekp4p == 0:
                # we measure Cartesian coords
                line_length = self.measureLine(self.linelayer, point_of, False) + self.offsetkp4p
                whole_line_length = self.measureWholeLine(self.linelayer, False)
                dcc_dist = self.closestPt(linetoMeasurekp4p, point_of, True, False)[1]  # second False means we use Cartesian distance
            elif self.geodetic_usekp4p == 2:
                # we measure geodetic coords
                line_length = self.measureLine(linetoMeasurekp4p, point_of, True) + self.offsetkp4p
                whole_line_length = self.measureWholeLine(self.linelayer, True)
                dcc_dist = self.closestPt(linetoMeasurekp4p, point_of, True, True)[1]
            if self.reversekp_statuskp4p == 0:
                kp_dist =  line_length
            elif self.reversekp_statuskp4p == 2:
                kp_dist = whole_line_length - line_length

            attributes_nf.append(round(kp_dist, self.kpdeckp4p)) # round with precision values
            attributes_nf.append(round(dcc_dist,self.dccdeckp4p)) # round with precision values
            attributes_nf.append(self.formatCoord(point_of)[0]) # insert latitude
            attributes_nf.append(self.formatCoord(point_of)[1]) # insert longitude
            new_feat.setAttributes(attributes_nf) # set new attributes to new feat
            provider_ptLayer.addFeatures([ new_feat ]) # add new feature to provider
            
        new_pt_layer.commitChanges() # finish editing layer
        return new_pt_layer # return the new layer with new features with new attributes
        
    def formatCoord(self, pt):
        canvasCRS = self.canvas.mapSettings().destinationCrs()
        epsg4326 = QgsCoordinateReferenceSystem('EPSG:4326')
        # convert point to wgs84 for conversion
        if canvasCRS == epsg4326:
            pt4326 = pt
        else:
            transform = QgsCoordinateTransform(canvasCRS, epsg4326, QgsProject.instance())
            pt4326 = transform.transform(pt.x(), pt.y())
        lat = self.convertDD2DM(pt4326.y(), True, 4)
        lon = self.convertDD2DM(pt4326.x(), False, 4)
        return [lat,lon]
                
    def convertDD2DM(self, coord, islat, prec):
        """Convert decimal degrees to DM - taken from latlontools plugin"""
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
        dmsSpace = ' ' # put some spaces in there
        zeroes = 1 # this will be used for padding
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
            s = '{:0{}.0f}\xB0{}{:0{}.0{prec}f}\'{}{}'.format(deg, zeroes*2, dmsSpace, dmin, prec+zeroes*3, dmsSpace, unit, prec=prec)
        else:
            s = '{:0{}.0f}\xB0{}{:0{}.0{prec}f}\'{}{}'.format(deg, zeroes*3, dmsSpace, dmin, prec+zeroes*3, dmsSpace, unit, prec=prec)
        return(s)
        
    def measureWholeLine(self, linetoMeasure, geodetic_measure=True):
        """This will take a line segment and return its geodetic length"""
        srcCRS = linetoMeasure.sourceCrs() #get CRS from line
        feature = linetoMeasure.getFeature(0) #get first feature from line
        wgs84=KPTool.epsg4326 #define EPSG 43226            
        if feature.geometry().isMultipart(): #get nodes out of data
            ptdata = [feature.geometry().asMultiPolyline()]
        else:
            feature = linetoMeasure.getFeature(1)
            ptdata = [[feature.geometry().asPolyline()]]
        for seg in ptdata:
            if len(seg) < 1: #should never happen
                self.iface.messageBar().pushMessage('Something is strange with your line file',level=Qgis.Critical, duration=2)
                continue 
            for pts in seg: # now we get all nodes from start to end of line
                numpoints = len(pts)
                
                if numpoints < 2: # should never happen
                    self.iface.messageBar().pushMessage('Something is strange with your line file',level=Qgis.Critical, duration=2)
                    continue                    
                ptStart_raw = QgsPointXY(pts[0].x(), pts[0].y()) # get initial point coords
                # Calculate the total distance of this line segment
                distance = 0.0
                for x in range(1,numpoints):
                    # from point one (since we preextend segment), cumulative
                    ptEnd_raw = QgsPointXY(pts[x].x(), pts[x].y()) # coordinates of next point
                    if geodetic_measure==True:
                        if srcCRS != wgs84:
                            # Convert to 4326 - just to be safe when using geodetic measure
                            geomTo4326 = QgsCoordinateTransform(srcCRS, wgs84, QgsProject.instance()) # convert if needed
                            ptStart = geomTo4326.transform(ptStart_raw)
                            ptEnd = geomTo4326.transform(ptEnd_raw)
                        length = KPTool.geod.Inverse(ptStart.y(), ptStart.x(), ptEnd.y(), ptEnd.x())['s12']  # geodetic distance between begin and end of subsegment                    
                    else:
                        length = QgsDistanceArea().measureLine(ptStart_raw, ptEnd_raw)
                    distance += length # add to cumulative
                    ptStart_raw = ptEnd_raw # make ready for next pair of points 
                out = distance/1000 # Distance converted KM        
        return out
    
    def putKPPointsAlongLine(self, source, maxseglen):
        "We travel along the line and test if each consecutive segment should contain KPs and how many"
        layercrs = source.sourceCrs()
        if layercrs != KPTool.epsg4326:
            transto4326 = QgsCoordinateTransform(layercrs, KPTool.epsg4326, QgsProject.instance())
            transfrom4326 = QgsCoordinateTransform(KPTool.epsg4326, layercrs, QgsProject.instance())
            
        new_pt_layer = QgsVectorLayer('Point', 'KP_points_'+ str(source.name()), 'memory') # define new layer we will return
        new_pt_layer.setCrs(source.crs()) # get crs from input layer
        new_pt_layer.startEditing() 
        provider_ptLayer = new_pt_layer.dataProvider() # provider for new layer to add the features to
        provider_ptLayer.addAttributes([QgsField('KP', QVariant.Double)])
        iterator = source.getFeatures()
        KP_label_count = 0
        
        for cnt, feature in enumerate(iterator):
            if feature.geometry().isMultipart():
                seg = feature.geometry().asMultiPolyline()
            else:
                seg = [feature.geometry().asPolyline()]
            numseg = len(seg)
            if numseg < 1 or len(seg[0]) < 2:
                self.iface.messageBar().pushMessage('Less than one segment in line layer',level=Qgis.Critical, duration=2)
                continue
            for line in seg:
                numpoints = len(line)
                if self.Reverse_KP_points == 2: #reverse point order for reverse KP
                    _ = line.reverse()
                ptStart = QgsPointXY(line[0][0], line[0][1])
                new_kp_point = self.createFeatureFromPoint(ptStart, KP_label_count)
                provider_ptLayer.addFeatures([ new_kp_point ])
                remaining_dist = 0.0 #remaining distance to next point
                if layercrs != KPTool.epsg4326:  # Convert to 4326
                    ptStart = transto4326.transform(ptStart)
                for x in range(1, numpoints):
                    ptEnd = QgsPointXY(line[x][0], line[x][1])
                    if layercrs != KPTool.epsg4326:  # Convert to 4326
                        ptEnd = transto4326.transform(ptEnd)
                    gline = KPTool.geod.InverseLine(ptStart.y(), ptStart.x(), ptEnd.y(), ptEnd.x())
                    if remaining_dist+gline.s13 > maxseglen: #we have to place at least one KP                       
                        s = maxseglen - remaining_dist
                        g = gline.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                        ptKP = QgsPointXY(g['lon2'], g['lat2'])
                        if layercrs != KPTool.epsg4326:
                            ptKP = transfrom4326.transform(ptKP)
                        KP_label_count = KP_label_count + maxseglen
                        remaining_dist = remaining_dist + gline.s13 - maxseglen
                        new_kp_point = self.createFeatureFromPoint(ptKP, KP_label_count)
                        provider_ptLayer.addFeatures([ new_kp_point ])
                        if remaining_dist > maxseglen: #we need to place more KP pts
                            extra_from_start = s
                            n = int(remaining_dist / maxseglen)                              
                            for i in range(0, n):
                                s = maxseglen * (i+1) + extra_from_start
                                g = gline.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                                ptKP = QgsPointXY(g['lon2'], g['lat2'])
                                if layercrs != KPTool.epsg4326:  # Convert each point back to the output CRS
                                    ptKP = transfrom4326.transform(ptKP)
                                KP_label_count = KP_label_count + maxseglen
                                remaining_dist = remaining_dist - maxseglen
                                new_kp_point = self.createFeatureFromPoint(ptKP, KP_label_count)
                                provider_ptLayer.addFeatures([ new_kp_point ])
                    else: #no KPs placed in this segment, keep the cumulative distance                         
                        remaining_dist = remaining_dist + gline.s13  
                    ptStart = ptEnd
        new_pt_layer.commitChanges()
        return new_pt_layer
    
    def putKPPointsAlongLineCart(self, source, maxseglen):
        "this uses the .interpolate method from the QgsGeometry class, which works well in Cartesian coordinates"
        new_pt_layer = QgsVectorLayer('Point', 'KP_points_'+ str(source.name()), 'memory') # define new layer we will return
        new_pt_layer.setCrs(source.crs()) # get crs from input layer
        new_pt_layer.startEditing() 
        provider_ptLayer = new_pt_layer.dataProvider() # provider for new layer to add the features to
        provider_ptLayer.addAttributes([QgsField('KP', QVariant.Double)])
        iterator = source.getFeatures()
        for cnt, feature in enumerate(iterator):
            if feature.geometry().isMultipart():
                seg = feature.geometry()
            else:
                seg = feature.geometry()
                seg.convertToMultiType()
            if self.Reverse_KP_points == 2: #reverse point order for reverse KP
                seg = self.reverseLineDirection(seg)
            n = int(seg.length() / maxseglen)
            KP_label_count = -maxseglen # so we get zero with first interp point
            for i in range (0,n+1):
                ptKP = seg.interpolate(i*maxseglen)
                ptKP = ptKP.asPoint()
                KP_label_count = KP_label_count + maxseglen
                new_kp_point = self.createFeatureFromPoint(ptKP, KP_label_count)
                provider_ptLayer.addFeatures([ new_kp_point ])
        new_pt_layer.commitChanges()
        return new_pt_layer
    
    def createFeatureFromPoint(self, point, attributes):
        "takes an QgsPointXY and makes it a feature, adds attributes"
        new_feat = QgsFeature()
        new_feat.setGeometry(QgsGeometry.fromPointXY(point)) # create and set geometry from old feature
        new_feat.setAttributes([round(attributes, 3)])
        return new_feat
        
    def reverseLineDirection(self, line_geom):
        mls1 = line_geom.get()
        mls2 = QgsMultiLineString()
        # For each reversed linestring, visited in reverse order
        for i in [QgsLineString([*i][::-1]) for i in [*mls1][::-1]]:
            _ = mls2.addGeometry(i) # add it to new geometry
        new_geometry = QgsGeometry(mls2)    
        return new_geometry
        