# KP_find
KP Find QGIS Plugin

	Measures geodetic distance along a line, i.e. the chainage (in KP - kilometer points) and offset. 
	Clicking copies the information to the clipboard. Can also make a copy of a point layer with new attributes for KP and DOL.
	Please note: the plugin only accepts layers with the same CRS. It takes the first feature in the line layer - 
	it is best if the line layer has only one feature, but the plugin will issue a warning if not the case.
	Dependent on geographiclib (included in plugin itself).
	Code borrowed from the Lat Lon Tools, Shapetools and Closest Point plugins.
	Reliable results only with projected layers!
	Tested on QGIS version 3.10.2.
