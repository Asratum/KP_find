# KP_find
KP Find QGIS Plugin

	Measures geodetic distance along a line, i.e. the chainage (in KP - kilometer points) and offset. 
	Clicking copies the information to the clipboard. Can also make a copy of a point layer with new attributes for KP and DOL, Lat and Lon.
	Please note: the plugin only accepts layers with the same CRS. It takes the first feature in the line layer - 
	it is best if the line layer has only one feature, but the plugin will issue a warning if not the case.
	Dependent on geographiclib (included in plugin itself).
	Code borrowed from the Lat Lon Tools, Shapetools and Closest Point plugins.
	Reliable results only with projected layers!

	3.0.4 - Added reverse KP an custom KP offset, proper warning when using geographical CRS
	3.0.3 - bug fix, lat/lon added as attribute to point layer measurements
	3.0.2 - KP measurements extend past the end points of the lines
	3.0.1 - Dialog windows stay on top, fixed a bug where you keep picking a line if you get error
	3.0.0 - Initial QGIS release
