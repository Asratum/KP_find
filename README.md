# KP_find
KP Find QGIS Plugin

	Measures geodetic distance along a line, i.e. the chainage (in KP - kilometer points) and offset. 
	Clicking copies the information to the clipboard. Can generate chainage points, geodetic and cartesian.
	Can also make a copy of a point layer with new attributes for KP and distance-off-line (DOL), Lat and Lon.
	Please note: perpendicularity changes with each CRS and the distortions it creates - the chainage points are
	those perpendicular in the CRS of the line layer.
	If more than one feature in the line layer is present, the user will be prompted to choose one by feature ID.
	
	3.2.0 - Bugfixes, handles differences in CRS in canvas and layers, option to select which feature is measured
	3.1.1 - Additional output option added
	3.1.0 - Code rewrite, added geodetic and Cartesian chainage along lines, bug fixes
	3.0.5 - Added support for outputting Cartesian distances instead of geodetic ones
	3.0.4 - Added reverse KP and custom KP offset, proper warning when using geographical CRS
	3.0.3 - bug fix, lat/lon added as attribute to point layer measurements
	3.0.2 - KP measurements extend past the end points of the lines
	3.0.1 - Dialog windows stay on top, fixed a bug where you keep picking a line if you get an error
	3.0.0 - Initial QGIS release
