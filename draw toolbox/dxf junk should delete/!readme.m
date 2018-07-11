
_____________________________dxf2coord 1.1 (freeware)

description:
the dxf2coord script reads the coordinates of points, lines, polylines
(LWPOLYLINES with and without elevation) 3d polylines, 3dfaces and circles 
out of a acad r2000 - r2004 ascii .dxf file. this simple script treads the 
.dxf file as a string and looks for the keywords for the interested entity 
(e.g. 'AcDbPolyline' for polylines) and reads the following coordinates. 
for details see the dxf reference.

caution: this script won't read additional information like layer color,
linewidth, linestyle etc.!!!

everybody can feel free to expand the functionality of this script but 
please send the improved file to the mathworks file exchange that other 
interested people can use it too!
this package contains 2 different .m scripts, dxf2coord_11_cellarray.m and 
dxf2coord_11_matrix.m and this !readme.m file
independend on which of the two scripts u use,the coordinates are saved in 
the form (object_id  x-coords  y-coords  z-coords) except for circles: 
(object_id  x-coords  y-coords  z-coords  radius)

if an entity does not exist in the .dxf file the corresponding output variable 
stays empty. this script is a temporary solution because real .dxf import is 
not easy and also depends on the version of the .dxf file.(by the way i am 
too stupid to program that and would have no time to do it anyway)

advantages:
+) different entities can be included in one .dxf file and will be read out seperately
+) closed polylines are detected, so the last pair of coordinates is the same as the first
+) the script can easily be adapted e.g. as function

disadvantages:
+) although i have tried to program this script to be as fast as possible
   (or how far my programming knowlegde reaches) the script is relatively slow 
   with bigger .dxf files (e.g. bigger than 10,11 Mb, depends on your machine) 
   the resulting boringness corresponds to the size of the .dxf file 
   
   
dxf2coord_11_cellarray.m_____________________________

in this script all output variables are cell arrays  


dxf2coord_11_matrix.m_____________________________

all output variables are saved as matrices


i hope you find these scripts useful!!
lukas
   
















   