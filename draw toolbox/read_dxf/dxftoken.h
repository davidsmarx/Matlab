////////////////////////////////////////////////////////////
//
// Name:    dxf.h
//
// Author:  Steven Michael
//          (smichael@ll.mit.edu)
//
// Date:    3/10/2005
//
// Description:
//
//    The following simply lists the DXF tokens for the 
//    DXF C++ class
//
////////////////////////////////////////////////////////////

#ifndef _DXTOKEN_H_
#define _DXTOKEN_H_

typedef enum {
  BAD_TOKEN,
  SECTION,
  HEADER,
  ENDSEC,
  SEQEND,
  BLOCKS,
  BLOCK,
  ENDBLK,
  VERTEX,
  TABLES,
  TABLE,
  LAYER,
  CONTIN,
  ENDTAB,
  ENTITIES,
  THREE_D_FACE,
  END_OF_FILE,
  POLYLINE,
  INTEGER,
  REAL,
  NAME,
  NTOKENS
} DXFEntity;

extern char *entityStrings[NTOKENS];


typedef enum {
  DXF_START=0,
  DXF_TEXTVAL=1,
  DXF_NAME=2,
  DXF_OTHERNAME1=3,
  DXF_OTHERNAME2=4,
  DXF_ENTITYHANDLE=5,
  DXF_LINETYPE=6,
  DXF_TEXTSTYLE=7,
  DXF_LAYERNAME=8,
  DXF_PRIMARY_X=10,
  DXF_OTHER_X_1=11,
  DXF_PRIMARY_Y=20,
  DXF_OTHER_Y_1=21,
  DXF_PRIMARY_Z=30,
  DXF_OTHER_Z=31,
  DXF_ELEVATION=38,
  DXF_THICKNESS=39,
  DXF_FLOATVALE=40,
  DXF_REPEAT=49,
  DXF_ANGLE1=50,
  DXF_ANGLE2=51,
  DXF_ANGLE3=52,
  DXF_ANGLE4=53,
  DXF_ANGLE5=54,
  DXF_ANGLE6=55,
  DXF_ANGLE7=56,
  DXF_ANGLE8=57,
  DXF_ANGLE9=58,
  DXF_COLORNUM=62,
  DXF_ENTITIES_FLG=66,
  DXF_ENT_IDENT=67,
  DXF_SEVENTYFLAG=70,
  DXF_SEVENTYONEFLAG=71,
  DXF_SEVENTYTWOFLAG=72,
  DXF_VIEWSTAT=69,
  DXF_COMMENT=999
} DXFValues;

#endif
