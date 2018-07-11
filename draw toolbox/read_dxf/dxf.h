////////////////////////////////////////////////////////////
//
// Name:    dxf.header
//
// Author:  Steven Michael
//          (smichael@ll.mit.edu)
//
// Date:    3/10/2005
//
// Description:
//
//   The following is a c++ class header of a 
//   DXF file reader
//
//   It does not handle all the DXF format specifications,
//   but is sufficient to read in all the ones that I am
//   interested in.  Note that it does not handle colors.
//
////////////////////////////////////////////////////////////

#ifndef _DXF_H_
#define _DXF_H_

#include "dxftoken.h"



typedef union {
  double     d;
  int        i;
  char      *s;
  DXFEntity  e;
} TokVal;


typedef struct {
  int    type;
  TokVal val;
} Token;


typedef float  Vertex[3];
typedef Vertex Facet[3];
typedef int    FaceIndex[4];


class DXF {
 public:
  DXF(const char *fname=(const char *)0);
  virtual ~DXF();

  int read_file(const char *);
  int load_facets();

  int    nFacet;
  Facet *facets;

 protected:
  int    nTok;
  Token *tokens;
  
  int free_tokens();
  int read_3d_face(Token *, Facet *, int &);
  int read_vertex(Token *);

  Vertex       *vertices;
  FaceIndex    *faceIndices;
  int           nVert;
  int           nFaceIndex;

  int add_vertex(Vertex);
  int add_face_index(FaceIndex);

  int reset_polyline();
  int indices_to_faces();

}; // end of DXF class definition

#endif
