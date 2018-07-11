////////////////////////////////////////////////////////////
//
// Name:    dxf.cpp
//
// Author:  Steven Michael
//          (smichael@ll.mit.edu)
//
// Date:    3/10/2005
//
// Description:
//
//   The following is a c++ class implentation of a 
//   DXF file reader
//
//   It does not handle all the DXF format specifications,
//   but is sufficient to read in all the ones that I am
//   interested in.  Note that it does not handle colors.
//
////////////////////////////////////////////////////////////

#include "dxf.h"

#include <string>
#include <stdio.h>
#include <fstream>

DXF::DXF(const char *fname)
{
  tokens = (Token *)NULL;
  nTok = 0;
  facets = (Facet *)NULL;
  if(fname) {
    if(read_file(fname)) return;
    if(load_facets()) return;
  }
  vertices = (Vertex *)NULL;
  faceIndices = (FaceIndex *)NULL;
  nVert = 0;
  nFaceIndex = 0;

} // end of constructor

DXF::~DXF()
{
  free_tokens();
  if(facets) delete[] facets;
  if(vertices) free(vertices);
  if(faceIndices) free(faceIndices);
  return;
} // end of destructor

int DXF::free_tokens()
{
  if(!tokens) return 0;
  delete[] tokens;
  tokens = (Token *)NULL;
  nTok = 0;
  return 0;
} // end of free_tokens

int DXF::load_facets()
{
  // Delete the old facets,
  // if they exist
  if(facets) delete[] facets;

  // First, count max the number of facets to allocate
  // Each 3-D Face could have 2 facets
  int nalloc=0;
  for(int i=0;i<nTok;i++) {
    if(tokens[i].type==0 && tokens[i].val.e == THREE_D_FACE)
      nalloc++;
    else if(tokens[i].type==0 && tokens[i].val.e == VERTEX) {
      i++;
      while(tokens[i].type!=0) {
	if(tokens[i].type == 70 && tokens[i].val.i == 128)
	  nalloc++;
	i++;
      }
      i--;
    }
  }
  nalloc*=2;
  facets = new Facet[nalloc];
  
  // Read in the facets
  nFacet = 0;
  for(int i=0;i<nTok;i++) {
    int nFacRead;
    // Look for a 3D face
    if((tokens[i].type==0) && (tokens[i].val.e == THREE_D_FACE)) {
      int ret = read_3d_face(tokens+i,facets+nFacet,nFacRead);
      if(ret < 0) {
	printf("Error reading 3d face.\n");
	return -1;
      }
      i += ret;
      nFacet += nFacRead;
    }
    // Look for a vertex
    else if((tokens[i].type==0) && (tokens[i].val.e == VERTEX)) {
      int ret = read_vertex(tokens+i);
      if(ret < 0) {
	printf("Error reading vertex.\n");
	return -1;
      }
      i += ret;
    }
    // Look for a polyline
    else if((tokens[i].type == 0) && (tokens[i].val.e == POLYLINE)) {
      indices_to_faces();
      reset_polyline();
    }
  }		   
  // Convert vertex indices to faces
  // If necessary
  indices_to_faces();

  return 0;
} // end of load_facets
  
int DXF::indices_to_faces(void)
{
  Facet *f  = facets+nFacet;

  for(int i=0;i<nFaceIndex;i++) {
    if(faceIndices[i][0] > nVert || faceIndices[i][1] > nVert || 
       faceIndices[i][2] > nVert || faceIndices[i][3] > nVert) {
      printf("Invalid face index: %d\n",i);
      return -1;
    }
    memcpy((*f)[0],vertices[faceIndices[i][0]-1],sizeof(Vertex));
    memcpy((*f)[1],vertices[faceIndices[i][1]-1],sizeof(Vertex));
    memcpy((*f)[2],vertices[faceIndices[i][2]-1],sizeof(Vertex));
    nFacet++;
    f++;
    if(faceIndices[i][3] != -1) {
      memcpy((*f)[0],vertices[faceIndices[i][0]-1],sizeof(Vertex));
      memcpy((*f)[1],vertices[faceIndices[i][2]-1],sizeof(Vertex));
      memcpy((*f)[2],vertices[faceIndices[i][3]-1],sizeof(Vertex));
      nFacet++;
      f++;
    } 
  }
  return 0;
}


int DXF::reset_polyline(void)
{
  if(vertices) {
    free(vertices);
    vertices = (Vertex *)NULL;
  }
  if(faceIndices) {
    free(faceIndices);
    faceIndices = (FaceIndex *)NULL;
  }
  nVert = 0;
  nFaceIndex = 0;
	return 0;
} // end of reset_polyline


#define NALLOC 1000
int DXF::add_vertex(Vertex v)
{
  // Allocate memory
  if(!vertices) {
    vertices = (Vertex *)malloc(sizeof(Vertex)*NALLOC);
    nVert = 0;
  }
  else if(nVert % NALLOC == 0) {
    vertices = (Vertex *)realloc(vertices,sizeof(Vertex)*(nVert+NALLOC));
  }
  memcpy(vertices[nVert],v,sizeof(Vertex));
  nVert++;

  return 0;
} 

int DXF::add_face_index(FaceIndex f)
{
  if(!faceIndices) {
    faceIndices = (FaceIndex *)malloc(sizeof(FaceIndex)*NALLOC);
    nFaceIndex = 0;
  }
  else if(nFaceIndex % NALLOC == 0) {
    faceIndices = (FaceIndex *)realloc(faceIndices,
				       sizeof(FaceIndex)*(nFaceIndex+NALLOC));
  }
  memcpy(faceIndices[nFaceIndex],f,sizeof(FaceIndex));
  nFaceIndex++;
  return 0;
} // end of add_face_index

int DXF::read_vertex(Token *t)
{
  Vertex v1,v2,v3,v4;
  FaceIndex f;
  int type = 0;
  int count=0;
  v1[0] = 0.0;v1[1] = 0.0;v1[2] = 0.0;
  f[0] = -1;f[1] = -1;f[2] = -1;f[3] = -1;

  t++;
  while(t->type != 0) {
    switch(t->type) {
    case 10:
      v1[0] = (float)t->val.d;
      break;
    case 20:
      v1[1] = (float)t->val.d;
      break;
    case 30:
      v1[2] = (float)t->val.d;
      break;
    case 11:
      v2[0] = (float)t->val.d;
      break;
    case 21:
      v2[1] = (float)t->val.d;
      break;
    case 31:
      v2[2] = (float)t->val.d;
      break;
    case 12:
      v3[0] = (float)t->val.d;
      break;
    case 22:
      v3[1] = (float)t->val.d;
      break;
    case 32:
      v3[2] = (float)t->val.d;
      break;
    case 13:
      v4[0] = (float)t->val.d;
      break;
    case 23:
      v4[1] = (float)t->val.d;
      break;
    case 33:
      v4[2] = (float)t->val.d;
      break;
    case 70:
      type = t->val.i;
      break;
    case 71:
      f[0] = t->val.i;
      break;
    case 72:
      f[1] = t->val.i;
      break;
    case 73:
      f[2] = t->val.i;
      break;
    case 74:
      f[3] = t->val.i;
      break;
    default: 
      break;
    }
    t++;
    count++;
  }
  if(type==128)
    add_face_index(f);
  else if(type==192)
    add_vertex(v1);
  
  return count;
} // end of read_vertex
      

int DXF::read_3d_face(Token *t, Facet *f, int &nFac)
{

  nFac=0;
  int count=0;
  if(t->type != 0 && t->val.e != THREE_D_FACE) 
    return -1;
  
  t++;

  // Read in 4 possible vertices
  // (That could define 2 facets)
  Vertex v1,v2,v3,v4;
  while(t->type != 0) {
    switch(t->type) {
    case 10:
      v1[0] = (float)t->val.d;
      break;
    case 20:
      v1[1] = (float)t->val.d;
      break;
    case 30:
      v1[2] = (float)t->val.d;
      break;
    case 11:
      v2[0] = (float)t->val.d;
      break;
    case 21:
      v2[1] = (float)t->val.d;
      break;
    case 31:
      v2[2] = (float)t->val.d;
      break;
    case 12:
      v3[0] = (float)t->val.d;
      break;
    case 22:
      v3[1] = (float)t->val.d;
      break;
    case 32:
      v3[2] = (float)t->val.d;
      break;
    case 13:
      v4[0] = (float)t->val.d;
      break;
    case 23:
      v4[1] = (float)t->val.d;
      break;
    case 33:
      v4[2] = (float)t->val.d;
      break;
    default: break;
    }
    t++;
    count++;
  }

  nFac=1;
  // Copy the vertices into the facet
  memcpy((*f)[0],v1,sizeof(Vertex));
  memcpy((*f)[1],v2,sizeof(Vertex));
  memcpy((*f)[2],v3,sizeof(Vertex));


  // See if there is a 4th vertex
  if(v3[0] != v4[0] || 
     v3[1] != v4[1] || 
     v3[2] != v4[2]) {
    nFac=2;
    f++;
    memcpy((*f)[0],v1,sizeof(Vertex));
    memcpy((*f)[1],v3,sizeof(Vertex));
    memcpy((*f)[2],v4,sizeof(Vertex));
  }

  return count;
} // end of read_3d_face

int DXF::read_file(const char *filename)
{
  // Read in the file to a buffer
	std::ifstream is(filename);
  if(is.is_open() == 0) {
    printf("file not open\n");
    return -1;
  }
 

	is.seekg(0,std::ios::end);
  int fsize = is.tellg();
	is.seekg(0,std::ios::beg);
  
  char *buf = new char[fsize+1];
  is.read(buf,fsize);
  buf[fsize] = '\0';
  is.close();

  int nlines=0;
  char *s,*s2;
  s = buf;s2=s;
  while(s <= buf+fsize) {
    s2 = strstr(s,"\n");
    if(s2==0)
      break;
    nlines++;
    s = s2+1;
  }
  tokens = new Token[nlines/2];
  memset(tokens,0,sizeof(Token)*nlines/2);

  s=buf;s2=s;
  
  nTok=-1;
	while(s<=buf+fsize-1) {

    // Increment the # of tokens
    nTok++;

    // Truncate the line
    s2 = strstr(s,"\n");
    if(s2) *s2='\0';
		else
			break;

    // Read in the token type
    tokens[nTok].type = atoi(s);
    //sscanf(s,"%d",&tokens[nTok].type);

    // Goto the next line
    s = s2+1;
    if(s >= buf+fsize) break;

    // Truncate the line
    s2 = strstr(s,"\n");
    if(s2) *s2 = '\0';

    // Read in the token values
    if(tokens[nTok].type == 0) {
      tokens[nTok].val.e = BAD_TOKEN;
      for(int i=0;i<NTOKENS;i++) {
	if(!strncmp(s,entityStrings[i],strlen(entityStrings[i]))) {
	  tokens[nTok].val.e = (DXFEntity) i;
	  break;
	}
      }
    }
    else if(tokens[nTok].type < 60)
      tokens[nTok].val.d = atof(s);
    else if(tokens[nTok].type < 100) 
      tokens[nTok].val.i = atoi(s);
    else if(tokens[nTok].type < 150)
      tokens[nTok].val.d = atof(s);
    else if(tokens[nTok].type < 180)
      tokens[nTok].val.i = atoi(s);

    // Goto the next line
    s = s2+1;
    if(s >= buf+fsize) break;
  }
	delete[] buf;

  return 0;
} // end of read_file

#ifdef _TEST_
int main(int argc, char *argv[])
{
  DXF *dxf = new DXF;
  dxf->read_file(argv[1]);
  delete dxf;
  return 0;
}
#endif
