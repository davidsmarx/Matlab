/* this file contains the actual definitions of */
/* the IIDs and CLSIDs */

/* link this file in with the server and any clients */


/* File created by MIDL compiler version 5.01.0164 */
/* at Tue Aug 09 12:28:40 2005
 */
/* Compiler settings for C:\My_Documents\internal_metrology\diffraction\int_met_comserver\server.idl:
    Os (OptLev=s), W1, Zp8, env=Win32, ms_ext, c_ext
    error checks: allocation ref bounds_check enum stub_data 
*/
//@@MIDL_FILE_HEADING(  )
#ifdef __cplusplus
extern "C"{
#endif 


#ifndef __IID_DEFINED__
#define __IID_DEFINED__

typedef struct _IID
{
    unsigned long x;
    unsigned short s1;
    unsigned short s2;
    unsigned char  c[8];
} IID;

#endif // __IID_DEFINED__

#ifndef CLSID_DEFINED
#define CLSID_DEFINED
typedef IID CLSID;
#endif // CLSID_DEFINED

const IID IID_IDiffraction = {0x44182096,0xCE96,0x4DC8,{0x8D,0x63,0xEC,0x2F,0xC6,0xFA,0x8D,0x85}};


const IID LIBID_DiffractionLib = {0x7FB676A4,0x8C3D,0x49A4,{0x80,0x4D,0x07,0x35,0x30,0xB6,0xB4,0x00}};


const CLSID CLSID_Diffraction = {0xE4BF1897,0x350C,0x488F,{0xB2,0x46,0xE1,0x4B,0xAF,0x2A,0x07,0x56}};


#ifdef __cplusplus
}
#endif

