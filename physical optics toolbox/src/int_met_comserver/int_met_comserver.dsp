# Microsoft Developer Studio Project File - Name="int_met_comserver" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=int_met_comserver - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "int_met_comserver.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "int_met_comserver.mak" CFG="int_met_comserver - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "int_met_comserver - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "int_met_comserver - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "int_met_comserver - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt /winapp
# ADD F90 /automatic /browser /compile_only /libs:dll /nologo /threads /warn:argument_checking /warn:nofileopt /winapp
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /win32 /cpp_cmdfpp /cpp_opt"/a /m /B /extend_source 132"
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /i "$(INTDIR)" /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 kernel32.lib /nologo /stack:0x8000000,0x1000000 /subsystem:windows /machine:I386
# Begin Custom Build - Performing Registration
OutDir=.\Release
TargetPath=.\Release\int_met_comserver.exe
InputPath=.\Release\int_met_comserver.exe
SOURCE="$(InputPath)"

"$(OutDir)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	"$(TargetPath)" /Regserver > "$(OutDir)\regsvr32.trg"

# End Custom Build

!ELSEIF  "$(CFG)" == "int_met_comserver - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt /winapp
# ADD F90 /automatic /browser /check:bounds /compile_only /dbglibs /debug:full /libs:dll /nologo /reentrancy:threaded /threads /traceback /warn:argument_checking /warn:nofileopt /winapp
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /win32 /cpp_cmdfpp /cpp_opt"/a /m /B /extend_source 132"
# SUBTRACT MTL /mktyplib203
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /i "$(INTDIR)" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /stack:0x8000000,0x1000000 /subsystem:windows /incremental:no /debug /machine:I386 /pdbtype:sept
# Begin Custom Build - Performing Registration
OutDir=.\Debug
TargetPath=.\Debug\int_met_comserver.exe
InputPath=.\Debug\int_met_comserver.exe
SOURCE="$(InputPath)"

"$(OutDir)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	"$(TargetPath)" /Regserver > "$(OutDir)\regsvr32.trg"

# End Custom Build

!ENDIF 

# Begin Target

# Name "int_met_comserver - Win32 Release"
# Name "int_met_comserver - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE="..\2003.10.14 update\Angle_Units.F90"
DEP_F90_ANGLE=\
	".\Release\Constants.mod"\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\Array_Routines_fftw.F90"
DEP_F90_ARRAY=\
	".\Release\Constants.mod"\
	".\Release\Discrete_Transforms.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Clipping_Routines.f90"
DEP_F90_CLIPP=\
	".\Release\Array_Routines.mod"\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	".\Release\Polygons.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Constants.F90"
DEP_F90_CONST=\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\Discrete_Transforms_fftw.F90"
DEP_F90_DISCR=\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\fftw_module.mod"\
	".\Release\Kinds.mod"\
	".\Release\Math_Routines.mod"\
	".\Release\Singleton_FFT.mod"\
	".\Release\Sorensen_FFT.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Error_Exit.F90"
# End Source File
# Begin Source File

SOURCE=.\fftw_module.f90
DEP_F90_FFTW_=\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\int_met_routines.f90
DEP_F90_INT_M=\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	".\Release\Optics_Routines.mod"\
	".\Release\Polygons.mod"\
	".\Release\SI_Units.mod"\
	".\Release\Utility_Routines.mod"\
	".\Release\Wavefronts.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Kinds.F90"
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Math_Routines.F90"
DEP_F90_MATH_=\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\MyCorner_Cubes.F90"
DEP_F90_MYCOR=\
	".\Release\Array_Routines.mod"\
	".\Release\Clipping_Routines.mod"\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	".\Release\Polygons.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\MyWavefronts_fftw_mt.F90"
DEP_F90_MYWAV=\
	".\Release\Array_Routines.mod"\
	".\Release\Clipping_Routines.mod"\
	".\Release\Constants.mod"\
	".\Release\Corner_Cubes.mod"\
	".\Release\Discrete_Transforms.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\fftw_module.mod"\
	".\Release\Kinds.mod"\
	".\Release\Math_Routines.mod"\
	".\Release\Polygons.mod"\
	".\Release\Quadrature_Routines.mod"\
	".\Release\Utility_Routines.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.11.19 update\Optics_Routines 031202.f90"
DEP_F90_OPTIC=\
	".\Release\Angle_Units.mod"\
	".\Release\Array_Routines.mod"\
	".\Release\Clipping_Routines.mod"\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	".\Release\Polygons.mod"\
	".\Release\SI_Units.mod"\
	".\Release\Wavefronts.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Polygons.F90"
DEP_F90_POLYG=\
	".\Release\Array_Routines.mod"\
	".\Release\Constants.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\Quadrature_Routines_#3F5803.f90"
DEP_F90_QUADR=\
	".\Release\Array_Routines.mod"\
	".\Release\Clipping_Routines.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	".\Release\Math_Routines.mod"\
	".\Release\Polygons.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\SI_Units.F90"
DEP_F90_SI_UN=\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\Singleton_FFT.F90"
DEP_F90_SINGL=\
	".\Release\Kinds.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update fftw\Sorensen_FFT.F90"
DEP_F90_SOREN=\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	".\Release\Math_Routines.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\threadwrappers.f90
DEP_F90_THREA=\
	".\Release\Corner_Cubes.mod"\
	".\Release\IntMetRoutines.mod"\
	".\Release\Kinds.mod"\
	".\Release\Optics_Routines.mod"\
	".\Release\Polygons.mod"\
	".\Release\Wavefronts.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\UDiffractionTY.f90
DEP_F90_UDIFF=\
	".\Release\Kinds.mod"\
	".\Release\SI_Units.mod"\
	".\Release\ThreadWrappers.mod"\
	".\Release\Wavefronts.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\UIDiffraction.f90
DEP_F90_UIDIF=\
	".\Release\Array_Routines.mod"\
	".\Release\Corner_Cubes.mod"\
	".\Release\Diffraction_Types.mod"\
	".\Release\Error_Exit.mod"\
	".\Release\IntMetRoutines.mod"\
	".\Release\Kinds.mod"\
	".\Release\Optics_Routines.mod"\
	".\Release\Polygons.mod"\
	".\Release\ThreadWrappers.mod"\
	".\Release\Wavefronts.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\2003.10.14 update\Utility_Routines.F90"
DEP_F90_UTILI=\
	".\Release\Error_Exit.mod"\
	".\Release\Kinds.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\resource.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\int_met_comserver.rc
# End Source File
# End Group
# Begin Group "Do Not Edit"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\clsfact.f90
DEP_F90_CLSFA=\
	".\Release\Diffraction_Types.mod"\
	".\Release\IClassFactory_Types.mod"\
	".\Release\IDiffraction_Methods.mod"\
	".\Release\int_met_comserver_global.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\clsfactty.f90
# End Source File
# Begin Source File

SOURCE=.\DiffractionTY.f90
DEP_F90_DIFFR=\
	".\Release\Diffraction_USE.mod"\
	".\Release\int_met_comserver_global.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\exemain.f90
DEP_F90_EXEMA=\
	".\Release\Diffraction_Types.mod"\
	".\Release\IClassFactory_Methods.mod"\
	".\Release\IClassFactory_Types.mod"\
	".\Release\int_met_comserver_global.mod"\
	".\Release\ServerHelper.mod"\
	{$(INCLUDE)}"ole32.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\IDiffraction.f90
DEP_F90_IDIFF=\
	".\Release\Diffraction_Types.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\int_met_comserverGlobal.f90
DEP_F90_INT_ME=\
	".\Release\IClassFactory_Types.mod"\
	{$(INCLUDE)}"ole32.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\server.idl
# End Source File
# Begin Source File

SOURCE=.\serverhelper.f90
DEP_F90_SERVE=\
	".\Release\Diffraction_Types.mod"\
	".\Release\int_met_comserver_global.mod"\
	
# End Source File
# End Group
# Begin Source File

SOURCE="..\fftw-3.0.1-w32-pl1\fftw3fortran.lib"
# End Source File
# Begin Source File

SOURCE="..\fftw-3.0.1-w32-pl1\fftw3.lib"
# End Source File
# End Target
# End Project
