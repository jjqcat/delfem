# Microsoft Developer Studio Project File - Name="stlib_cadmshview" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** 編集しないでください **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=stlib_cadmshview - Win32 Debug
!MESSAGE これは有効なﾒｲｸﾌｧｲﾙではありません。 このﾌﾟﾛｼﾞｪｸﾄをﾋﾞﾙﾄﾞするためには NMAKE を使用してください。
!MESSAGE [ﾒｲｸﾌｧｲﾙのｴｸｽﾎﾟｰﾄ] ｺﾏﾝﾄﾞを使用して実行してください
!MESSAGE 
!MESSAGE NMAKE /f "stlib_cadmshview.mak".
!MESSAGE 
!MESSAGE NMAKE の実行時に構成を指定できます
!MESSAGE ｺﾏﾝﾄﾞ ﾗｲﾝ上でﾏｸﾛの設定を定義します。例:
!MESSAGE 
!MESSAGE NMAKE /f "stlib_cadmshview.mak" CFG="stlib_cadmshview - Win32 Debug"
!MESSAGE 
!MESSAGE 選択可能なﾋﾞﾙﾄﾞ ﾓｰﾄﾞ:
!MESSAGE 
!MESSAGE "stlib_cadmshview - Win32 Release" ("Win32 (x86) Static Library" 用)
!MESSAGE "stlib_cadmshview - Win32 Debug" ("Win32 (x86) Static Library" 用)
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "stlib_cadmshview - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /D "_MBCS" /D "_LIB" /D "WIN32" /D "NDEBUG" /D "__VISUALC__" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x411 /d "NDEBUG"
# ADD RSC /l 0x411 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "stlib_cadmshview - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "_LIB" /D "__VISUALC__" /D "WIN32" /D "_DEBUG" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x411 /d "_DEBUG"
# ADD RSC /l 0x411 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "stlib_cadmshview - Win32 Release"
# Name "stlib_cadmshview - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\cad\brep.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cad\brep2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cad\cad_elem2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cad\cad_obj2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cad\drawer_cad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\drawer_msh.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\mesh3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\mesh3d_extrude.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\mesher2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\mesher3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\meshkernel2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\msh\meshkernel3d.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\include\delfem\cad\brep.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\cad\brep2d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\cad2d_interface.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\cad_com.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\cad\cad_elem2d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\cad_obj2d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\DelFEM\drawer_cad.h
# End Source File
# Begin Source File

SOURCE=..\..\include\DelFEM\drawer_msh.h
# End Source File
# Begin Source File

SOURCE=..\..\include\DelFEM\mesh3d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\DelFEM\mesh_interface.h
# End Source File
# Begin Source File

SOURCE=..\..\include\DelFEM\mesh_primitive.h
# End Source File
# Begin Source File

SOURCE=..\..\include\DelFEM\mesher2d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\msh\meshkernel2d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\msh\meshkernel3d.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\cad\objset_cad.h
# End Source File
# End Group
# End Target
# End Project
