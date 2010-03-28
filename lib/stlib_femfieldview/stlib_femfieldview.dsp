# Microsoft Developer Studio Project File - Name="stlib_femfieldview" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** 編集しないでください **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=stlib_femfieldview - Win32 Debug
!MESSAGE これは有効なﾒｲｸﾌｧｲﾙではありません。 このﾌﾟﾛｼﾞｪｸﾄをﾋﾞﾙﾄﾞするためには NMAKE を使用してください。
!MESSAGE [ﾒｲｸﾌｧｲﾙのｴｸｽﾎﾟｰﾄ] ｺﾏﾝﾄﾞを使用して実行してください
!MESSAGE 
!MESSAGE NMAKE /f "stlib_femfieldview.mak".
!MESSAGE 
!MESSAGE NMAKE の実行時に構成を指定できます
!MESSAGE ｺﾏﾝﾄﾞ ﾗｲﾝ上でﾏｸﾛの設定を定義します。例:
!MESSAGE 
!MESSAGE NMAKE /f "stlib_femfieldview.mak" CFG="stlib_femfieldview - Win32 Debug"
!MESSAGE 
!MESSAGE 選択可能なﾋﾞﾙﾄﾞ ﾓｰﾄﾞ:
!MESSAGE 
!MESSAGE "stlib_femfieldview - Win32 Release" ("Win32 (x86) Static Library" 用)
!MESSAGE "stlib_femfieldview - Win32 Debug" ("Win32 (x86) Static Library" 用)
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "stlib_femfieldview - Win32 Release"

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

!ELSEIF  "$(CFG)" == "stlib_femfieldview - Win32 Debug"

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
# ADD CPP /nologo /MDd /W2 /Gm /GX /ZI /Od /D "_LIB" /D "__VISUALC__" /D "WIN32" /D "_DEBUG" /D "_MBCS" /YX /FD /GZ /c
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

# Name "stlib_femfieldview - Win32 Release"
# Name "stlib_femfieldview - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\femfield\drawer_field.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\drawer_field_edge.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\drawer_field_face.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\drawer_field_image_based_flow_vis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\drawer_field_streamline.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\drawer_field_vector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\elem_ary.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\eval.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\field.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\field_world.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femeqn\ker_emat_tri.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femfield\node_ary.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\include\delfem\drawer_field.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\drawer_field_edge.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\drawer_field_face.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\drawer_field_image_based_flow_vis.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\drawer_field_streamline.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\drawer_field_vector.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\elem_ary.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\eval.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\field.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\field_world.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\node_ary.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\objset.h
# End Source File
# End Group
# End Target
# End Project
