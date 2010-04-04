# Microsoft Developer Studio Project File - Name="stlib_femmat" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** 編集しないでください **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=stlib_femmat - Win32 Debug
!MESSAGE これは有効なﾒｲｸﾌｧｲﾙではありません。 このﾌﾟﾛｼﾞｪｸﾄをﾋﾞﾙﾄﾞするためには NMAKE を使用してください。
!MESSAGE [ﾒｲｸﾌｧｲﾙのｴｸｽﾎﾟｰﾄ] ｺﾏﾝﾄﾞを使用して実行してください
!MESSAGE 
!MESSAGE NMAKE /f "stlib_femmat.mak".
!MESSAGE 
!MESSAGE NMAKE の実行時に構成を指定できます
!MESSAGE ｺﾏﾝﾄﾞ ﾗｲﾝ上でﾏｸﾛの設定を定義します。例:
!MESSAGE 
!MESSAGE NMAKE /f "stlib_femmat.mak" CFG="stlib_femmat - Win32 Debug"
!MESSAGE 
!MESSAGE 選択可能なﾋﾞﾙﾄﾞ ﾓｰﾄﾞ:
!MESSAGE 
!MESSAGE "stlib_femmat - Win32 Release" ("Win32 (x86) Static Library" 用)
!MESSAGE "stlib_femmat - Win32 Debug" ("Win32 (x86) Static Library" 用)
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "stlib_femmat - Win32 Release"

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
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /D "_LIB" /D "__VISUALC__" /D "_MBCS" /D "WIN32" /D "NDEBUG" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x411 /d "NDEBUG"
# ADD RSC /l 0x411 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "stlib_femmat - Win32 Debug"

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

# Name "stlib_femmat - Win32 Release"
# Name "stlib_femmat - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\ls\eigen_lanczos.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ls\linearsystem.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femls\linearsystem_field.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\femls\linearsystem_fieldsave.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\mat_blkcrs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\matdia_blkcrs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\matdiafrac_blkcrs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\matdiainv_blkdia.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\matfrac_blkcrs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\matprolong_blkcrs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\ordering_blk.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ls\preconditioner.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ls\solver_ls_iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\solver_mat_iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\solver_mg.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matvec\vector_blk.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\include\delfem\matvec\bcflag_blk.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\diamat_blk.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\ls\eigen_lanczos.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\ker_mat.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\ls\linearsystem.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\femls\linearsystem_field.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\femls\linearsystem_fieldsave.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\linearsystem_interface_eqnsys.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\ls\linearsystem_interface_solver.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\mat_blkcrs.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\matdia_blkcrs.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\matdiafrac_blkcrs.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\matdiainv_blkdia.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\matfrac_blkcrs.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\matprecond_blk.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\matprolong_blkcrs.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\ordering_blk.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\ls\preconditioner.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\ls\solver_interface.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\ls\solver_ls_iter.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\solver_mat_iter.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\solver_mg.h
# End Source File
# Begin Source File

SOURCE=..\..\include\delfem\matvec\vector_blk.h
# End Source File
# End Group
# End Target
# End Project
