
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/glut.h>
#else
#  include <GL/glut.h>
#endif


#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"
#include "delfem/drawer_msh.h"

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"

Com::View::CCamera mvp_trans;
double mov_begin_x, mov_begin_y;
int press_button;

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	if( press_button == GLUT_MIDDLE_BUTTON ){
		mvp_trans.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	else if( press_button == GLUT_RIGHT_BUTTON ){
		mvp_trans.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

void myGlutResize(int w, int h)
{
	mvp_trans.SetWindowAspect((double)w/h);

	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
}

Com::View::CDrawerArray drawer_ary;

void myGlutSpecial(int Key, int x, int y)
{
	switch(Key)
	{
	case GLUT_KEY_PAGE_UP:
		if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
			if( mvp_trans.IsPers() ){
				const double tmp_fov_y = mvp_trans.GetFovY() + 10.0;
				mvp_trans.SetFovY( tmp_fov_y );
			}
		}
		else{
			const double tmp_scale = mvp_trans.GetScale() * 0.9;
			mvp_trans.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_PAGE_DOWN:
		if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
			if( mvp_trans.IsPers() ){
				const double tmp_fov_y = mvp_trans.GetFovY() - 10.0;
				mvp_trans.SetFovY( tmp_fov_y );
			}
		}
		else{
			const double tmp_scale = mvp_trans.GetScale() * 1.111;
			mvp_trans.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_HOME :
		drawer_ary.InitTrans(mvp_trans);
		mvp_trans.Fit();
		break;
	case GLUT_KEY_END :
		if( mvp_trans.IsPers() ) mvp_trans.SetIsPers(false);
		else{ mvp_trans.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	::glutPostRedisplay();
}

void myGlutDisplay(void)
{
	::glClearColor(1.0, 1.0, 1.0, 1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1, 4.0 );
   
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(mvp_trans);

	drawer_ary.Draw();
	::glutSwapBuffers();
}






void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/(double)win_w;
	mov_begin_y = (win_h-2.0*y)/(double)win_h;
	press_button = button;
	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ){
		const unsigned int size_buffer = 128;
		unsigned int select_buffer[size_buffer];
		Com::View::PickPre(size_buffer,select_buffer,   x,y,  5,5,   mvp_trans);
		drawer_ary.DrawSelection();
		std::vector<Com::View::SSelectedObject> aSelecObj
			= Com::View::PickPost(select_buffer,   x,y,   mvp_trans );

		if( aSelecObj.size() > 0 ){
			drawer_ary.AddSelected( aSelecObj[0].name );	
		}
		else{
			drawer_ary.ClearSelected();
		}
	}
}

bool SetNewProblem()
{
	const unsigned int nprob = 14;
	static unsigned int iprob = 13;

	if( iprob == 0 )
	{
		Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.2) );
			vec_ary.push_back( Com::CVector2D(0.5,0.2) );
			vec_ary.push_back( Com::CVector2D(0.5,0.8) );
			vec_ary.push_back( Com::CVector2D(1.0,0.8) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.ConnectVertex_Line(6,3);
		}
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.05);
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 1 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 2 )
	{
		Msh::CMesher2D mesh_2d;
		std::cout << "hoge0" << std::endl;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		std::cout << "hoge1" << std::endl;
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		std::cout << "hoge2" << std::endl;
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 3 )
	{
		Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.5) );
			vec_ary.push_back( Com::CVector2D(0.5,0.5) );
			vec_ary.push_back( Com::CVector2D(0.5,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d, 0.1);
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 0.5, 0.1 );
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_3d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_3d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 4 )	// ２次元問題の設定
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 5.0, 0.5 );
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 5 )
	{
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_hex.msh");
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_3d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_3d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 6 )
	{
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_tet.msh");
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_3d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_3d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 7 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 5.0, 0.5 );
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 8 )
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l = 0;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary );
		}
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.6,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.4,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.4) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.6,0.4) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.4,0.4) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.4) );
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.02);
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 9 )
	{
		Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.resize(4);
			vec_ary[0] = Com::CVector2D(0.0,0.0);
			vec_ary[1] = Com::CVector2D(1.0,0.0);
			vec_ary[2] = Com::CVector2D(1.0,1.0);
			vec_ary[3] = Com::CVector2D(0.0,1.0);
			cad_2d.AddPolygon( vec_ary );
		}
		cad_2d.SetCurve_Arc(1,true, -0.5);
		cad_2d.SetCurve_Arc(2,false,-0.5);
		cad_2d.SetCurve_Arc(3,true, -0.5);
		cad_2d.SetCurve_Arc(4,false,-0.5);
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.05);
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み込み
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 10 )
	{
		Cad::CCadObj2D cad_2d;
		{	// 正方形に矩形の穴
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.2));
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.2));
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.8));
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.8));
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.05,1);
		{	// ファイル保存
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// ファイル読み出し
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 11 )	// カットされたメッシュ
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l;
		unsigned int id_e1, id_e2, id_e3, id_e4, id_e5;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(0.3,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary );
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.3,0.5) );
			id_e1 = cad_2d.ConnectVertex_Line(2,id_v1);
			unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.7,0.5) );
			unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.7,0.2) );
			unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.7,0.8) );
			unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.5) );
			unsigned int id_v6 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.9,0.5) );
			id_e2 = cad_2d.ConnectVertex_Line(id_v2,id_v3);
			id_e3 = cad_2d.ConnectVertex_Line(id_v2,id_v4);
			id_e4 = cad_2d.ConnectVertex_Line(id_v2,id_v5);
			id_e5 = cad_2d.ConnectVertex_Line(id_v2,id_v6);
		}
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.2);
//		mesh_2d.Tesselation(cad_2d);

		{
			std::vector<unsigned int> aIdMsh_Inc;
			aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l,Cad::LOOP) );
			std::vector<unsigned int> aIdMshBar_Cut;
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e1,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e2,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e5,Cad::EDGE) );
			////////////////
			std::vector< std::vector<int> > aLnods;
			std::vector<unsigned int> mapVal2Co;
			mesh_2d.GetClipedMesh(aLnods,mapVal2Co,   aIdMsh_Inc,aIdMshBar_Cut);
		}

		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 12 )	// カットされたメッシュ
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l1, id_l2;
		unsigned int id_e3, id_e4;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.5, 0.0) );
			vec_ary.push_back( Com::CVector2D(2.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(2.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(1.5, 1.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );
			unsigned int id_l0 = cad_2d.AddPolygon( vec_ary );
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l0, Com::CVector2D(0.5,0.5) );
			unsigned int id_e1 = cad_2d.ConnectVertex_Line(2,7);
			unsigned int id_e2 = cad_2d.ConnectVertex_Line(3,6);
			unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,id_e1, Com::CVector2D(1.0,0.5) );
			unsigned int id_v3 = cad_2d.AddVertex(Cad::EDGE,1,     Com::CVector2D(0.5,0.0) );
			id_e3 = cad_2d.ConnectVertex_Line(id_v1,id_v2);
			id_e4 = cad_2d.ConnectVertex_Line(id_v1,id_v3);
			id_l1 = 1;
			id_l2 = 2;
		}
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.2);
//		mesh_2d.Tesselation(cad_2d);
		{
			std::vector<unsigned int> aIdMsh_Inc;
			aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l1,Cad::LOOP) );
			aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l2,Cad::LOOP) );
			std::vector<unsigned int> aIdMshBar_Cut;
//			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
			////////////////
			std::vector< std::vector<int> > aLnods;
			std::vector<unsigned int> mapVal2Co;
			mesh_2d.GetClipedMesh(aLnods,mapVal2Co,   aIdMsh_Inc,aIdMshBar_Cut);
		}

		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 13 ){
		Cad::CCadObj2D cad_2d;
		unsigned int id_e3, id_e4;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );	// 1
			vec_ary.push_back( Com::CVector2D(1.5, 0.0) );	// 2
			vec_ary.push_back( Com::CVector2D(1.5, 0.4) );	// 3
			vec_ary.push_back( Com::CVector2D(1.0, 0.4) );	// 4
			vec_ary.push_back( Com::CVector2D(1.0, 0.5) );	// 5
			vec_ary.push_back( Com::CVector2D(2.0, 0.5) );	// 6
			vec_ary.push_back( Com::CVector2D(2.0, 1.0) );	// 7
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );	// 8
			vec_ary.push_back( Com::CVector2D(0.0, 0.5) );	// 9
			unsigned int id_l0 = cad_2d.AddPolygon( vec_ary );
			unsigned int id_e1 = cad_2d.ConnectVertex_Line(5,9);
			cad_2d.ShiftLayer(id_l0,true);
			const double col[3] = { 0.9, 0.4, 0.4 };
			cad_2d.SetColor_Loop(id_l0, col);
			cad_2d.AddVertex(Cad::EDGE,3, Com::CVector2D(1.3,0.5) );
		}
		Msh::CMesher2D mesh_2d;
		mesh_2d.Meshing_ElemLength(cad_2d,0.05);
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( mvp_trans );
	}

	iprob++;
	if( iprob == nprob ) iprob = 0;

	return true;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' は ESC の ASCII コード */
	  exit(0);
	  break;
  case ' ':
	  SetNewProblem();
	  ::glMatrixMode(GL_PROJECTION);
	  ::glLoadIdentity();
	  Com::View::SetProjectionTransform(mvp_trans);
	  break;
  default:
    break;
  }
}

int main(int argc,char* argv[])
{
	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 400);
	::glutInit(&argc, argv);
	::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("DelFEM demo");

    // Set call back function
	::glutMotionFunc(myGlutMotion);
	::glutMouseFunc(myGlutMouse);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutSpecialFunc(myGlutSpecial);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutIdleFunc(myGlutIdle);
	
	SetNewProblem();

    // Enter Main Loop
	::glutMainLoop();
	return 0;
}
