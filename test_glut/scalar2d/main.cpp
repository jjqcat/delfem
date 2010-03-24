
#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"		// カメラクラスCCamera
#include "delfem/drawer_gl_utility.h"	// GLのための便利関数群

#include "delfem/cad_obj2d.h"		// ２次元形状クラスCCadObj2D
#include "delfem/mesher2d.h"		// ２次元メッシュクラスCMesher2D

#include "delfem/field.h"	// 有限要素法離散場クラスCField
#include "delfem/field_world.h"		// 有限要素法離散場管理クラスCFieldWorld
#include "delfem/drawer_field.h"	// 有限要素法離散場可視化クラス
#include "delfem/drawer_field_face.h"	// 有限要素法離散場可視化クラス
#include "delfem/drawer_field_edge.h"	// 有限要素法離散場可視化クラス
#include "delfem/drawer_field_vector.h"	// 有限要素法離散場可視化クラス

#include "delfem/eqnsys_scalar.h"	// 有限要素法スカラー型方程式クラス

using namespace Fem::Field;
using namespace Fem::Ls;

Fem::Field::CFieldWorld world;
Fem::Eqn::CEqnSystem_Scalar2D eqn_scalar;
double dt = 0.001;
View::CDrawerArrayField drawer_ary;
Com::View::CCamera mvp_trans;
double mov_begin_x, mov_begin_y;
unsigned int id_base;

void RenderBitmapString(float x, float y, void *font,char *string)
{
  char *c;
  ::glRasterPos2f(x, y);
  for (c=string; *c != '\0'; c++) {
	  ::glutBitmapCharacter(font, *c);
  }
}


void ShowFPS()
{
	int* font=(int*)GLUT_BITMAP_8_BY_13;
	static char s_fps[32];
	{
		static int frame, timebase;
		int time;
		frame++;
		time=glutGet(GLUT_ELAPSED_TIME);
		if (time - timebase > 500) {
			sprintf(s_fps,"FPS:%4.2f",frame*1000.0/(time-timebase));
			timebase = time;
			frame = 0;
		}
	}
	char s_tmp[32];

	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];

	::glMatrixMode(GL_PROJECTION);
	::glPushMatrix();
	::glLoadIdentity();
	::gluOrtho2D(0, win_w, 0, win_h);
	::glMatrixMode(GL_MODELVIEW);
	::glPushMatrix();
	::glLoadIdentity();
	::glScalef(1, -1, 1);
	::glTranslatef(0, -win_h, 0);
//	::glDisable(GL_LIGHTING);
	::glDisable(GL_DEPTH_TEST);
	::glColor3d(1.0, 0.0, 0.0);
	strcpy(s_tmp,"DelFEM demo");
	RenderBitmapString(10,15, (void*)font, s_tmp);
	::glColor3d(0.0, 0.0, 1.0);
	strcpy(s_tmp,"Press \"space\" key!");
	RenderBitmapString(120,15, (void*)font, s_tmp);
	::glColor3d(0.0, 0.0, 0.0);
	RenderBitmapString(10,30, (void*)font, s_fps);
//	::glEnable(GL_LIGHTING);
	::glEnable(GL_DEPTH_TEST);
	::glPopMatrix();
	::glMatrixMode(GL_PROJECTION);
	::glPopMatrix();
	::glMatrixMode(GL_MODELVIEW);
}

bool SetNewProblem()
{
	const unsigned int nprob = 20;
	static unsigned int iprob = 0;

	static int id_val_bc0=0, id_val_bc1=0, id_val_bc2=0;
	
	if( iprob == 0 )	// ２次元問題の設定
	{
		world.Clear();
		drawer_ary.Clear();
		////////////////
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.5));
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.9));
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.9));
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.5));
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}
		// メッシュを作る
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		Fem::Field::CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		eqn_scalar.SetDomain_Field(id_base,world);    // 方程式の設定
		dt = 0.02;
		eqn_scalar.SetTimeIntegrationParameter(dt);
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetAxialSymmetry(false);
		// 全体の方程式の係数設定
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetCapacity(30.0);
		eqn_scalar.SetAdvection(0);

		id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::LOOP),world);
		{
			CField& field = world.GetField(id_val_bc0);
//			field.SetValue("cos(2*PI*t+0.1)",0,world,true);
			field.SetValue("floor(1+0.8*cos(2*PI*t+0.1))",0,Fem::Field::VALUE,world,true);
		}
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(1,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc1);
			field.SetValue(1.0,0,Fem::Field::VALUE,world,false);
		}
		id_val_bc2 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc2);
			field.SetValue(-1.0,0,Fem::Field::VALUE,world,false);
		}

		// 描画オブジェクトの登録
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
//		drawer_ary.PushBack( new View::CDrawerFaceContour(id_field_val,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(mvp_trans);	// 視線座標変換行列の初期化
	}
	else if( iprob == 1 ){
		eqn_scalar.SetCapacity(10);
	}
	else if( iprob == 2 ){
		eqn_scalar.SetCapacity(5);
	}
	else if( iprob == 3 ){
		eqn_scalar.SetCapacity(1);
	}
	else if( iprob == 4 ){
		eqn_scalar.SetStationary(true);
	}
	else if( iprob == 5 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 6 ){
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetCapacity(10);
		dt = 0.02;
		eqn_scalar.SetTimeIntegrationParameter(dt);
	}
	else if( iprob == 7 ){
		const unsigned int id_field_velo = world.MakeField_FieldElemDim(id_base,2,Fem::Field::VECTOR2,VELOCITY);
		std::cout << "Velo : " << id_field_velo << std::endl;
		{	// 流速場の設定
			CField& field = world.GetField(id_field_velo);
			field.SetValue(" (y-0.5)", 0,Fem::Field::VELOCITY,world,false);
			field.SetValue("-(x-0.5)" ,1,Fem::Field::VELOCITY,world,false);
		}
		{	// 固定境界条件の設定
			CField& field = world.GetField(id_val_bc0);
			field.SetValue("floor(1+0.8*cos(3*t))",0,Fem::Field::VALUE,world,true);
		}
		{	// 周囲の固定境界条件の設定
			const CIDConvEAMshCad conv = world.GetIDConverter(id_base);
			std::vector<unsigned int> m_aIDEA;
			m_aIDEA.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(3,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_val_bc1 = eqn_scalar.AddFixElemAry(m_aIDEA,world);
			CField& field = world.GetField(id_val_bc1);
			field.SetValue(-1.0,0,Fem::Field::VALUE,world,false);
		}
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetAxialSymmetry(false);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		// 方程式の係数の設定
		eqn_scalar.SetAlpha(0.00001);
		eqn_scalar.SetCapacity(1.0);
		eqn_scalar.SetAdvection(id_field_velo);
	}
	else if( iprob == 8 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 9 ){
		eqn_scalar.SetStationary(true);
	}
	else if( iprob == 10 ){
		eqn_scalar.SetSaveStiffMat(false);
	}
	else if( iprob == 11 )
	{
		////////////////
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			{
				std::vector<Com::CVector2D> vec_ary;
				vec_ary.push_back( Com::CVector2D(0.0,0.0) );
				vec_ary.push_back( Com::CVector2D(1.0,0.0) );
				vec_ary.push_back( Com::CVector2D(1.0,1.0) );
				vec_ary.push_back( Com::CVector2D(0.0,1.0) );
				cad_2d.AddPolygon( vec_ary );
			}
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(0.5,0.0));
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(0.5,1.0));
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
		}
		// メッシュを作る
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );	// メッシュで表される場のハンドルを得る
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		eqn_scalar.SetDomain_FieldElemAry(id_base,conv.GetIdEA_fromCad(2,Cad::LOOP),world);
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(true);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		// 方程式の設定
		dt = 0.02;
		eqn_scalar.SetTimeIntegrationParameter(dt);
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetCapacity(30.0);
		eqn_scalar.SetAdvection(false);
		// 境界条件の設定
		id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc0);
			field.SetValue("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world,true);
//			field.SetValue("floor(1+0.8*cos(2*PI*t+0.1))",0,world,true);
		}
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc1);
			field.SetValue(1.0,0,Fem::Field::VALUE,world,false);
		}

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_base,true,world) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 12 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 13 ){
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
	}
	else if( iprob == 14 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 15 ){		
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
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		eqn_scalar.SetDomain_Field(id_base,world);
		eqn_scalar.SetStationary(true);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		// 方程式の設定
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetCapacity(30.0);
		eqn_scalar.SetAdvection(false);
		eqn_scalar.SetSaveStiffMat(false);
		{
			Fem::Eqn::CEqn_Scalar2D eqn1 = eqn_scalar.GetEqnation(conv.GetIdEA_fromCad(2,Cad::LOOP));
			eqn1.SetAlpha(10.0);
			eqn_scalar.SetEquation(eqn1);
		}

		id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc0);
			field.SetValue("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world,true);
		}
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc1);
			field.SetValue(1.0,0,Fem::Field::VALUE,world,false);
		}

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 16 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 17 ){
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
	}
	else if( iprob == 18 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 19 ){
        Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.1) );
			vec_ary.push_back( Com::CVector2D(0.0,0.1) );
			cad_2d.AddPolygon(vec_ary);
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		eqn_scalar.SetDomain_Field(id_base,world);
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetSaveStiffMat(false);
		dt = 1;
		eqn_scalar.SetTimeIntegrationParameter(dt,0.5);
		eqn_scalar.SetAxialSymmetry(true);
		// 全体の方程式の係数設定
		eqn_scalar.SetAlpha(48.0);
		eqn_scalar.SetCapacity(480*7.86*1000);
		eqn_scalar.SetAdvection(0);
		eqn_scalar.SetSource(0);

		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		{
			CField& field = world.GetField(id_field_val);
			field.SetValue(500.0,0,Fem::Field::VALUE,world,false);
		}

		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		{
			CField& field = world.GetField(id_val_bc1);
			field.SetValue(0,0,Fem::Field::VALUE,world,false);
		}

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,0,500) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(mvp_trans);	// 視線座標変換行列の初期化
	}

	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}


////////////////////////////////////////////////////////////////

double cur_time = 0.0;
bool is_animation = true;
/*
void SetProjectionTransform()
{
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	if( mvp_trans.IsPers() ){	// 透視投影変換
		double fov_y,aspect,clip_near,clip_far;
		mvp_trans.GetPerspective(fov_y,aspect,clip_near,clip_far);
		::gluPerspective(fov_y,aspect,clip_near,clip_far);
	}
	else{	// 正規投影変換
		const double inv_scale = 1.0/mvp_trans.GetScale();
		const double asp = mvp_trans.GetWindowAspect();
		const double h_h = mvp_trans.GetHalfViewHeight()*inv_scale;
		const double h_w = mvp_trans.GetHalfViewHeight()*inv_scale*asp;
		::glOrtho(-h_w,h_w, -h_h, h_h, -10.0, 10.0);
	}
}

void SetModelViewTransform()
{
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	{	// 物体を平衡移動させる
		double x,y,z;
		mvp_trans.GetCenterPosition(x,y,z);
		::glTranslatef( x, y, z );
	}
	{	// 物体を回転させる
		double rot[16];
		mvp_trans.RotMatrix44Trans(rot);
		::glMultMatrixd(rot);
	}
	{	// 物体の中心を原点にする
		double x,y,z;
		mvp_trans.GetObjectCenter(x,y,z);
		::glTranslatef( -x, -y, -z );
	}
}
*/

// リサイズ時のコールバック関数
void myGlutResize(int w, int h)
{
	mvp_trans.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);

	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	::glutPostRedisplay();
}

// 描画時のコールバック関数
void myGlutDisplay(void)
{
//	::glClearColor(0.2, 0.7, 0.7 ,1.0);
	::glClearColor(1.0, 1.0, 1.0 ,1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(mvp_trans);

	if( is_animation ){
		cur_time += dt;
		world.FieldValueExec(cur_time);
		eqn_scalar.Solve(world);
		if( eqn_scalar.GetAry_ItrNormRes().size() > 0 ){
			std::cout << "Iter : " << eqn_scalar.GetAry_ItrNormRes()[0].first << " ";
			std::cout << "Res : " << eqn_scalar.GetAry_ItrNormRes()[0].second << std::endl;
		}
		drawer_ary.Update(world);
	}

	ShowFPS();
	drawer_ary.Draw();
	glutSwapBuffers();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	mvp_trans.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
	case 'q':
	case 'Q':
	case '\033':
		exit(0);  /* '\033' ? ESC ? ASCII ??? */
	case 'a':
		if( is_animation ){ is_animation = false; }
		else{ is_animation = true; }
		break;
	case ' ':
		SetNewProblem();
		break;
	}
	::glutPostRedisplay();
}

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
		mvp_trans.Fit();
		break;
	case GLUT_KEY_END :
		if( mvp_trans.IsPers() ) mvp_trans.SetIsPers(false);
		else{ mvp_trans.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	// glutの初期設定
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("DelFEM demo");

	// コールバック関数の設定
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);

    // 問題設定
	SetNewProblem();

	// メインループ
	glutMainLoop();
	return 0;
}
