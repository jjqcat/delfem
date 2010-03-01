
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

#include "delfem/camera.h"

#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"

#include "delfem/femls/zlinearsystem.h"
#include "delfem/femls/zsolver_ls_iter.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/femeqn/eqn_helmholtz.h"

using namespace Fem::Field;
using namespace Fem::Ls;

Fem::Field::CFieldWorld world;
double dt = 0.02;
View::CDrawerArrayField drawer_ary;
Com::View::CCamera mvp_trans;
double mov_begin_x, mov_begin_y;

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
	::glDisable(GL_LIGHTING);
//	::glDisable(GL_DEPTH_TEST);
	::glColor3d(1.0, 0.0, 0.0);
	strcpy(s_tmp,"DelFEM demo");
	RenderBitmapString(10,15, (void*)font, s_tmp);
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
	const unsigned int nprob = 1;
	static unsigned int iprob = 0;
	
	if( iprob == 0 )	// ２次元問題の設定
	{
		////////////////
		Cad::CCadObj2D cad_2d;
		unsigned int id_v;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,2.0) );
			vec_ary.push_back( Com::CVector2D(0.0,2.0) );
			const unsigned int id_l = cad_2d.AddPolygon(vec_ary);
			id_v = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.05) );
		}
		// メッシュを作る
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		Fem::Field::CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		// 方程式の設定
		const unsigned int id_field_val = world.MakeField_FieldElemDim(id_base,2,ZSCALAR,VALUE,CORNER);
//		unsigned int id_field_bc0 = world.GetPartialField(id_field_val,conv.GetIdEA_fromCad(2,1));
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> aEA;
			aEA.push_back( conv.GetIdEA_fromCad(1,1) );
			aEA.push_back( conv.GetIdEA_fromCad(2,1) );
			aEA.push_back( conv.GetIdEA_fromCad(3,1) );
			aEA.push_back( conv.GetIdEA_fromCad(4,1) );
			id_field_bc1 = world.GetPartialField(id_field_val,aEA);
			CField& field = world.GetField(id_field_bc1);
			field.SetValue(1,0,world,false);
		}

		CZLinearSystem ls;
		CZPreconditioner_ILU prec;
		ls.AddPattern_Field(id_field_val,world);
//		ls.SetFixedBoundaryCondition_Field(id_field_bc0,world);
//		ls.SetFixedBoundaryCondition_Field(id_field_bc1,world);
		prec.SetFillInLevel(1);
		prec.SetLinearSystem(ls);

		double wave_length = 0.4;
		ls.InitializeMarge();
		Fem::Eqn::AddLinSys_Helmholtz(ls,wave_length,world,id_field_val);
		Fem::Eqn::AddLinSys_SommerfeltRadiationBC(ls,wave_length,world,id_field_bc1);
		double res = ls.FinalizeMarge();
		prec.SetValue(ls);

		{
			const unsigned int id_ea_v = conv.GetIdEA_fromCad(id_v,0);
			std::cout << id_ea_v << std::endl;
			const CElemAry& ea = world.GetEA(id_ea_v);
			const CElemAry::CElemSeg& es = ea.GetSeg(1);
			assert( ea.ElemType() == Fem::Field::POINT );
			unsigned int noes[1];
			es.GetNodes(0,noes);
			std::cout << noes[0] << std::endl;
			ls.GetResidualPtr(id_field_val,CORNER,world)->AddValue(noes[0],0,Com::Complex(1,0));
		}

		std::cout << "Residual : " << res << std::endl;
		{
			double tol = 1.0e-6;
			unsigned int iter = 2000;
//			Fem::Ls::Solve_CG(tol,iter,ls);
//			Fem::Ls::Solve_PCG(tol,iter,ls,prec);
			Fem::Ls::Solve_PCOCG(tol,iter,ls,prec);
//			Fem::Ls::Solve_CGNR(tol,iter,ls);
//			Fem::Ls::Solve_BiCGSTAB(tol,iter,ls);
//			Fem::Ls::Solve_BiCGStabP(tol,iter,ls,prec);
			std::cout << iter << " " << tol << std::endl;
		}
		ls.UpdateValueOfField(id_field_val,world,VALUE);

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-0.05,0.05) );
//		drawer_ary.PushBack( new View::CDrawerFaceContour(id_field_val,world) );
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

// リサイズ時のコールバック関数
void myGlutResize(int w, int h)
{
	mvp_trans.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);
	SetProjectionTransform();
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

	SetModelViewTransform();

//	ShowFPS();
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
		::glPopMatrix();
		::glMatrixMode(GL_MODELVIEW);
		::SetProjectionTransform();
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
	
	SetProjectionTransform();
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{

	SetNewProblem();

	// glutの初期設定
	glutInitWindowPosition(200,200);
	glutInitWindowSize(250, 250);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// コールバック関数の設定
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);

	// メインループ
	glutMainLoop();
	return 0;
}
