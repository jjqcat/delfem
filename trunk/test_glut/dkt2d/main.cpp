

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

#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"

#include "delfem/femls/linearsystem_field.h"

#include "delfem/femeqn/eqn_dkt.h"

using namespace Fem::Ls;
using namespace Fem::Field;

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
	int* font = (int*)GLUT_BITMAP_8_BY_13;
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
	strcpy(s_tmp,"DelFEM domo");
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

void myGlutIdle(){	// アイドル時のコールバック関数
	glutPostRedisplay();
}

Com::View::CCamera mvp_trans;

void myGlutResize(int w, int h)
{	// ウィンドウサイズ変更時のコールバック関数
	mvp_trans.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	glutPostRedisplay();
}

View::CDrawerArrayField drawer_ary;
Fem::Field::CFieldWorld world;
double mov_begin_x, mov_begin_y;
bool is_animation = true;

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	mvp_trans.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
//	mvp_trans.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	::glutPostRedisplay();
}

void myGlutDisplay(void)
{	// 描画時のコールバック関数
	::glClearColor(1.0, 1.0, 1.0, 1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1, 4.0 );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(mvp_trans);

	drawer_ary.Draw();
	ShowFPS();

	glutSwapBuffers();
}


void SetNewProblem();
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
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(mvp_trans);
		break;
	}
	::glutPostRedisplay();
}


void SetNewProblem()
{
	const unsigned int nprob = 2;
	static unsigned int iprob = 0;

	static int id_val_bc0=0, id_val_bc1=0, id_val_bc2=0;
	
	if( iprob == 0 ){	// ２次元問題の設定
		Cad::CCadObj2D cad_2d;
		{	// 形を作る
            std::vector<Com::CVector2D> vec_ary;
            vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
            vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
            vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
            vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			cad_2d.AddPolygon( vec_ary );
            cad_2d.AddVertex(Cad::LOOP,1,Com::CVector2D(0.0,0.0));
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		unsigned int id_field_rot = world.MakeField_FieldElemDim(id_base,2,VECTOR2,VALUE,CORNER);
		unsigned int id_field_deflect = world.MakeField_FieldElemDim(id_base,2,SCALAR,VALUE,CORNER);

		Fem::Ls::CLinearSystem_Field ls;
		LsSol::CPreconditioner_ILU prec;
		ls.AddPattern_Field(id_field_deflect,world);
		ls.AddPattern_Field(id_field_rot,id_field_deflect,world);
		prec.SetLinearSystem(ls.m_ls);

		unsigned int id_field_rot_fix0 = world.GetPartialField(id_field_rot,conv.GetIdEA_fromCad(2,1));
		unsigned int id_field_rot_fix1 = world.GetPartialField(id_field_rot,conv.GetIdEA_fromCad(4,1));
		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix1,world);

		unsigned int id_field_def_fix0 = world.GetPartialField(id_field_deflect,conv.GetIdEA_fromCad(2,1));
		unsigned int id_field_def_fix1 = world.GetPartialField(id_field_deflect,conv.GetIdEA_fromCad(4,1));
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix1,world);

		{
			Fem::Field::CField& field = world.GetField(id_field_rot_fix1);
			field.SetValue(-1,1,world,false);
		}
		ls.InitializeMarge();
		Fem::Eqn::AddLinearSystem_DKT2D_Static(ls,world,id_field_deflect,id_field_rot);
		double res = ls.FinalizeMarge();
		prec.SetValue(ls.m_ls);

		std::cout << "Residual : " << res << std::endl;
		{
			double tol = 1.0e-6;
			unsigned int iter = 10000;
			Sol::Solve_CG(tol,iter,ls);
//			Fem::Sol::Solve_PCG(tol,iter,ls,prec);
			std::cout << iter << " " << tol << std::endl;
		}
		ls.UpdateValueOfField(id_field_deflect,world,VALUE);
		ls.UpdateValueOfField(id_field_rot,world,VALUE);


		// 描画オブジェクトの登録
		drawer_ary.Clear();
	//	drawer_ary.PushBack( new View::CDrawerVector(id_field_rot,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,true,world) );
		drawer_ary.InitTrans( mvp_trans );
	}
	else if( iprob == 1 ){
		Msh::CMesher2D msh;
		msh.ReadFromFile_GiDMsh("../input_file/rect_tri.msh");
		world.Clear();
		const unsigned int id_base = world.AddMesh( msh );

		unsigned int id_field_rot = world.MakeField_FieldElemDim(id_base,2,  VECTOR2,VALUE,CORNER);
		unsigned int id_field_deflect = world.MakeField_FieldElemDim(id_base,2,  SCALAR,VALUE,CORNER);

		unsigned int id_field_rot_fix0 = world.GetPartialField(id_field_rot,3);
		unsigned int id_field_rot_fix1 = world.GetPartialField(id_field_rot,4);
		unsigned int id_field_def_fix0 = world.GetPartialField(id_field_deflect,3);
		unsigned int id_field_def_fix1 = world.GetPartialField(id_field_deflect,4);
		{
			Fem::Field::CField& field = world.GetField(id_field_def_fix1);
			field.SetValue(-3,0,world,false);
		}

		Fem::Ls::CLinearSystem_Field ls;
		LsSol::CPreconditioner_ILU prec;
		ls.AddPattern_Field(id_field_deflect,world);
		ls.AddPattern_Field(id_field_rot,id_field_deflect,world);
		prec.SetLinearSystem(ls.m_ls);

		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix1,world);
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix1,world);
		ls.InitializeMarge();
		Fem::Eqn::AddLinearSystem_DKT2D_Static(ls,world,id_field_deflect,id_field_rot);
		double res = ls.FinalizeMarge();
		prec.SetValue(ls.m_ls);

		std::cout << "Residual : " << res << std::endl;
		{
			double tol = 1.0e-6;
			unsigned int iter = 10000;
			Sol::Solve_CG(tol,iter,ls);
//			Fem::Sol::Solve_PCG(tol,iter,ls,prec);
			std::cout << iter << " " << tol << std::endl;
		}
		ls.UpdateValueOfField(id_field_deflect,world,VALUE);
		ls.UpdateValueOfField(id_field_rot,world,VALUE);


		// 描画オブジェクトの登録
		drawer_ary.Clear();
	//	drawer_ary.PushBack( new View::CDrawerVector(id_field_rot,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,true,world) );
		drawer_ary.InitTrans( mvp_trans );
	}

	iprob++;
	if( iprob == nprob ) iprob=0;
}

int main(int argc,char* argv[])
{
	SetNewProblem();

	// GLUTの初期設定
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("DelFEM Demo");

	// コールバック関数の設定
	glutIdleFunc(myGlutIdle);
	glutKeyboardFunc(myGlutKeyboard);
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutSpecialFunc(myGlutSpecial);;
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	
	// メインループ
	glutMainLoop();
	return 0;
}
