////////////////////////////////////////////////////////////////
//                                                            //
//		Moving LinearSolid                                    //
//                                                            //
//          Copy Rights (c) Nobuyuki Umetani 2007             //
//          e-mail : numetani@gmail.com                       //
////////////////////////////////////////////////////////////////

#pragma warning( disable : 4786 )
#define for if(0); else for

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>
#include <time.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"

#include "delfem/matvec/mat_blkcrs.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/ls/preconditioner.h"

#include "delfem/femeqn/eqn_linear_solid2d.h"
#include "delfem/femeqn/eqn_st_venant.h"
#include "delfem/femls/linearsystem_field.h"

#include "delfem/field_world.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"


using namespace Fem::Ls;
using namespace MatVec;
using namespace Fem::Field;

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


void ShowFPS(){
	static int frame, timebase;
	static char s[30];
	int time;
	int* font=(int*)GLUT_BITMAP_8_BY_13;

	frame++;
	time=glutGet(GLUT_ELAPSED_TIME);
	if (time - timebase > 500) {
		sprintf(s,"FPS:%4.2f",frame*1000.0/(time-timebase));
		timebase = time;
		frame = 0;
	}

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
	::glColor3d(1.0, 1.0, 0.0);
	RenderBitmapString(10,15, (void*)font, "DelFEM demo");
	::glColor3d(1.0, 1.0, 1.0);
	RenderBitmapString(10,30,(void *)font,s);
//	::glEnable(GL_LIGHTING);
	::glEnable(GL_DEPTH_TEST);
	::glPopMatrix();
	::glMatrixMode(GL_PROJECTION);
	::glPopMatrix();
	::glMatrixMode(GL_MODELVIEW);
}

void myGlutResize(int w, int h)
{
	mvp_trans.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
//	mvp_trans.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	mvp_trans.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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


double cur_time = 0.0;
const double dt = 0.02;
bool is_animation = false;

void SetProblem();
void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
	case 'q':
	case 'Q':
	case '\033':
		exit(0);  /* '\033' ? ESC ? ASCII ??? */
	case 'a':
		is_animation = !is_animation;
		break;
	case ' ':
		SetProblem();
		break;
	}
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
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
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}




////////////////////////////////


Fem::Field::CFieldWorld world;
View::CDrawerArrayField drawer_ary;

const double gravity[2] = { 0.0, -10.0 };

unsigned int id_disp;
CLinearSystem_Field ls;
LsSol::CPreconditioner_ILU prec;

void Solve(){	
	////////////////
	MatVec::CVector_Blk velo_pre;
	{
		Fem::Field::CField& field = world.GetField(id_disp);
		unsigned int id_ns_v = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
		unsigned int id_na = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		Fem::Field::CNodeAry& na = world.GetNA(id_na);		
		Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true,world,VELOCITY);
		unsigned int nblk = ns_v.Size();		
		unsigned int len = ns_v.Length();
		velo_pre.Initialize(nblk,len);
		na.GetValueFromNodeSegment(id_ns_v,velo_pre);
	}	
	for(unsigned int itr=0;itr<4;itr++){
		double res;
		{
			ls.InitializeMarge();
//            Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_BackwardEular			
			Fem::Eqn::AddLinSys_StVenant2D_NonStatic_BackwardEular
			(dt,
			 ls,
			 0.00, 4000.0,
			 1.0, gravity[0], gravity[1],
			 world, id_disp,
			 velo_pre,
			 itr==0);
			res = ls.FinalizeMarge();
		}
    {
			double conv_ratio = 1.0e-9;
			unsigned int max_iter = 1000;
      prec.SetValue(ls.m_ls);
      LsSol::CLinearSystemPreconditioner lsp(ls.m_ls,prec);
			LsSol::Solve_PCG(conv_ratio,max_iter,lsp);
			std::cout << max_iter << " " << conv_ratio << std::endl;
    }
		ls.UpdateValueOfField_BackwardEular(dt, id_disp  ,world,itr==0);
		std::cout << "iter : " << itr << "  Res : " << res << std::endl;
		if( res < 1.0e-6 ) break;
	}
//    getchar();
}

// 描画時のコールバック関数
void myGlutDisplay(void)
{
	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(mvp_trans);

	if( is_animation  ){
		cur_time += dt;
//		world.FieldValueExec(cur_time);
		Solve();
//		Solve2();
		drawer_ary.Update(world);
	}

	drawer_ary.Draw();
	ShowFPS();
	::glutSwapBuffers();
}



void SetProblem()
{
	const unsigned int nprob = 1;	// number of problems
	static int iprob = 0;

	if( iprob == 0 )
	{
    cur_time = 0;
		Cad::CCadObj2D cad_2d;
 		{	// define shape
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0,1.0) );
      vec_ary.push_back( Com::CVector2D(5,1.0) );
      vec_ary.push_back( Com::CVector2D(5,2.0) );
      vec_ary.push_back( Com::CVector2D(0,2.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.2) );
    const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		////////////////		
		id_disp = world.MakeField_FieldElemDim(id_base,2, VECTOR2, VALUE|VELOCITY|ACCELERATION, CORNER);
//		world.FieldValueExec(cur_time);
        unsigned int id_disp_fix0;
    {   // 境界条件の設定
      std::vector<unsigned int> aIdEAFix;
      //		aIdEAFix.push_back(3);
			aIdEAFix.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
      id_disp_fix0 = world.GetPartialField(id_disp,aIdEAFix);
      Fem::Field::CField& field = world.GetField(id_disp_fix0);
      field.SetValue("sin(t)",    0,Fem::Field::VALUE,world,true);
      field.SetValue("sin(0.5*t)",1,Fem::Field::VALUE,world,true);
    }
		
		// set linear system
    ls.Clear();
    ls.AddPattern_Field(id_disp,world);
		ls.SetFixedBoundaryCondition_Field(id_disp_fix0,world);
    // set preconditioner
    prec.SetFillInLevel(0);
    prec.SetLinearSystem(ls.m_ls);

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_disp,false,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_disp,false,world) );
	//	drawer_ary.PushBack( new View::CDrawerEdge(id_disp,true ,world) );
	}

	iprob++;
	if( iprob == nprob ){
		iprob = 0;
	}

}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	// Initialize GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// Define callback functions
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);
	
	// 問題の設定
	SetProblem();
	drawer_ary.InitTrans(mvp_trans);
	
	// メインループ
	glutMainLoop();
	return 0;
}
