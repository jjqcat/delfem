////////////////////////////////////////////////////////////////
//                                                            //
//		Test of Hyper Elastic Solid 3D                        //
//                                                            //
//          Copy Rights (c) Nobuyuki Umetani 2008             //
//          e-mail : numetani@gmail.com                       //
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#endif
#define for if(0);else for

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>
#include <time.h>

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/cad_obj2d.h"
#include "delfem/mesh_primitive.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femeqn/eqn_linear_solid2d.h"
#include "delfem/femeqn/eqn_linear_solid3d.h"
#include "delfem/femeqn/eqn_hyper.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"

using namespace Fem::Ls;
using namespace Fem::Field;

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
int pressed_btn;
bool is_animation = true;

void RenderBitmapString(float x, float y, void *font,char *string)
{
  char *c;
  ::glRasterPos2f(x, y);
  for (c=string; *c != '\0'; c++) {
	  ::glutBitmapCharacter(font, *c);
  }
}

void ShowBackGround(){  
	::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
	::glLoadIdentity();
	::glMatrixMode(GL_PROJECTION);
    ::glPushMatrix();
	::glLoadIdentity();
    ::glDisable(GL_DEPTH_TEST);
    ::glBegin(GL_QUADS);
    ::glColor3d(0.2,0.7,0.7);
    ::glVertex3d(-1,-1,0);
    ::glVertex3d( 1,-1,0);
    ::glColor3d(1,1,1);
    ::glVertex3d( 1, 1,0);
    ::glVertex3d(-1, 1,0);
    ::glEnd();
    ::glEnable(GL_DEPTH_TEST);
	::glMatrixMode(GL_PROJECTION);
    ::glPopMatrix();
	::glMatrixMode(GL_MODELVIEW);
    ::glPopMatrix();
}

void ShowFPS(){
	static char s_fps[32];
	int* font=(int*)GLUT_BITMAP_8_BY_13;
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
	char s_tmp[30];

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
//	::glColor3d(1.0, 1.0, 0.0);
	::glColor3d(1.0, 0.0, 0.0);
	strcpy(s_tmp,"DelFEM demo");
	RenderBitmapString(10,15, (void*)font, s_tmp);
	::glColor3d(0.0, 0.0, 1.0);
	strcpy(s_tmp,"Press \"space\" key!");
	RenderBitmapString(120,15, (void*)font, s_tmp);
//	::glColor3d(1.0, 0.0, 0.0);
	::glColor3d(0.0, 0.0, 0.0);
	RenderBitmapString(10,30, (void*)font, s_fps);
//	::glEnable(GL_LIGHTING);
	::glEnable(GL_DEPTH_TEST);
	::glPopMatrix();
	::glMatrixMode(GL_PROJECTION);
	::glPopMatrix();
	::glMatrixMode(GL_MODELVIEW);
}

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
    if( pressed_btn == GLUT_MIDDLE_BUTTON ){
	    camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
    }
    if( pressed_btn == GLUT_RIGHT_BUTTON ){
	    camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
    }
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
    pressed_btn = button;
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
		Com::View::SetProjectionTransform(camera);
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
			if( camera.IsPers() ){
				const double tmp_fov_y = camera.GetFovY() + 10.0;
				camera.SetFovY( tmp_fov_y );
			}
		}
		else{
			const double tmp_scale = camera.GetScale() * 0.9;
			camera.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_PAGE_DOWN:
		if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
			if( camera.IsPers() ){
				const double tmp_fov_y = camera.GetFovY() - 10.0;
				camera.SetFovY( tmp_fov_y );
			}
		}
		else{
			const double tmp_scale = camera.GetScale() * 1.111;
			camera.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_HOME :
		camera.Fit();
		break;
	case GLUT_KEY_END :
		if( camera.IsPers() ) camera.SetIsPers(false);
		else{ camera.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

////////////////////////////////


Fem::Field::CFieldWorld world;
View::CDrawerArrayField drawer_ary;
double cur_time = 0.0;
const double dt = 0.06;
const double newmarkb_gamma = 0.59;
const double newmarkb_beta = 0.25*(0.5+newmarkb_gamma)*(0.5+newmarkb_gamma);
unsigned int id_field_disp;
unsigned int id_field_lambda;
Fem::Ls::CLinearSystem_Field ls;
LsSol::CPreconditioner_ILU prec;

void StepTime(){
    const double c1 = 200;
    const double c2 = 200;
    const double rho = 1.8;
    const double g[3] = { 0,0,0 };
    for(unsigned int iitr=0;iitr<2;iitr++){
        ls.InitializeMarge();
        Fem::Eqn::AddLinSys_Hyper3D_NonStatic_NewmarkBeta(
            dt, newmarkb_gamma, newmarkb_beta, ls,
		    c1, c2,
		    rho, g[0], g[1], g[2],    
		    id_field_disp, id_field_lambda, world, 
		    iitr==0);
/*        Fem::Eqn::AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta(
		    dt, newmarkb_gamma, newmarkb_beta, ls,
		    lambda, myu,
		    rho, g[0], g[1], g[2],
		    world,
		    id_field_disp );*/
        const double res = ls.FinalizeMarge();
        prec.SetValue(ls.m_ls);
        {
            double conv_ratio = 1.0e-6;
            unsigned int iteration = 400;
            LsSol::CLinearSystemPreconditioner lsp(ls.m_ls,prec);
//            Sol::Solve_PCG(conv_ratio,iteration, lsp);
            LsSol::Solve_PBiCGSTAB(conv_ratio,iteration, lsp);
            std::cout << "Iter NR : " << iitr << "   Res : " << res << "   Iter : " << iteration << "   Conv : " << conv_ratio << std::endl;
        }
        ls.UpdateValueOfField_NewmarkBeta(newmarkb_gamma,newmarkb_beta,dt, id_field_disp,  world,iitr==0);
        ls.UpdateValueOfField_NewmarkBeta(newmarkb_gamma,newmarkb_beta,dt, id_field_lambda,world,iitr==0);
    }
}

void myGlutDisplay(void)
{
	::glClearColor(0.2f, 0.7f, 0.7f,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);
    
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
    Com::View::SetProjectionTransform(camera);

	if( is_animation ){
		cur_time += dt;
		world.FieldValueExec(cur_time);
        StepTime();
		drawer_ary.Update(world);
	}

    ShowBackGround();
	drawer_ary.Draw();
	ShowFPS();



	::glutSwapBuffers();
}

void SetNewProblem()
{
	const unsigned int nprob = 1;	// –â‘è”
	static unsigned int iprob = 0;

	if( iprob == 0 ){
        cur_time = 0;
        Msh::CMesh_Primitive_Hexahedra mesh_3d(0.5, 4, 6,  1,8,8);
//        Msh::CMesh_Primitive_ThickCylinder mesh_3d(2,2.5,20, 1, 12, 6);
        world.Clear();
		const unsigned int id_base = world.AddMesh( mesh_3d );
        const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
        id_field_disp   = world.MakeField_FieldElemDim(id_base,3,VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
        id_field_lambda = world.MakeField_FieldElemDim(id_base,3,SCALAR, VALUE|VELOCITY|ACCELERATION,BUBBLE);
		unsigned int id_field_bc1 = world.GetPartialField(id_field_disp,conv.GetIdEA_fromMsh(2));
		{
			CField& bc1_field = world.GetField(id_field_bc1);
			bc1_field.SetValue("1*sin(2*t)", 0,Fem::Field::VALUE, world,true);	// bc1_field‚ÌyÀ•W‚É’PU“®‚ð’Ç‰Á
			bc1_field.SetValue("1*sin(t)"  , 1,Fem::Field::VALUE, world,true);	// bc1_field‚ÌyÀ•W‚É’PU“®‚ð’Ç‰Á
//			bc1_field.SetValue("-1*sin(t)" , 2,Fem::Field::VALUE, world,true);	// bc1_field‚ÌyÀ•W‚É’PU“®‚ð’Ç‰Á
		}
//		unsigned int id_field_bc2 = world.GetPartialField(id_field_disp,conv.GetIdEA_fromMsh(3));

        ls.Clear();
        ls.AddPattern_Field(id_field_disp,world);
        ls.AddPattern_Field(id_field_lambda,id_field_disp,world);
        ls.SetFixedBoundaryCondition_Field(id_field_bc1,world);
//        ls.SetFixedBoundaryCondition_Field(id_field_bc2,world);

        prec.Clear();
        prec.SetFillInLevel(0);
        prec.SetLinearSystem(ls.m_ls);

		// •`‰æƒIƒuƒWƒFƒNƒg‚Ì“o˜^
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);
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
	// Initialization of GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// Setting of call back function
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);
	
	SetNewProblem();	// set new problem
	glutMainLoop();		// enter main loop
	return 0;
}
