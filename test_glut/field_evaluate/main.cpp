
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
#include <stdlib.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

using namespace Fem::Field;

Com::View::CCamera camera;
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
	static char s_fps[30];
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
	::glColor3d(1.0, 0.0, 0.0);
	strcpy(s_tmp,"DelFEM demo");
	RenderBitmapString(10,15, (void*)font, s_tmp);
	::glColor3d(0.0, 0.0, 1.0);
	strcpy(s_tmp,"Press \"Space\" key!");
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

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
//	camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	glutPostRedisplay();
}

View::CDrawerArrayField drawer_ary;

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
		drawer_ary.InitTrans(camera);
		camera.Fit();
		break;
	case GLUT_KEY_END :
		if( camera.IsPers() ) camera.SetIsPers(false);
		else{ camera.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

unsigned int id_field_val;
unsigned int id_base;
Fem::Field::CFieldWorld world;

double cur_time = 0.0;
bool is_animation = true;

void myGlutDisplay(void)
{
	::glClearColor(1.0, 1.0, 1.0 ,1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

	if( is_animation ){
		cur_time += 0.1;
		world.FieldValueExec(cur_time);
		drawer_ary.Update(world);
	}
	drawer_ary.Draw();
	ShowFPS();

	glutSwapBuffers();
}

void myGlutIdle(){
	glutPostRedisplay();
}

bool SetNewProblem()
{
	const unsigned int nprob = 13;
	static unsigned int iprob = 0;

	static int id_val_bc0=0, id_val_bc1=0, id_val_bc2=0;
	
	if( iprob == 0 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
			assert( world.IsIdField(id_field_val) );
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(10*sqrt(x^2+y^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 1 ){
		CField& val_field = world.GetField(id_field_val);
		val_field.SetValue("sin(2*PI*x-t)*sin(2*PI*y-t)", 0,Fem::Field::VALUE, world,true);
	}
	else if( iprob == 2 )
	{
		Msh::CMesher2D msh_2d;
		msh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
//		msh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		world.Clear();
		id_base = world.AddMesh( msh_2d );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(0.5*sqrt((x+1)^2+y^2)-0.1*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera);
	}
	else if( iprob == 3 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.07);
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d,1.0,0.07);
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(10*sqrt(x^2+y^2+z^2)-PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 4 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d,5.0,0.5);
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(1.0*sqrt(x^2+y^2+z^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 5 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
			cad_2d.AddVertex(Cad::EDGE, 1, Com::CVector2D(0.0, -0.5) );
			cad_2d.AddVertex(Cad::EDGE, 3, Com::CVector2D(0.0,  0.5) );
			cad_2d.ConnectVertex_Line(5,6);
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
			assert( world.IsIdField(id_field_val) );
			unsigned int id_field0 = world.GetPartialField(id_field_val,conv.GetIdEA_fromCad(1,Cad::LOOP));
			assert( world.IsIdField(id_field0) );
			CField& val_field0 = world.GetField(id_field0);
			val_field0.SetValue("sin(10*sqrt((x+0.5)^2+y^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
			unsigned int id_field1 = world.GetPartialField(id_field_val,conv.GetIdEA_fromCad(2,Cad::LOOP));
			CField& val_field1 = world.GetField(id_field1);
			val_field1.SetValue("sin(10*sqrt((x-0.5)^2+y^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 6 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.AddVertex(Cad::EDGE, 1, Com::CVector2D(0.0, -0.5) );
			cad_2d.AddVertex(Cad::EDGE, 3, Com::CVector2D(0.0,  0.5) );
			cad_2d.ConnectVertex_Line(5,6);
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		{
			id_field_val = world.MakeField_FieldElemAry(id_base,conv.GetIdEA_fromCad(1,Cad::LOOP),Fem::Field::SCALAR);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(10*sqrt((x+0.5)^2+y^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 7 ){
		unsigned int id_field_grad = world.MakeField_FieldElemDim(id_field_val,2,VECTOR2,VALUE,BUBBLE);
		{
			CField& field_grad = world.GetField(id_field_grad);
//			field_grad.SetGradient(id_field_val,world);
			field_grad.SetValue("0.1*sin(t)", 0,Fem::Field::VALUE, world,true);
			field_grad.SetValue("0.1*cos(t)", 1,Fem::Field::VALUE, world,true);
		}
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_grad,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 8 ){
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_tet.msh");
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		id_field_val = world.MakeField_FieldElemDim(id_base,3,SCALAR);
		{
			id_field_val = world.MakeField_FieldElemAry(id_base,1,Fem::Field::SCALAR);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(t+0.5*x)", 0,Fem::Field::VALUE, world,true);
		}
		unsigned int id_field_grad = world.MakeField_FieldElemDim(id_field_val,3,VECTOR3,VALUE,BUBBLE);
		{
			CField& field_grad = world.GetField(id_field_grad);
			field_grad.SetGradient(id_field_val,world,true);
			field_grad.ExportFile_Inp("grad_tet.inp",world);
		}
		drawer_ary.Clear();
//		drawer_ary.PushBack( new View::CDrawerFaceContour(id_field_val,world,-1.0,1.0) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_grad,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_grad,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 9 ){
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_hex.msh");
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		id_field_val = world.MakeField_FieldElemDim(id_base,3,SCALAR);
		{
			id_field_val = world.MakeField_FieldElemAry(id_base,1,Fem::Field::SCALAR);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(t+0.5*x)", 0,Fem::Field::VALUE, world,true);
		}
//		unsigned int id_field_grad = world.MakeField_AllRegion(VECTOR3,VALUE,BUBBLE);
		unsigned int id_field_grad = world.MakeField_FieldElemDim(id_field_val,3,VECTOR3,VALUE,BUBBLE);
		{
			CField& field_grad = world.GetField(id_field_grad);
			field_grad.SetGradient(id_field_val,world,true);
			field_grad.ExportFile_Inp("grad_hex.inp",world);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_grad,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_grad,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 10 )
	{
		Msh::CMesher2D msh_2d;
//		msh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		msh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		world.Clear();
		id_base = world.AddMesh( msh_2d );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR,VALUE,BUBBLE);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(x+y-0.1*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world) );
		drawer_ary.InitTrans( camera);
	}
	else if( iprob == 11 )
	{
		Msh::CMesher2D msh_2d;
//		msh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		msh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(msh_2d,5.0,1);
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,3,Fem::Field::SCALAR,VALUE,BUBBLE);
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(0.5*sqrt(x^2+y^2+z^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world) );
		drawer_ary.InitTrans( camera);
	}
	else if( iprob == 12 )
	{	// ÉoÉuÉãêﬂì_Çä‹ÇÒÇæï‚ä‘
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		{
			id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR,VALUE,CORNER|BUBBLE);
            std::cout << id_field_val << std::endl;
			assert( world.IsIdField(id_field_val) );
			CField& val_field = world.GetField(id_field_val);
			val_field.SetValue("sin(10*sqrt(x^2+y^2)-2*PI*t)", 0,Fem::Field::VALUE, world,true);
		}
		unsigned int id_field_vec;
		{
			id_field_vec = world.MakeField_FieldElemDim(id_base,2,Fem::Field::VECTOR2,VALUE,CORNER|BUBBLE);
			assert( world.IsIdField(id_field_vec) );
			CField& vec_field = world.GetField(id_field_vec);
			vec_field.SetValue("0.05*sin(t)", 0,Fem::Field::VALUE, world,true);
			vec_field.SetValue("0.05*cos(t)", 1,Fem::Field::VALUE, world,true);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerVector(id_field_vec,world) );
		drawer_ary.InitTrans( camera );
	}

	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
	  exit(0);
	  break;
  case 'a' :
	  is_animation = true;
	  break;
  case 's' :
	  if( is_animation ){
		cur_time += 0.1;
		world.FieldValueExec(cur_time);
		drawer_ary.Update(world);
		glutPostRedisplay();
	  }
	  is_animation = true;
	  break;
  case ' ':
	  SetNewProblem();
	  ::glMatrixMode(GL_PROJECTION);
	  ::glLoadIdentity();
	  Com::View::SetProjectionTransform(camera);
  default:
    break;
  }
}

int main(int argc,char* argv[])
{
	SetNewProblem();

	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 400);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("Field View");

	glutKeyboardFunc(myGlutKeyboard);
	glutIdleFunc(myGlutIdle);
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutSpecialFunc(myGlutSpecial);
	
	glutMainLoop();
	return 0;
}
