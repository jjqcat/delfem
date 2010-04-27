
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
#include <stdio.h>
#include <cstdlib> //(exit)
#include <deque>

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/rigid/rigidbody.h"
#include "delfem/rigid/linearsystem_rigid.h"

Com::View::CCamera mvp_trans;
double mov_begin_x, mov_begin_y;
int press_button;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<Rigid::CRigidBody3D> aRB;
std::vector<Rigid::CConstraint*> apFix;

bool is_animation = false;

double cur_time = 0;
const double dt = 0.05;
const double newmark_gamma = 0.7;
const double newmark_beta = 0.25*(0.5+newmark_gamma)*(0.5+newmark_gamma);

const Com::CVector3D gravity( 0, 0, -1.0 );

void StepTime3()
{
    Ls::CLinearSystem_RigidBody_CRS2 ls(aRB,apFix);
    Ls::CPreconditioner_RigidBody_CRS2 prec;
    prec.SetLinearSystem(ls);
    ////////////////
    ls.InitializeMarge();
    ls.UpdateValueOfRigidSystem(aRB,apFix,   dt,newmark_gamma,newmark_beta,     true);
    double norm_res0;
	for(unsigned int itr=0;itr<10;itr++){
        ls.InitializeMarge();
        for(unsigned int irb=0;irb<aRB.size();irb++){
            aRB[irb].AddLinearSystem(ls,irb,      dt,newmark_gamma,newmark_beta,      gravity,itr==0 );
        }
        for(unsigned int ifix=0;ifix<apFix.size();ifix++){
            apFix[ifix]->AddLinearSystem(ls,ifix,         dt,newmark_gamma,newmark_beta,       aRB, itr==0);
        }
        ////////////////////////////////
        const double res = ls.FinalizeMarge();
        if( res < 1.0e-30 ) break;
        if( itr==0 ){ norm_res0 = res; }
        std::cout << "itr : " << itr << "     Residual : " << res << std::endl;
//        ls.ReSizeTmpVecSolver(1);
//        ls.COPY(-1,0);
//        ls.COPY(-1,-2);
        {
            ls.COPY(-1,-2);
            prec.SetValue(ls);
            prec.Solve( ls.GetVector(-2) );
        }
//        ls.MATVEC(-1.0, -2, 1, 0);
//        std::cout << "dot : " << ls.DOT(0,0) << std::endl;
//        prec.SetValue(ls);
//        ls2.Solve(prec);
/*        {
    		double conv_ratio = 1.0e-6;
		    unsigned int max_iter = 1000;
            prec.SetValue(ls);
            Ls::CLinearSystemPreconditioner_RigidBody_CRS2 lsp(ls,prec);
            Sol::Solve_PBiCGSTAB(conv_ratio,max_iter,lsp);
//            Sol::Solve_BiCGSTAB(conv_ratio,max_iter,ls);
//            std::cout << "       solver itr : " << max_iter << "  conv : " << conv_ratio << std::endl;
        }*/
        ////////////////////////////////
        ls.UpdateValueOfRigidSystem(aRB,apFix,   dt,newmark_gamma,newmark_beta,     false);
        if( res < norm_res0*1.0e-8 ) return;
	}
}

void RenderBitmapString(float x, float y, void *font,char *string)
{
  char *c;
  ::glRasterPos2f(x, y);
  for (c=string; *c != '\0'; c++) {
	  ::glutBitmapCharacter(font, *c);
  }
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

void DrawBackGround()
{
    // 背景描画
    ::glMatrixMode(GL_MODELVIEW);  // ModelView行列をクリア
    ::glPushMatrix();
    ::glLoadIdentity();
    ::glMatrixMode(GL_PROJECTION); // Projection行列をクリア
    ::glPushMatrix();
    ::glLoadIdentity();
    ::glDisable(GL_DEPTH_TEST);    // DepthTestを無効にする
    ::glBegin(GL_QUADS);
    ::glColor3d(0.2,0.7,0.7);   // 下の色を設定
    ::glVertex3d(-1,-1,0);  // 左下の点
    ::glVertex3d( 1,-1,0);  // 右下の点
    ::glColor3d(1,1,1);         // 上の色を設定
    ::glVertex3d( 1, 1,0);  // 右上の点
    ::glVertex3d(-1, 1,0);  // 左上の点
    ::glEnd();
    ::glEnable(GL_DEPTH_TEST);     // DepthTestを有効にする
    ::glMatrixMode(GL_PROJECTION); // Projection行列を元に戻す
    ::glPopMatrix();
    ::glMatrixMode(GL_MODELVIEW);  // ModelView行列を元に戻す
    ::glPopMatrix();
}

void myGlutResize(int w, int h)
{
	mvp_trans.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	glutPostRedisplay();
}

void DrawRigidBody(const Rigid::CRigidBody3D& r){
    const unsigned int imode = 1;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        ::glPushMatrix();
        ::glTranslated( 
            r.ini_pos_cg.x + r.disp_cg.x,
            r.ini_pos_cg.y + r.disp_cg.y,
            r.ini_pos_cg.z + r.disp_cg.z );
        {
            double rot0[16];
			r.GetInvRotMatrix44(rot0);
            ::glMultMatrixd(rot0);
        }
        ::glColor3d(0,0,0);
	    ::glutSolidCube(0.2);
    //    glutSolidTeapot(1.0);
    //    glutWireTeapot(0.2);
    //    glutWireOctahederon();
    //    glutWireDodecahedron();
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glColor3d(1,0,0);
        ::glVertex3d(0,0,0);
        ::glVertex3d(0.4,0,0);
        ::glColor3d(0,1,0);
        ::glVertex3d(0,0,0);
        ::glVertex3d(0,0.4,0);
        ::glColor3d(0,0,1);
        ::glVertex3d(0,0,0);
        ::glVertex3d(0,0,0.4);
        ::glEnd();
	    ::glPopMatrix();
    }
}


void DrawConstraint(const Rigid::CConstraint& c){
}

void myGlutDisplay(void)
{
	::glClearColor(0.2f, .7f, 0.7f, 1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(mvp_trans);

    if( is_animation ){        
        StepTime3();	// 剛体の運動を解く
        cur_time += dt;
		// 剛体エネルギー計算
        double eng = 0;
        for(unsigned int irb=0;irb<aRB.size();irb++){
            double e = 0;
            Com::CVector3D omg = aRB[irb].Omega;
            e += 0.5*(omg.x*omg.x*aRB[irb].mineatia[0]
                     +omg.y*omg.y*aRB[irb].mineatia[1]
                     +omg.z*omg.z*aRB[irb].mineatia[2]);
            Com::CVector3D velo = aRB[irb].velo_cg;
            e += 0.5*Com::Dot(velo,velo)*aRB[irb].mass;
            e -= Com::Dot(gravity,aRB[irb].disp_cg)*aRB[irb].mass;
            eng += e;
        }
        std::cout << "cur time " << cur_time << " " << eng << std::endl;
    }

    DrawBackGround();
/*  {
        ::glColor3d(0.7,0.7,0.7);
        ::glBegin(GL_QUADS);
        ::glVertex3d(-2,-2,0);	::glVertex3d( 2,-2,0);	::glVertex3d( 2, 2,0);	::glVertex3d(-2, 2,0);
        ::glEnd();
    } */
    {
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glColor3d(1,0,0);	::glVertex3d(0,0,0);	::glVertex3d(1,0,0);
        ::glColor3d(0,1,0);	::glVertex3d(0,0,0);	::glVertex3d(0,1,0);
        ::glColor3d(0,0,1);	::glVertex3d(0,0,0);	::glVertex3d(0,0,1);
        ::glEnd();
    }

    for(unsigned int irb=0;irb<aRB.size();irb++){
		DrawRigidBody(aRB[irb]);
//	    aRB[irb].Draw();
    }
    for(unsigned int ifix=0;ifix<apFix.size();ifix++){
//		apFix[ifix]->Draw(aRB);
    }

    ShowFPS();

	glutSwapBuffers();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

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

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
	press_button = button;
}

void SetNewProblem(){
    const unsigned int nprob = 9;
    static int iprob = 0;
  
    for(unsigned int icst=0;icst<apFix.size();icst++){ delete apFix[icst]; }
    apFix.clear(); 
    aRB.clear();
    cur_time = 0;

    if( iprob == 0 ){
        aRB.resize(1);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(0);
        pFix->SetIniPosFix(0,0,0);
        apFix.push_back( pFix );
    }
    else if( iprob == 1 ){
        double tot_len = 3.0;
        const unsigned int nRB = 6; 
        aRB.resize(nRB);
        const double div_len = tot_len / nRB;
        for(unsigned int irb=0;irb<nRB;irb++){
            aRB[irb].ini_pos_cg.SetVector(div_len*(irb+1),0,0);
            if( irb == 0 ){
                Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(irb);
                pFix->SetIniPosFix(0,0,0);
                apFix.push_back( pFix );
            }
            else{
                Rigid::CJoint_Spherical* pFix = new Rigid::CJoint_Spherical(irb-1,irb);
                pFix->SetIniPosJoint(div_len*(irb+0.5),0,0);
                apFix.push_back( pFix );
            }
        }
    }
    else if( iprob == 2 ){
        aRB.resize(3);  
        {
            aRB[0].ini_pos_cg.SetVector(1.0,0,0);
            aRB[1].ini_pos_cg.SetVector(1.0,1,0);
            aRB[2].ini_pos_cg.SetVector(2.0,1,0);
        }
        {
            Rigid::CFix_Hinge* pFix = new Rigid::CFix_Hinge(0);
            pFix->SetIniPosFix(0,0,0);
            pFix->SetAxis(0,1,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_Hinge* pFix = new Rigid::CJoint_Hinge(0,1);
            pFix->SetIniPosJoint(1,0.5,0);
            pFix->SetAxis(1,0,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_Spherical* pFix = new Rigid::CJoint_Spherical(1,2);
            pFix->SetIniPosJoint(1.5,1,0);
            apFix.push_back( pFix );
        }
    }
    else if( iprob == 3 ){
        double tot_len = 3.0;
        const unsigned int nRB = 3;
        aRB.resize(nRB);
        const double div_len = tot_len / nRB;
        for(unsigned int irb=0;irb<nRB;irb++){
            aRB[irb].ini_pos_cg.SetVector(div_len*(irb+1),0,0);
            if( irb == 0 ){
                Rigid::CFix_Hinge* pFix = new Rigid::CFix_Hinge(irb);
                pFix->SetIniPosFix(0,0,0);
                pFix->SetAxis(0.0,0.5,1);
                apFix.push_back( pFix );
            }
            else{
                Rigid::CJoint_Spherical* pFix = new Rigid::CJoint_Spherical(irb-1,irb);
                pFix->SetIniPosJoint(div_len*(irb+0.5),0,0);
                apFix.push_back( pFix );
            }
        }
    }
    else if( iprob == 4 ){
        aRB.resize(1);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        Rigid::CFix_HingeRange* pFix = new Rigid::CFix_HingeRange(0);
        pFix->SetIniPosFix(0,0,0);
        pFix->SetAxis(0,1,0);
        pFix->SetRange(-30,60);
        apFix.push_back( pFix );
    }
    else if( iprob == 5 ){
        aRB.resize(3);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        aRB[1].ini_pos_cg.SetVector(2.0,0,0);
        aRB[2].ini_pos_cg.SetVector(3.0,0,0);
        {
            Rigid::CFix_HingeRange* pFix = new Rigid::CFix_HingeRange(0);
            pFix->SetIniPosFix(0,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-30,60);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_Hinge* pFix = new Rigid::CJoint_Hinge(0,1);
            pFix->SetIniPosJoint(1.5,0,0);
            pFix->SetAxis(0,1,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_Hinge* pFix = new Rigid::CJoint_Hinge(1,2);
            pFix->SetIniPosJoint(2.5,0,0);
            pFix->SetAxis(0,1,0);
            apFix.push_back( pFix );
        }
    }
    else if( iprob == 6 ){
        aRB.resize(2);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        aRB[1].ini_pos_cg.SetVector(2.0,0,0);
        {
            Rigid::CFix_HingeRange* pFix = new Rigid::CFix_HingeRange(0);
            pFix->SetRange(-30,35);
//            Rigid::CFix_Hinge* pFix = new Rigid::CFix_Hinge(0);
            pFix->SetIniPosFix(0,0,0);
            pFix->SetAxis(0,1,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_HingeRange* pFix = new Rigid::CJoint_HingeRange(0,1);
            pFix->SetIniPosJoint(1.5,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-30,30);
            apFix.push_back( pFix );
        }
    }
    else if( iprob == 7 ){
        aRB.resize(4);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        aRB[1].ini_pos_cg.SetVector(2.0,0,0);
        aRB[2].ini_pos_cg.SetVector(3.0,0,0);
        aRB[3].ini_pos_cg.SetVector(4.0,0,0);
        {
            Rigid::CFix_HingeRange* pFix = new Rigid::CFix_HingeRange(0);
            pFix->SetRange(-30,35);
//            Rigid::CFix_Hinge* pFix = new Rigid::CFix_Hinge(0);
            pFix->SetIniPosFix(0,0,0);
            pFix->SetAxis(0,1,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_HingeRange* pFix = new Rigid::CJoint_HingeRange(0,1);
            pFix->SetIniPosJoint(1.5,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-30,30);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_HingeRange* pFix = new Rigid::CJoint_HingeRange(1,2);
            pFix->SetIniPosJoint(2.5,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-30,30);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_HingeRange* pFix = new Rigid::CJoint_HingeRange(2,3);
            pFix->SetIniPosJoint(3.5,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-30,30);
            apFix.push_back( pFix );
        }
    }
    if( iprob == 8 ){
        aRB.resize(1);
        aRB[0].ini_pos_cg.SetVector(0.3,0.0,1.0);
        aRB[0].Omega.z = 1.0;
        aRB[0].mineatia[2] = 10.0;
        aRB[0].mass = 0.001;
        Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(0);
        pFix->SetIniPosFix(0,0,0);
        apFix.push_back( pFix );
    }
    iprob++;
    if( iprob == nprob ) iprob = 0;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' は ESC の ASCII コード */
	  exit(0);
	  break;
  case 'a':
      is_animation = !is_animation;
	  break;
  case 's':
//      StepTime();
      break;
  case 'c':
      for(unsigned int irb=0;irb<aRB.size();irb++){
          aRB[irb].Clear();
      }
      for(unsigned int icst=0;icst<apFix.size();icst++){
          apFix[icst]->Clear();
      }
      break;
  case ' ':
      SetNewProblem();
  default:
    break;
  }
}

void myGlutSpecialFunc(int key, int x, int y){
    switch(key){
    case GLUT_KEY_PAGE_UP :     
        mvp_trans.SetScale( mvp_trans.GetScale()*0.9    );
        break;
    case GLUT_KEY_PAGE_DOWN :
        mvp_trans.SetScale( mvp_trans.GetScale()*1.1111 );
        break;
    }
	::glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
}

int main(int argc,char* argv[])
{
	// Initailze GLUT
	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 300);
	::glutInit(&argc, argv);	
    ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("Cad View");

	// Set callback function
	::glutMotionFunc(myGlutMotion);
	::glutMouseFunc(myGlutMouse);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutKeyboardFunc(myGlutKeyboard);
    ::glutSpecialFunc(myGlutSpecialFunc);
	::glutIdleFunc(myGlutIdle);

    mvp_trans.SetRotationMode(Com::View::ROT_3D);
//    mvp_trans.SetIsPers(true);
    {
        Com::CBoundingBox bb(-2,2,-2,2,-2,2);
        mvp_trans.SetObjectBoundingBox(bb);
        mvp_trans.Fit(bb);
    }
    SetNewProblem();

	// Enter main loop
	::glutMainLoop();
	return 0;
}
