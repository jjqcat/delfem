
#include <QtGui>
#include <QtOpenGL>

#include <math.h>

#include "glwidget.h"

#include "delfem/camera.h"
#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"	// 有限要素法離散場可視化クラス
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/drawer_field_image_based_flow_vis.h"
#include "delfem/drawer_field_streamline.h"

#include "delfem/eqnsys_fluid.h"		// 有限要素法流体ソルバクラスCEqnSystem_Fluid2D


using namespace Fem::Field;

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{
    cur_time = 0;
    dt = 0.05;
    SetNewProblem();

    QTimer *timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(StepTime()));
    timer->start(20);
}

GLWidget::~GLWidget()
{
    makeCurrent();
}

void GLWidget::initializeGL()
{

}

void GLWidget::paintGL()
{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    Com::View::SetModelViewTransform(camera);

    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    Com::View::SetProjectionTransform(camera);

    drawer_ary.Draw();
}

void GLWidget::resizeGL(int w, int h)
{
    camera.SetWindowAspect((double)w/h);
    ::glViewport(0, 0, w, h);
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    Com::View::SetProjectionTransform(camera);
    updateGL();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
}

void GLWidget::StepTime()
{
    cur_time += dt;
    world.FieldValueExec(cur_time);
    fluid.Solve(world);
    world.FieldValueDependExec();
    drawer_ary.Update(world);
    updateGL();
}

void GLWidget::SetNewProblem()
{
	const unsigned int nprob = 15;	// 問題数
	static unsigned int iprob = 0;
	
	if( iprob == 0 )	// キャビティフロー問題，定常ストークス流れ
	{
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			unsigned int id_l = cad_2d.AddPolygon( vec_ary );
			cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.0,0.0));
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		
		unsigned int id_field_press = fluid.GetIdField_Press();
		//		unsigned int id_field_press_bc0 = world.GetPartialField(id_field_press,10);
		//		fluid.AddFixField(id_field_press_bc0,world);
		
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.5*sin(0.05*t)", 0,Fem::Field::VELOCITY, world,true);
			//			bc0_field_velo.SetVelocity("0.1", 0, world,true);
		}
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		fluid.SetRho(0.1);
		fluid.SetMyu(0.0002);
		fluid.SetStokes();
		fluid.SetIsStationary(true);
		dt = 0.5;
		fluid.SetTimeIntegrationParameter(dt);
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
                drawer_ary.InitTrans( camera );
	}
	else if( iprob == 1 )	// キャビティーフロー問題，非定常ストークス流れ
	{
		fluid.SetIsStationary(false);
	}
	else if( iprob == 2 )	// キャビティーフロー問題，非定常ストークス流れ(Rhoを大きく)
	{
		fluid.SetRho(0.5);
	}
	else if( iprob == 3 )	// キャビティーフロー問題，非定常Naiver-Stokes流れ
	{
		fluid.SetRho(0.02);
		fluid.SetMyu(0.00001);
		fluid.SetNavierStokes();
	}
	else if( iprob == 4 )	// キャビティーフロー問題，バブル補間，定常ストークス流れ
	{
		fluid.Clear();
		fluid.SetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		fluid.SetStokes();
		fluid.SetIsStationary(true);
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		unsigned int id_field_press = fluid.GetIdField_Press();
        std::cout << "press : " << id_field_press << std::endl;
		unsigned int id_field_press_bc0 = world.GetPartialField(id_field_press,conv.GetIdEA_fromCad(1,Cad::VERTEX));
		fluid.AddFixField(id_field_press_bc0,world);
		
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
        std::cout << "velo : " << id_field_velo << std::endl;
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.5*sin(0.1*t)", 0,Fem::Field::VELOCITY, world,true);
			//			bc0_field_velo.SetVelocity("0.1", 0, world,true);
		}
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
                drawer_ary.InitTrans( camera );
	}
	else if( iprob == 5 )	// L字形
	{
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.5) );
			vec_ary.push_back( Com::CVector2D(0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(0.5, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		fluid.UnSetInterpolationBubble();
		fluid.Clear();
		fluid.UpdateDomain_Field(id_base,world);
		fluid.SetTimeIntegrationParameter(dt);
		
		unsigned int id_field_press = fluid.GetIdField_Press();
		//		unsigned int id_field_press_bc0 = world.GetPartialField(id_field_press,10);
		//		fluid.AddFixField(id_field_press_bc0,world);
		
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.1*sin(0.1*t)", 0,Fem::Field::VELOCITY, world,true);
			//			bc0_field_velo.SetVelocity("0.1", 0, world,true);
		}
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(6,Cad::EDGE) );
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		fluid.SetRho(0.1);
		fluid.SetMyu(0.0002);
		//		fluid.UnSetStationary(world);
		fluid.SetIsStationary(true);
		fluid.SetStokes();
		fluid.SetTimeIntegrationParameter(dt);
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
                drawer_ary.InitTrans( camera );
	}
	else if( iprob == 6 )
	{
		fluid.SetIsStationary(false);
		fluid.SetNavierStokes();
		fluid.SetRho(1.3);
		fluid.SetMyu(0.0002);
	}
	else if( iprob == 7 )
	{
		Cad::CCadObj2D cad_2d;
		{	// 正方形にが２つに分割
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(0.5,0.0));
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(0.5,1.0));
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
		}
		
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.03) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		fluid.Clear();
		fluid.UpdateDomain_FieldElemAry(id_base, conv.GetIdEA_fromCad(2,Cad::LOOP) ,world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		
		unsigned int id_field_bc1 = fluid.AddFixElemAry( conv.GetIdEA_fromCad(3,Cad::EDGE) ,world);
		{
			Fem::Field::CField& field = world.GetField(id_field_bc1);
			field.SetValue("0.3*sin(0.5*t)", 1,Fem::Field::VELOCITY, world,true);
		}
		unsigned int id_field_bc2 = fluid.AddFixElemAry( conv.GetIdEA_fromCad(5,Cad::EDGE) ,world);
		
		fluid.SetRho(0.1);
		fluid.SetMyu(0.0002);
		fluid.SetStokes();
		fluid.SetIsStationary(true);
		fluid.SetTimeIntegrationParameter(dt);
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_base,true,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
                drawer_ary.InitTrans( camera );
	}
	else if( iprob == 8 ){
		fluid.SetIsStationary(false);
	}
	else if( iprob == 9 ){
		fluid.SetNavierStokes();
	}
	else if( iprob == 10 ){	// カルマン渦列
		Cad::CCadObj2D cad_2d;
		{	// 正方形に矩形の穴
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.6) );
			vec_ary.push_back( Com::CVector2D(0.0,0.6) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.2));
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.2));
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.4));
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.4));
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}
		
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		
		fluid.Clear();
		fluid.UpdateDomain_FieldElemAry(id_base,conv.GetIdEA_fromCad(1,Cad::LOOP),world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		
		unsigned int id_field_bc1 = fluid.AddFixElemAry( conv.GetIdEA_fromCad(4,Cad::EDGE) ,world);
		{
			Fem::Field::CField& field = world.GetField(id_field_bc1);
			field.SetValue(0.1,0,Fem::Field::VELOCITY,world,true);
		}      
		{
			std::vector<unsigned int> aIdEAFixVelo;
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			unsigned int id_field_bc2 = fluid.AddFixElemAry(aIdEAFixVelo,world);
		}
		{
			std::vector<unsigned int> aIdEAFixVelo;
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(5,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(6,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(7,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(8,Cad::EDGE) );
			unsigned int id_field_bc3 = fluid.AddFixElemAry(aIdEAFixVelo,world);
		}
		dt = 0.13;
		fluid.SetRho(200.0);
		fluid.SetMyu(0.0001);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
                drawer_ary.InitTrans( camera );
	}
	else if( iprob == 11 )	// 上下で分離している問題
	{
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.5) );
			vec_ary.push_back( Com::CVector2D(1.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 0.5) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.ConnectVertex_Line(3,6);
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		
		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.1*sin(0.1*t)", 0,Fem::Field::VELOCITY, world,true);
		}
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(6,Cad::EDGE) );
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		dt = 0.8;
		fluid.SetRho(1);
		fluid.SetMyu(0.0002);
		fluid.SetIsStationary(false);
		//		fluid.SetStationary(world);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
                drawer_ary.InitTrans( camera );
	}
    else if( iprob == 12 ){ // バックステップ流れ
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D( 0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D( 1.4, 0.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 0.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 1.0) );
			vec_ary.push_back( Com::CVector2D( 1.4, 1.0) );
			vec_ary.push_back( Com::CVector2D(-0.5, 1.0) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.7) );
			vec_ary.push_back( Com::CVector2D( 0.0, 0.7) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.ConnectVertex_Line(2,5);	
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		
		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);
		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.2", 0,Fem::Field::VELOCITY, world,true);
		}
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(5,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(7,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(8,Cad::EDGE) );
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		fluid.SetRho(5);
		fluid.SetMyu(0.0002);
		fluid.SetIsStationary(false);
		//		fluid.SetStationary(world);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		{
         Fem::Eqn::CEqn_Fluid2D eqn_fluid = fluid.GetEquation( conv.GetIdEA_fromCad(2,Cad::LOOP) );
			eqn_fluid.SetMyu(0.01);
			fluid.SetEquation(eqn_fluid);
		}
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
                drawer_ary.InitTrans( camera );
	}
	if( iprob == 13 ){ // 流体の中に力が働いている領域
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D( 0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D( 3.0, 0.0) );
			vec_ary.push_back( Com::CVector2D( 3.0, 3.0) );
			vec_ary.push_back( Com::CVector2D( 0.0, 3.0) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary );
			vec_ary.clear();
			vec_ary.push_back( Com::CVector2D( 1.0, 1.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 1.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 1.5) );
			vec_ary.push_back( Com::CVector2D( 1.0, 1.5) );
			cad_2d.AddPolygon( vec_ary, id_l );
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		
		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0;
		{
			std::vector<unsigned int> id_ea_bc;
			id_ea_bc.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_field_bc0 = fluid.AddFixElemAry(id_ea_bc,world);
		}
		fluid.SetRho(5);
		fluid.SetMyu(0.005);
		fluid.SetIsStationary(false);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		{
         Fem::Eqn::CEqn_Fluid2D eqn_fluid = fluid.GetEquation( conv.GetIdEA_fromCad(2,Cad::LOOP) );
			eqn_fluid.SetBodyForce(0.0,0.5);
			fluid.SetEquation(eqn_fluid);
		}
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world);
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
                drawer_ary.InitTrans( camera );
	}
	else if( iprob == 14 )
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l;
		unsigned int id_e1, id_e2, id_e3,id_e4,id_e5,id_e6;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D( 0.0,0.0) );
			vec_ary.push_back( Com::CVector2D( 0.5,0.0) );
			vec_ary.push_back( Com::CVector2D( 2.0,0.0) );
			vec_ary.push_back( Com::CVector2D( 2.0,1.0) );
			vec_ary.push_back( Com::CVector2D( 0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary );
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.5) );
			id_e1 = cad_2d.ConnectVertex_Line(2,id_v1);
			unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.0,0.3) );
			unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.0,0.9) );
			id_e2 = cad_2d.ConnectVertex_Line(id_v2,id_v3);
			unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.5,0.4) );
			unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.5,0.1) );
			unsigned int id_v6 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.5,0.7) );
			unsigned int id_v7 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.2,0.4) );
			unsigned int id_v8 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.8,0.4) );
			id_e3 = cad_2d.ConnectVertex_Line(id_v4,id_v5);
			id_e4 = cad_2d.ConnectVertex_Line(id_v4,id_v6);
			id_e5 = cad_2d.ConnectVertex_Line(id_v4,id_v7);
			id_e6 = cad_2d.ConnectVertex_Line(id_v4,id_v8);
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.05);
		
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh(mesh_2d);
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // ID変換クラス
		unsigned int id_base2 = 0;
		{
			std::vector<unsigned int> mapVal2Co;
			std::vector< std::vector<int> > aLnods;
			{
				std::vector<unsigned int> aIdMsh_Inc;
				aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l,Cad::LOOP) );
				std::vector<unsigned int> aIdMshBar_Cut;
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e1,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e2,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e5,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e6,Cad::EDGE) );
				mesh_2d.GetClipedMesh(aLnods,mapVal2Co, aIdMsh_Inc,aIdMshBar_Cut);
			}
			std::vector<unsigned int> aIdEA_Inc;
			aIdEA_Inc.push_back( conv.GetIdEA_fromCad(1,Cad::LOOP) );
			id_base2 = world.SetCustomBaseField(id_base,aIdEA_Inc,aLnods,mapVal2Co);
		}
		
		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_FieldVeloPress(id_base,id_base2,world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0;
		{
			std::vector<unsigned int> id_ea_bc;
			id_ea_bc.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e1,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e2,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e3,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e4,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e5,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e6,Cad::EDGE) );
			id_field_bc0 = fluid.AddFixElemAry(id_ea_bc,world);
		}
		unsigned int id_field_bc1 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(5,Cad::EDGE),world);
		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc1);
			bc0_field_velo.SetValue("0.1*sin(t*PI*0.1+0.01)", 0,Fem::Field::VELOCITY, world,true);
		}
		
		dt = 0.8;
		fluid.SetRho(10);
		fluid.SetMyu(0.005);
		fluid.SetIsStationary(false);
		//		fluid.SetStokes();
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(     id_field_velo, world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
                drawer_ary.InitTrans( camera );
	}
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
        Com::View::SetProjectionTransform(camera);
	//	::glutPostRedisplay();
	
	iprob++;
	if( iprob == nprob ){ iprob = 0; }
}
