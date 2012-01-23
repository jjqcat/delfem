
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#pragma warning( disable : 4786 ) 
#define for if(0);else for

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>

#if defined(__APPLE__) && (__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/cad_obj2d.h"
#include "delfem/drawer_cad.h"

#include "delfem/mesher2d.h"
#include "delfem/drawer_msh.h"
#include "delfem/drawer_cad.h"

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d_move.h"

bool is_draw_msh = true;

double mov_begin_x, mov_begin_y;
//int press_button;
int imodifier;

unsigned int id_part_cad;
Cad::CAD_ELEM_TYPE itype_part_cad;
int ictrl_cad;

// 0:drag
// 1:drag curve
// 2:add point
// 3:delete point
// 4:smooth curve
int imode = 0;


Com::View::CCamera_2D camera;
Cad::CCadObj2D_Move cad_2d;
Cad::View::CDrawer_Cad2D drawerCAD;
bool is_updated_cad = false;

std::vector< Com::CVector2D > aVecStrok;

void myGlutPassiveMotion( int x, int y ){
  
  { // hilight cad element
    const unsigned int size_buffer = 128;
    unsigned int select_buffer[size_buffer];
    Com::View::PickPre(size_buffer,select_buffer,   x,y,  5,5,  camera);
    drawerCAD.DrawSelection(0);
    std::vector<Com::View::SSelectedObject> aSelecObj = Com::View::PickPost(select_buffer, x,y, camera);
    drawerCAD.ClearModeShow(1);
    if( aSelecObj.size() != 0 ){
      unsigned int id_part_cad0 = 0;
      Cad::CAD_ELEM_TYPE itype_part_cad0 = Cad::NOT_SET;
      int ictrl_cad0;
      drawerCAD.GetCadPartID(aSelecObj[0].name,itype_part_cad0,id_part_cad0,ictrl_cad0);
      if( drawerCAD.GetModeShow(itype_part_cad0,id_part_cad0) != 2 ){
        drawerCAD.SetModeShow(itype_part_cad0, id_part_cad0, 1 );
      }
    }
  }  
}

void myGlutMouse(int button, int state, int x, int y)
{
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/(double)win_w;
	mov_begin_y = (win_h-2.0*y)/(double)win_h;
	imodifier = glutGetModifiers();
  
  const unsigned int id_part_cad0 = id_part_cad;
  Cad::CAD_ELEM_TYPE itype_part_cad0 = itype_part_cad;
  { // hilight the picked element
    const unsigned int size_buffer = 128;
    unsigned int select_buffer[size_buffer];
    Com::View::PickPre(size_buffer,select_buffer, x,y, 5,5, camera);
    drawerCAD.DrawSelection(0);
    std::vector<Com::View::SSelectedObject> aSelecObj = Com::View::PickPost(select_buffer, x,y, camera );
    if( state == GLUT_DOWN ){
      drawerCAD.ClearModeShow(2);
    }
    id_part_cad = 0;
    if( aSelecObj.size() != 0 ){
      drawerCAD.GetCadPartID(aSelecObj[0].name,itype_part_cad,id_part_cad,ictrl_cad);
      drawerCAD.SetModeShow(itype_part_cad,id_part_cad,2);
    }
    else{
      itype_part_cad = Cad::NOT_SET;
      id_part_cad = 0;
    }
  }

	if( state == GLUT_DOWN ){
		if( imode == 1 ){ 
      if( itype_part_cad == Cad::EDGE && cad_2d.IsElemID(itype_part_cad,id_part_cad) ){
        const Com::CVector3D& oloc = camera.ProjectionOnPlane(mov_begin_x,mov_begin_y);  
        if( cad_2d.GetEdgeCurveType(id_part_cad) == Cad::CURVE_POLYLINE ){
          bool res = cad_2d.PreCompDragPolyline(id_part_cad,Com::CVector2D(oloc.x,oloc.y));
        }      
        return; 
      }
    }
		if( imode == 2 ){   // add point
			aVecStrok.clear();
			if( itype_part_cad == Cad::VERTEX ){ return; }
			const Com::CVector3D& oloc = camera.ProjectionOnPlane(mov_begin_x,mov_begin_y);
			const Com::CVector2D& v0 = Com::CVector2D(oloc.x,oloc.y);
			unsigned int id_v0 = cad_2d.AddVertex(itype_part_cad,id_part_cad, v0 ).id_v_add;
			itype_part_cad = Cad::VERTEX;
			id_part_cad = id_v0;
			if( !cad_2d.IsElemID(Cad::VERTEX,id_v0) ){ aVecStrok.clear(); return; }
      drawerCAD.UpdateCAD_TopologyGeometry(cad_2d);
		}
		if( imode == 3 ){   // delete point
			if( cad_2d.IsElemID(itype_part_cad,id_part_cad) ){
				bool iflag = cad_2d.RemoveElement(itype_part_cad,id_part_cad);
				if( iflag ){
          drawerCAD.UpdateCAD_TopologyGeometry(cad_2d);
				}
			}
			imode = 0;
		}
	}
	if( state == GLUT_UP ){
		if( imode == 2 ){
//			imode = 0;
			if( itype_part_cad0 != Cad::VERTEX || !cad_2d.IsElemID(Cad::VERTEX,id_part_cad0 ) ){ 
				aVecStrok.clear();
				return;
			}
			unsigned int id_v0 = id_part_cad0;
			unsigned int id_v1 = 0;
			if( itype_part_cad == Cad::VERTEX ){
				if( !cad_2d.IsElemID(Cad::VERTEX,id_part_cad) ){ aVecStrok.clear(); return; }
				id_v1 = id_part_cad;
			}
			else{
				const Com::CVector2D& v1 = aVecStrok[aVecStrok.size()-1];
				id_v1 = cad_2d.AddVertex(itype_part_cad, id_part_cad, v1 ).id_v_add;
				if( !cad_2d.IsElemID(Cad::VERTEX,id_v1) ){ aVecStrok.clear(); return; }
			}
			Cad::CEdge2D e(id_v0,id_v1);
			{
//        e.itype = Com::CURVE_POLYLINE;
				const unsigned int n = aVecStrok.size();
				const Com::CVector2D& pos = cad_2d.GetVertexCoord(e.GetIdVtx(true));
				const Com::CVector2D& poe = cad_2d.GetVertexCoord(e.GetIdVtx(false));
				e.SetVtxCoords(pos,poe);
				const double sqlen = Com::SquareLength(poe-pos);
				const Com::CVector2D& eh = (poe-pos)*(1/sqlen);
				const Com::CVector2D ev(-eh.y,eh.x);
        {
          std::vector<double> aRelCo;
          aRelCo.reserve(n*2);
          for(unsigned int i=0;i<n;i++){
            double x1 = Com::Dot(aVecStrok[i]-pos,eh);
            double y1 = Com::Dot(aVecStrok[i]-pos,ev);
            aRelCo.push_back(x1);
            aRelCo.push_back(y1);
          }          
          e.SetCurve_Polyline(aRelCo);
        }
				std::vector<Com::CVector2D> aVecDiv;
				const unsigned int ndiv = e.GetCurveLength()/(camera.GetHalfViewHeight()*0.04)+1;
				e.GetCurveAsPolyline(aVecDiv,ndiv);
        std::vector<double> aRelCo;
				aRelCo.reserve(aVecDiv.size()*2);
				for(unsigned int i=0;i<aVecDiv.size();i++){
					double x1 = Com::Dot(aVecDiv[i]-pos,eh);
					double y1 = Com::Dot(aVecDiv[i]-pos,ev);
					aRelCo.push_back(x1);
					aRelCo.push_back(y1);
				}
        e.SetCurve_Polyline(aRelCo);
			}
			cad_2d.ConnectVertex(e);
      drawerCAD.UpdateCAD_TopologyGeometry(cad_2d);
			aVecStrok.clear();
		}
	}
}

void myGlutMotion( int x, int y )
{
	if( imode == 2 ){
		{ // hilight
			const unsigned int size_buffer = 128;
			unsigned int select_buffer[size_buffer];
			Com::View::PickPre(size_buffer,select_buffer,   x,y,  5,5,  camera);
      drawerCAD.DrawSelection(0);
			std::vector<Com::View::SSelectedObject> aSelecObj = Com::View::PickPost(select_buffer,   x,y,   camera );
      drawerCAD.ClearModeShow(1);
			unsigned int id_part_cad0 = 0;
			Cad::CAD_ELEM_TYPE itype_part_cad0 = Cad::NOT_SET;
      int ictrl_cad0;
			if( aSelecObj.size() != 0 ){
				drawerCAD.GetCadPartID(aSelecObj[0].name,itype_part_cad0,id_part_cad0,ictrl_cad0);
				drawerCAD.SetModeShow(itype_part_cad0,id_part_cad0,1);
			}
			else{
				itype_part_cad0 = Cad::NOT_SET;
				id_part_cad0 = 0;
			}
		}
	}
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;

	if( imodifier == GLUT_ACTIVE_SHIFT ){
		camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	else if( imodifier == 0 ){
		if( imode == 0 ){	// MoveMode
			if( itype_part_cad == Cad::VERTEX && cad_2d.IsElemID(Cad::VERTEX,id_part_cad) ){
				const Com::CVector3D& oloc = camera.ProjectionOnPlane(mov_end_x,mov_end_y);
        cad_2d.MoveVertex(id_part_cad,Com::CVector2D(oloc.x,oloc.y));
				is_updated_cad = true;
			}
			if( itype_part_cad == Cad::EDGE && cad_2d.IsElemID(Cad::EDGE,id_part_cad) ){
				const Com::CVector3D& oloc0 = camera.ProjectionOnPlane(mov_begin_x, mov_begin_y);
				const Com::CVector3D& oloc1 = camera.ProjectionOnPlane(mov_end_x,   mov_end_y);
        if( ictrl_cad == -1 ){
          cad_2d.MoveEdge(id_part_cad, Com::CVector2D(oloc1.x-oloc0.x,oloc1.y-oloc0.y) );
        }
        else{
          cad_2d.MoveEdgeCtrl(id_part_cad,ictrl_cad,Com::CVector2D(oloc1.x-oloc0.x,oloc1.y-oloc0.y) );
        }
				is_updated_cad = true;
			}
			if( itype_part_cad == Cad::LOOP && cad_2d.IsElemID(Cad::LOOP,id_part_cad) ){
				const Com::CVector3D& oloc0 = camera.ProjectionOnPlane(mov_begin_x,mov_begin_y);
				const Com::CVector3D& oloc1 = camera.ProjectionOnPlane(mov_end_x, mov_end_y);
        cad_2d.MoveLoop(id_part_cad,  Com::CVector2D(oloc1.x-oloc0.x,oloc1.y-oloc0.y));
				is_updated_cad = true;
			}
		}
    if( imode == 1 ){
      if( itype_part_cad == Cad::EDGE && cad_2d.IsElemID(itype_part_cad,id_part_cad) ){
        const Com::CVector3D& oloc = camera.ProjectionOnPlane(mov_end_x,mov_end_y);  
        if(      cad_2d.GetEdgeCurveType(id_part_cad) == Cad::CURVE_ARC ){
          cad_2d.DragArc(id_part_cad,Com::CVector2D(oloc.x,oloc.y));
        }
        else if( cad_2d.GetEdgeCurveType(id_part_cad) == Cad::CURVE_POLYLINE ){
          cad_2d.DragPolyline(id_part_cad,Com::CVector2D(oloc.x,oloc.y));
        }
        is_updated_cad = true;
      }
    }
		if( imode == 2 ){	// sketch strok
			const Com::CVector3D& oloc = camera.ProjectionOnPlane(mov_end_x,mov_end_y);
			aVecStrok.push_back( Com::CVector2D(oloc.x,oloc.y) );		
    }
    if( imode == 4 ){ // smooth strok
      if( itype_part_cad == Cad::EDGE && cad_2d.IsElemID(itype_part_cad,id_part_cad) ){
        const Com::CVector3D& oloc = camera.ProjectionOnPlane(mov_end_x,mov_end_y);      
        cad_2d.SmoothingPolylineEdge(id_part_cad, 3,
                                     Com::CVector2D(oloc.x,oloc.y), 0.1); 
        is_updated_cad = true;
      }
    }
	}
	
	////////////////
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);

	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
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
    {
      double rot[9]; camera.RotMatrix33(rot);
      camera.Fit(drawerCAD.GetBoundingBox(rot));
    }        	
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


void myGlutDisplay(void)
{
  ::glClearColor(1.0f,1.0f,1.0f,1.0f);
  ::glClearColor(0.3f,0.8f,0.8f,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
   
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

	if( is_updated_cad ){
		is_updated_cad = false;
		drawerCAD.UpdateCAD_Geometry(cad_2d);
	}
  
  {
    ::glDisable(GL_DEPTH_TEST);
    ::glMatrixMode(GL_PROJECTION);
    ::glPushMatrix();            
    ::glLoadIdentity();    
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();                
    ::glLoadIdentity();            
    char str[256];
    if(      imode == 0 ){ sprintf(str,"drag"); }
    else if( imode == 1 ){ sprintf(str,"edit curve"); }
    else if( imode == 2 ){ sprintf(str,"sketch"); }
    else if( imode == 3 ){ sprintf(str,"delete"); }
    else if( imode == 4 ){ sprintf(str,"smooth"); }     
    void* font = GLUT_BITMAP_HELVETICA_18;    
    RenderBitmapString(-0.95, +0.8, font,str);        
    ::glMatrixMode(GL_PROJECTION);
    ::glPopMatrix();            
    ::glMatrixMode(GL_MODELVIEW);    
    ::glPopMatrix();    
    ::glEnable(GL_LIGHTING);    
    ::glEnable(GL_DEPTH_TEST);    
  }  

	{
		::glBegin(GL_LINE_STRIP);
		for(unsigned int ivec=0;ivec<aVecStrok.size();ivec++){
			::glVertex2d( aVecStrok[ivec].x, aVecStrok[ivec].y );
		}
		::glEnd();
	}

	drawerCAD.Draw();

	::glutSwapBuffers();
}

bool SetNewProblem()
{
	const unsigned int nprob = 9;
	static int iprob = 8;
		
	static int id_l0;
	if( iprob == 0 )
	{
		cad_2d.Clear();
 		{	// Make model
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(2.0,0.0) );
      vec_ary.push_back( Com::CVector2D(2.0,1.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(1.0,2.0) );
      vec_ary.push_back( Com::CVector2D(0.0,2.0) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
	}
	else if( iprob == 1 )
	{
		cad_2d.Clear();
 		{	// Make model
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.5,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		cad_2d.SetCurve_Arc(1,false,-0.2);
		cad_2d.SetCurve_Arc(2,false, 0.5);
		cad_2d.ConnectVertex_Line(1,2);
	}
	else if( iprob == 2 )
	{
		cad_2d.Clear();
 		{	// Make model
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.5,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		cad_2d.SetCurve_Arc(1,false,  0);
		cad_2d.SetCurve_Arc(2,true, -1.0);
		cad_2d.SetCurve_Arc(3,true,  -1.0);
		cad_2d.SetCurve_Arc(5,true,  -0.5);
	}
	else if( iprob == 3 ){
		cad_2d.Clear();
		{	// ê≥ï˚å`Ç…ãÈå`ÇÃåä
			unsigned int id_l = 0;
			{
        std::vector<Com::CVector2D> vec_ary;
        vec_ary.push_back( Com::CVector2D(0.0,0.0) );
        vec_ary.push_back( Com::CVector2D(1.0,0.0) );
        vec_ary.push_back( Com::CVector2D(1.0,1.0) );
        vec_ary.push_back( Com::CVector2D(0.0,1.0) );        
				id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			}
      const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.2)).id_v_add;
      const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.2)).id_v_add;
      const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.8)).id_v_add;
      const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.8)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}
	}
	else if( iprob == 4 ){
		cad_2d.Clear();
		{	// ê≥ï˚å`Ç…ãÈå`ÇÃåä
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );      
			cad_2d.AddPolygon( vec_ary );
			std::vector<Com::CVector2D> aRelCo;
			aRelCo.resize(3);
			{
				aRelCo[0].x = 0.25;	aRelCo[0].y = -0.1;
				aRelCo[1].x = 0.5;	aRelCo[1].y = -0.0;
				aRelCo[2].x= 0.75;	aRelCo[2].y = -0.1;
			}
			cad_2d.SetCurve_Polyline(1,aRelCo);
			cad_2d.SetCurve_Arc(2,true,-1);
		}
	}
	else if( iprob == 5 ){
		cad_2d.Clear();
		unsigned int id_l_l, id_l_r;
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
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_e1 = cad_2d.ConnectVertex_Line(5,9).id_e_add;
			cad_2d.ShiftLayer_Loop(id_l0,true);
			const double col[3] = { 0.9, 0.4, 0.4 };
			cad_2d.SetColor_Loop(id_l0, col);
			cad_2d.AddVertex(Cad::EDGE,3, Com::CVector2D(1.3,0.5) );
			cad_2d.GetIdLoop_Edge(id_l_l, id_l_r, id_e1);
		}
		////////////////
    drawerCAD.UpdateCAD_TopologyGeometry(cad_2d);    
    camera.Fit(drawerCAD.GetBoundingBox());
  }
	else if( iprob == 6 ){
		cad_2d.ShiftLayer_Loop(id_l0,false);
		cad_2d.ShiftLayer_Loop(id_l0,false);
	}
	else if( iprob == 7 )
	{
		cad_2d.Clear();
 		{	// Make model
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(2.0,0.0) );
      vec_ary.push_back( Com::CVector2D(3.0,0.0) );
      vec_ary.push_back( Com::CVector2D(3.0,1.0) );
      vec_ary.push_back( Com::CVector2D(2.0,1.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;      
		}
	}
	else if( iprob == 8 ){
		cad_2d.Clear();
 		{	// Make model
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.5,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		cad_2d.SetCurve_Bezier(1,0.2,+0.5, 0.2,-0.5);
		cad_2d.SetCurve_Bezier(2,0.5,+0.5, 0.2,-0.5);
		cad_2d.ConnectVertex_Line(1,2);
	}
  
  drawerCAD.Clear();
  drawerCAD.SetPointSize(5);
  drawerCAD.UpdateCAD_TopologyGeometry(cad_2d);    
  camera.Fit(drawerCAD.GetBoundingBox());

	::glMatrixMode(GL_PROJECTION);	
	::glLoadIdentity();	
	Com::View::SetProjectionTransform(camera);

	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
	switch (key) {     
    case '0':
      imode = 0;
      break;
    case '1':
      imode = 1;
      break;
    case '2':
      imode = 2;
      break;
    case '3':
      imode = 3;
      break;
    case '4':
      imode = 4;
      break;
    case 'q':
    case 'Q':
    case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
      exit(0);
      break;
    case 's':
      break;
    case ' ':
      SetNewProblem();
      break;
    default:
      break;      
	}
}

void myGlutMenu(int value)
{
  switch(value) {
    case 1:
      imode = 1;
      break;
    case 2:	// ChangeToLine
      if( itype_part_cad != Cad::EDGE ) break;
      std::cout << "ChangeToLine" << std::endl;
      is_updated_cad = cad_2d.SetCurve_Line(id_part_cad);
      break;
    case 3:	// ChangeToArc
      if( itype_part_cad != Cad::EDGE ) break;
      std::cout << "ChangeToArc" << std::endl;
		{
			unsigned int id_vs, id_ve;
			cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_part_cad);
      const Com::CVector2D& vs = cad_2d.GetVertexCoord(id_vs);
      const Com::CVector2D& ve = cad_2d.GetVertexCoord(id_ve);
      const double dist = sqrt( Com::SquareLength(vs,ve) )*10.0;
			is_updated_cad = cad_2d.SetCurve_Arc(id_part_cad,false,dist);
		}
      break;
    case 4:
      if( itype_part_cad != Cad::EDGE ) break;
      std::cout << "ChangeToPolyline" << std::endl;
      cad_2d.SetCurve_Polyline(id_part_cad);   
      break;
  } 
}


int main(int argc,char* argv[])
{

	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 300);
	::glutInit(&argc, argv);
	::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("MSH View");

	::glutMotionFunc(myGlutMotion);
	::glutPassiveMotionFunc(myGlutPassiveMotion);
	::glutMouseFunc(myGlutMouse);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutSpecialFunc(myGlutSpecial);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutIdleFunc(myGlutIdle);

	::glutCreateMenu(myGlutMenu);
  ::glutAddMenuEntry("EditCurve", 1);
  ::glutAddMenuEntry("ChangeToLine", 2);
  ::glutAddMenuEntry("ChangeToArc", 3);
  ::glutAddMenuEntry("ChangeToPolyline", 4);
  ::glutAttachMenu(GLUT_RIGHT_BUTTON);
	
	SetNewProblem();
	
	::glutMainLoop();
	return 0;
}