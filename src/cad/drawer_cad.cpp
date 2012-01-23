/*
DelFEM (Finite Element Analysis)
Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

////////////////////////////////////////////////////////////////
// DrawerCad.cpp : CADモデル描画クラス(CDrawerCad)の実装
// このクラスは何が起きても絶対落ちないように実装すること．
// assertionも原則しないこと
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#if defined(_WIN32)
#include <windows.h>
#if defined(__VISUALC__)
#pragma comment (lib, "winmm.lib")     /* link with Windows MultiMedia lib */
#pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */


#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <assert.h>
#include <iostream>

#include "delfem/drawer_cad.h"
#include "delfem/cad2d_interface.h"
#include "delfem/vector2d.h"
#include "delfem/mesher2d.h"

using namespace Cad::View;
using namespace Com;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void CDrawerRubberBand::Draw() const
{
	::glLineWidth(2);
	::glColor3d(0.0,0.0,0.0);
	::glBegin(GL_LINES);
	if( imode == 0 ){
		::glVertex3d(initial.x, initial.y, initial.z);
		::glVertex3d(mouse.x, mouse.y, mouse.z);
	}
	else if( imode == 1 ){
		for(unsigned int iline=0;iline<fix.size();iline++){
			::glVertex3d(fix[iline].x, fix[iline].y, fix[iline].z);
			::glVertex3d(mouse.x, mouse.y, mouse.z);
		}
	}
	else if( imode == 2 ){
		assert( fix.size() == 2 );
		::glVertex3d(
			fix[0].x + mouse.x - initial.x,  
			fix[0].y + mouse.y - initial.y,  
			fix[0].z + mouse.z - initial.z );
		::glVertex3d(
			fix[1].x + mouse.x - initial.x,
			fix[1].y + mouse.y - initial.y, 
			fix[1].z + mouse.z - initial.z );
		for(unsigned int iline=0;iline<fix_stat.size();iline++){
			::glVertex3d(
				fix_stat[iline].second.x, 
				fix_stat[iline].second.y, 
				fix_stat[iline].second.z);
			unsigned int ifix = fix_stat[iline].first;
			::glVertex3d(
				fix[ifix].x + mouse.x - initial.x, 
				fix[ifix].y + mouse.y - initial.y, 
				fix[ifix].z + mouse.z - initial.z );
		}
	}
	::glEnd();
}

Com::CBoundingBox3D CDrawerRubberBand::GetBoundingBox( double rot[] ) const
{
	////////////////
	double x_max,x_min,  y_max,y_min,  z_max,z_min;
	{
		const double x1 = mouse.x; const double y1 = mouse.y; const double z1 = mouse.z;
		const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
		const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
		const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
		x_max = x2;  x_min = x2;
		y_max = y2;  y_min = y2;
		z_max = z2;  z_min = z2;
	}
	for(unsigned int iline=0;iline<fix.size();iline++){
		const double x1 = fix[iline].x; const double y1 = fix[iline].y; const double z1 = fix[iline].z;
		const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
		const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
		const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
		x_max = ( x2 > x_max ) ? x2 : x_max;  x_min = ( x2 < x_min ) ? x2 : x_min;			
		y_max = ( y2 > y_max ) ? y2 : y_max;  y_min = ( y2 < y_min ) ? y2 : y_min;			
		z_max = ( z2 > z_max ) ? z2 : z_max;  z_min = ( z2 < z_min ) ? z2 : z_min;			
	}
	double x_size = x_max - x_min;
	double y_size = y_max - y_min;
	double z_size = z_max - z_min;
	double x_cent1 = ( x_max + x_min )*0.5;
	double y_cent1 = ( y_max + y_min )*0.5;
	double z_cent1 = ( z_max + z_min )*0.5;
//	double x_cent = x_cent1*rot[0]+y_cent1*rot[3]+z_cent1*rot[6];
//	double y_cent = x_cent1*rot[1]+y_cent1*rot[4]+z_cent1*rot[7];
//  double z_cent = x_cent1*rot[2]+y_cent1*rot[5]+z_cent1*rot[8];
	return CBoundingBox3D(
		x_cent1-x_size*0.5,x_cent1+x_size*0.5,
		y_cent1-y_size*0.5,y_cent1+y_size*0.5,
		z_cent1-z_size*0.5,z_cent1+z_size*0.5 );
}

CDrawerRubberBand::CDrawerRubberBand(const Cad::CCadObj2D& cad, unsigned int id_v)
{
	this->sutable_rot_mode = 1;	
  Cad::CBRepSurface::CItrVertex itrv = cad.GetItrVertex(id_v);
	for(;!itrv.IsEnd();itrv++){
		unsigned int id_e; bool is_same_dir;
		itrv.GetIdEdge_Behind(id_e,is_same_dir);
		unsigned int id_vs = cad.GetIdVertex_Edge(id_e,true );
    unsigned int id_ve = cad.GetIdVertex_Edge(id_e,false);
    const unsigned int id_v1 = (is_same_dir) ? id_ve : id_vs;
		assert( ((is_same_dir) ? id_vs : id_ve) == id_v );
    assert( cad.IsElemID(Cad::VERTEX,id_v1) );
		const Com::CVector2D& vec2d = cad.GetVertexCoord(id_v1);
		fix.push_back( CVector3D(vec2d.x, vec2d.y, 0.0) );
	}
	imode = 1;
}

CDrawerRubberBand::CDrawerRubberBand(const Cad::CCadObj2D& cad,
                                     unsigned int id_e, const Com::CVector3D& initial)
{
	this->sutable_rot_mode = 1;	
	if( !cad.IsElemID(Cad::EDGE,id_e) ) return;
	this->initial = initial;
	const unsigned int id_vs = cad.GetIdVertex_Edge(id_e,true );
  const unsigned int id_ve = cad.GetIdVertex_Edge(id_e,false);
	fix.clear();
	{
		{
      assert( cad.IsElemID(Cad::VERTEX,id_vs) );
			const Com::CVector2D& vec2d = cad.GetVertexCoord(id_vs);
			fix.push_back( CVector3D(vec2d.x,vec2d.y,0) );
		}
		{
      Cad::CBRepSurface::CItrVertex itrv = cad.GetItrVertex(id_vs);
			for(;!itrv.IsEnd();itrv++){
				unsigned int id_e0; bool is_same_dir0;
				itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
				if( id_e0 == id_e ) continue;
				unsigned int id_vs0 = cad.GetIdVertex_Edge(id_e0,true);
        unsigned int id_ve0 = cad.GetIdVertex_Edge(id_e0,false);
        const unsigned int id_v1 = (is_same_dir0) ? id_ve0 : id_vs0;
				assert( ((is_same_dir0) ? id_vs0 : id_ve0) == id_vs );
				assert( cad.IsElemID(Cad::VERTEX,id_v1) );
				const Com::CVector2D& vec2d = cad.GetVertexCoord(id_v1);
				fix_stat.push_back( std::make_pair(0,Com::CVector3D(vec2d.x, vec2d.y, 0.0)) );
			}
		}
		{
      assert( cad.IsElemID(Cad::VERTEX,id_ve) );
			const Com::CVector2D& vec2d = cad.GetVertexCoord(id_ve);
			fix.push_back( CVector3D(vec2d.x,vec2d.y,0) );
		}
		{
      Cad::CBRepSurface::CItrVertex itrv = cad.GetItrVertex(id_ve);
			for(;!itrv.IsEnd();itrv++){
				unsigned int id_e0; bool is_same_dir0;
				itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
				if( id_e0 == id_e ) continue;
				unsigned int id_vs0 =	cad.GetIdVertex_Edge(id_e,true );
        unsigned int id_ve0 = cad.GetIdVertex_Edge(id_e,false);
				const unsigned int id_v1 = (is_same_dir0) ? id_ve0 : id_vs0;
				assert( ((is_same_dir0) ? id_vs0 : id_ve0) == id_ve );
				assert( cad.IsElemID(Cad::VERTEX,id_v1) );
				const Com::CVector2D& vec2d = cad.GetVertexCoord(id_v1);
				fix_stat.push_back( std::make_pair(1,Com::CVector3D(vec2d.x, vec2d.y, 0.0)) );
			}
		}
	}
	imode = 2;
}


////////////////////////////////////////////////////////////////

bool CDrawer_Cad2D::CDrawPart::SetTriAry(const Msh::CTriAry2D& TriAry)
{
	this->id_cad = TriAry.id_l_cad;	assert( id_cad != 0 );
	this->id_msh = TriAry.id;		assert( id_msh != 0 );
	this->itype  = Cad::LOOP;
	////////////////
	npoel = 3;
	nelem = TriAry.m_aTri.size();
	if( pIndexArray != 0 ){ delete[] pIndexArray; }
	pIndexArray = new unsigned int [nelem*npoel];
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		for(unsigned int ipoel=0;ipoel<npoel;ipoel++){
			pIndexArray[ielem*npoel+ipoel] = TriAry.m_aTri[ielem].v[ipoel];
		}
	}
  ////////////////
  color[0]=0.8;   color[1]=0.8;   color[2]=0.8;
	return true;
}

bool CDrawer_Cad2D::CDrawPart::SetBarAry(const Msh::CBarAry& BarAry)
{
	this->id_cad = BarAry.id_e_cad;	assert( id_cad != 0 );
	this->id_msh = BarAry.id;		assert( id_msh != 0 );
	this->itype  = Cad::EDGE;
	////////////////
	npoel = 2;
	nelem = BarAry.m_aBar.size();
	if( pIndexArray != 0 ){ delete[] pIndexArray; }
	pIndexArray = new unsigned int [nelem*npoel];
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		for(unsigned int ipoel=0;ipoel<npoel;ipoel++){
			pIndexArray[ielem*npoel+ipoel] = BarAry.m_aBar[ielem].v[ipoel];
		}
	}
  ////////////////
  color[0]=0.0;   color[1]=0.0;   color[2]=0.0;
	return true;
}


bool CDrawer_Cad2D::CDrawPart::SetVertex(const Msh::SVertex& vtx)
{
	this->id_cad = vtx.id_v_cad;	assert( id_cad != 0 );
	this->id_msh = vtx.id;		assert( id_msh != 0 );
	this->itype  = Cad::VERTEX;
	////////////////
	npoel = 1;
	nelem = 1;
	if( pIndexArray != 0 ){ delete[] pIndexArray; }
	pIndexArray = new unsigned int [nelem*npoel];
  pIndexArray[0] = vtx.v;
  ////////////////
  color[0]=0.0;   color[1]=0.0;   color[2]=0.0;
	return true;
}

void CDrawer_Cad2D::CDrawPart::DrawElements() const
{
  if( npoel == 1 ){
    ::glDrawElements(GL_POINTS,nelem*npoel,GL_UNSIGNED_INT,pIndexArray); 
		return;
  }
	if( npoel == 2 ){ 
		::glDrawElements(GL_LINES    ,nelem*npoel,GL_UNSIGNED_INT,pIndexArray); 
		return;
	}
	if( npoel == 3 ){ 
		::glDrawElements(GL_TRIANGLES,nelem*npoel,GL_UNSIGNED_INT,pIndexArray); 
		return;
	}
}

////////////////////////////////////////////////////////////////

bool CDrawer_Cad2D::UpdateCAD_TopologyGeometry(const Cad::CCadObj2D &cad_2d)
{
	this->sutable_rot_mode = 1;	
  //! 今までのDrawerPartの配列を一端バッファにコピーして，必要な物だけを戻す
  std::vector<CDrawPart*> apDrawPart_old = m_apDrawPart;
  for(unsigned int idp=0;idp<apDrawPart_old.size();idp++){
    apDrawPart_old[idp]->id_msh = 0;
    apDrawPart_old[idp]->imode_show = 0;
  }
  m_apDrawPart.clear();
  
	int ilayer_min, ilayer_max;
	cad_2d.GetLayerMinMax(ilayer_min, ilayer_max);
	double layer_height = 1.0/(ilayer_max-ilayer_min+1); 

	{	// 面をセット
    const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(LOOP);
    for(unsigned int iil=0;iil<aIdL.size();iil++){
      const unsigned int id_l = aIdL[iil];
			double height = 0;
			{
				unsigned int ilayer = cad_2d.GetLayer(Cad::LOOP,id_l);
				height = (ilayer-ilayer_min)*layer_height;
			}
      unsigned int idp0 = 0;
      for(;idp0<apDrawPart_old.size();idp0++){
        CDrawPart& dpo = *apDrawPart_old[idp0];        
        if( dpo.itype == Cad::LOOP && dpo.id_cad == id_l )
        {          
          dpo.id_msh = 1;          
          dpo.SetHeight(height);
          double cd[3];
          cad_2d.GetColor_Loop(id_l,cd);
          dpo.color[0] = (float)cd[0];
          dpo.color[1] = (float)cd[1];
          dpo.color[2] = (float)cd[2];          
          this->m_apDrawPart.push_back( apDrawPart_old[idp0] );
          break;
        }
      }
      if( idp0 == apDrawPart_old.size() ){
        CDrawPart* dp = new CDrawPart;
        dp->id_cad = id_l;
        dp->itype = LOOP;        
        dp->SetHeight( height );
        double cd[3];
        cad_2d.GetColor_Loop(id_l,cd);
        dp->color[0] = (float)cd[0];
        dp->color[1] = (float)cd[1];
        dp->color[2] = (float)cd[2];              
        this->m_apDrawPart.push_back( dp );
      }      
		}
	}
  
	{	// set edge
    const std::vector<unsigned int>& aIdE = cad_2d.GetAryElemID(EDGE);
    for(unsigned int iie=0;iie<aIdE.size();iie++){
      const unsigned int id_e = aIdE[iie];
			double height = 0;
			{
				int ilayer = cad_2d.GetLayer(Cad::EDGE,id_e);
				height += (ilayer-ilayer_min+0.01)*layer_height;
			}
      unsigned int idp0 = 0;
      for(;idp0<apDrawPart_old.size();idp0++){
        CDrawPart& dpo = *apDrawPart_old[idp0];
        if( dpo.itype == Cad::EDGE && dpo.id_cad == id_e ){
          dpo.id_msh = 1;          
          dpo.SetHeight( height );
          this->m_apDrawPart.push_back( &dpo );
          break;
        }
      }
      if( idp0 == apDrawPart_old.size() ){
        CDrawPart* dp = new CDrawPart;
        dp->id_cad = id_e;
        dp->itype = EDGE;
        dp->SetHeight( height );
        this->m_apDrawPart.push_back( dp );
      }
      {
        CDrawPart& dp = *m_apDrawPart[m_apDrawPart.size()-1];      
        const CEdge2D& edge = cad_2d.GetEdge(id_e);
        dp.aCtrlPoint.clear();      
        dp.itype_curve = edge.GetCurveType();
        if( edge.GetCurveType() == CURVE_ARC ){
          Com::CVector2D po_c; double radius;
          edge.GetCenterRadius(po_c,radius);
          dp.aCtrlPoint.push_back(po_c);
        }
        else if( edge.GetCurveType() == CURVE_BEZIER ){
          const std::vector<Com::CVector2D>& aCo = edge.GetCurvePoint();
          dp.aCtrlPoint.push_back( aCo[0] );
          dp.aCtrlPoint.push_back( aCo[1] );
        }        
      }      
		}
	}
  
	{	// set vertex
    const std::vector<unsigned int>& aIdV = cad_2d.GetAryElemID(VERTEX);    
		for(unsigned int iiv=0;iiv<aIdV.size();iiv++){
			const unsigned int id_v_cad = aIdV[iiv];
			int ilayer = cad_2d.GetLayer(Cad::VERTEX,id_v_cad);      
			const double height = (ilayer-ilayer_min+0.1)*layer_height;
      unsigned int idp0 = 0;
      for(;idp0<apDrawPart_old.size();idp0++){
        CDrawPart& dpo = *apDrawPart_old[idp0];
        if( dpo.itype == Cad::VERTEX && dpo.id_cad == id_v_cad ){
          dpo.id_msh = 1;
          dpo.SetHeight( height );
          this->m_apDrawPart.push_back( &dpo );
          break;
        }
      }
      if( idp0 == apDrawPart_old.size() ){
        CDrawPart* dp = new CDrawPart;
        dp->id_cad = id_v_cad;
        dp->itype = VERTEX;
        dp->SetHeight( height );
        this->m_apDrawPart.push_back( dp );
      }      
		}
	}
  
  for(unsigned int idp=0;idp<apDrawPart_old.size();idp++){
    if( apDrawPart_old[idp]->id_msh == 0 ){ delete apDrawPart_old[idp]; }
  }
  apDrawPart_old.clear();  
  
  UpdateCAD_Geometry(cad_2d);
	return true;
}

void CDrawer_Cad2D::UpdateCAD_Geometry(const Cad::CCadObj2D& cad_2d)
{
	Msh::CMesher2D mesh(cad_2d);
	for(unsigned int idp=0;idp<m_apDrawPart.size();idp++){
    CDrawPart& dp = *m_apDrawPart[idp];
		dp.Clear();
		const unsigned int id_cad = dp.id_cad;
		Cad::CAD_ELEM_TYPE itype_cad = dp.itype;
		if( !cad_2d.IsElemID(itype_cad,id_cad) ){ continue; }
		const unsigned int id_msh = mesh.GetElemID_FromCadID(id_cad,itype_cad);
		if( id_msh == 0 ) continue;
		Msh::MSH_TYPE msh_type;
		unsigned int nelem, iloc, id_cad0;
		mesh.GetMshInfo(id_msh,  nelem,msh_type,iloc,id_cad0);
		assert( id_cad0 == id_cad );
		if( msh_type == Msh::TRI ){
			dp.SetTriAry( mesh.GetTriArySet()[iloc] );
      double cd[3];
      cad_2d.GetColor_Loop(id_cad0,cd);
      dp.color[0] = (float)cd[0];
      dp.color[1] = (float)cd[1];
      dp.color[2] = (float)cd[2];
    }
		else if( msh_type == Msh::BAR ){
			dp.SetBarAry( mesh.GetBarArySet()[iloc] );
      assert( itype_cad == EDGE );
      const CEdge2D& edge = cad_2d.GetEdge(id_cad);      
      dp.itype_curve = edge.GetCurveType();      
      dp.aCtrlPoint.clear();      
      if( edge.GetCurveType() == CURVE_ARC ){
        Com::CVector2D po_c; double radius;
        edge.GetCenterRadius(po_c,radius);
        dp.aCtrlPoint.push_back(po_c);
      }
      else if( edge.GetCurveType() == CURVE_BEZIER ){
        const std::vector<Com::CVector2D>& aCo = edge.GetCurvePoint();
        dp.aCtrlPoint.push_back( aCo[0] );
        dp.aCtrlPoint.push_back( aCo[1] );
      }
		}
    else if( msh_type == Msh::VERTEX ){
			dp.SetVertex( mesh.GetVertexAry()[iloc] );
    }
	}
	
	{	// 座標をセット
		const std::vector<CVector2D>& aVec2D = mesh.GetVectorAry();
		const unsigned int npoin = aVec2D.size();
		const unsigned int ndim = 2;
		m_vertex_ary.SetSize(npoin,ndim);
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			m_vertex_ary.pVertexArray[ipoin*ndim  ] = aVec2D[ipoin].x;
			m_vertex_ary.pVertexArray[ipoin*ndim+1] = aVec2D[ipoin].y;
		}
    if( m_vertex_ary.pUVArray != 0 ){      
      for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
        m_vertex_ary.pUVArray[ipoin*ndim  ] = aVec2D[ipoin].x*tex_scale;
        m_vertex_ary.pUVArray[ipoin*ndim+1] = aVec2D[ipoin].y*tex_scale;
      }      
    }    
	}
}

void CDrawer_Cad2D::GetCadPartID
(const int selec_flag[],
 Cad::CAD_ELEM_TYPE& part_type, unsigned int& part_id, int& ictrl)
{
  unsigned int idp = selec_flag[1];
  if( idp < m_apDrawPart.size() ){
    CDrawPart& dp = *m_apDrawPart[idp];
    part_type = dp.itype;
    part_id = dp.id_cad;
    ictrl = selec_flag[2];
    return;
  }
	part_type = Cad::NOT_SET;
	part_id = 0;
}

void CDrawer_Cad2D::SetModeShow(Cad::CAD_ELEM_TYPE itype, unsigned int id, int imode_show)
{
	bool iflag = false;
  for(unsigned int iea=0;iea<m_apDrawPart.size();iea++){
    if( m_apDrawPart[iea]->itype == itype && m_apDrawPart[iea]->id_cad == id ){
      m_apDrawPart[iea]->imode_show = imode_show;
      assert( iflag == false );
      iflag = true;
    }
  }
}

int CDrawer_Cad2D::GetModeShow(Cad::CAD_ELEM_TYPE itype, unsigned int id)
{
  for(unsigned int iea=0;iea<m_apDrawPart.size();iea++){
    if( m_apDrawPart[iea]->itype == itype && m_apDrawPart[iea]->id_cad == id ){
      return m_apDrawPart[iea]->imode_show;
    }
  }
}

void CDrawer_Cad2D::ClearModeShow(int imode_show)
{
	for(unsigned int iea=0;iea<m_apDrawPart.size();iea++){
    if( m_apDrawPart[iea]->imode_show == imode_show ){ 
      m_apDrawPart[iea]->imode_show = 0;
    }
	}
}
/*
void CDrawer_Cad2D::Hide(Cad::CAD_ELEM_TYPE part_type, unsigned int part_id)
{
    if( part_type == Cad::EDGE || part_type == Cad::LOOP ){
		for(unsigned int iea=0;iea<m_apDrawPart.size();iea++){
            if(    m_apDrawPart[iea]->id_cad == part_id
                && m_apDrawPart[iea]->itype  == part_type){
				m_apDrawPart[iea]->is_show = false;
			}
		}
	}
    else if( part_type == Cad::VERTEX ){
		for(unsigned int iv=0;iv<this->m_aIndexVertex.size();iv++){
            if( m_aIndexVertex[iv].id_cad == part_id ){
				m_aIndexVertex[iv].is_show = false;
			}
		}
	}
}
*/
void CDrawer_Cad2D::SetIsShow(bool is_show, Cad::CAD_ELEM_TYPE itype_part_cad, unsigned int id_part_cad)
{
  int imode = (is_show) ? 0 : -1;
//  if( itype_part_cad == Cad::EDGE || itype_part_cad == Cad::LOOP ){
		for(unsigned int iea=0;iea<m_apDrawPart.size();iea++){
      if(    m_apDrawPart[iea]->id_cad == id_part_cad
         && m_apDrawPart[iea]->itype  == itype_part_cad ){
				m_apDrawPart[iea]->imode_show = imode;
			}
		}
  /*
	}
  else if( itype_part_cad == Cad::VERTEX ){
		for(unsigned int iv=0;iv<this->m_aIndexVertex.size();iv++){
			if( m_aIndexVertex[iv].id_cad == id_part_cad ){
				m_aIndexVertex[iv].imode_show = imode;
			}
		}
	}
   */
}

void CDrawer_Cad2D::SetIsShow
(bool is_show, 
 Cad::CAD_ELEM_TYPE itype_part_cad, const std::vector<unsigned int>& aIdPart )
{
	for(unsigned int i=0;i<aIdPart.size();i++){
		unsigned int id = aIdPart[i];
		this->SetIsShow(is_show,itype_part_cad,id);
	}
}

void CDrawer_Cad2D::SetRigidDisp
(unsigned int id_l, double xdisp, double ydisp)
{	
	for(unsigned int iea=0;iea<m_apDrawPart.size();iea++){
    if(    m_apDrawPart[iea]->id_cad == id_l
       && m_apDrawPart[iea]->itype  == Cad::LOOP ){
			m_apDrawPart[iea]->xdisp = xdisp;
			m_apDrawPart[iea]->ydisp = ydisp;
		}
	}
}

void CDrawer_Cad2D::HideEffected(const Cad::CCadObj2D& cad_2d,
                                 Cad::CAD_ELEM_TYPE part_type, unsigned int part_id)
{
  if(      part_type == Cad::VERTEX ){
		const unsigned int id_v = part_id;
    for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_v);!itrv.IsEnd();itrv++){
			unsigned int id_e0; bool is_same_dir0;
			itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
			this->SetIsShow(false,Cad::EDGE,id_e0);
			const unsigned int id_l = itrv.GetIdLoop();
			this->SetIsShow(false,Cad::LOOP,id_l);
		}
	}
  else if( part_type == Cad::EDGE ){
		if( !cad_2d.IsElemID(Cad::EDGE,part_id) ) return;
		unsigned int id_vs = cad_2d.GetIdVertex_Edge(part_id,true );		
    unsigned int id_ve = cad_2d.GetIdVertex_Edge(part_id,false);		
    for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_vs);!itrv.IsEnd();itrv++){
			unsigned int id_e0; bool is_same_dir0;
			itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
			this->SetIsShow(false,Cad::EDGE,id_e0);
			const unsigned int id_l = itrv.GetIdLoop();
			this->SetIsShow(false,Cad::LOOP,id_l);
		}		
		for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_ve);!itrv.IsEnd();itrv++){
			unsigned int id_e0; bool is_same_dir0;
			itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
			this->SetIsShow(false,Cad::EDGE,id_e0);
			const unsigned int id_l = itrv.GetIdLoop();
			this->SetIsShow(false,Cad::LOOP,id_l);
		}
	}
}

void CDrawer_Cad2D::ShowEffected(const Cad::CCadObj2D& cad_2d,
                                 Cad::CAD_ELEM_TYPE part_type, unsigned int part_id)
{
	std::vector<unsigned int> aEdgeID, aLoopID;
	if( part_type == 1 ){
		const unsigned int id_v = part_id;
		for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_v);!itrv.IsEnd();itrv++){
			unsigned int id_e0; bool is_same_dir0;
			itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
			this->SetIsShow(true,Cad::EDGE,id_e0);
			const unsigned int id_l = itrv.GetIdLoop();
			this->SetIsShow(true,Cad::LOOP,id_l);
		}
	}
}

void CDrawer_Cad2D::EnableUVMap(bool is_uv_map){
  m_vertex_ary.EnableUVMap(is_uv_map);
  const unsigned int npoin = m_vertex_ary.NPoin();
  const unsigned int ndim = 2;
  if( m_vertex_ary.pUVArray != 0 ){      
    for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
      m_vertex_ary.pUVArray[ipoin*ndim  ] = m_vertex_ary.pVertexArray[ipoin*ndim  ]*3;
      m_vertex_ary.pUVArray[ipoin*ndim+1] = m_vertex_ary.pVertexArray[ipoin*ndim+1]*3;
    }      
  }  
}

void CDrawer_Cad2D::SetTextureScale(double tex_scale)
{
  this->tex_scale = tex_scale;
  const unsigned int npoin = m_vertex_ary.NPoin();
  const unsigned int ndim = 2;
  if( m_vertex_ary.pUVArray != 0 ){      
    for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
      m_vertex_ary.pUVArray[ipoin*ndim  ] = m_vertex_ary.pVertexArray[ipoin*ndim  ]*tex_scale;
      m_vertex_ary.pUVArray[ipoin*ndim+1] = m_vertex_ary.pVertexArray[ipoin*ndim+1]*tex_scale;
    }      
  }      
}

void myGlVertex2d(const Com::CVector2D& p){
  glVertex2d(p.x,p.y);
}

void CDrawer_Cad2D::Draw() const
{
	::glEnable(GL_DEPTH_TEST);
  ::glDisable(GL_CULL_FACE);    // to make the program simple...
	const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  const bool is_texture  = ::glIsEnabled(GL_TEXTURE_2D);  
  const bool is_blend    = ::glIsEnabled(GL_BLEND);
  ::glDisable(GL_LIGHTING);

	const unsigned int ndim = this->m_vertex_ary.NDim();

	////////////////////////////////////////////////////////////////
	// モデルの描画
  /*
  { // draw vertecies
    ::glDisable(GL_TEXTURE_2D);
    ////////////////
    ::glPointSize(m_pointsize);
    ::glBegin(GL_POINTS);
    for(unsigned int iver=0;iver<this->m_aIndexVertex.size();iver++){
      if( this->m_aIndexVertex[iver].imode_show == -1 ) continue;
      const double height = this->m_aIndexVertex[iver].height;
      if(      this->m_aIndexVertex[iver].imode_show == 1 ){ ::glColor3d(1.0,1.0,0.0); }
      else if( this->m_aIndexVertex[iver].imode_show == 2 ){ ::glColor3d(1.0,0.0,0.0); }
      
      else{ ::glColor3d(0.0,0.0,0.0);	}
      unsigned int ipo0 = this->m_aIndexVertex[iver].id_v;
      ::glVertex3d(m_vertex_ary.pVertexArray[ipo0*ndim+0], 
                   m_vertex_ary.pVertexArray[ipo0*ndim+1], 
                   height );
    }
    ::glEnd();      
    if( is_texture ){ glEnable(GL_TEXTURE_2D); }
  } 
   */
  ////
  // regist vertex array
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);
  if( is_texture && m_vertex_ary.pUVArray!=0 ){
    ::glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    ::glTexCoordPointer(2,GL_DOUBLE,0,m_vertex_ary.pUVArray);  
    ::glMatrixMode(GL_TEXTURE);
    ::glLoadIdentity();
    ::glTranslated(-tex_cent_x, -tex_cent_y, 0.0);          
  }
  ::glPointSize(m_pointsize);
  ::glLineWidth(m_linewidth);
	for(unsigned int idp=0;idp<m_apDrawPart.size();idp++){
    const CDrawPart& dp = *m_apDrawPart[idp];
    if( dp.imode_show == -1 ) continue;
		const double height = dp.height;
		const double xdisp = dp.xdisp;
		const double ydisp = dp.ydisp;
    if( dp.itype == Cad::VERTEX ){
      ::glDisable(GL_TEXTURE_2D);
      const double height = dp.height;
      if(      dp.imode_show == 1 ){ ::glColor3d(1.0,1.0,0.0); }
      else if( dp.imode_show == 2 ){ ::glColor3d(1.0,0.0,0.0); }
      else{ ::glColor3d(0.0,0.0,0.0);	}
			::glTranslated(0.0,0.0, height);            
      dp.DrawElements();
			::glTranslated(0.0,0.0,-height);            
      if( is_texture ){ glEnable(GL_TEXTURE_2D); }
    }
    if( dp.itype == Cad::EDGE )	// draw edge
    {
      ::glDisable(GL_TEXTURE_2D);
      ::glLineWidth(m_linewidth);
      if( this->m_is_anti_aliasing ){ // anti aliasing
        ::glEnable(GL_LINE_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        ::glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
      }
			if(      dp.imode_show == 1 ){ ::glColor3d(1.0,1.0,0.0);  }
			else if( dp.imode_show == 2 ){ ::glColor3d(1.0,0.0,0.0);  }
      else{ ::glColor3d(0,0,0); }
			::glTranslated(0.0,0.0, height);
			dp.DrawElements();
      if( dp.imode_show > 0 ){ // draw ctrl point        
        ::glBegin(GL_POINTS);
        for(unsigned int icp=0;icp<dp.aCtrlPoint.size();icp++){
          const Com::CVector2D& cp = dp.aCtrlPoint[icp];
          ::glVertex3d(cp.x,cp.y,0.0);
        }    
        ::glEnd();        
        // draw line between ctrl point and point
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1 , 0xF0F0);
        ::glPolygonStipple((const GLubyte*)m_mask);
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        unsigned int ivs = dp.pIndexArray[0];        
        unsigned int ive = dp.pIndexArray[dp.nelem*dp.npoel-1];
        const double* pva = m_vertex_ary.pVertexArray;
        Com::CVector2D ps(pva[ivs*2+0],pva[ivs*2+1]);
        Com::CVector2D pe(pva[ive*2+0],pva[ive*2+1]);        
        if( dp.itype_curve == CURVE_ARC ){
          ::myGlVertex2d(ps);
          ::myGlVertex2d(dp.aCtrlPoint[0]);
          ::myGlVertex2d(pe);    
          ::myGlVertex2d(dp.aCtrlPoint[0]);          
        }
        if( dp.itype_curve == CURVE_BEZIER ){
          ::myGlVertex2d(ps);
          ::myGlVertex2d(dp.aCtrlPoint[0]);
          ::myGlVertex2d(pe);    
          ::myGlVertex2d(dp.aCtrlPoint[1]);                    
        }
        ::glEnd();
        ::glDisable(GL_LINE_STIPPLE);              
      }
			::glTranslated(0.0,0.0,-height);
      ::glDisable(GL_LINE_SMOOTH);
      ::glDisable(GL_BLEND);
      if( is_texture ){ glEnable(GL_TEXTURE_2D); }      
		}
    else if( dp.itype == Cad::LOOP ) // draw loop
    {
      ::glDisable(GL_BLEND);
      if( dp.imode_show > 0 ){
				::glEnable(GL_POLYGON_STIPPLE);
				::glPolygonStipple((const GLubyte*)m_mask);
        if( dp.imode_show == 1 ){ ::glColor3d(1.0,1.0,0.0); }
        if( dp.imode_show == 2 ){ ::glColor3d(1.0,0.0,0.0); }
        ::glTranslated(0.0,0.0,+height+0.001);
				dp.DrawElements();
        ::glTranslated(0.0,0.0,-height-0.001);
				::glDisable(GL_POLYGON_STIPPLE);
			}
      if( dp.imode_show ) continue;
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, dp.color); 
      ::glColor3fv(dp.color);
			::glTranslated(+xdisp,+ydisp,+height);
			dp.DrawElements();
			::glTranslated(-xdisp,-ydisp,-height);
		}
	}
  ::glDisableClientState(GL_VERTEX_ARRAY);
  ::glDisableClientState(GL_TEXTURE_COORD_ARRAY);      
  if( is_lighting ){ ::glEnable(GL_LIGHTING);   } else{ ::glDisable(GL_LIGHTING);    }
  if( is_blend    ){ ::glEnable(GL_BLEND);      } else{ ::glDisable(GL_BLEND);      }
  if( is_texture  ){ ::glEnable(GL_TEXTURE_2D); } else{ ::glDisable(GL_TEXTURE_2D); }
	return;
}

void CDrawer_Cad2D::DrawSelection(unsigned int idraw) const
{
  ////////////////
  const bool is_blend       = ::glIsEnabled(GL_BLEND);
  const bool is_line_smooth = ::glIsEnabled(GL_LINE_SMOOTH);
  const bool is_texture     = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_BLEND);
  ::glDisable(GL_LINE_SMOOTH);
  ::glDisable(GL_TEXTURE_2D);
	const unsigned int ndim = this->m_vertex_ary.NDim();
	::glPushName(idraw);
	// モデルの描画
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);
	for(unsigned int idp=0;idp<m_apDrawPart.size();idp++){
		const CDrawPart& dp = *m_apDrawPart[idp];
		const double height = dp.height;
		::glPushName(idp);
    ::glTranslated(0.0,0.0,+height);
    dp.DrawElements();
    if( dp.itype == EDGE && dp.imode_show == 2){
      for(unsigned int icp=0;icp<dp.aCtrlPoint.size();icp++){
        const Com::CVector2D& cp = dp.aCtrlPoint[icp];
        ::glPushName(icp);
        ::glBegin(GL_POINTS);        
        ::myGlVertex2d(cp);
        ::glEnd();                      
        ::glPopName();
      }    
    }
    ::glTranslated(0.0,0.0,-height);
		::glPopName();
	}
	::glDisableClientState(GL_VERTEX_ARRAY);	
	::glPopName();
  
  if( is_blend       ){ ::glEnable(GL_BLEND);       }
  if( is_line_smooth ){ ::glEnable(GL_LINE_SMOOTH); }
  if( is_texture     ){ ::glEnable(GL_TEXTURE_2D);  }

	return;
}
