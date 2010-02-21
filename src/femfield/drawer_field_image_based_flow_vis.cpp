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
// DrawerField.cpp : 場可視化クラス(DrawerField)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
#endif

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

#include <assert.h>
#include <iostream>
#include <vector>
#include <stdio.h>

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include "delfem/drawer_field_image_based_flow_vis.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

// 三角形の面積を求める関数
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}

bool CEdgeTextureColor::Update(const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(this->id_field_velo) ){ return false; }
	if( !world.IsIdField(this->id_part_field_velo) ){ return false; }
	if( !world.IsIdEA(this->id_ea) ){ return false; }
	const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == Fem::Field::LINE );
	if( nelem != ea.Size() ){
		nelem  = ea.Size();
		if( m_aXYVeloElem != 0 ){ delete[] m_aXYVeloElem; }
		if( m_aXYElem     != 0 ){ delete[] m_aXYElem; }
		m_aXYVeloElem = new double [nelem*4];
		m_aXYElem     = new double [nelem*4];
	}
	const Fem::Field::CField& fv = world.GetField(this->id_part_field_velo);
//	const Fem::Field::CField::CNodeSegInNodeAry& nans = fv.GetNodeSegInNodeAry(Fem::Field::CORNER);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea,CORNER,true, world);
	const Fem::Field::CElemAry::CElemSeg& es_c = fv.GetElemSeg(id_ea,CORNER,false,world);
	assert( es_v.GetSizeNoes() == 2 );
	assert( es_c.GetSizeNoes() == 2 );
	for(unsigned int ielem=0;ielem<nelem;ielem++)
	{
		unsigned int noes[2];
		es_c.GetNodes(ielem,noes);
		double co[2][2];
		ns_c.GetValue(noes[0],co[0]);
		ns_c.GetValue(noes[1],co[1]);
		////////////////
		es_v.GetNodes(ielem,noes);
		double velo[2][2];
		ns_v.GetValue(noes[0],velo[0]);	
		ns_v.GetValue(noes[1],velo[1]);
		////////////////
		m_aXYElem[ielem*4+0] = co[0][0];
		m_aXYElem[ielem*4+1] = co[0][1];
		m_aXYElem[ielem*4+2] = co[1][0];
		m_aXYElem[ielem*4+3] = co[1][1];
		////////////////
		const double r = velo_scale;
		m_aXYVeloElem[ielem*4+0] = co[0][0] + r*velo[0][0];	
		m_aXYVeloElem[ielem*4+1] = co[0][1] + r*velo[0][1];
		m_aXYVeloElem[ielem*4+2] = co[1][0] + r*velo[1][0];
		m_aXYVeloElem[ielem*4+3] = co[1][1] + r*velo[1][1];
	}
	return true;
}

void CEdgeTextureColor::Draw() const
{
	::glColor3d(color[0], color[1], color[2]);
	::glBegin(GL_QUADS);
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		::glVertex2dv(m_aXYElem+ielem*4+2);
		::glVertex2dv(m_aXYElem+ielem*4+0);
		::glVertex2dv(m_aXYVeloElem+ielem*4+0);
		::glVertex2dv(m_aXYVeloElem+ielem*4+2);
	}
	::glEnd();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


CDrawerImageBasedFlowVis::CDrawerImageBasedFlowVis()
{
	nelem = 0;
	aXYVeloElem = 0;
	aXYElem = 0;
	////////////////
	iPtn = 0;
    alpha = 0;
	imode = 1;
    ////////////////
    velo_scale = 0.5;
	this->MakePattern();
}

CDrawerImageBasedFlowVis::CDrawerImageBasedFlowVis(const unsigned int id_field, 
												   const Fem::Field::CFieldWorld& world,
												   unsigned int imode)
{
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClear(GL_COLOR_BUFFER_BIT);
	////////////////
	nelem = 0;
	aXYVeloElem = 0;
	aXYElem = 0;
	////////////////
    this->imode = imode;
    if( imode == 0 ){ alpha = 0; }
	else{ alpha = 5; }
	iPtn = 0;
	m_nameDisplayList = 0;
    ////////////////
    velo_scale = 0.3;
    ////////////////
	this->MakePattern();
	this->Set( id_field, world );
}

void CDrawerImageBasedFlowVis::ClearDisplayList()
{
	if( m_nameDisplayList == 0 ) return;
	::glDeleteLists(m_nameDisplayList, m_nPattern);
	m_nameDisplayList = 0;
}

void CDrawerImageBasedFlowVis::MakePattern() 
{ 
	this->ClearDisplayList();
	const unsigned int NPN = 32;	// テクスチャのサイズ
	if( imode == 0 ){
		m_nPattern = 0;
		m_nameDisplayList = 0;
	}
	else if( imode == 1 )
	{
		m_nPattern = 32;
		m_nameDisplayList = glGenLists(m_nPattern);
		int lut[256];
		for(unsigned int i=0; i<256; i++){ lut[i] = (i < 127) ? 0 : 255; }
		int phase[NPN][NPN];
		for(unsigned int i=0; i<NPN; i++){
		for(unsigned int j=0; j<NPN; j++){ 
			phase[i][j] = rand() % 256; 
		}
		}
		GLubyte pat[NPN][NPN][4];
		for(unsigned int ipat=0; ipat<m_nPattern; ipat++){
			const unsigned int t = ipat*256/m_nPattern;
			for(unsigned int i=0; i<NPN; i++){
			for(unsigned int j=0; j<NPN; j++){
				GLubyte val = lut[(t + phase[i][j]) % 255];
				pat[i][j][0] = val;
				pat[i][j][1] = val;
				pat[i][j][2] = val;
				pat[i][j][3] = alpha;
			}
			}
			glNewList(m_nameDisplayList + ipat, GL_COMPILE);
			glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,   GL_RGBA, GL_UNSIGNED_BYTE, pat);
			glEndList();
		}
	}
	else if( imode == 2 )
	{
		m_nPattern = 1;
		m_nameDisplayList = glGenLists(m_nPattern);
		GLubyte pat[NPN][NPN][4];
		for(unsigned int i=0; i<NPN; i++){
		for(unsigned int j=0; j<NPN; j++){
			const int phase = ( abs((int)(i % 16)) < 3 || abs((int)(j % 16)) < 3 ) ?  0 : 255;
			pat[i][j][0] = phase;
			pat[i][j][1] = phase;
			pat[i][j][2] = phase;
            pat[i][j][3] = alpha;
		}
		}
		glNewList(m_nameDisplayList, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,   GL_RGBA, GL_UNSIGNED_BYTE, pat);
		glEndList();
	}
	else if( imode == 3 )
	{
		m_nPattern = 1;
		m_nameDisplayList = glGenLists(m_nPattern);
		GLubyte pat[NPN][NPN][4];
		for(unsigned int i=0; i<NPN; i++){
		for(unsigned int j=0; j<NPN; j++){			
			const int phase = ( abs((int)(i % 16)) < 3 && abs((int)(j % 16)) < 3 ) ?  0 : 255;
			pat[i][j][0] = phase;
			pat[i][j][1] = phase;
			pat[i][j][2] = phase;
			pat[i][j][3] = alpha;
		}
		}
		glNewList(m_nameDisplayList, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,   GL_RGBA, GL_UNSIGNED_BYTE, pat);
		glEndList();		
   }
}

bool CDrawerImageBasedFlowVis::Set(unsigned int id_field_velo, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field_velo) ){ return false; }
	this->m_IdFieldVelo = id_field_velo;
	const Fem::Field::CField& fv = world.GetField(id_field_velo);
	const std::vector<unsigned int>& aIdEA = fv.GetAry_IdElemAry();
	unsigned int nelem0 = 0;
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		assert( world.IsIdEA(id_ea) );
		const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
		nelem0 += ea.Size();
	}
	if( nelem != nelem0 ){
		nelem = nelem0;
		if( aXYVeloElem != 0 ){ delete[] aXYVeloElem; }
        if( aXYElem     != 0 ){ delete[] aXYElem; }
		aXYVeloElem = new double [nelem*6];
		aXYElem = new double [nelem*6];
	}
	return this->Update(world);
}
		
bool CDrawerImageBasedFlowVis::Update(const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(this->m_IdFieldVelo) ){ return false; }
	const Fem::Field::CField& fv = world.GetField(this->m_IdFieldVelo);
	const std::vector<unsigned int>& aIdEA = fv.GetAry_IdElemAry();
	{
		unsigned int nelem0 = 0;
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
			nelem0 += ea.Size();
		}
		if( nelem != nelem0 ){
			nelem = nelem0;
			if( aXYVeloElem != 0 ){ delete[] aXYVeloElem; }
			if( aXYElem     != 0 ){ delete[] aXYElem; }
			aXYVeloElem = new double [nelem*6];
			aXYElem     = new double [nelem*6];
		}
	}
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	unsigned int ielem_cur = 0;
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
	{
		const unsigned int id_ea = aIdEA[iiea];
		assert( world.IsIdEA(id_ea) );
		const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea,CORNER,true,world);
		const Fem::Field::CElemAry::CElemSeg& es_c = fv.GetElemSeg(id_ea,CORNER,false,world);
		for(unsigned int ielem=0;ielem<es_v.GetSizeElem();ielem++){
			unsigned int noes[3];
			es_c.GetNodes(ielem,noes);
			double co[3][2];
			ns_c.GetValue(noes[0],co[0]);
			ns_c.GetValue(noes[1],co[1]);
			ns_c.GetValue(noes[2],co[2]);
			////////////////
			es_v.GetNodes(ielem,noes);
			double velo[3][2];
			ns_v.GetValue(noes[0],velo[0]);
			ns_v.GetValue(noes[1],velo[1]);
			ns_v.GetValue(noes[2],velo[2]);
			////////////////
			aXYElem[ielem_cur*6+0] = co[0][0];
			aXYElem[ielem_cur*6+1] = co[0][1];
			aXYElem[ielem_cur*6+2] = co[1][0];
			aXYElem[ielem_cur*6+3] = co[1][1];
			aXYElem[ielem_cur*6+4] = co[2][0];
			aXYElem[ielem_cur*6+5] = co[2][1];
			////////////////
            const double r = velo_scale;
			aXYVeloElem[ielem_cur*6+0] = co[0][0] + r*velo[0][0];
			aXYVeloElem[ielem_cur*6+1] = co[0][1] + r*velo[0][1];
			aXYVeloElem[ielem_cur*6+2] = co[1][0] + r*velo[1][0];
			aXYVeloElem[ielem_cur*6+3] = co[1][1] + r*velo[1][1];
			aXYVeloElem[ielem_cur*6+4] = co[2][0] + r*velo[2][0];
			aXYVeloElem[ielem_cur*6+5] = co[2][1] + r*velo[2][1];
			ielem_cur++;
		}
	}
    assert( ielem_cur == nelem );
	for(unsigned int iec=0;iec<aEdgeColor.size();iec++){
		aEdgeColor[iec].Update(world);
	}
	return true;
}

Com::CBoundingBox CDrawerImageBasedFlowVis::GetBoundingBox( double rot[] ) const
{
	if( nelem == 0 ){ return Com::CBoundingBox(-1,1, -1,1, 0,0); }
	double x_min,x_max, y_min,y_max;
	x_min = x_max = aXYElem[0];
	y_min = y_max = aXYElem[1];
	for(unsigned int i=0;i<nelem*3;i++){
		const double x0 = aXYElem[i*2+0];
		const double y0 = aXYElem[i*2+1];
		x_min = ( x0 < x_min ) ? x0 : x_min;
		x_max = ( x0 > x_max ) ? x0 : x_max;
		y_min = ( y0 < y_min ) ? y0 : y_min;
		y_max = ( y0 > y_max ) ? y0 : y_max;
	}
	return Com::CBoundingBox(x_min,x_max, y_min,y_max, -1,1);
}

void CDrawerImageBasedFlowVis::Draw() const
{
	for(unsigned int iec=0;iec<aEdgeColor.size();iec++){ aEdgeColor[iec].Draw(); }

	////////////////////////////////
	// ModelView-Projection行列を取得
	GLint view[4];
	::glGetIntegerv(GL_VIEWPORT,view);
	const int win_w = view[2];
	const int win_h = view[3];
	const double invww = 1.0/win_w;
	const double invwh = 1.0/win_h;
	GLdouble model[16], prj[16];
	::glGetDoublev(GL_MODELVIEW_MATRIX,  model);
	::glGetDoublev(GL_PROJECTION_MATRIX, prj);
	////////////////
	::glEnable(GL_TEXTURE_2D);
	::glDisable(GL_DEPTH_TEST);
	if( imode != 0 ){	
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		::glDisable(GL_BLEND);
		::glColor4f(1.0f, 1.0f, 1.0f, 1.f);	
		::glBegin(GL_TRIANGLES);
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			const double* co_obj = &aXYElem[ielem*6];
			double co_win[3][2], dtmp1;
			::gluProject(co_obj[0*2+0],co_obj[0*2+1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
			::gluProject(co_obj[1*2+0],co_obj[1*2+1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
			::gluProject(co_obj[2*2+0],co_obj[2*2+1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);
			const double* velo_co = &aXYVeloElem[ielem*6];
			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2dv(velo_co+0);
			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2dv(velo_co+2);
			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2dv(velo_co+4);
		}
		::glEnd();
	}
	else{		
//		::glEnable(GL_BLEND); 
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		////////////////
//		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		::glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		::glBegin(GL_TRIANGLES);
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			const double* co_obj = &aXYElem[ielem*6];
			double co_win[3][2], dtmp1;
			::gluProject(co_obj[0*2+0],co_obj[0*2+1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
			::gluProject(co_obj[1*2+0],co_obj[1*2+1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
			::gluProject(co_obj[2*2+0],co_obj[2*2+1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);
			double r = +0.0;
			double diffuse[3][2];
			{
				double gc[2] = { (co_obj[0]+co_obj[2]+co_obj[4])/3.0, (co_obj[1]+co_obj[3]+co_obj[5])/3.0 };
				diffuse[0][0] = (co_obj[0] - gc[0])*r;
				diffuse[0][1] = (co_obj[1] - gc[1])*r;
				diffuse[1][0] = (co_obj[2] - gc[0])*r;
				diffuse[1][1] = (co_obj[3] - gc[1])*r;
				diffuse[2][0] = (co_obj[4] - gc[0])*r;
				diffuse[2][1] = (co_obj[5] - gc[1])*r;
			}
			const double* velo_co = &aXYVeloElem[ielem*6];
			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2d(velo_co[0]+diffuse[0][0], velo_co[1]+diffuse[0][1]);
			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2d(velo_co[2]+diffuse[1][0], velo_co[3]+diffuse[1][1]);
			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2d(velo_co[4]+diffuse[2][0], velo_co[5]+diffuse[2][1]);
		}
		::glEnd();
//		::glDisable(GL_BLEND); 

		/*
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND); 
		static const GLfloat blend[] = { 1.0, 1.0, 1.0, 0.0 };
		glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, blend);
		*/
		////////////////
		////////////////
		::glDisable(GL_DEPTH_TEST);
		::glEnable(GL_BLEND);		
		::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//		::glBlendFunc(GL_ONE, GL_ZERO);
//		::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//		::glBlendFunc(GL_ZERO,GL_ALPHA);
//		::glBlendFunc(GL_ALPHA,GL_ZERO);
//		::glBlendFunc(GL_ALPHA,GL_ALPHA);
//		::glBlendFunc(GL_SRC_ALPHA, GL_ZERO);
//		::glBlendFunc(GL_ZERO,GL_SRC_COLOR);
//		::glBlendFunc(GL_SRC_COLOR,GL_ZERO);
//		::glBlendFunc(GL_ONE_MINUS_SRC_COLOR,GL_SRC_COLOR);
//		::glBlendFunc(GL_DST_COLOR, GL_ZERO);
//		::glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ONE);
//		::glBlendFunc(GL_DST_COLOR,GL_ONE_MINUS_DST_COLOR);
		////////////////
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		::glColor4f(1.0f, 1.0f, 1.0f, 0.8f);
		::glBegin(GL_TRIANGLES);
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			const double* co_obj = &aXYElem[ielem*6];
			const double area = TriArea2D(co_obj,co_obj+2,co_obj+4);
			double co_win[3][2], dtmp1;
			::gluProject(co_obj[0*2+0],co_obj[0*2+1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
			::gluProject(co_obj[1*2+0],co_obj[1*2+1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
			::gluProject(co_obj[2*2+0],co_obj[2*2+1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);
			double diffuse[3][2];
/*			{
				double gc[2] = { (co_obj[0]+co_obj[2]+co_obj[4])/3.0, (co_obj[1]+co_obj[3]+co_obj[5])/3.0 };
				diffuse[0][0] = r*(1-2*(double)rand()/(RAND_MAX+1)); diffuse[0][1] = r*(1-2*(double)rand()/(RAND_MAX+1));
				diffuse[1][0] = r*(1-2*(double)rand()/(RAND_MAX+1)); diffuse[1][1] = r*(1-2*(double)rand()/(RAND_MAX+1));
				diffuse[2][0] = r*(1-2*(double)rand()/(RAND_MAX+1)); diffuse[2][1] = r*(1-2*(double)rand()/(RAND_MAX+1));
			}*/
/*			{
				double gc[2] = { (co_obj[0]+co_obj[2]+co_obj[4])/3.0, (co_obj[1]+co_obj[3]+co_obj[5])/3.0 };
				double tmpx = r*(1-2*(double)rand()/(RAND_MAX+1));
				double tmpy = r*(1-2*(double)rand()/(RAND_MAX+1));
				diffuse[0][0] = tmpx; diffuse[0][1] = tmpy;
				diffuse[1][0] = tmpx; diffuse[1][1] = tmpy;
				diffuse[2][0] = tmpx; diffuse[2][1] = tmpy;
			}*/
			{
				double r1 = sqrt(area)*0.5;
				double r2 = 1;
//				double r1 = 0;
//				double r2 = 0;
				double gc[2] = { (co_obj[0]+co_obj[2]+co_obj[4])/3.0, (co_obj[1]+co_obj[3]+co_obj[5])/3.0 };
				double tmpx = r1*(1-2*(double)rand()/(RAND_MAX+1));
				double tmpy = r1*(1-2*(double)rand()/(RAND_MAX+1));
				diffuse[0][0] = tmpx + (co_obj[0] - gc[0])*r2; 
				diffuse[0][1] = tmpy + (co_obj[1] - gc[1])*r2;
				diffuse[1][0] = tmpx + (co_obj[2] - gc[0])*r2; 
				diffuse[1][1] = tmpy + (co_obj[3] - gc[1])*r2;
				diffuse[2][0] = tmpx + (co_obj[4] - gc[0])*r2; 
				diffuse[2][1] = tmpy + (co_obj[5] - gc[1])*r2;
			}
//			const double* velo_co = &aXYVeloElem[ielem*6];
//			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2d(velo_co[0]+diffuse[0][0], velo_co[1]+diffuse[0][1]);
//			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2d(velo_co[2]+diffuse[1][0], velo_co[3]+diffuse[1][1]);
//			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2d(velo_co[4]+diffuse[2][0], velo_co[5]+diffuse[2][1]);			
			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2d(co_obj[0]+diffuse[0][0], co_obj[1]+diffuse[0][1]);
			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2d(co_obj[2]+diffuse[1][0], co_obj[3]+diffuse[1][1]);
			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2d(co_obj[4]+diffuse[2][0], co_obj[5]+diffuse[2][1]);
		}
		::glEnd();
		::glDisable(GL_BLEND);
	}
	iPtn = iPtn + 1;

	if( m_nPattern != 0 ){
		::glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND); 
		glCallList(iPtn % m_nPattern + m_nameDisplayList);
		::glBegin(GL_TRIANGLES);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			const double co[3][2] = {
				{ aXYElem[ielem*6+0], aXYElem[ielem*6+1] },
				{ aXYElem[ielem*6+2], aXYElem[ielem*6+3] },
				{ aXYElem[ielem*6+4], aXYElem[ielem*6+5] } };
/*			unsigned int r=4;
			::glTexCoord2d(r*co[0][0]-1,r*co[0][1]-1);		::glVertex2dv(co[0]);
			::glTexCoord2d(r*co[1][0]-1,r*co[1][1]-1);		::glVertex2dv(co[1]);
			::glTexCoord2d(r*co[2][0]-1,r*co[2][1]-1);		::glVertex2dv(co[2]);*/
			double r = 4.0;	// 流す画像の周期をきめる(大きければ細かい)
            ::glTexCoord2d(r*co[0][0],r*co[0][1]);		::glVertex2dv(co[0]);
            ::glTexCoord2d(r*co[1][0],r*co[1][1]);		::glVertex2dv(co[1]);
            ::glTexCoord2d(r*co[2][0],r*co[2][1]);		::glVertex2dv(co[2]);
		}
		::glEnd();
		::glDisable(GL_BLEND);
	}
	::glEnable(GL_DEPTH_TEST);
	
//	::glDisable(GL_BLEND);
	::glEnable(GL_BLEND); 
	::glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, win_w, win_h, 0);
	::glDisable(GL_BLEND);
}
