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

#include "delfem/drawer_field_vector.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

CDrawerVector::CDrawerVector(){
	m_paVer = 0;
	nline = 0;
}

CDrawerVector::CDrawerVector(unsigned int id_field, const Fem::Field::CFieldWorld& world){
	m_paVer = 0;
	nline = 0;
	this->Set(id_field, world);
}

CDrawerVector::~CDrawerVector(){
	if( m_paVer != 0 ){ delete m_paVer; }
}
Com::CBoundingBox CDrawerVector::GetBoundingBox( double rot[] ) const{
	return m_paVer->GetBoundingBox(rot);
}

bool CDrawerVector::Update(const Fem::Field::CFieldWorld& world){
	const Fem::Field::CField& field = world.GetField(id_field);

	if( field.IsPartial() ){
		std::cout  << "未実装" << std::endl;
		getchar();
		assert(0);
	}

	nline = 0;
	{
		{	// 頂点配列をセット
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			if( id_na_val !=  0 ){
				assert( world.IsIdNA(id_na_val) );
				const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val);
				nline += na_val.Size();
			}
		}	
		{	// 頂点配列をセット
			unsigned int id_na_val = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
			if( id_na_val != 0 ){
				assert( world.IsIdNA(id_na_val) );
				const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val);
				nline += na_val.Size();
			}
		}
	}

	const unsigned int ndim = field.GetNDimCoord();
	if( m_paVer == 0 ){	m_paVer = new Com::View::CVertexArray(2*nline,ndim); }
	else if( m_paVer->NDim() != ndim || m_paVer->NPoin() != 2*nline ){ 
		assert(0);
		return false; 
	}

	unsigned int icoun = 0;
	// コーナー節点についてベクトルを作る
	if( field.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
		const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(CORNER);
		assert( world.IsIdNA(nsna_c.id_na_va) );
		const Fem::Field::CNodeAry& na_c_val = world.GetNA(nsna_c.id_na_va);
		const unsigned int npoin_va = na_c_val.Size();
		unsigned int id_ns_c_v;
		{
			if(      nsna_c.id_ns_va != 0 ) id_ns_c_v = nsna_c.id_ns_va;
			else if( nsna_c.id_ns_ve != 0 ) id_ns_c_v = nsna_c.id_ns_ve;
			else if( nsna_c.id_ns_ac != 0 ) id_ns_c_v = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		assert( na_c_val.IsSegID(id_ns_c_v) );
		const Fem::Field::CNodeAry::CNodeSeg& ns_c_val = na_c_val.GetSeg(id_ns_c_v);
		assert( ndim == ns_c_val.GetLength() );
		assert( world.IsIdNA(nsna_c.id_na_co) );
		const Fem::Field::CNodeAry& na_c_co = world.GetNA(nsna_c.id_na_co);
		const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(nsna_c.id_ns_co);
		double coord[3],value[3];
		for(unsigned int ipoin=0;ipoin<npoin_va;ipoin++){
			unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
			ns_c_co.GetValue(ipoin_co,coord);
			ns_c_val.GetValue(ipoin,value);
			for(unsigned int idim=0;idim<ndim;idim++){	// 点０の座標をセット
				m_paVer->pVertexArray[ipoin*2*ndim     +idim] = coord[idim];
				m_paVer->pVertexArray[ipoin*2*ndim+ndim+idim] = coord[idim]+value[idim];
			}
		}
		icoun = npoin_va;
	}
	if( field.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 ){
		const CField::CNodeSegInNodeAry& nsna_b = field.GetNodeSegInNodeAry(BUBBLE);
		assert( world.IsIdNA(nsna_b.id_na_va) );
		const Fem::Field::CNodeAry& na_val = world.GetNA(nsna_b.id_na_va);
		const unsigned int npoin_va = na_val.Size();
		unsigned int id_ns_v;
		{
			if(      nsna_b.id_ns_va != 0 ) id_ns_v = nsna_b.id_ns_va;
			else if( nsna_b.id_ns_ve != 0 ) id_ns_v = nsna_b.id_ns_ve;
			else if( nsna_b.id_ns_ac != 0 ) id_ns_v = nsna_b.id_ns_ac;
			else{ assert(0); }
		}
		assert( na_val.IsSegID(id_ns_v) );
		const Fem::Field::CNodeAry::CNodeSeg& ns_val = na_val.GetSeg(id_ns_v);
		assert( ndim == ns_val.GetLength() );
		if( world.IsIdNA(nsna_b.id_na_co) ){
//			std::cout << "バブルで対応する座標節点がある場合" << std::endl;
			const Fem::Field::CNodeAry& na_co = world.GetNA(nsna_b.id_na_co);
			const CNodeAry::CNodeSeg& ns_co = na_co.GetSeg(nsna_b.id_ns_co);
			double coord[3],value[3];
			for(unsigned int ipoin=0;ipoin<npoin_va;ipoin++){
				unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
				ns_co.GetValue(ipoin_co,coord);
				ns_val.GetValue(ipoin,value);
				for(unsigned int idim=0;idim<ndim;idim++){	// 点０の座標をセット
					m_paVer->pVertexArray[(icoun+ipoin)*2*ndim     +idim] = coord[idim];
					m_paVer->pVertexArray[(icoun+ipoin)*2*ndim+ndim+idim] = coord[idim]+value[idim];
				}
			}
		}
		else{
//			std::cout << "バブルで座標が要素の中心にある場合" << std::endl;
			unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
			unsigned int id_ns_c_co = field.GetNodeSegInNodeAry(CORNER).id_ns_co;
			assert( world.IsIdNA(id_na_c_co) );
			const CNodeAry& na_c_co = world.GetNA(id_na_c_co);
			assert( na_c_co.IsSegID(id_ns_c_co) );
			const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(id_ns_c_co);
			assert( world.IsIdNA(nsna_b.id_na_va) );
			const CNodeAry& na_va = world.GetNA(nsna_b.id_na_va);
			assert( na_va.IsSegID(id_ns_v) );
			unsigned int noes[64];
			const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
			for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
				unsigned int id_ea = aIdEA[iiea];
				const CElemAry& ea = world.GetEA(id_ea);
				const CElemAry::CElemSeg& es_c_co = field.GetElemSeg(id_ea,CORNER,false,world);
				assert( es_c_co.GetIdNA() == id_na_c_co );
				const CElemAry::CElemSeg& es_b_va = field.GetElemSeg(id_ea,BUBBLE,true,world);
				assert( es_b_va.GetIdNA() == nsna_b.id_na_va );
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					double coord_cnt[3];
					{	// 要素の中心の座標を取得
						for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] = 0.0; }
						const unsigned int nnoes = es_c_co.GetSizeNoes();
						es_c_co.GetNodes(ielem,noes);
						double coord[3];
						for(unsigned int inoes=0;inoes<nnoes;inoes++){
							unsigned int ipoi0 = noes[inoes];
							assert( ipoi0 < na_c_co.Size() );
							ns_c_co.GetValue(ipoi0,coord);
							for(unsigned int idim=0;idim<ndim;idim++){
								coord_cnt[idim] += coord[idim];
							}
						}
						for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] /= nnoes; }
					}
					double value[3];
					{	// バブル節点の値を取得
						es_b_va.GetNodes(ielem,noes);
						unsigned int ipoi0 = noes[0];
						assert( ipoi0 < na_va.Size() );
						ns_val.GetValue(ipoi0,value);
					}
					for(unsigned int idim=0;idim<ndim;idim++){	// 点０の座標をセット
						m_paVer->pVertexArray[icoun*2*ndim     +idim] = coord_cnt[idim];
						m_paVer->pVertexArray[icoun*2*ndim+ndim+idim] = coord_cnt[idim]+value[idim];
					}
					icoun++;
				} // end ielem
			} // end iei
		} // end if
	}	
	return true;
}

bool CDrawerVector::Set(unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
	this->id_field = id_field;
	if( m_paVer != 0 ){ delete m_paVer; m_paVer=0; }
	this->Update(world);
	unsigned int ndim = m_paVer->NDim();
	if( ndim == 2 ){
		this->sutable_rot_mode = 1;
	}
	else if( ndim == 3 ){
		this->sutable_rot_mode = 3;
	}
	return true;
}

void CDrawerVector::Draw() const{
	if( nline == 0 ) return;
	assert( m_paVer != 0 );
	::glDisable(GL_TEXTURE_2D);
	::glColor3d(0.0,0.0,0.0);
	::glLineWidth(2);
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(m_paVer->NDim(),GL_DOUBLE,0,m_paVer->pVertexArray);
	::glDrawArrays(GL_LINES,0,nline*2);
	::glDisableClientState(GL_VERTEX_ARRAY);
}
