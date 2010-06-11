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
// DrawerField.cpp : èÍâ¬éãâªÉNÉâÉX(DrawerField)ÇÃé¿ëï
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
	if( this->m_paVer == 0 ) return Com::CBoundingBox();
	return m_paVer->GetBoundingBox(rot);
}

bool CDrawerVector::Update_VECTOR(const Fem::Field::CFieldWorld& world)	// ÉxÉNÉgÉãï`âÊ
{
	assert( world.IsIdField(id_field) );
	if( !world.IsIdField(id_field) ) return false;
	const Fem::Field::CField& field = world.GetField(id_field);
	if( field.IsPartial() ){
		std::cout  << "ñ¢é¿ëï" << std::endl;
		getchar();
		assert(0);
	}
	////////////////
	assert( field.GetFieldType() == VECTOR2 || field.GetFieldType() == VECTOR3 );
	nline = 0;
	{
		// í∏ì_îzóÒÇÉZÉbÉg(CORNER)
		unsigned int id_na_val_c = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		if( id_na_val_c !=  0 ){
			assert( world.IsIdNA(id_na_val_c) );
			const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val_c);
			nline += na_val.Size();
		}
		// í∏ì_îzóÒÇÉZÉbÉg(BUBBLE)
		unsigned int id_na_val_b = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		if( id_na_val_b != 0 ){
			assert( world.IsIdNA(id_na_val_b) );
			const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val_b);
			nline += na_val.Size();
		}
	}

	int ilayer_min, ilayer_max;
	{
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		if( aIdEA.size() > 0 ){
			ilayer_min = field.GetLayer(aIdEA[0]);
			ilayer_max = ilayer_min;
		}
		else{ ilayer_min=0; ilayer_max=0; }
		for(unsigned int iiea=1;iiea<aIdEA.size();iiea++){
			int ilayer = field.GetLayer(aIdEA[iiea]);
			ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
			ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
		}
	}

//	std::cout << ilayer_min << " " << ilayer_max << std::endl;

	const unsigned int ndim_co = field.GetNDimCoord();
	if( m_paVer == 0 ){
		if( ilayer_min == ilayer_max ){
			m_paVer = new Com::View::CVertexArray(2*nline,ndim_co); 
		}
		else{
			assert( ndim_co == 2 );
			m_paVer = new Com::View::CVertexArray(2*nline,3); 
		}
	}
	else{
		assert( m_paVer->NPoin() == 2*nline );
		if( ilayer_min == ilayer_max ){ assert( m_paVer->NDim() == ndim_co ); }
		else{ assert( ndim_co == 2 && m_paVer->NDim() == 3 ); }
	}

	const unsigned int ndim_va = m_paVer->NDim();

	unsigned int icoun = 0;
	// ÉRÅ[ÉiÅ[êﬂì_Ç…Ç¬Ç¢ÇƒÉxÉNÉgÉãÇçÏÇÈ
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
		assert( world.IsIdNA(nsna_c.id_na_co) );
		const Fem::Field::CNodeAry& na_c_co = world.GetNA(nsna_c.id_na_co);
		const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(nsna_c.id_ns_co);
		double coord[3],value[3];
		for(unsigned int ipoin=0;ipoin<npoin_va;ipoin++){
			unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
			ns_c_co.GetValue(ipoin_co,coord);
			ns_c_val.GetValue(ipoin,value);
			for(unsigned int idim=0;idim<ndim_co;idim++){	// ì_ÇOÇÃç¿ïWÇÉZÉbÉg
				m_paVer->pVertexArray[ipoin*2*ndim_va        +idim] = coord[idim];
				m_paVer->pVertexArray[ipoin*2*ndim_va+ndim_va+idim] = coord[idim]+value[idim];
			}
			if( ilayer_min != ilayer_max ){
				assert( ndim_co == 2 );
				assert( ndim_va == 3 );
            m_paVer->pVertexArray[ipoin*2*ndim_va        +2] = 0.01;
            m_paVer->pVertexArray[ipoin*2*ndim_va+ndim_va+2] = 0.01;
			}
		}
		if( ilayer_min != ilayer_max ){ // çÇÇ≥ÇçÏÇÈ
			const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
			for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
				const unsigned int id_ea = aIdEA[iiea];
				const int ilayer = field.GetLayer(id_ea);
				const double height = (ilayer+0.01-ilayer_min)/(1+ilayer_max-ilayer_min);
				const CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,CORNER,true,world);
				const unsigned int nnoes = es.GetSizeNoes();
				assert( nnoes < 16 );
				unsigned int noes[16];
				for(unsigned int ielem=0;ielem<es.GetSizeElem();ielem++){
					es.GetNodes(ielem,noes);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						unsigned int ipo0 = noes[inoes];
						m_paVer->pVertexArray[ipo0*2*ndim_va        +2] = height;
						m_paVer->pVertexArray[ipo0*2*ndim_va+ndim_va+2] = height;
					}
				}
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
		assert( ndim_co == ns_val.GetLength() );
		if( world.IsIdNA(nsna_b.id_na_co) ){
//			std::cout << "ÉoÉuÉãÇ≈ëŒâûÇ∑ÇÈç¿ïWêﬂì_Ç™Ç ÇÈèÍçá" << std::endl;
			const Fem::Field::CNodeAry& na_co = world.GetNA(nsna_b.id_na_co);
			const CNodeAry::CNodeSeg& ns_co = na_co.GetSeg(nsna_b.id_ns_co);
			double coord[3],value[3];
			for(unsigned int ipoin=0;ipoin<npoin_va;ipoin++){
				unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
				ns_co.GetValue(ipoin_co,coord);
				ns_val.GetValue(ipoin,value);
				for(unsigned int idim=0;idim<ndim_co;idim++){	// ì_ÇOÇÃç¿ïWÇÉZÉbÉg
					m_paVer->pVertexArray[(icoun+ipoin)*2*ndim_va       +idim] = coord[idim];
					m_paVer->pVertexArray[(icoun+ipoin)*2*ndim_va+ndim_va+idim] = coord[idim]+value[idim];
				}
			}
		}
		else{
//			std::cout << "ÉoÉuÉãÇ≈ç¿ïWÇ™óvëfÇÃíÜêSÇ…Ç ÇÈèÍçá" << std::endl;
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
					{	// óvëfÇÃíÜêSÇÃç¿ïWÇéÊìæ
						for(unsigned int idim=0;idim<ndim_co;idim++){ coord_cnt[idim] = 0.0; }
						const unsigned int nnoes = es_c_co.GetSizeNoes();
						es_c_co.GetNodes(ielem,noes);
						double coord[3];
						for(unsigned int inoes=0;inoes<nnoes;inoes++){
							unsigned int ipoi0 = noes[inoes];
							assert( ipoi0 < na_c_co.Size() );
							ns_c_co.GetValue(ipoi0,coord);
							for(unsigned int idim=0;idim<ndim_co;idim++){
								coord_cnt[idim] += coord[idim];
							}
						}
						for(unsigned int idim=0;idim<ndim_co;idim++){ coord_cnt[idim] /= nnoes; }
					}
					double value[3];
					{	// ÉoÉuÉãêﬂì_ÇÃílÇéÊìæ
						es_b_va.GetNodes(ielem,noes);
						unsigned int ipoi0 = noes[0];
						assert( ipoi0 < na_va.Size() );
						ns_val.GetValue(ipoi0,value);
					}
					for(unsigned int idim=0;idim<ndim_co;idim++){	// ì_ÇOÇÃç¿ïWÇÉZÉbÉg
						m_paVer->pVertexArray[icoun*2*ndim_va        +idim] = coord_cnt[idim];
						m_paVer->pVertexArray[icoun*2*ndim_va+ndim_va+idim] = coord_cnt[idim]+value[idim];
					}
					icoun++;
				} // end ielem
			} // end iei
		} // end if
	}
	return true;
}

bool CDrawerVector::Update_SSTR2(const Fem::Field::CFieldWorld& world)	// éÂâûóÕï\é¶
{
	assert( world.IsIdField(id_field) );
	if( !world.IsIdField(id_field) ) return false;
	const Fem::Field::CField& field = world.GetField(id_field);
	if( field.IsPartial() ){
		std::cout  << "ñ¢é¿ëï" << std::endl;
		getchar();
		assert(0);
	}
	////////////////
	assert( field.GetFieldType() == STSR2 );
	if( field.GetNodeSegInNodeAry(BUBBLE).id_na_va == 0 ) return false;
	nline = 0;
	{
		// í∏ì_îzóÒÇÉZÉbÉg(BUBBLE)
		unsigned int id_na_val_b = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		assert( world.IsIdNA(id_na_val_b) );
		const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val_b);
		nline = na_val.Size();
	}
	int ilayer_min, ilayer_max;
	{
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		if( aIdEA.size() > 0 ){
			ilayer_min = field.GetLayer(aIdEA[0]);
			ilayer_max = ilayer_min;
		}
		else{ ilayer_min=0; ilayer_max=0; }
		for(unsigned int iiea=1;iiea<aIdEA.size();iiea++){
			int ilayer = field.GetLayer(aIdEA[iiea]);
			ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
			ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
		}
	}
	assert( ilayer_min == ilayer_max );

	const unsigned int ndim = field.GetNDimCoord();
	if( m_paVer == 0 ){
		m_paVer = new Com::View::CVertexArray(2*nline,ndim); 
	}
	else{
		assert( m_paVer->NPoin() == 2*nline );
		assert( m_paVer->NDim()  == ndim );
	}

	unsigned int icoun = 0;

	assert( field.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 );
	{
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
		assert( ns_val.GetLength() == ndim*(ndim+1)/2 );
		{
//			std::cout << "ÉoÉuÉãÇ≈ç¿ïWÇ™óvëfÇÃíÜêSÇ…Ç ÇÈèÍçá" << std::endl;
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
					{	// óvëfÇÃíÜêSÇÃç¿ïWÇéÊìæ
						for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] = 0.0; }
						const unsigned int nnoes = es_c_co.GetSizeNoes();
						es_c_co.GetNodes(ielem,noes);
						double coord[3];
						for(unsigned int inoes=0;inoes<nnoes;inoes++){
							unsigned int ipoi0 = noes[inoes];
							assert( ipoi0 < na_c_co.Size() );
							ns_c_co.GetValue(ipoi0,coord);
							for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] += coord[idim]; }
						}
						for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] /= nnoes; }
					}
					double value[6];	assert( ns_val.GetLength() <= 6 );
					{	// ÉoÉuÉãêﬂì_ÇÃílÇéÊìæ
						es_b_va.GetNodes(ielem,noes);
						unsigned int ipoi0 = noes[0];
						assert( ipoi0 < na_va.Size() );
						ns_val.GetValue(ipoi0,value);
					}
					for(unsigned int idim=0;idim<ndim;idim++){	// ì_ÇOÇÃç¿ïWÇÉZÉbÉg
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

bool CDrawerVector::Update(const Fem::Field::CFieldWorld& world)
{
	const Fem::Field::CField& field = world.GetField(id_field);
	if( field.IsPartial() ){
		std::cout  << "ñ¢é¿ëï" << std::endl;
		getchar();
		assert(0);
	}
	////////////////
	if( field.GetFieldType() == VECTOR2 || 
	    field.GetFieldType() == VECTOR3    ){ return this->Update_VECTOR(world); }
	else if( field.GetFieldType() == STSR2 ){ return this->Update_SSTR2(world);  }
	return true;
}

bool CDrawerVector::Set(unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
	this->id_field = id_field;
	if( m_paVer != 0 ){ delete m_paVer; m_paVer=0; }
	this->Update(world);
	{
		const CField& field = world.GetField(id_field);
		unsigned int ndim = field.GetNDimCoord();
		if(      ndim == 2 ){ this->sutable_rot_mode = 1; }
		else if( ndim == 3 ){ this->sutable_rot_mode = 3; }
	}
	return true;
}

void CDrawerVector::Draw() const{
	if( nline == 0 ) return;
	assert( m_paVer != 0 );
   ::glEnable(GL_DEPTH_TEST);
	::glDisable(GL_TEXTURE_2D);
	::glColor3d(0.0,0.0,0.0);
	::glLineWidth(2);

   ::glEnableClientState(GL_VERTEX_ARRAY);
   ::glVertexPointer(m_paVer->NDim(),GL_DOUBLE,0,m_paVer->pVertexArray);
   if( m_paVer->NDim() == 2 ){ ::glTranslated(0,0,+0.01); }
   ::glDrawArrays(GL_LINES,0,nline*2);
   if( m_paVer->NDim() == 2 ){ ::glTranslated(0,0,-0.01); }
   ::glDisableClientState(GL_VERTEX_ARRAY);
}
