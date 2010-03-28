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

#include "delfem/drawer_field.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;


CIndexArrayElem::CIndexArrayElem(unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world)
{
//	std::cout << "CIndexArrayElem::CIndexArrayElem" << std::endl;
	itype = ELEM_TYPE_NOT_SET; 
	is_selected = false;
	this->id_ea=id_ea;
	this->id_es=id_es;
	color[0] = 0.8; color[1] = 0.8; color[2] = 0.8;
	nElem = 0;
	pIA_Elem = 0;
	pColor = 0;
	ilayer = 0;
	////////////////
	if( !world.IsIdEA(id_ea) ) return;
	const CElemAry& ea = world.GetEA(id_ea);
	if(      ea.ElemType() == Fem::Field::POINT){}
	else if( ea.ElemType() == Fem::Field::LINE ){ Set_Line(id_ea,id_es, world); }
	else if( ea.ElemType() == Fem::Field::TRI  ){ Set_Tri( id_ea,id_es, world); }
	else if( ea.ElemType() == Fem::Field::QUAD ){ Set_Quad(id_ea,id_es, world); }
	else if( ea.ElemType() == Fem::Field::TET  ){ Set_Tet( id_ea,id_es, world); }
	else if( ea.ElemType() == Fem::Field::HEX  ){ Set_Hex( id_ea,id_es, world); }
}

////////////////////////////////////////////////////////////////

/*
void View::CIndexArrayElem::DrawElements()
{
}
*/

bool View::CIndexArrayElem::Set_Line(unsigned int id_ea, unsigned int id_es, 
									 const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA(id_ea) ) return false;
	const CElemAry& ea = world.GetEA(id_ea);
	if( !ea.IsSegID(id_es) ) return false;
	if( ea.ElemType() != LINE ) return false;
	////////////////
	itype = Fem::Field::LINE;
	this->id_ea = id_ea;
	this->id_es = id_es;
	const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
	nElem = ea.Size();
	if( this->pIA_Elem != 0 ) delete[] pIA_Elem;
	pIA_Elem = new unsigned int [nElem*2];
	unsigned int nnoes = es.GetSizeNoes();
	assert( nnoes == 2 );
	for(unsigned int iedge=0;iedge<nElem;iedge++){
		es.GetNodes(iedge,pIA_Elem+iedge*2);
	}
	return true;
}

bool View::CIndexArrayElem::Set_Tri(unsigned int id_ea, unsigned int id_es, 
									const Fem::Field::CFieldWorld& world)
{	
//	std::cout << "View::CIndexArrayElem::Set_Tri" << std::endl;
	if( !world.IsIdEA(id_ea) ) return false;
	const CElemAry& ea = world.GetEA(id_ea);
	if( !ea.IsSegID(id_es) ) return false;
	if( ea.ElemType() != TRI ) return false;
	////////////////
	itype = Fem::Field::TRI;
	this->id_ea = id_ea;
	this->id_es = id_es;
	const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
	nElem = ea.Size();
	if( this->pIA_Elem != 0 ) delete[] pIA_Elem;
	pIA_Elem = new unsigned int [nElem*3];
	unsigned int nnoes = es.GetSizeNoes();
	assert( nnoes == 3 );
	for(unsigned int itri=0;itri<nElem;itri++){
		es.GetNodes(itri,pIA_Elem+itri*3);
	}
	return true;
}

bool View::CIndexArrayElem::Set_Quad(unsigned int id_ea, unsigned int id_es, 
									 const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA(id_ea) ) return false;
	const CElemAry& ea = world.GetEA(id_ea);
	if( !ea.IsSegID(id_es) ) return false;
	if( ea.ElemType() != QUAD ) return false;
	////////////////
	itype = Fem::Field::QUAD;
	this->id_ea = id_ea;
	this->id_es = id_es;
	const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
	nElem = ea.Size();
	if( this->pIA_Elem != 0 ) delete[] pIA_Elem;
	pIA_Elem = new unsigned int [nElem*4];
	const unsigned int nnoes = es.GetSizeNoes();
	assert( nnoes == 4 );
	for(unsigned int iquad=0;iquad<nElem;iquad++){
		es.GetNodes(iquad,pIA_Elem+iquad*4);
	}
	return true;
}

bool View::CIndexArrayElem::Set_Tet(unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA(id_ea) ) return false;
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.IsSegID(id_es) );
	if( !ea.IsSegID(id_es) ) return false;
	if( ea.ElemType() != TET ) return false;
	////////////////
	itype = Fem::Field::TET;
	this->id_ea = id_ea;
	this->id_es = id_es;
	unsigned int id_es_add = 0;
	CElemAry* pEA = ea.MakeBoundElemAry(id_es,id_es_add,aIndElem);
	const CElemAry::CElemSeg& es = pEA->GetSeg(id_es_add);
	nElem = pEA->Size();
	const unsigned int npofa = 3;
	pIA_Elem = new unsigned int [nElem*npofa];
	for(unsigned int iface=0;iface<nElem;iface++){
		es.GetNodes(iface,pIA_Elem+iface*npofa);
	}
	delete pEA;
	return true;
}

bool View::CIndexArrayElem::Set_Hex(unsigned int id_ea, unsigned int id_es, 
									const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA(id_ea) ) return false;
	const CElemAry& ea = world.GetEA(id_ea);
	if( !ea.IsSegID(id_es) ) return false;
	if( ea.ElemType() != HEX ) return false;
	////////////////
	itype = Fem::Field::HEX;
	this->id_ea = id_ea;
	this->id_es = id_es;
	unsigned int id_es_add = 0;
	CElemAry* pEA = ea.MakeBoundElemAry(id_es,id_es_add,aIndElem);
	nElem = pEA->Size();
	assert( aIndElem.size() == nElem );
	assert( pEA->IsSegID(id_es_add) );
	const CElemAry::CElemSeg& es = pEA->GetSeg(id_es_add);
	const unsigned int npofa = 4;
	pIA_Elem = new unsigned int [nElem*npofa];
	for(unsigned int iface=0;iface<nElem;iface++){
		es.GetNodes(iface,pIA_Elem+iface*npofa);
	}
	delete pEA;
	return true;
}

////////////////////////////////////////////////////////////////

bool View::CIndexArrayElem::SetColor(
		unsigned int id_es_v, unsigned int id_ns_v, const Fem::Field::CFieldWorld& world,
		const std::auto_ptr<CColorMap>& color_map )
{
	if( itype == Fem::Field::TRI )	// TRI
	{
		if( !world.IsIdEA(id_ea) ) return false;
		const CElemAry& ea = world.GetEA(id_ea);
		if( !ea.IsSegID(id_es) ) return false;
		if( ea.ElemType() != TRI ) return false;
		const CElemAry::CElemSeg& es_v = ea.GetSeg(id_es_v);
		assert( es_v.GetSizeNoes() == 1 );
		unsigned int id_na_v = es_v.GetIdNA();
		assert( world.IsIdNA(id_na_v) );
		const CNodeAry& na = world.GetNA(id_na_v);
		assert( es_v.GetMaxNoes() < na.Size() );
		assert( na.IsSegID(id_ns_v) );
		const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_v);
		assert( ns.GetLength() == 1 );
		if( this->pColor == 0 ){ this->pColor = new float [nElem*3]; }
		for(unsigned int itri=0;itri<nElem;itri++){
			unsigned int inode0;
			es_v.GetNodes(itri,&inode0);
			double val;
			ns.GetValue(inode0,&val);
			color_map->GetColor(pColor+itri*3,val);
		}
		return true;
	}
	else if( itype == Fem::Field::QUAD )
	{
		if( !world.IsIdEA(id_ea) ) return false;
		const CElemAry& ea = world.GetEA(id_ea);
		if( !ea.IsSegID(id_es) ) return false;
		if( ea.ElemType() != QUAD ) return false;
		const CElemAry::CElemSeg& es_v = ea.GetSeg(id_es_v);
		assert( es_v.GetSizeNoes() == 1 );
		unsigned int id_na_v = es_v.GetIdNA();
		assert( world.IsIdNA(id_na_v) );
		const CNodeAry& na = world.GetNA(id_na_v);
		assert( es_v.GetMaxNoes() < na.Size() );
		assert( na.IsSegID(id_ns_v) );
		const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_v);
		assert( ns.GetLength() == 1 );
		if( this->pColor == 0 ){ this->pColor = new float [nElem*3]; }
		for(unsigned int iquad=0;iquad<nElem;iquad++){
			unsigned int inode0;
			es_v.GetNodes(iquad,&inode0);
			double val;
			ns.GetValue(inode0,&val);
			color_map->GetColor(pColor+iquad*3,val);
		}
		return true;
	}
	else if( itype == Fem::Field::TET )
	{
		if( !world.IsIdEA(id_ea) ) return false;
		const CElemAry& ea = world.GetEA(id_ea);
		const unsigned int ntet = ea.Size();
		if( !ea.IsSegID(id_es) ) return false;
		if( ea.ElemType() != TET ) return false;
		const CElemAry::CElemSeg& es_v = ea.GetSeg(id_es_v);
		assert( es_v.GetSizeNoes() == 1 );
		unsigned int id_na_v = es_v.GetIdNA();
		assert( world.IsIdNA(id_na_v) );
		const CNodeAry& na = world.GetNA(id_na_v);
		assert( es_v.GetMaxNoes() < na.Size() );
		assert( na.IsSegID(id_ns_v) );
		const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_v);
		assert( ns.GetLength() == 1 );
		if( this->pColor == 0 ){ this->pColor = new float [nElem*3]; }
		for(unsigned int iface=0;iface<nElem;iface++){
			unsigned int itet = aIndElem[iface];
			assert( itet < ntet );
			unsigned int inode0;
			es_v.GetNodes(itet,&inode0);
			double val;
			ns.GetValue(inode0,&val);
			color_map->GetColor(pColor+iface*3,val);
		}
		return true;
	}
	else if( itype == Fem::Field::HEX )
	{	
		if( !world.IsIdEA(id_ea) ) return false;
		const CElemAry& ea = world.GetEA(id_ea);
		const unsigned int nhex = ea.Size();
		if( !ea.IsSegID(id_es) ) return false;
		if( ea.ElemType() != HEX ) return false;
		const CElemAry::CElemSeg& es_v = ea.GetSeg(id_es_v);
		assert( es_v.GetSizeNoes() == 1 );
		unsigned int id_na_v = es_v.GetIdNA();
		assert( world.IsIdNA(id_na_v) );
		const CNodeAry& na = world.GetNA(id_na_v);
		assert( es_v.GetMaxNoes() < na.Size() );
		assert( na.IsSegID(id_ns_v) );
		const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_v);
		assert( ns.GetLength() == 1 );
		if( this->pColor == 0 ){ this->pColor = new float [nElem*3]; }
		for(unsigned int iface=0;iface<nElem;iface++){
			unsigned int ihex = aIndElem[iface];
			assert( ihex < nhex );
			unsigned int inode0;
			es_v.GetNodes(ihex,&inode0);
			double val;
			ns.GetValue(inode0,&val);			
			color_map->GetColor(pColor+iface*3,val);
		}
		return true;
	}
	return true;
}

////////////////////////////////////////////////////////////////

void View::CIndexArrayElem::DrawElements()
{
	if( this->pColor == 0 ){ 
		::glColor3d(color[0],color[1],color[2]);
		if(      itype == Fem::Field::LINE )	// line
		{
			::glDrawElements(GL_LINES,    nElem*2,GL_UNSIGNED_INT,this->pIA_Elem);
		}
		else if( itype == Fem::Field::TRI  || itype == Fem::Field::TET )	// tri 
		{
			::glDrawElements(GL_TRIANGLES,nElem*3,GL_UNSIGNED_INT,this->pIA_Elem);
		}
		else if( itype == Fem::Field::QUAD || itype == Fem::Field::HEX )	// quad
		{	
			::glDrawElements(GL_QUADS,    nElem*4,GL_UNSIGNED_INT,this->pIA_Elem);
		}
		return;
	}
	if( itype == Fem::Field::QUAD || itype == Fem::Field::HEX ){
		::glBegin(GL_QUADS);
		for(unsigned int iface=0;iface<nElem;iface++){
			::glColor3fv( pColor+iface*3 );
			::glArrayElement( pIA_Elem[iface*4  ] );
			::glArrayElement( pIA_Elem[iface*4+1] );
			::glArrayElement( pIA_Elem[iface*4+2] );
			::glArrayElement( pIA_Elem[iface*4+3] );
		}
		::glEnd();
	}	
	else if( itype == Fem::Field::TRI || itype == Fem::Field::TET )
	{
		::glBegin(GL_TRIANGLES);
		for(unsigned int itri=0;itri<nElem;itri++){
			::glColor3fv( pColor+itri*3 );
			::glArrayElement( pIA_Elem[itri*3  ] );
			::glArrayElement( pIA_Elem[itri*3+1] );
			::glArrayElement( pIA_Elem[itri*3+2] );
		}
		::glEnd();
	}
}
