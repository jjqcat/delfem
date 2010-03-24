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


/*! @file
@brief 場可視化クラス(Fem::Field::View::CDrawerField, Fem::Field::View::CDrawerFace)のインターフェース
@author Nobuyuki Umetani
*/


#if !defined(DRAWER_FIELD_H)
#define DRAWER_FIELD_H

#include "delfem/field_world.h"
#include "delfem/drawer.h"
#include "delfem/drawer_gl_utility.h"
#include "delfem/elem_ary.h"

namespace Fem{
namespace Field{
namespace View{

//! OpenGL頂点配列用のIndexArray格納クラス
class CIndexArrayElem
{
public : 
	CIndexArrayElem(unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world);
	
	~CIndexArrayElem(){
		if( pIA_Elem != 0 ){ delete[] pIA_Elem; }
		if( pColor   != 0 ){ delete[] pColor; }
	}
	void DrawElements();

	unsigned int GetElemDim() const {
		if( itype == 2 ){ return 1; }
		if( itype == 3 || itype == 4 ){ return 2; }
		if( itype == 5 || itype == 6 ){ return 3; }
		return 0;
	}
	unsigned int GetIdEA() const { return id_ea; }
	void SetColor(double r, double g, double b){ color[0]=r; color[1]=g; color[2]=b; }
	bool SetColor(unsigned int id_es_v, unsigned int id_ns_v, const Fem::Field::CFieldWorld& world, double min, double max);
private:
	bool Set_Line(unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world);
	bool Set_Tri( unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world);
	bool Set_Quad(unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world);
	bool Set_Tet( unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world);
	bool Set_Hex( unsigned int id_ea, unsigned int id_es, const Fem::Field::CFieldWorld& world);
public:
	int ilayer;
private:
	// これはelem_ary.hで定義されているELEM_TYPEを使うべき
	Fem::Field::ELEM_TYPE itype; // 0:none 1:vertex 2:line 3:tri 4:quad 5:tet 6:hex	
	bool is_selected;
	unsigned int id_ea;
	unsigned int id_es;
    double color[3];
	////////////////
	unsigned int nElem;
	unsigned int* pIA_Elem;
	////////////////
	float* pColor;
	std::vector<unsigned int> aIndElem;	// tetやhexの面から，要素idexへのマップ
};

//! Vertexのインデックスを格納するクラス
class CIndexVertex{
public:
	CIndexVertex(unsigned int id_v, unsigned int id_ea, unsigned int id_es) : id_v(id_v), id_ea(id_ea), id_es(id_es){}
	unsigned int id_v;
	unsigned int id_ea;
	unsigned int id_es;
	bool is_selected;
};

////////////////////////////////////////////////////////////////

//! 場可視化クラスの抽象クラス
class CDrawerField : public Com::View::CDrawer 
{
public:
	virtual bool Update(const Fem::Field::CFieldWorld& world) = 0;
};

////////////////////////////////////////////////////////////////

//! DrawerField配列クラス。(C++のデモプログラムからのみ使われる)
class CDrawerArrayField : public Com::View::CDrawerArray
{
public:
    virtual ~CDrawerArrayField(){}
	bool Update(const Fem::Field::CFieldWorld& world){
		for(unsigned int idraw=0;idraw<m_drawer_ary.size();idraw++){
			CDrawerField* pDF = static_cast<CDrawerField*>(this->m_drawer_ary[idraw]);
			if( pDF != 0 ){ pDF->Update(world); }
		}
		return true;
	}
	virtual void PushBack(Com::View::CDrawer* pDrawer){
		assert( pDrawer != 0 );
		this->m_drawer_ary.push_back(pDrawer);
	}
};

}	// end namespace View
}	// end namespace Feild
}	// end namespace Fem

#endif
