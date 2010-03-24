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
@brief べクトルで速度場を可視化するクラス(Fem::Field::View::CDrawerVector)のインターフェース
@author Nobuyuki Umetani
*/


#if !defined(DRAWER_FIELD_VECTOR_H)
#define DRAWER_FIELD_VECTOR_H

#include "delfem/drawer_field.h"

namespace Fem{
namespace Field{
namespace View{
	
//! ベクトルの描画クラス
class CDrawerVector : public CDrawerField
{
public:
	CDrawerVector();
	CDrawerVector(unsigned int id_field, const Fem::Field::CFieldWorld& world );
	virtual ~CDrawerVector();
	Com::CBoundingBox GetBoundingBox( double rot[] ) const;
	virtual void DrawSelection(unsigned int idraw) const{};
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual void Draw() const;
	virtual bool Update(const Fem::Field::CFieldWorld& world);
private:
	bool Set(unsigned int id_field, const Fem::Field::CFieldWorld& world );
private:
	unsigned int id_field;
	unsigned int nline;	// 辺の数
	Com::View::CVertexArray* m_paVer;	// 頂点配列
};

}
}
}

#endif
