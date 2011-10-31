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
// Vector3D.cpp: CVector3D クラスのインプリメンテーション
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#endif

#include <cassert>
#include <math.h>
#include <iostream>
#include <stack>

#include "delfem/vector3d.h"

using namespace Com;

////////////////////////////////////////////////////////////////////
// メンバ関数のフレンド関数
////////////////////////////////////////////////////////////////////

namespace Com{

bool operator == (const CVector3D& lhs, const CVector3D& rhs){
	if( fabs(lhs.x - rhs.x) < NEARLY_ZERO
		&& fabs(lhs.y - rhs.y) < NEARLY_ZERO
		&& fabs(lhs.z - rhs.z) < NEARLY_ZERO )
		return true;
	else return false;
}

bool operator != (const CVector3D& lhs, const CVector3D& rhs){
	if( lhs == rhs )	return false;
	else return true;
}


}

//////////////////////////////////////////////////////////////////////
//	メンバ関数の非フレンド関数
//////////////////////////////////////////////////////////////////////

namespace Com{

void CVector3D::SetNormalizedVector()
{
	double mag;

	mag = Length();
	x /= mag;
	y /= mag;
	z /= mag;
}

void CVector3D::SetZero()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

}

