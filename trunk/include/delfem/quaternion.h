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
@brief 四元数クラス(Com::CQuaternion)のインターフェイス
@author Nobuyuki Umetani
*/


#if !defined(QUATERNION_H)
#define QUATERNION_H

#include <math.h>

#include "delfem/vector3d.h"

namespace Com{

class CQuaternion;

//! @{
CQuaternion operator+(const CQuaternion&, const CQuaternion&);	//!< 四元数同士の足し算
CQuaternion operator-(const CQuaternion&, const CQuaternion&);	//!< 四元数同士の引き算
CQuaternion operator*(double, const CQuaternion&);	//!< 四元数の実数倍
CQuaternion operator*(const CQuaternion&, double);	//!< 四元数の実数倍
CQuaternion operator/(const CQuaternion&, double);	//!< 四元数を実数で割る
CQuaternion operator*(const CQuaternion&, const CQuaternion&);	//!< 四元数同士の掛け算
CVector3D Rotate(const CQuaternion&, const CVector3D&);
CVector3D UnRotate(const CQuaternion&, const CVector3D&);
//! @}

//! 四元数(クオータニオン)クラス
class CQuaternion  
{
public:
	CQuaternion()
		:vector(){
		real = 1.0;
	}
	CQuaternion(const CQuaternion& rhs)
		:vector( rhs.vector ){
		real = rhs.real;
	}
	CQuaternion(const CVector3D& axis);
	CQuaternion(double real, const CVector3D& vector);
	CQuaternion(const CVector3D& a_vector, const CVector3D& b_vector);
	~CQuaternion();

	////////////////
	// Getメソッド

	CQuaternion GetConjugate() const;	//!< 共役クオータニオンを得る
	double GetReal() const{ return real; }	//!< 実部を得る
	CVector3D GetVector(){ return vector; }	//!< 虚部を得る

	//! 正規化
	void Normalize();
	//! 単位クオータニオンに設定
	void SetUnit(){ real = 1.0; vector.SetZero(); }
	//! 軸性ベクトルからの初期化
	void AxisToQuat(const CVector3D& axis);
	void VectorTrans(const CVector3D& a_vector, const CVector3D& b_vector);
    void RotMatrix33(double* m) const;

	friend bool operator==(const CQuaternion&, const CQuaternion&);
	friend bool operator!=(const CQuaternion&, const CQuaternion&);

	CQuaternion& operator=(const CQuaternion&);
	CQuaternion& operator+=(const CQuaternion&);
	CQuaternion& operator-=(const CQuaternion&);
	CQuaternion& operator*=(const CQuaternion&);
	CQuaternion& operator*=(double);
	CQuaternion& operator/=(const CQuaternion&);

	
private:
	CVector3D vector;	//!< 虚部
	double real;	//!< 実部
private:
	double Length();
};

}

#endif // !defined(QUATERNION_H)
