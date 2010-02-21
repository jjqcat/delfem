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
@brief St.Venant-Kirchhoff体の要素剛性作成部のインターフェース
@author Nobuyuki Umetani
@sa http://ums.futene.net/wiki/FEM/46454D20666F722053742E56656E616E742D4B69726368686F6666204D6174657269616C.html
*/



#if !defined(EQN_ST_VENSNT_H)
#define EQN_ST_VENANT_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include "delfem/linearsystem_interface_eqnsys.h"

namespace Fem{
namespace Ls{
	class CLinearSystem_Field;
}
namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn
{
////////////////////////////////
// ２Ｄの方程式

/*! @defgroup eqn_stvenant St.Venant-Kirchhoff体の連立一次方程式へのマージ関数群
@ingroup FemEqnMargeFunction
　
@f$ S = \frac{\partial W}{\partial E}@f$
*/
/*!@{*/

/*!
@brief ２次元静的St.Venant-Kirchhoff体のマージ
@param [in,out] ls マージされる連立一次方程式
@param lambda [in] ラメ第一定数 @f$ \lambda @f$
@param myu [in] ラメ第二定数 @f$ \mu @f$
@param rho [in] 質量密度
@
*/
bool AddLinSys_StVenant2D_Static(
		Fem::Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu,
		double rho, double f_x, double f_y,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp,
		unsigned int id_ea = 0);

/*!
@brief ２次元動的St.Venant-Kirchhoff体のマージ
@param [in,out] ls 連立一次方程式
@param lambda [in] ラメ第一定数 @f$ \lambda @f$
@param myu [in] ラメ第二定数 @f$ \mu @f$
@param rho [in] 質量密度
*/
bool AddLinSys_StVenant2D_NonStatic_NewmarkBeta(
		double dt, double gamma_newmark, double beta, 
        Fem::Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu,
		double rho, double g_x, double g_y,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp, 
		bool is_initial,
		unsigned int id_ea = 0 );
		
////////////////////////////////
// ３Ｄの方程式

/*!
@brief ３次元静的St.Venant-Kirchhoff体のマージ
@param [in,out] ls 連立一次方程式
@param lambda [in] ラメ第一定数 @f$ \lambda @f$
@param myu [in] ラメ第二定数 @f$ \mu @f$
@param rho [in] 質量密度
@param g_x [in] x方向の体積力
@param g_y [in] y方向の体積力
@param g_z [in] z方向の体積力
*/
bool AddLinSys_StVenant3D_Static(
        Fem::Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp,
		unsigned int id_ea = 0 );

/*!
@brief ３次元動的St.Venant-Kirchhoff体のマージ
@param [in,out] ls 連立一次方程式
@lambda [in] ラメ第一定数
@myu [in] ラメ第二定数
@rho [in] 質量密度
@g_x [in] x方向の体積力
@g_y [in] y方向の体積力
@g_z [in] z方向の体積力
*/
bool AddLinSys_StVenant3D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta,
        Fem::Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp,
		bool is_initial,
		unsigned int id_ea = 0);

/*!@}*/
}
}

#endif
