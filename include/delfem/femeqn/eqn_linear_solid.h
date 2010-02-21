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
@brief 線形弾性体の要素剛性作成部のインターフェース
@author Nobuyuki Umetani
@sa http://ums.futene.net/wiki/FEM/46454D20666F72204C696E65617220456C6173746963.html
*/

#if !defined(EQN_LINEAR_SOLID_H)
#define EQN_LINEAR_SOLID_H

#include <vector>

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

namespace Fem{
namespace Ls
{
class CLinearSystem_Field;
class CLinearSystem_Save;
class CLinearSystem_SaveDiaM_NewmarkBeta;
class CLinearSystem_Eigen;
class CPreconditioner;
}

namespace Field
{
class CField;
class CFieldWorld;
}

namespace Eqn{

class CLinearSystem_EqnInterface;

/*! @defgroup eqn_linear_solid 線形弾性体の方程式をマージする関数群
@ingroup FemEqnMargeFunction
　
*/
//! @{

////////////////////////////////////////////////////////////////
// ２次元問題

// 静的線形弾性体
bool AddLinSys_LinearSolid2D_Static(
		Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu, double rho, double g_x, double g_y,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp,
		unsigned int id_ea = 0 );

// 静的　熱応力　線形弾性体
bool AddLinSys_LinearSolidThermalStress2D_Static(
		Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu, double rho, double g_x, double g_y, 	double thermoelastic, 
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, unsigned int id_field_temp,
		unsigned int id_ea = 0 );

// 動的線形弾性体
bool AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta, Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu, double rho, double g_x, double g_y,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
		bool is_initial = true,
		unsigned int id_ea = 0 );

// 動的　熱応力　線形弾性体
bool AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta, Fem::Ls::CLinearSystem_Field& ls,
		double lambda, double myu, double rho, double g_x, double g_y, double thermoelastic, 
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, unsigned int id_field_temp, 
		bool is_inital = true,
		unsigned int id_ea = 0);

// 静的線形弾性体(剛性行列を保存)
bool AddLinSys_LinearSolid2D_Static_SaveStiffMat(
		Fem::Ls::CLinearSystem_Save& ls,
		double lambda, double myu, double rho, double g_x, double g_y,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
		unsigned int id_ea = 0 );

// 動的線形弾性体(剛性行列を保存)
bool AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta(
		Fem::Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls,
		double lambda, double myu, double rho, double g_x, double g_y,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
		unsigned int id_ea = 0 );

// 動的弾性体の固有値解析用の行列を作る
bool AddLinSys_LinearSolid2D_Eigen(
		Fem::Ls::CLinearSystem_Eigen& ls,
		double lambda, double myu, double rho,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
		unsigned int id_ea = 0 );

////////////////////////////////////////////////////////////////
// ３次元問題

// 静的線形弾性体
bool AddLinSys_LinearSolid3D_Static(
		Fem::Ls::CLinearSystem_Field& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp );

// 動的線形弾性体
bool AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta,
		Eqn::CLinearSystem_EqnInterface& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp );

// 静的線形弾性体(剛性行列を保存)
bool AddLinSys_LinearSolid3D_Static_SaveStiffMat(
		Fem::Ls::CLinearSystem_Save& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp );

// 動的弾性体の固有値解析用の行列を作る
bool AddLinSys_LinearSolid3D_Eigen(
		Fem::Ls::CLinearSystem_Eigen& ls,
		double lambda, double myu, double rho,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp );
//! @}
}
}

#endif
