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
// Eqn_LinearSolid.cpp : 線形弾性体の要素剛性作成部の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"

#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_hex.h"
#include "delfem/femeqn/eqn_linear_solid.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// ２Ｄの方程式

////////////////
// 静的問題

static bool AddLinSys_LinearSolid2D_Static_P1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea);

static bool AddLinSys_LinearSolid2D_Static_P1b(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int iei);

static bool AddLinSys_LinearSolid2D_Static_P1_SaveStiffMat(
		CLinearSystem_Save& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea);

////////////////
// 静的熱応力問題

static bool AddLinSys_LinearSolidThermalStress2D_Static_P1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double rho, double g_x, double g_y,
		double thermoelastic,
		unsigned int id_field_disp, unsigned int id_field_temp, const CFieldWorld& world,
		const unsigned int iei);

////////////////
// 動的問題

static bool AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta_P1(				
		double gamma, double beta, double dt,
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_val, const CFieldWorld& world,
		bool is_initial,
		unsigned int id_ea);

static bool AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta_P1(
		CLinearSystem_SaveDiaM_NewmarkBeta& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea);

static bool AddLinSys_LinearSolid2D_Eigen_P1(				
		CLinearSystem_Eigen& ls, 
		double lambda, double myu, double  rho,
		const unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea);

static bool AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta_P1(
		double dt, double gamma, double beta, ILinearSystem_Eqn& ls, 
		double lambda, double myu, double  rho, double g_x, double g_y, double thermoelastic, 
		const CFieldWorld& world, const unsigned int id_field_val, const unsigned int id_field_temp,
		const bool is_initial, 
		const unsigned int id_ea);


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// ３Ｄの方程式

////////////////
// 定常問題

static bool AddLinSys_LinearSolid3D_Static_P1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double rho, double g_x, double g_y, double g_z,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea);

static bool AddLinSys_LinearSolid3D_Static_Q1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double rho, double g_x, double g_y, double g_z,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea);

static bool AddLinSys_LinearSolid3D_Static_P1(
		CLinearSystem_Save& ls, 
		double lambda, double myu,
		double rho, double g_x, double g_y, double g_z,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea );

////////////////
// 非定常問題

static bool AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta_P1(				
		double gamma, double beta, double dt,
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea);

static bool AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta_Q1(				
		double gamma, double beta, double dt,
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea);

static bool AddLinSys_LinearSolid3D_Eigen_P1(
		CLinearSystem_Eigen& ls, 				
		double lambda, double myu, double  rho,
		const unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea);


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


bool Fem::Eqn::AddLinSys_LinearSolid2D_Static_SaveStiffMat(
		CLinearSystem_Save& ls,
		double lambda, double myu, double  rho, double g_x, double g_y,
		const CFieldWorld& world, unsigned int id_field_val, 
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != VECTOR2 ) return false;

	if( id_ea != 0 ){	
		if( field_val.GetInterpolationType(id_ea,world) == TRI11 ){
			return AddLinSys_LinearSolid2D_Static_P1_SaveStiffMat(
				ls,
				lambda, myu,
				rho, g_x, g_y,id_field_val,world,
				id_ea);
		}
		assert(0); 
		return false;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_LinearSolid2D_Static_SaveStiffMat(// 再帰文
				ls,
				lambda, myu, rho, g_x, g_y,
				world, id_field_val, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}


bool Fem::Eqn::AddLinSys_LinearSolid3D_Static_SaveStiffMat(
		CLinearSystem_Save& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const CFieldWorld& world,
		const unsigned int id_field_val )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != VECTOR3 ) return false;

	const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if( field_val.GetInterpolationType(id_ea,world) == TET11 ){
			AddLinSys_LinearSolid3D_Static_P1(
				ls,
				lambda, myu,
				rho, g_x, g_y, g_z,
                id_field_val,world,id_ea);
		}
		else{ assert(0); }
	}
	return true;
}


// 静的 熱応力 線形弾性体
bool Fem::Eqn::AddLinSys_LinearSolidThermalStress2D_Static(
		Eqn::ILinearSystem_Eqn& ls,
		double lambda, double myu, double rho, double g_x, double g_y, double thermoelastic, 
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, unsigned int id_field_temp, 
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_disp) ){ assert(0); return false; }
	const CField& field_val = world.GetField(id_field_disp);
	if( field_val.GetFieldType() != VECTOR2 ){ assert(0); return false; }

	if( !world.IsIdField(id_field_temp) ){ assert(0); return false; }
	const CField& field_temp = world.GetField(id_field_temp);
	if( field_temp.GetFieldType() != SCALAR ){ assert(0); return false; }

	if( id_ea != 0 ){
		bool res;
		if( field_val.GetInterpolationType(id_ea,world) == TRI11 ){
			res = AddLinSys_LinearSolidThermalStress2D_Static_P1(
				ls,
				lambda,myu,rho, g_x,g_y, thermoelastic,
                id_field_disp,id_field_temp,
				world,id_ea);
		}
		else{ 
			assert(0); 
			res = false;
		}
		assert( res );
		return res;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			// 再帰文
			bool res = Fem::Eqn::AddLinSys_LinearSolidThermalStress2D_Static(
				ls,
				lambda,myu,rho, g_x,g_y, thermoelastic, 
				world,id_field_disp,id_field_temp, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return false;
}


bool Fem::Eqn::AddLinSys_LinearSolid2D_Static(
        Eqn::ILinearSystem_Eqn& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const CFieldWorld& world, const unsigned int id_field_val, 
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ){
		assert(0);
		return false;
	}
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != VECTOR2 ){
		assert(0);
		return false;
	}

	if( id_ea != 0 ){
		if( field_val.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinSys_LinearSolid2D_Static_P1(
				ls,
				lambda, myu,
				rho, g_x, g_y,
                id_field_val,world,
                id_ea);
		}
		else if( field_val.GetInterpolationType(id_ea,world) == TRI1011 ){
			AddLinSys_LinearSolid2D_Static_P1b(
				ls,
				lambda, myu,
				rho, g_x, g_y,
                id_field_val,world,id_ea);
		}
		else{ assert(0); }
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_LinearSolid2D_Static(
				ls,
				lambda, myu, rho, g_x, g_y,
				world, id_field_val, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}


	return true;
}


// 動的弾性体の固有値解析用の行列を作る
bool Fem::Eqn::AddLinSys_LinearSolid2D_Eigen(
		Fem::Ls::CLinearSystem_Eigen& ls,
		double lambda, double myu, double rho,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_disp) ){
		assert(0);
		return false;
	}
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR2 ){
		assert(0);
		return false;
	}

	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinSys_LinearSolid2D_Eigen_P1(
				ls,
				lambda, myu, rho,
                id_field_disp,world,
				id_ea);
		}
		else{ assert(0); }
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_LinearSolid2D_Eigen(
				ls,
				lambda, myu, rho,
				world, id_field_disp, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}


	return true;
}

// 動的線形弾性体(剛性行列を保存)
bool Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta(
		Fem::Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls,
		double lambda, double myu, double rho, double g_x, double g_y,
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_disp) ){
		assert(0);
		return false;
	}
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR2 ){
		assert(0);
		return false;
	}

	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta_P1(
				ls,
				lambda, myu, rho, g_x, g_y,
                id_field_disp,world,
				id_ea);
		}
		else{ assert(0); }
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta(
				ls,
				lambda, myu, rho, g_x, g_y,
				world, id_field_disp, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}

// 動的弾性体の固有値解析用の行列を作る
bool Fem::Eqn::AddLinSys_LinearSolid3D_Eigen(
		Fem::Ls::CLinearSystem_Eigen& ls,
		double lambda, double myu, double rho,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_disp )
{
	if( !world.IsIdField(id_field_disp) ){
		assert(0);
		return false;
	}
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR3 ){
		assert(0);
		return false;
	}

	const std::vector<unsigned int>& aIdEA = field_disp.GetAry_IdElemAry();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if( field_disp.GetInterpolationType(id_ea,world) == TET11 ){
			AddLinSys_LinearSolid3D_Eigen_P1(
				ls,
				lambda, myu, rho,
                id_field_disp,world,
				id_ea);
		}
		else{ assert(0); }
	}

	return true;
}


bool Fem::Eqn::AddLinSys_LinearSolid3D_Static(
		CLinearSystem_Field& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const CFieldWorld& world,
		const unsigned int id_field_val)
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != VECTOR3 ) return false;

	const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if( field_val.GetInterpolationType(id_ea,world) == TET11 ){
			AddLinSys_LinearSolid3D_Static_P1(
				ls,
				lambda, myu,
				rho, g_x, g_y, g_z,
                id_field_val,world,id_ea);
		}
		else if( field_val.GetInterpolationType(id_ea,world) == HEX11 ){
			AddLinSys_LinearSolid3D_Static_Q1(
				ls,
				lambda, myu,
				rho, g_x, g_y, g_z,
                id_field_val,world,id_ea);
		}
		else{ assert(0); }
	}

	return true;
}

bool Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta,
		ILinearSystem_Eqn& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const CFieldWorld& world,
		unsigned int id_field_val, 
		bool is_initial,
		unsigned int id_ea )
{	
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != VECTOR2 ) return false;

	if( id_ea != 0 ){
		if( field_val.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta_P1(
				gamma, beta, dt,
				ls,
				lambda, myu, rho, g_x, g_y,
                id_field_val,world,
				is_initial,
				id_ea);
		}
		else{ assert(0); }
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta(
				dt, gamma, beta,
				ls,
				lambda, myu, rho, g_x, g_y,
				world, id_field_val, 
				is_initial,
				id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}


// 静的 熱応力 線形弾性体
bool Fem::Eqn::AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta, Fem::Ls::CLinearSystem_Field& ls,
		double lambda, double myu, double rho, double g_x, double g_y, double thermoelastic, 
		const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, unsigned int id_field_temp, 
		bool is_inital, 
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_disp) ){ assert(0); return false; }
	const CField& field_val = world.GetField(id_field_disp);
	if( field_val.GetFieldType() != VECTOR2 ){ assert(0); return false; }

	if( !world.IsIdField(id_field_temp) ){ assert(0); return false; }
	const CField& field_temp = world.GetField(id_field_temp);
	if( field_temp.GetFieldType() != SCALAR ){ assert(0); return false; }

	if( id_ea != 0 ){
		bool res;
		if( field_val.GetInterpolationType(id_ea,world) == TRI11 ){
			res = AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta_P1(
				dt,gamma,beta, ls,
				lambda,myu,rho, g_x,g_y, thermoelastic,
				world, id_field_disp,id_field_temp,
				is_inital, id_ea);
		}
		else{ 
			assert(0); 
			res = false;
		}
		assert( res );
		return res;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			// 再帰文
			bool res = Fem::Eqn::AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta(
				dt, gamma, beta,
				ls,
				lambda,myu,rho, g_x,g_y, thermoelastic, 
				world,id_field_disp,id_field_temp,
				is_inital,
				id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return false;
}

bool Fem::Eqn::AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta(
		double dt, double gamma, double beta,
		ILinearSystem_Eqn& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const CFieldWorld& world,
		unsigned int id_field_val)
{
	
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != VECTOR3 ) return false;

	const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if( field_val.GetInterpolationType(id_ea,world) == TET11 ){
			AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta_P1(
				gamma, beta, dt,
				ls,
				lambda, myu,
				rho, g_x, g_y, g_z,
                id_field_val,world,
                id_ea);
		}
		else if( field_val.GetInterpolationType(id_ea,world) == HEX11 ){
			AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta_Q1(
				gamma, beta, dt,
				ls,
				lambda, myu,
				rho, g_x, g_y, g_z,
                id_field_val,world,id_ea);
		}
		else{ assert(0); }
	}

	return true;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


static bool AddLinSys_LinearSolid2D_Static_P1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		unsigned int id_field_disp, const CFieldWorld& world,
		unsigned int id_ea)
{
//	std::cout << "LinearSolid2D Static Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va = field_disp.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;		assert( nno == es_co.GetSizeNoes() );
	const unsigned int ndim = 2;

	unsigned int noes[nno];	// 要素内の節点の節点番号

	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eres[nno][ndim];		// 要素内残差ベクトル

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_disp,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_disp = field_disp.GetNodeSeg(CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_co   = field_disp.GetNodeSeg(CORNER,false,world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{		
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);			
		// 節点の座標、値を取ってくる
		double coords[nno][ndim];		// 要素節点座標
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(noes[ino],coords[ino]);
		}
		// 要素の節点番号を取ってくる
		es_va.GetNodes(ielem,noes);			
		// 節点の座標、値を取ってくる
		double disp[nno][ndim];		// 要素節点変位	
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_disp.GetValue(noes[ino],disp[ino]);
		}

		////////////////////////////////

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);
		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項
		TriDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2]);
		// 要素剛性行列を求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno][0][0] = area*(lambda+myu)*dldx[ino][0]*dldx[jno][0];
			emat[ino][jno][0][1] = area*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
			emat[ino][jno][1][0] = area*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
			emat[ino][jno][1][1] = area*(lambda+myu)*dldx[ino][1]*dldx[jno][1];
			const double dtmp1 = (dldx[ino][1]*dldx[jno][1]+dldx[ino][0]*dldx[jno][0])*area*myu;
			emat[ino][jno][0][0] += dtmp1;
			emat[ino][jno][1][1] += dtmp1;
		}
		}
		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres[ino][0] = area*rho*g_x*0.33333333333333333;
			eres[ino][1] = area*rho*g_y*0.33333333333333333;
		}
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1];
			eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1];
		}
		}

		////////////////////////////////

		// 全体剛性行列に要素剛性行列をマージ
		mat_cc.Mearge(nno,noes, nno,noes,   ndim*ndim,   &emat[0][0][0][0]);
		// 要素内残差をマージ
		for(unsigned int ino=0;ino<nno;ino++){	
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
		}
	}

	return true;
}

static bool AddLinSys_LinearSolidThermalStress2D_Static_P1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double rho, double g_x, double g_y,
		double thermoelastic,
		unsigned int id_field_disp, unsigned int id_field_temp,
		const CFieldWorld& world,
		unsigned int id_ea)
{
//	std::cout << "LinearSolidThermalStress2D Static Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ){ assert(0); return false; }
	const CField& field_disp = world.GetField(id_field_disp);
	assert( field_disp.GetFieldType() == VECTOR2 );

	if( !world.IsIdField(id_field_temp) ){ assert(0); return false; }
	const CField& field_temp = world.GetField(id_field_temp);
	assert( field_temp.GetFieldType() == SCALAR );

	if( field_disp.GetIdElemSeg(id_ea,CORNER,true,world) == 0 ){ assert(0); return false; }
	if( field_temp.GetIdElemSeg(id_ea,CORNER,true,world) == 0 ){ assert(0); return false; }

	const CElemAry::CElemSeg& es_c_disp = field_disp.GetElemSeg(id_ea,CORNER,true,world);
	const CElemAry::CElemSeg& es_c_temp = field_temp.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;		assert( nno == es_c_disp.GetSizeNoes() );
	const unsigned int ndim = 2;

	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eres[nno][ndim];		// 要素内残差ベクトル
	
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_disp,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_disp.GetNodeSeg(CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_tmp = field_temp.GetNodeSeg(CORNER,true, world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int noes[nno];
		// 要素の節点番号を取ってくる
		es_c_temp.GetNodes(ielem,noes);
		// 節点の値を取ってくる
		double temp[nno];
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_tmp.GetValue(noes[ino], &temp[ino]);
		}
		// 要素の節点番号を取ってくる
		es_c_disp.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		double coords[nno][ndim];		// 要素節点座標
		double disp[  nno][ndim];		// 要素節点変位
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue( noes[ino],coords[ino]);
			ns_c_val.GetValue(noes[ino],  disp[ino]);
		}

		////////////////////////////////

		// 面積を求める
		const double area = TriArea(coords[0],coords[1],coords[2]);
		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項
		TriDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2]);
		// 要素剛性行列を求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno][0][0] = area*(lambda+myu)*dldx[ino][0]*dldx[jno][0];
			emat[ino][jno][0][1] = area*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
			emat[ino][jno][1][0] = area*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
			emat[ino][jno][1][1] = area*(lambda+myu)*dldx[ino][1]*dldx[jno][1];
			const double dtmp1 = (dldx[ino][1]*dldx[jno][1]+dldx[ino][0]*dldx[jno][0])*area*myu;
			emat[ino][jno][0][0] += dtmp1;
			emat[ino][jno][1][1] += dtmp1;
		}
		}
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres[ino][0] = area*rho*g_x*0.333333333333333333;
			eres[ino][1] = area*rho*g_y*0.333333333333333333;
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1];
			eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1];
		}
		}
		{	// 熱応力部分
			const double ave_temp = (temp[0]+temp[1]+temp[2])*0.33333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				eres[ino][0] += dldx[ino][0]*ave_temp*thermoelastic;
				eres[ino][1] += dldx[ino][1]*ave_temp*thermoelastic;
			}
		}

		////////////////////////////////

		// 全体剛性行列に要素剛性行列をマージ
		mat_cc.Mearge(nno,noes, nno,noes,   ndim*ndim,   &emat[0][0][0][0]);
		// 要素内残差をマージ
		for(unsigned int ino=0;ino<nno;ino++){	
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
		}
	}

	return true;
}

static bool AddLinSys_LinearSolid2D_Static_P1b(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_disp, 
		const CFieldWorld& world, 
		const unsigned int id_ea)
{
//	std::cout << "LinearSolid2D Static Tri P1b" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_val = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);
	const CElemAry::CElemSeg& es_b = field_val.GetElemSeg(id_ea,BUBBLE,true,world);

	const unsigned int nno_c = 3;	assert( nno_c == es_c.GetSizeNoes() );
	const unsigned int nno_b = 1;	assert( nno_b == es_b.GetSizeNoes() );
	const unsigned int ndim = 2;

	unsigned int noes_c[nno_c];	// 要素内の節点の節点番号
	unsigned int noes_b;

	double emat_cc[nno_c][nno_c][ndim][ndim];	// 要素剛性行列
	double emat_cb[nno_c][ndim][ndim];	// 要素剛性行列
	double emat_bc[nno_c][ndim][ndim];	// 要素剛性行列
	double emat_bb[ndim][ndim];	// 要素剛性行列
	
	double eres_c[nno_c][ndim];		// 要素内残差ベクトル
	double eres_b[ndim];		// 要素内残差ベクトル

	double coords[nno_c][ndim];		// 要素節点座標
	double disp_c[nno_c][ndim];		// 要素節点変位
	double disp_b[ndim];		// 要素節点変位

	double dldx[nno_c][ndim];		// 形状関数の空間微分
	double zero_order_term[nno_c];	// 形状関数の定数項
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(id_field_disp,CORNER,world);
	CMatDia_BlkCrs& mat_bb = ls.GetMatrix(id_field_disp,BUBBLE,world);
	CMat_BlkCrs&    mat_cb = ls.GetMatrix(id_field_disp,CORNER,id_field_disp,BUBBLE,world);
	CMat_BlkCrs&    mat_bc = ls.GetMatrix(id_field_disp,BUBBLE,id_field_disp,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_disp,CORNER,world);
	CVector_Blk&    res_b  = ls.GetResidual(id_field_disp,BUBBLE,world);

	const CNodeAry::CNodeSeg& ns_c_va = field_val.GetNodeSeg(CORNER,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_b_va = field_val.GetNodeSeg(BUBBLE,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		es_c.GetNodes(ielem,noes_c);  // 要素の節点番号を取ってくる
		for(unsigned int ino=0;ino<nno_c;ino++){ // 節点の座標、値を取ってくる
			ns_c_co.GetValue(noes_c[ino],coords[ino]);
			ns_c_va.GetValue(noes_c[ino],disp_c[ino]);
		}
		es_b.GetNodes(ielem,&noes_b);	// 要素の節点番号を取ってくる
		ns_b_va.GetValue(noes_b,disp_b);	// 節点の座標、値を取ってくる

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno_c*nno_c*ndim*ndim;i++){ *(&emat_cc[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno_c*ndim*ndim;      i++){ *(&emat_cb[0][0][0]   +i) = 0.0; }
		for(unsigned int i=0;i<nno_c*ndim*ndim;      i++){ *(&emat_bc[0][0][0]   +i) = 0.0; }
		for(unsigned int i=0;i<ndim*ndim;            i++){ *(&emat_bb[0][0]      +i) = 0.0; }

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);

		// 形状関数のｘｙ微分を求める
		TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

		double vc_b[4];
		vc_b[0] = 1.0/3.0; vc_b[1] = 1.0/3.0; vc_b[2] = 1.0/3.0; vc_b[3] = 27.0;
		{
			const double tmp_val1 = vc_b[3]*vc_b[3]*area/180.0*( 
				dldx[0][0]*dldx[0][0]+dldx[0][1]*dldx[0][1]+
				dldx[1][0]*dldx[1][0]+dldx[1][1]*dldx[1][1]+
				dldx[2][0]*dldx[2][0]+dldx[2][1]*dldx[2][1] );
            for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
            for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
				const double tmp1 
					= area*myu*(dldx[ino_c][0]*dldx[jno_c][0]+dldx[ino_c][1]*dldx[jno_c][1])
					+ myu*vc_b[ino_c]*vc_b[jno_c]*tmp_val1;
				emat_cc[ino_c][jno_c][0][0] = tmp1;
				emat_cc[ino_c][jno_c][1][1] = tmp1;
			}
			}
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				const double tmp1 = -myu*vc_b[ino_c]*tmp_val1;
				emat_cb[ino_c][0][0] = tmp1;
				emat_cb[ino_c][1][1] = tmp1;
				emat_bc[ino_c][0][0] = tmp1;
				emat_bc[ino_c][1][1] = tmp1;
			}
			emat_bb[0][0] = myu*tmp_val1;
			emat_bb[1][1] = myu*tmp_val1;
		}
		
		////////////////
		// 外力ベクトルを求める
		for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
			eres_c[ino_c][0] = g_x*area*11.0/60.0;
			eres_c[ino_c][1] = g_y*area*11.0/60.0;
		}
		eres_b[0] = g_x*area*27.0/60.0;
		eres_b[1] = g_y*area*27.0/60.0;

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno_c;ino++){
			for(unsigned int jno=0;jno<nno_c;jno++){
				eres_c[ino][0] -= emat_cc[ino][jno][0][0]*disp_c[jno][0]+emat_cc[ino][jno][0][1]*disp_c[jno][1];
				eres_c[ino][1] -= emat_cc[ino][jno][1][0]*disp_c[jno][0]+emat_cc[ino][jno][1][1]*disp_c[jno][1];
			}
			eres_c[ino][0] -= emat_cb[ino][0][0]*disp_b[0]+emat_cb[ino][0][1]*disp_b[1];
			eres_c[ino][1] -= emat_cb[ino][1][0]*disp_b[0]+emat_cb[ino][1][1]*disp_b[1];
		}

		for(unsigned int ino=0;ino<nno_c;ino++){
			eres_b[0] -= emat_bc[ino][0][0]*disp_c[ino][0]+emat_bc[ino][0][1]*disp_c[ino][1];
			eres_b[1] -= emat_bc[ino][1][0]*disp_c[ino][0]+emat_bc[ino][1][1]*disp_c[ino][1];
		}
		eres_b[0] -= emat_bb[0][0]*disp_b[0]+emat_bb[0][1]*disp_b[1];
		eres_b[1] -= emat_bb[1][0]*disp_b[0]+emat_bb[1][1]*disp_b[1];

		////////////////////////////////
		
		mat_cc.Mearge(nno_c,noes_c,  nno_c, noes_c, ndim*ndim, &emat_cc[0][0][0][0]);
		mat_cb.Mearge(nno_c,noes_c,  nno_b,&noes_b, ndim*ndim, &emat_cb[0][0][0]);
		mat_bc.Mearge(nno_b,&noes_b, nno_c, noes_c, ndim*ndim, &emat_bc[0][0][0]);
		mat_bb.Mearge(nno_b,&noes_b, nno_b,&noes_b, ndim*ndim, &emat_bb[0][0]);
		for(unsigned int ino=0;ino<nno_c;ino++){	// 要素内残差をマージ
			res_c.AddValue(noes_c[ino],0,eres_c[ino][0]);
			res_c.AddValue(noes_c[ino],1,eres_c[ino][1]);
		}
		res_b.AddValue(noes_b,0,eres_b[0]);
		res_b.AddValue(noes_b,1,eres_b[1]);
	}

	return true;
}

static bool AddLinSys_LinearSolid3D_Static_P1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea )
{
//	std::cout << "LinearSolid3D Static Tet P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 4;
	const unsigned int ndim = 3;

	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eres[nno][ndim];		// 要素内残差ベクトル

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
	    unsigned int noes[nno];	// 要素内の節点の節点番号
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
	    double coords[nno][ndim];		// 要素節点座標
	    double disp[  nno][ndim];		// 要素節点変位
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co .GetValue(noes[ino],coords[ino]);
			ns_c_val.GetValue(noes[ino],disp[  ino]);
		}

		////////////////////////////////

		// 面積を求める
		const double vol = TetVolume(coords[0],coords[1],coords[2],coords[3]);
		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項
		TetDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2],coords[3]);

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno][0][0] = vol*( lambda*dldx[ino][0]*dldx[jno][0]+myu*dldx[jno][0]*dldx[ino][0] );
			emat[ino][jno][0][1] = vol*( lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1] );
			emat[ino][jno][0][2] = vol*( lambda*dldx[ino][0]*dldx[jno][2]+myu*dldx[jno][0]*dldx[ino][2] );
			emat[ino][jno][1][0] = vol*( lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0] );
			emat[ino][jno][1][1] = vol*( lambda*dldx[ino][1]*dldx[jno][1]+myu*dldx[jno][1]*dldx[ino][1] );
			emat[ino][jno][1][2] = vol*( lambda*dldx[ino][1]*dldx[jno][2]+myu*dldx[jno][1]*dldx[ino][2] );
			emat[ino][jno][2][0] = vol*( lambda*dldx[ino][2]*dldx[jno][0]+myu*dldx[jno][2]*dldx[ino][0] );
			emat[ino][jno][2][1] = vol*( lambda*dldx[ino][2]*dldx[jno][1]+myu*dldx[jno][2]*dldx[ino][1] );
			emat[ino][jno][2][2] = vol*( lambda*dldx[ino][2]*dldx[jno][2]+myu*dldx[jno][2]*dldx[ino][2] );
            const double dtmp1 = dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2];
			emat[ino][jno][0][0] += vol*myu*dtmp1;
			emat[ino][jno][1][1] += vol*myu*dtmp1;
			emat[ino][jno][2][2] += vol*myu*dtmp1;
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres[ino][0] = vol*rho*g_x*0.25;
			eres[ino][1] = vol*rho*g_y*0.25;
			eres[ino][2] = vol*rho*g_z*0.25;
			for(unsigned int jno=0;jno<nno;jno++){
				eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
				eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
				eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
			}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nno,noes, nno,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nno;ino++){	// 要素内残差をマージ
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
			res_c.AddValue(noes[ino],2,eres[ino][2]);
		}
	}

	return true;
}

static bool AddLinSys_LinearSolid2D_Static_P1_SaveStiffMat(
		CLinearSystem_Save& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea)
{
//	std::cout << "LinearSolid2D Static Tri P1 save" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_co = field_val.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va = field_val.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	unsigned int noes[nno];	// 要素内の節点の節点番号

	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eforce[nno][ndim];		// 要素内外力ベクトル

	double coords[nno][ndim];		// 要素節点座標
	double disp[  nno][ndim];		// 要素節点変位

	CMatDia_BlkCrs& mat_cc  = ls.GetMatrix(id_field_val,CORNER,world);
	CVector_Blk&    force_c = ls.GetForce( id_field_val,CORNER,world);

	CMat_BlkCrs& mat_cc_bound = ls.GetMatrix_Boundary(id_field_val,CORNER,id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue( noes[ino],coords[ino]);
		}
		// 要素の節点番号を取ってくる
		es_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_val.GetValue(noes[ino],  disp[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<         nno*ndim;i++){ *(&eforce[0][0]    +i) = 0.0; }

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);
		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項			
		TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 emat[ino][jno][idim][jdim] 
					+= area*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  emat[ino][jno][idim][idim] += area*myu*dtmp1;
		   }
		}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eforce[ino][0] = area*rho*g_x*0.3333333333333333333;
			eforce[ino][1] = area*rho*g_y*0.3333333333333333333;
		}

		////////////////////////////////
		
		mat_cc.Mearge(      nno,noes, nno,noes,  ndim*ndim,&emat[0][0][0][0]);
		mat_cc_bound.Mearge(nno,noes, nno,noes,  ndim*ndim,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){	// 要素内残差をマージ
			force_c.AddValue(noes[ino],0,eforce[ino][0]);
			force_c.AddValue(noes[ino],1,eforce[ino][1]);
		}
	}

	return true;
}

static bool AddLinSys_LinearSolid3D_Eigen_P1(
		CLinearSystem_Eigen& ls, 
		double lambda, double myu, double  rho,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea)
{
	std::cout << "LinearSolid3D Eigen TetP1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nnoes = 4;	assert( nnoes == es_c.GetSizeNoes() );
	const unsigned int ndim = 3;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double eKmat[nnoes][nnoes][ndim][ndim];
	double eMmat[nnoes][ndim][ndim];
	double coords[nnoes][ndim];		// 要素節点座標

	CMatDia_BlkCrs& Kmat = ls.GetMatrix(       id_field_val,CORNER,world);
	CDiaMat_Blk&    Mmat = ls.GetDiaMassMatrix(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg( CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<      nnoes*ndim*ndim;i++){ *(&eMmat[0][0][0]   +i) = 0.0; }

		// 面積を求める
		const double vol = TetVolume(coords[0],coords[1],coords[2],coords[3]);
		// 形状関数のｘｙ微分を求める
		double dldx[nnoes][ndim];		// 形状関数の空間微分
		double zero_order_term[nnoes];	// 形状関数の定数項					
		TetDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2],coords[3]);

		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int jno=0;jno<nnoes;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 eKmat[ino][jno][idim][jdim] 
					+= vol*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  eKmat[ino][jno][idim][idim] += vol*myu*dtmp1;
		   }
		}
		}

		{
			const double dtmp1 = vol*rho*0.25;
			for(unsigned int ino=0;ino<nnoes;ino++){
				eMmat[ino][0][0] += dtmp1;
				eMmat[ino][1][1] += dtmp1;
				eMmat[ino][2][2] += dtmp1;
			}
		}

		////////////////////////////////
		
		Kmat.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&eKmat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){
			Mmat.Mearge(noes[ino],ndim*ndim,&eMmat[ino][0][0]);
		}
	}
	return true;
}

		
static bool AddLinSys_LinearSolid2D_Eigen_P1(
		CLinearSystem_Eigen& ls, 				
		double lambda, double myu, double  rho,
		const unsigned int id_field_disp, const CFieldWorld& world,
		unsigned int id_ea)
{
	std::cout << "LinearSolid2D Eigen P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_val = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nnoes = 3;	assert( nnoes == es_co.GetSizeNoes() );
	const unsigned int ndim = 2;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double eKmat[nnoes][nnoes][ndim][ndim];
	double eMmat[nnoes][ndim][ndim];

	double coords[nnoes][ndim];		// 要素節点座標

	CMatDia_BlkCrs& Kmat = ls.GetMatrix(       id_field_disp,CORNER,world);
	CDiaMat_Blk&    Mmat = ls.GetDiaMassMatrix(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);
		// 節点の座標を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<      nnoes*ndim*ndim;i++){ *(&eMmat[0][0][0]   +i) = 0.0; }

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);

		// 形状関数のｘｙ微分を求める
		double dldx[nnoes][ndim];		// 形状関数の空間微分
		double zero_order_term[nnoes];	// 形状関数の定数項
		TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int jno=0;jno<nnoes;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 eKmat[ino][jno][idim][jdim] 
					+= area*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  eKmat[ino][jno][idim][idim] += area*myu*dtmp1;
		   }
		}
		}

		{
			const double dtmp1 = area*rho*0.33333333333333333333;
			for(unsigned int ino=0;ino<nnoes;ino++){
				eMmat[ino][0][0] += dtmp1;
				eMmat[ino][1][1] += dtmp1;
			}
		}

		////////////////

		Kmat.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&eKmat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){
			Mmat.Mearge(noes[ino],ndim*ndim,&eMmat[ino][0][0]);
		}
	}
	return true;
}

static bool AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta_P1(
		CLinearSystem_SaveDiaM_NewmarkBeta& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea)
{
//	std::cout << "LinearSolid2D NonStatic SaveDiaM NewmarkBeta Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;	assert( nno == es_co.GetSizeNoes() );
	const unsigned int ndim = 2;

	unsigned int noes[nno];	// 要素内の節点の節点番号

	double eKmat[nno][nno][ndim][ndim];
	double eMmat[nno][ndim][ndim];
	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eqf_out[nno][ndim];	// 要素内外力ベクトル

	CMatDia_BlkCrs& mat_cc          = ls.GetMatrix(         id_field_val,CORNER,world);
	CVector_Blk&    force_c         = ls.GetForce(          id_field_val,CORNER,world);
	CDiaMat_Blk&    Mmat            = ls.GetDiaMassMatrix(  id_field_val,CORNER,world);
	CMat_BlkCrs&    mat_cc_boundary = ls.GetMatrix_Boundary(id_field_val,CORNER,id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);
		double coords[nno][ndim];		// 要素節点座標
		// 節点の座標を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<    nno*ndim*ndim;i++){ *(&eMmat[0][0][0]   +i) = 0.0; }
		for(unsigned int i=0;i<         nno*ndim;i++){ *(&eqf_out[0][0]    +i) = 0.0; }

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);

		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項
		TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 eKmat[ino][jno][idim][jdim] 
					+= area*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  eKmat[ino][jno][idim][idim] += area*myu*dtmp1;
		   }
		}
		}

		{
			const double dtmp1 = area*rho*0.33333333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				eMmat[ino][0][0] += dtmp1;
				eMmat[ino][1][1] += dtmp1;
			}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out[ino][0] = area*rho*g_x*0.33333333333333333333333333;
			eqf_out[ino][1] = area*rho*g_y*0.33333333333333333333333333;
		}

		////////////////

		const double beta = ls.GetBeta();
		const double dt = ls.GetDt();

		{
			double dtmp1 = beta*dt*dt;
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				if( ino == jno ){
					for(unsigned int i=0;i<ndim*ndim;i++){
						(&emat[ino][jno][0][0])[i] = (&eMmat[ino][0][0])[i]+dtmp1*(&eKmat[ino][jno][0][0])[i];
					}
				}
				else{
					for(unsigned int i=0;i<ndim*ndim;i++){
						(&emat[ino][jno][0][0])[i] = dtmp1*(&eKmat[ino][jno][0][0])[i];
					}
				}
			}
			}
		}
		// 全体剛性行列に要素剛性行列をマージ
		mat_cc.Mearge(nno,noes, nno,noes, ndim*ndim, &emat[0][0][0][0]);
		mat_cc_boundary.Mearge(nno,noes, nno,noes, ndim*ndim, &eKmat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			Mmat.Mearge(noes[ino],ndim*ndim,&eMmat[ino][0][0]);
		}
		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			force_c.AddValue( noes[ino],0,eqf_out[ino][idim]);
		}
		}
	}
	return true;
}


static bool AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta_P1(				
		double gamma, double beta, double dt, ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_val, const CFieldWorld& world,
		const bool is_initial, 
		const unsigned int id_ea)
{
//	std::cout << "LinearSolid2D NonStatic NewmarkBeta Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_co = field_val.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va = field_val.GetElemSeg(id_ea,CORNER,true, world);

	const unsigned int nno = 3;	assert( nno == es_co.GetSizeNoes() );
	const unsigned int ndim = 2;

	unsigned int noes[nno];	// 要素内の節点の節点番号

	double eKmat[nno][nno][ndim][ndim];
	double eMmat[nno][nno][ndim][ndim];
	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eqf_out[nno][ndim];	// 要素内外力ベクトル
	double eres[nno][ndim];		// 要素内残差ベクトル

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_val.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_acc = field_val.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);
		double coords[nno][ndim];		// 要素節点座標
		// 節点の座標を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
		}
		// 要素の節点番号を取ってくる
		es_va.GetNodes(ielem,noes);
		double disp[  nno][ndim];		// 要素節点変位
		double acc[   nno][ndim];
		double velo[  nno][ndim];
		// 節点の値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_val.GetValue( noes[ino],disp[  ino]);
			ns_c_velo.GetValue(noes[ino],velo[  ino]);
			ns_c_acc.GetValue( noes[ino],acc[   ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<         nno*ndim;i++){ *(&eqf_out[0][0]    +i) = 0.0; }

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);

		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項
		TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 eKmat[ino][jno][idim][jdim] 
					+= area*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  eKmat[ino][jno][idim][idim] += area*myu*dtmp1;
		   }
		}
		}

		{
			const double dtmp1 = area*rho*0.0833333333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat[ino][jno][0][0] += dtmp1;
					eMmat[ino][jno][1][1] += dtmp1;
				}
				eMmat[ino][ino][0][0] += dtmp1;
				eMmat[ino][ino][1][1] += dtmp1;
			}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out[ino][0] = area*rho*g_x*0.33333333333333333333333333;
			eqf_out[ino][1] = area*rho*g_y*0.33333333333333333333333333;
		}

		////////////////

		{
			double dtmp1 = beta*dt*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
			}
		}
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres[ino][0] = eqf_out[ino][0];
			eres[ino][1] = eqf_out[ino][1];
			for(unsigned int jno=0;jno<nno;jno++){
				eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1];
				eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1];
				eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1];
			}
		}
		if( is_initial ){
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eres[ino][0] -= dt*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]);
					eres[ino][1] -= dt*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]);
				}
				for(unsigned int jno=0;jno<nno;jno++){
					eres[ino][0] -= 0.5*dt*dt*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]);
					eres[ino][1] -= 0.5*dt*dt*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]);
				}
			}
		}

		////////////////////////////////

		// 全体剛性行列に要素剛性行列をマージ
		mat_cc.Mearge(nno,noes, nno,noes, ndim*ndim, &emat[0][0][0][0]);
		// 要素内残差をマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
		}
	}
	return true;
}

static bool AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta_P1(
		double dt, double gamma, double beta, ILinearSystem_Eqn& ls, 
		double lambda, double myu, double  rho, double g_x, double g_y, double thermoelastic, 
		const CFieldWorld& world, const unsigned int id_field_val, const unsigned int id_field_temp,
		const bool is_initial, 
		const unsigned int id_ea)
{
//	std::cout << "LinearSolid2D NonStatic NewmarkBeta Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_disp = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_temp) ) return false;
	const CField& field_temp = world.GetField(id_field_temp);

	const unsigned int nno = 3;	
	const unsigned int ndim = 2;

	const CElemAry::CElemSeg& es_co   = field_disp.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va   = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_temp = field_temp.GetElemSeg(id_ea,CORNER,true, world);

	assert( nno == es_co.GetSizeNoes() );

	unsigned int noes[nno];	// 要素内の節点の節点番号

	double eKmat[nno][nno][ndim][ndim];
	double eMmat[nno][nno][ndim][ndim];
	double emat[nno][nno][ndim][ndim];	// 要素剛性行列
	double eqf_out[nno][ndim];	// 要素内外力ベクトル
	double eres[nno][ndim];		// 要素内残差ベクトル

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_co   = field_disp.GetNodeSeg(CORNER,false,world);
	const CNodeAry::CNodeSeg& ns_c_val  = field_disp.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_acc  = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_temp = field_temp.GetNodeSeg(CORNER,true,world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_temp.GetNodes(ielem,noes);
		double temp[nno];		// 要素節点座標
		// 節点の座標を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_temp.GetValue(  noes[ino],&temp[ino]);
		}
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);
		double coords[nno][ndim];		// 要素節点座標
		// 節点の座標を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
		}
		// 要素の節点番号を取ってくる
		es_va.GetNodes(ielem,noes);
		double disp[  nno][ndim];		// 要素節点変位
		double acc[   nno][ndim];
		double velo[  nno][ndim];
		// 節点の値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_val.GetValue( noes[ino],disp[  ino]);
			ns_c_velo.GetValue(noes[ino],velo[  ino]);
			ns_c_acc.GetValue( noes[ino],acc[   ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<         nno*ndim;i++){ *(&eqf_out[0][0]    +i) = 0.0; }

		// 面積を求める
		double area = TriArea(coords[0],coords[1],coords[2]);

		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];		// 形状関数の空間微分
		double zero_order_term[nno];	// 形状関数の定数項
		TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 eKmat[ino][jno][idim][jdim] 
					+= area*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  eKmat[ino][jno][idim][idim] += area*myu*dtmp1;
		   }
		}
		}

		{
			const double dtmp1 = area*rho*0.0833333333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat[ino][jno][0][0] += dtmp1;
					eMmat[ino][jno][1][1] += dtmp1;
				}
				eMmat[ino][ino][0][0] += dtmp1;
				eMmat[ino][ino][1][1] += dtmp1;
			}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out[ino][0] = area*rho*g_x*0.33333333333333333333333333;
			eqf_out[ino][1] = area*rho*g_y*0.33333333333333333333333333;
		}

		////////////////

		{
			double dtmp1 = beta*dt*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
			}
		}
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres[ino][0] = eqf_out[ino][0];
			eres[ino][1] = eqf_out[ino][1];
			for(unsigned int jno=0;jno<nno;jno++){
				eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1];
				eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1];
				eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1];
			}
		}
		{	// 熱応力部分
			const double ave_temp = (temp[0]+temp[1]+temp[2])*0.33333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				eres[ino][0] += dldx[ino][0]*ave_temp*thermoelastic;
				eres[ino][1] += dldx[ino][1]*ave_temp*thermoelastic;
			}
		}
		if( is_initial ){
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eres[ino][0] -= dt*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]);
					eres[ino][1] -= dt*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]);
				}
				for(unsigned int jno=0;jno<nno;jno++){
					eres[ino][0] -= 0.5*dt*dt*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]);
					eres[ino][1] -= 0.5*dt*dt*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]);
				}
			}
		}

		////////////////////////////////

		// 全体剛性行列に要素剛性行列をマージ
		mat_cc.Mearge(nno,noes, nno,noes, ndim*ndim, &emat[0][0][0][0]);
		// 要素内残差をマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
		}
	}
	return true;
}



static bool AddLinSys_LinearSolid3D_Static_P1(
		CLinearSystem_Save& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea )
{
//	std::cout << "LinearSolid2D Tetrahedra 4-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nnoes = 4;
	const unsigned int ndim = 3;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eforce[nnoes][ndim];		// 要素内外力ベクトル

	double coords[nnoes][ndim];		// 要素節点座標
	double disp[  nnoes][ndim];		// 要素節点変位

	double dldx[nnoes][ndim];		// 形状関数の空間微分
	double zero_order_term[nnoes];	// 形状関数の定数項
				
	CMatDia_BlkCrs& mat_cc       = ls.GetMatrix(id_field_val,CORNER,world);
	CVector_Blk&    force_c      = ls.GetForce( id_field_val,CORNER,world);
	CMat_BlkCrs&    mat_cc_bound = ls.GetMatrix_Boundary(id_field_val,CORNER,id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(noes[ino],coords[ino]);
			ns_c_val.GetValue(noes[ino],disp[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce[0][0]    +i) = 0.0; }

		// 面積を求める
		const double vol = TetVolume(coords[0],coords[1],coords[2],coords[3]);

		// 形状関数のｘｙ微分を求める
		TetDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2],coords[3]);

		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int jno=0;jno<nnoes;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
		      for(unsigned int jdim=0;jdim<ndim;jdim++){
		         emat[ino][jno][idim][jdim]
				     += vol*( lambda*dldx[ino][idim]*dldx[jno][jdim]
				             +myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  emat[ino][jno][idim][idim] += vol*myu*dtmp1;
		   }
		}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
			eforce[ino][0] = vol*rho*g_x*0.25;
			eforce[ino][1] = vol*rho*g_y*0.25;
			eforce[ino][2] = vol*rho*g_z*0.25;
		}

		////////////////////////////////
		
		mat_cc      .Mearge(nnoes,noes, nnoes,noes,  ndim*ndim,&emat[0][0][0][0]);
		mat_cc_bound.Mearge(nnoes,noes, nnoes,noes,  ndim*ndim,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
			force_c.AddValue(noes[ino],0,eforce[ino][0]);
			force_c.AddValue(noes[ino],1,eforce[ino][1]);
			force_c.AddValue(noes[ino],2,eforce[ino][2]);
		}
	}

	return true;
}



static bool AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta_P1(
		double gamma, double beta, double dt,
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea)
{
//	std::cout << "LinearSolid NonStatic NewmarkBeta 3D Tetrahedra 4-point 1st order" << std::endl;
//	std::cout << gamma << " " << beta << " " << dt << " " << lambda << " " << myu << " " << rho << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nnoes = 4;	assert( nnoes == es_c.GetSizeNoes() );
	const unsigned int ndim = 3;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double eKmat[nnoes][nnoes][ndim][ndim];
	double eMmat[nnoes][nnoes][ndim][ndim];
	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eqf_out[nnoes][ndim];	// 要素内外力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	double coords[nnoes][ndim];		// 要素節点座標
	double disp[  nnoes][ndim];		// 要素節点変位
	double acc[   nnoes][ndim];
	double velo[  nnoes][ndim];

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg( CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_val.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_acc = field_val.GetNodeSeg( CORNER,true, world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg( CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
			ns_c_val.GetValue( noes[ino],disp[  ino]);
			ns_c_velo.GetValue(noes[ino],velo[  ino]);
			ns_c_acc.GetValue( noes[ino],acc[   ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eqf_out[0][0]    +i) = 0.0; }

		// 面積を求める
		const double vol = TetVolume(coords[0],coords[1],coords[2],coords[3]);
		// 形状関数のｘｙ微分を求める
		double dldx[nnoes][ndim];		// 形状関数の空間微分
		double zero_order_term[nnoes];	// 形状関数の定数項					
		TetDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2],coords[3]);

		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int jno=0;jno<nnoes;jno++){
		   double dtmp1 = 0.0;
		   for(unsigned int idim=0;idim<ndim;idim++){
			  for(unsigned int jdim=0;jdim<ndim;jdim++){
				 eKmat[ino][jno][idim][jdim] 
					+= vol*( lambda*dldx[ino][idim]*dldx[jno][jdim]+myu*dldx[jno][idim]*dldx[ino][jdim] );
			  }
			  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
		   }
		   for(unsigned int idim=0;idim<ndim;idim++){
			  eKmat[ino][jno][idim][idim] += vol*myu*dtmp1;
		   }
		}
		}

		{
			const double dtmp1 = vol*rho*0.05;
			for(unsigned int ino=0;ino<nnoes;ino++){
				for(unsigned int jno=0;jno<nnoes;jno++){
					eMmat[ino][jno][0][0] += dtmp1;
					eMmat[ino][jno][1][1] += dtmp1;
					eMmat[ino][jno][2][2] += dtmp1;
				}
				eMmat[ino][ino][0][0] += dtmp1;
				eMmat[ino][ino][1][1] += dtmp1;
				eMmat[ino][ino][2][2] += dtmp1;
			}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
			eqf_out[ino][0] = vol*rho*g_x*0.25;
			eqf_out[ino][1] = vol*rho*g_y*0.25;
			eqf_out[ino][2] = vol*rho*g_z*0.25;
		}

		////////////////

		{
			double dtmp1 = beta*dt*dt;
			for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){
				(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i] + dtmp1*(&eKmat[0][0][0][0])[i];
			}
		}
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
			eres[ino][0] = eqf_out[ino][0];
			eres[ino][1] = eqf_out[ino][1];
			eres[ino][2] = eqf_out[ino][2];
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= 
					 eKmat[ino][jno][idim][0]*disp[jno][0]
					+eKmat[ino][jno][idim][1]*disp[jno][1]
					+eKmat[ino][jno][idim][2]*disp[jno][2];
			}
			}
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= 
					 eMmat[ino][jno][idim][0]*acc[jno][0]
					+eMmat[ino][jno][idim][1]*acc[jno][1]
					+eMmat[ino][jno][idim][2]*acc[jno][2];
			}
			}
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= dt*(
					 eKmat[ino][jno][idim][0]*velo[jno][0]
					+eKmat[ino][jno][idim][1]*velo[jno][1]
					+eKmat[ino][jno][idim][2]*velo[jno][2]);
			}
			}
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= 0.5*dt*dt*(
					 eKmat[ino][jno][idim][0]*acc[jno][0]
					+eKmat[ino][jno][idim][1]*acc[jno][1]
					+eKmat[ino][jno][idim][2]*acc[jno][2]);
			}
			}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
			res_c.AddValue(noes[ino],2,eres[ino][2]);
		}
	}
	return true;
}

static bool AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta_Q1(
		double gamma, double beta, double dt,
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea )
{
	std::cout << "LinearSolid NonStatic NewmarkBeta 3D Hexahedra 8-point 1st order" << std::endl;
	std::cout << gamma << " " << beta << " " << dt << " " << lambda << " " << myu << " " << rho << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == HEX );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);

	unsigned int num_integral = 1;
	const unsigned int nInt = NIntLineGauss[num_integral];
	const double (*Gauss)[2] = LineGauss[num_integral];
	double detjac, detwei;

	const unsigned int nnoes = 8;	assert( nnoes == es_c.GetSizeNoes() );
	const unsigned int ndim = 3;

	double eKmat[nnoes][nnoes][ndim][ndim];
	double eMmat[nnoes][nnoes][ndim][ndim];
	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eqf_out[nnoes][ndim];	// 要素内外力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	double dndx[nnoes][ndim];		// 形状関数の空間微分
	double an[nnoes];				// 形状関数の値
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg( CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_val.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_acc = field_val.GetNodeSeg( CORNER,true, world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg( CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
	    unsigned int noes[nnoes];	// 要素内の節点の節点番号
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
	    double coords[nnoes][ndim];		// 要素節点座標
	    double disp[  nnoes][ndim];		// 要素節点変位
	    double acc[   nnoes][ndim];
	    double velo[  nnoes][ndim];
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(  noes[ino],coords[ino]);
			ns_c_val.GetValue( noes[ino],disp[  ino]);
			ns_c_velo.GetValue(noes[ino],velo[  ino]);
			ns_c_acc.GetValue( noes[ino],acc[   ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eqf_out[0][0]    +i) = 0.0; }

		double vol = 0.0;
		for(unsigned int ir1=0;ir1<nInt;ir1++){
		for(unsigned int ir2=0;ir2<nInt;ir2++){
		for(unsigned int ir3=0;ir3<nInt;ir3++){
			const double r1 = Gauss[ir1][0];
			const double r2 = Gauss[ir2][0];
			const double r3 = Gauss[ir3][0];
			ShapeFunc_Hex8(r1,r2,r3,coords,detjac,dndx,an);
			detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
			vol += detwei;
            for(unsigned int ino=0;ino<nnoes;ino++){
            for(unsigned int jno=0;jno<nnoes;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
					for(unsigned int jdim=0;jdim<ndim;jdim++){
						eKmat[ino][jno][idim][jdim]
							+= detwei*( lambda*dndx[ino][idim]*dndx[jno][jdim]
							              +myu*dndx[jno][idim]*dndx[ino][jdim] );
					}
					dtmp1 += dndx[ino][idim]*dndx[jno][idim];
				}
				for(unsigned int idim=0;idim<ndim;idim++){
					eKmat[ino][jno][idim][idim] += detwei*myu*dtmp1;
				}
            }
            }
            for(unsigned int ino=0;ino<nnoes;ino++){
            for(unsigned int jno=0;jno<nnoes;jno++){
				eMmat[ino][jno][0][0] += rho*detwei*an[ino]*an[jno];
				eMmat[ino][jno][1][1] += rho*detwei*an[ino]*an[jno];
				eMmat[ino][jno][2][2] += rho*detwei*an[ino]*an[jno];
			}
			}
			// 要素節点等価外力ベクトルを積n分毎に足し合わせる
			for(unsigned int ino=0;ino<nnoes;ino++){
				eqf_out[ino][0] += detwei*rho*g_x*an[ino];
				eqf_out[ino][1] += detwei*rho*g_y*an[ino];
				eqf_out[ino][2] += detwei*rho*g_z*an[ino];
			}
		}
		}
		}

		////////////////

		{
			double dtmp1 = beta*dt*dt;
			for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){
				(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i] + dtmp1*(&eKmat[0][0][0][0])[i];
			}
		}
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
			eres[ino][0] = eqf_out[ino][0];
			eres[ino][1] = eqf_out[ino][1];
			eres[ino][2] = eqf_out[ino][2];
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= 
					 eKmat[ino][jno][idim][0]*disp[jno][0]
					+eKmat[ino][jno][idim][1]*disp[jno][1]
					+eKmat[ino][jno][idim][2]*disp[jno][2];
			}
			}
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= 
					 eMmat[ino][jno][idim][0]*acc[jno][0]
					+eMmat[ino][jno][idim][1]*acc[jno][1]
					+eMmat[ino][jno][idim][2]*acc[jno][2];
			}
			}
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= dt*(
					 eKmat[ino][jno][idim][0]*velo[jno][0]
					+eKmat[ino][jno][idim][1]*velo[jno][1]
					+eKmat[ino][jno][idim][2]*velo[jno][2]);
			}
			}
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres[ino][idim] -= 0.5*dt*dt*(
					 eKmat[ino][jno][idim][0]*acc[jno][0]
					+eKmat[ino][jno][idim][1]*acc[jno][1]
					+eKmat[ino][jno][idim][2]*acc[jno][2]);
			}
			}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
			res_c.AddValue(noes[ino],2,eres[ino][2]);
		}
	}
	return true;
}


static bool AddLinSys_LinearSolid3D_Static_Q1(
		ILinearSystem_Eqn& ls, 
		double lambda, double myu,
		double  rho, double g_x, double g_y, double g_z,
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea)
{
//	std::cout << "LinearSolid3D Static Hex Q1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == HEX );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,false,world);

	unsigned int num_integral = 1;
	const unsigned int nInt = NIntLineGauss[num_integral];
	const double (*Gauss)[2] = LineGauss[num_integral];
	double detjac, detwei;

	const unsigned int nnoes = 8;
	const unsigned int ndim = 3;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eforce[nnoes][ndim];		// 要素内外力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	double coords[nnoes][ndim];		// 要素節点座標
	double disp[  nnoes][ndim];		// 要素節点変位

	double dndx[nnoes][ndim];		// 形状関数の空間微分
	double an[nnoes];				// 形状関数の値
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(noes[ino],coords[ino]);
			ns_c_val.GetValue(noes[ino],disp[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce[0][0]    +i) = 0.0; }

		double vol = 0.0;
		for(unsigned int ir1=0;ir1<nInt;ir1++){
		for(unsigned int ir2=0;ir2<nInt;ir2++){
		for(unsigned int ir3=0;ir3<nInt;ir3++){
			const double r1 = Gauss[ir1][0];
			const double r2 = Gauss[ir2][0];
			const double r3 = Gauss[ir3][0];
			ShapeFunc_Hex8(r1,r2,r3,coords,detjac,dndx,an);
			detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
			vol += detwei;
            for(unsigned int ino=0;ino<nnoes;ino++){
            for(unsigned int jno=0;jno<nnoes;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
					for(unsigned int jdim=0;jdim<ndim;jdim++){
						emat[ino][jno][idim][jdim]
							+= detwei*( lambda*dndx[ino][idim]*dndx[jno][jdim]
							              +myu*dndx[jno][idim]*dndx[ino][jdim] );
					}
					dtmp1 += dndx[ino][idim]*dndx[jno][idim];
				}
				for(unsigned int idim=0;idim<ndim;idim++){
					emat[ino][jno][idim][idim] += detwei*myu*dtmp1;
				}
			}
			}
			// 要素節点等価外力ベクトルを積n分毎に足し合わせる
			for(unsigned int ino=0;ino<nnoes;ino++){
				eforce[ino][0] += detwei*rho*g_x*an[ino];
				eforce[ino][1] += detwei*rho*g_y*an[ino];
				eforce[ino][2] += detwei*rho*g_z*an[ino];
			}
		}
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
			eres[ino][0] = eforce[ino][0];
			eres[ino][1] = eforce[ino][1];
			eres[ino][2] = eforce[ino][2];
			for(unsigned int jno=0;jno<nnoes;jno++){
				eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
				eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
				eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
			}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
			res_c.AddValue(noes[ino],0,eres[ino][0]);
			res_c.AddValue(noes[ino],1,eres[ino][1]);
			res_c.AddValue(noes[ino],2,eres[ino][2]);
		}
	}

	return true;
}


