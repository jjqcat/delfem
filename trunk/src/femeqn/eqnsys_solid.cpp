
#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <math.h>

#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/eqn_linear_solid.h"
#include "delfem/femeqn/eqn_st_venant.h"

#include "delfem/eqnsys_solid.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// �R�c�̕�����

CEqn_Solid3D_Linear::CEqn_Solid3D_Linear(unsigned int id_field, Fem::Field::CFieldWorld& world) 
: m_IsGeomNonlin(false), m_IsSaveStiffMat(false), m_IsStationary(false)
{
	m_lambda = 0.0;
	m_myu = 0.0;
	m_rho = 1.0;
	
	m_g_x = 0.0;
	m_g_y = 0.0;
	m_g_z = 0.0;

	this->SetDomain_Field(id_field, world);
}

CEqn_Solid3D_Linear::CEqn_Solid3D_Linear()
: m_IsGeomNonlin(false), m_IsSaveStiffMat(false), m_IsStationary(false)
{
	m_lambda = 0.0;
	m_myu = 0.0;
	m_rho = 1.0;
	
	m_g_x = 0.0;
	m_g_y = 0.0;
	m_g_z = 0.0;
}

double CEqn_Solid3D_Linear::MakeLinearSystem(const Fem::Field::CFieldWorld& world, bool is_initial)
{	
	if( pLS==0 || pPrec==0 ) this->InitializeLinearSystem(world);

	// �A���ꎟ�����������
	pLS->InitializeMarge();	// �A���ꎟ������������������
	if( this->m_IsGeomNonlin ){
		assert( !this->m_IsSaveStiffMat );
		if( this->m_IsStationary ){
			// �ÓISt.Venant-Kichhoff��
			Fem::Eqn::AddLinSys_StVenant3D_Static(
				*pLS,
				m_lambda, m_myu, m_rho,   m_g_x, m_g_y, m_g_z,
				world, m_IdFieldDisp);
		}
		else{ // ���ISt.Vennat-Kirchhoff��
			Fem::Eqn::AddLinSys_StVenant3D_NonStatic_NewmarkBeta(
				m_dt, m_gamma_newmark, m_beta_newmark,
				*pLS,
				m_lambda, m_myu, m_rho,   m_g_x, m_g_y, m_g_z,
				world, m_IdFieldDisp,
				is_initial);
		}
	}
	else{
		if( this->m_IsStationary ){
			if( this->m_IsSaveStiffMat ){		
				// �ÓI���`�e����(�����s��ۑ�)
				Fem::Eqn::AddLinSys_LinearSolid3D_Static_SaveStiffMat(
					*((CLinearSystem_Save*)pLS),
					m_lambda, m_myu, m_rho,   m_g_x, m_g_y, m_g_z,
					world,
					m_IdFieldDisp);
			}
			else{ // �ÓI���`�e����
				Fem::Eqn::AddLinSys_LinearSolid3D_Static(
					*pLS,
					m_lambda, m_myu, m_rho,   m_g_x, m_g_y, m_g_z,
					world, 
					m_IdFieldDisp);
			}
		}
		else{ // ���I���`�e����
			assert( !this->m_IsSaveStiffMat );
			Fem::Eqn::AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta(
				m_dt, m_gamma_newmark, m_beta_newmark,
				*pLS,
				m_lambda, m_myu, m_rho,   m_g_x, m_g_y, m_g_z,
				world,
				m_IdFieldDisp);
		}
	}
	const double norm_res = pLS->FinalizeMarge();

	// �O�����s������
	pPrec->SetValue( (*pLS).m_ls );

	return norm_res;
}

bool CEqn_Solid3D_Linear::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{
	if( pLS!=0 || pPrec!=0 ) ClearLinearSystemPreconditioner();

	// �A���ꎟ�������N���X�̍쐬
	if( this->m_IsGeomNonlin ){ assert( !this->m_IsSaveStiffMat ); }
	 
	if( this->m_IsSaveStiffMat ){ pLS = new CLinearSystem_Save; }
	else{ pLS = new CLinearSystem_Field; }

	// �A���ꎟ�������N���X�̐ݒ�
	pLS->AddPattern_Field(m_IdFieldDisp,world);	// val_field����ł���S�̍����s���ǉ�����
	for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
		const unsigned int id_field = m_aIdFixField[idf].first; 
		const int idof = m_aIdFixField[idf].second;
		if( idof == -1 ){
			pLS->SetFixedBoundaryCondition_Field( id_field, 0, world ); // bc0_field���Œ苫�E�����ɐݒ�	
			pLS->SetFixedBoundaryCondition_Field( id_field, 1, world ); // bc0_field���Œ苫�E�����ɐݒ�	
			pLS->SetFixedBoundaryCondition_Field( id_field, 2, world ); // bc0_field���Œ苫�E�����ɐݒ�
		}
		else{
			pLS->SetFixedBoundaryCondition_Field( id_field, idof, world ); // bc0_field���Œ苫�E�����ɐݒ�	
		}
	}
		
	// �O�����N���X�̍쐬
	assert( pPrec == 0 );
    pPrec = new LsSol::CPreconditioner_ILU( (*pLS).m_ls, 1 );

	return true;
}

bool CEqn_Solid3D_Linear::Solve(Fem::Field::CFieldWorld& world)
{
	this->m_aItrNormRes.clear();
	if( this->m_IsGeomNonlin ){
		assert( !this->m_IsSaveStiffMat );
		if( pLS == 0 || pPrec == 0 ){ this->InitializeLinearSystem(world); }
        double ini_norm_res = 0;
		for(unsigned int iitr=0;iitr<40;iitr++){
//			std::cout << iitr << std::endl;
			double norm_res = this->MakeLinearSystem(world,iitr==0);
			if( iitr==0 ){
				if( norm_res < 1.0e-20 ) break;	// �����c������������Δ�����
				ini_norm_res = norm_res;
			}
			if( norm_res < ini_norm_res * 1.0e-6 ) break;
//			std::cout << iitr << " " << norm_res << " " << norm_res / ini_norm_res << std::endl;
			{	// �s�������
				double conv_ratio = 1.0e-6;
				unsigned int max_iter = 1000;
				// Solve with Preconditioned Conjugate Gradient
                LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
				Sol::Solve_PCG(conv_ratio,max_iter,lsp);
				// Solve with Conjugate Gradient
			//	Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);
//				std::cout << max_iter << " " << conv_ratio << std::endl;
				this->m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
			}
			if( this->m_IsStationary ){
				pLS->UpdateValueOfField(m_IdFieldDisp,world,VALUE);
			}
			else{
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
					m_IdFieldDisp,world, iitr==0 );
			}
		}
	}
	else{ 
		if( pLS == 0 || pPrec == 0 ){
			this->InitializeLinearSystem(world);
			this->MakeLinearSystem(world,true);
		}
		else{ if( !this->m_IsSaveStiffMat ){ this->MakeLinearSystem(world,true); } }

		if( this->m_IsSaveStiffMat ){ pLS->MakeResidual(world); }

		{	// �s�������
			double conv_ratio = 1.0e-6;
			unsigned int max_iter = 1000;
			// Solve with Preconditioned Conjugate Gradient
            LsSol::CLinearSystemPreconditioner lsp( (*pLS).m_ls, *pPrec );
			Sol::Solve_PCG(conv_ratio,max_iter,lsp);	
			// Solve with Conjugate Gradient
		//	Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);
//			std::cout << max_iter << " " << conv_ratio << std::endl;
			this->m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
		}
		// �����X�V����
		if( this->m_IsStationary ){
			if( this->m_IsSaveStiffMat ){ 
				((Fem::Ls::CLinearSystem_Save*)pLS)->UpdateValueOfField(m_IdFieldDisp,world,VALUE); 
			}
			else{
				pLS->UpdateValueOfField(m_IdFieldDisp,world,VALUE);
			}
		}
		else{
			assert( !this->m_IsSaveStiffMat );
			pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
				m_IdFieldDisp,world,true);
		}
	}
	return true;
}

bool CEqn_Solid3D_Linear::SetDomain_Field(unsigned int id_field_base, Fem::Field::CFieldWorld& world){
	{	// ���̓t�B�[���h�̍��W�ߓ_�Z�O�����g��dof��3���ǂ����`�F�b�N����
//		unsigned int id_field_base = world.GetFieldBaseID();
		assert( world.IsIdField(id_field_base) );
		const CField& field_base = world.GetField(id_field_base);
		assert( field_base.GetNDimCoord() == 3 );
	}
	m_IdFieldDisp = world.MakeField_FieldElemDim(id_field_base, 3, 
        VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
	this->ClearLinearSystemPreconditioner();
	return true;
}


bool CEqn_Solid3D_Linear::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return true;
}

unsigned int CEqn_Solid3D_Linear::AddFixElemAry( 
		unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world, idof );
}

bool CEqn_Solid3D_Linear::ClearFixElemAry(
		unsigned int id_ea, Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA( id_ea ) ) return false;
	for(unsigned int ifix=0;ifix<m_aIdFixField.size();ifix++){
		const unsigned int id_field_fix = m_aIdFixField[ifix].first;
		const Fem::Field::CField& field = world.GetField(id_field_fix);
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		if( aIdEA.size() != 1 ){
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
		}
		if( aIdEA[0] == id_ea ){
			m_aIdFixField.erase( m_aIdFixField.begin()+ifix );
			this->ClearLinearSystem();
		}
	}
	return true;
}


void CEqn_Solid3D_Linear::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}

unsigned int CEqn_Solid3D_Linear::AddFixElemAry( 
		const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
	const unsigned int id_field = world.GetPartialField(m_IdFieldDisp, aIdEA);
	if( id_field == 0 ) return 0;
	assert( world.IsIdField(id_field) );
	{
		CField& field = world.GetField(id_field);
		unsigned int nlen_val = field.GetNLenValue();
		for(unsigned int ilen=0;ilen<nlen_val;ilen++){
			field.SetValue(       0.0,ilen,world,false);
			field.SetVelocity(    0.0,ilen,world,false);
			field.SetAcceleration(0.0,ilen,world,false);
		}
	}
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return id_field;
}




////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// �Q�c�̕�����




bool CEqn_Solid2D::AddLinSys(Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world )
{
	if( this->m_IsGeomNonlin ){
		// �ÓISt.Venant-Kirchhoff��
		return Fem::Eqn::AddLinSys_StVenant2D_Static(
			ls,
			m_lambda, m_myu, m_rho,   m_g_x, m_g_y,
			world,m_IdFieldDisp,
			m_id_ea );
	}
	else{
		if( !world.IsIdField(m_IdFieldTemperature) ){
			// �ÓI���`�e����
			return Fem::Eqn::AddLinSys_LinearSolid2D_Static(
				ls,
				m_lambda, m_myu, m_rho,   m_g_x, m_g_y,
				world, m_IdFieldDisp,
				m_id_ea);
		}
		else{ 
			// �ÓI���`�e����(�M����)
			return Fem::Eqn::AddLinSys_LinearSolidThermalStress2D_Static(
				ls,
				m_lambda, m_myu, m_rho,   m_g_x, m_g_y, 0.02,
				world, m_IdFieldDisp, m_IdFieldTemperature,
				m_id_ea);
		}
	}
	return false;
}

// �A���ꎟ�������}�[�W���\�b�h
bool CEqn_Solid2D::AddLinSys_NewmarkBetaAPrime( double dt, double gamma, double beta, bool is_initial, 
	Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world )
{
	if( this->m_IsGeomNonlin ){
		// ���ISt.Venant-Kirchhoff��
		return Fem::Eqn::AddLinSys_StVenant2D_NonStatic_NewmarkBeta(
			dt, gamma, beta, ls,
			m_lambda, m_myu, m_rho,   m_g_x, m_g_y,
			world, m_IdFieldDisp,
			is_initial,
			m_id_ea);
	}
	else{
		if( !world.IsIdField(m_IdFieldTemperature) ){
			// ���I���`�e����
			return Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta(
				dt, gamma, beta, ls,
				m_lambda, m_myu, m_rho,   m_g_x, m_g_y,
				world, m_IdFieldDisp,
				is_initial,
				m_id_ea);
		}
		else{
			// ���I���`�e����
			return Fem::Eqn::AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta(
				dt, gamma, beta, ls,
				m_lambda, m_myu, m_rho, m_g_x, m_g_y, 0.02,
				world, m_IdFieldDisp, m_IdFieldTemperature,
				is_initial,
				m_id_ea);
		}
	}

	return false;
}

// �A���ꎟ�������}�[�W���\�b�h
bool CEqn_Solid2D::AddLinSys_NewmarkBetaAPrime_Save(
	Fem::Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls, const Fem::Field::CFieldWorld& world )
{
	assert( !this->m_IsGeomNonlin );
	assert( !world.IsIdField(m_IdFieldTemperature) );
	// ���I���`�e����
	return Fem::Eqn::AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta(
		ls,
		m_lambda, m_myu, m_rho,   m_g_x, m_g_y,
		world, m_IdFieldDisp,
		m_id_ea);

	return false;
}

bool CEqn_Solid2D::AddLinSys_Save( Fem::Ls::CLinearSystem_Save& ls, const Fem::Field::CFieldWorld& world )
{
	return Fem::Eqn::AddLinSys_LinearSolid2D_Static_SaveStiffMat(
		ls,
		m_lambda, m_myu, m_rho,   m_g_x, m_g_y,
		world, m_IdFieldDisp);
	return false;
}


CEqnSystem_Solid2D::CEqnSystem_Solid2D(unsigned int id_field, Fem::Field::CFieldWorld& world) 
: m_IsSaveStiffMat(false)
{
	m_rho_back = 1;
	m_young_back = 1;
	m_poisson_back = 0;
	m_is_plane_stress_back = true;
	this->UpdateDomain_Field(id_field, world);
}

CEqnSystem_Solid2D::CEqnSystem_Solid2D()
: m_IsSaveStiffMat(false)
{
	m_rho_back = 1;
	m_young_back = 1;
	m_poisson_back = 0;
	m_is_plane_stress_back = true;
}

bool CEqnSystem_Solid2D::EqnationProperty(bool& is_nonlin){
	is_nonlin = false;
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if(  m_aEqn[ieqn].IsGeometricalNonlinear() )	is_nonlin = true;
	}
	return true;
}

double CEqnSystem_Solid2D::MakeLinearSystem(const Fem::Field::CFieldWorld& world, bool is_initial)
{	
//	std::cout << "CEqnSystem_Solid2D::MakeLinearSystem" << std::endl;
	if( pLS==0 ){ this->InitializeLinearSystem(world); }
	// �A���ꎟ�����������
	assert( pLS != 0 );
	pLS->InitializeMarge();	// �A���ꎟ������������������
	{	// �A���������Ƀ}�[�W
		bool is_nonlin;
		this->EqnationProperty(is_nonlin); 
		if( this->m_IsSaveStiffMat ){
//			std::cout << "Save Stiff Mat " << std::endl;
			assert( !is_nonlin ); // �s��ۑ����ł���̂͐��`�������̂Ƃ��̂�
			if( m_IsStationary ){
				for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
					m_aEqn[ieqn].AddLinSys_Save( *(CLinearSystem_Save*)pLS, world);
				}
			}
			else{
				((CLinearSystem_SaveDiaM_NewmarkBeta*)pLS)->SetNewmarkParameter(m_beta_newmark,m_gamma_newmark,m_dt);
				for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
					m_aEqn[ieqn].AddLinSys_NewmarkBetaAPrime_Save( *(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS, world);
				}
			}
		}
		else{
//			std::cout << "NotSave Stiff Mat " << std::endl;
			if( m_IsStationary ){
//				std::cout << "Add LinSys " << std::endl;
				for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
//					std::cout << "Add LinSys " << ieqn << " " << m_aEqn[ieqn].GetIdField_Disp() << std::endl;
					m_aEqn[ieqn].AddLinSys(*pLS, world);
				}
			}
			else{
				for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
					m_aEqn[ieqn].AddLinSys_NewmarkBetaAPrime(m_dt,m_gamma_newmark,m_beta_newmark,is_initial,  
						*pLS, world);
				}
			}
		}
	}
/*
	{	// �׏d���E������g�ݓ����
		const Fem::Field::CField& field_disp = world.GetField(m_id_val);
		unsigned int id_na_val = field_disp.GetNodeSegInNodeAry(CORNER).id_na_va;
		unsigned int id_na_co = field_disp.GetNodeSegInNodeAry(CORNER).id_na_co;
		MatVec::CVector_Blk* res = pLS->GetResidualPtr(m_id_val,CORNER,world);
		for(unsigned int iload=0;iload<m_aLoad.size();iload++){
			unsigned int id_ea = m_aLoad[iload].first;
			double load = m_aLoad[iload].second;
			const CElemAry& ea = world.GetEA(id_ea);
			std::vector<unsigned int> aIdES = ea.GetAry_SegID();
			unsigned int id_es_va = 0;
			unsigned int id_es_co = 0;
			for(unsigned int iies=0;iies<aIdES.size();iies++){	// id_na_val������id_es������
				unsigned int id_es0 = aIdES[iies];
				if( ea.GetSeg(id_es0).GetIdNA() == id_na_val ){ 
					id_es_va = id_es0;
				}
				if( ea.GetSeg(id_es0).GetIdNA() == id_na_co ){ 
					id_es_co = id_es0;
				}
			}
			if( id_es_va == 0 ) continue;
			double tot_len = 0;
			{
				const CElemAry::CElemSeg& es_co = ea.GetSeg(id_es_co);
				unsigned int nnoes = es_co.GetSizeNoes();
				if( nnoes != 2 ) continue;	// �ӗv�f�łȂ���΂Ȃ�Ȃ�
				unsigned int noes_co[2];
				const CNodeAry::CNodeSeg& ns_co = field_disp.GetNodeSeg(CORNER,false,world);
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es_co.GetNodes(ielem,noes_co);
					const unsigned int inode0 = noes_co[0];
					const unsigned int inode1 = noes_co[1];
					double coord0[2], coord1[2];
					ns_co.GetValue(inode0,coord0);
					ns_co.GetValue(inode1,coord1);
					double len = sqrt( (coord1[0]-coord0[0])*(coord1[0]-coord0[0])+(coord1[1]-coord0[1])*(coord1[1]-coord0[1]) );
					tot_len += len;
				}
			}
			{
				const CElemAry::CElemSeg& es_va = ea.GetSeg(id_es_va);
				const CElemAry::CElemSeg& es_co = ea.GetSeg(id_es_va);
				unsigned int nnoes = es_va.GetSizeNoes();
				if( nnoes != 2 ) continue;	// �ӗv�f�łȂ���΂Ȃ�Ȃ�
				unsigned int noes_co[2];
				const CNodeAry::CNodeSeg& ns_co = field_disp.GetNodeSeg(CORNER,false,world);
				unsigned int noes_va[2];
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es_co.GetNodes(ielem,noes_co);
					const unsigned int inode0 = noes_co[0];
					const unsigned int inode1 = noes_co[1];
					double coord0[2], coord1[2];
					ns_co.GetValue(inode0,coord0);
					ns_co.GetValue(inode1,coord1);
					double len = sqrt( (coord1[0]-coord0[0])*(coord1[0]-coord0[0])+(coord1[1]-coord0[1])*(coord1[1]-coord0[1]) );
					es_va.GetNodes(ielem,noes_va);					
					res->AddValue(noes_va[0],1,load*len*0.5/tot_len);
					res->AddValue(noes_va[1],1,load*len*0.5/tot_len);
				}
			}
		}
	}*/
	const double norm_res = pLS->FinalizeMarge();
	this->m_is_cleared_value_ls = false;

	return norm_res;
}

bool CEqnSystem_Solid2D::MakePreconditioner(){
//	std::cout << "Value Preconditioner" << std::endl;
	if( pPrec==0 ){ this->InitializePreconditioner(); }
	// �O�����s������
	pPrec->SetValue( (*pLS).m_ls );
	this->m_is_cleared_value_prec = false;
	return true;
}

bool CEqnSystem_Solid2D::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{
//	std::cout << "Initialize LinearSystem" << std::endl;
	if( pLS  !=0 ) this->ClearLinearSystem();
	assert( pLS == 0 );
	// �A���ꎟ�������N���X�̍쐬
	if( this->m_IsSaveStiffMat ){
		bool is_nonlin;
		this->EqnationProperty(is_nonlin); 
		assert(!is_nonlin ); // ����`�������͍s��ۑ��ł��Ȃ�
		if( m_IsStationary ){ pLS = new CLinearSystem_Save; }
		else{ pLS = new CLinearSystem_SaveDiaM_NewmarkBeta; } // ���I�������̏ꍇ�̍s��ۑ��͎������Ă��Ȃ�
	}
	else{ pLS = new CLinearSystem_Field; }

	// �A���ꎟ�������N���X�̐ݒ�
	pLS->AddPattern_Field(m_IdFieldDisp,world);	// val_field����ł���S�̍����s���ǉ�����
	for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
		const unsigned int id_field = m_aIdFixField[idf].first; 
		const int idof = m_aIdFixField[idf].second; 
		if( idof == -1 ){
			pLS->SetFixedBoundaryCondition_Field( id_field, 0, world ); // bc0_field���Œ苫�E�����ɐݒ�	
			pLS->SetFixedBoundaryCondition_Field( id_field, 1, world ); // bc0_field���Œ苫�E�����ɐݒ�	
		}
		else{
			assert( idof < 2 );
			pLS->SetFixedBoundaryCondition_Field( id_field, idof, world ); // bc0_field���Œ苫�E�����ɐݒ�	
		}
	}
	return true;
}

bool CEqnSystem_Solid2D::InitializePreconditioner()
{		
//	std::cout << "Initialize Preconditioner" << std::endl;
	// �O�����N���X�̍쐬
	if( pPrec!=0 ) this->ClearPreconditioner();
	assert( pPrec == 0 );
    pPrec = new LsSol::CPreconditioner_ILU( (*pLS).m_ls, 1 );
/*	pPrec = new CPreconditioner_ILU;
	(*(CPreconditioner_ILU*)pPrec).SetOrdering(false);
	(*(CPreconditioner_ILU*)pPrec).SetFillInLevel(-1);
	pPrec->SetLinearSystem(*pLS);*/
	return true;
}

bool CEqnSystem_Solid2D::SetEquation( const CEqn_Solid2D& eqn )
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if( m_aEqn[ieqn].GetIdEA() == eqn.GetIdEA() ){
			m_aEqn[ieqn] = eqn;
			this->ClearLinearSystemPreconditioner();
			return true;
		}
	}
	return false;
}

CEqn_Solid2D CEqnSystem_Solid2D::GetEqnation(unsigned int id_ea) const
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if( m_aEqn[ieqn].GetIdEA() == id_ea ){
			return m_aEqn[ieqn];
		}
	}
	return CEqn_Solid2D(0,0);
}

bool CEqnSystem_Solid2D::Solve(Fem::Field::CFieldWorld& world)
{
//	std::cout << "CEqnSystem_Solid2D::Solve" << std::endl;
	bool is_nonlin;
	this->EqnationProperty(is_nonlin); 
	this->m_aItrNormRes.clear();
	if( is_nonlin ){	// ����`�̏ꍇ��Newton���������Ȃ���΂����Ȃ��̂ŕ�����
		assert( !this->m_IsSaveStiffMat );
		if( pLS   == 0 ){ this->InitializeLinearSystem(world); }
		if( pPrec == 0 ){ this->InitializePreconditioner();    }
        double ini_norm_res = 0;
		for(unsigned int iitr=0;iitr<10;iitr++){
			double norm_res = this->MakeLinearSystem(world,iitr==0);
			this->MakePreconditioner();
			if( iitr==0 ){
				if( norm_res < 1.0e-20 ) break;
				ini_norm_res = norm_res;
			}
			if( norm_res < ini_norm_res * 1.0e-6 ) break;
//			std::cout << iitr << " " << norm_res << " " << ini_norm_res << " " << norm_res / ini_norm_res << std::endl;
			{	// �s�������
				double conv_ratio = 1.0e-6;
				unsigned int max_iter = 1000;
                LsSol::CLinearSystemPreconditioner lsp( (*pLS).m_ls, *pPrec );
				Sol::Solve_PCG(conv_ratio,max_iter,lsp);
				this->m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
//				std::cout << max_iter << " " << conv_ratio << std::endl;
			}
			if( m_IsStationary ){
				pLS->UpdateValueOfField(m_IdFieldDisp,world,VALUE);
			}
			else{
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
						m_IdFieldDisp,world,iitr==0);
			}
		}
	}
	else{ 
//		std::cout << "MakeLinearSystem LinearSolid " << std::endl;
		if( pLS == 0 ){
			this->InitializeLinearSystem(world);
			this->MakeLinearSystem(world,true);
		}	
		else{
			if( !this->m_IsSaveStiffMat ){ this->MakeLinearSystem(world,true); }
		}
		if( this->m_IsSaveStiffMat ){ 
			if( this->m_is_cleared_value_ls ){
				this->MakeLinearSystem(world,true);
			}
			pLS->MakeResidual(world);
		}
		assert( this->m_is_cleared_value_ls   == false );
		////////////////
		if( pPrec == 0 ){
			this->InitializePreconditioner();
			this->MakePreconditioner();
		}
		else if( !this->m_IsSaveStiffMat || this->m_is_cleared_value_prec ){
			this->MakePreconditioner();
		}
		assert( this->m_is_cleared_value_prec == false );

		{	// �s�������
			double conv_ratio = 1.0e-6;
			unsigned int max_iter = 1000;
            LsSol::CLinearSystemPreconditioner lsp( (*pLS).m_ls, *pPrec );
			Sol::Solve_PCG(conv_ratio,max_iter,lsp);
//			Fem::Sol::Solve_CG(conv_ratio,max_iter,*pLS);
//			std::cout << max_iter << " " << conv_ratio << std::endl;
//			m_num_iter = max_iter;
//			m_conv_ratio = conv_ratio;
			this->m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
		}

		// �����X�V����
		if( m_IsStationary ){
			pLS->UpdateValueOfField(m_IdFieldDisp,world,VALUE); 
		}
		else{
			if( !this->m_IsSaveStiffMat ){
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
					m_IdFieldDisp,world,true);
			}
			else{
				pLS->UpdateValueOfField(m_IdFieldDisp,world,ACCELERATION);
			}
		}
	}
	return true;
}

bool CEqnSystem_Solid2D::UpdateDomain_Field(unsigned int id_base, Fem::Field::CFieldWorld& world)
{
	m_IdFieldDisp  = world.MakeField_FieldElemDim(id_base,2,
		Fem::Field::VECTOR2,VALUE|VELOCITY|ACCELERATION,CORNER);
    assert( world.IsIdField(m_IdFieldDisp) );

	{	// �����v�f�z��ID��������������������C������g���D�Ȃ���ΐV�K�ɒǉ�
		std::vector<CEqn_Solid2D> aEqn_old = m_aEqn;
		m_aEqn.clear();
		const CField& field = world.GetField(m_IdFieldDisp);
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			unsigned int ieqn0=0;
			for(;ieqn0<aEqn_old.size();ieqn0++){
				if( aEqn_old[ieqn0].GetIdEA() == id_ea ){ 
					const unsigned int ieqn1 = m_aEqn.size();
					m_aEqn.push_back( aEqn_old[ieqn0] );
					m_aEqn[ieqn1].SetIdFieldDisp(m_IdFieldDisp);
					break; 
				}
			}
			if( ieqn0 != aEqn_old.size() ){ continue; }
            CEqn_Solid2D eqn1(id_ea,m_IdFieldDisp);
			eqn1.SetYoungPoisson(m_young_back, m_poisson_back, m_is_plane_stress_back);
			eqn1.SetRho(m_rho_back);
            m_aEqn.push_back( eqn1 );
		}
	}

    std::cout << "Size Eqn : " << m_aEqn.size() << std::endl;
    for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
        const CEqn_Solid2D& eqn = m_aEqn[ieqn];
		double young, poisson;
		eqn.GetYoungPoisson(young,poisson);
//        std::cout << ieqn << " " << eqn.GetIdEA() << " " << young << " " << poisson << " " << std::endl;
		m_aEqn[ieqn].SetGravitation(0,-0.1);
    }
	

	if( !world.IsIdField(m_IdFieldDisp) ) return false;
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}
/*
bool CEqnSystem_Solid2D::SetDomain_Field(unsigned int id_field_base, Fem::Field::CFieldWorld& world)
{
	{	// ���̓t�B�[���h�̍��W�ߓ_�Z�O�����g��dof��2���ǂ����`�F�b�N����
//		unsigned int id_field_base = world.GetFieldBaseID();
		assert( world.IsIdField(id_field_base) );
		const CField& field_base = world.GetField(id_field_base);
		assert( field_base.GetNDimCoord() == 2 );
	}
	m_IdFieldDisp = world.MakeField_FieldElemDim(id_field_base, 2,
        VECTOR2,VALUE|VELOCITY|ACCELERATION,CORNER);
	const CIDConvEAMshCad conv = world.GetIDConverter(id_field_base);
	{
		m_aEqn.clear();
		const CField& field = world.GetField(m_IdFieldDisp);
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			const CElemAry& ea = world.GetEA(id_ea);
			if( ea.ElemType() != TRI && ea.ElemType() != QUAD ) continue;
			m_aEqn.push_back( CEqn_Solid2D(id_ea,m_IdFieldDisp) );
		}
	}
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}
*/
bool CEqnSystem_Solid2D::SetDomain_FieldEA(unsigned int id_field_base, unsigned int id_ea, 
                                           Fem::Field::CFieldWorld& world)
{
	{	// ���̓t�B�[���h�̍��W�ߓ_�Z�O�����g��dof��2���ǂ����`�F�b�N����
//		unsigned int id_field_base = world.GetFieldBaseID();
		assert( world.IsIdField(id_field_base) );
		const CField& field_base = world.GetField(id_field_base);
		assert( field_base.GetNDimCoord() == 2 );
	}
	m_IdFieldDisp = world.MakeField_FieldElemAry(id_field_base, id_ea,
        VECTOR2,VALUE|VELOCITY|ACCELERATION,CORNER);
	{
		m_aEqn.clear();
		const CField& field = world.GetField(m_IdFieldDisp);
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			const CElemAry& ea = world.GetEA(id_ea);
			if( ea.ElemType() != TRI && ea.ElemType() != QUAD ) continue;
			m_aEqn.push_back( CEqn_Solid2D(id_ea,m_IdFieldDisp) );
		}
	}
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}


/*
bool CEqnSystem_Solid2D::ToplogicalChangeCad_InsertLoop(Fem::Field::CFieldWorld& world, 
	unsigned int id_l_back, unsigned id_l_ins)
{
	{	// ���̓t�B�[���h�̍��W�ߓ_�Z�O�����g��dof��2���ǂ����`�F�b�N����
		unsigned int id_field_base = world.GetFieldBaseID();
		assert( world.IsIdField(id_field_base) );
		const CField& field_base = world.GetField(id_field_base);
		assert( field_base.GetNDimCoord() == 2 );
	}
	m_id_val = world.MakeField_AllRegion(VECTOR2,VALUE|VELOCITY|ACCELERATION,CORNER);
	{
		m_aEqn.clear();
		const CField& field = world.GetField(m_id_val);
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			m_aEqn.push_back( CEqn_Solid2D(id_ea,m_id_val) );
		}
	}
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}
*/


void CEqnSystem_Solid2D::SetLoad(double load, unsigned int iea, Fem::Field::CFieldWorld& world)
{
	unsigned int iload0;
	for(iload0=0;iload0<m_aLoad.size();iload0++){
		const unsigned int iea0 = m_aLoad[iload0].first;
		if( iea0 == iea ){
			m_aLoad[iload0].second = load;
			break;
		}
	}
	if( iload0 == m_aLoad.size() ){
		m_aLoad.push_back( std::make_pair(iea,load) );
	}
}


void CEqnSystem_Solid2D::ClearLoad(unsigned int id_ea)
{
	unsigned int iload0;
	for(iload0=0;iload0<m_aLoad.size();iload0++){
		const unsigned int iea0 = m_aLoad[iload0].first;
		if( iea0 == id_ea ){
			m_aLoad.erase( m_aLoad.begin()+iload0);
			break;
		}
	}
}




void CEqnSystem_Solid2D::SetYoungPoisson( double young, double poisson, bool is_plane_stress )
{
	this->m_young_back = young;
	this->m_poisson_back = poisson;
	this->m_is_plane_stress_back = is_plane_stress;
	////////////////
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetYoungPoisson(young,poisson,is_plane_stress);
	}
	////////////////
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Solid2D::SetRho( double rho )
{
	this->m_rho_back = rho;
	////////////////
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetRho(rho);
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Solid2D::SetGravitation( double g_x, double g_y )
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetGravitation(g_x,g_y);
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Solid2D::SetGeometricalNonlinear( bool is_nonlin )
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetGeometricalNonlinear(is_nonlin);
	}
	if( is_nonlin ){
		this->m_IsSaveStiffMat = false;
		this->ClearLinearSystem();
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Solid2D::SetThermalStress( unsigned int id_field_temperature )
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetThermalStress( id_field_temperature );
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Solid2D::SetStationary( bool is_stat )
{
	m_IsStationary = is_stat;
	this->ClearLinearSystemPreconditioner();
}

void CEqnSystem_Solid2D::SetSaveStiffMat( bool is_save )
{
	bool is_nonlin;
	this->EqnationProperty(is_nonlin);
	if( is_nonlin ){ this->m_IsSaveStiffMat = false; }
	else{ this->m_IsSaveStiffMat = is_save; }
	this->ClearLinearSystemPreconditioner();
}


// �X�J���[�̉��͑����l��ǉ�����D���̂������[�h������
// mode = 0 : �~�[�[�X
// mode = 1 : �ő剞��
// �􉽊w�I�Ȕ���`�����l�����邩�ǂ����́CEqnObj�����Č��߂�悤�ɂ���
bool CEqnSystem_Solid2D::SetScalarStressValue(unsigned int id_field_str, Fem::Field::CFieldWorld& world)
{	
	if( !world.IsIdField(id_field_str) ) return false;
	Fem::Field::CField& field_str = world.GetField(id_field_str);

	if( !world.IsIdField(m_IdFieldDisp) ) return false;
	Fem::Field::CField& field_dis = world.GetField(m_IdFieldDisp);

	const std::vector<unsigned int>& aIdEA_from = field_dis.GetAry_IdElemAry();
	const std::vector<unsigned int>& aIdEA_to   = field_str.GetAry_IdElemAry();
	if( aIdEA_from.size() != aIdEA_to.size() ) return false;

	const unsigned int niea = aIdEA_from.size();
	for(unsigned int iiea=0;iiea<niea;iiea++)
	{
		Fem::Field::INTERPOLATION_TYPE type_from, type_to;
		{
			if( aIdEA_from[iiea] != aIdEA_to[iiea] ){
				assert(0);
				return false;
			}
			const unsigned int id_ea = aIdEA_from[iiea];
			type_from = field_dis.GetInterpolationType(id_ea,world);
			type_to   = field_str.GetInterpolationType(id_ea,world);
		}

		unsigned int nnoes, ndim;
		if( type_from==TRI11 && type_to==TRI1001 ){
			nnoes = 3; ndim = 2;
		}
		else{
			std::cout << "Error!-->Not Implimented!" << std::endl;
			std::cout << type_from << " " << type_to << std::endl;
			assert(0);
			getchar();
		}

		unsigned int id_ea = aIdEA_to[iiea];
		const CElemAry& ea = world.GetEA(id_ea);
		const CElemAry::CElemSeg& es_c_co = field_dis.GetElemSeg(id_ea,CORNER,false,world);
		const CElemAry::CElemSeg& es_c_va = field_dis.GetElemSeg(id_ea,CORNER,true, world);
		const CElemAry::CElemSeg& es_b_va = field_str.GetElemSeg(id_ea,BUBBLE,true, world);

		Fem::Field::CField::CNodeSegInNodeAry nans_c = field_dis.GetNodeSegInNodeAry(CORNER);
		assert( world.IsIdNA(nans_c.id_na_co) );
		assert( world.IsIdNA(nans_c.id_na_va) );
		const CNodeAry& na_c_co = world.GetNA( nans_c.id_na_co);
		const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg( nans_c.id_ns_co);
		const CNodeAry& na_c_va = world.GetNA( nans_c.id_na_va);
		const CNodeAry::CNodeSeg& ns_c_va = na_c_va.GetSeg( nans_c.id_ns_va);

		Fem::Field::CField::CNodeSegInNodeAry nans_b = field_str.GetNodeSegInNodeAry(BUBBLE);
		unsigned int id_na_b_va = nans_b.id_na_va;
		unsigned int id_ns_b_va = nans_b.id_ns_va;

		assert( world.IsIdNA(id_na_b_va) );
		CNodeAry& na_b_va = world.GetNA(id_na_b_va);
		CNodeAry::CNodeSeg& ns_b_va = na_b_va.GetSeg(id_ns_b_va);

        const unsigned int nnoes_c = es_c_co.GetSizeNoes();
        assert( nnoes_c < 64 );
		unsigned int noes[64];

		for(unsigned int ielem=0;ielem<ea.Size();ielem++){
			////////////////
			double coord[16][3];
			// ���W(coord)�ƒl(value)�����
			es_c_co.GetNodes(ielem,noes);
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				unsigned int ipoi0 = noes[inoes];
				assert( ipoi0 < na_c_co.Size() );
				ns_c_co.GetValue(ipoi0,coord[inoes]);
			}
			////////////////
			double disp[16][3];	// �ߓ_�ψ�
			es_c_va.GetNodes(ielem,noes);
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				unsigned int ipoi0 = noes[inoes];
				assert( ipoi0 < na_c_va.Size() );
				for(unsigned int idim=0;idim<ndim;idim++){
					ns_c_va.GetValue(ipoi0,disp[inoes]);
				}
			}
			////////////////
			double dudx[2][2];	// �ό`���z
			if( type_from == TRI11 ){
				double dldx[3][2];
				double const_term[3];
				TriDlDx(dldx,const_term, coord[0],coord[1],coord[2]);
				for(unsigned int i=0;i<ndim*ndim;i++){ (&dudx[0][0])[i] = 0.0; }
				for(unsigned int knoes=0;knoes<nnoes;knoes++){
					dudx[0][0] += disp[knoes][0]*dldx[knoes][0];
					dudx[0][1] += disp[knoes][0]*dldx[knoes][1];
					dudx[1][0] += disp[knoes][1]*dldx[knoes][0];
					dudx[1][1] += disp[knoes][1]*dldx[knoes][1];
				}
			}
			double myu, lambda;
			{
				const CEqn_Solid2D& eqn = this->GetEqnation(id_ea);
				eqn.GetLambdaMyu(lambda,myu);
			}
			double strain[2][2];	// �c
			// �􉽊w�I����`����
//			strain[0][0] = 0.5*(dudx[0][0]+dudx[0][0]+dudx[0][0]*dudx[0][0]+dudx[1][0]*dudx[1][0]);
//			strain[0][1] = 0.5*(dudx[0][1]+dudx[1][0]+dudx[0][0]*dudx[0][1]+dudx[1][1]*dudx[1][0]);
//			strain[1][0] = 0.5*(dudx[1][0]+dudx[0][1]+dudx[0][1]*dudx[0][0]+dudx[1][0]*dudx[1][1]);
//			strain[1][1] = 0.5*(dudx[1][1]+dudx[1][1]+dudx[0][1]*dudx[0][1]+dudx[1][1]*dudx[1][1]);
			// �􉽊w�I����`�Ȃ�
			strain[0][0] = 0.5*(dudx[0][0]+dudx[0][0]);
			strain[0][1] = 0.5*(dudx[0][1]+dudx[1][0]);
			strain[1][0] = 0.5*(dudx[1][0]+dudx[0][1]);
			strain[1][1] = 0.5*(dudx[1][1]+dudx[1][1]);
			////////////////
			double stress[2][2];
			{
				stress[0][0] = myu*strain[0][0];
				stress[0][1] = myu*strain[0][1];
				stress[1][0] = myu*strain[1][0];
				stress[1][1] = myu*strain[1][1];
				const double dtmp1 = lambda*(strain[0][0]+strain[1][1]);
				stress[0][0] += dtmp1;
				stress[1][1] += dtmp1;
			}
			double mises;
			{
				const double dtmp1 = 0.5*stress[0][0]*stress[0][0]+0.5*stress[1][1]*stress[1][1]
				+0.5*(stress[1][1]-stress[0][0])*(stress[1][1]-stress[0][0])
				+3*stress[0][1]*stress[0][1];
				mises = sqrt(dtmp1);
			}
			////////////////
			double maxprinciple;
			{
				const double d1 = stress[0][0]+stress[1][1];
				const double d2 = stress[0][0]*stress[1][1]-stress[0][1]*stress[1][0];
				const double d3 = d1*d1-4*d2;
				assert( d3 >= 0 );
				maxprinciple = 0.5*(d1+sqrt(d3));
			}
			{
				unsigned int noes[16];
				es_b_va.GetNodes(ielem,noes);
				unsigned int ipoi0 = noes[0];
				assert( ipoi0 < na_b_va.Size() );
				ns_b_va.SetValue(ipoi0,0,mises);
			}
		}
	}
	return true;
}


bool CEqnSystem_Solid2D::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return true;
}

unsigned int CEqnSystem_Solid2D::AddFixElemAry( 
		unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world, idof );
}

bool CEqnSystem_Solid2D::ClearFixElemAry(
		unsigned int id_ea, Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA( id_ea ) ) return false;
	for(unsigned int ifix=0;ifix<m_aIdFixField.size();ifix++){
		const unsigned int id_field_fix = m_aIdFixField[ifix].first;
		const Fem::Field::CField& field = world.GetField(id_field_fix);
		const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
		if( aIdEA.size() != 1 ){
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
		}
		if( aIdEA[0] == id_ea ){
			m_aIdFixField.erase( m_aIdFixField.begin()+ifix );
			this->ClearLinearSystem();
		}
	}
	return true;
}


void CEqnSystem_Solid2D::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}

unsigned int CEqnSystem_Solid2D::AddFixElemAry( 
		const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
	const unsigned int id_field = world.GetPartialField(m_IdFieldDisp, aIdEA );
	if( id_field == 0 ) return 0;
	assert( world.IsIdField(id_field) );
	{
		CField& field = world.GetField(id_field);
		unsigned int nlen_val = field.GetNLenValue();
		for(unsigned int ilen=0;ilen<nlen_val;ilen++){
			field.SetValue(       0.0,ilen,world,false);
			field.SetVelocity(    0.0,ilen,world,false);
			field.SetAcceleration(0.0,ilen,world,false);
		}
	}
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return id_field;
}