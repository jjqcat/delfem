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
// FieldWorld.cpp：場管理クラス(CFieldWorld)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
    #pragma warning ( disable : 4996 )
#endif

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <assert.h>

#include "delfem/field_world.h"
#include "delfem/elem_ary.h"

#include "delfem/mesh_interface.h"

#include "delfem/field.h"

#include "delfem/matvec/vector_blk.h"

using namespace Fem::Field;

////////////////////////////////////////////////////////////////
// 生成/消滅
////////////////////////////////////////////////////////////////

CFieldWorld::CFieldWorld(){
	std::cout << "CFieldWorld::CFieldWorld" << std::endl;
}

CFieldWorld::~CFieldWorld(){
	std::cout << "CFieldWorld::~CFieldWorld" << std::endl;
	this->Clear();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void CFieldWorld::Clear()
{
	// 要素配列の削除
	std::vector<unsigned int> id_ary_pea = m_apEA.GetAry_ObjID();
	unsigned int iid_pea;
	for(iid_pea=0;iid_pea<id_ary_pea.size();iid_pea++){ 
		unsigned int id_pEA = id_ary_pea[iid_pea];
		assert( m_apEA.IsObjID(id_pEA) );
		CElemAry* pEA = m_apEA.GetObj(id_pEA);
		delete pEA;
	}
	m_apEA.Clear();

	// 節点配列の削除
	std::vector<unsigned int> id_ary_pna = m_apNA.GetAry_ObjID();
	unsigned int iid_pna;
	for(iid_pna=0;iid_pna<id_ary_pna.size();iid_pna++){ 
		unsigned int id_pNA = id_ary_pna[iid_pna];
		assert( m_apNA.IsObjID(id_pNA) );
		CNodeAry* pna = m_apNA.GetObj(id_pNA);
		delete pna;
	}
	m_apNA.Clear();

	// 場配列の削除
	std::vector<unsigned int> id_ary_pField = m_apField.GetAry_ObjID();
	unsigned int iid_pField;
	for(iid_pField=0;iid_pField<id_ary_pField.size();iid_pField++){ 
		unsigned int id_pField = id_ary_pField[iid_pField];
		assert( m_apField.IsObjID(id_pField) );
		CField* pField = m_apField.GetObj(id_pField);
		delete pField;
	}
	m_apField.Clear();

	m_map_field_conv.clear();
}

bool CFieldWorld::UpdateMeshCoord(const unsigned int id_base, const Msh::CMesh_Interface& mesh)
{
    assert( this->IsIdField(id_base) );
    if( !this->IsIdField(id_base) ) return false;
	CField& field_base = this->GetField(id_base);
	const unsigned int ndim = field_base.GetNDimCoord();
    assert( ndim == mesh.GetDimention() );
    if( ndim != mesh.GetDimention() ) return false;

	unsigned int id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	Field::CNodeAry& na = this->GetNA( id_na_co );
	Field::CNodeAry::CNodeSeg& ns_coord = field_base.GetNodeSeg(CORNER,false,*this,VALUE);

	std::vector<double> coord;
    mesh.GetCoord(coord);

    assert( coord.size() == ndim*na.Size() );
    if( coord.size() != ndim*na.Size() ) return false;
		
	if( ndim == 2 ){
		for(unsigned int inode=0;inode<na.Size();inode++){
			double x_value = coord[inode*2+0];
			ns_coord.SetValue(inode, 0, x_value );
			double y_value = coord[inode*2+1];
			ns_coord.SetValue(inode, 1, y_value );
		}
	}
	else if( ndim == 3 ){ 
		for(unsigned int inode=0;inode<na.Size();inode++){
			double x_value = coord[inode*3+0];
			ns_coord.SetValue(inode, 0, x_value );
			double y_value = coord[inode*3+1];
			ns_coord.SetValue(inode, 1, y_value );
			double z_value = coord[inode*3+2];
			ns_coord.SetValue(inode, 2, z_value );
		}
	}
	
	return true;
}

// 今の所三角形要素のリコネクティングにしか対応していない．
// 将来ConvをCADから切り離すためにid_cad系は使わずに，mshだけで対応する
bool CFieldWorld::UpdateConnectivity( const unsigned int id_base, const Msh::CMesh_Interface& mesh )
{
	assert( this->IsIdField(id_base) );
    if( !this->IsIdField(id_base) ) return false;
	CField& field_base = this->GetField(id_base);
	unsigned int id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	const CNodeAry& na = this->GetNA(id_na_co);
	const std::vector< std::pair<unsigned int, unsigned int> >& aEaEs = na.GetAryEaEs();
    const std::vector<unsigned int>& aIdMsh = mesh.GetAry_ID();
    CIDConvEAMshCad conv = this->GetIDConverter(id_base);
    std::vector<int> lnods;

	for(unsigned int iid_msh=0;iid_msh<aIdMsh.size();iid_msh++){
        const unsigned int id_msh = aIdMsh[iid_msh];
		const unsigned int id_ea  = conv.GetIdEA_fromMsh(id_msh);
		unsigned int id_es = 0;
		for(unsigned int iEaEs=0;iEaEs<aEaEs.size();iEaEs++){
			if( aEaEs[iEaEs].first == id_ea ){
				id_es = aEaEs[iEaEs].second;
				break;
			}
		}
        assert( this->IsIdEA(id_ea) );
		CElemAry& ea = this->GetEA(id_ea);
        assert( ea.IsSegID(id_es) );
		const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
        if( ea.ElemType() == TRI ){
            Msh::MSH_TYPE type = mesh.GetConnectivity(lnods,id_msh);
            assert( type == Msh::TRI );
		    for(unsigned int itri=0;itri<ea.Size();itri++){
			    for(unsigned int inotri=0;inotri<3;inotri++){
			        int ino0 = lnods[itri*3+inotri];
				    es.SetNodes(itri,inotri,ino0);
			    }
		    }
        }
	}
	return true;
}

unsigned int CFieldWorld::AddMesh(const Msh::CMesh_Interface& mesh)
{
	unsigned int id_na, id_ns_co;
	{
		std::vector<double> coord;
		mesh.GetCoord(coord);
		const unsigned int ndim = mesh.GetDimention();
		assert( coord.size() % ndim == 0 );
		const unsigned int nnode = coord.size() / ndim;
		CNodeAry* pna = new CNodeAry(nnode);
		id_ns_co = pna->GetFreeSegID();
		CNodeAry::CNodeSeg ns_co(ndim,"COORD");
		std::vector< std::pair<unsigned int,CNodeAry::CNodeSeg> > ns_input_ary;
		ns_input_ary.push_back( std::make_pair(id_ns_co,ns_co) );
		pna->AddSegment( ns_input_ary, coord );
		id_na = this->m_apNA.GetFreeObjID();
		const unsigned int tmp_id = this->m_apNA.AddObj( std::make_pair(id_na,pna) );
		assert( tmp_id == id_na );
	}

	unsigned int max_id_msh;
	{	// 要素IDの最大値を求める
		max_id_msh = 0;
		const std::vector<unsigned int>& aID = mesh.GetAry_ID();
		for(unsigned int iid=0;iid<aID.size();iid++){
			if( max_id_msh < aID[iid] ) max_id_msh = aID[iid];
		}
	}

    CIDConvEAMshCad conv;
	conv.m_aIdAry.reserve(max_id_msh+1);

	std::vector< std::pair<unsigned int, unsigned int> > aEaEs;

	{
		const std::vector<unsigned int>& aID = mesh.GetAry_ID();
		for(unsigned int iid=0;iid<aID.size();iid++){
			const unsigned int id_msh = aID[iid];
			std::vector<int> lnods;
			Msh::MSH_TYPE msh_type = mesh.GetConnectivity(lnods,id_msh);
			unsigned int id_cad_part, id_msh_before_ext, inum_ext;
            mesh.GetInfo(id_msh, id_cad_part, id_msh_before_ext, inum_ext);
			unsigned int itype_cad_part;
			unsigned int nnoel = 0;
			ELEM_TYPE elem_type;
			if(      msh_type == Msh::HEX     ){ 
                nnoel = 8; elem_type = HEX;   
                if( inum_ext == 0 ){ itype_cad_part = 3; }
                else{ 
                    assert( inum_ext % 2 == 0 );
                    itype_cad_part = 2; 
                }
            }
			else if( msh_type == Msh::TET     ){ 
                nnoel = 4; elem_type = TET;
                if( inum_ext == 0 ){ itype_cad_part = 3; }
                else{ 
                    assert( inum_ext%2 == 0 );
                    itype_cad_part = 2; 
                }   
            }
			else if( msh_type == Msh::QUAD    ){ 
                nnoel = 4; elem_type = QUAD;  
                if( inum_ext == 0 ){        itype_cad_part = 2; }
                else if( inum_ext%2 == 0 ){ itype_cad_part = 1; }   
                else{                       itype_cad_part = 2; }
            }
			else if( msh_type == Msh::TRI     ){ 
                nnoel = 3; elem_type = TRI;   
                if( inum_ext == 0 ){        itype_cad_part = 2; }
                else if( inum_ext%2 == 0 ){ itype_cad_part = 1; } 
                else{                       itype_cad_part = 2; }
            }
			else if( msh_type == Msh::BAR     ){ 
                nnoel = 2; elem_type = LINE;  
                if( inum_ext == 0 ){        itype_cad_part = 1; }
                else if( inum_ext%2 == 0 ){ itype_cad_part = 0; }   
                else{                       itype_cad_part = 1; }
            }
			else if( msh_type == Msh::VERTEX  ){ 
                nnoel = 1; elem_type = POINT; 
                if( inum_ext != 0 ){ assert( inum_ext % 2 == 1 ); }
                itype_cad_part = 0; 
            }
			else{ assert(0); }
			assert( lnods.size() % nnoel == 0 );
			const unsigned int nelem = lnods.size()/nnoel;
			CElemAry* pea = new CElemAry(nelem,elem_type);
			unsigned int id_es = pea->GetFreeSegID();
			std::vector<CElemAry::CElemSeg> es_ary;
			es_ary.push_back( CElemAry::CElemSeg(id_es,id_na,CORNER) );
			const std::vector<int>& res = pea->AddSegment(es_ary,lnods);
			assert( res.size() == 1 );
            assert( res[0] == (int)id_es );
			const unsigned int id_ea = this->m_apEA.GetFreeObjID();
			const unsigned int tmp_id = this->m_apEA.AddObj( std::make_pair(id_ea,pea)  );
			assert( tmp_id == id_ea );
			aEaEs.push_back( std::make_pair(id_ea,id_es) );
			CIDConvEAMshCad::CInfoCadMshEA info;
			{
				info.id_ea = id_ea;
				info.id_part_msh = id_msh;
				info.id_part_cad = id_cad_part;
				info.itype_part_cad = itype_cad_part;
                info.id_part_msh_before_extrude = id_msh_before_ext;
                info.inum_extrude = inum_ext;
			}
			conv.m_aIdAry.push_back( info );
		}
	}

	{	// 包含関係(include relation)を作る
		CNodeAry& na = this->GetNA(id_na);
		for(unsigned int ieaes=0;ieaes<aEaEs.size();ieaes++){
			na.AddEaEs( aEaEs[ieaes] );
		}
		const std::vector<unsigned int>& id_ary_msh = mesh.GetAry_ID();
		for(unsigned int iid_ary_msh=0;iid_ary_msh<id_ary_msh.size();iid_ary_msh++){
			const unsigned int id_msh = id_ary_msh[iid_ary_msh];
			const unsigned int id_ea = conv.GetIdEA_fromMsh(id_msh); assert( this->IsIdEA(id_ea) );
			unsigned int iEaEs = 0;
			for(;iEaEs<aEaEs.size();iEaEs++){
				if( id_ea == aEaEs[iEaEs].first ) break;
			}
			assert( iEaEs != aEaEs.size() );
			std::vector<unsigned int> id_ary_msh_inc = mesh.GetIncludeElemIDAry(id_msh);
			for(unsigned int iid_ary_msh_inc=0;iid_ary_msh_inc<id_ary_msh_inc.size();iid_ary_msh_inc++){
				unsigned int id_msh_inc = id_ary_msh_inc[iid_ary_msh_inc];
				unsigned int id_ea_inc = conv.GetIdEA_fromMsh(id_msh_inc);
				assert( this->IsIdEA(id_ea_inc) );
				unsigned int iEaEs_inc;
				for(iEaEs_inc=0;iEaEs_inc<aEaEs.size();iEaEs_inc++){
					if( id_ea_inc == aEaEs[iEaEs_inc].first ) break;
				}
				assert( iEaEs != aEaEs.size() );
				na.SetIncludeEaEs_InEaEs( aEaEs[iEaEs_inc], aEaEs[iEaEs] );
			}
		}
	}

	const unsigned int id_field_base = this->SetBaseField(id_na,id_ns_co,aEaEs);
	if( id_field_base == 0 ){ return 0; }
    m_map_field_conv.insert( std::make_pair(id_field_base, conv) );
	return id_field_base;
}


unsigned int CFieldWorld::SetCustomBaseField(unsigned int id_base,
	std::vector<unsigned int> aIdEA_Inc,
	std::vector< std::vector<int> >& aLnods,
	std::vector<unsigned int>& mapVal2Co)
{
	Fem::Field::CField& field_base = this->GetField(id_base);
	const unsigned int nno_v = mapVal2Co.size();
	unsigned int id_na_v = this->m_apNA.GetFreeObjID();
	id_na_v = m_apNA.AddObj( std::make_pair(id_na_v,new CNodeAry(nno_v)) );
	unsigned int id_na_c = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_ns_c = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;

	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	for(unsigned int iidea=0;iidea<aIdEA_Inc.size();iidea++){
		const unsigned int id_ea0 = aIdEA_Inc[iidea];
		unsigned int id_es_c = field_base.GetIdElemSeg(id_ea0,CORNER,false,*this);
		unsigned int id_es_v;
		{
			Fem::Field::CElemAry& ea = this->GetEA(id_ea0);
			id_es_v = ea.GetFreeSegID();
			std::vector<Fem::Field::CElemAry::CElemSeg> aEs;
			aEs.push_back( Fem::Field::CElemAry::CElemSeg(id_es_v,id_na_v,CORNER) );
			ea.AddSegment( aEs, aLnods[iidea] );
		}
		aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea0, id_es_v,id_es_c, 0,0, 0,0) );
	}

	CField* pField = new CField(
		0,
		aElemIntp,
		CField::CNodeSegInNodeAry(id_na_c,false,id_ns_c,  id_na_v,false,0,0,0), 
		CField::CNodeSegInNodeAry(),
		*this
		);
	const unsigned int id_base_new = this->m_apField.AddObj( std::make_pair(0,pField) );
	return id_base_new;
}


bool CFieldWorld::UpdateConnectivity_CustomBaseField(const unsigned int id_base,
	const std::vector<unsigned int>& aIdEA_Inc, 
	const std::vector< std::vector<int> >& aLnods,
	const std::vector<unsigned int>& mapVal2Co)
{
	Fem::Field::CField& fb = this->GetField(id_base);
	const unsigned int id_na = fb.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va;
	const std::vector<unsigned int>& aIdEA = fb.GetAry_IdElemAry();
	assert( aIdEA_Inc.size() == aIdEA.size() );
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		unsigned int id_ea = aIdEA[iiea];
		assert( aIdEA_Inc[iiea] == id_ea );
        assert( this->IsIdEA(id_ea) );
		CElemAry& ea = this->GetEA(id_ea);
		assert( ea.ElemType() == Fem::Field::TRI );
        unsigned int id_es = 0;
		{
			const std::vector<unsigned int>& aIdES = ea.GetAry_SegID();
			for(unsigned int iies=0;iies<aIdES.size();iies++){
				unsigned int id_es0 = aIdES[iies];
				const Fem::Field::CElemAry::CElemSeg& es = ea.GetSeg(id_es0);
				if( id_na == es.GetIdNA() ){ 
					id_es = id_es0; 
					break;
				}
			}
		}
        assert( ea.IsSegID(id_es) );
		const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
		const std::vector<int>& lnods = aLnods[iiea];
		assert( lnods.size() == ea.Size() * 3 );
		for(unsigned int itri=0;itri<ea.Size();itri++){
			for(unsigned int inotri=0;inotri<3;inotri++){
				int ino0 = lnods[itri*3+inotri];
				es.SetNodes(itri,inotri,ino0);
			}
        }
	}
	return true;
}

int CFieldWorld::InitializeFromFile(const std::string& file_name, long& offset)
{
	while(1){
		////////////////////////////////
		std::string strElemType;
		unsigned int ID;
		{
			FILE *fp;
			const unsigned int buff_size = 512;
			char stmp1[buff_size];

			if( (fp = fopen(file_name.c_str(),"r"))== NULL ){
				std::cout << "CFieldWorld::InitializeFromFile" << std::endl;
				std::cout << "Error!-->Cannot Open File : " << " " << file_name.c_str() << std::endl;
				return 1;
			}

			fseek(fp,offset,SEEK_SET);

			while( 1 ){
				if( fgets(stmp1,buff_size,fp) == NULL ){ return 1; }
				if( stmp1[0] == '#' ) continue;
//				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
				break;
			}
			assert( strncmp(stmp1,"$$$$$$$$",6)==0 );

			while( 1 ){
				if( fgets(stmp1,buff_size,fp) == NULL ){ return -1; }
				if( stmp1[0] == '#' ) continue;
//				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
				break;
			}
			strElemType = stmp1;
			{
				unsigned int ipos = strElemType.find("\n");
				strElemType.erase(ipos);
			}

			{
				while( fgets(stmp1,buff_size,fp) != NULL ){
					if( stmp1[0] == '#' ) continue;
//					if( strspn(stmp1," \n") != strlen(stmp1) ) break;
					break;
				}
				int tmp_id;
				sscanf(stmp1,"%d",&tmp_id);
				assert( tmp_id > 0 );
				ID = tmp_id;
			}
			fclose(fp);
		}
		////////////////////////////////
//		std::cout << strElemType << " " << ID << std::endl;

		if( strElemType == "ELEM_ARY" ){
			CElemAry* pea = new CElemAry;
			int ires =  pea->InitializeFromFile(file_name,offset);
			if( ires !=0 ){
				delete pea;
				assert( ires != 1 );	// Not EOF
				return ires;
			}
			this->m_apEA.AddObj( std::make_pair(ID,pea) );
		}
		else if( strElemType == "NODE" ){
			CNodeAry* pna = new CNodeAry;
			int ires =  pna->InitializeFromFile(file_name,offset);
			if( ires !=0 ){
				delete pna;
				assert( ires != 1 );	// Not EOF
				return ires;
			}
			this->m_apNA.AddObj( std::make_pair(ID,pna) );
		}
		else if( strElemType == "FIELD" ){
			CField* pfld = new CField;
			int ires = pfld->InitializeFromFile(file_name,offset);
			if( ires != 0 ){
				std::cout << "CFieldWorld fail field " << std::endl;
				delete pfld;
				assert( ires != 1 );
				return ires;
			}
			this->m_apField.AddObj( std::make_pair(ID,pfld) );
		}
		else{
			std::cout << "Skip Input : " << strElemType << std::endl;
			// 読み飛ばしルーティンをここに付け加える
			FILE *fp;
			const unsigned int buff_size = 512;
			char stmp1[buff_size];

			if( (fp = fopen(file_name.c_str(),"r"))== NULL ){
				std::cout << "Error!-->Cannot Open File" << std::endl;
				assert(0);
				return false;
			}

			fseek(fp,offset,SEEK_SET);
			fgets(stmp1,buff_size,fp);
			for(;;){
				offset = ftell(fp);
				if( fgets(stmp1,buff_size,fp)==NULL ) break;
				if( strncmp(stmp1,"$$$$$$",6) == 0 ){ break; }
			}
			fclose(fp);
			if( feof(fp) == 0 ) break;
		}
		{	// ファイルの終わりだったらBreakする部分
			FILE *fp;
			const unsigned int buff_size = 512;
			char stmp1[buff_size];
			if( (fp = fopen(file_name.c_str(),"r"))== NULL ){
				std::cout << "Error!-->Cannot Open File" << std::endl;
				assert(0);
				return false;
			}
			fseek(fp,offset,SEEK_SET);
			fgets(stmp1,buff_size,fp);
			int ires = feof(fp);
			fclose(fp);
			if( ires != 0 ){ break; }
		}
	}
/*	std::vector< std::pair<unsigned int, unsigned int> > aEaEs;
	this->SetBaseField(1,1,aEaEs);
	CNodeAry& na = this->GetNA(1);
	na.SetIncludeEaEs_InEaEs( std::make_pair(1,1), std::make_pair(7,1) );
	na.SetIncludeEaEs_InEaEs( std::make_pair(1,1), std::make_pair(8,1) );*/
	return 0;
}


int CFieldWorld::WriteToFile(const std::string& file_name, long& offset) const {
	{	// 節点配列の出力
		std::vector<unsigned int> id_ary_pna = m_apNA.GetAry_ObjID();
		unsigned int iid_pna;
		for(iid_pna=0;iid_pna<id_ary_pna.size();iid_pna++){ 
			unsigned int id_pna = id_ary_pna[iid_pna];
			std::cout << "id_pna " << id_pna << std::endl;
			assert( m_apNA.IsObjID(id_pna) );
			const CNodeAry* pna = m_apNA.GetObj(id_pna);
			assert( pna != 0 );
			int ires =  pna->WriteToFile(file_name,offset,id_pna);
			if( ires != 0 ){ return ires; }
		}
	}
	{	// 要素配列の出力
		std::vector<unsigned int> id_ary_pea = m_apEA.GetAry_ObjID();
		unsigned int iid_pea;
		for(iid_pea=0;iid_pea<id_ary_pea.size();iid_pea++){ 
			unsigned int id_pea = id_ary_pea[iid_pea];
			assert( m_apEA.IsObjID(id_pea) );
			const CElemAry* pea = m_apEA.GetObj(id_pea);
			assert( pea != 0 );
			int ires = pea->WriteToFile(file_name,offset,id_pea);
			if( ires != 0 ){ return ires; }
		}
	}
	{	// 場の出力
		std::vector<unsigned int> id_ary_pfld = this->m_apField.GetAry_ObjID();
		unsigned int iid_pfld;
		for(iid_pfld=0;iid_pfld<id_ary_pfld.size();iid_pfld++){ 
			unsigned int id_pfld = id_ary_pfld[iid_pfld];
			assert( m_apField.IsObjID(id_pfld) );
			const CField* pfld = m_apField.GetObj(id_pfld);
			assert( pfld != 0 );
			int ires = pfld->WriteToFile(file_name,offset,id_pfld);
			if( ires != 0 ){ return ires; }
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// 節点配列に関係する関数群

bool CFieldWorld::IsIdNA( unsigned int id_na ) const
{
	if( this->m_apNA.IsObjID(id_na) ){ return true; }
	return false;
}

const std::vector<unsigned int>& CFieldWorld::GetAry_IdNA() const 
{ 
	return this->m_apNA.GetAry_ObjID(); 
}

const CNodeAry& CFieldWorld::GetNA(unsigned int id_na) const
{
	assert( this->m_apNA.IsObjID(id_na) );
	if( !this->m_apNA.IsObjID(id_na) ) throw;
	return *m_apNA.GetObj(id_na);
}

CNodeAry& CFieldWorld::GetNA(unsigned int id_na)
{
	assert( this->m_apNA.IsObjID(id_na) );
	if( !this->m_apNA.IsObjID(id_na) ) throw;
	return *m_apNA.GetObj(id_na);
}

unsigned int CFieldWorld::AddNodeAry(unsigned int size)
{
	CNodeAry* pNA = new CNodeAry( size );
	return this->m_apNA.AddObj( std::make_pair(0,pNA) );
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// 要素配列に関係する関数群

bool CFieldWorld::IsIdEA( unsigned int id_ea ) const{
	if( this->m_apEA.IsObjID(id_ea) ){ return true; }
	return false;
}

const std::vector<unsigned int>& CFieldWorld::GetAry_IdEA() const { return this->m_apEA.GetAry_ObjID(); }

const CElemAry& CFieldWorld::GetEA(unsigned int id_ea) const{
	assert( this->m_apEA.IsObjID(id_ea) );
	if( !this->m_apEA.IsObjID(id_ea) ) throw;
	return *m_apEA.GetObj(id_ea);
}

CElemAry& CFieldWorld::GetEA(unsigned int id_ea){
	assert( this->m_apEA.IsObjID(id_ea) );
	if( !this->m_apEA.IsObjID(id_ea) ) throw;
	return *m_apEA.GetObj(id_ea);
}

unsigned int CFieldWorld::AddElemAry(unsigned int size, ELEM_TYPE elem_type){
	CElemAry* pEA = new CElemAry( size, elem_type );
	unsigned int id_ea = this->m_apEA.AddObj( std::make_pair(0,pEA) );
	return id_ea;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// Fieldに関係する関数群

bool CFieldWorld::IsIdField( unsigned int id_field ) const{
	if( this->m_apField.IsObjID(id_field) ){ return true; }
	return false;
}

const std::vector<unsigned int>& CFieldWorld::GetAry_IdField() const { return this->m_apField.GetAry_ObjID(); }

// Fieldの参照を得る関数（ポインタじゃなくて参照なのは、勝手にdeleteされないため）
const CField& CFieldWorld::GetField(unsigned int id_field) const{
	assert( this->m_apField.IsObjID(id_field) );
	if( !this->m_apField.IsObjID(id_field) ) throw "invalid field id";
	return *m_apField.GetObj(id_field);
}

CField& CFieldWorld::GetField(unsigned int id_field){
	assert( this->m_apField.IsObjID(id_field) );
	if( !this->m_apField.IsObjID(id_field) ) throw "invalid field id";
	return *m_apField.GetObj(id_field);
}

unsigned int CFieldWorld::SetBaseField(
	unsigned int id_na, unsigned int id_ns_co,
	const std::vector< std::pair< unsigned int, unsigned int> >& aEaEs )
{
	std::vector< Field::CField::CElemInterpolation > aElemIntp;
	for(unsigned int ieaes=0;ieaes<aEaEs.size();ieaes++){
		aElemIntp.push_back( 
			CField::CElemInterpolation( aEaEs[ieaes].first, 
				0,aEaEs[ieaes].second,
				0,0,
				0,0 ) 
		);
	}
	unsigned int ndim = 0;
	{
		assert( this->IsIdNA(id_na) );
		const CNodeAry& na = this->GetNA(id_na);
		assert( na.IsSegID(id_ns_co) );
		const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_co);
		ndim = ns.GetLength();
	}

	CField* pField = new CField(
		0,
		aElemIntp,
		CField::CNodeSegInNodeAry(id_na,false,id_ns_co,  0,false,0,0,0), 
		CField::CNodeSegInNodeAry(),
		*this
		);

	const unsigned int m_id_field_base = this->m_apField.AddObj( std::make_pair(0,pField) );
	return m_id_field_base;
}

unsigned int CFieldWorld::MakeField_FieldElemAry(
    unsigned int id_field_base,                                               
	unsigned int id_ea, 
	Fem::Field::FIELD_TYPE field_type, const int fdt, const int nct )
{
	if( !this->IsIdEA(id_ea) ) return 0;
	assert( nct == CORNER );

	unsigned int id_es_val = 0;
	unsigned int id_na_val = 0;
	{   // 既存のNAで対応できるか調べる
		CElemAry& ea = this->GetEA(id_ea);
		const std::vector<unsigned int> aIdEs = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<aIdEs.size();iid_es++){
			unsigned int id_es0 = aIdEs[iid_es];
			assert( ea.IsSegID(id_es0) );
			const CElemAry::CElemSeg es0 = ea.GetSeg(id_es0);
			if( es0.GetElSegType() != CORNER ){ continue; }
			const unsigned int id_na = es0.GetIdNA();
			assert( this->IsIdNA(id_na) );
			const CNodeAry& na = this->GetNA(id_na);
			assert( es0.GetMaxNoes() < na.Size() );
            if( es0.GetMaxNoes()+1 != na.Size() ){ continue; }
            ////////////////
			const unsigned int nnode = na.Size();
			std::vector<unsigned int> flg_vec;
			flg_vec.resize(nnode,0);
			unsigned int noes[10];
			unsigned int nnoes = es0.GetSizeNoes();
			for(unsigned int ielem=0;ielem<ea.Size();ielem++){
				es0.GetNodes(ielem,noes);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){ flg_vec[ noes[inoes] ] = 1; }
			}
			bool iflg = true;
			for(unsigned int inode=0;inode<nnode;inode++){
				if( flg_vec[inode] == 0 ){ iflg = false; break; }
			}
			if( iflg == true ){	// 既存のセグメントで対応できる場合
				std::cout << " Use Not New " << std::endl;
				id_es_val = id_es0;
				id_na_val = id_na;
				break;
			}
		}
		if( id_na_val == 0 && id_es_val == 0 ){	// 新しい値セグメントを追加
			unsigned int id_es1 = 0;
			for(unsigned int iid_es=0;iid_es<aIdEs.size();iid_es++){
				unsigned int id_es0 = aIdEs[iid_es];
				assert( ea.IsSegID(id_es0) );
				const CElemAry::CElemSeg es0 = ea.GetSeg(id_es0);
				if( es0.GetElSegType() == CORNER ){ 
					id_es1 = id_es0;
					break; 
				}
			}
			const CElemAry::CElemSeg es1 = ea.GetSeg(id_es1);
			unsigned int nnode;
			std::vector<int> flg_vec;
			{	// coからvaのテーブルであるflg_vecを作る
				const unsigned int mx_noes = es1.GetMaxNoes();
				flg_vec.resize(mx_noes+1,-2);
				unsigned int noes[10];
				unsigned int nnoes = es1.GetSizeNoes();
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es1.GetNodes(ielem,noes);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						flg_vec[ noes[inoes] ] = -1;
					}
				}
				nnode = 0;
				for(unsigned int i=0;i<mx_noes+1;i++){
					if( flg_vec[i] == -1 ){
						flg_vec[i] = nnode;
						nnode++;
					}
				}
			}
			std::cout << "New Nnode for Val : " << nnode << std::endl;
			id_na_val = this->AddNodeAry(nnode);
			std::vector<int> lnods;
			{	// lnodsを作る
				const unsigned int nnoes = es1.GetSizeNoes();
				const unsigned int nelem = ea.Size();
				lnods.resize( nnoes*nelem );;
				unsigned int noes[16];
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es1.GetNodes(ielem,noes);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						int inode0 = noes[inoes];
						int inode1 = flg_vec[inode0];
						assert( inode1 >= 0 );
						lnods[ielem*nnoes+inoes] = inode1;
					}
				}
			}
			id_es_val = ea.GetFreeSegID();
			std::vector<CElemAry::CElemSeg> es_ary;
			es_ary.push_back( CElemAry::CElemSeg(id_es_val,id_na_val,CORNER) );
			const std::vector<int>& res = ea.AddSegment(es_ary,lnods);
            assert( res.size() == 1 && res[0] == (int)id_es_val );
			{	// 包含関係を入れる
				CNodeAry& na_val = this->GetNA(id_na_val);
				na_val.AddEaEs( std::make_pair(id_ea,id_es_val) );
			}
		}	// 新しい要素/節点セグメントを追加
	}

	std::cout << "VAL   EA:" << id_ea << "  ES:" << id_es_val << "   NA:" << id_na_val << std::endl;

    assert( id_field_base != 0 );

	unsigned int ndim_coord;
	unsigned int id_na_co = 0;
	unsigned int id_ns_co = 0;
	bool is_part_co;
	unsigned int id_es_co = 0;
	{
		assert( this->IsIdField(id_field_base) );
		const CField& field_base = this->GetField(id_field_base);
		ndim_coord = field_base.GetNDimCoord();
		id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
		id_ns_co = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;
		CElemAry& ea = this->GetEA(id_ea);
		const std::vector<unsigned int> aIdEs = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<aIdEs.size();iid_es++){
			unsigned int id_es0 = aIdEs[iid_es];
			assert( ea.IsSegID(id_es0) );
			const CElemAry::CElemSeg es0 = ea.GetSeg(id_es0);
			if( es0.GetElSegType() != CORNER ){ continue; }
			if( id_na_co == es0.GetIdNA() ){
				id_es_co = id_es0;
			}
		}
		assert( id_es_co != 0 );
		////////////////
		assert( this->IsIdNA(id_na_co) );
		const CNodeAry& na = this->GetNA(id_na_co);
		const std::vector< std::pair<unsigned int,unsigned int> >& aEaEs = na.GetAry_EaEs_Min();
		if( aEaEs.size() == 1 && aEaEs[0].first == id_ea ){ is_part_co = false; }
		else{ is_part_co = true; }
	}

	Field::CField::CNodeSegInNodeAry na_c;
	{
		na_c.id_na_co = id_na_co;
		na_c.id_ns_co = id_ns_co;
		na_c.id_na_va = id_na_val;
		na_c.is_part_va = false;
		na_c.is_part_co = is_part_co;
	}

	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	{
		aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea, id_es_val,id_es_co, 0,0, 0,0) );
	}

	CField* pField = new CField(
		0,
		aElemIntp,
		na_c, Field::CField::CNodeSegInNodeAry(),
		*this );

	const unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	assert( id_field != 0 );

	pField->SetValueType(field_type,fdt,*this);
	{
		assert( this->IsIdField(id_field) );
		const CField& field = this->GetField(id_field);
		assert( field.AssertValid(*this) );
	}
	return id_field;
}


// 入力された形状の座標をもつ場を全体に追加(値、速度、加速度を作るかどうか指定する。全て偽の場合は何も作らない)
// 宣言の部分でデフォルト引数を使っていて、デフォルトで値だけをCORNER節点に持つ関数となっている．
// 節点全体を含む最少の要素配列IDで補間する
/*
unsigned int CFieldWorld::MakeField_AllRegion(unsigned int id_field_base,
                                              Fem::Field::FIELD_TYPE field_type, 
											  const int fdt, const int nct )
{
	std::cout << "CFieldWorld::MakeField_AllRegion" << std::endl;

	if( !(fdt&VALUE) && !(fdt&VELOCITY) && !(fdt&ACCELERATION) ) return 0;
	if( this->m_apNA.GetAry_ObjID().size() == 0 ) return 0;
	if( this->m_apEA.GetAry_ObjID().size() == 0 ) return 0;

	assert( this->IsIdField(m_id_field_base) );
	const CField& field_base = this->GetField(m_id_field_base);

	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	{
		const unsigned int id_na_c = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
		assert( this->IsIdNA(id_na_c) );
		const CNodeAry& na = this->GetNA(id_na_c);
		const std::vector< std::pair<unsigned int, unsigned int> >& aEaEsMin = na.GetAry_EaEs_Min();
		for(unsigned int iieaes=0;iieaes<aEaEsMin.size();iieaes++){
			const unsigned int id_ea = aEaEsMin[iieaes].first;
			const unsigned int id_es = aEaEsMin[iieaes].second;
			Fem::Field::CField::CElemInterpolation ei(id_ea,   0,id_es,   0,0,   0,0);
			aElemIntp.push_back(ei);
		}
	}

	Field::CField::CNodeSegInNodeAry na_c = field_base.GetNodeSegInNodeAry(CORNER);
	Field::CField::CNodeSegInNodeAry na_b;

	if( nct & CORNER ){
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			Fem::Field::CField::CElemInterpolation& ei = aElemIntp[iei];
			ei.id_es_c_va = ei.id_es_c_co;
		}
		na_c.id_na_va = na_c.id_na_co;
	}
	if( nct & BUBBLE ){
		unsigned int nbubble = 0;
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			unsigned int id_ea = aElemIntp[iei].id_ea;
			assert( this->IsIdEA(id_ea) );
			Fem::Field::CElemAry& ea = this->GetEA(id_ea);
			nbubble += ea.Size();
		}
		unsigned int id_na_b = this->AddNodeAry(nbubble);
		assert( this->IsIdNA(id_na_b) );
		na_b.id_na_va = id_na_b;

		if( aElemIntp.size() != 1 ){
			std::cout << "Error!-->未実装" << std::endl;
			assert(0);
		}
		{
			CField::CElemInterpolation& ei = aElemIntp[0];	// TODO:複数領域に対応していない
			unsigned int id_ea = ei.id_ea;
			CElemAry& ea = this->GetEA(id_ea);
			std::vector<int> lnods;
			{
				lnods.resize(ea.Size());
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					lnods[ielem] = ielem;
				}
			}
			std::vector<CElemAry::CElemSeg> es_ary;
			unsigned int id_es_b_va = ea.GetFreeSegID();
			es_ary.push_back( CElemAry::CElemSeg(id_es_b_va,id_na_b,BUBBLE) );
			const std::vector<int>& res = ea.AddSegment(es_ary,lnods);
			assert( res.size() == 1 && res[0] == id_es_b_va );
			ei.id_es_b_co = 0;
			ei.id_es_b_va = id_es_b_va;
		}
	}

	CField* pField = new CField(0,   aElemIntp,   na_c,na_b,   *this );

	const unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	assert( id_field != 0 );
	pField->SetValueType(field_type,fdt,*this);
	{
		assert( this->IsIdField(id_field) );
		const CField& field = this->GetField(id_field);
		assert( field.AssertValid(*this) );
	}
	return id_field;
}
*/

unsigned int CFieldWorld::MakeField_FieldElemDim(
		unsigned int id_field_base, 
        int idim_elem,  // 要素の次元（－１なら全ての要素）
		Fem::Field::FIELD_TYPE field_type, 
		const int fdt, const int nct )
{
	std::cout << "CFieldWorld::MakeField_Field" << std::endl;
	if( !(fdt&VALUE) && !(fdt&VELOCITY) && !(fdt&ACCELERATION) ) return 0;
	if( this->m_apNA.GetAry_ObjID().size() == 0 ) return 0;
	if( this->m_apEA.GetAry_ObjID().size() == 0 ) return 0;

	assert( this->IsIdField(id_field_base) );
	const CField& field_base = this->GetField(id_field_base);

	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	Fem::Field::CField::CNodeSegInNodeAry na_c, na_b;
    
	na_c = field_base.GetNodeSegInNodeAry(CORNER);
	na_c.id_na_va = 0;
    
    {
	    const std::vector<unsigned int>& aIdEA = field_base.GetAry_IdElemAry();
	    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		    const unsigned int id_ea = aIdEA[iiea];
            if( idim_elem != -1 ){
                const Fem::Field::CElemAry& ea = this->GetEA(id_ea);
                if( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) != idim_elem ) continue;
            }
		    const unsigned int id_es_c_co = field_base.GetIdElemSeg(id_ea,CORNER,false,*this);
            aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea, 0,id_es_c_co, 0,0, 0,0 ) );
        }
    }

	if( nct & CORNER ){
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			const unsigned int id_ea = aElemIntp[iei].id_ea;
            {
    			Fem::Field::CElemAry& ea = this->GetEA(id_ea);
                assert( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) == idim_elem );
            }
            if( field_base.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
				aElemIntp[iei].id_es_c_va = field_base.GetIdElemSeg(id_ea,CORNER,true, *this);
			}
            else{
				aElemIntp[iei].id_es_c_va = field_base.GetIdElemSeg(id_ea,CORNER,false,*this);
            }
		}
		na_c.id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
		na_c.id_ns_co = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;
        na_c.is_part_va = false;
        if( field_base.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
            na_c.id_na_va = field_base.GetNodeSegInNodeAry(CORNER).id_na_va;
        }
        else{
            na_c.id_na_va = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
        }
	}
	if( nct &  BUBBLE ){
        // 節点配列の追加
        unsigned int nbubble = 0;
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			const unsigned int id_ea = aElemIntp[iei].id_ea;
			Fem::Field::CElemAry& ea = this->GetEA(id_ea);
            assert( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) == idim_elem );
			nbubble += ea.Size();
		}
		unsigned int id_na = this->AddNodeAry(nbubble);
		assert( this->IsIdNA(id_na) );
		na_b.id_na_va = id_na;
        // 要素配列の追加
		unsigned int ibubble = 0;
		std::vector<int> lnods;
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			const unsigned int id_ea = aElemIntp[iei].id_ea;
			CElemAry& ea = this->GetEA(id_ea);
            assert( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) == idim_elem );
			{	// コネクティビティ(lnods)を作る
				lnods.resize(ea.Size());
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					lnods[ielem] = ibubble;
					ibubble++;
				}
			}
			std::vector<CElemAry::CElemSeg> es_ary;
			const unsigned int id_es_b_va = ea.GetFreeSegID();
			es_ary.push_back( CElemAry::CElemSeg(id_es_b_va,id_na,BUBBLE) );
            const std::vector<int>& ares = ea.AddSegment(es_ary,lnods);
            assert( ares.size() == 1 && ares[0] == (int)id_es_b_va );
            aElemIntp[iei].id_es_b_va = id_es_b_va;
		}
		assert( ibubble == nbubble );
	}
    
	CField* pField = new CField(
		0,
		aElemIntp,
		na_c, na_b,
		*this );

	const unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	assert( id_field != 0 );
	pField->SetValueType(field_type,fdt,*this);
	{
		assert( this->IsIdField(id_field) );
		const CField& field = this->GetField(id_field);
		assert( field.AssertValid(*this) );
	}
	return id_field;
}

unsigned int CFieldWorld::GetPartialField(unsigned int id_field_val, unsigned int id_ea)
{
    // いずれはid_eaが配列なバージョンを用いる．
	std::cout << "GetPartialField" << std::endl;

	if( !this->IsIdEA(id_ea) ){ return 0; }
	if( !this->IsIdField(id_field_val) ){ return 0; }

	const CField& field_val = this->GetField(id_field_val);
	if( !field_val.IsValid() ){ return 0; }
	if( !(field_val.GetIDFieldParent()==0) ){ return 0; }	// 親フィールドでなければならない．

	// 既に作られているPartialFieldで適合するものを探す
	{
		const std::vector<unsigned int> aIdField = this->m_apField.GetAry_ObjID();
		for(unsigned int iid_f=0;iid_f<aIdField.size();iid_f++){
			unsigned int id_f = aIdField[iid_f];
			const Fem::Field::CField& field = *m_apField.GetObj(id_f);
			if( field.GetIDFieldParent() != id_field_val ) continue;
			const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
			if( aIdEA.size() != 1 ) continue;
			if( aIdEA[0] == id_ea ){ return id_f; }
		}
	}

	{
		bool iflag = false;
		const unsigned int id_na_va_c = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		assert( this->IsIdEA(id_na_va_c) );
		const CNodeAry& na_va = this->GetNA(id_na_va_c);
		const unsigned int id_na_co_c = field_val.GetNodeSegInNodeAry(CORNER).id_na_co;
		assert( this->IsIdEA(id_na_co_c) );
		const CNodeAry& na_co = this->GetNA(id_na_co_c);
		const std::vector<unsigned int>& aIdEA = field_val.GetAry_IdElemAry();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea0 = aIdEA[iiea];
			unsigned int id_es0 = field_val.GetIdElemSeg(id_ea0,CORNER,true,*this);
			std::cout << "Try : " << id_ea0 << " " << id_es0 << " " << id_ea << std::endl;
			const unsigned int id_es1 = na_va.IsContainEa_InEaEs( std::make_pair(id_ea0,id_es0), id_ea );
			if( id_es1 != 0 ){
				iflag = true;
				break;
			}
			unsigned int id_es2 = field_val.GetIdElemSeg(id_ea0,CORNER,false,*this);
			std::cout << "Try : " << id_ea0 << " " << id_es2 << " " << id_ea << std::endl;
			const unsigned int id_es3 = na_co.IsContainEa_InEaEs( std::make_pair(id_ea0,id_es2), id_ea );
			if( id_es3 != 0 ){
				iflag = true;
				break;
			}
		}
		if( !iflag ){ 
            std::cout << " not partial field " << std::endl;
            return 0; 
        }
	}

	Field::CField::CElemInterpolation ei(id_ea, 0,0, 0,0, 0,0);
	{	// corner節点の追加
		unsigned int id_na_c_co  = field_val.GetNodeSegInNodeAry(CORNER).id_na_co;
		unsigned int id_na_c_val = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		assert( this->IsIdEA(id_ea) );
		const CElemAry& ea =this->GetEA(id_ea);
		const std::vector<unsigned int>& id_es_ary = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<id_es_ary.size();iid_es++){
			const unsigned int id_es = id_es_ary[iid_es];
			assert( ea.IsSegID(id_es) );
			const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
			const unsigned int id_na0 = es.GetIdNA();
			if( id_na0 == id_na_c_val ){ ei.id_es_c_va = id_es; }
			if( id_na0 == id_na_c_co  ){ ei.id_es_c_co = id_es; }
		}
	}
	assert( ei.id_es_c_co != 0 );
	if( ei.id_es_c_va == 0 ){	// この値場は新しく作られたものなので、部分場に新しい要素配列を追加する
		unsigned int id_es_co = ei.id_es_c_co;
		assert( this->IsIdEA(id_ea) );
		CElemAry& ea =this->GetEA(id_ea);
		assert( ea.IsSegID(id_es_co) );
		const CElemAry::CElemSeg& es_co = ea.GetSeg(id_es_co);
		std::vector<int> lnods;
		{
			std::vector<int> map_co2va;
			{
				const CNodeAry& na_co = this->GetNA( field_val.GetNodeSegInNodeAry(CORNER).id_na_co );
				unsigned int nnode_co = na_co.Size();
				map_co2va.resize( nnode_co, -1);
				const CNodeAry& na_va = this->GetNA( field_val.GetNodeSegInNodeAry(CORNER).id_na_va );
				unsigned int nnode_va = na_va.Size();
				for(unsigned int inode_va=0;inode_va<nnode_va;inode_va++){
					unsigned int inode_co = field_val.GetMapVal2Co(inode_va);
					map_co2va[inode_co] = inode_va;
				}
			}
			unsigned int nelem = ea.Size();
			unsigned int nnoes = es_co.GetSizeNoes();
			lnods.resize( nelem*nnoes );
			unsigned int noes[16];
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				es_co.GetNodes(ielem,noes);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){
					unsigned int inode_co = noes[inoes];
					int inode_va = map_co2va[inode_co];
					assert( inode_va != -1 );   // これが失敗した場合は，このEAは場の一部ではない
					lnods[ielem*nnoes+inoes] = inode_va;
				}
			}
		}
		unsigned int id_es_va = ea.GetFreeSegID();
		std::vector<CElemAry::CElemSeg> es_ary;
		es_ary.push_back( CElemAry::CElemSeg(id_es_va,
			field_val.GetNodeSegInNodeAry(CORNER).id_na_va,
			CORNER) );
		ea.AddSegment(es_ary,lnods);
		ei.id_es_c_va = id_es_va;
	}

	std::vector< Field::CField::CElemInterpolation > aElemIntp;
	{
		aElemIntp.push_back(ei);
	}

	
	Field::CField::CNodeSegInNodeAry na_c;
	{
		na_c = field_val.GetNodeSegInNodeAry(CORNER);
		na_c.is_part_co = true;
		na_c.is_part_va = true;
	}

	Field::CField::CNodeSegInNodeAry na_b;
	// TODO:部分場のBUBBLE節点が含まれるかどうかを判断しなければならない．
	// とりあえず今は境界条件のためのPartialFieldと思って保留
/*	{
		na_b = field_val.GetNodeSegInNodeAry(BUBBLE);
		na_b.is_part_co = true;
		na_b.is_part_va = true;
	}*/

	// ここに辺節点の追加に関する物を書く
	CField* pField = new CField(
		id_field_val,
		aElemIntp,
		na_c, na_b,
		*this );

	unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	return id_field;
}

unsigned int CFieldWorld::GetPartialField(unsigned int id_field_val, 
										  std::vector<unsigned int> id_ea_ary)
{
    if( id_ea_ary.size() == 0 ) return 0;
	std::cout << "GetPartialField" << std::endl;
	for(unsigned int iid_ea=0;iid_ea<id_ea_ary.size();iid_ea++){
		unsigned int id_ea = id_ea_ary[iid_ea];
		if( !this->IsIdEA(id_ea) ){	return 0; }
	}
	if( !this->IsIdField(id_field_val) ) return 0;
	const Fem::Field::CField& field_val = this->GetField(id_field_val);
	if( !field_val.IsValid() ){ return 0; }
	assert( field_val.AssertValid(*this) );

	std::vector< Field::CField::CElemInterpolation > aElemIntp;
	unsigned int id_na_c_co = field_val.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_na_c_val = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;

	std::vector<int> map_co2val;
	if( id_na_c_co != id_na_c_val ){
		const CNodeAry& na_co = this->GetNA(id_na_c_co);
		map_co2val.resize( na_co.Size(), -1 );
		const CNodeAry& na_va = this->GetNA(id_na_c_val);
		const unsigned int npoin_val = na_va.Size();
		for(unsigned int ipoin=0;ipoin<npoin_val;ipoin++){
			unsigned int ipoin_co0 = field_val.GetMapVal2Co(ipoin);
            assert( ipoin_co0 < map_co2val.size() );
			map_co2val[ipoin_co0] = ipoin;
		}
	}
	for(unsigned int iid_ea=0;iid_ea<id_ea_ary.size();iid_ea++){
		unsigned int id_ea = id_ea_ary[iid_ea];
		Field::CField::CElemInterpolation ei(id_ea, 0,0, 0,0, 0,0);
		assert( this->IsIdEA(id_ea) );
		CElemAry& ea =this->GetEA(id_ea);
		const std::vector<unsigned int>& id_es_ary = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<id_es_ary.size();iid_es++){
			const unsigned int id_es = id_es_ary[iid_es];
			assert( ea.IsSegID(id_es) );
			const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
			if( es.GetElSegType() != CORNER ) continue;
			const unsigned int id_na0 = es.GetIdNA();
			if( id_na0 == id_na_c_val ){ ei.id_es_c_va = id_es; }
			if( id_na0 == id_na_c_co  ){ ei.id_es_c_co = id_es; }
		}
		assert( ei.id_es_c_co != 0 );
		if( ei.id_es_c_va == 0 ){
			assert( map_co2val.size() != 0 );
			const CElemAry::CElemSeg& es_co = ea.GetSeg( ei.id_es_c_co );
			std::vector<int> lnods;
			const unsigned int nelem = ea.Size();
			const unsigned int nnoes = es_co.GetSizeNoes();
			lnods.resize(nelem*nnoes);
			unsigned int noes_co[32];
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				es_co.GetNodes(ielem,noes_co);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){
					unsigned int inode_co = noes_co[inoes];
                    assert( inode_co < map_co2val.size() );
					int inode_va = map_co2val[inode_co];
					assert( inode_va >= 0 );
					lnods[ielem*nnoes+inoes] = inode_va;
				}
			}
			std::vector< CElemAry::CElemSeg > es_ary;
			unsigned int id_es_val = ea.GetFreeSegID();
			es_ary.push_back( CElemAry::CElemSeg( id_es_val, id_na_c_val, CORNER) );
			ea.AddSegment( es_ary, lnods);
			ei.id_es_c_va = id_es_val;
		}
		aElemIntp.push_back( ei );
	}

	Field::CField::CNodeSegInNodeAry na_c;
	{
		na_c = field_val.GetNodeSegInNodeAry(CORNER);
		na_c.is_part_co = true;
		na_c.is_part_va = true;
	}

	Field::CField::CNodeSegInNodeAry na_b;
	// TODO:部分場のBUBBLE節点が含まれるかどうかを判断しなければならない．
	// とりあえず今は境界条件のためのPartialFieldと思って保留
/*	{
		na_b = field_val.GetNodeSegInNodeAry(BUBBLE);
		na_b.
		na_b.is_part_co = true;
		na_b.is_part_va = true;
	}*/

	// ここに辺節点の追加に関する物を書く
	CField* pField = new CField(
		id_field_val,
		aElemIntp,
		na_c,  na_b,
		*this
	);

	unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	return id_field;
}

void CFieldWorld::FieldValueExec(double time){
	const std::vector<unsigned int>& id_field_ary = this->m_apField.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_field_ary.size();iid++){
		unsigned int id_field = id_field_ary[iid];
		CField* pField = m_apField.GetObj(id_field);
		if( !pField->IsDepend() ){
			pField->ExecuteValue(time,*this);
		}
	}
}


void CFieldWorld::FieldValueDependExec(){
	const std::vector<unsigned int>& id_field_ary = this->m_apField.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_field_ary.size();iid++){
		unsigned int id_field = id_field_ary[iid];
		CField* pField = m_apField.GetObj(id_field);
		if( pField->IsDepend() ){
			pField->ExecuteValue(0.0,*this);
		}
	}
}

