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
@brief 場管理クラス(Fem::Field::CFieldWorld)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(FIELD_WORLD_H)
#define FIELD_WORLD_H

#include <map>

#include "delfem/elem_ary.h"	// CElemAryの参照返しをしているので必要
#include "delfem/node_ary.h"	// CNodeAryの参照返しをしているので必要
#include "delfem/objset.h"		// ID管理テンプレートクラス
#include "delfem/cad_com.h"		// ID管理テンプレートクラス

namespace Msh{
	class IMesh;	
}

namespace Fem{
namespace Field{

//! 場の種類
enum FIELD_TYPE{ 
	NO_VALUE,	//!< 設定されていない
	SCALAR,		//!< 実スカラー
	VECTOR2,	//!< ２次ベクトル
	VECTOR3,	//!< ３次ベクトル
	ZSCALAR		//!< 複素スカラー
};

class CField;	

/*! 
@brief 要素配列，メッシュ，CADのIDを解決するクラス
@ingroup Fem
*/
class CIDConvEAMshCad
{
	friend class CFieldWorld;
public:
	void Clear(){
		m_aIdAry.clear();
	}
    bool IsIdEA(unsigned int id_ea) const{
        for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
            if( m_aIdAry[iid].id_ea == id_ea ) return true;
        }
        return false;
    }
	unsigned int GetIdEA_fromMsh(unsigned int id_part_msh) const {
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
			if( m_aIdAry[iid].id_part_msh == id_part_msh ){
				return m_aIdAry[iid].id_ea;
			}
		}
		return 0;
	}
	unsigned int GetIdEA_fromMshExtrude(unsigned int id_part_msh, unsigned int inum_ext) const {
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
			if(    m_aIdAry[iid].id_part_msh_before_extrude == id_part_msh
                && m_aIdAry[iid].inum_extrude == inum_ext ){
				return m_aIdAry[iid].id_ea;
			}
		}
		return 0;
	}
	// itype_cad_part : Vertex(0), Edge(1), Loop(2)
    unsigned int GetIdEA_fromCad(unsigned int id_part_cad, Cad::CAD_ELEM_TYPE itype_cad_part, unsigned int inum_ext = 0) const {
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
			if(    m_aIdAry[iid].id_part_cad    == id_part_cad 
				&& m_aIdAry[iid].itype_part_cad == itype_cad_part 
                && m_aIdAry[iid].inum_extrude   == inum_ext ){
				return m_aIdAry[iid].id_ea;
			}
		}
		return 0;
	}
	void GetIdCad_fromIdEA(unsigned int id_ea, unsigned int& id_part_cad, Cad::CAD_ELEM_TYPE& itype_part_cad ) const{
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
			if( m_aIdAry[iid].id_ea == id_ea ){
				id_part_cad    = m_aIdAry[iid].id_part_cad;
				itype_part_cad = m_aIdAry[iid].itype_part_cad;
				return;
			}
		}
		id_part_cad    = 0;
		itype_part_cad = Cad::NOT_SET;
	}
private:
	class CInfoCadMshEA{
	public:
		unsigned int id_ea;
		unsigned int id_part_msh;
		unsigned int id_part_cad;
		Cad::CAD_ELEM_TYPE itype_part_cad;
		unsigned int id_part_msh_before_extrude;    // 突き出す前のメッシュID
		unsigned int inum_extrude;  // 突き出されてない(0), 底面(1), 側面(2), 上面(3)
	};
	std::vector<CInfoCadMshEA> m_aIdAry;
};

/*! 
@brief 場管理クラス
@ingroup Fem
*/
class CFieldWorld{
public:
	// デフォルトコンストラクタ
	CFieldWorld();
	// デストラクタ
	~CFieldWorld();
    /*!
	@brief ３次元メッシュから有限要素法補間場を構築
	@remarks 複数メッシュにも対応できるようにMeshのIDを返したい
	*/
	unsigned int AddMesh(const Msh::IMesh& mesh);

	unsigned int SetCustomBaseField(unsigned int id_base,
		std::vector<unsigned int> aIdEA_Inc,
		std::vector< std::vector<int> >& aLnods,
		std::vector<unsigned int>& mapVal2Co);

	CIDConvEAMshCad GetIDConverter(unsigned int id_field_base) const {
        std::map<unsigned int,CIDConvEAMshCad>::const_iterator itr = m_map_field_conv.find(id_field_base);
        if( itr == m_map_field_conv.end() ){
            CIDConvEAMshCad conv;
            return conv;
        }
        return itr->second;
    }

	//! 値を全て削除して初期化する
	void Clear();

	bool UpdateMeshCoord(    const unsigned int id_base, const Msh::IMesh& mesh);
	bool UpdateConnectivity( const unsigned int id_base, const Msh::IMesh& mesh );
	bool UpdateConnectivity_CustomBaseField(const unsigned int id_base,
		const std::vector<unsigned int>& aIdEA_Inc, 
		const std::vector< std::vector<int> >& aLnods,
		const std::vector<unsigned int>& mapVal2Co);

	////////////////////////////////////////////////////////////////
	// 要素配列に関係する関数群

	//! 要素配列のIDかどうか調べる
	bool IsIdEA( unsigned int id_ea ) const;
	//! 要素配列のIDを全て配列で得る
	const std::vector<unsigned int>& GetAry_IdEA() const;
	//! 要素配列の参照を得る関数(const)
	const CElemAry& GetEA(unsigned int id_ea) const;
	//! 要素配列の参照を得る関数(非const)
	CElemAry& GetEA(unsigned int id_ea);
	//! 要素配列を外部から追加する関数(.netからの追加に使う)
	unsigned int AddElemAry(unsigned int size, ELEM_TYPE elem_type);
	bool AddIncludeRelation(unsigned int id_ea, unsigned int id_ea_inc);	 // AddElemAryと統一したい．

	////////////////////////////////////////////////////////////////
	// 節点配列に関係する関数群

	//! 節点配列のIDかどうか調べる
	bool IsIdNA( unsigned int id_na ) const;
	//! 節点配列のIDを全て配列で得る
	const std::vector<unsigned int>& GetAry_IdNA() const;
	//! 節点配列の参照を得る関数(const)
	const CNodeAry& GetNA(unsigned int id_na) const;
	//! 節点配列の参照を得る関数(非const)
	CNodeAry& GetNA(unsigned int id_na);
	//! 節点配列を外部から追加する関数(.netからの追加に使う)
	unsigned int AddNodeAry(unsigned int size);

	////////////////////////////////////////////////////////////////
	// Fieldに関係する関数群

	//! FieldのIDかどうかを返す
	bool IsIdField( unsigned int id_field ) const;
	//! FieldのIDを全て配列で得る
	const std::vector<unsigned int>& GetAry_IdField() const;
	//! 場クラス(Fem::Field::CField)の参照を得る関数(const)
	const CField& GetField(unsigned int id_field) const;
	//! 場クラス(Fem::Field::CField)の参照を得る関数(非const)
	CField& GetField(unsigned int id_field);

	////////////////////////////////////////////////////////////////
	// Field追加関数

	/*! 
	@brief 入力された形状の座標をもつ場を全体に追加
	@param[in] field_type 値
	@param[in] derivative_type 速度、加速度を作るかどうか指定する。０の場合は何も作らない)
	@param[in] node_configuration_type 補間の種類
	*/
//	unsigned int MakeField_AllRegion(Field::FIELD_TYPE field_type = NO_VALUE, const int derivative_type = 1, const int node_configuration_type = 1 );
	unsigned int MakeField_FieldElemAry(unsigned int id_field, unsigned int id_ea, Field::FIELD_TYPE field_type = NO_VALUE, const int derivative_type = 1, const int node_configuration_type = 1 );
	unsigned int MakeField_FieldElemDim(unsigned int id_field, int idim_elem,      Field::FIELD_TYPE field_type = NO_VALUE, const int derivative_type = 1, const int node_configuration_type = 1 );
	//! 要素ID(id_ea)から成る部分場の取得
	unsigned int GetPartialField(unsigned int id_field, unsigned int IdEA );
	//! 要素ID配列(id_ea)から成る部分場の取得
	unsigned int GetPartialField(unsigned int id_field, std::vector<unsigned int> aIdEA);

	void FieldValueExec(double time);
	void FieldValueDependExec();

	////////////////////////////////////////////////////////////////
	// ファイルに関係する関数群

	int InitializeFromFile(const std::string& file_name, long& offset);
	int WriteToFile(const std::string& file_name, long& offset) const;
private:        
/*	unsigned int SetBaseField(
		unsigned int id_na, unsigned int id_ns_co,
		const std::vector< std::pair< unsigned int, unsigned int> >& pEaEs );*/
private:
	Com::CObjSet<CElemAry*> m_apEA;		//!< 要素配列集合
	Com::CObjSet<CNodeAry*> m_apNA;		//!< 節点配列集合
	Com::CObjSet<CField*> m_apField;	//!< 場集合

    std::map<unsigned int,CIDConvEAMshCad> m_map_field_conv;
};

}	// end namespace Field
}	// end namespace Fem

#endif
