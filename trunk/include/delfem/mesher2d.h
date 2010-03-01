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
@brief ２次元メッシュクラス(Msh::CMesher2D)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(MSHER_2D_H)
#define MSHER_2D_H

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <set>
#include <map>

#include "delfem/vector2d.h"
#include "delfem/serialize.h"
#include "delfem/mesh_interface.h"
#include "delfem/cad2d_interface.h"

#include "delfem/msh/meshkernel2d.h"

////////////////////////////////////////////////

namespace Msh{

/*! 
@addtogroup Msh2D
*/
// @{

//! 頂点構造体
struct SVertex{
public:
	SVertex() : id(0), id_v_cad(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_v_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int v;	//!< 点のID
};

//! 線要素配列
class CBarAry{
public:
	CBarAry() : id(0), id_e_cad(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_e_cad;	//!< CADの辺ID（CADに関連されてなければ０）
	unsigned int id_se[2];
	unsigned int id_lr[2];
	std::vector<SBar> m_aBar;	//!< 辺要素の配列
};

//! ２次元３角形要素配列
class CTriAry2D{
public:
	CTriAry2D() : id(0), id_l_cad(0), ilayer(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_l_cad;	//!< CADの面ID（CADに関連されてなければ０）
	int ilayer;
	std::vector<STri2D> m_aTri;	//!< ３角形要素の配列
};

//! ２次元４角形要素配列
class CQuadAry2D{
public:
	CQuadAry2D() : id(0), id_l_cad(0), ilayer(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_l_cad;	//!< CADの面ID(CADに関連されてなければ０)
	int ilayer;
	std::vector<SQuad2D> m_aQuad;	//!< ４角形要素の配列
};

////////////////////////////////////////////////

/*!
@brief ２次元メッシュクラス
@ingroup Msh2D

要素の種類に依存せずに通しでIDが振られている．
辺要素に関しては，CADの辺と同じ向きに要素番号順に並んでいる．（荷重境界条件からの制約）
*/
class CMesher2D : public IMesh
{
public:
	//! できるだけ少ない要素数でメッシュを切る
	CMesher2D(const Cad::ICad2D& cad_2d){ this->Tesselation(cad_2d); }
	//! 要素の長さがelenとなるようにメッシュを切る
	CMesher2D(const Cad::ICad2D& cad_2d, double elen){ this->Meshing_ElemLength(cad_2d,elen); }
	//! デフォルトコンストラクタ
	CMesher2D(){}
    virtual ~CMesher2D(){}

	////////////////////////////////
	
	//! デフォルトコンストラクタで初期化された後や，Clearされた後でメッシュを切る(elen : メッシュ幅)
	bool Meshing_ElemLength(const Cad::ICad2D& cad_2d,double elen, unsigned int id_l=0);
	bool Meshing_ElemLength(const Cad::ICad2D& cad_2d,double elen,  const std::vector<unsigned int>& aIdLoop);
	//! デフォルトコンストラクタで初期化された後や，Clearされた後でメッシュを切る(esize : 要素数)
	bool Meshing_ElemSize( const Cad::ICad2D& cad_2d, unsigned int esize, unsigned int id_l = 0);
	bool Meshing_ElemSize( const Cad::ICad2D& cad_2d, unsigned int esize, const std::vector<unsigned int>& aIdLoop);
	//! 出来るだけ点を追加しないように、メッシュを切る
	bool Tesselation(const Cad::ICad2D& cad_2d, unsigned int id_l = 0);
	bool Tesselation(const Cad::ICad2D& cad_2d, const std::vector<unsigned int>& aIdLoop);

	////////////////////////////////

	virtual unsigned int GetDimention() const{ return 2; }	//!< 座標の次元（２）を返す
	virtual void GetInfo(unsigned int id_msh,
        unsigned int& id_cad, unsigned int& id_msh_before_ext, unsigned int& inum_ext,
		int& ilayer) const
    { 
		const int itype = m_ElemType[id_msh];
		const int iloc = m_ElemLoc[id_msh];
		if(      itype == 0 ){ id_cad = m_aVertex[ iloc].id_v_cad; }
		else if( itype == 1 ){ id_cad = m_aBarAry[ iloc].id_e_cad; }
		else if( itype == 2 ){ id_cad = m_aTriAry[ iloc].id_l_cad; }
        else if( itype == 3 ){ id_cad = m_aQuadAry[iloc].id_l_cad; }
		else{ assert(0); }
        id_msh_before_ext = 0;
        inum_ext = 0;
	}
	virtual void GetCoord(std::vector<double>& coord) const{
		unsigned int nnode = aVec2D.size();
		coord.resize( nnode*2 );
		for(unsigned int inode=0;inode<nnode;inode++){
			coord[inode*2  ] = aVec2D[inode].x;
			coord[inode*2+1] = aVec2D[inode].y;
		}
	}
	virtual std::vector<unsigned int> GetAry_ID() const{
		std::vector<unsigned int> id_ary;
		for(unsigned int id=1;id<m_ElemLoc.size();id++){
			if( m_ElemLoc[id] == -1 ) continue;
			id_ary.push_back(id);
		}
		return id_ary;
	}
	virtual std::vector<unsigned int> GetIncludeElemIDAry(unsigned int id_msh) const{
		{	// エラーの場合は空の配列を返す
			std::vector<unsigned int> id_ary;
			if( id_msh >= m_ElemLoc.size() ) return id_ary;
			if( m_ElemLoc[id_msh] == -1 ) return id_ary;
			if( id_msh >= this->m_include_relation.size() ) return id_ary;
		}
		return m_include_relation[id_msh];
	}
	/*! メッシュの接続関係を取得する
	コピーを無くしたデータの受け渡しを実装する予定
	*/
	MSH_TYPE GetConnectivity(unsigned int id_msh,std::vector<int>& lnods) const;

	bool GetClipedMesh(
		std::vector< std::vector<int> >& lnods_tri, 
		std::vector<unsigned int>& mapVal2Co, 
		const std::vector<unsigned int>& aIdMsh_Ind,
		const std::vector<unsigned int>& aIdMshBar_Cut );

	//! データをすべてクリアする
	virtual void Clear();

	void SmoothingMesh_Laplace(unsigned int num_iter);
	void SmoothingMesh_Delaunay(unsigned int& num_reconnect);

	////////////////////////////////
	// const関数
	
	unsigned int GetElemID_FromCadID(unsigned int id_cad,  Cad::CAD_ELEM_TYPE type_cad) const;
	bool IsID(unsigned int id) const;
	bool GetMshInfo(unsigned int id, unsigned int& nelem, MSH_TYPE& msh_type, unsigned int& iloc,  unsigned int& id_cad ) const;

	////////////////////////////////
	// 要素配列、節点配列に関するGetメソッド
    
	const std::vector<CTriAry2D>& GetTriArySet() const { return m_aTriAry; }
	const std::vector<CQuadAry2D>& GetQuadArySet() const { return m_aQuadAry; }
	const std::vector<CBarAry>& GetBarArySet() const { return m_aBarAry; }
	const std::vector<SVertex>& GetVertexAry() const { return m_aVertex; }
	const std::vector<Com::CVector2D>& GetVectorAry() const { return aVec2D; }

	////////////////////////////////////////////////////////////////
	// IOメソッド
	
	//! 読み込み書き出し
	bool Serialize( Com::CSerializer& serialize );	
	//! GiDメッシュの読み込み
	bool ReadFromFile_GiDMsh(const std::string& file_name);

protected:
	unsigned int FindMaxID() const;
	int CheckMesh();	// 異常がなければ０を返す
	unsigned int GetFreeObjID();
	void MakeElemLocationType(); // 要素の場所と種類をハッシュ（IDが引数）している配列を初期化
	void MakeIncludeRelation(const Cad::ICad2D& cad);	// include_relationを作る

	////////////////
	bool MakeBoundary_SplitBarAry(std::vector<SBar>& aBar, unsigned int id_inner);
	// もしもis_inverted == trueならmax_aspectは無効な値が入る
	void CheckMeshQuality(bool& is_inverted, double& max_aspect, const double ave_edge_len);

	////////////////
	// メッシュ切り関係のルーティン

	bool MakeMesh_Edge(const Cad::ICad2D& cad_2d, unsigned int id_e, const double len);
	bool MakeMesh_Loop(const Cad::ICad2D& cad_2d, unsigned int id_l, const double len);

	bool Tesselate_LoopAround(const Cad::ICad2D& cad_2d, const unsigned int id_l);
	bool Tessalate_Edge(const Cad::ICad2D& cad_2d, const unsigned int id_e );
	bool Tesselate_Loop(const Cad::ICad2D& cad_2d, const unsigned int id_l);

	////////////////////////////////
	// トポロジー取得

	bool GetMesh_FrontBackStartEnd(unsigned int id_msh_bar, 
		unsigned int& id_msh_f, unsigned int& id_msh_b, 
		unsigned int& id_msh_s, unsigned int& id_msh_e) const;	
	
	class CTriAround
	{
	public:
		CTriAround(unsigned int id_msh, unsigned int ielem, unsigned int inoel )
			: id_msh(id_msh), ielem(ielem), inoel(inoel){}
		unsigned int id_msh;
		unsigned int ielem;
		unsigned int inoel;
	};

	bool GetTriMesh_Around(
		std::vector< CTriAround >& aTriAround,
		unsigned int ibarary, 
		unsigned int ibar, unsigned int inobar, bool is_left,
		const std::vector<unsigned int>& aFlgMshCut );

	////////////////////////////////
	// CADとの接続関係のルーティン

	// もしここにメッシュが切られていないとfalseを返す
	bool FindElemLocType_CadIDType(
		unsigned int& iloc, unsigned int& itype, 
		unsigned int id_cad_part, Cad::CAD_ELEM_TYPE itype_cad );

	// このループの面積と，現在切られているメッシュから最適な辺の長さを決定する
	double GetAverageEdgeLength(const Cad::ICad2D& cad_2d, 
        const std::set<unsigned int>& aIdL);

protected:
	std::vector<int> m_ElemType;	// vertex(0) bar(1) tri(2) quad(3)	いつでも有効(なければ-1が返る)
	std::vector<int> m_ElemLoc;		// index of elem_ary				いつでも有効(なければ-1が返る)　
	std::vector< std::vector<unsigned int> > m_include_relation;	// どの要素配列がどの要素配列を含んでるか

	std::vector<SVertex> m_aVertex;		// type(0)
	std::vector<CBarAry> m_aBarAry;		// type(1)
	std::vector<CTriAry2D> m_aTriAry;	// type(2)
	std::vector<CQuadAry2D> m_aQuadAry;	// type(3)

	std::vector<Com::CVector2D> aVec2D;
};

// @}
}

#endif
