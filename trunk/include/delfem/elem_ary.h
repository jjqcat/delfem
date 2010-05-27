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
@brief 要素配列クラス(Fem::Field::CElemAry)のインターフェース
@author Nobuyuki Umetani
*/

#if defined(__VISUALC__)
    #pragma warning( disable : 4786 )
#endif

#if !defined(ELEM_ARY_H)
#define ELEM_ARY_H

#include <vector>
#include <assert.h>
#include <cstdlib> //(abs)
#include <cstring> //(strspn, strlen, strncmp, strtok)

#include "delfem/objset.h"
#include "delfem/indexed_array.h"

namespace Fem{
namespace Field{

//! 要素セグメントの種類
enum ELSEG_TYPE{ 
	CORNER=1,	//!< コーナー節点
	EDGE=2,		//!< 辺上の節点
	BUBBLE=4	//!< 要素内の節点
};

// この整数のIndexはelem_ary.cppで使われるので，むやみに変更しないこと
// ここだけは別のヘッダファイルに移したほうがいいかも．このIndexだけ欲しいクラス(drawer_field.hとか)のために
//! 要素の種類
enum ELEM_TYPE{	
	ELEM_TYPE_NOT_SET=0,
	POINT=1,	//!< 点要素
	LINE=2,		//!< 線要素
	TRI=3, 		//!< 三角形要素
	QUAD=4, 	//!< 四角形要素
	TET=5, 		//!< 四面体要素
	HEX=6 		//!< 六面体要素
};

/*! 
@brief 要素配列クラス
@ingroup Fem
*/
class CElemAry
{
	
public:
	//! 要素セグメント．要素の場所(Corner,Bubble..etc)ごとの節点番号を取得するインターフェース関数
	class CElemSeg{
		friend class CElemAry;
	public:
//		// そのうちコンストラクタを作る予定
		CElemSeg(unsigned int id, unsigned int id_na, ELSEG_TYPE elseg_type)
			: m_id(id), m_id_na(id_na), m_elseg_type(elseg_type){}

		// Getメソッド
		unsigned int GetMaxNoes() const { return max_noes; }	//!< ノード番号の一番大きなものを得る（このnoesを格納するためには一つ大きな配列が必要なので注意）
		unsigned int GetSizeNoes() const { return m_nnoes; }	//!< 要素辺りの節点数を返す
		unsigned int GetSizeElem() const { return nelem; }	//!< 要素数を得る
		unsigned int GetIdNA() const { return m_id_na; }	//!< この要素セグメントが参照する節点配列IDを得る
		unsigned int GetID() const { return m_id; }	//!< 要素セグメントIDを得る
		ELSEG_TYPE GetElSegType() const { return m_elseg_type; }	//!< 要素セグメントタイプ(Fem::Field::CORNER,Fem::Field::BUBBLE,Fem::Field::EDGE)を得る

		/*
		// 辺は向きを持っているのでコネクティビティが符号を持つようにしたい
		void GetNodes(const unsigned int& ielem, int* noes ) const {
			for(unsigned int inoes=0;inoes<m_nnoes;inoes++){ 
				noes[inoes] = pLnods[ielem*npoel+begin+inoes]; 
			}
		}*/
		//! 節点番号を取得
		void GetNodes(const unsigned int& ielem, unsigned int* noes ) const {
			for(unsigned int inoes=0;inoes<m_nnoes;inoes++){ 
				assert( ielem*npoel+begin+inoes < nelem*npoel );
				noes[inoes] = abs(pLnods[ielem*npoel+begin+inoes]); 
			}
		}
		//! 節点番号を設定
		void SetNodes(unsigned int ielem, unsigned int idofes, int ino ) const {
			pLnods[ielem*npoel+begin+idofes] = ino;
		}
		//! 各場所(Corner,Bubble)に定義されている要素節点の数を出す
		static unsigned GetLength(ELSEG_TYPE elseg_type, ELEM_TYPE elem_type)
		{	// スタティックな関数なので実体は持たなくて良い
			if( elem_type == POINT ){
				if( elseg_type == CORNER ){ return 1; }
				if( elseg_type == EDGE   ){ return 2; }
				if( elseg_type == BUBBLE ){ return 2; }	// 点と点を結ぶLagrange未定乗数を定義するときに使うから０じゃダメ
				assert(0);
			}
			else if( elem_type == LINE ){
				if( elseg_type == CORNER ){ return 2; }
				if( elseg_type == EDGE   ){ return 1; }
				if( elseg_type == BUBBLE ){ return 0; }
				assert(0);
			}
			else if( elem_type == TRI ){
				if( elseg_type == CORNER ){ return 3; }
				if( elseg_type == EDGE   ){ return 3; }
				if( elseg_type == BUBBLE ){ return 1; }
				assert(0);
			}
			else if( elem_type == QUAD ){
				if( elseg_type == CORNER ){ return 4; }
				if( elseg_type == BUBBLE ){ return 1; }
				if( elseg_type == EDGE   ){ return 4; }
				assert(0);
			}
			else if( elem_type == TET ){
				if( elseg_type == CORNER ){ return 4; }
				if( elseg_type == BUBBLE ){ return 1; }
				if( elseg_type == EDGE   ){ return 6; }
				assert(0);
			}
			else if( elem_type == HEX ){
				if( elseg_type == CORNER ){ return 8; }
				if( elseg_type == BUBBLE ){ return 1; }
				if( elseg_type == EDGE   ){ return 12; }
				assert(0);
			}
			else{ assert(0); }
			return 0;
		}
		//! 要素の次元を取得する
        static unsigned GetElemDim(ELEM_TYPE elem_type){
            if(      elem_type == POINT ){ return 0; }
            else if( elem_type == LINE  ){ return 1; }
            else if( elem_type == TRI   ){ return 2; }
            else if( elem_type == QUAD  ){ return 2; }
            else if( elem_type == TET   ){ return 3; }
            else if( elem_type == HEX   ){ return 3; }
			else{ assert(0); }
			return 0;
        }
	private: // 初期化で必要な関数
		unsigned int m_id;	// <- そのうち廃止
		unsigned int m_id_na; // <- 復活、ないと大変なことになった。
		enum ELSEG_TYPE m_elseg_type;
		// int id_es_corner <- そのうち追加予定
	private: // 初期化にいらない変数
		unsigned int max_noes;	// このElementSegmentに属する節点番号の最大のもの
		unsigned int begin;		// m_aLnodsのどこから始まるか (<npoel)
		unsigned int m_nnoes;		// ElementSegmentの長さ
	private: // CElemAryによって随時与えられる変数
		mutable int* pLnods;
		mutable unsigned int npoel;
		mutable unsigned int nelem;
	};
public:
	//! デフォルト・コンストラクタ
	CElemAry(){
		m_nElem = 0; npoel = 0; m_pLnods = 0;
	}
	/*! 
	@brief コンストラクタ
	@param [in] nelem 要素数
	@param [in] elem_type 要素の種類(TRIとかHEXとか)
	*/
	CElemAry(unsigned int nelem, ELEM_TYPE elem_type) : m_nElem(nelem), m_ElemType(elem_type){
		npoel = 0; m_pLnods = 0;
	}
	//! デストラクタ
	virtual ~CElemAry(){
		if( this->m_pLnods != 0 ) delete[] m_pLnods;
	}

	bool IsSegID( unsigned int id_es ) const { return m_aSeg.IsObjID(id_es); }	//!< 引数が使用されている要素セグメントのIDか調べる関数
	const std::vector<unsigned int>& GetAry_SegID() const { return this->m_aSeg.GetAry_ObjID(); } //!< 要素セグメントIDの配列を得る関数
	unsigned int GetFreeSegID() const{ return m_aSeg.GetFreeObjID(); }	//!< 使われていない要素セグメントIDを取得する関数
	//! 要素セグメントを取得する関数
	const CElemSeg& GetSeg(unsigned int id_es) const{			
		assert( this->m_aSeg.IsObjID(id_es) );
		if( !m_aSeg.IsObjID(id_es) ) throw;
		const CElemSeg& es = m_aSeg.GetObj(id_es);
		es.pLnods = this->m_pLnods;
		es.npoel = this->npoel;
		es.nelem = this->m_nElem;
		return es;
	}

	virtual ELEM_TYPE ElemType() const { return m_ElemType; }	//!< 要素の種類を取得
	virtual unsigned int Size() const { return m_nElem; }	//!< 要素数を取得

	//! CRSデータを作る(非対角用)
	virtual bool MakePattern_FEM(	
		const unsigned int& id_es0, const unsigned int& id_es1, 
        Com::CIndexedArray& crs ) const;

	//! CRSデータを作る(対角用)
	virtual bool MakePattern_FEM(	
		const unsigned int& id_es0, 
        Com::CIndexedArray& crs ) const;

	/*!
	２次節点のための関数
	*/
	virtual bool MakeEdge(const unsigned int& id_es_co, unsigned int& nedge, std::vector<unsigned int>& edge_ary) const;
	
	/*!
	２次節点のための関数
	*/
	virtual bool MakeElemToEdge(const unsigned int& id_es_corner, 
		const unsigned int& nedge, const std::vector<unsigned int>& edge_ary,
		std::vector<int>& el2ed ) const;
	/*!
	@brief 境界要素を作る(可視化のための関数)
	*/
	virtual CElemAry* MakeBoundElemAry(unsigned int id_es_corner, unsigned int& id_es_add, std::vector<unsigned int>& aIndElemFace) const;
	
	//! ファイルIO関数
	int InitializeFromFile(const std::string& file_name, long& offset);
	int WriteToFile(       const std::string& file_name, long& offset, unsigned int id) const;

	// lnodsはunsigned int にすべきでは？
	//! 要素セグメントの追加
	virtual std::vector<int> AddSegment(std::vector<CElemSeg>& es_ary, 
		const std::vector<int>& lnods );

	// 要素を囲む要素を作る．内部でelsuelのメモリ領域確保はしないので，最初から確保しておく
	bool MakeElemSurElem( const unsigned int& id_es_corner, int* elsuel) const;
private:

	bool MakePointSurElem( 
		const unsigned int id_es, 
        Com::CIndexedArray& elsup ) const;

	bool MakePointSurPoint( 
        const unsigned int id_es, const Com::CIndexedArray& elsup, bool isnt_self,
        Com::CIndexedArray& psup ) const;

    bool MakeElemSurElem( const unsigned int& id_es_corner, Com::CIndexedArray& elsup, int* elsuel) const;
private:


protected:
    unsigned int m_nElem;	//!< 要素の数
    ELEM_TYPE m_ElemType;	//!< 要素タイプ
	unsigned int npoel;		//!< 一要素辺りの説点数
	// コネクティビティはunsigned intにすべきでは？（混合要素への対応のため？）
	int * m_pLnods;			//!< コネクティビティ格納配列
	
	Com::CObjSet<CElemSeg> m_aSeg;	//!< 要素セグメントの集合
};

}
}

#endif
