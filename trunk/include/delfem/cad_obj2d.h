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
@brief ２次元ＣＡＤクラス(Cad::CCadObj2Dm)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(CAD_OBJ_2D_M_H)
#define CAD_OBJ_2D_M_H

#include "delfem/vector2d.h"
#include "delfem/serialize.h"
#include "delfem/cad2d_interface.h"
#include "delfem/objset.h"
#include "delfem/cad/brep2d.h"
#include "delfem/cad/cad_elem2d.h"

namespace Cad{

class CVertex2D;
class CLoop2D;
class CEdge2D;
class CTopology;

/*! 
@brief ２次元ＣＡＤモデルクラス
@ingroup CAD
*/
class CCadObj2D : public Cad::ICad2D
{
public:
	//! Iterator which go around loop (面の周りを回るイテレータ)
	class CItrLoop : public Cad::ICad2D::CItrLoop
	{
	public:
		CItrLoop(const CCadObj2D* pCadObj2D, unsigned int id_l)
			: itrl(pCadObj2D->m_BRep.GetItrLoop(id_l)){}
		void Begin(){ itrl.Begin(); }	//! ループの初めの辺に戻す
		void operator++(){ itrl++; } //!< move to next edge
		void operator++(int n){ itrl++; }	//!< dummy operator ( work same as ++)
		//! return current id of edge and return if this edge is counter-clock wise
		bool GetIdEdge(unsigned int& id_e, bool& is_same_dir) const { return itrl.GetIdEdge(id_e,is_same_dir); }
		//! Shift next child loop. if this is end of loop, return false
		bool ShiftChildLoop(){ return itrl.ShiftChildLoop(); }
		bool IsEndChild() const { return itrl.IsEndChild(); }	//!< return true if iterator go around
		unsigned int GetIdVertex() const {        return itrl.GetIdVertex(); }	//!< 現在の頂点を返す
		unsigned int GetIdVertex_Ahead()  const { return itrl.GetIdVertex_Ahead(); }	//!< 一つ先の頂点を返す
        unsigned int GetIdVertex_Behind() const { return itrl.GetIdVertex_Behind(); }	//!< 一つ前の頂点を返す		
		bool IsEnd() const { return itrl.IsEnd(); }	//!< return true if went around this loop
	private:
		CBRep2D::CItrLoop itrl;
	};
	//! iterator go around vertex
	class CItrVertex : public Cad::ICad2D::CItrVertex
	{
	public:		
		CItrVertex(const CCadObj2D* pCadObj2D, unsigned int id_v)
			: itrv(pCadObj2D->m_BRep.GetItrVertex(id_v)){}
		void operator++(){ itrv++; }	//!< go around (cc-wise) loop around vertex 
		void operator++(int n){ itrv++; }	//!< dummy operator (for ++)		
		//! cc-wise ahead  edge-id and its direction(true root of edge is this vertex)
		bool GetIdEdge_Ahead( unsigned int& id_e, bool& is_same_dir) const { return itrv.GetIdEdge_Ahead( id_e,is_same_dir); }
		//! cc-wise behind edge-id and its direction(true root of edge is this vertex)
		bool GetIdEdge_Behind(unsigned int& id_e, bool& is_same_dir) const { return itrv.GetIdEdge_Behind(id_e,is_same_dir); }
		unsigned int GetIdLoop() const{ return itrv.GetIdLoop(); } //!< get loop-id		
		bool IsEnd() const { return itrv.IsEnd(); }	//!< return true if iterator go around
	private:
		CBRep2D::CItrVertex itrv;
	};
    friend class CItrLoop;
    friend class CItrVertex;
	//! デフォルトコンストラクタ
	CCadObj2D();
	//! デストラクタ
	virtual ~CCadObj2D();
	//! 全ての要素を削除して初期化
	void Clear();

	////////////////////////////////
	// Getメソッド

	//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
    virtual std::auto_ptr<Cad::ICad2D::CItrLoop> GetItrLoop(unsigned int id_l) const{
        return std::auto_ptr<Cad::ICad2D::CItrLoop>( new CItrLoop(this,id_l) );	// インスタンスを生成
	}
	//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
    virtual std::auto_ptr<Cad::ICad2D::CItrVertex> GetItrVertex(unsigned int id_v) const{
        return std::auto_ptr<Cad::ICad2D::CItrVertex>( new CItrVertex(this,id_v) );	// インスタンスを生成
	}
	////////////////////////////////////////////////
	// 構成要素へのアクセス
	virtual bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const;
	virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE itype) const;

	////////////////////////////////
	// レイヤ関係の関数
    virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const;
	virtual void GetLayerMinMax(int& layer_min, int& layer_max) const;

	////////////////////////////////
	// ループのメンバ関数

	//! @{
    //! ID:id_lのループの色を返す
    virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const;
    //! ID:id_lのループの色を設定する
    virtual bool SetColor_Loop(unsigned int id_l, const double color[3] );
	//! ID:id_lのループの面積を返えす
	virtual double GetArea_Loop(unsigned int id_l) const;
	bool ShiftLayer_Loop(unsigned int id_l, bool is_up);
	//! @}

	////////////////////////////////
	// 辺のメンバ関数
	// 　　円弧について
	// 　　　　is_left_side==trueならこの円弧は始点終点を結ぶ線の始点から見て左側にある．is_left_side==falseなら右側
	// 　　　　distは円の中心から始点終点を結ぶ線の中点が左側にどれだけ離れているかという値

	//! @{
	//! 辺の始点と終点の頂点のIDを返す
	virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const;
    //! 辺の両側のループのIDを返す
    virtual bool GetIdLoop_Edge(unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const;
	/*!
	@brief 辺の形状タイプを返す
    @retval 0 線分
    @retval 1 円弧
    @retval 2 メッシュ
	*/
	virtual int GetEdgeCurveType(const unsigned int& id_e) const;
	//! ID:id_eの辺の情報を得る
	virtual bool GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const;
	virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const;
    virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen = -1) const;    
    //! ID:id_eの辺をndiv個に分割したものを得る. 但し始点，終点はそれぞれps,peとする
    virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
		unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const;
	//! @}

    bool GetPointOnCurve_OnCircle(unsigned int id_e,
                                  const Com::CVector2D& v0, double len, bool is_front,
                                  bool& is_exceed, Com::CVector2D& out) const;
    //! 入力点から最も近い辺上の点と距離を返す
    Com::CVector2D GetNearestPoint(unsigned int id_e, const Com::CVector2D& po_in) const;

	////////////////////////////////
	// 頂点のメンバ関数

	//! @{
	/*!
	@brief 頂点の座標を得る
	@param[in] id_v 頂点のID
	@retval 頂点の座標
	*/
	virtual const Com::CVector2D& GetVertexCoord(unsigned int id_v) const;
	//! @}

	////////////////////////////////////////////////
	// 形状操作関数(トポロジーを変える)
	
	//! @{ 
	/*!
	@brief 頂点集合から面を作る
	@param[in] vec_ary 左回りの頂点の配列
	@param[in] id_l 多角形を追加する面のID, 省略された場合は，領域の外部に面が作られる．
	@retval 新しくできた面ＩＤ 成功した場合
	@retval ０ 失敗した場合
	*/
	unsigned int AddPolygon( const std::vector<Com::CVector2D>& vec_ary, unsigned int id_l = 0);

	/*!
	@brief 頂点を加える関数
	@param[in] itype 頂点を加える先の要素のタイプ
	@param[in] id 頂点を加える先の要素のID
	@param[in] vec 頂点の座標
	@retval 新しくできた頂点のID 成功した場合
	@retval ０ 失敗した場合
	@remarks itype == Cad::NOT_SETだったら，形状の外側に頂点を追加
	*/
	unsigned int AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id, const Com::CVector2D& vec);
	
	/*!
	@brief 要素を消去する関数
	@param[in] itype 要素のタイプ
	@param[in] id 要素のID
	@retval 成功したかどうか
	@pre 辺を消去する時は辺の両側が面でなければならない
	*/
	bool RemoveElement(Cad::CAD_ELEM_TYPE itype, unsigned int id);
	/*!
	@brief ２つの頂点(ID:id_v1,id_v2)からEdgeを作る、
	@retval 新しくできた辺のＩＤ 成功した場合
	@retval ０ 失敗した場合
	*/
	unsigned int ConnectVertex(CEdge2D edge);
	/*!
	@brief ２つの頂点(ID:id_v1,id_v2)からEdgeを作る、
	@retval 新しくできた辺のＩＤ 成功した場合
	@retval ０ 失敗した場合
	*/
	unsigned int ConnectVertex_Line(unsigned int id_v1, unsigned int id_v2){
		Cad::CEdge2D e(id_v1,id_v2, 0, 0, 0 );
		return this->ConnectVertex(e);
	}
	//! @}

	////////////////////////////////////////////////
	// 形状操作関数（トポロジーを変えない）

    bool SetCurve_Polyline(const unsigned int id_e);
	//! ID:id_eの辺をメッシュにセットする
	bool SetCurve_Polyline(const unsigned int id_e, const std::vector<Com::CVector2D>& aVec);
	//! ID:id_eの辺を円弧にセットする
    bool SetCurve_Arc(const unsigned int id_e, bool is_left_side=false, double rdist=10.0);
    //! ID:id_eの辺を直線にセットする
	bool SetCurve_Line(const unsigned int id_e);

	////////////////////////////////////////////////
	// IOルーティン

	//! DXFへ書き出し
	bool WriteToFile_dxf(const std::string& file_name) const;
	//! 読み込み書き出し
	bool Serialize( Com::CSerializer& serialize );	
protected:
	// id_vs, id_ve, po_s, po_eが代入されたEdgeを返す
	const CEdge2D& GetEdge(unsigned int id_e) const;
	// id_vs, id_ve, po_s, po_eが代入されたEdgeを返す
	CEdge2D& GetEdge(unsigned int id_e);
	int AssertValid() const;
	bool CheckIsPointInside_ItrLoop(ICad2D::CItrLoop& itrl, const Com::CVector2D& p1) const;

//    返り値が０ならid_ul1の中の点は全てid_ul2の中
//    返り値が１ならid_ul1の中の点の一部はid_ul2の中で一部は外
//    返り値が２ならid_ul1の中の点は全てid_ul2の外
	unsigned int CheckInOut_ItrLoopPoint_ItrLoop(ICad2D::CItrLoop& itrl1, ICad2D::CItrLoop& itrl2) const;
	
	bool CheckIsPointInsideLoop(unsigned int id_l1, const Com::CVector2D& point) const;
//    ループが健全かどうか調べる．０なら健全
    int CheckLoop(unsigned int id_l) const;

	// 返り値は半辺のIDでid_vとpointを結ぶ半直線からid_vを中心に時計周りに最初に出会う半辺
	// この半辺の属するループはこの半直線と重なっている。
	// id_vが浮遊点の場合は浮遊点周りのループが帰る
	CBRep2D::CItrVertex FindCorner_HalfLine(unsigned int id_v, const Com::CVector2D& point) const;
	bool CheckLoopIntersection(unsigned int id_l=0) const;
	double GetArea_ItrLoop(ICad2D::CItrLoop& itrl) const;
protected:
	////////////////
	Com::CObjSet<CLoop2D>   m_LoopSet;
	Com::CObjSet<CEdge2D>   m_EdgeSet;
	Com::CObjSet<CVertex2D> m_VertexSet;
	////////////////
	CBRep2D m_BRep;	// 位相情報を持っているクラス
};

}	// end namespace CAD

#endif
