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
@brief 抽象２次元ＣＡＤクラス(Cad::CCadObj2D)のインターフェース
@remarks このインターフェースさえ継承できれば、Meshを切ったりCDrawerCADで描画できる
@author Nobuyuki Umetani
*/

#if !defined(CAD_2D_INTERFACE_H)
#define CAD_2D_INTERFACE_H

#include <memory>   // autoptrのために必要

#include "delfem/vector2d.h"
#include "delfem/cad_com.h"

namespace Cad{

/*! 
@ingroup CAD
@brief 2D CAD model class (２次元ＣＡＤモデルクラス)
*/
class ICad2D
{
public:
    //! iterator go around loop
	class CItrLoop
	{
	public:
		virtual void Begin() = 0;	//!< back to initial point of current use-loop
		virtual void operator++() = 0; //!< move to next edge
		virtual void operator++(int n)= 0;	//!< dummy operator (for ++)
		//! return current edge id and whether if this edge is same dirrection as loop
		virtual bool GetIdEdge(unsigned int& id_e, bool& is_same_dir) const = 0;	
		virtual bool ShiftChildLoop() = 0;	//!< move to next use-loop in this loop
		virtual bool IsEndChild() const = 0;	//!< return true if iterator go around
		virtual unsigned int GetIdVertex() const = 0;	//!< return current vertex id
		virtual unsigned int GetIdVertex_Ahead()  const = 0;	//!< return next vertex
		virtual unsigned int GetIdVertex_Behind() const = 0;	//!< return previous vertex
		virtual bool IsEnd() const = 0;	//!< return true if iterator go around
	};
	//! iterator go around vertex
	class CItrVertex
	{ 
	public:		
		virtual void operator++() = 0;	//!< go around (cc-wise) loop around vertex 
		virtual void operator++(int n) = 0;	//!< dummy operator (for ++)		
		//! cc-wise ahead  edge-id and its direction(true root of edge is this vertex)
		virtual bool GetIdEdge_Ahead(unsigned int& id_e, bool& is_same_dir) const = 0;
		//! cc-wise behind edge-id and its direction(true root of edge is this vertex)
		virtual bool GetIdEdge_Behind(unsigned int& id_e, bool& is_same_dir) const = 0;
		virtual unsigned int GetIdLoop() const = 0; //!< get loop-id		
		virtual bool IsEnd() const = 0;	//!< return true if iterator go around
	};
public:	
    ICad2D(){}	//!< needs defalut constructor
    virtual ~ICad2D(){}	//!< virtual destructor is must for interface class

	////////////////////////////////
	// Getメソッド

  //! 辺の両側のループのIDを返す
  virtual bool GetIdLoop_Edge(unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const = 0;
	//! 辺の始点と終点の頂点のIDを返す
	virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const = 0;  

	////////////////////////////////////////////////
	// 構成要素へのアクセス

	//!@{
	//! idが使われているかどうかを調べる関数
	virtual bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const = 0;
	//! すべてのIDを配列にして返す関数
	virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE) const = 0;
	//!@}
	
	// レイヤ関係の関数
    virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const = 0;
	virtual void GetLayerMinMax(int& layer_min, int& layer_max) const = 0;

	////////////////////////////////
	// ループのメンバ関数

	//! @{
    //! ID:id_lのループの色を返す(本来このクラスは位相と幾何情報以外を持つべきではないかもしれないので暫定的)
    virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const = 0;
	//! ID:id_lのループの面積を返えす
	virtual double GetArea_Loop(unsigned int id_l) const = 0;
	//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
    virtual std::auto_ptr<CItrLoop> GetItrLoop(unsigned int id_l) const = 0;
	//! ID:id_vの頂点を構成する辺やループをめぐるイテレータを返す関数
    virtual std::auto_ptr<CItrVertex> GetItrVertex(unsigned int id_v) const = 0;
	//! @}

	////////////////////////////////
	// 辺のメンバ関数
	// 　　円弧について
	// 　　　　is_left_side==trueならこの円弧は始点終点を結ぶ線の始点から見て左側にある．is_left_side==falseなら右側
	// 　　　　distは円の中心から始点終点を結ぶ線の中点が左側にどれだけ離れているかという値
	//! @{
    //ここから省略できそう？(このinterfaceは描画とメッシュ切のためだとして)
    /*!
	@brief 辺の形状タイプを返す
	@retval ０ 線分
    @retval １ 円弧
    @retval ２ メッシュ
	*/
	virtual int GetEdgeCurveType(const unsigned int& id_e) const = 0;
	//! ID:id_eの辺の情報を得る
	virtual bool GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const = 0;
	virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const = 0;
    //ここまで省略できそう？

    //! ID:id_eの辺のメッシュ分割を得る(elen<=0ならできるだけ詳細にメッシュを切ろうとする)
    virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen) const = 0;
    //! ID:id_eの辺をndiv個に分割したものを得る. 但し始点，終点はそれぞれps,peとする
    virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
                               unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const = 0;
	//! @}

	////////////////////////////////
	// 頂点のメンバ関数

	//! @{
	/*!
	@brief 頂点の座標を得る
	@param[in] id_v 頂点のID
	@retval 頂点の座標
	*/
	virtual const Com::CVector2D& GetVertexCoord(unsigned int id_v) const = 0;
	//! @}
};

}	// end namespace CAD

#endif
