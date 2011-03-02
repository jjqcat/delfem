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
@brief interface of 2D cad class (Cad::CCadObj2Dm)
@author Nobuyuki Umetani
*/

#if !defined(CAD_OBJ_2D_H)
#define CAD_OBJ_2D_H

#include <vector>

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
@brief 2 dimentional cad model class
@ingroup CAD
*/
class CCadObj2D : public Cad::ICad2D_Msh
{
public:     
  class CResAddVertex{
  public:
    CResAddVertex(){
      id_v_add = 0;
      id_e_add = 0;
    }
  public:
    unsigned int id_v_add;
    unsigned int id_e_add;
  };
  class CResAddPolygon{
  public:
    CResAddPolygon(){
      id_l_add = 0;
    }
    CResAddPolygon(const CResAddPolygon& lhs){
      this->id_l_add = lhs.id_l_add;
      this->aIdV = lhs.aIdV;
      this->aIdE = lhs.aIdE;
    }
  public:
    unsigned int id_l_add;
    std::vector<unsigned int> aIdV;
    std::vector<unsigned int> aIdE;
  };
  
  
  ////////////////////////////////
	// constructor & destructor
  
	//! default constructor
	CCadObj2D();
  //! copy constructor
  CCadObj2D(const CCadObj2D& cad);
	//! destructor
	virtual ~CCadObj2D();
	//! initialization clear all element
	void Clear();

	////////////////////////////////
	// Get method

	//! function gives iterator which travel vtx and edge inside the loop (ID:id_l)
  virtual std::auto_ptr<Cad::IItrLoop> GetPtrItrLoop(unsigned int id_l) const{
    return std::auto_ptr<Cad::IItrLoop>( new CBRepSurface::CItrLoop(m_BRep,id_l) );	// instance
	}
	virtual bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const;
	virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE itype) const;  
  virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const;
  //! 辺の両側のループのIDを返す
  virtual bool GetIdLoop_Edge(unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const;    
  CBRepSurface::CItrVertex GetItrVertex(unsigned int id_v) const { return CBRepSurface::CItrVertex(m_BRep,id_v); }
  CBRepSurface::CItrLoop GetItrLoop(  unsigned int id_l) const { return CBRepSurface::CItrLoop(m_BRep,id_l); }
    
	// functions related to layer
  virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const;
	virtual void GetLayerMinMax(int& layer_min, int& layer_max) const;
  bool ShiftLayer_Loop(unsigned int id_l, bool is_up);
  
  double GetMinClearance() const { return min_clearance; }
  
	// loop function  
	//! @{
	bool CheckIsPointInsideLoop(unsigned int id_l1, const Com::CVector2D& point) const;    
  double SignedDistPointLoop(unsigned int id_l1, const Com::CVector2D& point, unsigned int id_v_ignore=0) const;  
  //! get color(double[3]) of loop(ID:id_l), return false if there is no loop(ID:id_l)
  virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const;
  //! ID:id_l set color of loop
  virtual bool SetColor_Loop(unsigned int id_l, const double color[3] );
	//! ID:id_l return are of the loop
	virtual double GetArea_Loop(unsigned int id_l) const;
	//! @}
  	////////////////////////////////
	// 辺のメンバ関数
	// 　　円弧について
	// 　　　　is_left_side==trueならこの円弧は始点終点を結ぶ線の始点から見て左側にある．is_left_side==falseなら右側
	// 　　　　distは円の中心から始点終点を結ぶ線の中点が左側にどれだけ離れているかという値

	//! 辺の始点と終点の頂点のIDを返す
	// return edge with vertex id(id_vs, id_ve) and vertex coord
	const CEdge2D& GetEdge(unsigned int id_e) const;  
	/*!
	@brief Get Geometric Type of the Curve
    @retval 0 Line
    @retval 1 Arc
    @retval 2 Polyline
	*/
	virtual int GetEdgeCurveType(const unsigned int& id_e) const;
	//! get information of edge(ID:id_e)
	virtual bool GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const;
	virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const;
  
  virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen = -1) const;    
  //! ID:id_eの辺をndiv個に分割したものを得る. 但し始点，終点はそれぞれps,peとする
  virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
                                 unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const;

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
	virtual Com::CVector2D GetVertexCoord(unsigned int id_v) const;
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
	CResAddPolygon AddPolygon(const std::vector<Com::CVector2D>& vec_ary, unsigned int id_l = 0);

	/*!
	@brief 頂点を加える関数
	@param[in] itype 頂点を加える先の要素のタイプ
	@param[in] id 頂点を加える先の要素のID
	@param[in] vec 頂点の座標
	@retval 新しくできた頂点のID 成功した場合
	@retval ０ 失敗した場合
	@remarks itype == Cad::NOT_SETだったら，形状の外側に頂点を追加
	*/
  CResAddVertex AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id, const Com::CVector2D& vec);
	
	/*!
	@brief 要素を消去する関数
	@param[in] itype 要素のタイプ
	@param[in] id 要素のID
	@retval 成功したかどうか
	@pre 辺を消去する時は辺の両側が面でなければならない
	*/
	bool RemoveElement(Cad::CAD_ELEM_TYPE itype, unsigned int id);

	CBRepSurface::CResConnectVertex ConnectVertex(CEdge2D edge);
  CBRepSurface::CResConnectVertex ConnectVertex_Line(unsigned int id_v1, unsigned int id_v2){
		Cad::CEdge2D e(id_v1,id_v2, 0, 0, 0 );
		return this->ConnectVertex(e);
	}
	//! @}

	////////////////////////////////////////////////
	// 形状操作関数（トポロジーを変えない）

  bool SetCurve_Polyline(const unsigned int id_e);
	//! set edge (ID:id_e) mesh
	bool SetCurve_Polyline(unsigned int id_e, const std::vector<Com::CVector2D>& aVec);
	//! set edge (ID:id_e) arc 
  bool SetCurve_Arc(const unsigned int id_e, bool is_left_side=false, double rdist=10.0);
  //! set edge (ID:id_e) straight line
	bool SetCurve_Line(const unsigned int id_e);

	////////////////////////////////////////////////
	// IO routines
	
	virtual bool WriteToFile_dxf(const std::string& file_name, double scale) const;
	bool Serialize( Com::CSerializer& serialize );		
protected:
	// return edge with vertex id(id_vs, id_ve) and vertex coord  
	CEdge2D& GetEdgeRef(unsigned int id_e);  
	int AssertValid() const;
	bool CheckIsPointInside_ItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& p1) const;
  double DistPointItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& point) const;
  
  
  // ret:0 all the points in id_ul1 are inside id_ul2
  // ret:1 there are both inside/outside id_ul1 points in id_ul1 or ambiguous
  // ret:2 all the points in id_ul1 are outside id_ul2
	unsigned int CheckInOut_ItrLoopPoint_ItrLoop(CBRepSurface::CItrLoop& itrl1, CBRepSurface::CItrLoop& itrl2) const;
	
	// assert loop : return 0 if loop is OK
  int CheckLoop(unsigned int id_l) const;

	// 返り値は半辺のIDでid_vとpointを結ぶ半直線からid_vを中心に時計周りに最初に出会う半辺
	// この半辺の属するループはこの半直線と重なっている。
	// id_vが浮遊点の場合は浮遊点周りのループが帰る
	CBRepSurface::CItrVertex FindCorner_HalfLine(unsigned int id_v, const Com::CVector2D& point) const;
	bool CheckIntersection_Loop(unsigned int id_l=0) const;
  bool CheckIntersection_EdgeAgainstLoop(const CEdge2D& edge,unsigned int id_l=0) const;
	double GetArea_ItrLoop(CBRepSurface::CItrLoop& itrl) const;
protected:
	////////////////
	Com::CObjSet<CLoop2D>   m_LoopSet;
	Com::CObjSet<CEdge2D>   m_EdgeSet;
	Com::CObjSet<CVertex2D> m_VertexSet;
	////////////////
	CBRepSurface m_BRep;	// class which have topology
  double min_clearance;
};

}	// end namespace CAD

#endif
