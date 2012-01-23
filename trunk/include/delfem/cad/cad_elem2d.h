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
@brief Interfaces define the geometry of 2d cad elements
@author Nobuyuki Umetani
*/

#if !defined(CAD_ELEM_2D_H)
#define CAD_ELEM_2D_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <assert.h>
#include <iostream> // needed only in debug

#include "delfem/vector2d.h"

////////////////////////////////////////////////////////////////

namespace Cad{
  
enum CURVE_TYPE {
  CURVE_END_POINT,
  CURVE_LINE,
  CURVE_ARC,
  CURVE_POLYLINE,
  CURVE_BEZIER
};
  
  /*!
@addtogroup CAD
*/
//!@{

//! 2dim loop class
class CLoop2D{
public:
	CLoop2D(const CLoop2D& rhs){
    m_color[0]=rhs.m_color[0];  m_color[1]=rhs.m_color[1];  m_color[2]=rhs.m_color[2];
		ilayer = rhs.ilayer;
  }
	CLoop2D(){
    m_color[0]=0.8; m_color[1]=0.8; m_color[2]=0.8;
		ilayer = 0;
  }  
public:
  double m_color[3];
	unsigned int ilayer;
};
  
/*  
double GetDist_LineSeg_Point(const Com::CVector2D& po_c,
                             const Com::CVector2D& po_s, const Com::CVector2D& po_e);
double GetDist_LineSeg_LineSeg(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,
                               const Com::CVector2D& po_s1, const Com::CVector2D& po_e1);
// line-line intersection detection
bool IsCross_LineSeg_LineSeg(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,
                             const Com::CVector2D& po_s1, const Com::CVector2D& po_e1 );
*/  
//! circle-circle interseciton detection
bool IsCross_Circle_Circle(const Com::CVector2D& po_c0, double radius0,
                           const Com::CVector2D& po_c1, double radius1,
                           Com::CVector2D& po0, Com::CVector2D& po1 );  
/*!
 @brief 円弧と直線の交点を求める
 交点がある場合は２つの交点のposからpoeへのパラメータがt1,t2に入る．
 @retval true 交点がある場合
 @retval false 交点が無い場合
 */  
bool IsCross_Line_Circle(const Com::CVector2D& po_c, const double radius, 
                         const Com::CVector2D& po_s, const Com::CVector2D& po_e, double& t0, double& t1);
//! 点と直線の一番近い点を探す
double FindNearestPointParameter_Line_Point(const Com::CVector2D& po_c,
                                            const Com::CVector2D& po_s, const Com::CVector2D& po_e);    
  
Com::CVector2D GetProjectedPointOnCircle(const Com::CVector2D& c, double r, 
                                         const Com::CVector2D& v);

//! 2dim edge
class CEdge2D{
public:
  CEdge2D(const CEdge2D& rhs) :
  itype(rhs.itype),
  is_left_side(rhs.is_left_side), dist(rhs.dist), 
  aRelCo(rhs.aRelCo),
  aCo(rhs.aCo),  
  id_v_s(rhs.id_v_s), id_v_e(rhs.id_v_e), po_s(rhs.po_s), po_e(rhs.po_e)
  {  
  }
  CEdge2D() : id_v_s(0), id_v_e(0), itype(CURVE_LINE){}
  CEdge2D(unsigned int id_v_s, unsigned int id_v_e) : id_v_s(id_v_s), id_v_e(id_v_e), itype(CURVE_LINE){}
  
  void SetVtxCoords(const Com::CVector2D& ps, const Com::CVector2D& pe) const {
    po_s = ps;
    po_e = pe;
    bb_.isnt_empty = false;
    unsigned int n = aRelCo.size()/2;
    aCo.resize(n);
    if( n > 0 ){
      const Com::CVector2D& gh = po_e - po_s;
      const Com::CVector2D gv(-gh.y, gh.x);
      for(unsigned int i=0;i<n;i++){
        const Com::CVector2D posc = po_s + gh*aRelCo[i*2+0] + gv*aRelCo[i*2+1];
        aCo[i] = posc;
      }      
    }    
  }
  void SetIdVtx(unsigned int id_vs, unsigned int id_ve) const{
    this->id_v_s = id_vs;
    this->id_v_e = id_ve;
  }
  
  // return minimum distance between this edge and e1
  // if the edge obviouly intersects, this function return 0 or -1
  double Distance(const CEdge2D& e1) const;

  // get bounding box of edge
  // lazy evaluation
  // make sure the value is set in po_s, po_e
  const Com::CBoundingBox2D& GetBoundingBox() const{
    if( bb_.isnt_empty ){ return bb_; }
    double xmin,xmax, ymin,ymax;
    this->GetBoundingBox(xmin,xmax, ymin,ymax);
    bb_ = Com::CBoundingBox2D(xmin,xmax, ymin,ymax);
    return bb_;
  }
	
	bool IsCrossEdgeSelf() const;	// check self intersection
	bool IsCrossEdge(const CEdge2D& e1) const;	// intersection between me and e1
	//! 一端が共有された辺同士の交差判定
	bool IsCrossEdge_ShareOnePoint(const CEdge2D& e1, bool is_share_s0, bool is_share_s1) const;
	//! 両端が共有された辺同士の交差判定
	bool IsCrossEdge_ShareBothPoints(const CEdge2D& e1, bool is_share_s1s0) const;

  /*!
   @brief カーブと辺の２頂点を結ぶ直線で囲まれる面積を計算(直線の右側にあれば＋)
   @remarks ループの面積を得るのに使う
   */  
	double AreaEdge() const;

	//! 辺の始点/終点における接線を計算する
	Com::CVector2D GetTangentEdge(bool is_s) const;
	//! 入力点から最も近い辺上の点と距離を返す
	Com::CVector2D GetNearestPoint(const Com::CVector2D& po_in) const;

  // get number of intersection between half line direction (=dir) from point (=org)
  // this function is used for in-out detection
	int NumIntersect_AgainstHalfLine(const Com::CVector2D& org, const Com::CVector2D& dir) const;
  bool GetNearestIntersectionPoint_AgainstHalfLine(Com::CVector2D& sec, const Com::CVector2D& org, const Com::CVector2D& dir) const;
  bool GetCurveAsPolyline(std::vector<Com::CVector2D>& aCo, int ndiv) const;
  double GetCurveLength() const;

  
	////////////////////////////////
	

	// 現在の辺が２つに分割されて，一端がedge_aに入る
	bool Split(Cad::CEdge2D& edge_a, const Com::CVector2D& pa);
	// is_add_aheadはe1がこの辺の前にあるか，is_same_dirはe1がこの辺と同じ向きか
	bool ConnectEdge(const Cad::CEdge2D& e1, bool is_add_ahead, bool is_same_dir);
  
  // get vertex on edge with distance (len) from point v0 along the edge
  // is_front==true:same direction is_front==false:opposite direciton
  bool GetPointOnCurve_OnCircle(const Com::CVector2D& v0, double len, bool is_front,
                                bool& is_exceed, Com::CVector2D& out) const;    
  
  CURVE_TYPE GetCurveType() const { return itype; }
  void SetCurve_Line(){ itype = CURVE_LINE; }
  void SetCurve_Arc(bool is_left_side,double dist){ 
    this->itype = CURVE_ARC;
    this->dist = dist;
    this->is_left_side = is_left_side;
  }
  void SetCurve_Polyline(const std::vector<double>& aRelCo){
    this->itype = CURVE_POLYLINE;
    this->aRelCo = aRelCo;
    aCo.clear();
    const unsigned int n = aRelCo.size()/2;
		Com::CVector2D v0 = po_e-po_s;
		Com::CVector2D v1(-v0.y,v0.x);
		for(unsigned int i=0;i<n;i++){
			Com::CVector2D po0 = po_s + v0*aRelCo[i*2+0] + v1*aRelCo[i*2+1];
      aCo.push_back(po0);
    }      
  }
  void SetCurve_Bezier(double cx0, double cy0, double cx1, double cy1){
    this->itype = CURVE_BEZIER;
    aRelCo.resize(4);
    this->aRelCo[0] = cx0;
    this->aRelCo[1] = cy0;
    this->aRelCo[2] = cx1;
    this->aRelCo[3] = cy1;
    Com::CVector2D v0 = po_e-po_s;
		Com::CVector2D v1(-v0.y,v0.x);
		for(unsigned int i=0;i<2;i++){
			Com::CVector2D po0 = po_s + v0*aRelCo[i*2+0] + v1*aRelCo[i*2+1];
      aCo.push_back(po0);
    }          
  }
  
  
  ////////////////////////////////
  
  // Get Arc property if this curve is not arc, this returns false
  // if is_left_side is true, this arc lies left side from the line connect start and end point of this edge (ID:id_e).
  // The dist means how far is the center of the circle from the line connect start and end point of edge (ID:id_e).  
  void GetCurve_Arc(bool& is_left_side, double& dist) const{
    is_left_side = this->is_left_side;
    dist = this->dist;
  }  
  /*!
   @brief 円が円弧の時、円の中心と半径を計算する
   @remarks 円弧じゃなかったらfalseを返す
   */
	bool GetCenterRadius(Com::CVector2D& po_c, double& radius) const;
  bool GetCenterRadiusThetaLXY(Com::CVector2D& pc, double& radius,
                               double& theta, Com::CVector2D& lx, Com::CVector2D& ly) const;
      
  std::vector<double> GetCurveRelPoint() const { 
    return this->aRelCo;
  }
  void SetCurveRelPoint(const std::vector<double>& aRelCo0){ 
    this->aRelCo = aRelCo0;
    const Com::CVector2D& h = po_e-po_s;
		const Com::CVector2D v(-h.y,h.x);
    unsigned int n = aRelCo.size()/2;
    this->aCo.resize(n);
    for(unsigned int i=0;i<n;i++){
      Com::CVector2D p = po_s + h*aRelCo[i*2+0] + v*aRelCo[i*2+1];
      aCo[i] = p;
    }
  }  
  std::vector<Com::CVector2D> GetCurvePoint() const{
    return this->aCo;
  }
  
  inline unsigned int GetIdVtx(bool is_root) const {
    return is_root ? id_v_s : id_v_e;
  }
  inline Com::CVector2D GetVtxCoord(bool is_root) const {
    return is_root ? po_s : po_e;
  }

private:
  void GetBoundingBox(double& x_min, double& x_max, double& y_min, double& y_max) const;            
	//! 線分と円弧の交錯を判定する
	int NumCross_Arc_LineSeg(const Com::CVector2D& po_s1, const Com::CVector2D& po_e1) const;
	//! 弦と弧で張られる領域内部に点poが入っているかを調べる
	int IsInsideArcSegment(const Com::CVector2D& po) const;
	//! 円弧の中心からみて，点poと円弧が同じ方向に重なっているか？
	int IsDirectionArc(const Com::CVector2D& po) const;
private:
  CURVE_TYPE itype;		//!< 0:Line, 1:Arc, 2:Mesh 3:Bezier
private:
	// type Arc
  bool is_left_side;      //!< is the arc is formed left side of the line po_s, po_e
  double dist;            //!< the minimum distance between line and arc center point
  std::vector<double> aRelCo;
private:
//public:
  //! 干渉チェックの時にだけ一時的に辺の２頂点の座標が代入される
  mutable unsigned int id_v_s, id_v_e;	//!< start vertex
  mutable Com::CVector2D po_s, po_e;
  mutable Com::CBoundingBox2D bb_;
  mutable std::vector<Com::CVector2D> aCo;
};

//! ２次元幾何頂点クラス
class CVertex2D{
public:
	CVertex2D(const Com::CVector2D& point) : point(point){}	
	CVertex2D(const CVertex2D& rhs)
		: point(rhs.point){}	
public:
  Com::CVector2D point;   //!< coordinate
};

/*!
干渉チェックを行う
そのうち交錯位置の情報も返したい
*/
int CheckEdgeIntersection(const std::vector<CEdge2D>& aEdge);

//! @}
}

#endif
