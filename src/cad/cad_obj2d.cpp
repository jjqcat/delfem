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
// CadObj2D.cpp : ２次元ＣＡＤモデルクラス(CCadObj2D)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif

#define for if(0);else for

#include <iostream>
#include <set>
#include <map>
#include <vector>	
#include <cassert>	
#include <math.h>	
#include <cstring>	// strlen

#include "delfem/cad_obj2d.h"
#include "delfem/cad/cad_elem2d.h"

using namespace Cad;
using namespace Com;

CCadObj2D::CCadObj2D()
{
}

CCadObj2D::~CCadObj2D()
{
	this->Clear();
}

void CCadObj2D::Clear()
{
	this->m_LoopSet.Clear();
	this->m_EdgeSet.Clear();
	this->m_VertexSet.Clear();
	this->m_BRep.Clear();
}

////////////////////////////////////////////////////////////////
 
bool CCadObj2D::IsElemID(Cad::CAD_ELEM_TYPE itype,unsigned int id) const
{
    if(      itype == Cad::NOT_SET ){ return false; }
    else if( itype == Cad::VERTEX  ){ return m_VertexSet.IsObjID(id); }
    else if( itype == Cad::EDGE    ){ return m_EdgeSet.IsObjID(id);   }
    else if( itype == Cad::LOOP    ){ return m_LoopSet.IsObjID(id);   }
	else{ assert(0); }
	return false;
}

const std::vector<unsigned int> CCadObj2D::GetAryElemID(Cad::CAD_ELEM_TYPE itype) const
{
	if(      itype == Cad::VERTEX ){ return m_VertexSet.GetAry_ObjID(); }
	else if( itype == Cad::EDGE   ){ return m_EdgeSet.GetAry_ObjID();   }
	else if( itype == Cad::LOOP   ){ return m_LoopSet.GetAry_ObjID();   }
	assert(0);
	std::vector<unsigned int> null_vec;
	return null_vec;
}

const CVector2D& CCadObj2D::GetVertexCoord(unsigned int id_v) const
{
	assert( m_VertexSet.IsObjID(id_v) );
	const CVertex2D& v = m_VertexSet.GetObj(id_v);
	return v.point;
}

const CEdge2D& CCadObj2D::GetEdge(unsigned int id_e) const
{
	assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	assert( this->m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	unsigned int id_vs, id_ve;
	m_BRep.GetIdVertex_Edge(id_e, id_vs, id_ve);
	e.id_v_s = id_vs;
	e.id_v_e = id_ve;
	assert( m_BRep.IsElemID(Cad::VERTEX,id_vs) );
	assert( m_BRep.IsElemID(Cad::VERTEX,id_ve) );
	assert( m_VertexSet.IsObjID(id_vs) );
	assert( m_VertexSet.IsObjID(id_ve) );
	e.po_s = this->GetVertexCoord(id_vs);
	e.po_e = this->GetVertexCoord(id_ve);
	return e;
}

CEdge2D& CCadObj2D::GetEdge(unsigned int id_e)
{
	assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	assert( this->m_EdgeSet.IsObjID(id_e) );
	CEdge2D& e = m_EdgeSet.GetObj(id_e);
	unsigned int id_vs, id_ve;
	m_BRep.GetIdVertex_Edge(id_e, id_vs, id_ve);
	e.id_v_s = id_vs;
	e.id_v_e = id_ve;
	assert( m_BRep.IsElemID(Cad::VERTEX,id_vs) );
	assert( m_BRep.IsElemID(Cad::VERTEX,id_ve) );
	assert( m_VertexSet.IsObjID(id_vs) );
	assert( m_VertexSet.IsObjID(id_ve) );
	e.po_s = this->GetVertexCoord(id_vs);
	e.po_e = this->GetVertexCoord(id_ve);
	return e;
}

bool CCadObj2D::GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const
{
	m_BRep.IsElemID(Cad::EDGE,id_e);
	return m_BRep.GetIdVertex_Edge(id_e,id_v_s,id_v_e);
}

bool CCadObj2D::GetIdLoop_Edge(unsigned int& id_l_l, unsigned int& id_l_r, unsigned int id_e) const
{
	return this->m_BRep.GetIdLoop_Edge(id_e,id_l_l,id_l_r);
}

int CCadObj2D::GetEdgeCurveType(const unsigned int& id_e) const
{
	assert( m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	return e.itype;
}

bool CCadObj2D::GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const
{
	if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	is_left_side = e.is_left_side;
	dist = e.dist;
	return true;
}

bool CCadObj2D::GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const
{
	if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	aRelCoMesh = e.aRelCoMesh;
	return true;
}

bool CCadObj2D::GetCurve_Polyline(
        unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen) const
{
    if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
        return false;
    }
    const CEdge2D& e = this->GetEdge(id_e);
    double len = e.GetCurveLength();
    if( elen > 0 ){
        const unsigned int ndiv = (unsigned int)(len/elen)+1;
        return e.GetCurve_Mesh(aCo, ndiv);
    }
    return e.GetCurve_Mesh(aCo,-1);
}

//! ID:id_eの辺をndiv個に分割したものを得る. 但し始点，終点はそれぞれps,peとする
bool CCadObj2D::GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
        unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const
{
    if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
        return false;
    }
    const CEdge2D& e = m_EdgeSet.GetObj(id_e);
    e.po_s = ps;    e.po_e = pe;
    return e.GetCurve_Mesh(aCo,ndiv);
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

double CCadObj2D::GetArea_ItrLoop(ICad2D::CItrLoop& itrl) const
{
	double area = 0.0;
	for(itrl.Begin();!itrl.IsEnd();itrl++){
		unsigned int id_e;   bool is_same_dir;
		itrl.GetIdEdge(id_e,is_same_dir);
		assert( this->IsElemID(Cad::EDGE,id_e) );
		const CEdge2D& e = this->GetEdge(id_e);
		assert( ((is_same_dir) ? e.id_v_s : e.id_v_e) == itrl.GetIdVertex() );
		assert( ((is_same_dir) ? e.id_v_e : e.id_v_s) == itrl.GetIdVertex_Ahead() );
		// ここから辺の面積を足し合わせ
		const double earea = TriArea(e.po_s, e.po_e, CVector2D(0.0,0.0)) + e.AreaEdge();
		if( is_same_dir ){ area += earea; }
		else{              area -= earea; }
	}
	return area;
}

// ループの面積を返す
double CCadObj2D::GetArea_Loop(unsigned int id_l) const
{
	assert( m_LoopSet.IsObjID(id_l) );
	double area = 0.0;
    for(CBRep2D::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
		area += this->GetArea_ItrLoop(itrl);
	}
	return area;
}

// ID:id_lのループの色を返す
bool CCadObj2D::GetColor_Loop(unsigned int id_l, double color[3] ) const
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    const CLoop2D& l = m_LoopSet.GetObj(id_l);
    color[0] = l.m_color[0];
    color[1] = l.m_color[1];
    color[2] = l.m_color[2];
    return true;
}

int CCadObj2D::GetLayer_Loop(unsigned int id_l ) const
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    const CLoop2D& l = m_LoopSet.GetObj(id_l);
	return l.ilayer;
}

// ID:id_lのループの色を設定する
bool CCadObj2D::SetColor_Loop(unsigned int id_l, const double color[3] )
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    CLoop2D& l = m_LoopSet.GetObj(id_l);
    l.m_color[0] = color[0];
    l.m_color[1] = color[1];
    l.m_color[2] = color[2];
    return true;
}

bool CCadObj2D::ShiftLayer(unsigned int id_l, bool is_up)
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    CLoop2D& l = m_LoopSet.GetObj(id_l);
	if( is_up ){
		l.ilayer++;
	}
	else{
		l.ilayer--;
	}
	return true;
}

bool CCadObj2D::CheckIsPointInside_ItrLoop(ICad2D::CItrLoop& itrl, const CVector2D& point) const
{
//	std::cout << "CCadObj2D::CheckIsPointInside_ItrLoop" << std::endl;
	for(unsigned int i=1;i<29;i++){	// 29はできるだけ手ごろな素数
		unsigned int cross_counter = 0;
		bool iflg = true;
		CVector2D dir(sin(3.14*i/29.0),cos(3.14*i/29.0));
		for(itrl.Begin();!itrl.IsEnd();itrl++){
			unsigned int id_e;   bool is_same_dir;
			itrl.GetIdEdge(id_e,is_same_dir);
			if( id_e == 0 ) return false;
			assert( m_EdgeSet.IsObjID(id_e) );
			const CEdge2D& e = this->GetEdge(id_e);
			const int ires = e.NumIntersect_AgainstHalfLine( point,dir );
			// -1なら微妙なのでやり直し
			if( ires == -1 ){ iflg = false; break; }
			cross_counter += ires;
		}
		if( iflg == true ){ 
			if( cross_counter % 2 == 0 ) return false;
			return true;
		}
	}
	assert(0);	// ここまで来ることがありませんように！
	return false;
}


// 返り値が０ならid_ul1の中の点は全てid_ul2の中
// 返り値が１ならid_ul1の中の点の一部はid_ul2の中で一部は外
// 返り値が２ならid_ul1の中の点は全てid_ul2の外
unsigned int CCadObj2D::CheckInOut_ItrLoopPoint_ItrLoop(
	ICad2D::CItrLoop& itrl1, ICad2D::CItrLoop& itrl2) const
{
//	std::cout << "CheckInOut_ItrLoopPoint_ItrLoop " << std::endl;
	unsigned int count_out=0, count_in=0;
	for(itrl1.Begin();!itrl1.IsEnd();itrl1++){
		const unsigned int id_v = itrl1.GetIdVertex();
		const CVertex2D& v = this->m_VertexSet.GetObj(id_v);
		if( this->CheckIsPointInside_ItrLoop(itrl2, v.point) ){
			if( count_out!=0 ) return 1;
			count_in++;
		}
		else{
			if( count_in!=0 ) return 1;
			count_out++;
		}
	}
	if( count_in == 0 ){
		assert( count_out != 0 );
		return 2;
	}	
	assert( count_out == 0 );
	assert( count_in != 0 );
	return 0;
}

bool CCadObj2D::CheckIsPointInsideLoop(unsigned int id_l1, const CVector2D& point) const
{
	assert( m_LoopSet.IsObjID(id_l1) );
	CBRep2D::CItrLoop itrl = m_BRep.GetItrLoop(id_l1);
	if( !CheckIsPointInside_ItrLoop(itrl,point) ) return false;
	for(;itrl.IsEndChild();itrl.ShiftChildLoop()){
		if( CheckIsPointInside_ItrLoop(itrl,point) ) return false;
	}
	return true;
}

int CCadObj2D::CheckLoop(unsigned int id_l) const
{
//    std::cout << "Check Loop " << id_l << std::endl;
    {	// ループの自己干渉を調べる
		std::vector<Cad::CEdge2D> aEdge;
		if( this->CheckLoopIntersection(aEdge,id_l) ){ return 1; }
	}
	{	// 親ループの面積が正で子ループの面積が負であることを調べる
		for(CBRep2D::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			if( itrl.IsParent() ){
				if( itrl.GetType() != 2 ) return 2;
				if( this->GetArea_ItrLoop(itrl) < 0 ) return 2;
			}
			else if( itrl.GetType() == 2 ){
				if( this->GetArea_ItrLoop(itrl) > 0 ) return 2;
			}
		}
	}
	{	// 子ループ内の点が親ループの中に入っているかどうか調べる
        CBRep2D::CItrLoop itrl_p = m_BRep.GetItrLoop(id_l);
		for(CBRep2D::CItrLoop itrl_c=m_BRep.GetItrLoop(id_l);!itrl_c.IsEndChild();itrl_c.ShiftChildLoop()){
			if( itrl_c.IsParent() ) continue;
            if( this->CheckInOut_ItrLoopPoint_ItrLoop(itrl_c,itrl_p)!=0 ) return 3;
		}
	}	
	{   // 子ループ同士が互いにお互いの頂点を含まないかを調べる
		for(CBRep2D::CItrLoop itrl1=m_BRep.GetItrLoop(id_l);!itrl1.IsEndChild();itrl1.ShiftChildLoop()){		
			if( itrl1.IsParent() ){ continue; }
			for(CBRep2D::CItrLoop itrl2=m_BRep.GetItrLoop(id_l);!itrl2.IsEndChild();itrl2.ShiftChildLoop()){
				if( itrl2.IsParent() ){ continue; }
				if( itrl1.IsSameUseLoop(itrl2) ){ continue; }
				if( this->CheckInOut_ItrLoopPoint_ItrLoop(itrl1,itrl2)!=2 ) return 4;
			}
		}
	}
	{	// ループのの中の各角についてその角の中に他の点が入っていないかを確認する
		Com::CVector2D vec_zero(0,0);
		for(CBRep2D::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			for(itrl.Begin();!itrl.IsEnd();itrl++){	// ループの中の点をめぐる
				unsigned int id_vm, id_vf;
				Com::CVector2D dir;
				{
					id_vm = itrl.GetIdVertex();
					id_vf = itrl.GetIdVertex_Ahead();
					unsigned int id_e0;   bool is_same_dir0;
					itrl.GetIdEdge(id_e0,is_same_dir0);
					if( !this->m_EdgeSet.IsObjID(id_e0) ) continue;
					const CEdge2D& e0 = this->GetEdge(id_e0);
					dir = e0.GetTangentEdge(is_same_dir0);
					assert( ((is_same_dir0) ? e0.id_v_s : e0.id_v_e) == id_vm );
					assert( ((is_same_dir0) ? e0.id_v_e : e0.id_v_s) == id_vf );
				}
				for(CBRep2D::CItrVertex itrv=m_BRep.GetItrVertex(id_vm);!itrv.IsEnd();itrv++){	// 点周りの辺をめぐる
					unsigned int id_e0;   bool is_same_dir0;
					itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
					if( !this->m_EdgeSet.IsObjID(id_e0) ){ continue; }
					const Cad::CEdge2D& e0 = this->GetEdge(id_e0);
					assert( ((is_same_dir0) ? e0.id_v_s : e0.id_v_e) == id_vm );
					if( ((is_same_dir0) ? e0.id_v_e : e0.id_v_s) == id_vf ){ continue; }
					const CVector2D& tan0 = e0.GetTangentEdge(is_same_dir0);
					////////////////
					unsigned int id_e1;   bool is_same_dir1;
					itrv.GetIdEdge_Ahead(id_e1,is_same_dir1);
					if( !this->m_EdgeSet.IsObjID(id_e1) ){ continue; }
					const Cad::CEdge2D& e1 = this->GetEdge(id_e1);
					assert( ((is_same_dir1) ? e1.id_v_s : e1.id_v_e) == id_vm );
					if( ((is_same_dir1) ? e1.id_v_e : e1.id_v_s) == id_vf ){ continue; }
					const CVector2D& tan1 = e1.GetTangentEdge(is_same_dir1);
					////////////////
					const double area0 = TriArea(tan1,vec_zero,tan0);
					const double area1 = TriArea(tan1,vec_zero,dir);	
					const double area2 = TriArea(dir, vec_zero,tan0);
					if( (area0 > 0.0 && area1 > 0.0 && area2 > 0.0) || 
						(area0 < 0.0 && (area1 > 0.0 || area2 > 0.0) ) ){ return 5; }
				}
			}
		}
	}
    return 0;
}

int CCadObj2D::AssertValid() const 
{
	{	// ループについてチェックする
		std::vector<unsigned int> id_l_ary = m_LoopSet.GetAry_ObjID();
		for(unsigned iid_l=0;iid_l<id_l_ary.size();iid_l++){
			const unsigned int id_l = id_l_ary[iid_l];
            int res = CheckLoop(id_l);
            if( res!=0 ){
                if(      res == 1 ){ std::cout << "Intersectoin in the loop" << std::endl; }
				else if( res == 2 ){ std::cout << "Check area parent plus, childe minus" << std::endl; }
				else if( res == 3 ){ std::cout << "Check whether childe loop included in parent loop" << std::endl; }
				else if( res == 4 ){ std::cout << "Check childe loop excluded from other child loop" << std::endl; }
				else if( res == 5 ){ std::cout << "Check positive angle around vertex on the loop" << std::endl; }
				return res;
            }
		}
	}
    if( !m_BRep.AssertValid() ){ return 6; }
	return 0;
}

// Private関数
// 返り値は半辺のIDでid_vを基点としたdir1の方向に伸びる半直線からid_vを中心に時計周りに最初に出会う半辺
// この半辺の属するループはこの半直線と重なっている。
// id_vが浮遊点の場合は浮遊点周りのループが帰る
CBRep2D::CItrVertex CCadObj2D::FindCorner_HalfLine(unsigned int id_v, const CVector2D& dir1) const
{
	assert( m_VertexSet.IsObjID(id_v) );
	Com::CVector2D dir = dir1;
	dir.Normalize();
	const Com::CVector2D vec_zero(0,0);
	Cad::CBRep2D::CItrVertex itrv = this->m_BRep.GetItrVertex(id_v);
	if( itrv.CountEdge() < 2 ){ return itrv; }
	for(;!itrv.IsEnd();itrv++){
		unsigned int id_e0;   bool is_same_dir0;
		itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
		assert( this->m_EdgeSet.IsObjID(id_e0) );	// id_vが浮遊点の場合は省かれている
		const Cad::CEdge2D& e0 = this->GetEdge(id_e0);
		assert( ((is_same_dir0) ? e0.id_v_s : e0.id_v_e) == id_v );
		const CVector2D& tan0 = e0.GetTangentEdge(is_same_dir0);
		////////////////
		unsigned int id_e1;   bool is_same_dir1;
		itrv.GetIdEdge_Ahead(id_e1,is_same_dir1);
		assert( this->m_EdgeSet.IsObjID(id_e1) );	// id_vが浮遊点の場合は省かれている
		const Cad::CEdge2D& e1 = this->GetEdge(id_e1);
		assert( ((is_same_dir1) ? e1.id_v_s : e1.id_v_e) == id_v );
		const CVector2D& tan1 = e1.GetTangentEdge(is_same_dir1);
		////////////////
		assert( id_e0 != id_e1 );	// id_vが端点の場合は省かれているはず
		const double area0 = TriArea(tan1,vec_zero,tan0);
		const double area1 = TriArea(tan1,vec_zero,dir);	
		const double area2 = TriArea(dir, vec_zero,tan0);
		if( area0 > 0.0 ){	if( area1 > 0.0 && area2 > 0.0 ){ return itrv; } }
		else{				if( area1 > 0.0 || area2 > 0.0 ){ return itrv; } }
	}
	return itrv;
}

unsigned int CCadObj2D::ConnectVertex(CEdge2D edge)
{
	const unsigned int id_v1 = edge.id_v_s;
	const unsigned int id_v2 = edge.id_v_e;
	if( !m_VertexSet.IsObjID(id_v1) ){ return 0; }	
	if( !m_VertexSet.IsObjID(id_v2) ){ return 0; }
	if( id_v1 == id_v2 ){ return 0; }
	assert( m_VertexSet.IsObjID(id_v1) );
	assert( m_VertexSet.IsObjID(id_v2) );

	if( edge.itype == 0 )	// v1,v2を辺に持つ直線の辺が存在するかどうかを調べる
	{
		const std::vector<unsigned int>& id_ary_e = m_EdgeSet.GetAry_ObjID();
		for(unsigned int iid_e=0;iid_e<id_ary_e.size();iid_e++){
			const unsigned int id_e = id_ary_e[iid_e]; 	assert( m_EdgeSet.IsObjID(id_e) );
			const CEdge2D& e = this->GetEdge(id_e);
			if( e.itype != 0 ) continue;
			const unsigned int id_v_e = e.id_v_e;
			const unsigned int id_v_s = e.id_v_s;
			if( (id_v_s-id_v1)*(id_v_s-id_v2) != 0 ) continue;
			if( (id_v_e-id_v1)*(id_v_e-id_v2) != 0 ) continue;
			return 0;
		}
	}
	edge.po_s = m_VertexSet.GetObj(id_v1).point;
	edge.po_e = m_VertexSet.GetObj(id_v2).point;
	if( edge.IsCrossEdgeSelf() ){ return 0; }	// 自己交差があるなら勿論辺の生成は不可
	////////////////
	unsigned int id_l, id_e_add;
	bool is_left_ladd = false;
	{
		const CBRep2D::CItrVertex& itrv1 = this->FindCorner_HalfLine(id_v1, edge.GetTangentEdge(true)  );
		const CBRep2D::CItrVertex& itrv2 = this->FindCorner_HalfLine(id_v2, edge.GetTangentEdge(false) );
		if( itrv1.GetIdLoop() != itrv2.GetIdLoop() ) return 0;
		id_l = itrv1.GetIdLoop();
		// もしも，V1→V2と繋いだときにV1-V2の左側のループを返す
		if( itrv1.IsSameUseLoop(itrv2) && ( !itrv1.IsParent() || id_l == 0 ) ){
			const std::vector< std::pair<unsigned int,bool> >& aIdEDir = m_BRep.GetItrLoop_ConnectVertex(itrv1,itrv2);
			assert( !aIdEDir.empty() );
			double area = TriArea(edge.po_s, edge.po_e, CVector2D(0,0)) + edge.AreaEdge();
			for(unsigned int ie=0;ie<aIdEDir.size();ie++){
				unsigned int id_e = aIdEDir[ie].first;
				assert( this->IsElemID(Cad::EDGE,id_e) );
				const CEdge2D& e = this->GetEdge(id_e);
				const double earea = e.AreaEdge() + TriArea(e.po_s, e.po_e, CVector2D(0,0));
				if( aIdEDir[ie].second ){ area += earea; }
				else{                     area -= earea; }
			}
			is_left_ladd = ( area > 0 );
		}
		id_e_add = m_BRep.ConnectVertex(itrv1,itrv2, is_left_ladd);
	}
	{
		int tmp_id = this->m_EdgeSet.AddObj( std::make_pair(id_e_add,edge) );
        assert( tmp_id == (int)id_e_add );
	}

	unsigned int id_l_add;
	{
		unsigned int id_l_l, id_l_r;
		m_BRep.GetIdLoop_Edge(id_e_add, id_l_l,id_l_r);
		if( id_l_l == id_l_r ){ // 二つに分離しなかった場合
			assert( this->AssertValid()==0 );
			return id_e_add; 
		}	
		id_l_add = (is_left_ladd) ? id_l_l : id_l_r; 
		assert( id_l_add != id_l );
		assert( ((!is_left_ladd) ? id_l_l : id_l_r) == id_l );
	}

	if( m_BRep.IsElemID(Cad::LOOP,id_l) ) {
		for(;;){
			bool iflg = true;
			CBRep2D::CItrLoop itrl_add_inner = m_BRep.GetItrLoop_SideEdge(id_e_add,  is_left_ladd);
			CBRep2D::CItrLoop itrl_add_outer = m_BRep.GetItrLoop_SideEdge(id_e_add, !is_left_ladd);
			for(CBRep2D::CItrLoop itrl_c=m_BRep.GetItrLoop(id_l);!itrl_c.IsEndChild();itrl_c.ShiftChildLoop()){
				if( itrl_c.IsParent() ) continue;
				if( itrl_c.IsSameUseLoop(itrl_add_outer) ) continue;
				unsigned int ires = this->CheckInOut_ItrLoopPoint_ItrLoop(itrl_c,itrl_add_inner);
//				std::cout << ires << std::endl;
				assert( ires == 0 || ires == 2 );
				if( ires == 0 ){ 
					m_BRep.SwapItrLoop(itrl_c,id_l_add); 
//					std::cout << "Swap " << std::endl;
					iflg = false;
					break;
				}
			}
			if( iflg ) break;
		}
	}
    else{	// 外の領域に切った場合に，新しいループが自己干渉している場合はこのループは無し
		std::vector<Cad::CEdge2D> aEdge;
        if( this->CheckLoopIntersection(aEdge,id_l_add) ){
			m_BRep.SetHoleLoop(id_l_add);
			assert( this->AssertValid()==0 );
			return id_e_add;
		}
	}

	if( m_LoopSet.IsObjID(id_l) ){
		CLoop2D loop_add = this->m_LoopSet.GetObj(id_l);
		this->m_LoopSet.AddObj( std::make_pair(id_l_add,loop_add) );
	}
	else{
		this->m_LoopSet.AddObj( std::make_pair(id_l_add,CLoop2D()) );
	}
	assert( this->AssertValid()==0 );
	return id_e_add;
}


bool CCadObj2D::RemoveElement(Cad::CAD_ELEM_TYPE itype, unsigned int id)
{
	if( !this->IsElemID(itype,id) ) return false;
	if(      itype == Cad::EDGE   ){
		CBRep2D::CItrLoop itrl_l = m_BRep.GetItrLoop_SideEdge(id,true );
		CBRep2D::CItrLoop itrl_r = m_BRep.GetItrLoop_SideEdge(id,false);
		unsigned int id_l_l = itrl_l.GetIdLoop();
		unsigned int id_l_r = itrl_r.GetIdLoop();
		unsigned int id_v1, id_v2;
		m_BRep.GetIdVertex_Edge(id, id_v1,id_v2);
		CBRep2D::CItrVertex itrv1 = m_BRep.GetItrVertex(id_v1);
		CBRep2D::CItrVertex itrv2 = m_BRep.GetItrVertex(id_v2);
		bool is_del_cp = false;
		if( itrl_l.IsSameUseLoop(itrl_r) && itrl_l.IsParent() && itrl_l.GetIdLoop()!=0  
			&& itrv1.CountEdge() > 1 && itrv2.CountEdge() > 1 )	// 親ループと子ループを繋ぐ辺を削除した場合
		{
			const std::vector< std::pair<unsigned int,bool> >& aIdEDir = m_BRep.GetItrLoop_RemoveEdge(id);
			assert( !aIdEDir.empty() );
			{	// 明らかに面積が無いかどうか調べる
				unsigned int ie=0;
				for(;ie<aIdEDir.size();ie++){
					unsigned int je = 0;
					for(;je<aIdEDir.size();je++){
						if( ie == je ) continue;
						if( aIdEDir[ie].first == aIdEDir[je].first ){ 
							assert( aIdEDir[ie].second != aIdEDir[je].second );
							break;
						}
					}
					if( je == aIdEDir.size() ){ break; }
				}
				is_del_cp = ( ie == aIdEDir.size() );
			}
			if( !is_del_cp ){	// 面積を実際に調べてみる
				double area = 0.0;
				for(unsigned int ie=0;ie<aIdEDir.size();ie++){
					unsigned int id_e = aIdEDir[ie].first;
					assert( this->IsElemID(Cad::EDGE,id_e) );
					const CEdge2D& e = this->GetEdge(id_e);
					const double earea = e.AreaEdge() + TriArea(e.po_s, e.po_e, CVector2D(0,0));
					if( aIdEDir[ie].second ){ area += earea; }
					else{                     area -= earea; }
				}
				if( area < 0 ){ is_del_cp = true; }
			}
		}
		if( !m_BRep.RemoveEdge(id,is_del_cp) ){ return false; }
		m_EdgeSet.DeleteObj(id);
		if( !m_BRep.IsElemID(Cad::LOOP,id_l_l) ){ m_LoopSet.DeleteObj(id_l_l); }
		if( !m_BRep.IsElemID(Cad::LOOP,id_l_r) ){ m_LoopSet.DeleteObj(id_l_r); }
		assert( this->AssertValid()==0 );
		return true;
	}
	else if( itype == Cad::VERTEX ){ 
		CBRep2D::CItrVertex itrv = m_BRep.GetItrVertex(id);
		if( itrv.CountEdge() == 2 ){
			unsigned int id_e1,id_e2;   bool btmp1, btmp2;
			itrv.GetIdEdge_Ahead( id_e1,btmp1); 
			itrv.GetIdEdge_Behind(id_e2,btmp2);
			{			
				const unsigned int id_v1 = m_BRep.GetIdVertex_Edge(id_e1,!btmp1);
				const unsigned int id_v2 = m_BRep.GetIdVertex_Edge(id_e2,!btmp2);		
				assert( m_BRep.GetIdVertex_Edge(id_e1,btmp1) == id );
				assert( m_BRep.GetIdVertex_Edge(id_e2,btmp2) == id );
				std::vector<CEdge2D> aEdge;
				const std::vector<unsigned int>& aIdE = m_BRep.GetAryElemID(Cad::EDGE);
				for(unsigned int iie=0;iie<aIdE.size();iie++){
					const unsigned int id_e = aIdE[iie];
					if( id_e == id_e2 ) continue;
					const CEdge2D& e = this->GetEdge(id_e);
					if( id_e == id_e1 ){
						const CVector2D& v2 = this->GetVertexCoord(id_v2);
						if( btmp1 ){ e.id_v_s = id_v2; e.po_s = v2; }
						else{        e.id_v_e = id_v2; e.po_e = v2; }
					}
					aEdge.push_back(e);
				}				
				if( CheckEdgeIntersection(aEdge) == 1 ){
					return false;
				}
			}
			if( !m_BRep.RemoveVertex(id) ){ return false; }
			assert( m_BRep.IsElemID(Cad::EDGE,id_e1) );
			assert( !m_BRep.IsElemID(Cad::EDGE,id_e2) );
			m_EdgeSet.DeleteObj(id_e2);
			m_VertexSet.DeleteObj(id);
			assert( this->AssertValid()==0 );
			return true;
		}
		else if( itrv.CountEdge() == 0 ){
			if( !m_BRep.RemoveVertex(id) ){ return false; }
			m_VertexSet.DeleteObj(id);
			assert( this->AssertValid()==0 );
			return true;
		}
	}
	return false;
}

unsigned int CCadObj2D::AddPolygon( const std::vector<CVector2D>& aPoint_input, unsigned int id_l )
{
	const unsigned int npoint = aPoint_input.size();
	if( npoint < 3 ) return 0;

    std::vector<Com::CVector2D> aPoint = aPoint_input;
	{	// 交線判断
		const unsigned int n = aPoint_input.size();
		std::vector<CEdge2D> aEdge;
		for(unsigned int i=0;i<n-1;i++){
			CEdge2D e(i,i+1, 0,true,0.0);
			e.po_s = aPoint_input[i  ];
			e.po_e = aPoint_input[i+1];
			aEdge.push_back( e );
		}
		{
			CEdge2D e(n-1,0, 0,true,0.0);
			e.po_s = aPoint_input[n-1];
			e.po_e = aPoint_input[0];
			aEdge.push_back( e );
		}
		if( Cad::CheckEdgeIntersection(aEdge) != 0 ) return 0;
	}
	// (注)時計回りになってる場合でも以下は動きます
	std::vector<unsigned int> aIdV;
	std::vector<unsigned int> aIdE;
	// 点の追加
	aIdV.reserve(npoint);
	for(unsigned int ipoint=0;ipoint<npoint;ipoint++){
		unsigned int id_v0 = this->AddVertex(Cad::LOOP, id_l, aPoint[ipoint] );
		if( id_v0 == 0 ) goto FAIL_ADD_POLYGON_INSIDE_LOOP;
		aIdV.push_back(id_v0);
	}
	// 辺の追加
	aIdE.reserve(npoint);
	for(unsigned int iedge=0;iedge<npoint-1;iedge++){
		unsigned int id_e0 = this->ConnectVertex_Line( aIdV[iedge], aIdV[iedge+1] );
		if( id_e0 == 0 ) goto FAIL_ADD_POLYGON_INSIDE_LOOP;
		aIdE.push_back(id_e0);
	}
	{
		unsigned int id_e0 = this->ConnectVertex_Line( aIdV[npoint-1], aIdV[0] );
		if( id_e0 == 0 ) goto FAIL_ADD_POLYGON_INSIDE_LOOP;
		aIdE.push_back(id_e0);
	}
	assert( this->AssertValid() == 0 );
	// 新しく出来たループのIDを取得
	{	// 辺の両側のループを調べる
		unsigned int id_e0 = aIdE[ npoint-1 ];
		unsigned int id_l0, id_l1;
		m_BRep.GetIdLoop_Edge(id_e0, id_l0,id_l1);
		return ( id_l0 == id_l ) ? id_l1 : id_l0;
	}
	////////////////////////////////
	// 失敗処理ラベル
FAIL_ADD_POLYGON_INSIDE_LOOP :
	for(unsigned int iie=0;iie<aIdE.size();iie++){
		unsigned int id_e0 = aIdE[iie];
		this->RemoveElement(Cad::EDGE,id_e0);
	}
	for(unsigned int iiv=0;iiv<aIdV.size();iiv++){
		unsigned int id_v0 = aIdV[iiv];
		this->RemoveElement(Cad::VERTEX,id_v0);
	}
	assert( this->AssertValid()==0 );
	return 0;	
}

unsigned int CCadObj2D::AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id, const  Com::CVector2D& vec)
{
	if(      itype == Cad::NOT_SET || id == 0 )
	{
		unsigned int id_v_add = m_BRep.AddVertex_Loop(0);
		const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex2D(vec)) );
		assert( tmp_id ==(int)id_v_add );
		return id_v_add;
	}
	else if( itype == Cad::LOOP )
	{
		unsigned int id_l = id;
		assert( m_LoopSet.IsObjID(id_l) );
		if( !m_LoopSet.IsObjID(id_l) ) return 0;
		// この点がループの中に入っているかどうかを調べる
		if( !this->CheckIsPointInsideLoop(id_l,vec) ){ return 0; }
		unsigned int id_v_add = m_BRep.AddVertex_Loop(id_l);
		const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex2D(vec)) );
		assert( tmp_id == (int)id_v_add );
		assert( this->AssertValid()==0 );
		return id_v_add;
	}
	else if( itype == Cad::EDGE )
	{
		unsigned int id_e = id;
		assert( m_EdgeSet.IsObjID(id_e) );
		if( !m_EdgeSet.IsObjID(id_e) ) return 0;
		CEdge2D edge_old = this->GetEdge(id_e);
		// 入力点を辺の上に射影できるのか調べる
		Com::CVector2D vec_add = edge_old.GetNearestPoint(vec);
		if( SquareLength(vec_add-edge_old.po_e) < 1.0e-20 || SquareLength(vec_add-edge_old.po_s) < 1.0e-20 ){
			return 0;
		}
		////////////////////////////////
		// Leave Input Check Section

		unsigned int id_v_add = m_BRep.AddVertex_Edge(id_e);
		unsigned int id_e_add;
		{
			CBRep2D::CItrVertex itrv = m_BRep.GetItrVertex(id_v_add);
			bool is_same_dir0;
			unsigned int id_e_b,id_e_a;
			itrv.GetIdEdge_Behind(id_e_b,is_same_dir0);
			itrv.GetIdEdge_Ahead( id_e_a,is_same_dir0);
			id_e_add = ( id_e_b == id_e ) ? id_e_a : id_e_b;
		}		
		{	// Vertex追加
			const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex2D(vec_add)) );
			assert( tmp_id == (int)id_v_add );
		}
		{	// Edge追加
			const int tmp_id = m_EdgeSet.AddObj( std::make_pair(id_e_add,
				CEdge2D(id_v_add,edge_old.id_v_e, 0,false,0.0)) );
			assert( tmp_id == (int)id_e_add );
		}
		{
			CEdge2D& edge_a = this->GetEdge(id_e_add);
			edge_old.Split(edge_a, vec_add);
			CEdge2D& edge = this->GetEdge(id_e);
			assert( edge.id_v_e == id_v_add );
			edge = edge_old;
		}
		assert( this->AssertValid() == 0 );
		return id_v_add;
	}
	return 0;
}

bool CCadObj2D::SetCurve_Line(const unsigned int id_e)
{
//    std::cout << "CCadObj2D::SetCurve_Line" << id_e << std::endl;
	if( !m_EdgeSet.IsObjID(id_e) ){
		assert(0);
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
    CEdge2D& e = m_EdgeSet.GetObj(id_e);
    CEdge2D e_old = e;
	////////////////
	e.itype = 0;
	////////////////
	const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
		unsigned int id_l = aID_Loop[iid_l];
        if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
	}
	return true;
FAILURE:
	e = e_old;
	return false;
}

bool CCadObj2D::SetCurve_Arc(const unsigned int id_e, bool is_left_side, double rdist)
{
	if( !m_EdgeSet.IsObjID(id_e) ){
		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		assert(0);
		return false;
    }
	assert( m_EdgeSet.IsObjID(id_e) );
    CEdge2D& e = this->GetEdge(id_e);
	CEdge2D e_old = e;
    ////////////////////////////////
    // ここからを現在のCurveTypeによって決める,次の設定は直線の場合
	e.itype = 1;
	e.is_left_side = is_left_side;
	e.dist = rdist;
    // ここまで
    ////////////////////////////////
	const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
		unsigned int id_l = aID_Loop[iid_l];
        if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
	}
	return true;
FAILURE:
	e = e_old;
	return false;
}

bool CCadObj2D::SetCurve_Polyline(const unsigned int id_e)
{
    if( !m_EdgeSet.IsObjID(id_e) ){
        std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
        assert(0);
        return false;
    }
    assert( m_EdgeSet.IsObjID(id_e) );
    CEdge2D& e = this->GetEdge(id_e);
    CEdge2D e_old = e;
    ////////////////////////////////
    // ここからを現在のCurveTypeによって決める,次の設定は直線の場合
    const CVector2D& pos = e_old.po_s;
    const CVector2D& poe = e_old.po_e;
    std::vector<Com::CVector2D> aCo;
    e_old.GetCurve_Mesh(aCo,20);
    const double sqlen = Com::SquareLength(poe-pos);
    const Com::CVector2D& eh = (poe-pos)*(1/sqlen);
    const Com::CVector2D ev(-eh.y,eh.x);
    e.itype = 2;
    e.aRelCoMesh.clear();
    for(unsigned int ico=0;ico<aCo.size();ico++){
        double x1 = Com::Dot(aCo[ico]-pos,eh);
        double y1 = Com::Dot(aCo[ico]-pos,ev);
        e.aRelCoMesh.push_back(x1);
        e.aRelCoMesh.push_back(y1);
    }
    // ここまで
    ////////////////////////////////
    const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
    for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
        unsigned int id_l = aID_Loop[iid_l];
        if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
    }
    return true;
FAILURE:
    e = e_old;
    return false;
}

bool CCadObj2D::SetCurve_Polyline(const unsigned int id_e, const std::vector<Com::CVector2D>& aCo)
{
	if( !m_EdgeSet.IsObjID(id_e) ){
		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		assert(0);
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
	CEdge2D& e = this->GetEdge(id_e);
	CEdge2D e_old = e;
	////////////////
	e.itype = 2;
	{	// 相対座標を作る
		const unsigned int n = aCo.size();
		e.aRelCoMesh.resize(0);
		e.aRelCoMesh.reserve(n*2);
		const Com::CVector2D& pos = e.po_s;
		const Com::CVector2D& poe = e.po_e;
		const double sqlen = Com::SquareLength(poe-pos);
		const Com::CVector2D& eh = (poe-pos)*(1/sqlen);
		const Com::CVector2D ev(-eh.y,eh.x);
		for(unsigned int i=0;i<n;i++){
			double x0 = Com::Dot(aCo[i]-pos,eh);
			double y0 = Com::Dot(aCo[i]-pos,ev);
			e.aRelCoMesh.push_back(x0);
			e.aRelCoMesh.push_back(y0);
		}
	}
	////////////////
	const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
		unsigned int id_l = aID_Loop[iid_l];
        if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
	}
	return true;
FAILURE:
	e = e_old;
	return false;
}

// id_lが存在しない(0)場合は全ての辺について干渉チェックを行う
bool CCadObj2D::CheckLoopIntersection(const std::vector<Cad::CEdge2D>& aEdge, unsigned int id_l) const
{
	std::vector<unsigned int> aIdEdge;
    if( m_BRep.IsElemID(Cad::LOOP,id_l) ){
		for(CBRep2D::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			for(itrl.Begin();!itrl.IsEnd();itrl++){
				unsigned int id_e;   bool is_same_dir;
				if( !itrl.GetIdEdge(id_e,is_same_dir) ) continue;	// 浮遊点は飛ばす
				if( itrl.IsEdge_BothSideSameLoop() && !is_same_dir ){ continue; }	// 両側が一緒なら一回しかリストに加えない
				aIdEdge.push_back(id_e);
			}
		}
	}
	else{ aIdEdge = this->GetAryElemID(Cad::EDGE); }

	for(unsigned int ie=0;ie<aEdge.size();ie++){
		const unsigned int id_v_s = aEdge[ie].id_v_s;
		const unsigned int id_v_e = aEdge[ie].id_v_e;
		aEdge[ie].po_s = this->GetVertexCoord(id_v_s);
		aEdge[ie].po_e = this->GetVertexCoord(id_v_e);
	}
	////////////////
	const unsigned int nide = aIdEdge.size();
	const unsigned int ne = nide + aEdge.size();
	for(unsigned int ie=0;ie<ne;ie++){
        const CEdge2D& e_i = ( ie < nide ) ? this->GetEdge( aIdEdge[ie] ) : aEdge[ie-nide];
		if( e_i.IsCrossEdgeSelf() ){ return true; }
		const unsigned int ipo0 = e_i.id_v_s; 
		const unsigned int ipo1 = e_i.id_v_e;
		// edge_iのバウンディングボックスを取得
		double x_min_i, x_max_i, y_min_i, y_max_i;
		e_i.GetBoundingBox(x_min_i,x_max_i, y_min_i,y_max_i);
		////////////////
		for(unsigned int je=ie+1;je<ne;je++){
            const CEdge2D& e_j = ( je < nide ) ? this->GetEdge( aIdEdge[je] ) : aEdge[je-nide];
			const unsigned int jpo0 = e_j.id_v_s;
			const unsigned int jpo1 = e_j.id_v_e;
			// 共有点が無い場合
			if( (ipo0-jpo0)*(ipo0-jpo1)*(ipo1-jpo0)*(ipo1-jpo1) != 0 ){
				// BoundingBoxを用いて，交錯しないパターンを除外
				// edge_jのバウンディングボックスを取得
				double x_min_j, x_max_j, y_min_j, y_max_j;
				e_j.GetBoundingBox(x_min_j,x_max_j, y_min_j,y_max_j);
				if( x_min_j > x_max_i || x_max_j < x_min_i ) continue;	// 交錯がありえないパターンを除外
				if( y_min_j > y_max_i || y_max_j < y_min_i ) continue;	// 上に同じ
				// 交点が無いか判断する
				if( !e_i.IsCrossEdge(e_j) ){ continue; }
                return true;
			}
			else if( ipo0 == jpo0 && ipo1 == jpo1 ){
				if( e_i.IsCrossEdge_ShareBothPoints(e_j,true) == 1 ){ return true; }
			}
			else if( ipo0 == jpo1 && ipo1 == jpo0 ){
				if( e_i.IsCrossEdge_ShareBothPoints(e_j,false) == 1 ){ return true; }
			}
			else if( ipo0 == jpo0 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j, true, true)==1 ){ return true; }
			}
			else if( ipo0 == jpo1 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j, true,false)==1 ){ return true; }
			}
			else if( ipo1 == jpo0 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j,false, true)==1 ){ return true; }
			}
			else if( ipo1 == jpo1 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j,false,false)==1 ){ return true; }
			}
		}
	}
	return false;
}

// DXFファイルへの書き出し
bool CCadObj2D::WriteToFile_dxf(const std::string& file_name) const
{
	FILE *fp;
	if( (fp = ::fopen(file_name.c_str(),"w"))== NULL ){
		fclose(fp);
		assert(0);
		return false;
	}

	// オブジェクトのBoundingBoxを得る
	double x_min,x_max,  y_min,y_max;
	{
		const std::vector<unsigned int>& aIdEdge = this->m_EdgeSet.GetAry_ObjID();
		assert( aIdEdge.size() > 0 );
		{
			const unsigned int id_e = aIdEdge[0];
			const CEdge2D& edge = this->GetEdge(id_e);
			edge.GetBoundingBox(x_min,x_max,y_min,y_max);
		}
		for(unsigned int iid_e=1;iid_e<aIdEdge.size();iid_e++){
			const unsigned int id_e = aIdEdge[iid_e];
			const CEdge2D& edge = this->GetEdge(id_e);
			double x_min0,x_max0,  y_min0,y_max0;
			edge.GetBoundingBox(x_min0,x_max0,y_min0,y_max0);
			x_min = ( x_min < x_min0 ) ? x_min : x_min0;
			x_max = ( x_max > x_max0 ) ? x_max : x_max0;
			y_min = ( y_min < y_min0 ) ? y_min : y_min0;
			y_max = ( y_max > y_max0 ) ? y_max : y_max0;
		}
	}

	// ヘッダーセクション
	fprintf(fp, "  0\nSECTION\n");
	fprintf(fp, "  2\nHEADER\n");
	fprintf(fp, "  9\n$ACADVER\n  1\nAC1009\n");
	fprintf(fp, "  9\n$EXTMIN\n  10\n%lf\n  20\n%lf\n",x_min,y_min);
	fprintf(fp, "  9\n$EXTMAX\n  10\n%lf\n  20\n%lf\n",x_max,y_max);
	fprintf(fp, "  0\nENDSEC\n");

	// テーブルセクション
	fprintf(fp, "  0\nSECTION\n");
	fprintf(fp, "  2\nTABLES\n");
	fprintf(fp, "  0\nENDSEC\n");

	// ブロックセクション
	fprintf(fp, "  0\nSECTION\n");
	fprintf(fp, "  2\nBLOCKS\n");
	fprintf(fp, "  0\nENDSEC\n");

	// エンティティセクション
	fprintf(fp,"  0\nSECTION\n");
	fprintf(fp,"  2\nENTITIES\n");
	const std::vector<unsigned int>& aIdLoop = this->m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++){
		const unsigned int id_l = aIdLoop[iid_l];
		assert( m_LoopSet.IsObjID(id_l) );
        std::auto_ptr<Cad::ICad2D::CItrLoop> pItr = this->GetItrLoop(id_l);
		for(;;){
            for(;!pItr->IsEnd();(*pItr)++){
				bool is_same_dir;
				unsigned int id_e;
				if( !pItr->GetIdEdge(id_e,is_same_dir) ){ assert(0); fclose(fp); return false; }
				unsigned int id_vs, id_ve;
				this->GetIdVertex_Edge(id_vs, id_ve,id_e);
				assert( this->IsElemID(Cad::VERTEX,id_vs) );
				assert( this->IsElemID(Cad::VERTEX,id_ve) );
				const CVector2D& ps = this->GetVertexCoord(id_vs);
				const CVector2D& pe = this->GetVertexCoord(id_ve);
				if( this->GetEdgeCurveType(id_e) == 0 ){
					fprintf(fp,"  0\nLINE\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n",id_l);
					fprintf(fp,"  10\n%lf\n",ps.x);
					fprintf(fp,"  20\n%lf\n",ps.y);
					fprintf(fp,"  11\n%lf\n",pe.x);
					fprintf(fp,"  21\n%lf\n",pe.y);
				}
				else if( this->GetEdgeCurveType(id_e) == 1 ){
					const CEdge2D& edge = this->m_EdgeSet.GetObj(id_e);
					edge.po_s = ps;
					edge.po_e = pe;
					CVector2D pc;
					double r;
					double d1, d2;
					{
						edge.GetCenterRadius(pc,r);
						CVector2D vs = ps - pc;
						CVector2D ve = pe - pc;
						double ds = atan2(vs.y,vs.x); ds = ds * 180.0 / 3.14159265; if( ds < 0.0 ) ds += 360;
						double de = atan2(ve.y,ve.x); de = de * 180.0 / 3.14159265; if( de < 0.0 ) de += 360;
						if( edge.is_left_side ){ d1 = de; d2 = ds; }
						else{                    d1 = ds; d2 = de; }
					}
					fprintf(fp,"  0\nARC\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n  100\nAcDbCircle\n",id_l);
					fprintf(fp,"  10\n%lf\n",pc.x);	// 中心のｘ座標指定
					fprintf(fp,"  20\n%lf\n",pc.y);	// 中心のｙ座標指定
					fprintf(fp,"  40\n%lf\n",r);	// 半径指定
					fprintf(fp,"  100\nAcDbArc\n");
					fprintf(fp,"  50\n%lf\n",d1);
					fprintf(fp,"  51\n%lf\n",d2);
				}
			}
			if( !pItr->ShiftChildLoop() ) break;
		}
/*		
		fprintf(fp,"0\nPOLYLINE\n10\n0\n20\n0\n66\n1\n70\n0\n");
		for(;!itr.IsEnd();itr++){
			const unsigned int id_v = itr.GetIdVertex();
			assert( this->IsID_Vertex(id_v) );
			const CVector2D& v = this->GetVertexCoord(id_v);
			fprintf(fp,"0\nVERTEX\n10\n%lf\n20\n%lf\n",v.x,v.y);
		}
		fprintf(fp,"0\nENDSEQ\n");
*/
	}
	fprintf(fp, "  0\nENDSEC\n  0\nEOF\n");
	fclose(fp);
	return true;
}

bool CCadObj2D::Serialize( Com::CSerializer& arch )
{
	if( arch.IsLoading() ){	// 読み込み時の処理
		this->Clear();
		char stmp1[256];
		arch.Get("%s",stmp1);
		assert( strncmp(stmp1,"$$$$$$$$",8)==0 );
		arch.Get("%s",stmp1);
		if( strncmp(stmp1,"CadObj2D",8)!=0 ) return true;
		int nv, ne, nl;
		{
			arch.Get("%d%d%d",&nv, &ne, &nl);
			assert( nv>0 ); assert( ne>0 ); assert( nl>0 );
			m_VertexSet.Reserve(nv*2);
			m_EdgeSet.Reserve(ne*2);
			m_LoopSet.Reserve(nl*2);
		}
		////////////////////////////////////////////////
        for(int iv=0;iv<nv;iv++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$",4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"V2D",3) ==0 );
			int id;		arch.Get("%d",&id);		assert( id>0 );
			double x,y;	arch.Get("%lf%lf",&x,&y);
			int tmp_id = m_VertexSet.AddObj( std::make_pair(id,CVertex2D(CVector2D(x,y) )) );
			assert(tmp_id==id);
		}
        for(int ie=0;ie<ne;ie++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$",4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"E2D",3) ==0 );
			int id;				arch.Get("%d",&id);					assert( id>0 );
			int id_v_s,id_v_e;	arch.Get("%d%d",&id_v_s,&id_v_e);	assert( id_v_s>0 && id_v_e>0 );
            int itype;			arch.Get("%d",&itype);				assert( itype == 0 || itype == 1 || itype == 2);
			int i_is_left_side;
			double dist;
			arch.Get("%d%lf",&i_is_left_side,&dist);
			assert( i_is_left_side == 0 || i_is_left_side == 1 );
			const bool is_left_side = (i_is_left_side != 0 );
            std::vector<double> aRelCo;
            int npo;
            arch.Get("%d",&npo);
            assert( npo >= 0 );
            for(unsigned int ipo=0;ipo<(unsigned int)npo;ipo++){
                double x,y;
                arch.Get("%lf%lf",&x,&y);
                aRelCo.push_back(x);
                aRelCo.push_back(y);
            }
            Cad::CEdge2D e(id_v_s,id_v_e,itype,is_left_side,dist);
            e.aRelCoMesh = aRelCo;
            const int tmp_id = m_EdgeSet.AddObj( std::make_pair(id,e) );
			assert( tmp_id == id );            
		}
        for(int il=0;il<nl;il++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$",4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"L2D",3) ==0 );
            int id;         arch.Get("%d",&id);		assert( id>0 );
			int ilayer;		arch.Get("%d",&ilayer);
            double c[3];    arch.Get("%lf%lf%lf",&c[0],&c[1],&c[2]);
            CLoop2D l;
			l.ilayer = ilayer;
            l.m_color[0]=c[0];  l.m_color[1]=c[1];  l.m_color[2]=c[2];
            const int tmp_id = m_LoopSet.AddObj( std::make_pair(id,l) );
			assert( tmp_id == id );
		}
		m_BRep.Serialize(arch);
		this->AssertValid();
		return true;
	}
	else{ // 書き込み時の処理
        // クラスの名前の指定，サイズの指定
        {
			arch.Out("$$$$$$$$\n");
			arch.Out("CadObj2D\n");
			arch.Out("%d %d %d\n",m_VertexSet.GetAry_ObjID().size(), m_EdgeSet.GetAry_ObjID().size(),m_LoopSet.GetAry_ObjID().size());
		}
        // Vertex2Dの出力
        {
			const std::vector<unsigned int>& id_ary = m_VertexSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_v = id_ary[iid];
				assert( m_VertexSet.IsObjID(id_v) );
				const CVertex2D& v = m_VertexSet.GetObj(id_v);
				arch.Out("$$$$\n");
				arch.Out("V2D\n");
				arch.Out("%d\n",id_v);
				arch.Out("%lf %lf\n",v.point.x,v.point.y);
			}
		}
        // Edge2Dの出力
        {
			const std::vector<unsigned int> id_ary = m_EdgeSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_e = id_ary[iid];
				assert( m_EdgeSet.IsObjID(id_e) );
				const CEdge2D& e = m_EdgeSet.GetObj(id_e);
				arch.Out("$$$$\n");
				arch.Out("E2D\n");
				arch.Out("%d\n",id_e);
				arch.Out("%d %d\n",e.id_v_s,e.id_v_e);
				arch.Out("%d\n",e.itype);
				arch.Out("%d %lf\n",e.is_left_side,e.dist);
                const unsigned int n = e.aRelCoMesh.size()/2;
                arch.Out("%d\n",n);
                for(unsigned int i=0;i<n;i++){
                    arch.Out("%lf %lf\n",e.aRelCoMesh[i*2+0], e.aRelCoMesh[i*2+1]);
                }
			}
		}
        // Loop2Dの出力
        {
            const std::vector<unsigned int> id_ary= m_LoopSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_l = id_ary[iid];
				assert( m_LoopSet.IsObjID(id_l) );
				const CLoop2D& l = m_LoopSet.GetObj(id_l);
				arch.Out("$$$$\n");
				arch.Out("L2D\n");
				arch.Out("%d\n",id_l);
				arch.Out("%d\n",l.ilayer);
                arch.Out("%lf %lf %lf\n",l.m_color[0],l.m_color[1],l.m_color[2]);
			}
		}
		m_BRep.Serialize(arch);
	}
	return true;
}
