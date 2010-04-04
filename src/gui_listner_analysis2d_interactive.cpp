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

#if defined(__VISUALC__)
#pragma warning ( disable : 4996 )
#pragma warning ( disable : 4786 )
#endif

#include "delfem/gui_listner_analysis2d_interactive.h"

CGuiListner_Analysis2D_Interactive::CCadAgent::CCadAgent(CGuiListner_Analysis2D_Interactive* pAnalysis){
    this->m_pAnalysis = pAnalysis;
}

//! idが使われているかどうかを調べる関数
bool CGuiListner_Analysis2D_Interactive::CCadAgent::IsElemID(Cad::CAD_ELEM_TYPE itype,unsigned int id) const
{
	return m_pAnalysis->cad_2d.IsElemID(itype,id);
}

//! すべてのIDを配列にして返す関数
const std::vector<unsigned int> CGuiListner_Analysis2D_Interactive::CCadAgent::GetAryElemID(Cad::CAD_ELEM_TYPE itype) const
{
	return m_pAnalysis->cad_2d.GetAryElemID(itype);
}

//! ID:id_lのループの色を返す(本来このクラスは位相と幾何情報以外を持つべきではないかもしれないので暫定的)
bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetColor_Loop(unsigned int id_l, double color[3] ) const
{
	return m_pAnalysis->cad_2d.GetColor_Loop(id_l,color);
}

//! ID:id_lのループの面積を返えす
double CGuiListner_Analysis2D_Interactive::CCadAgent::GetArea_Loop(unsigned int id_l) const
{
	return m_pAnalysis->cad_2d.GetArea_Loop(id_l);
}

//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
std::auto_ptr<Cad::ICad2D::CItrLoop> CGuiListner_Analysis2D_Interactive::CCadAgent::GetItrLoop(
		unsigned int id_l) const
{
	return m_pAnalysis->cad_2d.GetItrLoop(id_l);
}

//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
std::auto_ptr<Cad::ICad2D::CItrVertex> CGuiListner_Analysis2D_Interactive::CCadAgent::GetItrVertex(
        unsigned int id_v) const
{
    return m_pAnalysis->cad_2d.GetItrVertex(id_v);
}

//! 辺の始点と終点の頂点のIDを返す
bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetIdVertex_Edge(
		unsigned int &id_vs, unsigned int& id_ve, unsigned int id_e) const
{
	return m_pAnalysis->cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
}

//! 辺の形状タイプを返す (０ 線分, １ 円弧);
int CGuiListner_Analysis2D_Interactive::CCadAgent::GetEdgeCurveType(const unsigned int& id_e) const
{
	return m_pAnalysis->cad_2d.GetEdgeCurveType(id_e);
}

//! ID:id_eの辺の情報を得る
bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetCurve_Arc(
		unsigned int id_e, bool& is_left_side, double& dist) const
{
	return m_pAnalysis->cad_2d.GetCurve_Arc(id_e,is_left_side,dist);
}

bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetCurve_Polyline(
		unsigned int id_e, std::vector<double>& aRelCoMesh) const
{
    return m_pAnalysis->cad_2d.GetCurve_Polyline(id_e,aRelCoMesh);
}

bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetCurve_Polyline(
        unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen) const
{
    return m_pAnalysis->cad_2d.GetCurve_Polyline(id_e,aCo,elen);
}


//! ID:id_eの辺をndiv個に分割したものを得る. 但し始点，終点はそれぞれps,peとする
bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
        unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const
{
    return m_pAnalysis->cad_2d.GetCurve_Polyline(id_e,aCo,ndiv,ps,pe);
}

//! 頂点の座標を得る 
const Com::CVector2D& CGuiListner_Analysis2D_Interactive::CCadAgent::GetVertexCoord(
		unsigned int id_v) const
{
	return m_pAnalysis->cad_2d.GetVertexCoord(id_v);
}

//! 辺の両側のループのIDを返す
bool CGuiListner_Analysis2D_Interactive::CCadAgent::GetIdLoop_Edge(
        unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const
{
    return m_pAnalysis->cad_2d.GetIdLoop_Edge(id_l_l,id_l_r,id_e);
}

int CGuiListner_Analysis2D_Interactive::CCadAgent::GetLayer(Cad::CAD_ELEM_TYPE type, unsigned int id_l ) const
{
    return m_pAnalysis->cad_2d.GetLayer(type,id_l);
}

void CGuiListner_Analysis2D_Interactive::CCadAgent::GetLayerMinMax(int& ilayer_min, int& ilayer_max) const
{
	m_pAnalysis->cad_2d.GetLayerMinMax(ilayer_min,ilayer_max);
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void CGuiListner_Analysis2D_Interactive::DrawSelection()
{
	pDrawerCAD->DrawSelection(0);
}

Com::CBoundingBox CGuiListner_Analysis2D_Interactive::GetBoundingBox(double rot[]) const
{
	return pDrawerCAD->GetBoundingBox(rot);
}

void CGuiListner_Analysis2D_Interactive::FollowMshToCad_ifNeeded()
{
    if( setIdVCad_NeedFollow.empty() ){ return; }
    std::cout << "FollowMshToCad_ifNeeded()" << std::endl;
    std::vector<unsigned int> aIdV_del;
    std::set<unsigned int>::iterator itr = setIdVCad_NeedFollow.begin();
    for(;itr!=setIdVCad_NeedFollow.end();itr++){
        const unsigned int id_v = *itr;
        assert( cad_2d.IsElemID(Cad::VERTEX,id_v) );
        unsigned int itype_ope;
        if( mesh_2d.FitMeshToCad_Vertex(cad_2d, id_v, itype_ope) ){
            aIdV_del.push_back(id_v);
        }
        if( itype_ope == 0 ){
            this->Solve_fromCad();
            is_updated_coord = false;
            is_updated_edge = false;
            setIdVCad_NeedFollow.clear();
            break;
        }
        is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
        is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
    }
    for(unsigned int iid_v=0;iid_v<aIdV_del.size();iid_v++){
        setIdVCad_NeedFollow.erase( aIdV_del[iid_v] );
    }
}

void CGuiListner_Analysis2D_Interactive::SetSelection(const std::vector<Com::View::SSelectedObject>& aSelecObj)
{
	pDrawerCAD->ClearSelected();
	m_id_cad_part = 0;  m_itype_cad_part = Cad::NOT_SET;
	if( aSelecObj.size() == 0 ) return;
	pDrawerCAD->GetCadPartID(aSelecObj[0].name,m_itype_cad_part,m_id_cad_part);
	pDrawerCAD->AddSelected( aSelecObj[0].name );
	m_picked_x = aSelecObj[0].picked_pos.x;
	m_picked_y = aSelecObj[0].picked_pos.y;
}

void CGuiListner_Analysis2D_Interactive::Cad_GetPicked(Cad::CAD_ELEM_TYPE& itype, unsigned int& id,
                                      double& x, double& y) const {
    id = m_id_cad_part;
    itype = m_itype_cad_part;
    x = m_picked_x;
    y = m_picked_y;
}

void CGuiListner_Analysis2D_Interactive::Cad_SetColor_Loop(unsigned int id_l, double color[3])
{
    cad_2d.SetColor_Loop(id_l,color);
    pDrawerCAD->UpdateCAD_Geometry(cad_2d);
}


////////////////////////////////////////////////////////////////
// 幾何変化(位相は同じ)

bool CGuiListner_Analysis2D_Interactive::Cad_SetCurve_Polyline(
        unsigned int id_e, const std::vector<Com::CVector2D>& aVec)
{
    if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
    if( cad_2d.GetEdgeCurveType(id_e) != 2 ) return true;
    {
        unsigned int id_msh = mesh_2d.GetElemID_FromCadID(id_e,Cad::EDGE);
        Msh::MSH_TYPE msh_type;
        unsigned int nelem,iloc,id_cad;
        mesh_2d.GetMshInfo(id_msh,nelem,msh_type,iloc,id_cad);
    }
    if( !cad_2d.SetCurve_Polyline(id_e,aVec) ) return false;
    is_updated_cad = true;
    if( m_is_solve_cad_change ){
        unsigned int itype_ope;
        if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
            unsigned int id_vs, id_ve;
            cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
            setIdVCad_NeedFollow.insert( id_vs );
            setIdVCad_NeedFollow.insert( id_ve );
        }
        is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
        is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
        is_updated_cad = true;
    }
    return true;
}

bool CGuiListner_Analysis2D_Interactive::Cad_DragArc(unsigned int id_e, const Com::CVector2D& pos_obj, double tol)
{
	if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
	if( cad_2d.GetEdgeCurveType(id_e) != 1 ) return true;
    if( !cad_2d.DragArc(id_e,pos_obj,tol) ) return false;
    is_updated_cad = true;
    if( m_is_solve_cad_change ){
        unsigned int itype_ope;
        if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
            unsigned int id_vs, id_ve;
            cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
            setIdVCad_NeedFollow.insert( id_vs );
            setIdVCad_NeedFollow.insert( id_ve );
        }
        is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
        is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
        is_updated_cad = true;
    }
    return true;
}

bool CGuiListner_Analysis2D_Interactive::Cad_SmoothingPolylineEdge(unsigned int id_e,
        unsigned int niter, const Com::CVector2D& pos, double radius)
{
    if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
    if( cad_2d.GetEdgeCurveType(id_e) != 2 ) return true;
    if( !cad_2d.SmoothingPolylineEdge(id_e, niter,pos,radius) ) return false;
    is_updated_cad = true;
    if( m_is_solve_cad_change ){
        unsigned int itype_ope;
        if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
            unsigned int id_vs, id_ve;
            cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
            setIdVCad_NeedFollow.insert( id_vs );
            setIdVCad_NeedFollow.insert( id_ve );
        }
        is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
        is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
        is_updated_cad = true;
    }
    return true;
}

bool CGuiListner_Analysis2D_Interactive::Cad_SetCurveType(unsigned int id_e, unsigned int itype )
{
//    std::cout << "CGuiListner_Analysis2D_Interactive::Cad_SetCurveType" << std::endl;
    if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
    if(      itype == 0 ){ is_updated_cad = cad_2d.SetCurve_Line(id_e); }
    else if( itype == 1 ){ is_updated_cad = cad_2d.SetCurve_Arc( id_e); }
    else if( itype == 2 ){ is_updated_cad = cad_2d.SetCurve_Polyline(id_e); }
	if( !is_updated_cad ) return false;
    ////////////////
	if(  m_is_solve_cad_change ){	// メッシュを動かす
		unsigned int itype_ope;
		if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
			unsigned int id_vs, id_ve;
			cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
			setIdVCad_NeedFollow.insert( id_vs );
			setIdVCad_NeedFollow.insert( id_ve );
		}
		is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
		is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
	}
	return true;
}

bool CGuiListner_Analysis2D_Interactive::Cad_Move(
	Cad::CAD_ELEM_TYPE itype_cad_part, unsigned int id_cad_part, 
	const Com::CVector2D& pos_pre, const Com::CVector2D& pos_cur, double tor)
{
	if( itype_cad_part == Cad::VERTEX && cad_2d.IsElemID(Cad::VERTEX,id_cad_part) ){
		if( !cad_2d.MoveVertex(id_cad_part,pos_cur,tor) ){ // CADオブジェクトの頂点を動かす
			std::cout << "Cad false" << std::endl;
            return false;
		}
		is_updated_cad = true;
        if(  !m_is_solve_cad_change  ) return true;
		if( mesh_2d.GetElemID_FromCadID(id_cad_part,Cad::VERTEX ) == 0 ) return true;
		// ここからはFEMメッシュを動かす
		unsigned int itype_ope = 0;
		if( !mesh_2d.FitMeshToCad_Vertex(cad_2d, id_cad_part, itype_ope) ){ // Meshオブジェクトの頂点を動かす
			if( itype_ope == 0 ){ 
				this->Solve_fromCad();
                return true;
			}
			else{ setIdVCad_NeedFollow.insert( id_cad_part ); }
		}
		else{ setIdVCad_NeedFollow.erase( id_cad_part ); }
		is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
		is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
	}
	if( itype_cad_part == Cad::EDGE && cad_2d.IsElemID(Cad::EDGE,id_cad_part) ){
		if( !cad_2d.MoveEdge(id_cad_part, pos_cur-pos_pre) ){ // CADオブジェクトの辺を動かす
			std::cout << "Cad false" << std::endl;
            return false;
		}
		is_updated_cad = true;
        if(  !m_is_solve_cad_change  ) return true;
		if( mesh_2d.GetElemID_FromCadID(id_cad_part,Cad::EDGE ) == 0 ) return true;
		// ここからはFEMメッシュを動かす
		unsigned int itype_ope;
		if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_cad_part,   itype_ope) ){ // Meshオブジェクトの辺を動かす
			unsigned int id_vs, id_ve;
			cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_cad_part);
			setIdVCad_NeedFollow.insert( id_vs );
			setIdVCad_NeedFollow.insert( id_ve );
		}
		is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
		is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
	}
	if( itype_cad_part == Cad::LOOP && cad_2d.IsElemID(Cad::LOOP,id_cad_part) ){
		if( !cad_2d.MoveLoop(id_cad_part, pos_cur-pos_pre, false) ){ // CADオブジェクトの面を動かす
			std::cout << "Cad false" << std::endl;
            return false;
		}
		is_updated_cad = true;
        if(  !m_is_solve_cad_change  ) return true;
//		if( mesh_2d.GetElemID_FromCadID(id_cad_part,Cad::LOOP ) == 0 ) return true;	// このループの周りにメッシュが切ってある場合
		// ここからはFEMメッシュを動かす
		unsigned int itype_ope = 0;
		if( !mesh_2d.FitMeshToCad_Loop(cad_2d,id_cad_part,   itype_ope) ){ // Meshオブジェクトの頂点を動かす
            std::auto_ptr<Cad::ICad2D::CItrLoop> pItr = cad_2d.GetItrLoop(id_cad_part);
            Cad::ICad2D::CItrLoop& itr = *pItr;
			for(;!itr.IsEnd();itr++){
				const unsigned int id_v = itr.GetIdVertex();
				setIdVCad_NeedFollow.insert( id_v );
			}
		}
		is_updated_coord = is_updated_coord || (itype_ope&1);	// 節点の移動
		is_updated_edge  = is_updated_edge  || (itype_ope&2);	// 要素の切り替え
	}
    return true;
}

////////////////////////////////////////////////////////////////
// 位相変化

unsigned int CGuiListner_Analysis2D_Interactive::Cad_MakeRectRegion(unsigned int id_l, 
		const Com::CVector2D& pos_obj0, const Com::CVector2D& pos_obj1 )
{
	const double x0 = pos_obj0.x;
	const double y0 = pos_obj0.y;
	const double x1 = pos_obj1.x;
    const double y1 = pos_obj1.y;
    unsigned int id_l0;
    {
        std::vector< Com::CVector2D > aPoint;
        aPoint.push_back( Com::CVector2D(x0,y0) );
        aPoint.push_back( Com::CVector2D(x1,y0) );
        aPoint.push_back( Com::CVector2D(x1,y1) );
        aPoint.push_back( Com::CVector2D(x0,y1) );
        id_l0 = cad_2d.AddPolygon( aPoint, id_l );
    }
	if( id_l0 == 0 ) return false;
    this->CreatedLoop(id_l0,id_l);
    mesh_2d.setIdLCad_CutMesh.insert(id_l0);
    this->Solve_fromCad();
	return id_l0;
}

unsigned int CGuiListner_Analysis2D_Interactive::Cad_ConnectVertex_Line(
        unsigned int id_v0, unsigned int id_v1)
{
	const unsigned int id_e = cad_2d.ConnectVertex_Line(id_v0,id_v1);
	if( id_e == 0 ) return 0;
	this->Solve_fromCad();
	return id_e;
}

unsigned int CGuiListner_Analysis2D_Interactive::Cad_ConnectVertex(const Cad::CEdge2D& e)
{
    const unsigned int id_e = cad_2d.ConnectVertex(e);
    if( id_e == 0 ) return 0;
    this->Solve_fromCad();
    return id_e;
}

unsigned int CGuiListner_Analysis2D_Interactive::Cad_AddVertex(
		Cad::CAD_ELEM_TYPE itype, unsigned int id, const Com::CVector2D& vec){
	unsigned int id_v = cad_2d.AddVertex(itype,id,vec);
	if( id_v != 0 ){ this->Solve_fromCad(); }
	return id_v;
}

bool CGuiListner_Analysis2D_Interactive::Cad_Remove(Cad::CAD_ELEM_TYPE itype, unsigned int id)
{
	if( itype == Cad::LOOP || itype == Cad::VERTEX ){        
		if( !cad_2d.RemoveElement(itype,id) ) return false;
	}
    if( itype == Cad::EDGE ){
		unsigned int id_l_l, id_l_r;
		cad_2d.GetIdLoop_Edge(id_l_l,id_l_r, id);
		const bool is_l_l_cut = ( mesh_2d.setIdLCad_CutMesh.find(id_l_l) != mesh_2d.setIdLCad_CutMesh.end() );
		const bool is_l_r_cut = ( mesh_2d.setIdLCad_CutMesh.find(id_l_r) != mesh_2d.setIdLCad_CutMesh.end() );
		if( !cad_2d.RemoveElement(itype,id) ) return false;
		if( is_l_l_cut != is_l_r_cut ){	// 一方がメッシュがあり，一方にメッシュがない場合は，できるだけメッシュを作る
			if( cad_2d.IsElemID(Cad::LOOP,id_l_l) ){ mesh_2d.setIdLCad_CutMesh.insert(id_l_l); }
			if( cad_2d.IsElemID(Cad::LOOP,id_l_r) ){ mesh_2d.setIdLCad_CutMesh.insert(id_l_r); }
		}
		// メッシュを切るループが存在しなければ消す
		for(std::set<unsigned int>::iterator itr=mesh_2d.setIdLCad_CutMesh.begin();itr!=mesh_2d.setIdLCad_CutMesh.end();itr++){
			unsigned int id_l = *itr;
			if( cad_2d.IsElemID(Cad::LOOP,id_l) ){ continue; }
			mesh_2d.setIdLCad_CutMesh.erase(itr);	// 辺を消すことで消えるループは１つ
			break;	// ここでbreakしないとエラーになる
		}
    }
    this->Solve_fromCad();
    return true;
}

void CGuiListner_Analysis2D_Interactive::Draw(unsigned int idraw_object)
{
	if( this->is_updated_cad ){
		this->is_updated_cad = false;
		pDrawerCAD->UpdateCAD_TopologyGeometry(cad_2d);
	}
	if(      idraw_object == 0 ){
		pAnalysis->Draw();
        pDrawerCAD->SetIsShow(true,Cad::LOOP, this->cad_2d.GetAryElemID(Cad::LOOP) );
        std::set<unsigned int>::iterator itr = mesh_2d.setIdLCad_CutMesh.begin();
        for(;itr!=mesh_2d.setIdLCad_CutMesh.end();itr++){
            const unsigned int id_l = *itr;
            pDrawerCAD->SetIsShow(false,Cad::LOOP,id_l);
        }
        pDrawerCAD->Draw();
	}
	else if( idraw_object == 1 ){
        pDrawerCAD->SetIsShow(true,Cad::LOOP, this->cad_2d.GetAryElemID(Cad::LOOP) );
		pDrawerCAD->Draw();
	}
	else if( idraw_object == 2 ){
        pDrawerMsh->Draw();
	}
	else if( idraw_object == 3 ){
        pDrawerMsh->Draw();
        pDrawerCAD->SetIsShow(true,Cad::LOOP, this->cad_2d.GetAryElemID(Cad::LOOP) );
        std::set<unsigned int>::iterator itr = mesh_2d.setIdLCad_CutMesh.begin();
        for(;itr!=mesh_2d.setIdLCad_CutMesh.end();itr++){
            const unsigned int id_l = *itr;
            pDrawerCAD->SetIsShow(false,Cad::LOOP,id_l);
        }
        pDrawerCAD->Draw();
	}
}

void CGuiListner_Analysis2D_Interactive::InitDrawer()
{
	if( pDrawerCAD != 0 ){ delete pDrawerCAD; }
    pDrawerCAD = new Cad::View::CDrawer_Cad2D(cad_2d);
	if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
	pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
	pAnalysis->InitDrawer();
}


void CGuiListner_Analysis2D_Interactive::SetAnalysisInitialize(IAnalysis2D* pAnalysis)
{
	this->pAnalysis = pAnalysis;
	pAnalysis->SetModelProblem(cad_2d,mesh_2d);
	if( pDrawerCAD != 0 ){ delete pDrawerCAD; }
    pDrawerCAD = new Cad::View::CDrawer_Cad2D(cad_2d);
	if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
	pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
}

void CGuiListner_Analysis2D_Interactive::Solve_fromCad()
{
	mesh_2d.Meshing(cad_2d);
	pAnalysis->SolveInitial(cad_2d,mesh_2d);
	this->InitDrawer();
}

void CGuiListner_Analysis2D_Interactive::Solve_ifNeeded(){	
	if( is_updated_cad && !(is_updated_coord||is_updated_edge) ){
		is_updated_cad = false;
		mesh_2d.Meshing(cad_2d);
		pAnalysis->SolveInitial(cad_2d,mesh_2d);
		return;
	}
	if( !m_is_solve_cad_change ){ return; }
	SOLVER_FLAG solver_flg = this->pAnalysis->Solve(mesh_2d,is_updated_coord,is_updated_edge);
	if( solver_flg != SUCCESS ){
		mesh_2d.Meshing(cad_2d);
		this->pAnalysis->SolveInitial(cad_2d,mesh_2d);
	}
	is_updated_coord = false;
	is_updated_edge  = false;
}

void CGuiListner_Analysis2D_Interactive::Serialize( Com::CSerializer& arch ){
	if( arch.IsLoading() ){
		const unsigned int buff_size = 256;
		char class_name[buff_size];
		arch.ReadDepthClassName(class_name,buff_size);
		assert( strcmp(class_name,"CGuiListner_Analysis2D_Interactive\n") == 0 );        
        arch.ShiftDepth(true);
		cad_2d.Serialize(arch);
		mesh_2d.Serialize(arch,true);
		arch.ShiftDepth(false);
		mesh_2d.Meshing(cad_2d);
	}
	else{
		arch.WriteDepthClassName("CGuiListner_Analysis2D_Interactive");
        arch.ShiftDepth(true);
		cad_2d.Serialize(arch);
		mesh_2d.Serialize(arch,true);
        arch.ShiftDepth(false);
	}
}