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
@brief リアルタイム有限要素法解析クラス
@author Nobuyuki Umetani
*/

#if !defined(GUI_LISTENER_INTERACTIVE_ANALYSIS_2D_H)
#define GUI_LISTENER_INTERACTIVE_ANALYSIS_2D_H

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

#include <GL/gl.h>
#include <GL/glu.h>

#include <set>
#include <map>
#include <vector>

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"

#include "delfem/cad_obj2d_move.h"
#include "delfem/drawer_cad.h"
#include "delfem/mesher2d_edit.h"
#include "delfem/drawer_msh.h"

#include "delfem/field.h"

#include "delfem/analysis2d_interface.h"

//! リアルタイム有限要素法解析クラス
class CGuiListner_Analysis2D_Interactive
{
public:
	class CCadAgent : public Cad::ICad2D
	{
	public:
        CCadAgent(CGuiListner_Analysis2D_Interactive* pAnalysis);
        virtual ~CCadAgent(){}

        virtual int GetLayer(Cad::CAD_ELEM_TYPE type,unsigned int id) const;
		virtual void GetLayerMinMax(int& ilayer_min, int& ilayer_max) const;

        //! 辺の両側のループのIDを返す
        virtual bool GetIdLoop_Edge(unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const;

		//! idが使われているかどうかを調べる関数
		virtual bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const;
		//! すべてのIDを配列にして返す関数
		virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE) const;
		//! ID:id_lのループの色を返す(本来このクラスは位相と幾何情報以外を持つべきではないかもしれないので暫定的)
		virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const;
		//! ID:id_lのループの面積を返えす
		virtual double GetArea_Loop(unsigned int id_l) const;
		//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
        virtual std::auto_ptr<ICad2D::CItrLoop>   GetItrLoop(  unsigned int id_l) const;
        //! ID:id_vのループを構成する頂点や辺をめぐるイテレータを返す関数
        virtual std::auto_ptr<ICad2D::CItrVertex> GetItrVertex(unsigned int id_v) const;

		//! 辺の始点と終点の頂点のIDを返す
		virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const;
		//! 辺の形状タイプを返す (０ 線分, １ 円弧);
		virtual int GetEdgeCurveType(const unsigned int& id_e) const;
		//! ID:id_eの辺の情報を得る
		virtual bool GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const;
        virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const;
        virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen = -1) const;
        //! ID:id_eの辺をndiv個に分割したものを得る. 但し始点，終点はそれぞれps,peとする
        virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
                                   unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const;
		//! 頂点の座標を得る 
		virtual const Com::CVector2D& GetVertexCoord(unsigned int id_v) const;

		////////////////
        //! ID:id_lのループの色を設定する(本来このクラスは位相と幾何情報以外を持つべきではないかもしれないので暫定的)
        bool SetColor_Loop(unsigned int id_l, double color[3] );
        bool Move(Cad::CAD_ELEM_TYPE itype_cad, unsigned int id_part_cad,
				const Com::CVector2D& pos_obj_pre, const Com::CVector2D& pos_obj_cur, double tor=-1.0){
            return m_pAnalysis->Cad_Move(itype_cad,id_part_cad, pos_obj_pre,pos_obj_cur,tor);
		}

	private:
		CGuiListner_Analysis2D_Interactive* m_pAnalysis;
	};
    friend class CCadAgent;
public:
	//! コンストラクタ
	CGuiListner_Analysis2D_Interactive(){
		is_updated_coord = false;
		is_updated_edge = false;
		m_is_solve_cad_change = true;
		pDrawerCAD = 0;
		pDrawerMsh = 0;
	}
	//! デストラクタ
    virtual ~CGuiListner_Analysis2D_Interactive(){
		if( pDrawerCAD != 0 ) delete pDrawerCAD;
		if( pDrawerMsh != 0 ) delete pDrawerMsh;
	}
	////////////////
	// publicな純粋仮想関数

	void SetAnalysisInitialize(IAnalysis2D* pAnalysis);
	void Draw(unsigned int idaw_obj = 0);
	void DrawSelection();
	void Solve_fromCad();	//!< 形状データから解きなおす関数
	void Solve_ifNeeded();	//!< 必要に応じて解く関数
	void Serialize( Com::CSerializer& arch );	//!< ファイルへの(保存/読み出し)関数

	//! 形状が変わったときに，方程式を解くかどうか
	void Enable_SolveCadChange(bool flg){ m_is_solve_cad_change = flg; }

	////////////////////////////////
	// GUIに関するクラス

    std::auto_ptr<Cad::ICad2D> Cad()
	{
        std::auto_ptr<Cad::ICad2D> p( new CCadAgent(this) );
		return p;
	}

	// バウンディング・ボックス取得
	// rot : ３×３の回転行列
	Com::CBoundingBox GetBoundingBox(double rot[]) const;

	// CADに追従していないメッシュを追従させる(Idle時に実行する)
	void FollowMshToCad_ifNeeded();

	////////////////////////////////
	// ピック関係の関数

	void SetSelection(const std::vector<Com::View::SSelectedObject>& aSelecObj);
	void Cad_GetPicked(Cad::CAD_ELEM_TYPE& itype, unsigned int& id, double& x, double& y) const ;

	////////////////////////////////
	// CADに関する関数

    void HilightCadTypeID(Cad::CAD_ELEM_TYPE itype, unsigned int id){
        if( !this->cad_2d.IsElemID(itype,id) ) return;
		pDrawerCAD->ClearSelected();
        pDrawerCAD->AddSelected(itype,id);
	}

	////////////////////////////////
	// Mshに関する関数

	// 要素数を設定
	void Msh_SetElemSize(unsigned int esize){
		if( esize == 0 ) return;
		mesh_2d.m_esize = esize;
		this->Solve_fromCad();
	}

	// このCADのループが解析メッシュが切られているかかどうかを調べる
    bool Msh_IsMesh_CadLoop(unsigned int id_l) const {
		return mesh_2d.setIdLCad_CutMesh.find(id_l) != mesh_2d.setIdLCad_CutMesh.end();
    }

	// ループを穴に設定する
    void Msh_DelMesh_CadLoop(unsigned int id_l){
		if( !cad_2d.IsElemID(Cad::LOOP,id_l) ) return;
		std::set<unsigned int>::iterator itr = mesh_2d.setIdLCad_CutMesh.find(id_l);
		if( itr == mesh_2d.setIdLCad_CutMesh.end() ) return;
        mesh_2d.setIdLCad_CutMesh.erase( itr );
		this->Solve_fromCad();
	}
	
	// 穴になっているループを埋める
	void Msh_SetMesh_CadLoop(unsigned int id_l){
		if( !cad_2d.IsElemID(Cad::LOOP,id_l) ) return;
        if( Msh_IsMesh_CadLoop(id_l) ) return;
		mesh_2d.setIdLCad_CutMesh.insert( id_l );
		this->Solve_fromCad();
	}
	void Msh_SetMeshingMode_Length(double elen){
		mesh_2d.m_imode_meshing = 2;
		mesh_2d.m_elen = elen;
	}
	void Msh_SetMeshingMode_ElemSize(unsigned int esize){
		mesh_2d.m_imode_meshing = 1;
		mesh_2d.m_esize = esize;
	}

	//! CADモデルをDXFに保存する
	void File_WriteDXF(const std::string& fname) const {
		cad_2d.WriteToFile_dxf(fname);
    }
public:
	unsigned int Cad_MakeRectRegion(unsigned int id_l, const Com::CVector2D& pos_obj0, const Com::CVector2D& pos_obj1 );
    unsigned int Cad_ConnectVertex_Line(unsigned int id_v0, unsigned int id_v1);
    unsigned int Cad_ConnectVertex(const Cad::CEdge2D& e);
    bool Cad_Remove(Cad::CAD_ELEM_TYPE itype, unsigned int id);
    unsigned int Cad_AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id, const Com::CVector2D& pos_obj0);
    void Cad_SetColor_Loop(unsigned int id_l, double color[3]);	
	bool Cad_SetCurveType(unsigned int id_e, unsigned int itype);	//! 既存のカーブに出来るだけ沿うように新しいカーブを設定する
    bool Cad_SetCurve_Polyline(unsigned int id_e, const std::vector<Com::CVector2D>& aVec);
    bool Cad_DragArc(unsigned int id_part_cad,const Com::CVector2D& pos_obj0, double tol = -1);
    bool Cad_Move(Cad::CAD_ELEM_TYPE itype_cad, unsigned int id_part_cad,
		const Com::CVector2D& pos_obj_pre, const Com::CVector2D& pos_obj_cur, double tor = -1.0);    
    //! id_eがメッシュ辺なら，スムージングをかけて滑らかにする(radiusが負なら全体をスムージングする)
    bool Cad_SmoothingPolylineEdge(unsigned int id_e, unsigned int niter, const Com::CVector2D& pos, double radius);
protected:
	void CreatedLoop(unsigned int id_l_created, unsigned int id_l_from){}
	void InitDrawer();
protected:
	// Cad関係の変数
    Cad::CCadObj2D_Move cad_2d;
	std::set<unsigned int> setIdVCad_NeedFollow;
	Cad::View::CDrawer_Cad2D* pDrawerCAD;

	// メッシュ関係の変数
	Msh::CMesher2D_Edit mesh_2d;
	Msh::View::CDrawerMsh2D* pDrawerMsh;
	
    // ピックされた物
	double m_picked_x, m_picked_y;
	unsigned int m_id_cad_part;
	Cad::CAD_ELEM_TYPE m_itype_cad_part;

	bool m_is_solve_cad_change;
	bool is_updated_cad;
	bool is_updated_coord;	// 要素の座標がUpdateされたかどうか
	bool is_updated_edge;	// 要素の辺がUpdateされたかどうか

    IAnalysis2D* pAnalysis;
};

#endif
