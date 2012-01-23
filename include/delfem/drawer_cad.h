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
@brief CADを可視化するためのクラス(Cad::View::CDrawer_Cad2D)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_CAD_H)
#define DRAWER_CAD_H

#include <vector>

#include "delfem/cad_obj2d.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace Msh{
	class CBarAry;
	class CTriAry2D;
  class SVertex;
}

namespace Cad{
namespace View{

/*! 
@brief CadをOpenGLで可視化するクラス
@ingroup CAD
*/
class CDrawer_Cad2D : public Com::View::CDrawer
{
public : 
	/*!
   @brief コンストラクタ
   @param[in] cad ＣＡＤの形状
   @param[in] imode 表示モード{0:面を描画する,1:面を描画しない(セレクションでは面を描画)}(デフォルトで0)
   @remaks 面を描画しないモードに設定してもセレクションはされるので，ピックは可能
	*/
  CDrawer_Cad2D(const Cad::CCadObj2D& cad)
  {    
    tex_scale = 1;
    tex_cent_x=0, tex_cent_y=0;    
		this->m_is_anti_aliasing = false;
    m_linewidth = 3;
    m_pointsize = 5;
		// create check pattern with 2pt interval
		for(unsigned int j=0;j<4;j++){
			unsigned int i;
			for(i=0; i<8 ;i++) m_mask[j*32   +i] = 0x33;
			for(i=0; i<8 ;i++) m_mask[j*32+ 8+i] = 0xcc;
			for(i=0; i<8 ;i++) m_mask[j*32+16+i] = 0x33;
			for(i=0; i<8 ;i++) m_mask[j*32+24+i] = 0xcc;
		}	
    this->UpdateCAD_TopologyGeometry(cad);
	}  
  CDrawer_Cad2D(){
    tex_scale = 1;  
    tex_cent_x=0, tex_cent_y=0;        
		this->m_is_anti_aliasing = false;
    m_linewidth = 3;
    m_pointsize = 5;   
 		// create check pattern with 2pt interval
    for(unsigned int j=0;j<4;j++){
      unsigned int i;
      for(i=0; i<8 ;i++) m_mask[j*32   +i] = 0x33;
      for(i=0; i<8 ;i++) m_mask[j*32+ 8+i] = 0xcc;
      for(i=0; i<8 ;i++) m_mask[j*32+16+i] = 0x33;
      for(i=0; i<8 ;i++) m_mask[j*32+24+i] = 0xcc;
    }
  }
	virtual ~CDrawer_Cad2D(){
    Clear();
	}
	// virtual関数
  void Clear(){
		for(unsigned int i=0;i<m_apDrawPart.size();i++){ delete m_apDrawPart[i]; }    
    m_apDrawPart.resize(0);
  }

  //! 幾何形状を更新する(トポロジーの変化は存在しないとして)
  void UpdateCAD_Geometry(const Cad::CCadObj2D&);
  //! トポロジーと幾何を更新する
  bool UpdateCAD_TopologyGeometry(const Cad::CCadObj2D&);

	//! @{
	//! 描画
	virtual void Draw() const;
	//! ピッキングのための描画をする
	virtual void DrawSelection(unsigned int idraw) const;
  //! 線の太さを設定
  void SetLineWidth(unsigned int linewidth){ m_linewidth = linewidth; }
  //! 点の大きさを設定
  void SetPointSize(unsigned int pointsize){ m_pointsize = pointsize; }
	//! get bounding box. z direction is -1 to 1.
  // if rot==0 or no argment return object axis arrigned bounding box
  // if rot!=0 rot should be 3x3 maxtirx and return rotated bounding box
	virtual Com::CBoundingBox3D GetBoundingBox( const double* rot=0 ) const {
		return m_vertex_ary.GetBoundingBox(	rot );
	}
  // delete someday
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
  void SetIsShow(bool is_show, Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
  void SetIsShow(bool is_show, Cad::CAD_ELEM_TYPE part_type, const std::vector<unsigned int>& aIdPart );  
  // end delete someday

	void SetModeShow(Cad::CAD_ELEM_TYPE itype, unsigned int id, int imode_show);
  int GetModeShow(Cad::CAD_ELEM_TYPE itype, unsigned int id);
	void ClearModeShow(int imode_show);
    //! @}

	void GetCadPartID(const int selec_flag[], Cad::CAD_ELEM_TYPE& part_type, unsigned int& part_id, int& ictrl);
  void HideEffected(const Cad::CCadObj2D& cad_2d, Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
  void ShowEffected(const Cad::CCadObj2D& cad_2d, Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
  void SetRigidDisp(unsigned int id_l, double xdisp, double ydisp);
  
  // texture
  void SetTextureScale(double tex_scale);
  void SetTexCenter(double cent_x, double cent_y){
    this->tex_cent_x = cent_x;
    this->tex_cent_y = cent_y;    
  }  
  void GetTexCenter(double& cent_x, double& cent_y){
    cent_x = this->tex_cent_x;
    cent_y = this->tex_cent_y;
  }       
  void EnableUVMap(bool is_uv_map);
private:
	////////////////
	class CDrawPart{
	public : 
		CDrawPart(){
			nelem = 0; npoel = 0; pIndexArray = 0; 
			itype = Cad::VERTEX; imode_show = 0;
      color[0]=0; color[1]=0; color[2]=0;
			height = 0;
			xdisp=0; ydisp=0;
		}
		virtual ~CDrawPart(){
			this->Clear();
		}
		void Clear(){
			delete[] pIndexArray;
			pIndexArray = 0;
			nelem = 0;
			npoel = 0;
		}
		void DrawElements() const;
		bool SetBarAry(const Msh::CBarAry& BarAry);
		bool SetTriAry(const Msh::CTriAry2D& TriAry);
    bool SetVertex(const Msh::SVertex& vtx);    
		void SetHeight(double h){ height = h; }
	public:
    int imode_show;  // -1:not_show 0:default 1:hilight 2:selected
    float color[3];
    ////////////////
		Cad::CAD_ELEM_TYPE itype;
		unsigned int id_cad;
    unsigned int id_msh;  // the mesh id in the tesselation
    ////////////////
		unsigned int nelem;
		unsigned int npoel;
    unsigned int* pIndexArray;
		////////////////
		double height;
		double xdisp, ydisp;
    ////////////////
    CURVE_TYPE itype_curve;		//!< 0:Line, 1:Arc, 2:Mesh 3:Bezier
    std::vector<Com::CVector2D> aCtrlPoint;
	};
  /*
	////////////////
	class CDrawPart_CadVertex{
	public:
    int imode_show; // -1:not_show 0:default 1:hilight 2:selected
    unsigned int id_cad;
		unsigned int id_msh;
		unsigned int id_v;
		////////////////
		double height;
	};	
   */
private:
  unsigned char m_mask[128];	//!< mask of hilight reagion
  std::vector<CDrawPart*> m_apDrawPart;				//! 辺や面など
//  std::vector<CDrawPart_CadVertex> m_aIndexVertex;	//! The Vtx
  Com::View::CVertexArray m_vertex_ary;	//! vertex array
  unsigned int m_linewidth;   //!< the width of the line
  unsigned int m_pointsize;   //!< the size of vertices
  double tex_cent_x, tex_cent_y;
  double tex_scale;
};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*!
@brief ラバーバンドクラス
@remark ここはCVector3Dじゃなくて，CVector2Dを使うべき．
*/
class CDrawerRubberBand : public Com::View::CDrawer
{
public:
	CDrawerRubberBand(const Com::CVector3D& initial_position){
		this->initial = initial_position;
		imode = 0;
	}
	CDrawerRubberBand(const Cad::CCadObj2D& cad, unsigned int id_v_cad);
	CDrawerRubberBand(const Cad::CCadObj2D& cad, unsigned int id_v, const Com::CVector3D& initial);
	virtual void Draw() const;
	virtual void DrawSelection(unsigned int idraw) const{}
	virtual Com::CBoundingBox3D GetBoundingBox( double rot[] ) const;
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual unsigned int WhatKindOfYou() const { return 0; }

	void SetMousePosition(const Com::CVector3D& mouse){
		this->mouse = mouse;
	}
private:
	// 0 : connect inital and mouse;
	// 1 : connect fix[] and mouse
	// 2 : connect fix_stat[].second and fix[ fix_stat[].first ]
	int imode;
	Com::CVector3D mouse;
	Com::CVector3D initial;
	std::vector<Com::CVector3D> fix;
	std::vector< std::pair<unsigned int,Com::CVector3D> > fix_stat;
};


} // end namespace View
} // end namespace Cad

#endif
