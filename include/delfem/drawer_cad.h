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

#include "delfem/cad2d_interface.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace Msh{
	class CBarAry;
	class CTriAry2D;
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
    CDrawer_Cad2D(const Cad::ICad2D& cad)
    {    
		this->m_is_anti_aliasing = false;
        m_linewidth = 3;
        m_pointsize = 5;
		// 2ptのチェック柄を作る
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
		this->m_is_anti_aliasing = false;
        m_linewidth = 3;
        m_pointsize = 5;
        // 2ptのチェック柄を作る
        for(unsigned int j=0;j<4;j++){
            unsigned int i;
            for(i=0; i<8 ;i++) m_mask[j*32   +i] = 0x33;
            for(i=0; i<8 ;i++) m_mask[j*32+ 8+i] = 0xcc;
            for(i=0; i<8 ;i++) m_mask[j*32+16+i] = 0x33;
            for(i=0; i<8 ;i++) m_mask[j*32+24+i] = 0xcc;
        }
    }
	virtual ~CDrawer_Cad2D(){
		for(unsigned int i=0;i<m_apIndexAry.size();i++){ delete m_apIndexAry[i]; }
	}
	// virtual関数

    //! 幾何形状を更新する(トポロジーの変化は存在しないとして)
    void UpdateCAD_Geometry(const Cad::ICad2D&);
    //! トポロジーと幾何を更新する
    bool UpdateCAD_TopologyGeometry(const Cad::ICad2D&);

	//! @{
	//! 描画
	virtual void Draw() const;
	//! ピッキングのための描画をする
	virtual void DrawSelection(unsigned int idraw) const;
    //! 線の太さを設定
    void SetLineWidth(unsigned int linewidth){ m_linewidth = linewidth; }
    //! 点の大きさを設定
    void SetPointSize(unsigned int pointsize){ m_pointsize = pointsize; }
	//! バウンディング・ボックスを得る
	virtual Com::CBoundingBox GetBoundingBox( double rot[] ) const {
		return m_vertex_ary.GetBoundingBox(	rot );
	}
	//! Hilight表示する要素に加える(selection_flagから)
	virtual void AddSelected(const int selec_flag[]);
	//! Hilight表示する要素に加える
	virtual void AddSelected(Cad::CAD_ELEM_TYPE itype, unsigned int id);
	//! Hilight表示をやめる
	virtual void ClearSelected();
    //! @}

	void GetCadPartID(const int selec_flag[], Cad::CAD_ELEM_TYPE& part_type, unsigned int& part_id);
    void HideEffected(const Cad::ICad2D& cad_2d, Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
    void ShowEffected(const Cad::ICad2D& cad_2d, Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
    void Hide(Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
    void Show(Cad::CAD_ELEM_TYPE part_type, unsigned int part_id);
private:	
	////////////////
	class CDrawPart{
	public : 
		CDrawPart(){
			nelem = 0; npoel = 0; pIndexArray = 0; 
			is_selected = false; itype = Cad::VERTEX; is_show=true;
            color[0]=0; color[1]=0; color[2]=0;
			height = 0;
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
		bool Set(const Msh::CBarAry& BarAry);
		bool Set(const Msh::CTriAry2D& TriAry);
		void SetHeight(double h){ height = h; }
	public:
		bool is_show;
		bool is_selected;
        double color[3];
        ////////////////
		Cad::CAD_ELEM_TYPE itype;
		unsigned int id_cad;
        unsigned int id_msh;
        ////////////////
		unsigned int nelem;
		unsigned int npoel;
        unsigned int* pIndexArray;
		////////////////
		double height;
	};
	////////////////
	class CDrawPart_CadVertex{
	public:
		unsigned int id_v;
		unsigned int id_cad;
		unsigned int id_msh;
        ////////////////
		bool is_selected;
		bool is_show;
		////////////////
		double height;
	};	
private:
    bool is_front_and_back;	//!< 両面書くかを指定
    unsigned char m_mask[128];	//!< 選択領域をドットで見せるのためのマスク
    std::vector<CDrawPart*> m_apIndexAry;				//! 辺や面など
    std::vector<CDrawPart_CadVertex> m_aIndexVertex;	//! 頂点
    Com::View::CVertexArray m_vertex_ary;	//! 点配列
    unsigned int m_linewidth;   //!< 線の太さ
    unsigned int m_pointsize;   //!< 点の大きさ
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
	CDrawerRubberBand(const Cad::ICad2D& cad, unsigned int id_v_cad);
	CDrawerRubberBand(const Cad::ICad2D& cad, unsigned int id_v, const Com::CVector3D& initial);
	virtual void Draw() const;
	virtual void DrawSelection(unsigned int idraw) const{}
	virtual Com::CBoundingBox GetBoundingBox( double rot[] ) const;
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
