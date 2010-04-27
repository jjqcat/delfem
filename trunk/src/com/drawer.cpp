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

#include "delfem/drawer.h"
#include "delfem/camera.h"
#include "delfem/vector3d.h"
#include "delfem/quaternion.h"

using namespace Com::View;

Com::CBoundingBox CDrawerArray::GetBoundingBox( double rm[] ) const
{
	if( m_drawer_ary.empty() ){
		return CBoundingBox(-0.5,0.5, -0.5,0.5, -0.5,0.5);
	}	
    CBoundingBox bb = m_drawer_ary[0]->GetBoundingBox(rm);
	for(unsigned int idraw=1;idraw<m_drawer_ary.size();idraw++){
        bb += m_drawer_ary[idraw]->GetBoundingBox(rm);
	}
	return bb;
}

Com::CBoundingBox CVertexArray::GetBoundingBox( double rot[] ) const
{
	if( pVertexArray == 0 ){ return Com::CBoundingBox(); }
	if( ndim == 2 ){		
		Com::CBoundingBox bb;
		{
			const double x1 = pVertexArray[0];
			const double y1 = pVertexArray[1];
			const double z1 = 0.0;
			const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
			const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
			const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
			bb = Com::CBoundingBox(x2,x2, y2,y2, z2,z2);
		}
		for(unsigned int ipoin=1;ipoin<npoin;ipoin++){
			const double x1 = pVertexArray[ipoin*2  ];
			const double y1 = pVertexArray[ipoin*2+1];
			const double z1 = 0.0;
			const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
			const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
			const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
			bb.x_max = ( x2 > bb.x_max ) ? x2 : bb.x_max;  bb.x_min = ( x2 < bb.x_min ) ? x2 : bb.x_min;			
			bb.y_max = ( y2 > bb.y_max ) ? y2 : bb.y_max;  bb.y_min = ( y2 < bb.y_min ) ? y2 : bb.y_min;			
			bb.z_max = ( z2 > bb.z_max ) ? z2 : bb.z_max;  bb.z_min = ( z2 < bb.z_min ) ? z2 : bb.z_min;			
		}
		return bb;
	}
	if( ndim == 3 ){
		Com::CBoundingBox bb;
		{
			const double x1 = pVertexArray[0];
			const double y1 = pVertexArray[1];
			const double z1 = 0.0;
			const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
			const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
			const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
			bb = Com::CBoundingBox(x2,x2, y2,y2, z2,z2);
		}
		for(unsigned int ipoin=1;ipoin<npoin;ipoin++){
			const double x1 = pVertexArray[ipoin*3  ];
			const double y1 = pVertexArray[ipoin*3+1];
			const double z1 = pVertexArray[ipoin*3+2];
			const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
			const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
			const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
			bb.x_max = ( x2 > bb.x_max ) ? x2 : bb.x_max;  bb.x_min = ( x2 < bb.x_min ) ? x2 : bb.x_min;			
			bb.y_max = ( y2 > bb.y_max ) ? y2 : bb.y_max;  bb.y_min = ( y2 < bb.y_min ) ? y2 : bb.y_min;			
			bb.z_max = ( z2 > bb.z_max ) ? z2 : bb.z_max;  bb.z_min = ( z2 < bb.z_min ) ? z2 : bb.z_min;			
		}
		return bb;
	}
	return Com::CBoundingBox();
}

void CDrawerArray::InitTrans(Com::View::CCamera& mvp_trans){
	{	// ñ?g]???[?h????f?
		unsigned int irot_mode = 0;
		for(unsigned int idraw=0;idraw<m_drawer_ary.size();idraw++){
			unsigned int irot_mode0 = m_drawer_ary[idraw]->GetSutableRotMode();
			irot_mode = (irot_mode0>irot_mode) ? irot_mode0 : irot_mode;
		}
		if(      irot_mode == 1 ){ mvp_trans.SetRotationMode(ROT_2D);  }
		else if( irot_mode == 2 ){ mvp_trans.SetRotationMode(ROT_2DH); }
		else if( irot_mode == 3 ){ mvp_trans.SetRotationMode(ROT_3D);  }
	}
	{	// ?o?E?g?f?B?g?O?{?b?N?X????f?
		double rot[9];
		mvp_trans.RotMatrix33(rot);
		Com::CBoundingBox bb = this->GetBoundingBox( rot );
		double x_cent = ( bb.x_max + bb.x_min )*0.5;
		double y_cent = ( bb.y_max + bb.y_min )*0.5;
		double z_cent = ( bb.z_max + bb.z_min )*0.5;
		double x_size = bb.x_max - bb.x_min;
		double y_size = bb.y_max - bb.y_min;
		double z_size = bb.z_max - bb.z_min;
		mvp_trans.SetObjectCenter(x_cent,y_cent,z_cent);
		mvp_trans.SetObjectSize(x_size,y_size,z_size);
	}
	mvp_trans.Fit();
}
