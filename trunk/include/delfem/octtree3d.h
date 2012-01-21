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
 @brief 3 dimentional vector clas (Com::CVector3D)
 @author Nobuyuki Umetani
 */

#if !defined(OCT_TREE_3D_H)
#define OCT_TREE_3D_H

namespace Com
{
  class COctTree
  {
  private:
    class CCell{
    public:
      CCell(unsigned int iparent){ 
        for(unsigned int i=0;i<8;i++){ data[i]=0; }
        size = 0;
        this->iparent = iparent;
      }
      int data[8];
      int size;
      unsigned int iparent;
    };
  public:
    COctTree();
    void SetBoundingBox( const CBoundingBox3D& bb );
    int InsertPoint( unsigned int ipo_ins, const CVector3D& vec_ins );
    
    bool Check() const;
    int GetIndexCell_IncludePoint( const CVector3D& vec ) const;
    void GetAllPointInCell( unsigned int icell0, std::vector<unsigned int>& ipoins ) const;
    void GetBoundaryOfCell(unsigned int icell0, CBoundingBox3D& bb ) const;
    bool IsPointInSphere( double radius, const CVector3D& vec ) const;
  private:
    std::vector< CCell > m_aCell;
    std::vector< std::pair<unsigned int,CVector3D> > m_aVec;
    double x_min,x_max, y_min,y_max, z_min,z_max;
  };  
}

#endif
