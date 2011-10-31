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
#pragma warning ( disable : 4786 ) 
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif  /* _WIN32 */
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

//#include <GL/glut.h>



#include <iostream>
#include "delfem/rigid/linearsystem_rigid.h"
#include "delfem/rigid/rigidbody.h"
#include "delfem/indexed_array.h"

static void CalcInvMat(double* a, const unsigned int& n, int& info )
{
	double tmp1;
  
	info = 0;
	unsigned int i,j,k;
	for(i=0;i<n;i++){
		if( fabs(a[i+i*n]) < 1.0e-30 ){
			info = 1;
			return;
		}
		if( a[i+i*n] < 0.0 ){
			info--;
		}
		tmp1 = 1.0 / a[i+i*n];
		a[i+i*n] = 1.0;
		for(k=0;k<n;k++){
			a[i+k*n] *= tmp1;
		}
		for(j=0;j<n;j++){
			if( j!=i ){
				tmp1 = a[j+i*n];
				a[j+i*n] = 0.0;
				for(k=0;k<n;k++){
					a[j+k*n] -= tmp1*a[i+k*n];
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void Ls::CLinearSystem_RigidBody_CRS2::SetRigidSystem
(const std::vector<Rigid::CRigidBody3D>& aRB,
 const std::vector<Rigid::CConstraint*>& aConst)
{
  this->Clear();
  nRB = aRB.size();
  nConst = aConst.size();
  const unsigned int nblk = nRB+nConst;
  {
    m_aBlkSize.resize(nblk);
    for(unsigned int irb=0;irb<nRB;irb++){
      m_aBlkSize[irb] = aRB[irb].GetDOF();
    }
    for(unsigned int icst=0;icst<nConst;icst++){
      m_aBlkSize[nRB+icst] = aConst[icst]->GetDOF();
    }
    m_mat.Initialize(nblk,m_aBlkSize, nblk,m_aBlkSize);
    m_residual.Initialize(nblk,m_aBlkSize);
    m_update.Initialize(nblk,m_aBlkSize);
  }
  {
    Com::CIndexedArray crs;
    crs.InitializeSize(nblk);
    crs.index.resize(nblk+1);
    crs.index[0] = 0;
    for(unsigned int i=0;i<nRB;i++){ crs.index[i+1]=0; }
    for(unsigned int i=nRB;i<nRB+nConst;i++){ crs.index[i+1]=0; }
    for(unsigned int icst=0;icst<nConst;icst++){
      const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
      for(unsigned int i=0;i<aIndRB.size();i++){
        const unsigned int irb0 = aIndRB[i];
        crs.index[irb0+1] += 1;
      }
      crs.index[icst+nRB+1] += aIndRB.size();
    }
    for(unsigned int i=0;i<nblk;i++){ 
      crs.index[i+1] = crs.index[i+1] + crs.index[i];
    }
    const unsigned int ncrs = crs.index[nblk];
    crs.array.resize(ncrs);
    for(unsigned int icst=0;icst<nConst;icst++){
      const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
      for(unsigned int i=0;i<aIndRB.size();i++){
        const unsigned int irb0 = aIndRB[i];
        const unsigned int icrs0 = crs.index[icst+nRB];
        crs.array[icrs0] = irb0;
        crs.index[icst+nRB] += 1;
        const unsigned int icrs1 = crs.index[irb0];
        crs.array[icrs1] = icst+nRB;
        crs.index[irb0] += 1;
      }
    }
    for(int i=nblk;i>0;i--){ 
      crs.index[i] = crs.index[i-1];
    }
    crs.index[0] = 0;
    crs.Sort();
    assert( crs.index[nRB+nConst] == ncrs );
    m_mat.AddPattern(crs);      
  }  
}


double Ls::CLinearSystem_RigidBody_CRS2::DOT(int iv1,int iv2)
{
  MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
  MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
  return vec1*vec2;
}

bool Ls::CLinearSystem_RigidBody_CRS2::COPY(int iv1,int iv2)  // {v2} := {v1}
{
  if( iv1 == iv2 ) return true;
  
  MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
  MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
  vec2 = vec1;
  return true;
}

bool Ls::CLinearSystem_RigidBody_CRS2::SCAL(double d,int iv)  // {v1} := alpha * {v1}
{
  MatVec::CVector_Blk& vec = this->GetVector(iv);
  vec *= d;
  return true;
}

bool Ls::CLinearSystem_RigidBody_CRS2::AXPY(double d,int iv1,int iv2) // {v2} := alpha*{v1} + {v2}
{ 
  MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
  MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
  vec2.AXPY(d,vec1);
  return true;
}

bool Ls::CLinearSystem_RigidBody_CRS2::MATVEC(double a,int iv1,double b,int iv2)  //  {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
{
  MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
  MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
  m_mat.MatVec(a,vec1,b,vec2);
  return true; 
}



bool Ls::CLinearSystem_RigidBody_CRS2::UpdateValueOfRigidSystem_NewmarkBetaAPrime
(std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
 double dt, double newmark_gamma, double newmark_beta, 
 bool is_first) const
{
  double tmp[6];
  for(unsigned int irb=0;irb<aRB.size();irb++){
    const unsigned int nlen = m_update.Len(irb);
    assert( nlen <= 6 );
    for(unsigned int ilen=0;ilen<nlen;ilen++){
      tmp[ilen] = m_update.GetValue(irb,ilen);
    }
    aRB[irb].UpdateSolution_NewmarkBetaAPrime( tmp, 
                                              dt, newmark_gamma, newmark_beta, 
                                              is_first);
  }
  for(unsigned int icst=0;icst<aConst.size();icst++){
    const unsigned int nlen = m_update.Len(icst+nRB);
    assert( nlen <= 6 );
    for(unsigned int ilen=0;ilen<nlen;ilen++){
      tmp[ilen] = m_update.GetValue(icst+nRB,ilen);
    }
    aConst[icst]->UpdateLambda_NewmarkBetaAPrime( tmp, dt, newmark_gamma, newmark_beta );
  }
  return true;
}

bool Ls::CLinearSystem_RigidBody_CRS2::UpdateValueOfRigidSystem_BackwardEular
(std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
 double dt, bool is_first) const
{
  double tmp[6];
  for(unsigned int irb=0;irb<aRB.size();irb++){
    const unsigned int nlen = m_update.Len(irb);
    assert( nlen <= 6 );
    for(unsigned int ilen=0;ilen<nlen;ilen++){
      tmp[ilen] = m_update.GetValue(irb,ilen);
    }
    aRB[irb].UpdateSolution_BackwardEular( tmp, dt, is_first );
  }
  for(unsigned int icst=0;icst<aConst.size();icst++){
    const unsigned int nlen = m_update.Len(icst+nRB);
    assert( nlen <= 6 );
    for(unsigned int ilen=0;ilen<nlen;ilen++){
      tmp[ilen] = m_update.GetValue(icst+nRB,ilen);
    }
    aConst[icst]->UpdateLambda_BackwardEular( tmp, dt );  
  }
  return true;
}

