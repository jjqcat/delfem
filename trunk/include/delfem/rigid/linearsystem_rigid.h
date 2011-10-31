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
 @brief 剛体クラス(Com::CRigidBody3D)の実装
 @author Nobuyuki Umetani
 */

#if !defined(LINEAR_SYSTEM_RIGID_H)
#define LINEAR_SYSTEM_RIGID_H

#include <vector>
#include <cassert>
#include <math.h>

#include "delfem/vector3d.h"
#include "delfem/matrix3d.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/ls/linearsystem_interface_solver.h"
#include "delfem/matvec/matdiafrac_blkcrs.h"


////////////////////////////////////////////////////////////////

namespace Rigid{
  class CRigidBody3D;
  class CConstraint;
}

namespace Ls
{  
  class CLinearSystem_RigidBody // Abstruct Class
  {
  public:
    virtual unsigned int GetSizeRigidBody() const = 0;
    
    ////////////////
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
                             const Com::CVector3D& vres, double d ) = 0;
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, unsigned int size,
                             const double* eres, double d ) = 0;
    virtual void SubResidual(const unsigned int ind, bool is_rb, 
                             const double* res) = 0;
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans) = 0;
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                                  unsigned int indl, bool is_rb_l, unsigned int offsetl,
                                  const Com::CVector3D& vec, double d, bool is_colum) = 0;
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr, unsigned int sizer,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl, unsigned int sizel,
                           const double* emat, double d) = 0;
    virtual bool UpdateValueOfRigidSystem_NewmarkBetaAPrime(
                                          std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
                                          double dt, double newmark_gamma, double newmark_beta, 
                                          bool is_first) const = 0;
  };
  
  
  class CLinearSystem_RigidBody_CRS2 : public CLinearSystem_RigidBody, public LsSol::ILinearSystem_Sol
  {
    friend class CPreconditioner_RigidBody_CRS;
  public:
    CLinearSystem_RigidBody_CRS2(){}
    CLinearSystem_RigidBody_CRS2(const std::vector<Rigid::CRigidBody3D>& aRB,
                                 const std::vector<Rigid::CConstraint*>& aConst){
      this->SetRigidSystem(aRB,aConst);
    }
    virtual ~CLinearSystem_RigidBody_CRS2(){
      this->Clear();
    }
    void Clear(){   // データを全て削除
    }
    void SetRigidSystem(const std::vector<Rigid::CRigidBody3D>& aRB,
                        const std::vector<Rigid::CConstraint*>& aConst);
    
    unsigned int GetSizeRigidBody() const { return nRB; }
    
    ////////////////////////////////////////////////////////////////
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
                             const Com::CVector3D& vres, double d )
    {
      unsigned int iblk = (is_rb) ? ind : ind+nRB;
      assert( offset+3 <= m_residual.Len(iblk) );
      m_residual.AddValue(iblk,offset+0,vres.x*d);
      m_residual.AddValue(iblk,offset+1,vres.y*d);
      m_residual.AddValue(iblk,offset+2,vres.z*d);
    }
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, unsigned int size,
                             const double* eres, double d )
    {
      unsigned int iblk = (is_rb) ? ind : ind+nRB;
      assert( offset+size <= m_residual.Len(iblk) );
      for(unsigned int i=0;i<size;i++){
        m_residual.AddValue(iblk,offset+i,eres[i]*d);
      }
    }
    virtual void SubResidual(const unsigned int ind, bool is_rb, const double* res)
    {
      unsigned int iblk = (is_rb) ? ind : ind+nRB;
      unsigned int len = m_residual.Len(iblk);
      for(unsigned int i=0;i<len;i++){
        m_residual.AddValue(iblk,i,res[i]*-1);
      }
    }
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans )
    {
      const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
      const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
      const unsigned int lencol = m_mat.LenBlkCol(iblk0);
      const unsigned int lenrow = m_mat.LenBlkRow(iblk1);
      assert( lencol*lenrow <= 36 );
      double tmp[36];
      for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
      if( isnt_trans ){
        for(unsigned int i=0;i<3;i++){
          for(unsigned int j=0;j<3;j++){
            tmp[(i+offsetr)*lenrow + (j+offsetl)] += m.mat[i*3+j]*d;
          }
        }
      }
      else{
        for(unsigned int i=0;i<3;i++){
          for(unsigned int j=0;j<3;j++){
            tmp[(i+offsetr)*lenrow + (j+offsetl)] += m.mat[j*3+i]*d;
          }
        }
      }
      m_mat.Mearge(1,&iblk0, 1, &iblk1, lencol*lenrow, tmp );
    }
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr, unsigned int sizer,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl, unsigned int sizel,
                           const double* emat, double d){
      const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
      const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
      const unsigned int lencol = m_mat.LenBlkCol(iblk0);
      const unsigned int lenrow = m_mat.LenBlkRow(iblk1);
      assert( lencol*lenrow <= 36 );
      double tmp[36];
      for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
      for(unsigned int i=0;i<sizer;i++){
        for(unsigned int j=0;j<sizel;j++){
          tmp[(i+offsetr)*lenrow + (j+offsetl)] += emat[i*sizel+j]*d;
        }
      }
      m_mat.Mearge(1,&iblk0, 1,&iblk1, lencol*lenrow, tmp );
    }
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                                  unsigned int indl, bool is_rb_l, unsigned int offsetl,
                                  const Com::CVector3D& vec, double d, bool is_column)
    {
      const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
      const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
      const unsigned int lencol = m_mat.LenBlkCol(iblk0);
      const unsigned int lenrow = m_mat.LenBlkRow(iblk1);
      assert( lencol*lenrow <= 36 );
      double tmp[36];
      for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
      if( is_column ){
        tmp[ (0+offsetr)*lenrow + (offsetl)] += vec.x*d;
        tmp[ (1+offsetr)*lenrow + (offsetl)] += vec.y*d;
        tmp[ (2+offsetr)*lenrow + (offsetl)] += vec.z*d;
      }
      else{
        tmp[ (offsetr)*lenrow + (0+offsetl)] += vec.x*d;
        tmp[ (offsetr)*lenrow + (1+offsetl)] += vec.y*d;
        tmp[ (offsetr)*lenrow + (2+offsetl)] += vec.z*d;
      }
      m_mat.Mearge(1,&iblk0, 1, &iblk1, lencol*lenrow, tmp );
    }
    virtual void InitializeMarge(){
      m_mat.SetZero();
      m_update.SetVectorZero();
      m_residual.SetVectorZero();
    }
    virtual double FinalizeMarge(){
      double norm_res = m_residual.GetSquaredVectorNorm();
      return sqrt( norm_res );
    }
    ////////////////////////////////////////////////////////////////
    virtual unsigned int GetTmpVectorArySize() const{ return m_aTmpVec.size(); }
    virtual bool ReSizeTmpVecSolver(unsigned int isize){
      if( m_aTmpVec.size() > isize ){ m_aTmpVec.resize(isize); }
      else if( m_aTmpVec.size() < isize ){
        const unsigned int isize0 = m_aTmpVec.size();
        m_aTmpVec.resize(isize);
        for(unsigned int i=isize0;i<isize;i++){
          m_aTmpVec[i].Initialize(nRB+nConst,m_aBlkSize);
        }
      }
      return true; 
    }
    virtual double DOT(int iv1,int iv2); // {v2}*{v1}
    virtual bool COPY(int iv1,int iv2);  // {v2} := {v1}
    virtual bool SCAL(double d,int iv);  // {v1} := alpha * {v1}
    virtual bool AXPY(double d,int iv1,int iv2); // {v2} := alpha*{v1} + {v2}
    virtual bool MATVEC(double a,int iv1,double b,int iv2);  //  {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
    
    // -1:res  -2:upd
    MatVec::CVector_Blk& GetVector(unsigned int iv){
	    if( iv >= 0 && iv < (int)this->GetTmpVectorArySize() ) return m_aTmpVec[iv];
	    else if( iv == -1 ) return this->m_residual;
      assert( iv == -2 );
      return this->m_update;
    }
    const MatVec::CMatDia_BlkCrs& GetMatrix() const { return m_mat; }
    
    virtual bool UpdateValueOfRigidSystem_NewmarkBetaAPrime(std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
                                                            double dt, double newmark_gamma, double newmark_beta, 
                                                            bool is_first) const;
    virtual bool UpdateValueOfRigidSystem_BackwardEular(std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
                                                        double dt, 
                                                        bool is_first) const;  
    ////////////////////////////////////////////////////////////////
  private:
    
  private:
    unsigned int nRB;
    unsigned int nConst;
    // nRB+nConstが行や列のブロックの数
    
    std::vector<unsigned int> m_aBlkSize;
    MatVec::CMatDia_BlkCrs m_mat;
    MatVec::CVector_Blk m_residual;
    MatVec::CVector_Blk m_update;
    std::vector< MatVec::CVector_Blk > m_aTmpVec;
  };
  
  
  class CPreconditioner_RigidBody_CRS2
  {
  public:
    CPreconditioner_RigidBody_CRS2(){}
    virtual ~CPreconditioner_RigidBody_CRS2(){}
    void SetLinearSystem(const CLinearSystem_RigidBody_CRS2& ls, int ilev = 0){
      m_mat.MakePattern_Initialize( ls.GetMatrix() );
      m_mat.AddFracPtn(-1);
      m_mat.MakePatternFinalize();
    }
    void SetValue(const CLinearSystem_RigidBody_CRS2& ls){
      m_mat.SetValue( ls.GetMatrix() );
    }
    bool Solve( MatVec::CVector_Blk& vec ) const{
      return m_mat.Solve( vec );
    }
  private:
    MatVec::CMatDiaFrac_BlkCrs m_mat;
  };
  
  //! 前処理行列クラスの抽象クラス
  class CLinearSystemPreconditioner_RigidBody_CRS2: public LsSol::ILinearSystemPreconditioner_Sol
  {
  public:
    CLinearSystemPreconditioner_RigidBody_CRS2(CLinearSystem_RigidBody_CRS2& ls, 
                                               CPreconditioner_RigidBody_CRS2& prec) : ls(ls), prec(prec){}
  public:
    //! ソルバに必要な作業ベクトルの数を得る
    virtual unsigned int GetTmpVectorArySize() const{ return ls.GetTmpVectorArySize(); }
    //! ソルバに必要な作業ベクトルの数を設定d
    virtual bool ReSizeTmpVecSolver(unsigned int size_new){ return ls.ReSizeTmpVecSolver(size_new); }
    
    virtual double DOT(int iv1, int iv2){ return ls.DOT(iv1,iv2); }
    virtual bool COPY(int iv1, int iv2){ return ls.COPY(iv1,iv2); }
    virtual bool SCAL(double alpha, int iv1){ return ls.SCAL(alpha,iv1); }
    virtual bool AXPY(double alpha, int iv1, int iv2){ return ls.AXPY(alpha,iv1,iv2); }
    virtual bool MATVEC(double alpha, int iv1, double beta, int iv2){ return ls.MATVEC(alpha,iv1,beta,iv2); }
    
    virtual bool SolvePrecond(int iv){
      assert( iv == -1 || iv == -2 || (iv>=0 && iv<ls.GetTmpVectorArySize()) );
      MatVec::CVector_Blk& vec = ls.GetVector(iv);
      return prec.Solve( vec );
      return true;
    }
  private:
    CPreconditioner_RigidBody_CRS2& prec;
    CLinearSystem_RigidBody_CRS2& ls;
  };
  
  
}

#endif