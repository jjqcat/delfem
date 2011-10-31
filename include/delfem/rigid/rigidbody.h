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
 @brief çÑëÃÉNÉâÉX(Com::CRigidBody3D)ÇÃé¿ëï
 @author Nobuyuki Umetani
 */

#if !defined(RIGID_BODY_H)
#define RIGID_BODY_H

#include <vector>
#include <cassert>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/vector3d.h"
#include "delfem/matrix3d.h"

namespace Ls{
  class CLinearSystem_RigidBody;
}

namespace Rigid
{
  
  class CRigidBody3D
  {
  public:
    CRigidBody3D();
    CRigidBody3D(const CRigidBody3D& rb);
    unsigned int GetDOF() const { return 6; }
    ////
    const Com::CVector3D& GetIniPosCG() const { return ini_pos_cg; }    
    void SetIniPosCG( const Com::CVector3D& pos){ ini_pos_cg = pos; }
    const Com::CVector3D& GetDispCG() const { return disp_cg; }
    void AddDispCG(const Com::CVector3D& d){ disp_cg += d; }
    const Com::CVector3D& GetVeloCG() const { return velo_cg; }
    void SetVeloCG( const Com::CVector3D& vcg ){ velo_cg = vcg; }    
    ////
    Com::CMatrix3 GetRotMatrix() const;
    Com::CMatrix3 GetRotMatrix_Inc() const;
    Com::CMatrix3 GetRotMatrix_Ini() const;
    void AddRotation_IncRight( const Com::CMatrix3& rot );    
    void AddRotation_IniLeft( const Com::CMatrix3& rot );        
    void SetIniRotation( const Com::CMatrix3& rot );
    void SetIncRotation( const Com::CMatrix3& rot );    
    const Com::CVector3D& GetOmega() const { return Omega; }
    void SetOmega( const Com::CVector3D& omg ){ Omega = omg; }
    const Com::CVector3D& GetdOmega() const { return dOmega; }    
    ////
    double GetMass() const { return mass; }
    void SetMass(double mass, const double mineatia[3]){
      this->mass = mass;
      this->mineatia[0] = mineatia[0];
      this->mineatia[1] = mineatia[1];
      this->mineatia[2] = mineatia[2];
    }
    
    ///
    void SetModeBox(double xlen, double ylen, double zlen, double rho ){
      imode = 1;
      this->xlen = xlen;
      this->ylen = ylen;
      this->zlen = zlen;
      this->mass = xlen*ylen*zlen*rho;
      this->mineatia[0] = 1.0/12.0*(ylen*ylen+zlen*zlen)*mass;
      this->mineatia[1] = 1.0/12.0*(zlen*zlen+xlen*xlen)*mass;
      this->mineatia[2] = 1.0/12.0*(xlen*xlen+ylen*ylen)*mass;
    }
    
    void Clear(){
      crv_inc[0] = 0; crv_inc[1] = 0; crv_inc[2] = 0;
      Omega.SetVector(0,0,0);
      dOmega.SetVector(0,0,0);
      velo_cg.SetVector(0,0,0);
      disp_cg.SetVector(0,0,0);
      acc_cg.SetVector(0,0,0);
    }
    Com::CVector3D GetCurrentPositionFromInitial(const Com::CVector3D& vec) const
    {
      const Com::CMatrix3& rot = this->GetRotMatrix_Inc();
      return rot.MatVec(vec-ini_pos_cg)+ini_pos_cg+disp_cg;
    }
    Com::CVector3D GetCurrentPosition(const Com::CVector3D& vec) const  // input is a relative position
    {
      return this->GetRotMatrix().MatVec(vec)+ini_pos_cg+disp_cg;
    }
    Com::CVector3D GetInitialPosition(const Com::CVector3D& vec) const  // input is a relative position
    {
      return this->GetRotMatrix_Ini().MatVec(vec)+ini_pos_cg;
    }    
    Com::CVector3D GetCurrentVelocity(const Com::CVector3D& vec) const
    {
      const Com::CVector3D& Qa = this->GetRotMatrix_Ini().MatVec(vec);
      return velo_cg + this->GetRotMatrix_Inc().MatVec( Cross(this->Omega,Qa) );
    }
    void UpdateSolution_NewmarkBetaAPrime(const double* upd,
                                          const double dt, const double newmark_gamma, const double newmark_beta, 
                                          bool is_first_iter);
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int irb,
                                           const double dt, const double newmark_gamma, const double newmark_beta,
                                           const Com::CVector3D& gravity, 
                                           bool is_first);
    
    void UpdateSolution_BackwardEular(const double* upd,
                                      const double dt, 
                                      bool is_first_iter);
    void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int irb,
                                       const double dt, 
                                       const Com::CVector3D& gravity, 
                                       bool is_first);    
    double GetKineticEnergy() const{
      double e = 0;
      e += 0.5*( Omega.x*Omega.x*mineatia[0]
                +Omega.y*Omega.y*mineatia[1]
                +Omega.z*Omega.z*mineatia[2]);
      e += 0.5*Com::Dot(velo_cg,velo_cg)*mass;      
      return e;
    }
    
    void Draw(bool is_initial) const;
  public:
    double crv_ini[3]; // initial_rotation    
    double crv_inc[3]; // incremental_rotation
    Com::CVector3D Omega;
    Com::CVector3D dOmega;
    
    // center of gravity
    Com::CVector3D ini_pos_cg;
    Com::CVector3D disp_cg;    
    Com::CVector3D velo_cg;
    Com::CVector3D acc_cg;

    double mass;
  public:
    unsigned int imode;
    double xlen, ylen, zlen;
    double mineatia[3];    
  };
  
  ////////////////////////////////////////////////////////////////
  
  class CConstraint
  {
  public:
    virtual unsigned int GetDOF() const = 0;
    virtual void Clear() = 0;
    virtual void UpdateLambda_NewmarkBetaAPrime(const double* upd, const double dt, const double newmark_gamma, const double newmark_beta) = 0;
    virtual void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                                   const double dt, const double newmark_gamma, const double newmark_beta,
                                                   const std::vector<CRigidBody3D>& aRB, 
                                                   bool is_initial ) const = 0;
    
    virtual void UpdateLambda_BackwardEular(const double* upd, double dt) = 0;
    virtual void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                               const double dt, 
                                               const std::vector<CRigidBody3D>& aRB, 
                                               bool is_initial ) const = 0;    
    const std::vector<unsigned int>& GetAry_IndexRB() const { return aIndRB; }
  public:
    std::vector<unsigned int> aIndRB;
  };
  
  
  
  class CFix_Spherical : public CConstraint
  {
  public:
    CFix_Spherical(unsigned int irb){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
      ini_pos_fix.SetVector(0,0,0);
      aIndRB.push_back(irb);
    }
    virtual unsigned int GetDOF() const { return 3; }
    virtual void Clear(){ 
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
    }
    void SetIniPosFix(double x, double y, double z){ ini_pos_fix.SetVector(x,y,z); }
    void UpdateLambda_NewmarkBetaAPrime(const double* upd, const double dt, const double newmark_gamma, const double newmark_beta);           
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                         const double dt, const double newmark_gamma, const double newmark_beta,
                         const std::vector<CRigidBody3D>& aRB, bool is_initial = false ) const;
    
    void UpdateLambda_BackwardEular(const double* upd, double dt){}
    virtual void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                               const double dt, 
                                               const std::vector<CRigidBody3D>& aRB, 
                                               bool is_initial ) const{}
  public:
    double lambda[3];
    Com::CVector3D ini_pos_fix;
  private:
  };
  
  class CFix_Hinge : public CConstraint
  {
  public:
    CFix_Hinge(unsigned int irb){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0;
      axis.SetVector(1,0,0);
      loc_coord[0].SetVector(0,0,0);
      loc_coord[1].SetVector(0,0,0);
      ini_pos_fix.SetVector(0,0,0);
      aIndRB.push_back(irb);
    }
    virtual unsigned int GetDOF() const { return 5; }
    virtual void Clear(){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;  lambda[3] = 0;  lambda[4] = 0;
    }
    void SetIniPosFix(double x, double y, double z){
      this->ini_pos_fix.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void UpdateLambda_NewmarkBetaAPrime(const double* upd, const double dt, const double newmark_gamma, const double newmark_beta);           
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                         const double dt, const double newmark_gamma, const double newmark_beta,
                         const std::vector<CRigidBody3D>& aRB, bool is_initial = false ) const;
    
    void UpdateLambda_BackwardEular(const double* upd, double dt){}  
    virtual void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                               const double dt, 
                                               const std::vector<CRigidBody3D>& aRB, 
                                               bool is_initial ) const{}
  public:
    double lambda[5];
    Com::CVector3D ini_pos_fix;
    Com::CVector3D loc_coord[2];
  private:
    Com::CVector3D axis;
  };
  /*
  class CFix_HingeRange : public CConstraint
  {
  public:
    CFix_HingeRange(unsigned int irb){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0; lambda[5] = 0;
      axis.SetVector(1,0,0);
      loc_coord[0].SetVector(0,0,0);
      loc_coord[1].SetVector(0,0,0);
      ini_pos_fix.SetVector(0,0,0);
      aIndRB.push_back(irb);
      min_t = -180;
      max_t =  180;
    }
    virtual unsigned int GetDOF() const { return 6; }
    virtual void Clear(){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;  lambda[3] = 0;  lambda[4] = 0; lambda[5] = 0;
    }
    void SetIniPosFix(double x, double y, double z){
      this->ini_pos_fix.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void SetRange(double min_t, double max_t){
      this->min_t = min_t;
      this->max_t = max_t;
    }
    void UpdateLambda(const double* upd, double scale);
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                         const double dt, const double newmark_gamma, const double newmark_beta,
                         const std::vector<CRigidBody3D>& aRB, 
                         bool is_initial = false ) const;
  public:
    double lambda[6];
    Com::CVector3D ini_pos_fix;
    Com::CVector3D loc_coord[2];
    double min_t;
    double max_t;
  private:
    Com::CVector3D axis;
  };*/
  
  class CJoint_Spherical : public CConstraint
  {
  public:
    CJoint_Spherical(unsigned int irb0, unsigned int irb1){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
      ini_pos_joint.SetVector(0,0,0);
      aIndRB.push_back(irb0);
      aIndRB.push_back(irb1);
    }
    virtual unsigned int GetDOF() const { return 3; }
    virtual void Clear(){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
    }
    void SetIniPosJoint(double x, double y, double z){
      this->ini_pos_joint.SetVector(x,y,z);
    }
    void UpdateLambda_NewmarkBetaAPrime(const double* upd, const double dt, const double newmark_gamma, const double newmark_beta);           
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                         const double dt, const double newmark_gamma, const double newmark_beta,
                         const std::vector<CRigidBody3D>& aRB, 
                         bool is_initial = false) const;
    
    void UpdateLambda_BackwardEular(const double* upd, double dt){}    
    virtual void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                               const double dt, 
                                               const std::vector<CRigidBody3D>& aRB, 
                                               bool is_initial ) const{}
    //    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
  public:
    double lambda[3];
    Com::CVector3D ini_pos_joint;
  private:
  };
  
  class CJoint_Hinge : public CConstraint
  {
  public:
    CJoint_Hinge(unsigned int irb0, unsigned int irb1){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0;
      ini_pos_joint.SetVector(0,0,0);
      aIndRB.push_back(irb0);
      aIndRB.push_back(irb1);
    }
    virtual unsigned int GetDOF() const { return 5; }
    virtual void Clear(){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0;
    }
    void SetIniPosJoint(double x, double y, double z){
      this->ini_pos_joint.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void UpdateLambda_NewmarkBetaAPrime(const double* upd, const double dt, const double newmark_gamma, const double newmark_beta);
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                         const double dt, const double newmark_gamma, const double newmark_beta,
                         const std::vector<CRigidBody3D>& aRB,
                         bool is_initial = false ) const;
    
    void UpdateLambda_BackwardEular(const double* upd, double dt){}
    virtual void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                               const double dt, 
                                               const std::vector<CRigidBody3D>& aRB, 
                                               bool is_initial ) const{}    
  public:    
    double lambda[5];
    Com::CVector3D ini_pos_joint;
    Com::CVector3D loc_coord[2];
  private:
    Com::CVector3D axis;
  };
  
  
  
  
  
  
  class CJoint_HingeRange : public CConstraint
  {
  public:
    CJoint_HingeRange(unsigned int irb0, unsigned int irb1){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0; lambda[5] = 0;
      ini_pos_joint.SetVector(0,0,0);
      aIndRB.push_back(irb0);
      aIndRB.push_back(irb1);
    }
    virtual unsigned int GetDOF() const { return 6; }
    virtual void Clear(){
      lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0; lambda[5] = 0;
    }
    void SetIniPosJoint(double x, double y, double z){
      this->ini_pos_joint.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void SetRange(double min_t, double max_t){
      this->min_t = min_t;
      this->max_t = max_t;
    }
   
    void UpdateLambda_NewmarkBetaAPrime(const double* upd, const double dt, const double newmark_gamma, const double newmark_beta);        
    void AddLinearSystem_NewmarkBetaAPrime(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                         const double dt, const double newmark_gamma, const double newmark_beta,
                         const std::vector<CRigidBody3D>& aRB,
                         bool is_initial = false ) const;    
    
    void UpdateLambda_BackwardEular(const double* upd, double dt){}
    virtual void AddLinearSystem_BackwardEular(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                                               const double dt, 
                                               const std::vector<CRigidBody3D>& aRB, 
                                               bool is_initial ) const{}    
  public:
    double lambda[6];
    Com::CVector3D ini_pos_joint;
    Com::CVector3D loc_coord[2];
    double min_t;
    double max_t;
  private:
    Com::CVector3D axis;
  };
  
  
  
  
}

#endif