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

#include <iostream>

#include "delfem/matrix3d.h"
#include "delfem/rigid/rigidbody.h"
#include "delfem/rigid/linearsystem_rigid.h"

Rigid::CRigidBody3D::CRigidBody3D(const Rigid::CRigidBody3D& rb){
  this->crv_ini[0] = rb.crv_ini[0];
  this->crv_ini[1] = rb.crv_ini[1];
  this->crv_ini[2] = rb.crv_ini[2];      
  
  this->crv_inc[0] = rb.crv_inc[0];
  this->crv_inc[1] = rb.crv_inc[1];
  this->crv_inc[2] = rb.crv_inc[2];
  
  this->Omega = rb.Omega;
  this->dOmega = rb.dOmega;
  
  this->ini_pos_cg = rb.ini_pos_cg;
  this->disp_cg = rb.disp_cg;
  this->velo_cg = rb.velo_cg;
  this->acc_cg = rb.acc_cg;
  
  this->mass = rb.mass;
  this->mineatia[0] = rb.mineatia[0];
  this->mineatia[1] = rb.mineatia[1];
  this->mineatia[2] = rb.mineatia[2];
  
  this->imode = rb.imode;
  this->xlen = rb.xlen;
  this->ylen = rb.ylen;
  this->zlen = rb.zlen;
}

Rigid::CRigidBody3D::CRigidBody3D(){
  this->ini_pos_cg.SetVector(0,0,0);
  this->disp_cg.SetVector(0,0,0);
  this->velo_cg.SetVector(0,0,0);
  this->acc_cg.SetVector(0,0,0);
  
  crv_ini[0]=0;  crv_ini[1]=0;  crv_ini[2]=0;
  crv_inc[0]=0;  crv_inc[1]=0;  crv_inc[2]=0;      
  Omega.SetVector(0,0,0);
  dOmega.SetVector(0,0,0);
  
  this->mineatia[0] = 1.0; 
  this->mineatia[1] = 1.0;
  this->mineatia[2] = 1.0;
  this->mass = 1.0;
  
  imode = 0;  
}



void Rigid::CRigidBody3D::SetIniRotation(const Com::CMatrix3& m ){
  m.GetCRV_RotMatrix(crv_ini);  
}

void Rigid::CRigidBody3D::SetIncRotation( const Com::CMatrix3& m ){
  m.GetCRV_RotMatrix(crv_inc);
}


Com::CMatrix3 Rigid::CRigidBody3D::GetRotMatrix() const
{
  Com::CMatrix3 m0, m1;
  m0.SetRotMatrix_CRV(crv_ini);
  m1.SetRotMatrix_CRV(crv_inc);
  Com::CMatrix3 m = m1.MatMat(m0);
  return m;
}

Com::CMatrix3 Rigid::CRigidBody3D::GetRotMatrix_Inc() const
{
  Com::CMatrix3 m;
  m.SetRotMatrix_CRV(crv_inc);
  return m;
}

Com::CMatrix3 Rigid::CRigidBody3D::GetRotMatrix_Ini() const
{
  Com::CMatrix3 m;
  m.SetRotMatrix_CRV(crv_ini);
  return m;
}


void Rigid::CRigidBody3D::AddRotation_IncRight( const Com::CMatrix3& rot )
{
  Com::CMatrix3 m; m.SetRotMatrix_CRV(crv_inc);
  m = m.MatMat(rot);
  m.GetCRV_RotMatrix(crv_inc);
}

void Rigid::CRigidBody3D::AddRotation_IniLeft( const Com::CMatrix3& rot )
{
  Com::CMatrix3 m; m.SetRotMatrix_CRV(crv_ini);
  m = rot.MatMat(m);
  m.GetCRV_RotMatrix(crv_ini);
}

void Rigid::CRigidBody3D::UpdateSolution_NewmarkBetaAPrime
(const double* upd, 
 double dt, double newmark_gamma, double newmark_beta, 
 bool is_first_iter )
{
  if( is_first_iter ){
    const unsigned int imode = 1;
    if( imode == 0 ){
      disp_cg += velo_cg*dt + acc_cg*0.5*dt*dt;
      velo_cg += acc_cg*dt;
      ////////////////
      {
        Com::CVector3D rot_update = dt*Omega + 0.5*dt*dt*dOmega; 
        Com::CMatrix3 m; m.SetRotMatrix_Cartesian(rot_update);
        AddRotation_IncRight( m );
      }
      Omega += dOmega*dt;
    }
    else{
      disp_cg += velo_cg*dt + acc_cg*(0.5-newmark_beta)*dt*dt;
      velo_cg += acc_cg*(1-newmark_gamma)*dt;
      acc_cg.SetZero();
      ////////////////
      {
        Com::CVector3D rot_update = dt*Omega + (0.5-newmark_beta)*dt*dt*dOmega; 
        Com::CMatrix3 m; m.SetRotMatrix_Cartesian(rot_update);       
        AddRotation_IncRight( m );        
      }
      Omega += dOmega*(1-newmark_gamma)*dt;
      dOmega.SetZero();
    }
    return;
  }
  const Com::CVector3D acc_disp_update( upd[0], upd[1], upd[2] );
  const Com::CVector3D acc_rot_update(  upd[3], upd[4], upd[5] );
  const double dtmp1 = dt*dt*newmark_beta;
  const double dtmp2 = dt*newmark_gamma;
  
	// 変位のアップデート
  disp_cg += dtmp1*acc_disp_update;
  // 速度のアップデート
  velo_cg += dtmp2*acc_disp_update;
  // 加速度のアップデート
  acc_cg += acc_disp_update;
  
  ////////////////////////////////////////////////
	{	// 回転のアップデート
    Com::CVector3D rot_update = dtmp1*acc_rot_update;
    Com::CMatrix3 m; m.SetRotMatrix_Cartesian(rot_update);           
    AddRotation_IncRight( m );    
	}
  const Com::CVector3D oo = Omega;
  const Com::CMatrix3 woo(oo);
  {   // 角速度のアップデート
    Omega += dtmp2*acc_rot_update;
    Omega += dtmp1*woo.MatVec(acc_rot_update);
  }
  {   // 角加速度のアップデート
    const Com::CMatrix3 wdoo(dOmega);
    const Com::CMatrix3 woowoo = woo.MatMat(woo);
    dOmega += acc_rot_update;
    dOmega += dtmp1*( 0.5*wdoo.MatVec(acc_rot_update)-1.0/6.0*woowoo.MatVec(acc_rot_update) );
  }
}


void Rigid::CRigidBody3D::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int irb,
 const double dt, const double newmark_gamma, const double newmark_beta,
 const Com::CVector3D& gravity, 
 bool is_initial)
{
	// 並進の残差
  ls.AddResidual( irb,true,0,  gravity - acc_cg, mass );
	{   // 慣性力から来る回転の残差
    double res[6] = {0,0,0, 0,0,0};
		res[3] = mineatia[0]*dOmega.x;
    res[4] = mineatia[1]*dOmega.y;
    res[5] = mineatia[2]*dOmega.z;
		const double wOmg[9] = { 0,-Omega.z,+Omega.y, +Omega.z,0,-Omega.x,  -Omega.y,+Omega.x, 0 };
    const double jo[3] = { mineatia[0]*Omega.x, mineatia[1]*Omega.y, mineatia[2]*Omega.z };
		for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){
        res[i+3] += wOmg[i*3+j]*jo[j];
      }
    }
    ls.SubResidual(irb,true, res);
  }
	{   // 剛性行列を作る
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity(mass);
    ls.AddMatrix(irb,true,0,  irb,true,0,  mtmp,1, true);
    mtmp.mat[0] = mineatia[0]; mtmp.mat[1] = 0;           mtmp.mat[2] = 0;
    mtmp.mat[3] = 0;           mtmp.mat[4] = mineatia[1]; mtmp.mat[5] = 0;
    mtmp.mat[6] = 0;           mtmp.mat[7] = 0;           mtmp.mat[8] = mineatia[2];
    ls.AddMatrix(irb,true,3,  irb,true,3,  mtmp,1, true);
	}
}



void Rigid::CRigidBody3D::UpdateSolution_BackwardEular
(const double* upd, 
 double dt,
 bool is_first_iter )
{
  assert( is_first_iter );
  
  const Com::CVector3D vel_update( upd[0], upd[1], upd[2] );
  const Com::CVector3D omg_update( upd[3], upd[4], upd[5] );

  velo_cg += vel_update; // update velocity
  disp_cg += dt*velo_cg; // update displacemnet
  
  ////////////////////////////////////////////////
  Omega += omg_update;  
	{	// 回転のアップデート
    Com::CMatrix3 m; m.SetRotMatrix_Cartesian(dt*Omega);
    AddRotation_IncRight( m );
	}
}


void Rigid::CRigidBody3D::AddLinearSystem_BackwardEular
(Ls::CLinearSystem_RigidBody& ls, unsigned int irb,
 const double dt, 
 const Com::CVector3D& gravity, 
 bool is_initial)
{
	// 並進の残差
  ls.AddResidual( irb,true,0,  dt*gravity, mass );
	{   // 慣性力から来る回転の残差
    double res[6] = {0,0,0, 0,0,0};
		res[3] = mineatia[0]*dOmega.x;
    res[4] = mineatia[1]*dOmega.y;
    res[5] = mineatia[2]*dOmega.z;
		const double wOmg[9] = { 0,-Omega.z,+Omega.y, +Omega.z,0,-Omega.x,  -Omega.y,+Omega.x, 0 };
    const double jo[3] = { mineatia[0]*Omega.x, mineatia[1]*Omega.y, mineatia[2]*Omega.z };
		for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){
        res[i+3] += wOmg[i*3+j]*jo[j];
      }
    }
    for(unsigned int i=0;i<6;i++){ res[i] *= dt; }
    ls.SubResidual(irb,true, res);    
  }
	{   // 剛性行列を作る
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity(mass);
    ls.AddMatrix(irb,true,0,  irb,true,0,  mtmp,1, true);
    mtmp.mat[0] = mineatia[0]; mtmp.mat[1] = 0;           mtmp.mat[2] = 0;
    mtmp.mat[3] = 0;           mtmp.mat[4] = mineatia[1]; mtmp.mat[5] = 0;
    mtmp.mat[6] = 0;           mtmp.mat[7] = 0;           mtmp.mat[8] = mineatia[2];
    ls.AddMatrix(irb,true,3,  irb,true,3,  mtmp,1, true);
	}
}

void Rigid::CRigidBody3D::Draw(bool is_initial) const
{  
  ::glPushMatrix();
  if( !is_initial ){
    const Com::CVector3D& v = this->GetIniPosCG() + this->GetDispCG();
    ::glTranslated(v.x,v.y,v.z);
  }
  else{
    const Com::CVector3D& v = this->GetIniPosCG();
    ::glTranslated(v.x,v.y,v.z);
  }
  if( !is_initial ){
    Com::CMatrix3 m = this->GetRotMatrix();
    double rot0[16];
    m.GetAffineTransMatElements(rot0);
    ::glMultMatrixd(rot0);
  }
  else{
    Com::CMatrix3 m = this->GetRotMatrix_Ini();
    double rot0[16];
    m.GetAffineTransMatElements(rot0);
    ::glMultMatrixd(rot0);        
  }
  if( imode == 0 ){        
    ::glScaled(0.1,0.1,0.1);
  }
  else if( imode == 1 ){
    ::glScaled(xlen,ylen,zlen);
  }
  ::glLineWidth(1);      
  {
    ::glColor3d(0,0,0);
    ::glutWireCube(1.0);        
  }
  ////
  ::glScaled(1,1,1);
  ::glBegin(GL_LINES);
  ::glColor3d(1,0,0);
  ::glVertex3d(0,0,0);
  ::glVertex3d(0.4,0,0);
  ::glColor3d(0,1,0);
  ::glVertex3d(0,0,0);
  ::glVertex3d(0,0.4,0);
  ::glColor3d(0,0,1);
  ::glVertex3d(0,0,0);
  ::glVertex3d(0,0,0.4);
  ::glEnd();
  ////
  ::glPopMatrix();      
}


////////////////////////////////////////////////////////////////

void Rigid::CFix_Spherical::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
 const double dt, const double newmark_gamma, const double newmark_beta, 
 const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
  const unsigned int irb = aIndRB[0];
  assert( irb < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb = aRB[irb];
  
  ////////////////////////////////////////////////
  const Com::CMatrix3& mrot = rb.GetRotMatrix_Inc();
  const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
  const Com::CVector3D& Xdistfix = rb.GetIniPosCG() - ini_pos_fix;
  
  // 拘束力残差
	// 並進の残差
  ls.AddResidual(irb,true,0, vlambda, 1);
	{   // 拘束力から来る回転の残差
    Com::CVector3D tmp_v = mrot.MatVecTrans(vlambda);
    Com::CMatrix3 wX(Xdistfix);
    const Com::CVector3D& vec_q = wX.MatVec(tmp_v);
    ls.AddResidual(irb,true,3, vec_q, -1);
	}
	{   // 拘束条件の変分
    const Com::CVector3D& rot_pos_cg = mrot.MatVec(Xdistfix) + ini_pos_fix;
    ls.AddResidual(icst,false,0,  rb.GetIniPosCG()+rb.GetDispCG()-rot_pos_cg, 1 );
	}
  ////////////////////////////////
  Com::CMatrix3 RotwX = mrot.MatMat( Com::CMatrix3(Xdistfix) );
  Com::CMatrix3 wXwRtL;
	{
    const Com::CVector3D RtL = mrot.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix);
    wXwRtL = wX.MatMat( Com::CMatrix3(RtL) );
	}   
	{   // 剛性行列を作る
		const double dtmp1 = dt*dt*newmark_beta;
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity(-dtmp1);
    ls.AddMatrix(irb, true, 0, icst,false,0, mtmp,   1,     true );
    ls.AddMatrix(icst,false,0, irb, true, 0, mtmp,   1,     true );
    ls.AddMatrix(irb, true, 3, irb, true, 3, wXwRtL, dtmp1, true );
    ls.AddMatrix(irb, true, 3, icst,false,0, RotwX, -dtmp1, false);
    ls.AddMatrix(icst,false,0, irb, true, 3, RotwX, -dtmp1, true );
	}
}

void Rigid::CFix_Spherical::UpdateLambda_NewmarkBetaAPrime
(const double* upd, double dt, double newmark_gamma, double newmark_beta ){
  double tmp = dt*dt*newmark_beta;
  lambda[0] += tmp*upd[0];
  lambda[1] += tmp*upd[1];
  lambda[2] += tmp*upd[2];
}

/*
 void Rigid::CFix_Spherical::Draw(const std::vector<CRigidBody3D>& aRB) const
 {
 const unsigned int imode = 2;
 if( imode == 0 ){
 }
 else if( imode == 1 ){
 ::glColor3d(1,1,0);
 ::glBegin(GL_POINTS);
 ::glVertex3d(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
 ::glEnd();
 }
 else if( imode == 2 ){
 const unsigned int irb = aIndRB[0];
 const CRigidBody3D& rb = aRB[irb];
 const Com::CVector3D& vec_j = rb.GetPositionFromInital(ini_pos_fix);
 const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.GetDispCG();
 
 ::glPushMatrix();
 ::glTranslated(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
 ::glColor3d(1,1,0);
 ::glutSolidSphere(0.1,10,5);
 ::glPopMatrix();
 
 ::glColor3d(1,1,1);
 ::glLineWidth(2);
 ::glBegin(GL_LINES);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
 ::glVertex3d(vec_cg.x,vec_cg.y,vec_cg.z);
 ::glEnd();
 }
 }
 */

////////////////////////////////////////////////////////////////

void Rigid::CFix_Hinge::SetAxis(double ax, double ay, double az){
  axis = Com::CVector3D(ax,ay,az);
  if( axis.Length() < 1.0e-20 ){ assert(0); return; }
  axis.SetNormalizedVector();
  GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

void Rigid::CFix_Hinge::UpdateLambda_NewmarkBetaAPrime
(const double* upd, double dt, double newmark_gamma, double newmark_beta ){
  double tmp = dt*dt*newmark_beta;
  lambda[0] += tmp*upd[0];
  lambda[1] += tmp*upd[1];
  lambda[2] += tmp*upd[2];
  lambda[3] += tmp*upd[3];
  lambda[4] += tmp*upd[4];
}

void Rigid::CFix_Hinge::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
 const double dt, const double newmark_gamma, const double newmark_beta,
 const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
  const unsigned int irb = aIndRB[0];
  assert( irb < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb = aRB[irb];
  
  const Com::CMatrix3& mrot = rb.GetRotMatrix();
  const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
  const Com::CVector3D& Xdistfix = rb.GetIniPosCG() - ini_pos_fix;
  
  ////////////////////////////////////////////////
  // 拘束力残差
	// 並進の残差
  ls.AddResidual(irb,true,0, vlambda, 1);
	{   // 位置拘束力から来る物体回転の残差
    const Com::CVector3D& tmp_vec = mrot.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix);
    ls.AddResidual(irb,true,3,  wX.MatVec(tmp_vec), -1 );
	}
	{   // 拘束所条件の変分
    double res[5] = {0,0,0,0,0};
    // 変位拘束条件の変分
    Com::CVector3D rot_pos_cg = mrot.MatVec(Xdistfix) + ini_pos_fix;
    Com::CVector3D res_disp = - rb.GetIniPosCG() - rb.GetDispCG() + rot_pos_cg;
    res[0] = res_disp.x; 
    res[1] = res_disp.y; 
    res[2] = res_disp.z;
    // 回転拘束条件の変分
    const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
    res[3] = Com::Dot(Rta,loc_coord[0]);
    res[4] = Com::Dot(Rta,loc_coord[1]);
    ls.SubResidual(icst,false,res);
	}
  {   // 回転高速から来る物体回転の残差
    const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
    const Com::CMatrix3 wTmp(lambda[3]*loc_coord[0] + lambda[4]*loc_coord[1]);
    ls.AddResidual(irb,true,3, wTmp.MatVec(Rta),-1 );
  }
  ////////////////////////////////
  Com::CMatrix3 RotwX = mrot.MatMat( Com::CMatrix3(Xdistfix) );
  Com::CMatrix3 wXwRtL;
	{
    const Com::CVector3D RtL = mrot.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix);
    wXwRtL = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  Com::CVector3D wm0Rta, wm1Rta;
  Com::CMatrix3 lwLCwRta;
  {
    const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
    const Com::CMatrix3 wm0(loc_coord[0]);
    const Com::CMatrix3 wm1(loc_coord[1]);
    wm0Rta = wm0.MatVec(Rta);
    wm1Rta = wm1.MatVec(Rta);
    Com::CMatrix3 wTmp(lambda[3]*loc_coord[0]+lambda[4]*loc_coord[1]);
    lwLCwRta = wTmp.MatMat( Com::CMatrix3(Rta) );
  }
	{   // 剛性行列を作る
		const double dtmp1 = dt*dt*newmark_beta;
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity();
    ls.AddMatrix(irb, true, 0, icst,false,0, mtmp,   -dtmp1, true);
    ls.AddMatrix(icst,false,0, irb, true, 0, mtmp,   -dtmp1, true);
    
    ls.AddMatrix(irb, true, 3, irb, true, 3, wXwRtL,  dtmp1, true);
    ls.AddMatrix(irb, true, 3, irb, true, 3, lwLCwRta,dtmp1, true);
    
    ls.AddMatrix(irb, true, 3, icst,false,0, RotwX,  -dtmp1, false);
    ls.AddMatrix(icst,false,0, irb, true, 3, RotwX,  -dtmp1, true );
    
    ls.AddMatrix_Vector(irb, true, 3, icst,false,3, wm0Rta,dtmp1, true );
    ls.AddMatrix_Vector(icst,false,3, irb, true, 3, wm0Rta,dtmp1, false);
    ls.AddMatrix_Vector(irb, true, 3, icst,false,4, wm1Rta,dtmp1, true );
    ls.AddMatrix_Vector(icst,false,4, irb, true, 3, wm1Rta,dtmp1, false);
	}
}
/*
 void Rigid::CFix_Hinge::Draw(const std::vector<CRigidBody3D>& aRB) const
 {
 const unsigned int imode = 1;
 if( imode == 0 ){
 }
 else if( imode == 1 ){
 const unsigned int irb = aIndRB[0];
 const CRigidBody3D& rb = aRB[irb];
 const Com::CVector3D& vec_j = rb.GetPositionFromInital(ini_pos_fix);
 const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.GetDispCG();
 
 ::glPushMatrix();
 ::glTranslated(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
 ::glColor3d(1,1,0);
 ::glutSolidSphere(0.1,10,5);
 ::glPopMatrix();
 
 ::glColor3d(1,1,1);
 ::glLineWidth(2);
 ::glBegin(GL_LINES);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 ::glVertex3d(vec_cg.x,vec_cg.y,vec_cg.z);
 ::glEnd();
 
 const Com::CVector3D& lcb0 = this->loc_coord[0];
 const Com::CVector3D& lcb1 = this->loc_coord[1];
 unsigned int ndiv = 16;
 const double dtheta = 2*3.1416/ndiv;
 const double radius = 0.4;
 ::glColor3d(0,1,1);
 ::glBegin(GL_TRIANGLE_FAN);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 for(unsigned int idiv=0;idiv<ndiv+1;idiv++){
 const Com::CVector3D v0 = vec_j + sin(idiv*dtheta  )*lcb0*radius + cos(idiv*dtheta  )*lcb1*radius;
 ::glVertex3d(v0.x,  v0.y,  v0.z);
 }
 ::glEnd();
 }
 }
 */

////////////////////////////////////////////////////////////////
/*
void Rigid::CFix_HingeRange::SetAxis(double ax, double ay, double az){
  axis = Com::CVector3D(ax,ay,az);
  if( axis.Length() < 1.0e-20 ){ assert(0); return; }
  axis.Normalize();
  GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

void Rigid::CFix_HingeRange::UpdateLambda
(const double* upd, double scale )
{
  lambda[0] += scale*upd[0];
  lambda[1] += scale*upd[1];
  lambda[2] += scale*upd[2];
  lambda[3] += scale*upd[3];
  lambda[4] += scale*upd[4];
  const double tmp2 = dt*newmark_gamma;
  lambda[5] += tmp2*upd[5];
}

void Rigid::CFix_HingeRange::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
 const double dt, double newmark_gamma, double newmark_beta,
 const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
  const unsigned int irb = aIndRB[0];
  assert( irb < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb = aRB[irb];
  
  const Com::CMatrix3& mrot = rb.GetRotMatrix();
  const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
  const Com::CVector3D& Xdistfix = rb.GetIniPosCG() - ini_pos_fix;
  
  ////////////////////////////////////////////////
  // 拘束力残差
	// 並進の残差
  ls.AddResidual(irb,true,0, vlambda, 1);
	{   // 位置拘束力から来る物体回転の残差
    const Com::CVector3D& tmp_vec = mrot.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix);
    ls.AddResidual(irb,true,3,  wX.MatVec(tmp_vec), -1 );
	}
	{   // 拘束所条件の変分
    double res[6] = {0,0,0,0,0,0};
    // 変位拘束条件の変分
    Com::CVector3D rot_pos_cg = mrot.MatVec(Xdistfix) + ini_pos_fix;
    Com::CVector3D res_disp = - rb.GetIniPosCG() - rb.GetDispCG() + rot_pos_cg;
    res[0] = res_disp.x; 
    res[1] = res_disp.y; 
    res[2] = res_disp.z;
    // 回転拘束条件の変分
    const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
    res[3] = Com::Dot(Rta,loc_coord[0]);
    res[4] = Com::Dot(Rta,loc_coord[1]);
    res[5] = 0;
    ls.SubResidual(icst,false,res);
	}
  {   // 回転高速から来る物体回転の残差
    const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
    const Com::CMatrix3 wTmp(lambda[3]*loc_coord[0] + lambda[4]*loc_coord[1]);
    ls.AddResidual(irb,true,3, wTmp.MatVec(Rta),-1 );
  }
  ////////////////////////////////
  Com::CMatrix3 RotwX = mrot.MatMat( Com::CMatrix3(Xdistfix) );
  Com::CMatrix3 wXwRtL;
	{
    const Com::CVector3D RtL = mrot.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix);
    wXwRtL = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  Com::CVector3D wm0Rta, wm1Rta;
  Com::CMatrix3 lwLCwRta;
  {
    const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
    const Com::CMatrix3 wm0(loc_coord[0]);
    const Com::CMatrix3 wm1(loc_coord[1]);
    wm0Rta = wm0.MatVec(Rta);
    wm1Rta = wm1.MatVec(Rta);
    Com::CMatrix3 wTmp(lambda[3]*loc_coord[0]+lambda[4]*loc_coord[1]);
    lwLCwRta = wTmp.MatMat( Com::CMatrix3(Rta) );
  }
  {   // 角度についてのKKT条件を入れる
    double m0Rm0 = Com::Dot(loc_coord[0],mrot.MatVec(loc_coord[0]));
    double m1Rm0 = Com::Dot(loc_coord[1],mrot.MatVec(loc_coord[0]));
    const double PI = 3.14159265358979323846;
    double u_min[2] = { cos(PI*min_t/180.0), sin(PI*min_t/180.0) };
    double u_max[2] = { cos(PI*max_t/180.0), sin(PI*max_t/180.0) };
    double v_min[2] = { -u_min[1],  u_min[0] };
    double v_max[2] = {  u_max[1], -u_max[0] };
    bool flg = false;
    double v[2];
    double f_value;
    const double delta_f = 1.0e-5;
    {
      double dmin = v_min[0]*m0Rm0 + v_min[1]*m1Rm0;
      double dmax = v_max[0]*m0Rm0 + v_max[1]*m1Rm0;
      //            std::cout << dmin << " " << dmax << " " << lambda[5] << std::endl;
      if( dmin < -delta_f || dmax < -delta_f ) flg = true;
      const double emin = u_min[0]*m0Rm0 + u_min[1]*m1Rm0;
      const double emax = u_max[0]*m0Rm0 + u_max[1]*m1Rm0;
      if( emin > emax ){ 
        v[0] = v_min[0]; v[1] = v_min[1];
        f_value = dmin;
      }
      else{                              
        v[0] = v_max[0]; v[1] = v_max[1];
        f_value = dmax;
      }
    }
    //		const double dtmp1 = dt*dt*newmark_beta;
    const double dtmp2 = dt*newmark_gamma;
    if( flg == true && is_initial ){
      std::cout << "consider contact" << std::endl;
      double df;
      Com::CVector3D Rtm = mrot.MatVecTrans( v[0]*loc_coord[0] + v[1]*loc_coord[1] );
      Com::CMatrix3 wm0( loc_coord[0] );
      Com::CVector3D wm0Rtm = wm0.MatVec( Rtm );
      Com::CVector3D res_t = lambda[5]*wm0Rtm;
      df = Com::Dot(wm0Rtm,rb.GetOmega());
      //            std::cout << " velo f : " << df << std::endl;
      ls.AddResidual(irb, true, 3,   res_t,-1 );
      ls.AddResidual(icst,false,5,1, &df,-1 );
      ls.AddMatrix_Vector(irb, true, 3, icst,false,5, wm0Rtm,dtmp2, true );
      ls.AddMatrix_Vector(icst,false,5, irb, true, 3, wm0Rtm,dtmp2, false);
    }
    else if( !is_initial && lambda[5]<-1.0e-5 ){
      //            std::cout << "after impact" << std::endl;
      Com::CVector3D Rtm = mrot.MatVecTrans( v[0]*loc_coord[0] + v[1]*loc_coord[1] );
      Com::CMatrix3 wm0( loc_coord[0] );
      Com::CVector3D wm0Rtm = wm0.MatVec( Rtm );
      Com::CVector3D res_t = lambda[5]*wm0Rtm;
      ls.AddResidual(irb, true, 3,   res_t,-1 );
      ////////////////
      double val0 = 0.0;
      ls.AddResidual(icst,false,5,1, &val0,-1 );
      double val1 = 1.0;
      ls.AddMatrix(icst,false,5,1,  icst,false,5,1, &val1, dtmp2);
    }
    else{
      //            std::cout << "not contact" << std::endl;
      const double val0 = 1;
      const double val1 = lambda[5];
      ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp2);
      ls.AddResidual(icst,false,5,1, &val1,-1 );
    }
  }
	{   // 剛性行列を作る
		const double dtmp1 = dt*dt*newmark_beta;
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity();
    ls.AddMatrix(irb, true, 0, icst,false,0, mtmp,   -dtmp1, true);
    ls.AddMatrix(icst,false,0, irb, true, 0, mtmp,   -dtmp1, true);
    
    ls.AddMatrix(irb, true, 3, irb, true, 3, wXwRtL,  dtmp1, true);
    ls.AddMatrix(irb, true, 3, irb, true, 3, lwLCwRta,dtmp1, true);
    
    ls.AddMatrix(irb, true, 3, icst,false,0, RotwX,  -dtmp1, false);
    ls.AddMatrix(icst,false,0, irb, true, 3, RotwX,  -dtmp1, true );
    
    ls.AddMatrix_Vector(irb, true, 3, icst,false,3, wm0Rta,dtmp1, true );
    ls.AddMatrix_Vector(icst,false,3, irb, true, 3, wm0Rta,dtmp1, false);
    ls.AddMatrix_Vector(irb, true, 3, icst,false,4, wm1Rta,dtmp1, true );
    ls.AddMatrix_Vector(icst,false,4, irb, true, 3, wm1Rta,dtmp1, false);
    
    //        ls.AddMatrix_Vector(icst,false,5, irb, true, 5, 1,dtmp1, false);
	}
}
 */
/*
 void Rigid::CFix_HingeRange::Draw(const std::vector<CRigidBody3D>& aRB) const
 {
 const unsigned int imode = 1;
 if( imode == 0 ){
 }
 else if( imode == 1 ){
 const unsigned int irb = aIndRB[0];
 const CRigidBody3D& rb = aRB[irb];
 const Com::CVector3D& vec_j = rb.GetPositionFromInital(ini_pos_fix);
 const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.GetDispCG();
 
 ::glPushMatrix();
 ::glTranslated(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
 ::glColor3d(1,1,0);
 ::glutSolidSphere(0.1,10,5);
 ::glPopMatrix();
 
 ::glColor3d(1,1,1);
 ::glLineWidth(2);
 ::glBegin(GL_LINES);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 ::glVertex3d(vec_cg.x,vec_cg.y,vec_cg.z);
 ::glEnd();
 
 const Com::CVector3D& lcb0 = this->loc_coord[0];
 const Com::CVector3D& lcb1 = this->loc_coord[1];
 unsigned int ndiv_t = 32;
 unsigned int ndiv0 = ndiv_t*(max_t-min_t      )/360.0 + 1;
 const double dtheta0 = 2*3.1416/ndiv0*(max_t-min_t)/360.0;
 const double radius = 1;
 ::glColor3d(0,1,1);
 ::glBegin(GL_TRIANGLE_FAN);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 for(unsigned int idiv=0;idiv<ndiv0+1;idiv++){
 const Com::CVector3D v0 = vec_j + 
 sin(idiv*dtheta0 - max_t*3.14/180)*lcb0*radius + 
 cos(idiv*dtheta0 - max_t*3.14/180)*lcb1*radius;
 ::glVertex3d(v0.x,  v0.y,  v0.z);
 }
 ::glEnd();
 }
 }
 */
////////////////////////////////////////////////////////////////

void Rigid::CJoint_Spherical::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
 const double dt, const double newmark_gamma, const double newmark_beta, 
 const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
  const unsigned int irb0 = aIndRB[0];
  assert( irb0 < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb0 = aRB[irb0];
  const Com::CMatrix3& mrot0 = rb0.GetRotMatrix();
  const Com::CVector3D& Xdistfix0 = rb0.GetIniPosCG()-ini_pos_joint;
  
  const unsigned int irb1 = aIndRB[1];
  assert( irb1 < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb1 = aRB[irb1];
  const Com::CMatrix3& mrot1 = rb1.GetRotMatrix();
  const Com::CVector3D& Xdistfix1 = rb1.GetIniPosCG()-ini_pos_joint;
  
  const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);  
  
  ////////////////////////////////////////////////
  // 拘束力残差（物体０）
	// 並進の残差
  ls.AddResidual(irb0,true,0,  vlambda,1);
	{   // 拘束力から来る回転の残差（物体０）
    Com::CVector3D tmp_vec = mrot0.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix0);
    ls.AddResidual(irb0,true,3, wX.MatVec(tmp_vec),-1);
	}
	{   // 拘束条件の変分（物体０）
    const Com::CVector3D& rot_pos_cg = mrot0.MatVec(Xdistfix0);
    ls.AddResidual( icst,false,0,   rb0.GetIniPosCG()+rb0.GetDispCG()-rot_pos_cg,1 );
	}
  Com::CMatrix3 RotwX0 = mrot0.MatMat( Com::CMatrix3(Xdistfix0) );
  Com::CMatrix3 wXwRtL0;
	{
    const Com::CVector3D& RtL = mrot0.MatVec(vlambda);
    const Com::CMatrix3 wX(Xdistfix0);
    wXwRtL0 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  ////////////////////////////////////////////////
  // 拘束力残差（物体１）
	// 並進の残差
  ls.AddResidual(irb1,true,0,  vlambda,-1);
	{   // 拘束力から来る回転の残差（物体１）
    const Com::CVector3D& tmp_vec = mrot1.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix1);
    ls.AddResidual(irb1,true,3, wX.MatVec(tmp_vec),1 ); 
	}
	{   // 拘束条件の変分（物体１）
    const Com::CVector3D& rot_pos_cg = mrot1.MatVec(Xdistfix1);
    ls.AddResidual(icst,false,0,  rb1.GetIniPosCG()+rb1.GetDispCG()-rot_pos_cg,-1 );
	}
  Com::CMatrix3 RotwX1 = mrot1.MatMat( Com::CMatrix3(Xdistfix1) );
  Com::CMatrix3 wXwRtL1;
	{
    const Com::CVector3D& RtL = mrot1.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix1);
    wXwRtL1 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
	{   // 剛性行列を作る
		const double dtmp1 = dt*dt*newmark_beta;
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity();
    ls.AddMatrix(irb0,true, 0, icst,false,0, mtmp,   -dtmp1, true);
    ls.AddMatrix(icst,false,0, irb0,true, 0, mtmp,   -dtmp1, true);
    ls.AddMatrix(irb1,true, 0, icst,false,0, mtmp,    dtmp1, true);
    ls.AddMatrix(icst,false,0, irb1,true, 0, mtmp,    dtmp1, true);
    
    ls.AddMatrix(irb0,true,3, irb0,true,3, wXwRtL0,  dtmp1, true);
    ls.AddMatrix(irb1,true,3, irb1,true,3, wXwRtL1, -dtmp1, true);
    
    ls.AddMatrix(irb0,true, 3,  icst,false,0,  RotwX0, -dtmp1, false);
    ls.AddMatrix(icst,false,0,  irb0,true, 3,  RotwX0, -dtmp1, true );
    ls.AddMatrix(irb1,true, 3,  icst,false,0,  RotwX1, +dtmp1, false);
    ls.AddMatrix(icst,false,0,  irb1,true, 3,  RotwX1, +dtmp1, true );    
	}
}

void Rigid::CJoint_Spherical::UpdateLambda_NewmarkBetaAPrime
(const double* upd, double dt, double newmark_gamma, double newmark_beta ){
  double tmp = dt*dt*newmark_beta;
  lambda[0] += tmp*upd[0];
  lambda[1] += tmp*upd[1];
  lambda[2] += tmp*upd[2];
}
/*
 void Rigid::CJoint_Spherical::Draw(const std::vector<CRigidBody3D>& aRB) const
 {
 const unsigned int imode = 1;
 if( imode == 0 ){
 }
 else if( imode == 1 ){
 const unsigned int irb0 = aIndRB[0];
 const CRigidBody3D& rb0 = aRB[irb0];
 const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.GetDispCG();
 
 const unsigned int irb1 = aIndRB[1];
 const CRigidBody3D& rb1 = aRB[irb1];
 const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.GetDispCG();
 
 const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
 
 ::glPushMatrix();
 ::glTranslated( vec_j.x, vec_j.y, vec_j.z );
 ::glColor3d(1,1,0);
 ::glutSolidSphere(0.1,10,5);
 ::glPopMatrix();
 ::glColor3d(1,1,1);
 ::glLineWidth(2);
 ::glBegin(GL_LINES);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
 ::glVertex3d(vec_cg0.x,vec_cg0.y,vec_cg0.z);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 ::glVertex3d(vec_cg1.x,vec_cg1.y,vec_cg1.z);
 ::glEnd();
 }
 }
 */
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////

void Rigid::CJoint_Hinge::SetAxis(double ax, double ay, double az){
  axis = Com::CVector3D(ax,ay,az);
  if( axis.Length() < 1.0e-20 ){ assert(0); return; }
  axis.SetNormalizedVector();
  GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

/*
 void Rigid::CJoint_Hinge::Draw(const std::vector<CRigidBody3D>& aRB) const
 {
 const unsigned int imode = 2;
 if( imode == 0 ){
 }
 else if( imode == 1 ){
 const unsigned int irb0 = aIndRB[0];
 const CRigidBody3D& rb0 = aRB[irb0];
 const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
 ::glColor3d(1,1,0);
 ::glBegin(GL_POINTS);
 ::glVertex3d(vec_j.x, vec_j.y, vec_j.z );
 ::glEnd();
 }
 else if( imode == 2 ){
 const unsigned int irb0 = aIndRB[0];
 const unsigned int irb1 = aIndRB[1];
 const CRigidBody3D& rb0 = aRB[irb0];
 const CRigidBody3D& rb1 = aRB[irb1];
 const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
 const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.GetDispCG();
 const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.GetDispCG();
 
 ::glPushMatrix();
 ::glTranslated( vec_j.x, vec_j.y, vec_j.z );
 ::glColor3d(1,1,0);
 ::glutSolidSphere(0.1,10,5);
 ::glPopMatrix();
 ::glColor3d(1,1,1);
 ::glLineWidth(2);
 ::glBegin(GL_LINES);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
 ::glVertex3d(vec_cg0.x,vec_cg0.y,vec_cg0.z);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 ::glVertex3d(vec_cg1.x,vec_cg1.y,vec_cg1.z);
 ::glEnd();
 
 const Com::CMatrix3& mrot = rb0.GetRotMatrix();
 const Com::CVector3D& lcb0 = mrot.MatVec(this->loc_coord[0]);
 const Com::CVector3D& lcb1 = mrot.MatVec(this->loc_coord[1]);
 unsigned int ndiv = 16;
 const double dtheta = 2*3.1416/ndiv;
 const double radius = 0.4;
 ::glColor3d(0,1,1);
 ::glBegin(GL_TRIANGLE_FAN);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 for(unsigned int idiv=0;idiv<ndiv+1;idiv++){
 const Com::CVector3D v0 = vec_j + sin(idiv*dtheta  )*lcb0*radius + cos(idiv*dtheta  )*lcb1*radius;
 ::glVertex3d(v0.x,  v0.y,  v0.z);
 }
 ::glEnd();
 }
 }
 */

void Rigid::CJoint_Hinge::UpdateLambda_NewmarkBetaAPrime
(const double* upd, double dt, double newmark_gamma, double newmark_beta ){
  double tmp = dt*dt*newmark_beta;
  lambda[0] += tmp*upd[0];
  lambda[1] += tmp*upd[1];
  lambda[2] += tmp*upd[2];
  lambda[3] += tmp*upd[3];
  lambda[4] += tmp*upd[4];
}

void Rigid::CJoint_Hinge::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
 const double dt, const double newmark_gamma, const double newmark_beta,
 const std::vector<CRigidBody3D>& aRB, bool is_initial ) const
{
  const unsigned int irb0 = aIndRB[0];
  assert( irb0 < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb0 = aRB[irb0];
  
  const unsigned int irb1 = aIndRB[1];
  assert( irb1 < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb1 = aRB[irb1];
  
  const Com::CMatrix3& mrot0 = rb0.GetRotMatrix();
  const Com::CMatrix3& mrot1 = rb1.GetRotMatrix();
  const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
  const Com::CVector3D Xdistfix0 = rb0.GetIniPosCG() - ini_pos_joint;
  const Com::CVector3D Xdistfix1 = rb1.GetIniPosCG() - ini_pos_joint;
  
  ////////////////////////////////////////////////
  // 拘束力残差（物体０）
	// 並進の残差
  ls.AddResidual(irb0,true,0,  vlambda,1);
	{   // 変位拘束力から来る回転の残差（物体０）
    const Com::CVector3D& tmp_vec = mrot0.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix0);
    ls.AddResidual(irb0,true,3,  wX.MatVec(tmp_vec),-1);
	}
	{   // 変位拘束条件の変分（物体０）
    Com::CVector3D rot_pos_cg = mrot0.MatVec(Xdistfix0);
    ls.AddResidual(icst,false,0,   rb0.GetIniPosCG()+rb0.GetDispCG()-rot_pos_cg,1 );
	}
  Com::CMatrix3 RotwX0 = mrot0.MatMat( Com::CMatrix3(Xdistfix0) );
  Com::CMatrix3 wXwRtL0;
	{
    const Com::CVector3D& RtL = mrot0.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix0);
    wXwRtL0 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  ////////////////////////////////////////////////
  // 拘束力残差（物体１）
	// 並進の残差
  ls.AddResidual(irb1,true,0,  vlambda,-1);
	{   // 拘束力から来る回転の残差（物体１）
    const Com::CVector3D& tmp_vec = mrot1.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix1);
    ls.AddResidual(irb1,true,3,  wX.MatVec(tmp_vec),1 );
	}
	{   // 拘束条件の変分（物体１）
    const Com::CVector3D& rot_pos_cg = mrot1.MatVec(Xdistfix1);
    ls.AddResidual(icst,false,0,  rb1.GetIniPosCG()+rb1.GetDispCG()-rot_pos_cg,-1);
	}
  Com::CMatrix3 RotwX1 = mrot1.MatMat( Com::CMatrix3(Xdistfix1) );
  Com::CMatrix3 wXwRtL1;
	{
    const Com::CVector3D& RtL = mrot1.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix1);
    wXwRtL1 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  
  ////////////////////////////////
  
  Com::CVector3D Rtwma[2][2];
  {   // 回転拘束から来る物体回転の残差
    const Com::CVector3D& R1a = mrot1.MatVec(axis);
    const Com::CVector3D& R0m0 = mrot0.MatVec(loc_coord[0]);
    const Com::CVector3D& R0m1 = mrot0.MatVec(loc_coord[1]);
    const Com::CMatrix3 wR0m0(R0m0);
    const Com::CMatrix3 wR0m1(R0m1);
    const Com::CVector3D& wR0m0R1a = wR0m0.MatVec(R1a);
    const Com::CVector3D& wR0m1R1a = wR0m1.MatVec(R1a);
    Rtwma[0][0] = mrot0.MatVecTrans(wR0m0R1a);
    Rtwma[0][1] = mrot0.MatVecTrans(wR0m1R1a);
    Rtwma[1][0] = mrot1.MatVecTrans(wR0m0R1a);
    Rtwma[1][1] = mrot1.MatVecTrans(wR0m1R1a);
    ls.AddResidual(irb0,true,3, lambda[3]*Rtwma[0][0] + lambda[4]*Rtwma[0][1],-1 );
    ls.AddResidual(irb1,true,3, lambda[3]*Rtwma[1][0] + lambda[4]*Rtwma[1][1], 1 );
    double res[5] = {0,0,0,0,0};
    res[3] = Com::Dot(R1a,R0m0);
    res[4] = Com::Dot(R1a,R0m1);
    ls.SubResidual(icst,false,res);
  }
  ////////////////////////////////////////////////////////////////
  
	{   // 剛性行列を作る
		const double dtmp1 = dt*dt*newmark_beta;
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity();
    ls.AddMatrix(irb0,true, 0, icst,false,0, mtmp,   -dtmp1, true);
    ls.AddMatrix(icst,false,0, irb0,true, 0, mtmp,   -dtmp1, true);
    ls.AddMatrix(irb1,true, 0, icst,false,0, mtmp,    dtmp1, true);
    ls.AddMatrix(icst,false,0, irb1,true, 0, mtmp,    dtmp1, true);
    
    ls.AddMatrix(irb0,true, 3, irb0,true, 3, wXwRtL0, dtmp1, true);
    ls.AddMatrix(irb1,true, 3, irb1,true, 3, wXwRtL1,-dtmp1, true);
    
    ls.AddMatrix(irb0,true, 3, icst,false,0,  RotwX0,-dtmp1, false);
    ls.AddMatrix(icst,false,0, irb0,true, 3,  RotwX0,-dtmp1, true );
    ls.AddMatrix(irb1,true, 3, icst,false,0,  RotwX1, dtmp1, false);
    ls.AddMatrix(icst,false,0, irb1,true, 3,  RotwX1, dtmp1, true );
    
    ls.AddMatrix_Vector(irb0,true, 3,  icst,false,3,  Rtwma[0][0], dtmp1, true );
    ls.AddMatrix_Vector(icst,false,3,  irb0,true, 3,  Rtwma[0][0], dtmp1, false);
    ls.AddMatrix_Vector(irb0,true, 3,  icst,false,4,  Rtwma[0][1], dtmp1, true );
    ls.AddMatrix_Vector(icst,false,4,  irb0,true, 3,  Rtwma[0][1], dtmp1, false);
    ls.AddMatrix_Vector(irb1,true, 3,  icst,false,3,  Rtwma[1][0],-dtmp1, true );
    ls.AddMatrix_Vector(icst,false,3,  irb1,true, 3,  Rtwma[1][0],-dtmp1, false);
    ls.AddMatrix_Vector(irb1,true, 3,  icst,false,4,  Rtwma[1][1],-dtmp1, true );
    ls.AddMatrix_Vector(icst,false,4,  irb1,true, 3,  Rtwma[1][1],-dtmp1, false);
	}
}






////////////////////////////////////////////////////////////////
/*
void Rigid::CJoint_HingeRange::SetAxis(double ax, double ay, double az){
  axis = Com::CVector3D(ax,ay,az);
  if( axis.Length() < 1.0e-20 ){ assert(0); return; }
  axis.Normalize();
  GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

 void Rigid::CJoint_HingeRange::Draw(const std::vector<CRigidBody3D>& aRB) const
 {
 const unsigned int imode = 2;
 if( imode == 0 ){
 }
 else if( imode == 1 ){
 const unsigned int irb0 = aIndRB[0];
 const CRigidBody3D& rb0 = aRB[irb0];
 const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
 ::glColor3d(1,1,0);
 ::glBegin(GL_POINTS);
 ::glVertex3d(vec_j.x, vec_j.y, vec_j.z );
 ::glEnd();
 }
 else if( imode == 2 ){
 const unsigned int irb0 = aIndRB[0];
 const unsigned int irb1 = aIndRB[1];
 const CRigidBody3D& rb0 = aRB[irb0];
 const CRigidBody3D& rb1 = aRB[irb1];
 const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
 const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.GetDispCG();
 const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.GetDispCG();
 
 ::glPushMatrix();
 ::glTranslated( vec_j.x, vec_j.y, vec_j.z );
 ::glColor3d(1,1,0);
 ::glutSolidSphere(0.1,10,5);
 ::glPopMatrix();
 ::glColor3d(1,1,1);
 ::glLineWidth(2);
 ::glBegin(GL_LINES);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
 ::glVertex3d(vec_cg0.x,vec_cg0.y,vec_cg0.z);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 ::glVertex3d(vec_cg1.x,vec_cg1.y,vec_cg1.z);
 ::glEnd();
 
 const Com::CMatrix3& mrot = rb0.GetRotMatrix();
 const Com::CVector3D& lcb0 = mrot.MatVec(this->loc_coord[0]);
 const Com::CVector3D& lcb1 = mrot.MatVec(this->loc_coord[1]);
 unsigned int ndiv_t = 32;
 unsigned int ndiv0 = ndiv_t*(max_t-min_t)/360.0 + 1;
 const double dtheta = 2*3.1416/ndiv0*(max_t-min_t)/360.0;
 const double radius = 1;
 ::glColor3d(0,1,1);
 ::glBegin(GL_TRIANGLE_FAN);
 ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
 for(unsigned int idiv=0;idiv<ndiv0+1;idiv++){
 const Com::CVector3D v0 = vec_j 
 + sin(idiv*dtheta - max_t*3.1416/180.0 )*lcb0*radius 
 + cos(idiv*dtheta - max_t*3.1416/180.0 )*lcb1*radius;
 ::glVertex3d(v0.x,  v0.y,  v0.z);
 }
 ::glEnd();
 }
 }
 */
/*
void Rigid::CJoint_HingeRange::UpdateLambda
(const double* upd, double scale)
{
  lambda[0] += tmp1*upd[0];
  lambda[1] += tmp1*upd[1];
  lambda[2] += tmp1*upd[2];
  lambda[3] += tmp1*upd[3];
  lambda[4] += tmp1*upd[4];
  //    lambda[5] += tmp1*upd[5];
  const double tmp2 = dt*newmark_gamma;
  lambda[5] += tmp2*upd[5];
}

void Rigid::CJoint_HingeRange::AddLinearSystem_NewmarkBetaAPrime
(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
 const double dt, const double newmark_gamma, const double newmark_beta,
 const std::vector<CRigidBody3D>& aRB, bool is_initial ) const
{
  const unsigned int irb0 = aIndRB[0];
  assert( irb0 < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb0 = aRB[irb0];
  
  const unsigned int irb1 = aIndRB[1];
  assert( irb1 < ls.GetSizeRigidBody() );
  const CRigidBody3D& rb1 = aRB[irb1];
  
  const Com::CMatrix3& mrot0 = rb0.GetRotMatrix();
  const Com::CMatrix3& mrot1 = rb1.GetRotMatrix();
  const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
  const Com::CVector3D Xdistfix0 = rb0.GetIniPosCG() - ini_pos_joint;
  const Com::CVector3D Xdistfix1 = rb1.GetIniPosCG() - ini_pos_joint;
  
  ////////////////////////////////////////////////
  // 拘束力残差（物体０）
	// 並進の残差
  ls.AddResidual(irb0,true,0,  vlambda,1);
	{   // 変位拘束力から来る回転の残差（物体０）
    const Com::CVector3D& tmp_vec = mrot0.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix0);
    ls.AddResidual(irb0,true,3,  wX.MatVec(tmp_vec),-1);
	}
	{   // 変位拘束条件の変分（物体０）
    Com::CVector3D rot_pos_cg = mrot0.MatVec(Xdistfix0);
    ls.AddResidual(icst,false,0,   rb0.GetIniPosCG()+rb0.GetDispCG()-rot_pos_cg,1 );
	}
  Com::CMatrix3 RotwX0 = mrot0.MatMat( Com::CMatrix3(Xdistfix0) );
  Com::CMatrix3 wXwRtL0;
	{
    const Com::CVector3D& RtL = mrot0.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix0);
    wXwRtL0 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  ////////////////////////////////////////////////
  // 拘束力残差（物体１）
	// 並進の残差
  ls.AddResidual(irb1,true,0,  vlambda,-1);
	{   // 拘束力から来る回転の残差（物体１）
    const Com::CVector3D& tmp_vec = mrot1.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix1);
    ls.AddResidual(irb1,true,3,  wX.MatVec(tmp_vec),1 );
	}
	{   // 拘束条件の変分（物体１）
    const Com::CVector3D& rot_pos_cg = mrot1.MatVec(Xdistfix1);
    ls.AddResidual(icst,false,0,  rb1.GetIniPosCG()+rb1.GetDispCG()-rot_pos_cg,-1);
	}
  Com::CMatrix3 RotwX1 = mrot1.MatMat( Com::CMatrix3(Xdistfix1) );
  Com::CMatrix3 wXwRtL1;
	{
    const Com::CVector3D& RtL = mrot1.MatVecTrans(vlambda);
    const Com::CMatrix3 wX(Xdistfix1);
    wXwRtL1 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
  
  ////////////////////////////////
  
  Com::CVector3D Rtwma[2][2];
  {   // 回転拘束から来る物体回転の残差
    const Com::CVector3D& R1a = mrot1.MatVec(axis);
    const Com::CVector3D& R0m0 = mrot0.MatVec(loc_coord[0]);
    const Com::CVector3D& R0m1 = mrot0.MatVec(loc_coord[1]);
    const Com::CMatrix3 wR0m0(R0m0);
    const Com::CMatrix3 wR0m1(R0m1);
    const Com::CVector3D& wR0m0R1a = wR0m0.MatVec(R1a);
    const Com::CVector3D& wR0m1R1a = wR0m1.MatVec(R1a);
    Rtwma[0][0] = mrot0.MatVecTrans(wR0m0R1a);
    Rtwma[0][1] = mrot0.MatVecTrans(wR0m1R1a);
    Rtwma[1][0] = mrot1.MatVecTrans(wR0m0R1a);
    Rtwma[1][1] = mrot1.MatVecTrans(wR0m1R1a);
    ls.AddResidual(irb0,true,3, lambda[3]*Rtwma[0][0] + lambda[4]*Rtwma[0][1],-1 );
    ls.AddResidual(irb1,true,3, lambda[3]*Rtwma[1][0] + lambda[4]*Rtwma[1][1], 1 );
    double res[6] = {0,0,0,0,0,0};
    res[3] = Com::Dot(R1a,R0m0);
    res[4] = Com::Dot(R1a,R0m1);
    ls.SubResidual(icst,false,res);
  }
  {   // 角度についてのKKT条件を入れる
    bool flg = false;
    double v[2];
    double f_value;
    const double delta_f = 1.0e-8;
    {
      const double mint = -max_t;
      const double maxt = -min_t;
      double R0m0R1m0 = Com::Dot(mrot0.MatVec(loc_coord[0]), mrot1.MatVec(loc_coord[0]));
      double R0m0R1m1 = Com::Dot(mrot0.MatVec(loc_coord[0]), mrot1.MatVec(loc_coord[1]));
      const double PI = 3.14159265358979323846;
      double u_min[2] = { cos(PI*mint/180.0), sin(PI*mint/180.0) };
      double u_max[2] = { cos(PI*maxt/180.0), sin(PI*maxt/180.0) };
      double v_min[2] = { -u_min[1],  u_min[0] };
      double v_max[2] = {  u_max[1], -u_max[0] };
      double dmin = v_min[0]*R0m0R1m0 + v_min[1]*R0m0R1m1;
      double dmax = v_max[0]*R0m0R1m0 + v_max[1]*R0m0R1m1;
      //            std::cout << dmin << " " << dmax << " " << lambda[5] << std::endl;
      if( dmin < -delta_f || dmax < -delta_f ) flg = true;
      const double emin = u_min[0]*R0m0R1m0 + u_min[1]*R0m0R1m1;
      const double emax = u_max[0]*R0m0R1m0 + u_max[1]*R0m0R1m1;
      if( emin > emax ){ 
        v[0] = v_min[0]; v[1] = v_min[1];
        f_value = dmin;
      }
      else{                              
        v[0] = v_max[0]; v[1] = v_max[1];
        f_value = dmax;
      }
    }
    const double dtmp2 = dt*newmark_gamma;
    if( flg == true && is_initial )
    {
      //            std::cout << "consider contact" << std::endl;
      Com::CVector3D m = v[0]*loc_coord[0] + v[1]*loc_coord[1];
      ////////////////
      Com::CVector3D R0tR1m = mrot0.MatVecTrans( mrot1.MatVec( m ) );
      Com::CMatrix3 wm0( loc_coord[0] );
      Com::CVector3D wm0R0tR1m = wm0.MatVec( R0tR1m );
      Com::CVector3D res_t0 = lambda[5]*wm0R0tR1m;
      ls.AddResidual(irb0, true, 3,   res_t0,-1 );
      ls.AddMatrix_Vector(irb0,true, 3, icst,false,5, wm0R0tR1m,dtmp2, true );
      ls.AddMatrix_Vector(icst,false,5, irb0,true, 3, wm0R0tR1m,dtmp2, false);
      ////////////////
      Com::CVector3D R1tR0m0 = mrot1.MatVecTrans( mrot0.MatVec( loc_coord[0] ) );
      Com::CMatrix3 wm( m );
      Com::CVector3D wmR1tR0m0 = wm.MatVec( R1tR0m0 );
      Com::CVector3D res_t1 = lambda[5]*wmR1tR0m0;
      ls.AddResidual(irb1, true, 3,   res_t1,-1 );
      ls.AddMatrix_Vector(irb1,true, 3, icst,false,5, wmR1tR0m0,dtmp2, true );
      ls.AddMatrix_Vector(icst,false,5, irb1,true, 3, wmR1tR0m0,dtmp2, false);
      ////////////////
      const double df = Com::Dot(wm0R0tR1m,rb0.GetOmega()) + Com::Dot(wmR1tR0m0,rb1.GetOmega());
      ls.AddResidual(icst,false,5,1, &df,-1 );
    }
    else if( !is_initial && lambda[5]<-1.0e-5 )
    {
      //            std::cout << "after impact" << std::endl;
      Com::CVector3D m = v[0]*loc_coord[0] + v[1]*loc_coord[1];
      ////////////////
      Com::CVector3D R0tR1m = mrot0.MatVecTrans( mrot1.MatVec( m ) );
      Com::CMatrix3 wm0( loc_coord[0] );
      Com::CVector3D wm0R0tR1m = wm0.MatVec( R0tR1m );
      Com::CVector3D res_t0 = lambda[5]*wm0R0tR1m;
      ls.AddResidual(irb0, true, 3,   res_t0,-1 );
      ////////////////
      Com::CVector3D R1tR0m0 = mrot1.MatVecTrans( mrot0.MatVec( loc_coord[0] ) );
      Com::CMatrix3 wm( m );
      Com::CVector3D wmR1tR0m0 = wm.MatVec( R1tR0m0 );
      Com::CVector3D res_t1 = lambda[5]*wmR1tR0m0;
      ls.AddResidual(irb1, true, 3,   res_t1,-1 );
      ////////////////
      double val0 = 0.0;
      ls.AddResidual(icst,false,5,1, &val0,-1 );
      double val1 = 1.0;
      ls.AddMatrix(icst,false,5,1,  icst,false,5,1, &val1, dtmp2);
    }
    else{
      //            std::cout << "not contact" << std::endl;
      const double val0 = 1;
      ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp2);
      const double val1 = lambda[5];
      ls.AddResidual(icst,false,5,1, &val1,-1 );
    }
  }
  ////////////////////////////////////////////////////////////////
  
	{   // 剛性行列を作る
		const double dtmp1 = dt*dt*newmark_beta;
    Com::CMatrix3 mtmp;
    mtmp.SetIdentity();
    ls.AddMatrix(irb0,true, 0, icst,false,0, mtmp,   -dtmp1, true);
    ls.AddMatrix(icst,false,0, irb0,true, 0, mtmp,   -dtmp1, true);
    ls.AddMatrix(irb1,true, 0, icst,false,0, mtmp,    dtmp1, true);
    ls.AddMatrix(icst,false,0, irb1,true, 0, mtmp,    dtmp1, true);
    
    ls.AddMatrix(irb0,true, 3, irb0,true, 3, wXwRtL0, dtmp1, true);
    ls.AddMatrix(irb1,true, 3, irb1,true, 3, wXwRtL1,-dtmp1, true);
    
    ls.AddMatrix(irb0,true, 3, icst,false,0,  RotwX0,-dtmp1, false);
    ls.AddMatrix(icst,false,0, irb0,true, 3,  RotwX0,-dtmp1, true );
    ls.AddMatrix(irb1,true, 3, icst,false,0,  RotwX1, dtmp1, false);
    ls.AddMatrix(icst,false,0, irb1,true, 3,  RotwX1, dtmp1, true );
    
    ls.AddMatrix_Vector(irb0,true, 3,  icst,false,3,  Rtwma[0][0], dtmp1, true );
    ls.AddMatrix_Vector(icst,false,3,  irb0,true, 3,  Rtwma[0][0], dtmp1, false);
    ls.AddMatrix_Vector(irb0,true, 3,  icst,false,4,  Rtwma[0][1], dtmp1, true );
    ls.AddMatrix_Vector(icst,false,4,  irb0,true, 3,  Rtwma[0][1], dtmp1, false);
    ls.AddMatrix_Vector(irb1,true, 3,  icst,false,3,  Rtwma[1][0],-dtmp1, true );
    ls.AddMatrix_Vector(icst,false,3,  irb1,true, 3,  Rtwma[1][0],-dtmp1, false);
    ls.AddMatrix_Vector(irb1,true, 3,  icst,false,4,  Rtwma[1][1],-dtmp1, true );
    ls.AddMatrix_Vector(icst,false,4,  irb1,true, 3,  Rtwma[1][1],-dtmp1, false);
	}
}
*/


