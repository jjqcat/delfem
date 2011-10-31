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

#if !defined(VECTOR_3D_H)
#define VECTOR_3D_H

#include <vector>
#include <cassert>
#include <math.h>
#include <iostream>

#define NEARLY_ZERO 1.e-16


namespace Com
{

// rule about naming, the method starts "Set" change it self (not const)
class CVector3D;

//! @{
CVector3D operator+(const CVector3D&, const CVector3D&);
CVector3D operator-(const CVector3D&, const CVector3D&);
CVector3D operator*(double, const CVector3D&);
CVector3D operator*(const CVector3D&, double);
//! @}

//! 3 dimentional vector class
class CVector3D  
{
public:
	CVector3D(double vx, double vy, double vz) : x(vx), y(vy), z(vz){}
	CVector3D(): x(0.0), y(0.0), z(0.0){}
	CVector3D(const CVector3D& rhs){
		x = rhs.x; y = rhs.y; z = rhs.z;
	}
	virtual ~CVector3D(){}

	void SetVector(double vx, double vy, double vz){ x = vx; y = vy; z = vz; }

	inline const CVector3D operator-() const{ return -1.0*(*this); }
	inline const CVector3D operator+() const{ return *this; }  
	inline CVector3D& operator=(const CVector3D& rhs){
		if( this != &rhs ){ x = rhs.x; y = rhs.y; z = rhs.z; }
		return *this;
	}
	inline CVector3D& operator+=(const CVector3D& rhs){
		x += rhs.x; y += rhs.y; z += rhs.z;
		return *this;
	}
	inline CVector3D& operator-=(const CVector3D& rhs){
		x -= rhs.x; y -= rhs.y; z -= rhs.z;
		return *this;
	}
	inline CVector3D& operator*=(double d){
		x *= d; y *= d; z *= d;
		return *this;
	}
	inline CVector3D& operator/=(double d){
		if( fabs(d) < NEARLY_ZERO ){ return *this; }
		x /= d; y /= d; z /= d;
		return *this;
	}
  inline double operator[](int i) const{
    if( i == 0 ) return x;
    if( i == 1 ) return y;
    if( i == 2 ) return z;
    return 0;
  }
  inline double& operator[](int i){
    if( i == 0 ) return x;
    if( i == 1 ) return y;
    if( i == 2 ) return z;
    assert(0);
    return x;
  }  
	inline CVector3D operator+(){ return *this; }
	inline CVector3D operator-(){ return CVector3D(-x,-y,-z); }

	friend bool operator==(const CVector3D&, const CVector3D&);
	friend bool operator!=(const CVector3D&, const CVector3D&);

	friend Com::CVector3D Cross(const Com::CVector3D&, const Com::CVector3D&);
	friend double Dot(const Com::CVector3D&, const Com::CVector3D&);

	inline double Length()  const{ return sqrt( x*x+y*y+z*z ); }
	inline double DLength() const{ return x*x+y*y+z*z; }
	void SetNormalizedVector();
	void SetZero();
public:
	double x;	//!< x axis coordinate
	double y;	//!< y axis coordinate
	double z;	//!< z axis coordinate 
};

//! @{
	
/*! 
@brief inner product
*/
inline double Dot(const Com::CVector3D &arg1, const Com::CVector3D &arg2)
{
	return arg1.x*arg2.x + arg1.y*arg2.y + arg1.z*arg2.z;
}

/*! 
@brief cross product
*/
inline Com::CVector3D Cross(const Com::CVector3D& arg1, const Com::CVector3D& arg2)
{
	CVector3D temp;
	temp.x = arg1.y*arg2.z - arg1.z*arg2.y;
	temp.y = arg1.z*arg2.x - arg1.x*arg2.z;
	temp.z = arg1.x*arg2.y - arg1.y*arg2.x;
	return temp;
}

//! add
inline CVector3D operator + (const CVector3D& lhs, const CVector3D& rhs){
	CVector3D temp = lhs;
	temp += rhs;
	return temp;
}

//! subtract
inline CVector3D operator - (const CVector3D& lhs, const CVector3D& rhs){
	CVector3D temp = lhs;
	temp -= rhs;
	return temp;
}

//! scale
inline CVector3D operator * (double d, const CVector3D& rhs){
	CVector3D temp = rhs;
	temp *= d;
	return temp;
}
  
//! scale
inline CVector3D operator * (const CVector3D& vec, double d){
  CVector3D temp = vec;
  temp *= d;
  return temp;
}  

//! divide by real number
inline CVector3D operator / (const CVector3D& vec, double d){
	CVector3D temp = vec;
	temp /= d;
	return temp;
}
  
inline std::ostream &operator<<(std::ostream &output, const CVector3D& v)
{
  output.setf(std::ios::scientific);
  output << v.x << " " << v.y << " " << v.z; 
  return output;
}
    

inline double ScalarTripleProduct3D(const double a[], const double b[], const double c[]){
	return a[0]*(b[1]*c[2] - b[2]*c[1]) 
        +a[1]*(b[2]*c[0] - b[0]*c[2])
        +a[2]*(b[0]*c[1] - b[1]*c[0]);
}

inline double Dot3D(const double a[], const double b[]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
  
inline void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}
  
inline double Length3D(const double v[3]){
    return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}  

static inline double SquareDistance3D(const double p0[3], const double p1[3]){
  return (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
}
  
  
static inline double Distance3D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

static inline double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
  double x, y, z;
  x = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  y = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return 0.5*sqrt( x*x + y*y + z*z );
}

inline void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
    n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
    n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
    n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
    a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
    const double invlen = 0.5/a;
    n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}
  
inline void NormalTri3D(double n[3], const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
}
  
static inline double TetVolume3D(const double v1[3],
                                 const double v2[3], 
                                 const double v3[3], 
                                 const double v4[3] )
{
  return	
  (   ( v2[0] - v1[0] )*( ( v3[1] - v1[1] )*( v4[2] - v1[2] ) - ( v4[1] - v1[1] )*( v3[2] - v1[2] ) )		
   -	( v2[1] - v1[1] )*( ( v3[0] - v1[0] )*( v4[2] - v1[2] ) - ( v4[0] - v1[0] )*( v3[2] - v1[2] ) )		
   +	( v2[2] - v1[2] )*( ( v3[0] - v1[0] )*( v4[1] - v1[1] ) - ( v4[0] - v1[0] )*( v3[1] - v1[1] ) )
   ) * 0.16666666666666666666666666666667;		
}

  

static void GetVertical2Vector3D
(const double vec_n[3], 
 double vec_x[3], double vec_y[3])
{
  const double vec_s[3] = {0,1,0};
  Cross3D(vec_x,vec_s,vec_n);
  const double len = Length3D(vec_x);
  if( len < 1.0e-10 ){
    const double vec_t[3] = {1,0,0};
    Cross3D(vec_x,vec_t,vec_n);  // z????
    Cross3D(vec_y,vec_n,vec_x);  // x????
  }
  else{
    const double invlen = 1.0/len;
    vec_x[0] *= invlen;
    vec_x[1] *= invlen;
    vec_x[2] *= invlen;
    Cross3D(vec_y,vec_n,vec_x);
  }
}  
  
  

//! 3D bounding box class
class CBoundingBox3D
{
public:
	CBoundingBox3D(){
		x_min=0;	x_max=0;
		y_min=0;	y_max=0;
		z_min=0;	z_max=0;
		isnt_empty = false;
	}
	CBoundingBox3D(double x_min0,double x_max0,  double y_min0,double y_max0,  double z_min0,double z_max0)
		: x_min(x_min0),x_max(x_max0),  
		  y_min(y_min0),y_max(y_max0),  
		  z_min(z_min0),z_max(z_max0)
	{
		assert( x_min <= x_max );
		assert( y_min <= y_max );
		assert( z_min <= z_max );
		isnt_empty = true;
	}
	CBoundingBox3D( const CBoundingBox3D& bb )
		: x_min(bb.x_min),x_max(bb.x_max), 
		  y_min(bb.y_min),y_max(bb.y_max),  
		  z_min(bb.z_min),z_max(bb.z_max), 
		  isnt_empty(bb.isnt_empty){}

	CBoundingBox3D& operator+=(const CBoundingBox3D& bb)
	{
		if( !bb.isnt_empty ) return *this;
		if( !isnt_empty ){
			x_max = bb.x_max;	x_min = bb.x_min;
			y_max = bb.y_max;	y_min = bb.y_min;
			z_max = bb.z_max;	z_min = bb.z_min;
      this->isnt_empty = bb.isnt_empty;
			return *this;
		}
		x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
		x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
		y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
		y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
		z_max = ( z_max > bb.z_max ) ? z_max : bb.z_max;
		z_min = ( z_min < bb.z_min ) ? z_min : bb.z_min;
		return *this;
	}
  CBoundingBox3D& operator+=(const CVector3D& v)
	{
		if( !isnt_empty ){
			x_max = v.x;	x_min = v.x;
			y_max = v.y;	y_min = v.y;
			z_max = v.z;	z_min = v.z;
      this->isnt_empty = true;
			return *this;
		}
		x_max = ( x_max > v.x ) ? x_max : v.x;
		x_min = ( x_min < v.x ) ? x_min : v.x;
		y_max = ( y_max > v.y ) ? y_max : v.y;
		y_min = ( y_min < v.y ) ? y_min : v.y;
		z_max = ( z_max > v.z ) ? z_max : v.z;
		z_min = ( z_min < v.z ) ? z_min : v.z;
		return *this;
	}  
	bool IsInside(const CVector3D& vec)
	{
		if( !isnt_empty ) return false; // âΩÇ‡Ç»Ç¢èÍçáÇÕèÌÇ…ãU
		if(  vec.x >= x_min && vec.x <= x_max
			&& vec.y >= y_min && vec.y <= y_max 
			&& vec.z >= z_min && vec.z <= z_max ) return true;
		return false;
	}
	bool IsPossibilityIntersectSphere(const CVector3D& vec, const double radius ) const
	{
		if( !isnt_empty ) return false; // âΩÇ‡Ç»Ç¢èÍçáÇÕèÌÇ…ãU
		if(vec.x < x_min-radius || vec.x > x_max+radius ||
       vec.y < y_min-radius || vec.y > y_max+radius || 
       vec.z < z_min-radius || vec.z > z_max+radius ) return false;
		// ä±è¬ÇµÇ»Ç¢Ç‚Ç¬Ç™ä‹Ç‹ÇÍÇƒÇ¢ÇÈÇ™ÅAÉAÉoÉEÉgÇ»Ç‹Ç‹Ç≈ÇÊÇ¢ÅD
		return true;
	}
  bool AddPoint(const Com::CVector3D& vec, double eps){
    if( eps <= 0 ){ return false; }
    if( isnt_empty ){
      x_min = ( x_min < vec.x-eps ) ? x_min : vec.x-eps;
      y_min = ( y_min < vec.y-eps ) ? y_min : vec.y-eps;
      z_min = ( z_min < vec.z-eps ) ? z_min : vec.z-eps;
      x_max = ( x_max > vec.x+eps ) ? x_max : vec.x+eps;
      y_max = ( y_max > vec.y+eps ) ? y_max : vec.y+eps;
      z_max = ( z_max > vec.z+eps ) ? z_max : vec.z+eps;                    
    }
    else{
      isnt_empty = true;
      x_min = vec.x-eps;
      y_min = vec.y-eps;
      z_min = vec.z-eps;      
      x_max = vec.x+eps;      
      y_max = vec.y+eps;      
      z_max = vec.z+eps;            
    }
    return true;
  }
  void SetValueToArray(double bb[8]) const{
    bb[0] = x_min;  bb[2] = y_min;  bb[4] = z_min;
    bb[1] = x_max;  bb[3] = y_max;  bb[5] = z_max;    
  }
  void ProjectOnLine(double& min_r, double& max_r, 
                   const Com::CVector3D& org, const Com::CVector3D& dir) const
  {
    const double d[8] = {
      Dot(dir,Com::CVector3D(x_min,y_min,z_min)-org),
      Dot(dir,Com::CVector3D(x_max,y_min,z_min)-org),
      Dot(dir,Com::CVector3D(x_min,y_max,z_min)-org),
      Dot(dir,Com::CVector3D(x_max,y_max,z_min)-org),
      Dot(dir,Com::CVector3D(x_min,y_min,z_max)-org),
      Dot(dir,Com::CVector3D(x_max,y_min,z_max)-org),
      Dot(dir,Com::CVector3D(x_min,y_max,z_max)-org),
      Dot(dir,Com::CVector3D(x_max,y_max,z_max)-org) };
    min_r = max_r = d[0];
    for(unsigned int i=1;i<8;i++){
      min_r = ( d[i] < min_r ) ? d[i] : min_r;
      max_r = ( d[i] > max_r ) ? d[i] : max_r;
    }
  }
public:
	double x_min,x_max,  y_min,y_max,  z_min,z_max;
	bool isnt_empty;	//!< false if there is nothing inside
};



//! Hight of a tetrahedra
inline double Height(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3, const CVector3D& v4){
	// get normal vector
	double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
	double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
	double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);	

	// normalize normal vector
	const double dtmp1 = 1.0 / sqrt( dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z );
	dtmp_x *= dtmp1;
	dtmp_y *= dtmp1;
	dtmp_z *= dtmp1;

	return (v4.x-v1.x)*dtmp_x+(v4.y-v1.y)*dtmp_y+(v4.z-v1.z)*dtmp_z;
}


////////////////////////////////////////////////////////////////
  
static void GetVertical2Vector
(const CVector3D& vec_n, 
 CVector3D& vec_x, CVector3D& vec_y)
{
  vec_x = Cross(CVector3D(0,1,0),vec_n);
  const double len = vec_x.Length();
  if( len < 1.0e-10 ){
    vec_x = Cross(CVector3D(1,0,0),vec_n);  // z????
    vec_x.SetNormalizedVector();
    vec_y = Cross(vec_n,vec_x);  // x????
  }
  else{
    const double invlen = 1.0/len;
    vec_x *= invlen;
    vec_y = Cross(vec_n,vec_x);
  }
}    
  
  
  

//! Volume of a tetrahedra
inline double TetVolume(const CVector3D& v1,
						const CVector3D& v2, 
						const CVector3D& v3, 
						const CVector3D& v4 )
{
	return	
		(   ( v2.x - v1.x )*( ( v3.y - v1.y )*( v4.z - v1.z ) - ( v4.y - v1.y )*( v3.z - v1.z ) )		
		-	( v2.y - v1.y )*( ( v3.x - v1.x )*( v4.z - v1.z ) - ( v4.x - v1.x )*( v3.z - v1.z ) )		
		+	( v2.z - v1.z )*( ( v3.x - v1.x )*( v4.y - v1.y ) - ( v4.x - v1.x )*( v3.y - v1.y ) )
		) * 0.16666666666666666666666666666667;		
}

//! élñ ëÃÇÃëÃêœ
inline double TetVolume( int iv1, int iv2, int iv3, int iv4, 
						const std::vector<CVector3D>& node)
{
	return TetVolume(node[iv1],node[iv2],node[iv3],node[iv4]);
}

////////////////////////////////////////////////

//! äOê⁄ÉxÉNÉgÉã
inline void Cross( CVector3D& lhs, const CVector3D& v1, const CVector3D& v2 ){
	lhs.x = v1.y*v2.z - v2.y*v1.z;
	lhs.y = v1.z*v2.x - v2.z*v1.x;
	lhs.z = v1.x*v2.y - v2.x*v1.y;
}

//! ÇRéüå≥ÇRäpå`ÇÃñ êœ
inline double TriArea(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3)
{
	double x, y, z;
	x = ( v2.y - v1.y )*( v3.z - v1.z ) - ( v3.y - v1.y )*( v2.z - v1.z );
	y = ( v2.z - v1.z )*( v3.x - v1.x ) - ( v3.z - v1.z )*( v2.x - v1.x );
	z = ( v2.x - v1.x )*( v3.y - v1.y ) - ( v3.x - v1.x )*( v2.y - v1.y );
	return 0.5*sqrt( x*x + y*y + z*z );
}

//! ÇRéüå≥ÇRäpå`ÇÃñ êœ
inline double TriArea(const int iv1, const int iv2, const int iv3, 
                      const std::vector<CVector3D>& node )
{
	return TriArea(node[iv1],node[iv2],node[iv3]);
}

//! ÇRéüå≥ÇRäpå`ÇÃñ êœÇÃÇQèÊ
inline double SquareTriArea(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3)
{
	double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
	double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
	double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
	return (dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z)*0.25;
}

////////////////////////////////////////////////

//! í∑Ç≥ÇÃÇQèÊ
inline double SquareDistance(const CVector3D& ipo0, const CVector3D& ipo1)
{
	return	( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y ) + ( ipo1.z - ipo0.z )*( ipo1.z - ipo0.z );
}

//! í∑Ç≥ÇÃÇQèÊ
inline double SquareLength(const CVector3D& point)
{
	return	point.x*point.x + point.y*point.y + point.z*point.z;
}

////////////////////////////////////////////////

//! length of vector
inline double Length(const CVector3D& point)
{
	return	sqrt( point.x*point.x + point.y*point.y + point.z*point.z );
}

//! distance between two points
inline double Distance(const CVector3D& ipo0, const CVector3D& ipo1)
{
	return	sqrt( SquareDistance(ipo0,ipo1) );
}

////////////////////////////////////////////////

//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡í∑Ç¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
inline double SqareLongestEdgeLength(
		const CVector3D& ipo0,
		const CVector3D& ipo1,
		const CVector3D& ipo2,
		const CVector3D& ipo3 )
{
	double edge1, edge2;
	edge1 = SquareDistance( ipo0, ipo1 );
	edge2 = SquareDistance( ipo0, ipo2 );
	if( edge2 > edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo0, ipo3 );
	if( edge2 > edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo1, ipo2 );
	if( edge2 > edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo1, ipo3 );
	if( edge2 > edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo2, ipo3 );
	if( edge2 > edge1 ) edge1 = edge2;
	return edge1;
}

////////////////////////////////////////////////

//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡í∑Ç¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
inline double LongestEdgeLength(
		const CVector3D& ipo0,
		const CVector3D& ipo1,
		const CVector3D& ipo2,
		const CVector3D& ipo3 )
{
	return sqrt( SqareLongestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

////////////////////////////////////////////////

//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡íZÇ¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
inline double SqareShortestEdgeLength(const CVector3D& ipo0,
						  const CVector3D& ipo1,
						  const CVector3D& ipo2,
						  const CVector3D& ipo3 )
{
	double edge1, edge2;
	edge1 = SquareDistance( ipo0, ipo1 );
	edge2 = SquareDistance( ipo0, ipo2 );
	if( edge2 < edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo0, ipo3 );
	if( edge2 < edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo1, ipo2 );
	if( edge2 < edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo1, ipo3 );
	if( edge2 < edge1 ) edge1 = edge2;
	edge2 = SquareDistance( ipo2, ipo3 );
	if( edge2 < edge1 ) edge1 = edge2;
	return edge1;
}

////////////////////////////////////////////////


//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡íZÇ¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
inline double ShortestEdgeLength(
		const CVector3D& ipo0,
		const CVector3D& ipo1,
		const CVector3D& ipo2,
		const CVector3D& ipo3 )
{
	return sqrt( SqareShortestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

////////////////////////////////////////////////

//! ñ@ê¸ÉxÉNÉgÉã
inline void Normal(
		CVector3D& vnorm, 
		const CVector3D& v1, 
		const CVector3D& v2, 
		const CVector3D& v3)
{
	vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
	vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
	vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);	
}

////////////////////////////////////////////////

//! íPà ñ@ê¸ÉxÉNÉgÉã
inline void UnitNormal(
		CVector3D& vnorm, 
		const CVector3D& v1, 
		const CVector3D& v2, 
		const CVector3D& v3)
{
	vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
	vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
	vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);	
	const double dtmp1 = 1.0 / Length(vnorm);
	vnorm.x *= dtmp1;
	vnorm.y *= dtmp1;
	vnorm.z *= dtmp1;
}

////////////////////////////////////////////////

/*! 
äOê⁄ãÖÇÃîºåa
*/
inline double SquareCircumradius(
		const CVector3D& ipo0, 
		const CVector3D& ipo1, 
		const CVector3D& ipo2, 
		const CVector3D& ipo3)
{
	double base[3][3] = {
		{ ipo1.x-ipo0.x, ipo1.y-ipo0.y, ipo1.z-ipo0.z },
		{ ipo2.x-ipo0.x, ipo2.y-ipo0.y, ipo2.z-ipo0.z },
		{ ipo3.x-ipo0.x, ipo3.y-ipo0.y, ipo3.z-ipo0.z }
	};
	double s[6] = {
		base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
		base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
		base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
		base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
		base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
		base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
	};
	const double vol = TetVolume(ipo0,ipo1,ipo2,ipo3)*6.0;
	if( vol < 1.0e-20 ){ assert(0); }
	const double inv_det = 1.0 / (vol*vol);
	double t[6] = {
		(s[1]*s[2]-s[3]*s[3])*0.5*inv_det,
		(s[2]*s[0]-s[4]*s[4])*0.5*inv_det,
		(s[0]*s[1]-s[5]*s[5])*0.5*inv_det,
		(s[4]*s[5]-s[0]*s[3])*0.5*inv_det,
		(s[5]*s[3]-s[1]*s[4])*0.5*inv_det,
		(s[3]*s[4]-s[2]*s[5])*0.5*inv_det,
	};
	double u[3] = {
		t[0]*s[0]+t[5]*s[1]+t[4]*s[2],
		t[5]*s[0]+t[1]*s[1]+t[3]*s[2],
		t[4]*s[0]+t[3]*s[1]+t[2]*s[2],
	};
	return  0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
	/*
	const double square_radius = 0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
	CVector3D vec1;
	vec1.x = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.x;
	vec1.y = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.y;
	vec1.z = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.z;
	std::cout << square_radius << " ";
	std::cout << SquareLength(vec1,ipo0) << " ";
	std::cout << SquareLength(vec1,ipo1) << " ";
	std::cout << SquareLength(vec1,ipo2) << " ";
	std::cout << SquareLength(vec1,ipo3) << std::endl;;
	return square_radius;
	*/
}

////////////////////////////////////////////////

/*! 
curcumradius of a tetrahedra
*/
inline double Circumradius(const CVector3D& ipo0, 
						   const CVector3D& ipo1, 
						   const CVector3D& ipo2, 
						   const CVector3D& ipo3){
	return sqrt( SquareCircumradius(ipo0,ipo1,ipo2,ipo3) );
}


inline CVector3D RotateVector(const CVector3D& vec0, const CVector3D& rot )
{
	const double theta = rot.Length();
	if( theta < 1.0e-30 ){
		return vec0;
	}
	Com::CVector3D e0 = rot;
	e0.SetNormalizedVector();
	Com::CVector3D e2 = Com::Cross(e0,vec0);
	if( e2.Length() < 1.0e-30 ){
		return vec0;
	}
	e2.SetNormalizedVector();
	Com::CVector3D e1 = Com::Cross(e2,e0);
	assert( fabs( e1.Length() - 1 ) < 1.0e-10 );
//	assert( e2.x*vec_0.x + e2.y*vec_0.y + e2.z*vec_0.z < 1.0e-10 );
	const double dot00 = Dot(vec0,e0);
	const double dot01 = Dot(vec0,e1);
	const double cost = cos(theta);
	const double sint = sin(theta);
	CVector3D vec1;
	vec1.x = dot00*e0.x + dot01*cost*e1.x + dot01*sint*e2.x;
	vec1.y = dot00*e0.y + dot01*cost*e1.y + dot01*sint*e2.y;
	vec1.z = dot00*e0.z + dot01*cost*e1.z + dot01*sint*e2.z;
	return vec1;
}



}	// end namespace Com

#endif // !defined(VECTOR_3D_H)
