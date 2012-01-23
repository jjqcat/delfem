/*
 *  vector2d.cpp
 *  dfm_core
 *
 *  Created by Nobuyuki Umetani on 9/19/11.
 *  Copyright 2011 The University of Tokyo. All rights reserved.
 *
 */

#include "delfem/vector2d.h"

//using namespace Com;

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
double Com::FindNearestPointParameter_Line_Point
(const Com::CVector2D& po_c,
 const Com::CVector2D& po_s, const Com::CVector2D& po_e)
{
  const Com::CVector2D& es = po_e-po_s;
  const Com::CVector2D& sc = po_s-po_c;
  const double a = Com::SquareLength(es);
  const double b = Com::Dot(es,sc);
  return - b/a;
}

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
double Com::GetDist_LineSeg_Point
(const Com::CVector2D& po_c,
 const Com::CVector2D& po_s, const Com::CVector2D& po_e)
{
  double t = FindNearestPointParameter_Line_Point(po_c, po_s, po_e);
  if( t < 0 ){ return Distance(po_s,po_c); }
  if( t > 1 ){ return Distance(po_e,po_c); }
  Com::CVector2D p = po_s + t*(po_e-po_s);
  return Distance(p,po_c);
}

bool Com::IsCross_LineSeg_LineSeg(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,
                                  const Com::CVector2D& po_s1, const Com::CVector2D& po_e1 )
{
  {
    const double min0x = ( po_s0.x < po_e0.x ) ? po_s0.x : po_e0.x;
    const double max0x = ( po_s0.x > po_e0.x ) ? po_s0.x : po_e0.x;
    const double max1x = ( po_s1.x > po_e1.x ) ? po_s1.x : po_e1.x;
    const double min1x = ( po_s1.x < po_e1.x ) ? po_s1.x : po_e1.x;
    const double min0y = ( po_s0.y < po_e0.y ) ? po_s0.y : po_e0.y;
    const double max0y = ( po_s0.y > po_e0.y ) ? po_s0.y : po_e0.y;
    const double max1y = ( po_s1.y > po_e1.y ) ? po_s1.y : po_e1.y;
    const double min1y = ( po_s1.y < po_e1.y ) ? po_s1.y : po_e1.y;
    const double len = ((max0x-min0x)+(max0y-min0y)+(max1x-min1x)+(max1y-min1y))*0.0001;
    //		std::cout << len << std::endl;
    if( max1x+len < min0x ) return false;
    if( max0x+len < min1x ) return false;
    if( max1y+len < min0y ) return false;
    if( max0y+len < min1y ) return false;
  }
  const double area1 = Com::TriArea(po_s0,po_e0,po_s1);
  const double area2 = Com::TriArea(po_s0,po_e0,po_e1);
  const double area3 = Com::TriArea(po_s1,po_e1,po_s0);
  const double area4 = Com::TriArea(po_s1,po_e1,po_e0);  
  //	std::cout << area1 << " " << area2 << " " << area3 << " " << area4 << std::endl;
  const double a12 = area1*area2; if( a12 > 0 ) return false;
  const double a34 = area3*area4; if( a34 > 0 ) return false;
  return true;
}

double Com::GetDist_LineSeg_LineSeg(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,
                                    const Com::CVector2D& po_s1, const Com::CVector2D& po_e1)
{
  if( Com::IsCross_LineSeg_LineSeg(po_s0,po_e0, po_s1,po_e1) ) return -1;
  const double ds1 = Com::GetDist_LineSeg_Point(po_s0,po_s1,po_e1);
  const double de1 = Com::GetDist_LineSeg_Point(po_e0,po_s1,po_e1);
  const double ds0 = Com::GetDist_LineSeg_Point(po_s1,po_s0,po_e0);
  const double de0 = Com::GetDist_LineSeg_Point(po_e1,po_s0,po_e0);
  double min_dist = ds1;
  min_dist = ( de1 < min_dist ) ? de1 : min_dist;
  min_dist = ( ds0 < min_dist ) ? ds0 : min_dist;
  min_dist = ( de0 < min_dist ) ? de0 : min_dist;    
  return min_dist;  
}



//! square root of circumradius
double Com::SquareCircumradius(const Com::CVector2D& p0, 
                               const Com::CVector2D& p1, 
                               const Com::CVector2D& p2 )
{
	const double area = TriArea(p0,p1,p2);
  
	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);
  
	return dtmp0*dtmp1*dtmp2/(16.0*area*area);
}

//! center of the circumcircle
bool Com::CenterCircumcircle(const Com::CVector2D& p0, 
                             const Com::CVector2D& p1, 
                             const Com::CVector2D& p2, 
                             Com::CVector2D& center)
{
  const double area = TriArea(p0,p1,p2);
  if( fabs(area) < 1.0e-10 ){ return false; }
  const double tmp_val = 1.0/(area*area*16.0);
  
  const double dtmp0 = SquareLength(p1,p2);
  const double dtmp1 = SquareLength(p0,p2);
  const double dtmp2 = SquareLength(p0,p1);
  
  const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
  const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
  const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
  center.x = etmp0*p0.x + etmp1*p1.x + etmp2*p2.x;
  center.y = etmp0*p0.y + etmp1*p1.y + etmp2*p2.y;
  return true;
}


////////////////////////////////

//! check if Delaunay condition satisfied
// 0 : p3 is inside circum circle on the p0,p1,p2
// 1 :       on         
// 2 :       outsdie 
int Com::DetDelaunay(const Com::CVector2D& p0, 
                     const Com::CVector2D& p1, 
                     const Com::CVector2D& p2, 
                     const Com::CVector2D& p3)
{
	const double area = TriArea(p0,p1,p2);
	if( fabs(area) < 1.0e-10 ){
		return 3;
	}
	const double tmp_val = 1.0/(area*area*16.0);
  
	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);
  
	const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
	const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
	const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
	const CVector2D out_center(
                             etmp0*p0.x + etmp1*p1.x + etmp2*p2.x,
                             etmp0*p0.y + etmp1*p1.y + etmp2*p2.y );
  
	const double qradius = SquareLength(out_center,p0);
	const double qdistance = SquareLength(out_center,p3);
  
  //	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
  //	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );
  
	const double tol = 1.0e-20;
	if( qdistance > qradius*(1.0+tol) ){ return 2; }	// outside the circumcircle
	else{
		if( qdistance < qradius*(1.0-tol) ){ return 0; }	// inside the circumcircle
		else{ return 1;	}	// on the circumcircle
	}
	return 0;
}

////////////////////////////////////////////////

