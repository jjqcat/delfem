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

#if !defined(CAD_SVG_H)
#define CAD_SVG_H

#include <string>

#include "delfem/cad_obj2d.h"

namespace Cad
{

std::vector<unsigned int> ReadSVG_AddLoopCad(const std::string& fname, CCadObj2D& cad_2d, double scale=1);
  
}

#endif

