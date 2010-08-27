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
@brief the basical definition refered from whole system including mesh and fem
@author Nobuyuki Umetani
*/

#if !defined(CAD_COM_H)
#define CAD_COM_H

#ifdef __VISUALC__
	#pragma warning( disable : 4786 )
#endif

namespace Cad{


/*!
@ingroup CAD
@brief the type of component 
*/
enum CAD_ELEM_TYPE{
  NOT_SET,  //!< not setted(for the error handling)
  VERTEX,   //!< vertex
  EDGE,     //!< edge
  LOOP,     //!< loop
	SOLID,		//!< solid
};


}

#endif
