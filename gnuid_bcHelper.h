// Gnuid, An Incompressible Navier-Stokes Solver for Hemodynamics 
// Copyright (C) 2010 Lorenzo Alessio Botti

/* This Incompressible Navier-Stokes Solver is free software; */
/* you can redistribute it and/or modify it under the terms of the */
/* GNU Lesser General Public  License as published by the Free Software Foundation */
/* either version 2.1 of the License, or (at your option) any later version. */

/* This software is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this software; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#ifndef __gnuid_bc_helper__
#define __gnuid_bc_helper__

#include <iostream>
#include "math.h"
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"

using namespace libMesh;

class GnuidBCHelper 
{
  public:
  GnuidBCHelper() {} 
  ~GnuidBCHelper() {}
  
  void init_dirichletprofile_bc(EquationSystems& es, const int& lid);
  void compute_dirichletprofile_bc(const Point& point, const Real& t, RealVectorValue& U_bc);
  static void bc_type_and_lid(const EquationSystems& es,  const unsigned int boundary_id, std::string& , int& lid);
  
  private:
  Real _t_period;
  Real _u_mean;
  Real _alpha;
  Real _vel_bc_Radius;
  std::string _vel_bc_Type;
  VectorValue<Real> _vel_bc_Normal;
  Point _vel_bc_Center;
  Real _scaling;
  std::vector<Complex> _f_modes;
};

#endif
