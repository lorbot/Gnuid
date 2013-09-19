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

#ifndef __ns_dg_solver__
#define __ns_dg_solver__

#include <iostream>
#include "math.h"
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"
#include "libmesh/vector_value.h"
#include "libmesh/equation_systems.h"

using namespace libMesh;

class GnuidSolver 
{
public:

  GnuidSolver(EquationSystems& es): _system(es) {} 
  ~GnuidSolver() {}
  
  void init();
  void solve();
  EquationSystems & system () const { return _system; }

  std::string working_directory;

protected:
  static void assemble_adv_diff(EquationSystems& es, const std::string& system_name);
  static void assemble_p_proj(EquationSystems& es, const std::string& system_name);

  EquationSystems& _system;
};

#endif
