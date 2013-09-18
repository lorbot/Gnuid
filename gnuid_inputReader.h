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

#ifndef __gnuid_input_reader__
#define __gnuid_input_reader__

#include <iostream>
#include "math.h"
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"

using namespace libMesh;

struct GnuidInputReader 
{
  GnuidInputReader() {} 
  ~GnuidInputReader() {}
  
  static void read_mesh(Mesh& mesh, const std::string& workingDir);
  static void read_es(EquationSystems& es, const std::string& workingDir);
};

#endif
