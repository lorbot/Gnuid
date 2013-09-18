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

#include <iostream>
#include <algorithm>
#include <math.h>

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/perf_log.h"
#include "libmesh/equation_systems.h"

#include "ns_dg_solver.h"
#include "gnuid_inputReader.h"

using namespace libMesh;
int main (int argc, char** argv)
{
  LibMeshInit init (argc, argv);
  {    
    PerfLog perf_log("gnuid");

    std::string working_directory = ".";
    if (argc >= 2)
    {
      working_directory = argv[1];
    }
 
    const unsigned int dim = 3;     
    Mesh mesh (dim);
    GnuidInputReader::read_mesh(mesh, working_directory);
    
    EquationSystems es (mesh);
    GnuidInputReader::read_es(es, working_directory);

    NS_DG_Solver ns_dg_Solver(es);
    ns_dg_Solver.init();
    ns_dg_Solver.working_directory = working_directory;

    perf_log.start_event("gnuid");
    ns_dg_Solver.solve();
    perf_log.stop_event("gnuid");
    
    std::cout<<es.parameters<<std::endl;
  }
  return EXIT_SUCCESS;
}
