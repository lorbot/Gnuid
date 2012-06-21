
#ifndef __VTKDiscontinuousWriter_h__
#define __VTKDiscontinuousWriter_h__

// C++ includes
#include <map>
#include <fstream>

// Local includes
#include "mesh_base.h"
#include "elem.h"
#include "mesh.h"
#include "dof_map.h"
#include "equation_systems.h"
#include "system.h"
#include "parallel.h"
#include "fe.h"

#include "fe_interface.h"


class VTKDiscontinuousWriter 
{
public:
  VTKDiscontinuousWriter () 
  {
    ;
  }

  ~VTKDiscontinuousWriter () {}

  void write_ascii(const std::string& dir, const unsigned int& step, const MeshBase& mesh, const std::vector<double>& soln);

  void build_discontinuous_solution_vector(MeshBase& mesh, libMesh::EquationSystems& systems, std::vector<Real>& solution);
private:

  void cell_connectivity(const Elem* elem, std::vector<unsigned int>& vtk_cell_connectivity);
  
  unsigned int cell_type(const Elem* elem);

  unsigned int cell_offset(const Elem* elem);
};

#endif 
