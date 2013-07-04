
#ifndef __VTKWriter_h__
#define __VTKWriter_h__

// C++ includes
#include <map>
#include <fstream>

// Local includes
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"
#include "libmesh/parallel.h"
#include "libmesh/fe.h"

#include "libmesh/fe_interface.h"

using namespace libMesh;

class VTKWriter 
{
public:
  VTKWriter () 
  {
    ;
  }

  ~VTKWriter () {}

  void write_ascii_discontinuous(const std::string& dir, const unsigned int& step, const MeshBase& mesh, const std::vector<double>& soln);
  
  void write_ascii_continuous(const std::string& dir, const unsigned int& step, const MeshBase& mesh, const std::vector<double>& soln);

  void build_discontinuous_solution_vector(MeshBase& mesh, libMesh::EquationSystems& systems, std::vector<Real>& solution);
  
  void build_continuous_solution_vector(MeshBase& mesh, libMesh::EquationSystems& systems, std::vector<Real>& solution);
private:

  void cell_connectivity(const Elem* elem, std::vector<unsigned int>& vtk_cell_connectivity);
  
  unsigned int cell_type(const Elem* elem);

  unsigned int cell_offset(const Elem* elem);
};

#endif 
