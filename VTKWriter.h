
#ifndef __VTKWriter_h__
#define __VTKWriter_h__

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
