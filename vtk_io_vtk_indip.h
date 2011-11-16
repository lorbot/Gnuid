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

#ifndef __vtk_io_vtk_indip_h__
#define __vtk_io_vtk_indip_h__

// C++ includes
#include <map>

// Local includes
#include "libmesh_common.h"
#include "mesh_input.h"
#include "mesh_output.h"
#include "elem.h"

using namespace libMesh;

// Forward declarations
class libMesh::MeshBase;
class libMesh::MeshData;

class VTKIO_NOVTK : public MeshInput<MeshBase>,
	      public MeshOutput<MeshBase>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  VTKIO_NOVTK (MeshBase& mesh, MeshData* mesh_data=NULL);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  VTKIO_NOVTK (const MeshBase& mesh, MeshData* mesh_data=NULL);

 /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
//  virtual void write_nodal_data (const std::string&,
//             const std::vector<Number>&,
//             const std::vector<std::string>&);

  /**
   * Overloads writing equation systems, this is done because when overloading
   * write_nodal_data there would be no way to export cell centered data
   */
  virtual void write_nodal_data(const std::string& fname, const std::vector<Number>& soln, const std::vector<std::string>& names); 

  virtual void read (const std::string& );
  
  void write_refined_mesh(bool refined_mesh);

  virtual void write (const std::string& );  

private:
  void write_ascii(const std::string& fname, const std::vector<Number>* soln, const std::vector<std::string>* names); 

  void cell_connectivity(const Elem* elem, std::vector<unsigned int>& vtk_cell_connectivity);
  
  unsigned int cell_type(const Elem* elem);

  unsigned int cell_offset(const Elem* elem);

  bool _write_refined_mesh;

  MeshData* _mesh_data;
};

inline
VTKIO_NOVTK::VTKIO_NOVTK (MeshBase& mesh, MeshData* mesh_data) :
	MeshInput<MeshBase> (mesh),
	MeshOutput<MeshBase>(mesh),
        _write_refined_mesh(true),
	_mesh_data(mesh_data)
{
//  untested();
}



inline
VTKIO_NOVTK::VTKIO_NOVTK (const MeshBase& mesh, MeshData* mesh_data) :
	MeshOutput<MeshBase>(mesh),
        _write_refined_mesh(true),
	_mesh_data(mesh_data)
{
//  untested();
}
#endif // #define __vtk_io_vtk_indip_h__
