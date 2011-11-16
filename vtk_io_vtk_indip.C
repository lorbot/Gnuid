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

#include <fstream>

// Local includes
#include "vtk_io_vtk_indip.h"
#include "mesh_base.h"
#include "mesh.h"
#include "elem.h"
#include "equation_systems.h"
#include "numeric_vector.h"

using namespace libMesh;

void VTKIO_NOVTK::cell_connectivity (const Elem* elem, std::vector<unsigned int>& vtk_cell_connectivity)
{
  switch(elem->type())
  {
  case TRI6:
    vtk_cell_connectivity.resize(6);
    vtk_cell_connectivity[0] = elem->node(0);
    vtk_cell_connectivity[1] = elem->node(1);
    vtk_cell_connectivity[2] = elem->node(2);
    vtk_cell_connectivity[3] = elem->node(3);
    vtk_cell_connectivity[4] = elem->node(4);
    vtk_cell_connectivity[5] = elem->node(5);
    break;
  case QUAD8:
    vtk_cell_connectivity.resize(8);
    vtk_cell_connectivity[0] = elem->node(0);
    vtk_cell_connectivity[1] = elem->node(1);
    vtk_cell_connectivity[2] = elem->node(2);
    vtk_cell_connectivity[3] = elem->node(3);
    vtk_cell_connectivity[4] = elem->node(4);
    vtk_cell_connectivity[5] = elem->node(5);
    vtk_cell_connectivity[6] = elem->node(6);
    vtk_cell_connectivity[7] = elem->node(7);
    break;
  case QUAD9:
    vtk_cell_connectivity.resize(8);
    vtk_cell_connectivity[0] = elem->node(0);
    vtk_cell_connectivity[1] = elem->node(1);
    vtk_cell_connectivity[2] = elem->node(2);
    vtk_cell_connectivity[3] = elem->node(3);
    vtk_cell_connectivity[4] = elem->node(4);
    vtk_cell_connectivity[5] = elem->node(5);
    vtk_cell_connectivity[6] = elem->node(6);
    vtk_cell_connectivity[7] = elem->node(7);
    vtk_cell_connectivity[8] = elem->node(8);
    break;
  case PRISM15:
    vtk_cell_connectivity.resize(15);
    vtk_cell_connectivity[0] = elem->node(0);
    vtk_cell_connectivity[1] = elem->node(1);
    vtk_cell_connectivity[2] = elem->node(2);
    vtk_cell_connectivity[3] = elem->node(3);
    vtk_cell_connectivity[4] = elem->node(4);
    vtk_cell_connectivity[5] = elem->node(5);
    vtk_cell_connectivity[6] = elem->node(6);
    vtk_cell_connectivity[7] = elem->node(7);
    vtk_cell_connectivity[8] = elem->node(8);
    vtk_cell_connectivity[9] = elem->node(12);
    vtk_cell_connectivity[10] = elem->node(13);
    vtk_cell_connectivity[11] = elem->node(14);
    vtk_cell_connectivity[12] = elem->node(9);
    vtk_cell_connectivity[13] = elem->node(10);
    vtk_cell_connectivity[14] = elem->node(11);
    break;
  case PRISM18:
    vtk_cell_connectivity.resize(18);
    vtk_cell_connectivity[0] = elem->node(0);
    vtk_cell_connectivity[1] = elem->node(1);
    vtk_cell_connectivity[2] = elem->node(2);
    vtk_cell_connectivity[3] = elem->node(3);
    vtk_cell_connectivity[4] = elem->node(4);
    vtk_cell_connectivity[5] = elem->node(5);
    vtk_cell_connectivity[6] = elem->node(6);
    vtk_cell_connectivity[7] = elem->node(7);
    vtk_cell_connectivity[8] = elem->node(8);
    vtk_cell_connectivity[9] = elem->node(12);
    vtk_cell_connectivity[10] = elem->node(13);
    vtk_cell_connectivity[11] = elem->node(14);
    vtk_cell_connectivity[12] = elem->node(9);
    vtk_cell_connectivity[13] = elem->node(10);
    vtk_cell_connectivity[14] = elem->node(11);
    vtk_cell_connectivity[15] = elem->node(15);
    vtk_cell_connectivity[16] = elem->node(16);
    vtk_cell_connectivity[17] = elem->node(17);
    break;
  case HEX27:
    vtk_cell_connectivity.resize(27);
    vtk_cell_connectivity[0] = elem->node(0);
    vtk_cell_connectivity[1] = elem->node(1);
    vtk_cell_connectivity[2] = elem->node(2);
    vtk_cell_connectivity[3] = elem->node(3);
    vtk_cell_connectivity[4] = elem->node(4);
    vtk_cell_connectivity[5] = elem->node(5);
    vtk_cell_connectivity[6] = elem->node(6);
    vtk_cell_connectivity[7] = elem->node(7);
    vtk_cell_connectivity[8] = elem->node(8);
    vtk_cell_connectivity[9] = elem->node(9);
    vtk_cell_connectivity[10] = elem->node(10);
    vtk_cell_connectivity[11] = elem->node(11);
    vtk_cell_connectivity[12] = elem->node(16);
    vtk_cell_connectivity[13] = elem->node(17);
    vtk_cell_connectivity[14] = elem->node(18);
    vtk_cell_connectivity[15] = elem->node(19);
    vtk_cell_connectivity[16] = elem->node(12);
    vtk_cell_connectivity[17] = elem->node(13);
    vtk_cell_connectivity[18] = elem->node(14);
    vtk_cell_connectivity[19] = elem->node(15);
    vtk_cell_connectivity[20] = elem->node(24);
    vtk_cell_connectivity[21] = elem->node(22);
    vtk_cell_connectivity[22] = elem->node(21);
    vtk_cell_connectivity[23] = elem->node(23);
    vtk_cell_connectivity[24] = elem->node(20);
    vtk_cell_connectivity[25] = elem->node(25);
    vtk_cell_connectivity[26] = elem->node(26);
    break;
  default:
    elem->connectivity(0,VTK,vtk_cell_connectivity);
  }
}

unsigned int VTKIO_NOVTK::cell_offset (const Elem* elem)
{
  std::vector<unsigned int> conn;
  this->cell_connectivity(elem,conn);
  return conn.size();
}

unsigned int VTKIO_NOVTK::cell_type(const Elem* elem)
{
  unsigned int celltype = 0; // initialize to something to avoid compiler warning
  
  switch(elem->type())
  {
  case EDGE2:				
    celltype = 3; //VTK_LINE
    break;
  case EDGE3:      
    celltype = 21; //VTK_QUADRATIC_EDGE
    break;// 1
  case TRI3:       
    celltype = 5;  //VTK_TRIANGLE
    break;// 3
  case TRI6:       
    celltype = 22;  //VTK_QUADRATIC_TRIANGLE
    break;// 4
  case QUAD4:      
    celltype = 9;  //VTK_QUAD
    break;// 5
  case QUAD8:      
    celltype = 23;  //VTK_QUADRATIC_QUAD
    break;// 6
  case QUAD9: 
    celltype = 28;   //VTK_BIQUADRATIC_QUAD
    break;//7
  case TET4:      
    celltype = 10;  //VTK_TETRA
    break;// 8
  case TET10:      
    celltype = 24; //VTK_QUADRATIC_TETRA
    break;// 9
  case HEX8:    
    celltype = 12;  //VTK_HEXAHEDRON
    break;// 10
  case HEX20:      
    celltype = 25;  //VTK_QUADRATIC_HEXAHEDRON
    break;// 11
  case HEX27: 
    celltype = 29;  //VTK_TRIQUADRATIC_QUADRATIC_HEXAHEDRON
    break;// 12
  case PRISM6:     
    celltype = 13;   //VTK_WEDGE
    break;// 13
  case PRISM15:
    celltype = 26;   //VTK_QUADRATIC_WEDGE
    break;// 14
  case PRISM18:
    celltype = 32;   //VTK_BIQUADRATIC_QUADRATIC_WEDGE
    break;// 15
  case PYRAMID5:
    celltype = 14;   //VTK_PYRAMID
    break;// 16
  case EDGE4:      
  case INFEDGE2:   
  case INFQUAD4:   
  case INFQUAD6:   
  case INFHEX8:    
  case INFHEX16:   
  case INFHEX18:   
  case INFPRISM6:  
  case INFPRISM12: 
  case NODEELEM:   
  case INVALID_ELEM:
  default:
    {
  	  std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
    }
  }
  return celltype;
}

void VTKIO_NOVTK::write_refined_mesh(bool refined_mesh)
{
   _write_refined_mesh = refined_mesh;
}

void VTKIO_NOVTK::write_nodal_data (const std::string& fname, const std::vector<Number>& soln, const std::vector<std::string>& names)
{
  if (libMesh::processor_id() == 0)
     this->write_ascii (fname, &soln, &names);
}


void VTKIO_NOVTK::write_ascii(const std::string& fname, const std::vector<Number>* soln, const std::vector<std::string>* names)
{ 
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();
  unsigned int dim = mesh.spatial_dimension();
  MeshBase::const_element_iterator* p_it  = NULL;
  MeshBase::const_element_iterator* p_end = NULL; 
  MeshBase::const_element_iterator ait  = mesh.active_elements_begin();
  MeshBase::const_element_iterator aend = mesh.active_elements_end(); 
  MeshBase::const_element_iterator lit  = mesh.level_elements_begin(0);
  MeshBase::const_element_iterator lend = mesh.level_elements_end(0); 
  if (_write_refined_mesh) 
  {
     p_it  = & ait;
     p_end = & aend; 
  }
  else 
  {
     p_it  = & lit;
     p_end = & lend; 
  }
  MeshBase::const_element_iterator it = (*p_it);
  MeshBase::const_element_iterator end = (*p_end);

  const unsigned int n_nodes = mesh.n_nodes();
  std::vector<bool> nodes_bool;
  nodes_bool.assign(n_nodes,false);
  unsigned int n_elems = 0;
  for (it = (*p_it) ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    n_elems++;
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
        nodes_bool[elem->node(i)] = true;
    }
  }
  std::vector<int> nodes_vtk_id;
  nodes_vtk_id.assign(n_nodes, -1);
  unsigned int vtk_id = 0;
  for (unsigned int i = 0; i < nodes_bool.size(); i++)
  {
    if (nodes_bool[i] == true)
    {
       nodes_vtk_id[i] = vtk_id;
       vtk_id++;
    }
  }
  
  // write header
  FILE* pFile=fopen(fname.c_str(),"w");
  fprintf (pFile,"%s","<VTKFile type=\"UnstructuredGrid\"");
  fprintf (pFile,"%s"," version=\"0.1\"");
  fprintf (pFile,"%s\n"," byte_order=\"LittleEndian\">");
  fprintf (pFile,"%s\n","   <UnstructuredGrid>");
  fprintf (pFile,"%s","      <Piece  ");
  fprintf (pFile,"%s%d%s","NumberOfPoints=\"",vtk_id,"\"  ");
  fprintf (pFile,"%s%d%s\n","NumberOfCells=\"",n_elems,"\">");

// write mesh nodes  
  MeshBase::const_node_iterator node_it = mesh.nodes_begin();
  const MeshBase::const_node_iterator node_end = mesh.nodes_end();
  fprintf (pFile,"%s\n","         <Points>");
  fprintf (pFile,"%s%d%s\n","            <DataArray type=\"Float32\" NumberOfComponents=\"",dim,"\" format=\"ascii\">");
  for(node_it=mesh.nodes_begin();node_it!=node_end;++node_it)
  {
    if (nodes_bool[(*node_it)->id()] == true)
    {
      for(unsigned int i = 0; i < dim; i++)
      {
        fprintf (pFile,"%f ",(*node_it)->operator()(i));
      }
      fprintf (pFile,"%s\n","");
    }
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","         </Points>");
  
////  write solutions 
  fprintf (pFile,"%s\n","         <PointData>");
  const unsigned int n_es_vars = names->size();  
  for(unsigned int j=0;j<n_es_vars;++j)
  {
    const std::string& varname = (*names)[j];
    fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" Name=\"",varname.c_str(),"\" NumberOfComponents=\"1\" format=\"ascii\">");
    for(node_it = mesh.nodes_begin();node_it!=node_end;++node_it)
    {
        unsigned int nd_id = (*node_it)->id();
        if (nodes_bool[nd_id] == true)
        {
          fprintf (pFile,"%f ",(*soln)[(nd_id*n_es_vars)+j]);
        }
    }
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");
  }
  fprintf (pFile,"%s\n","         </PointData>");

// write cell data
  fprintf (pFile,"%s\n","         <CellData>");
  fprintf (pFile,"%s\n","         </CellData>");

// write elements connectivity
  fprintf (pFile,"%s\n","         <Cells>");
  fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
  for (it = (*p_it) ; it != end; ++it)
  {
    std::vector<unsigned int> vtk_cell_connectivity;  
    const Elem *elem  = (*it);
    this->cell_connectivity(elem,vtk_cell_connectivity);
    for (unsigned int i = 0; i < vtk_cell_connectivity.size(); i++)
    {
       fprintf (pFile, "%d ",nodes_vtk_id[vtk_cell_connectivity[i]]);
    }
    fprintf (pFile,"%s\n","");
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
  unsigned int offset = 0;
  for (it = (*p_it) ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    offset += this->cell_offset(elem);
    fprintf (pFile, "%d ",offset);
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">");
  for (it = (*p_it) ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    fprintf (pFile, "%d ",this->cell_type(elem));
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","         </Cells>");
  fprintf (pFile,"%s\n","      </Piece>  ");        
  fprintf (pFile,"%s\n","   </UnstructuredGrid>");
  fprintf (pFile,"%s\n","</VTKFile>");
  fclose (pFile);       
} 
  
void VTKIO_NOVTK::read (const std::string& name)
{
}
  
void VTKIO_NOVTK::write (const std::string& name)
{
}	
  
  
  
