// Local includes
#include "VTKDiscontinuousWriter.h"

void VTKDiscontinuousWriter::cell_connectivity (const Elem* elem, std::vector<unsigned int>& vtk_cell_connectivity)
{
  switch(elem->type())
  {
  case TRI3:
    vtk_cell_connectivity.resize(3);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    break;
  case TRI6:
    vtk_cell_connectivity.resize(6);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    break;
  case QUAD4:
    vtk_cell_connectivity.resize(4);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    break;
  case QUAD8:
    vtk_cell_connectivity.resize(8);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    break;
  case QUAD9:
    vtk_cell_connectivity.resize(8);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    vtk_cell_connectivity[8] = 8;
    break;
  case TET4:
    vtk_cell_connectivity.resize(4);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    break;
  case TET10:
    vtk_cell_connectivity.resize(10);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    vtk_cell_connectivity[8] = 8;
    vtk_cell_connectivity[9] = 9;
    break;
  case PRISM6:
    vtk_cell_connectivity.resize(6);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 2;
    vtk_cell_connectivity[2] = 1;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 5;
    vtk_cell_connectivity[5] = 4;
    break;
  case PRISM15:
    vtk_cell_connectivity.resize(15);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    vtk_cell_connectivity[8] = 8;
    vtk_cell_connectivity[9] = 12;
    vtk_cell_connectivity[10] = 13;
    vtk_cell_connectivity[11] = 14;
    vtk_cell_connectivity[12] = 9;
    vtk_cell_connectivity[13] = 10;
    vtk_cell_connectivity[14] = 11;
    break;
  case PRISM18:
    vtk_cell_connectivity.resize(18);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    vtk_cell_connectivity[8] = 8;
    vtk_cell_connectivity[9] = 12;
    vtk_cell_connectivity[10] = 13;
    vtk_cell_connectivity[11] = 14;
    vtk_cell_connectivity[12] = 9;
    vtk_cell_connectivity[13] = 10;
    vtk_cell_connectivity[14] = 11;
    vtk_cell_connectivity[15] = 15;
    vtk_cell_connectivity[16] = 16;
    vtk_cell_connectivity[17] = 17;
    break;
  case HEX8:
    vtk_cell_connectivity.resize(8);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    break;
  case HEX20:
    vtk_cell_connectivity.resize(20);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    vtk_cell_connectivity[8] = 8;
    vtk_cell_connectivity[9] = 9;
    vtk_cell_connectivity[10] = 10;
    vtk_cell_connectivity[11] = 11;
    vtk_cell_connectivity[12] = 16;
    vtk_cell_connectivity[13] = 17;
    vtk_cell_connectivity[14] = 18;
    vtk_cell_connectivity[15] = 19;
    vtk_cell_connectivity[16] = 12;
    vtk_cell_connectivity[17] = 13;
    vtk_cell_connectivity[18] = 14;
    vtk_cell_connectivity[19] = 15;
    break;
  case HEX27:
    vtk_cell_connectivity.resize(27);
    vtk_cell_connectivity[0] = 0;
    vtk_cell_connectivity[1] = 1;
    vtk_cell_connectivity[2] = 2;
    vtk_cell_connectivity[3] = 3;
    vtk_cell_connectivity[4] = 4;
    vtk_cell_connectivity[5] = 5;
    vtk_cell_connectivity[6] = 6;
    vtk_cell_connectivity[7] = 7;
    vtk_cell_connectivity[8] = 8;
    vtk_cell_connectivity[9] = 9;
    vtk_cell_connectivity[10] = 10;
    vtk_cell_connectivity[11] = 11;
    vtk_cell_connectivity[12] = 16;
    vtk_cell_connectivity[13] = 17;
    vtk_cell_connectivity[14] = 18;
    vtk_cell_connectivity[15] = 19;
    vtk_cell_connectivity[16] = 12;
    vtk_cell_connectivity[17] = 13;
    vtk_cell_connectivity[18] = 14;
    vtk_cell_connectivity[19] = 15;
    vtk_cell_connectivity[20] = 24;
    vtk_cell_connectivity[21] = 22;
    vtk_cell_connectivity[22] = 21;
    vtk_cell_connectivity[23] = 23;
    vtk_cell_connectivity[24] = 20;
    vtk_cell_connectivity[25] = 25;
    vtk_cell_connectivity[26] = 26;
    break;
  case PYRAMID5:
    vtk_cell_connectivity.resize(5);
    vtk_cell_connectivity[0] = 3;
    vtk_cell_connectivity[1] = 2;
    vtk_cell_connectivity[2] = 1;
    vtk_cell_connectivity[3] = 0;
    vtk_cell_connectivity[4] = 4;
    break;
  default:
    {
      std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
      libmesh_error();
    }
  }
}

unsigned int VTKDiscontinuousWriter::cell_offset (const Elem* elem)
{
    std::vector<unsigned int> conn;
    elem->connectivity(0,VTK,conn);
    return conn.size();
}

unsigned int VTKDiscontinuousWriter::cell_type(const Elem* elem)
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
  case PRISM6:     
    celltype = 13;   //VTK_WEDGE
    break;// 12
  case PRISM15:   
    celltype = 26;   //VTK_HIGHER_ORDER_WEDGE
    break;// 13
  case PYRAMID5:
    celltype = 14;   //VTK_PYRAMID
    break;// 14
  case QUAD9:      
    celltype = 28;   //VTK_BIQUADRATIC_QUAD
    break;// 15
  case HEX27:  
    celltype = 29;   //VTK_TRIQUADRATIC_HEXAHEDRON
    break;// 16
  case PRISM18:   
    celltype = 32;   //VTK_HIGHER_ORDER_WEDGE
    break;// 17
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
      libmesh_error();
    }
  }
  return celltype;
}

void VTKDiscontinuousWriter::write_ascii(const std::string& work_dir, const unsigned int& step, const MeshBase& mesh, const std::vector<double>& soln)
{ 
  unsigned int dim = mesh.spatial_dimension();
  std::vector<std::string> names;
  names.push_back("u");
  names.push_back("v");
  names.push_back("w");
  names.push_back("p");
  names.push_back("taux");
  names.push_back("tauy");
  names.push_back("tauz");
  MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator pit = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator end = mesh.active_local_elements_end();

  const unsigned int n_nodes = soln.size()/names.size();
  unsigned int n_elems = 0;
  unsigned int vtk_id = 0;
  for (it = (pit) ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    n_elems++;
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
      vtk_id++;
    }
  }
  
  //write pvtu header
  if (libMesh::processor_id() == 0)
  {
    char fpname[1024];
    sprintf(fpname,"gnuid_%06d",step);
    char fname[1024];
    sprintf(fname,"g_%06d",step);
    FILE* parFile=fopen((work_dir+"/"+fpname+".pvtu").c_str(),"w");
    fprintf (parFile,"%s","<VTKFile type=\"PUnstructuredGrid\"");
    fprintf (parFile,"%s"," version=\"0.1\"");
    fprintf (parFile,"%s\n"," byte_order=\"LittleEndian\">");
    fprintf (parFile,"%s\n","   <PUnstructuredGrid GhostLevel=\"0\">");
    fprintf (parFile,"%s\n","         <PPoints>");
    fprintf (parFile,"%s%s%s\n","            <PDataArray type=\"Float32\" NumberOfComponents=\"","3","\">");
    fprintf (parFile,"%s\n","            </PDataArray>");
    fprintf (parFile,"%s\n","         </PPoints>");
    fprintf (parFile,"%s\n","         <PPointData>");
    for(unsigned int v = 0; v < names.size(); v++)
    {
      fprintf (parFile,"%s%s%s\n","            <PDataArray type=\"Float32\" Name=\"",names[v].c_str(),"\" NumberOfComponents=\"1\">");
      fprintf (parFile,"%s\n","            </PDataArray>");
    }
    fprintf (parFile,"%s\n","         </PPointData>");
    for (unsigned int i = 0; i < libMesh::n_processors(); i++)
    {
      fprintf (parFile,"%s%s%s%04d%s\n","      <Piece  Source=\"",(work_dir+"/"+fname).c_str(),"_",i,"p.vtu\">");
      fprintf (parFile,"%s\n","      </Piece>  ");        
    }
    fprintf (parFile,"%s\n","   </PUnstructuredGrid>");
    fprintf (parFile,"%s\n","</VTKFile>");
    fclose (parFile);       
  }

  char fname[1024];
  sprintf(fname,"g_%06d_%04dp",step,libMesh::processor_id());
  // write header
  FILE* pFile=fopen((work_dir+"/"+fname+".vtu").c_str(),"w");
  fprintf (pFile,"%s","<VTKFile type=\"UnstructuredGrid\"");
  fprintf (pFile,"%s"," version=\"0.1\"");
  fprintf (pFile,"%s\n"," byte_order=\"LittleEndian\">");
  fprintf (pFile,"%s\n","   <UnstructuredGrid>");
  fprintf (pFile,"%s","      <Piece  ");
  fprintf (pFile,"%s%d%s","NumberOfPoints=\"",vtk_id,"\"  ");
  fprintf (pFile,"%s%d%s\n","NumberOfCells=\"",n_elems,"\">");

// write mesh nodes  
  fprintf (pFile,"%s\n","         <Points>");
  dim = 3;
  fprintf (pFile,"%s%d%s\n","            <DataArray type=\"Float32\" NumberOfComponents=\"",dim,"\" format=\"ascii\">");
  for (it = (pit) ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    for (unsigned int n=0; n<elem->n_nodes(); n++)
    {
      for(unsigned int i = 0; i < 3; i++)
      {
        fprintf (pFile,"%f ",elem->point(n).operator()(i));
      }
      fprintf (pFile,"%s\n","");
    }
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","         </Points>");
  
////  write solutions 
  fprintf (pFile,"%s\n","         <PointData>");
  const unsigned int n_sys_vars = 7;
  for(unsigned int j=0;j<n_sys_vars;++j)
  {
    const std::string& varname = names[j];
    fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" Name=\"",varname.c_str(),"\" NumberOfComponents=\"1\" format=\"ascii\">");
    for(unsigned int i = 0; i < n_nodes; i++)
    {
      fprintf (pFile,"%f ",soln[(i*n_sys_vars)+j]);
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
  unsigned int nodes_counter = 0;
  for (it = (pit) ; it != end; ++it)
  {
    std::vector<unsigned int> vtk_cell_connectivity;  
    const Elem *elem  = (*it);
    this->cell_connectivity(elem,vtk_cell_connectivity);
    for (unsigned int i = 0; i < vtk_cell_connectivity.size(); i++)
    {
       fprintf (pFile, "%d ",vtk_cell_connectivity[i]+nodes_counter);
    }
    nodes_counter += vtk_cell_connectivity.size();
    fprintf (pFile,"%s\n","");
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
  unsigned int offset = 0;
  for (it = (pit) ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    offset += this->cell_offset(elem);
    fprintf (pFile, "%d ",offset);
  } 
  fprintf (pFile,"%s\n","");
  fprintf (pFile,"%s\n","            </DataArray>");
  fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">");
  for (it = (pit) ; it != end; ++it)
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

void VTKDiscontinuousWriter::build_discontinuous_solution_vector(MeshBase& mesh, libMesh::EquationSystems& systems, std::vector<Real>& soln)
{
  parallel_only();
  
  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nv  = 7;
  unsigned int tw = 0;
  const Real viscosity = systems.parameters.get<Real>("viscosity");
  {
    MeshBase::element_iterator       it  = mesh.active_local_elements_begin();
    const MeshBase::element_iterator end = mesh.active_local_elements_end(); 
    for ( ; it != end; ++it)
      tw += (*it)->n_nodes();
  }
  soln.resize(tw*nv);

  std::vector<Number> advdiff_soln;
  const System& systemAdvDiff = systems.get_system("AdvDiff");
  const unsigned int u_var = systemAdvDiff.variable_number ("u");
  const unsigned int v_var = systemAdvDiff.variable_number ("v");
  const unsigned int w_var = systemAdvDiff.variable_number ("w");
  systemAdvDiff.update_global_solution (advdiff_soln);
  
  std::vector<Number> pproj_soln;
  const System& systemPProj = systems.get_system("PProj");
  const unsigned int p_var = systemPProj.variable_number ("p");
  systemPProj.update_global_solution (pproj_soln);


  std::vector<Real> elem_soln_u;   
  std::vector<Real> elem_soln_v;   
  std::vector<Real> elem_soln_w;   
  std::vector<Real> elem_soln_p;   
  std::vector<unsigned int> dof_indices_u; 
  std::vector<unsigned int> dof_indices_v; 
  std::vector<unsigned int> dof_indices_w; 
  std::vector<unsigned int> dof_indices_p; 
  std::vector<Real> nodal_soln_u;  
  std::vector<Real> nodal_soln_v;  
  std::vector<Real> nodal_soln_w;  
  std::vector<Real> nodal_soln_p;  
  const FEType& fe_type  = systemAdvDiff.variable_type(u_var);
  const FEType& fe_type_p  = systemPProj.variable_type(p_var);
  const DofMap &dof_map  = systemAdvDiff.get_dof_map();
  const DofMap &dof_map_p  = systemPProj.get_dof_map();
  AutoPtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  const std::vector<std::vector<Real> >&  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe_elem_face->get_dphi();
  const std::vector<Point>& qface_normals = fe_elem_face->get_normals();
  
  MeshBase::element_iterator       it  = mesh.active_local_elements_begin();
  const MeshBase::element_iterator end = mesh.active_local_elements_end(); 
  unsigned int nn=0;

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;      
  
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    dof_map.dof_indices (elem, dof_indices_w, w_var);
    dof_map_p.dof_indices (elem, dof_indices_p, p_var);
    elem_soln_u.resize(dof_indices_u.size());
    elem_soln_v.resize(dof_indices_v.size());
    elem_soln_w.resize(dof_indices_w.size());
    elem_soln_p.resize(dof_indices_p.size());
    for (unsigned int i=0; i<dof_indices_u.size(); i++)
    {
      elem_soln_u[i] = advdiff_soln[dof_indices_u[i]];
      elem_soln_v[i] = advdiff_soln[dof_indices_v[i]];
      elem_soln_w[i] = advdiff_soln[dof_indices_w[i]];
    }
    for (unsigned int i=0; i<dof_indices_p.size(); i++)
    {
      elem_soln_p[i] = pproj_soln[dof_indices_p[i]];
    }
    
    FEInterface::nodal_soln (dim, fe_type, elem, elem_soln_u, nodal_soln_u);
    FEInterface::nodal_soln (dim, fe_type, elem, elem_soln_v, nodal_soln_v);
    FEInterface::nodal_soln (dim, fe_type, elem, elem_soln_w, nodal_soln_w);
    FEInterface::nodal_soln (dim, fe_type_p, elem, elem_soln_p, nodal_soln_p);
    libmesh_assert (nodal_soln_u.size() == elem->n_nodes());
    libmesh_assert (nodal_soln_v.size() == elem->n_nodes());
    libmesh_assert (nodal_soln_w.size() == elem->n_nodes());
    libmesh_assert (nodal_soln_p.size() == elem->n_nodes());
    for (unsigned int n = 0; n<nodal_soln_u.size(); n++)
    { 
      soln[nv*(nn+n)]     += nodal_soln_u[n];
      soln[nv*(nn+n) + 1] += nodal_soln_v[n];
      soln[nv*(nn+n) + 2] += nodal_soln_w[n];
      soln[nv*(nn+n) + 3] += nodal_soln_p[n];
    }
    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      if (elem->neighbor(side) == NULL)
      {
        AutoPtr<Elem> elem_side (elem->build_side(side));
        std::vector<Point> refspace_nodes;
        FEBase::get_refspace_nodes(elem_side->type(), refspace_nodes);
        fe_elem_face->reinit(elem, side, 1.e-10, &refspace_nodes);
        for (unsigned int qp=0; qp<refspace_nodes.size(); qp++)
	{
          Real du_dx = 0., du_dy = 0., du_dz = 0.;
          Real dv_dx = 0., dv_dy = 0., dv_dz = 0.;
          Real dw_dx = 0., dw_dy = 0., dw_dz = 0.;
          for (unsigned int i=0; i<dof_indices_u.size(); i++)
          {
            du_dx += dphi_face[i][qp](0) * elem_soln_u[i];
            du_dy += dphi_face[i][qp](1) * elem_soln_u[i];
            du_dz += dphi_face[i][qp](2) * elem_soln_u[i];
            dv_dx += dphi_face[i][qp](0) * elem_soln_v[i];
            dv_dy += dphi_face[i][qp](1) * elem_soln_v[i];
            dv_dz += dphi_face[i][qp](2) * elem_soln_v[i];
            dw_dx += dphi_face[i][qp](0) * elem_soln_w[i];
            dw_dy += dphi_face[i][qp](1) * elem_soln_w[i];
            dw_dz += dphi_face[i][qp](2) * elem_soln_w[i];
	  }
          Real sx = 2. * viscosity * du_dx;
          Real sy = 2. * viscosity * dv_dy;
          Real sz = 2. * viscosity * dw_dz;
          Real txy = viscosity * (dv_dx + du_dy);
          Real txz = viscosity * (dw_dx + du_dz);
          Real tyz = viscosity * (dw_dy + dv_dz);
	  Real taux = - sx  * qface_normals[qp](0) - txy * qface_normals[qp](1) - txz * qface_normals[qp](2); 
	  Real tauy = - txy * qface_normals[qp](0) - sy  * qface_normals[qp](1) - tyz * qface_normals[qp](2); 
	  Real tauz = - txz * qface_normals[qp](0) - tyz * qface_normals[qp](1) - sz  * qface_normals[qp](2); 
	  Real tn = taux * qface_normals[qp](0) + tauy * qface_normals[qp](1) + tauz * qface_normals[qp](2);
	  taux -= tn * qface_normals[qp](0);
	  tauy -= tn * qface_normals[qp](1);
	  tauz -= tn * qface_normals[qp](2);
          for (unsigned int i = 0; i<elem->n_nodes(); i++)
	  {
            if (elem->node(i) == elem_side->node(qp))
	    {
              soln[nv*(nn+i) + 4] += taux;
              soln[nv*(nn+i) + 5] += tauy;
              soln[nv*(nn+i) + 6] += tauz;
	    }
	  }
        }
      }
    }
    nn += nodal_soln_u.size();
  }
}
  
