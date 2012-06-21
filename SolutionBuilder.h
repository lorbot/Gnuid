#ifndef __solution_builder__h
#define __solution_builder__h

#include "mesh.h"
#include "dof_map.h"
#include "equation_systems.h"
#include "system.h"
#include "parallel.h"
#include "fe.h"

#include "fe_interface.h"

class SolutionBuilder
{
public:
  SolutionBuilder()
  {
    ;
  }

  ~SolutionBuilder() {}

  void build_discontinuous_solution_vector(Mesh& mesh, libMesh::EquationSystems& systems, std::vector<Real>& solution);
  
private:   
};


void SolutionBuilder::build_discontinuous_solution_vector(Mesh& mesh, libMesh::EquationSystems& systems, std::vector<Real>& soln)
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
        fe_elem_face->reinit(elem, side, 1, &refspace_nodes);
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
	  Real taux = sx  * qface_normals[qp](0) + txy * qface_normals[qp](1) + txz * qface_normals[qp](2); 
	  Real tauy = txy * qface_normals[qp](0) + sy  * qface_normals[qp](1) + tyz * qface_normals[qp](2); 
	  Real tauz = txz * qface_normals[qp](0) + tyz * qface_normals[qp](1) + sz  * qface_normals[qp](2); 
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

#endif
