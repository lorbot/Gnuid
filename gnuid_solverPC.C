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

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/boundary_info.h"
#include "libmesh/coupling_matrix.h"
#include "libmesh/fe_interface.h"
#include "libmesh/utility.h"
#include "libmesh/error_vector.h"
#include "libmesh/petsc_linear_solver.h"

#include "gnuid_solverINS.h"
#include "gnuid_bcHelper.h"

void GnuidSolver::assemble_p_proj(EquationSystems& es, const std::string& system_name)
{
  std::cout<<"assembling pressure projection system... ";
  std::cout.flush();

  libmesh_assert (system_name == "PProj");
  GnuidBCHelper bcHelper;
  
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = 3;
  const Real density = es.parameters.get<Real>("density");
  const Real dt= es.parameters.get<Real> ("dt");   
  const Real t = es.parameters.get<Real>("time");
  const Real step = es.parameters.get<unsigned int>("step");
  const unsigned int n_washin_steps = es.parameters.get<unsigned int>("n washin steps");
  const unsigned int n_second_order_start_step = es.parameters.get<unsigned int>("second order start step");
  
  TransientLinearImplicitSystem & systemPProj = es.get_system<TransientLinearImplicitSystem> ("PProj");
  const unsigned int p_var = systemPProj.variable_number ("p");
  
  TransientLinearImplicitSystem & systemAdvDiff = es.get_system<TransientLinearImplicitSystem>("AdvDiff");
  const unsigned int u_var = systemAdvDiff.variable_number("u");
  const unsigned int v_var = systemAdvDiff.variable_number("v");
  const unsigned int w_var = systemAdvDiff.variable_number("w");
  
  FEType fe_type = systemAdvDiff.variable_type(u_var);
  FEType fe_type_p = systemPProj.variable_type(p_var);

  AutoPtr<FEBase> fe  (FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_p  (FEBase::build(dim, fe_type_p));
  AutoPtr<FEBase> fe_elem_face_p(FEBase::build(dim, fe_type_p));
  AutoPtr<FEBase> fe_neighbor_face_p(FEBase::build(dim, fe_type_p));

#ifdef QORDER 
  QGauss qrule (dim, QORDER);
#else
  QGauss qrule (dim, fe_type_p.default_quadrature_order());
#endif
  fe->attach_quadrature_rule (&qrule);
  fe_p->attach_quadrature_rule (&qrule);
  
#ifdef QORDER 
  QGauss qface(dim-1, QORDER);
#else
  QGauss qface(dim-1, fe_type_p.default_quadrature_order());
#endif
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);
  fe_elem_face_p->attach_quadrature_rule(&qface);
  fe_neighbor_face_p->attach_quadrature_rule(&qface);

  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >& psi = fe_p->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_p->get_dphi();
  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >&  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<Real> >&  phi_neighbor_face = fe_neighbor_face->get_phi();
  const std::vector<std::vector<Real> >&  psi_face = fe_elem_face_p->get_phi();
  const std::vector<std::vector<Real> >&  psi_neighbor_face = fe_neighbor_face_p->get_phi();
  const std::vector<Real>& JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point>& qface_normals = fe_elem_face->get_normals();
  const std::vector<Point >& qface_points = fe_elem_face->get_xyz();
  std::vector<Point> qface_neighbor_point;
  
  const DofMap & dof_map = systemPProj.get_dof_map();
  const DofMap & dof_map_u = systemAdvDiff.get_dof_map();

  DenseMatrix<Real> Ke_p;
  DenseVector<Real> Fe_p;
  DenseVector<Real> Fee_p;
  DenseVector<Real> Fnn_p;
  
  std::vector<unsigned int> dof_indices_p;
  std::vector<unsigned int> neighbor_dof_indices_p;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> neighbor_dof_indices_u;
  std::vector<unsigned int> neighbor_dof_indices_v;
  std::vector<unsigned int> neighbor_dof_indices_w;
  
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

  unsigned int n_flow_bc = es.parameters.get<unsigned int>("n flow bound");
  std::vector<Real> surfaceArea(n_flow_bc,0.);
  std::vector<Real> surfaceFlow(n_flow_bc,0.);
  for (el = mesh.active_local_elements_begin() ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    for (unsigned int side=0; side<elem->n_sides(); side++)
      if (elem->neighbor(side) == NULL)
      {
        const unsigned short boundary_id = mesh.boundary_info->boundary_id(elem,side);
        for (unsigned int lid = 0; lid < n_flow_bc; lid++)
        {
          char flow_bc_id[1024];
          sprintf(flow_bc_id,"flow_bc_%02d_id",lid);
          if (boundary_id == es.parameters.get<unsigned short>(flow_bc_id))
          {
            dof_map_u.dof_indices (elem, dof_indices_u, u_var);
            dof_map_u.dof_indices (elem, dof_indices_v, v_var);
            dof_map_u.dof_indices (elem, dof_indices_w, w_var);
            const unsigned int n_dofs = dof_indices_u.size();
            fe_elem_face->reinit(elem, side);
            for (unsigned int qp=0; qp<qface.n_points(); qp++)
            {
              surfaceArea[lid] += JxW_face[qp];
              Real u = 0;
              Real v = 0;
              Real w = 0;
              for (unsigned int i = 0; i < n_dofs; i ++)
              {
                u += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_u[i]);
                v += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_v[i]);
                w += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_w[i]);
              }
              surfaceFlow[lid] += JxW_face[qp] * (u * qface_normals[qp](0) + v * qface_normals[qp](1) + w * qface_normals[qp](2));
            }
          }
	}
      }
  }
  std::vector<Real> u_mean_n(n_flow_bc);
  for (unsigned int lid = 0; lid < n_flow_bc; lid++)
    u_mean_n[lid] = surfaceFlow[lid]/surfaceArea[lid];
  

  for (el = mesh.active_local_elements_begin() ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    fe->reinit(elem);
    fe_p->reinit(elem);

    dof_map.dof_indices (elem, dof_indices_p, p_var);
    dof_map_u.dof_indices (elem, dof_indices_u, u_var);
    dof_map_u.dof_indices (elem, dof_indices_v, v_var);
    dof_map_u.dof_indices (elem, dof_indices_w, w_var);
    const unsigned int n_uvw_e_dofs = dof_indices_u.size();
    const unsigned int n_p_e_dofs = dof_indices_p.size();

    Ke_p.resize(n_p_e_dofs, n_p_e_dofs);
    Fe_p.resize(n_p_e_dofs);
 
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      RealVectorValue  grad_p_old (0.,0.,0.);
      for (unsigned int i=0; i<n_p_e_dofs; i++)
      {
         grad_p_old(0) += dpsi[i][qp](0) * systemPProj.old_solution(dof_indices_p[i]);
         grad_p_old(1) += dpsi[i][qp](1) * systemPProj.old_solution(dof_indices_p[i]);
         grad_p_old(2) += dpsi[i][qp](2) * systemPProj.old_solution(dof_indices_p[i]);
      }
      Real divU = 0.0, divU_old = 0.0, divU_older = 0.0;
      Real divU_extr = 0.0;
      for (unsigned int i=0; i<n_uvw_e_dofs; i++)
      {
        RealVectorValue U( systemAdvDiff.current_solution(dof_indices_u[i]),
                           systemAdvDiff.current_solution(dof_indices_v[i]),
                           systemAdvDiff.current_solution(dof_indices_w[i]));
        divU += U * dphi[i][qp];
      }
      
      if (step > n_second_order_start_step)
        divU_extr = (3./2.) * divU;  
      else
        divU_extr = divU;  
        
      for (unsigned int i=0; i<n_p_e_dofs; i++)
      {
        for (unsigned int j=0; j<n_p_e_dofs; j++)
        {
          Ke_p(i,j) += JxW[qp]*(dpsi[i][qp]*dpsi[j][qp]);
        }
        Fe_p(i) += JxW[qp] * (grad_p_old * dpsi[i][qp]);
        Fe_p(i) -= JxW[qp] * ((density/dt) * psi[i][qp] * divU_extr);
      }
    }
    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      if (elem->neighbor(side) == NULL)
      {
        const unsigned short boundary_id = mesh.boundary_info->boundary_id(elem,side);
	std::string bc_type;
        bcHelper.init_bcData(es, boundary_id, bc_type);
	int dirichletProfile = bc_type.compare("dirichlet_profile");
	int dirichletWall    = bc_type.compare("dirichlet_wall");
        if (dirichletProfile  == 0 || dirichletWall == 0)
	{
          fe_elem_face->reinit(elem, side);
          fe_elem_face_p->reinit(elem, side);
          RealVectorValue U_bc;
	  if (dirichletWall == 0)
	    U_bc.zero();
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            if (dirichletProfile == 0)
              bcHelper.compute_dirichletIOProfile(qface_points[qp],t,U_bc);
            const Real u_bc_n = U_bc * qface_normals[qp];
            RealVectorValue U (0.,0.,0.); 
            for (unsigned int i=0; i<n_uvw_e_dofs; i++)
            {
              U(0) += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_u[i]);
              U(1) += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_v[i]);
              U(2) += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_w[i]);
            }
            const Real u_n = U * qface_normals[qp];
            Real u_n_jump = 0.;
            if (step > n_second_order_start_step)
              u_n_jump = (3./2.) * (u_n - u_bc_n);
            else
              u_n_jump = (u_n - u_bc_n);
            for (unsigned int i=0; i<n_p_e_dofs; i++)
              Fe_p(i) += JxW_face[qp] * (density/dt) * u_n_jump * psi_face[i][qp];
          }
	}
        else if (bc_type.compare("dirichlet_defective") == 0)
	{
	  unsigned int lid;
	  RealVectorValue U_bc;
          bcHelper.compute_dirichletDefectiveData(t, U_bc, lid);
          fe_elem_face->reinit(elem, side);
          fe_elem_face_p->reinit(elem, side);
	  Real u_bc_n = U_bc * qface_normals[0];
          Real u_n_jump = 0;
          if (step > n_second_order_start_step)
            u_n_jump = (3./2.) * (u_mean_n[lid] - u_bc_n);
          else
            u_n_jump = (u_mean_n[lid] - u_bc_n);
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
            for (unsigned int i=0; i<n_p_e_dofs; i++)
               Fe_p(i) += JxW_face[qp] * (density/dt) * u_n_jump * psi_face[i][qp];
	}
	else if (bc_type.compare("neumann") == 0)
	{
          fe_elem_face->reinit(elem, side);
          fe_elem_face_p->reinit(elem, side);
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            Real penalty_bc = 1e10;
            Real pressure_outlet = 0.;
            for (unsigned int i=0; i<n_p_e_dofs; i++)
            {
              for (unsigned int j=0; j<n_p_e_dofs; j++)
                Ke_p(i,j) += JxW_face[qp] * penalty_bc * psi_face[i][qp] * psi_face[j][qp];
              Fe_p(i) += JxW_face[qp] * pressure_outlet * penalty_bc * psi_face[i][qp];
            }          
          }
        }
	else
	{
	  std::cerr<<"Unknown boundary condtion: "<<bc_type<<std::endl;
	  libmesh_error();
	}
      } 
      else
      {
        const Elem* neighbor = elem->neighbor(side);
        const unsigned int elem_id = elem->id();
        const unsigned int neighbor_id = neighbor->id();
        if ((neighbor->active() && (neighbor->level() == elem->level()) && (elem_id < neighbor_id)) || (neighbor->level() < elem->level()))
        {
          AutoPtr<Elem> elem_side (elem->build_side(side));
          fe_elem_face->reinit(elem, side);
          fe_elem_face_p->reinit(elem, side);
	  unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
          fe_neighbor_face->side_map (neighbor, elem_side.get(), side_neighbor, qface.get_points(), qface_neighbor_point);
          fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);
          fe_neighbor_face_p->reinit(neighbor, &qface_neighbor_point);

          dof_map.dof_indices (neighbor, neighbor_dof_indices_p, p_var);
          dof_map_u.dof_indices (neighbor, neighbor_dof_indices_u, u_var);
          dof_map_u.dof_indices (neighbor, neighbor_dof_indices_v, v_var);
          dof_map_u.dof_indices (neighbor, neighbor_dof_indices_w, w_var);
          const unsigned int n_uvw_n_dofs = neighbor_dof_indices_u.size();
          const unsigned int n_p_n_dofs = neighbor_dof_indices_p.size();
   
          Fee_p.resize(n_p_e_dofs);
          Fnn_p.resize(n_p_n_dofs);
          
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            RealVectorValue U (0.,0.,0.); 
            for (unsigned int i=0; i<n_uvw_e_dofs; i++)
            {
               U(0) += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_u[i]);
               U(1) += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_v[i]);
               U(2) += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_w[i]);
            }
            const Real u_n = U * qface_normals[qp];
            RealVectorValue U_neighbor (0.,0.,0.); 
            for (unsigned int i=0; i<n_uvw_n_dofs; i++)
            {
               U_neighbor(0) += phi_neighbor_face[i][qp] * systemAdvDiff.current_solution(neighbor_dof_indices_u[i]);
               U_neighbor(1) += phi_neighbor_face[i][qp] * systemAdvDiff.current_solution(neighbor_dof_indices_v[i]);
               U_neighbor(2) += phi_neighbor_face[i][qp] * systemAdvDiff.current_solution(neighbor_dof_indices_w[i]);
            }
            const Real u_n_neighbor = U_neighbor * qface_normals[qp];
            Real u_n_jump = 0.;
            if (step > n_second_order_start_step)
              u_n_jump = (3./2.) * (u_n - u_n_neighbor);
            else
              u_n_jump = (u_n - u_n_neighbor);
            for (unsigned int i=0; i<n_p_e_dofs; i++)          // elem elem matrix
              Fee_p(i) += 0.5 * (density/dt) * JxW_face[qp] * u_n_jump * psi_face[i][qp];

            for (unsigned int i=0; i<n_p_n_dofs; i++) // neighbor neighbor matrix
              Fnn_p(i) += 0.5  * (density/dt) * JxW_face[qp] * u_n_jump * psi_neighbor_face[i][qp];
          }
          systemPProj.rhs->add_vector(Fee_p, dof_indices_p);
          systemPProj.rhs->add_vector(Fnn_p, neighbor_dof_indices_p);
        }
      }
    }
    systemPProj.matrix->add_matrix(Ke_p, dof_indices_p);
    systemPProj.rhs->add_vector(Fe_p, dof_indices_p);
  }
  std::cout<<"done"<<std::endl;
}

