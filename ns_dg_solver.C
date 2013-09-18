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

#include "ns_dg_solver.h"
#include "VTKWriter.h"
#include "gnuid_bcHelper.h"

//#define QORDER TENTH 

using namespace libMesh;

void NS_DG_Solver::init()
{
  EquationSystems& es = this->system();
  if (es.parameters.get<bool>("read es file"))
  {
     TransientLinearImplicitSystem & systemAdvDiff = es.get_system<TransientLinearImplicitSystem> ("AdvDiff");
     systemAdvDiff.update();
     systemAdvDiff.attach_assemble_function (NS_DG_Solver::assemble_adv_diff);
     TransientLinearImplicitSystem & systemPProj = es.get_system<TransientLinearImplicitSystem> ("PProj");
     systemPProj.update();
     systemPProj.attach_assemble_function (NS_DG_Solver::assemble_p_proj);
  }
  else
  {
     const unsigned int v_order = es.parameters.get<unsigned int>("velocity approx order");
     const unsigned int p_order = es.parameters.get<unsigned int>("pressure approx order");
     Order v_libmesh = static_cast<Order>(v_order);
     Order p_libmesh = static_cast<Order>(p_order);
     TransientLinearImplicitSystem & systemAdvDiff = es.add_system<TransientLinearImplicitSystem> ("AdvDiff");
     systemAdvDiff.add_variable ("u", v_libmesh, MONOMIAL);
     systemAdvDiff.add_variable ("v", v_libmesh, MONOMIAL);
     systemAdvDiff.add_variable ("w", v_libmesh, MONOMIAL);
     systemAdvDiff.attach_assemble_function (NS_DG_Solver::assemble_adv_diff);
     TransientLinearImplicitSystem & systemPProj = es.add_system<TransientLinearImplicitSystem> ("PProj");
     systemPProj.add_variable ("p", p_libmesh, LAGRANGE);
     systemPProj.attach_assemble_function (NS_DG_Solver::assemble_p_proj);
     es.init();
  }
}

void NS_DG_Solver::solve()
{
  EquationSystems& es = this->system();

  Real time = es.parameters.get<Real>("time");
  const unsigned int n_timesteps = es.parameters.get<unsigned int>("n timesteps");
  const unsigned int n_washin_steps = es.parameters.get<unsigned int>("n washin steps");
  const unsigned int max_nonlinear_steps = es.parameters.get<unsigned int>("nonlinear solver maximum iterations");
  const Real nonlinear_tolerance = es.parameters.get<Real>("nonlinear solver tolerance");
  bool write_continuous = es.parameters.get<bool>("write continuous");
  bool write_discontinuous = es.parameters.get<bool>("write discontinuous");
  const unsigned int write_solution_interval = es.parameters.get<unsigned int>("write output interval");
  const unsigned int write_es_interval = es.parameters.get<unsigned int>("write es interval");
  unsigned int n_nonlinear_steps = 0;

  TransientLinearImplicitSystem & systemAdvDiff = es.get_system<TransientLinearImplicitSystem> ("AdvDiff");
  systemAdvDiff.linear_solver.get()->set_solver_type(GMRES);
//  systemAdvDiff.linear_solver.get()->set_preconditioner_type(ASM_PRECOND);
  TransientLinearImplicitSystem & systemPProj = es.get_system<TransientLinearImplicitSystem> ("PProj");
  systemPProj.linear_solver.get()->set_solver_type(BICGSTAB);
//  systemPProj.linear_solver.get()->set_preconditioner_type(ASM_PRECOND);

  MeshBase& mesh = es.get_mesh();
  VTKWriter vtk_io;
  
  unsigned int step = es.parameters.get<unsigned int>("step");
  step++;
  if (!es.parameters.get<bool>("read es file"))
  {
     std::cout << "################ Initializing pressure field "<< std::endl <<std::endl;
     systemPProj.solve();
  
     std::cout<<"number of linear iterations = "<<systemPProj.n_linear_iterations()<<std::endl;
     std::cout<<"initial linear residual = "<<(dynamic_cast<PetscLinearSolver<Number>*>(systemPProj.linear_solver.get()))->get_initial_residual()<<std::endl;
     std::cout<<"final linear residual = "<<systemPProj.final_linear_residual()<<std::endl;
  }

  for (; step<n_timesteps; ++step)
  {
    const Real dt = es.parameters.get<Real>("dt");
    time += dt;
    
    es.parameters.set<Real>("time") = time;
    es.parameters.set<unsigned int>("step") = step;

    std::cout<<"################ Solving time step "<<step<<" of "<<n_timesteps<<", time is "<<time<<std::endl;

    *systemAdvDiff.older_local_solution = *systemAdvDiff.old_local_solution;
    *systemAdvDiff.old_local_solution = *systemAdvDiff.current_local_solution;
   
    if (step > n_washin_steps)
      n_nonlinear_steps = max_nonlinear_steps;
    else
      n_nonlinear_steps = 1;
 
    for (unsigned int l=0; l<n_nonlinear_steps; ++l)
    {
      AutoPtr<NumericVector<Number> > last_nonlinear_soln (systemAdvDiff.solution->clone());
      last_nonlinear_soln->zero();
      last_nonlinear_soln->add(*systemAdvDiff.solution);
      
      systemAdvDiff.solve();

      last_nonlinear_soln->add (-1., *systemAdvDiff.solution);

      last_nonlinear_soln->close();
      const Real norm_delta = last_nonlinear_soln->l2_norm();
      const unsigned int n_linear_iterations = systemAdvDiff.n_linear_iterations();
      const Real final_linear_residual = systemAdvDiff.final_linear_residual();
      
      std::cout<<"number of linear iterations = "<<n_linear_iterations<<std::endl;
      std::cout<<"initial linear residual = "<<(dynamic_cast<PetscLinearSolver<Number>*>(systemAdvDiff.linear_solver.get()))->get_initial_residual()<<std::endl;
      std::cout<<"final linear residual = "<<final_linear_residual<<std::endl;
      std::cout<<"Nonlinear convergence: ||u - u_old|| = "<< norm_delta<< std::endl<<std::endl;

      if ((norm_delta < nonlinear_tolerance) &&
          (systemAdvDiff.final_linear_residual() < nonlinear_tolerance))
        {
          std::cout << " Nonlinear solver converged at step "
                    << l
                    << std::endl;
          break;
        }
    } // end nonlinear loop
    
    *systemPProj.older_local_solution = *systemPProj.old_local_solution;
    *systemPProj.old_local_solution = *systemPProj.current_local_solution;
    systemPProj.solve();
    
    std::cout<<"number of linear iterations = "<<systemPProj.n_linear_iterations()<<std::endl;
    std::cout<<"initial linear residual = "<<(dynamic_cast<PetscLinearSolver<Number>*>(systemPProj.linear_solver.get()))->get_initial_residual()<<std::endl;
    std::cout<<"final linear residual = "<<systemPProj.final_linear_residual()<<std::endl;

    if (step%write_solution_interval == 0)
    {
      std::vector<Real> local_solution;
      if (write_discontinuous)
      {
        vtk_io.build_discontinuous_solution_vector(mesh, es, local_solution);
        vtk_io.write_ascii_discontinuous(working_directory, step, mesh, local_solution);
      }
      if (write_continuous) 
      {
        vtk_io.build_continuous_solution_vector(mesh, es, local_solution);
        vtk_io.write_ascii_continuous(working_directory, step, mesh, local_solution);
      }
    }
    
    if (step%write_es_interval == 0)
    {
      es.allgather();    
      char filenamees[1024];
      sprintf(filenamees,"gnuid_es_%06d.xdr",step);
      es.write(working_directory + "/" + filenamees, libMeshEnums::ENCODE, 1 | 2);
      char filenamemesh[1024];
      sprintf(filenamemesh,"gnuid_mesh_%06d.xdr",step);
      mesh.write(working_directory + "/" + filenamemesh);
    }
    std::cout<<"################ End of time step"<<std::endl<<std::endl;
  }
}

void NS_DG_Solver::assemble_adv_diff(EquationSystems& es, const std::string& system_name)
{
  std::cout<<"assembling adv diff system... ";
  std::cout.flush();

  libmesh_assert (system_name == "AdvDiff");
  
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  GnuidBCHelper bcHelper;
  
  TransientLinearImplicitSystem & systemAdvDiff = es.get_system<TransientLinearImplicitSystem> ("AdvDiff");
  const unsigned int u_var = systemAdvDiff.variable_number ("u");
  const unsigned int v_var = systemAdvDiff.variable_number ("v");
  const unsigned int w_var = systemAdvDiff.variable_number ("w");
  TransientLinearImplicitSystem & systemPProj = es.get_system<TransientLinearImplicitSystem> ("PProj");
  const unsigned int p_var = systemPProj.variable_number ("p");

  const Real penalty = es.parameters.get<Real> ("penalty");
  const Real density = es.parameters.get<Real>("density");
  const Real viscosity = es.parameters.get<Real>("viscosity");
  const Real t = es.parameters.get<Real>("time");
  const Real step = es.parameters.get<unsigned int>("step");
  const Real dt = es.parameters.get<Real>("dt");
  const unsigned int n_washin_steps = es.parameters.get<unsigned int>("n washin steps");
  const unsigned int n_second_order_start_step = es.parameters.get<unsigned int>("second order start step");

  FEType fe_type = systemAdvDiff.variable_type(u_var);
  FEType fe_type_p = systemPProj.variable_type(p_var);

  AutoPtr<FEBase> fe  (FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_p  (FEBase::build(dim, fe_type_p));

#ifdef QORDER 
  QGauss qrule (dim, QORDER);
#else
  QGauss qrule (dim, fe_type.default_quadrature_order());
#endif
  fe->attach_quadrature_rule (&qrule);
  fe_p->attach_quadrature_rule (&qrule);
  
#ifdef QORDER 
  QGauss qface(dim-1, QORDER);
#else
  QGauss qface(dim-1, fe_type.default_quadrature_order());
#endif
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >& psi = fe_p->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_p->get_dphi();

  const std::vector<std::vector<Real> >&  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe_elem_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point>& qface_normals = fe_elem_face->get_normals();
  const std::vector<Point >& qface_points = fe_elem_face->get_xyz();
  std::vector<Point> qface_neighbor_point;
    
  const std::vector<std::vector<Real> >&  phi_neighbor_face = fe_neighbor_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_neighbor_face = fe_neighbor_face->get_dphi();
  
  const DofMap & dof_map = systemAdvDiff.get_dof_map();
  const DofMap & dof_map_p = systemPProj.get_dof_map();

  DenseMatrix<Real> Ke_uu, Ke_uu_n, Ke_vv_n, Ke_ww_n, Ke_uv_n, Ke_uw_n, Ke_vu_n, Ke_vw_n, Ke_wu_n, Ke_wv_n;
  DenseVector<Real> Fe_u, Fe_v, Fe_w, Fee_u, Fnn_u, Fee_v, Fnn_v, Fee_w, Fnn_w;
  
  DenseMatrix<Real> Kne_uu, Kne_uu_n, Kne_vv_n, Kne_ww_n, Kne_uv_n, Kne_uw_n, Kne_vu_n, Kne_vw_n, Kne_wu_n, Kne_wv_n;
  DenseMatrix<Real> Ken_uu, Ken_uu_n, Ken_vv_n, Ken_ww_n, Ken_uv_n, Ken_uw_n, Ken_vu_n, Ken_vw_n, Ken_wu_n, Ken_wv_n;
  DenseMatrix<Real> Kee_uu, Kee_uu_n, Kee_vv_n, Kee_ww_n, Kee_uv_n, Kee_uw_n, Kee_vu_n, Kee_vw_n, Kee_wu_n, Kee_wv_n;
  DenseMatrix<Real> Knn_uu, Knn_uu_n, Knn_vv_n, Knn_ww_n, Knn_uv_n, Knn_uw_n, Knn_vu_n, Knn_vw_n, Knn_wu_n, Knn_wv_n;
  
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;
  std::vector<unsigned int> neighbor_dof_indices_u;
  std::vector<unsigned int> neighbor_dof_indices_v;
  std::vector<unsigned int> neighbor_dof_indices_w;
  std::vector<unsigned int> neighbor_dof_indices_p;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    dof_map.dof_indices (elem, dof_indices_w, w_var);
    dof_map_p.dof_indices (elem, dof_indices_p, p_var);
    const unsigned int n_uvwp_e_dofs = dof_indices_u.size();
    const unsigned int n_p_e_dofs = dof_indices_p.size();

    fe->reinit(elem);
    fe_p->reinit(elem);
    Real elemVolume = 0.;
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      elemVolume += JxW[qp];
    }
    
    Ke_uu.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_uu_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_vv_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_ww_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_uv_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_uw_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_vu_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_vw_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_wu_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
    Ke_wv_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);

    Fe_u.resize(n_uvwp_e_dofs);
    Fe_v.resize(n_uvwp_e_dofs);
    Fe_w.resize(n_uvwp_e_dofs);
 
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      RealVectorValue grad_p_current (0.,0.,0.);
      RealVectorValue grad_p_old (0.,0.,0.);
      RealVectorValue grad_p_older (0.,0.,0.);
      for (unsigned int i=0; i<n_p_e_dofs; i++)
      {
        grad_p_current(0) += systemPProj.current_solution(dof_indices_p[i]) * dpsi[i][qp](0);
        grad_p_current(1) += systemPProj.current_solution(dof_indices_p[i]) * dpsi[i][qp](1);
        grad_p_current(2) += systemPProj.current_solution(dof_indices_p[i]) * dpsi[i][qp](2);
        grad_p_old(0) += systemPProj.old_solution(dof_indices_p[i]) * dpsi[i][qp](0);
        grad_p_old(1) += systemPProj.old_solution(dof_indices_p[i]) * dpsi[i][qp](1);
        grad_p_old(2) += systemPProj.old_solution(dof_indices_p[i]) * dpsi[i][qp](2);
        grad_p_older(0) += systemPProj.older_solution(dof_indices_p[i]) * dpsi[i][qp](0);
        grad_p_older(1) += systemPProj.older_solution(dof_indices_p[i]) * dpsi[i][qp](1);
        grad_p_older(2) += systemPProj.older_solution(dof_indices_p[i]) * dpsi[i][qp](2);
      }
      Real divU = 0.0;
      for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
      {
        RealVectorValue U( systemAdvDiff.current_solution(dof_indices_u[i]),
                           systemAdvDiff.current_solution(dof_indices_v[i]),
                           systemAdvDiff.current_solution(dof_indices_w[i]));
        divU += U * dphi[i][qp];
      }
      
      Real u = 0., u_old = 0.,u_older = 0.;
      Real v = 0., v_old = 0.,v_older = 0.;
      Real w = 0., w_old = 0.,w_older = 0.;
      Real du_dx = 0., du_dy = 0., du_dz = 0.;
      Real dv_dx = 0., dv_dy = 0., dv_dz = 0.;
      Real dw_dx = 0., dw_dy = 0., dw_dz = 0.;
      for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
      {
         u += phi[i][qp] * systemAdvDiff.current_solution(dof_indices_u[i]);
         v += phi[i][qp] * systemAdvDiff.current_solution(dof_indices_v[i]);
         w += phi[i][qp] * systemAdvDiff.current_solution(dof_indices_w[i]);
         u_old += phi[i][qp] * systemAdvDiff.old_solution(dof_indices_u[i]);
         v_old += phi[i][qp] * systemAdvDiff.old_solution(dof_indices_v[i]);
         w_old += phi[i][qp] * systemAdvDiff.old_solution(dof_indices_w[i]);
         u_older += phi[i][qp] * systemAdvDiff.older_solution(dof_indices_u[i]);
         v_older += phi[i][qp] * systemAdvDiff.older_solution(dof_indices_v[i]);
         w_older += phi[i][qp] * systemAdvDiff.older_solution(dof_indices_w[i]);
         du_dx += dphi[i][qp](0) * systemAdvDiff.current_solution(dof_indices_u[i]);
         du_dy += dphi[i][qp](1) * systemAdvDiff.current_solution(dof_indices_u[i]);
         du_dz += dphi[i][qp](2) * systemAdvDiff.current_solution(dof_indices_u[i]);
         dv_dx += dphi[i][qp](0) * systemAdvDiff.current_solution(dof_indices_v[i]);
         dv_dy += dphi[i][qp](1) * systemAdvDiff.current_solution(dof_indices_v[i]);
         dv_dz += dphi[i][qp](2) * systemAdvDiff.current_solution(dof_indices_v[i]);
         dw_dx += dphi[i][qp](0) * systemAdvDiff.current_solution(dof_indices_w[i]);
         dw_dy += dphi[i][qp](1) * systemAdvDiff.current_solution(dof_indices_w[i]);
         dw_dz += dphi[i][qp](2) * systemAdvDiff.current_solution(dof_indices_w[i]);
      }
      const RealVectorValue U (u, v, w);
      for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
      {
        for (unsigned int j=0; j<n_uvwp_e_dofs; j++)
        {
           Ke_uu(i,j) += JxW[qp]*(viscosity/density)*(dphi[i][qp]*dphi[j][qp]);
           Ke_uu(i,j) += JxW[qp]*phi[i][qp]*(U*dphi[j][qp]);
           Ke_uu(i,j) += JxW[qp]*0.5*(divU)*phi[i][qp]*phi[j][qp];
           if (step > n_washin_steps)
           {
             Ke_uu_n(i,j) += JxW[qp] * ((du_dx * phi[i][qp]*phi[j][qp]) + (0.5 * u * phi[i][qp]*dphi[j][qp](0)));
             Ke_uv_n(i,j) += JxW[qp] * ((du_dy * phi[i][qp]*phi[j][qp]) + (0.5 * u * phi[i][qp]*dphi[j][qp](1)));
             Ke_uw_n(i,j) += JxW[qp] * ((du_dz * phi[i][qp]*phi[j][qp]) + (0.5 * u * phi[i][qp]*dphi[j][qp](2)));
             Ke_vu_n(i,j) += JxW[qp] * ((dv_dx * phi[i][qp]*phi[j][qp]) + (0.5 * v * phi[i][qp]*dphi[j][qp](0)));
             Ke_vv_n(i,j) += JxW[qp] * ((dv_dy * phi[i][qp]*phi[j][qp]) + (0.5 * v * phi[i][qp]*dphi[j][qp](1)));
             Ke_vw_n(i,j) += JxW[qp] * ((dv_dz * phi[i][qp]*phi[j][qp]) + (0.5 * v * phi[i][qp]*dphi[j][qp](2)));
             Ke_wu_n(i,j) += JxW[qp] * ((dw_dx * phi[i][qp]*phi[j][qp]) + (0.5 * w * phi[i][qp]*dphi[j][qp](0)));
             Ke_wv_n(i,j) += JxW[qp] * ((dw_dy * phi[i][qp]*phi[j][qp]) + (0.5 * w * phi[i][qp]*dphi[j][qp](1)));
             Ke_ww_n(i,j) += JxW[qp] * ((dw_dz * phi[i][qp]*phi[j][qp]) + (0.5 * w * phi[i][qp]*dphi[j][qp](2)));
           }
           if (step > n_second_order_start_step)
           {
             Ke_uu(i,j) += JxW[qp]*(1.5/dt)*phi[i][qp]*phi[j][qp];
           }
           else
           {
             Ke_uu(i,j) += JxW[qp]*(1./dt)*phi[i][qp]*phi[j][qp];
           }
        }  
        if (step > n_second_order_start_step)
        {
          Fe_u(i) += JxW[qp]*(1./dt)*phi[i][qp]*(2 * u_old - 0.5 * u_older);
          Fe_v(i) += JxW[qp]*(1./dt)*phi[i][qp]*(2 * v_old - 0.5 * v_older);
          Fe_w(i) += JxW[qp]*(1./dt)*phi[i][qp]*(2 * w_old - 0.5 * w_older);
          Fe_u(i) -= JxW[qp] * phi[i][qp] * (1./density) * ((7./3.) * grad_p_current(0) - (5./3.) * grad_p_old(0) + (1./3.) * grad_p_older(0)); 
          Fe_v(i) -= JxW[qp] * phi[i][qp] * (1./density) * ((7./3.) * grad_p_current(1) - (5./3.) * grad_p_old(1) + (1./3.) * grad_p_older(1)); 
          Fe_w(i) -= JxW[qp] * phi[i][qp] * (1./density) * ((7./3.) * grad_p_current(2) - (5./3.) * grad_p_old(2) + (1./3.) * grad_p_older(2)); 
        }
        else
        {
          Fe_u(i) += JxW[qp]*(1./dt)*phi[i][qp]*(u_old);
          Fe_v(i) += JxW[qp]*(1./dt)*phi[i][qp]*(v_old);
          Fe_w(i) += JxW[qp]*(1./dt)*phi[i][qp]*(w_old);
          Fe_u(i) -= JxW[qp] * phi[i][qp] * (1./density) * (2. * grad_p_current(0) - grad_p_old(0)); 
          Fe_v(i) -= JxW[qp] * phi[i][qp] * (1./density) * (2. * grad_p_current(1) - grad_p_old(1)); 
          Fe_w(i) -= JxW[qp] * phi[i][qp] * (1./density) * (2. * grad_p_current(2) - grad_p_old(2)); 
        }
        if (step > n_washin_steps)
        {
          Fe_u(i) += JxW[qp] * ((u*du_dx+v*du_dy+w*du_dz) + (0.5*(du_dx+dv_dy+dw_dz)*u)) * phi[i][qp];
          Fe_v(i) += JxW[qp] * ((u*dv_dx+v*dv_dy+w*dv_dz) + (0.5*(du_dx+dv_dy+dw_dz)*v)) * phi[i][qp];
          Fe_w(i) += JxW[qp] * ((u*dw_dx+v*dw_dy+w*dw_dz) + (0.5*(du_dx+dv_dy+dw_dz)*w)) * phi[i][qp];
        }
      }
    }
    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      if (elem->neighbor(side) == NULL)
      {
        const unsigned short boundary_id = mesh.boundary_info->boundary_id(elem,side);
	std::string bc_type;
	int lid;
        bcHelper.bc_type_and_lid(es, boundary_id, bc_type, lid);
	int dirichletProfile = bc_type.compare("dirichlet_profile");
	int dirichletWall    = bc_type.compare("dirichlet_wall");
        if (dirichletProfile == 0 || dirichletWall == 0)
	{
          RealVectorValue U_bc;
          if (dirichletProfile == 0)
	    bcHelper.init_dirichletprofile_bc(es, lid);
	  else if (dirichletWall == 0)
	    U_bc.zero();
          
          fe_elem_face->reinit(elem, side);
          Real sideVolume = 0.;
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            sideVolume += JxW_face[qp];
          }
          const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
          double h = elemVolume/(sideVolume*elem->n_sides()*pow(elem_b_order,2.));
          
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
	    if (dirichletProfile == 0)
              bcHelper.compute_dirichletprofile_bc(qface_points[qp],t,U_bc);
            const Real u_bc_n = U_bc(0) * qface_normals[qp](0) + U_bc(1) * qface_normals[qp](1) + U_bc(2) * qface_normals[qp](2);
            Real u = 0.;
            Real v = 0.;
            Real w = 0.;
            for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
            {
               u += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_u[i]);
               v += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_v[i]);
               w += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_w[i]);
            }
            const Real u_n = u * qface_normals[qp](0) + v * qface_normals[qp](1) + w * qface_normals[qp](2);
            for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
            {
               for (unsigned int j=0; j<n_uvwp_e_dofs; j++)
               { 
                   Ke_uu(i,j) -= JxW_face[qp] * (viscosity/density) * (phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));
                   Ke_uu(i,j) -= JxW_face[qp] * (viscosity/density) * (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]));
                   Ke_uu(i,j) += JxW_face[qp] * (viscosity/density) * (penalty/h) * (phi_face[i][qp]*phi_face[j][qp]);
                   Ke_uu(i,j) -= JxW_face[qp] * 0.5 * u_n * (phi_face[i][qp] * phi_face[j][qp]);
                   if (step > n_washin_steps)
                   {
                     Ke_uu_n(i,j) -= JxW_face[qp] * (0.5 * u * qface_normals[qp](0) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_uv_n(i,j) -= JxW_face[qp] * (0.5 * u * qface_normals[qp](1) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_uw_n(i,j) -= JxW_face[qp] * (0.5 * u * qface_normals[qp](2) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_vu_n(i,j) -= JxW_face[qp] * (0.5 * v * qface_normals[qp](0) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_vv_n(i,j) -= JxW_face[qp] * (0.5 * v * qface_normals[qp](1) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_vw_n(i,j) -= JxW_face[qp] * (0.5 * v * qface_normals[qp](2) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_wu_n(i,j) -= JxW_face[qp] * (0.5 * w * qface_normals[qp](0) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_wv_n(i,j) -= JxW_face[qp] * (0.5 * w * qface_normals[qp](1) * phi_face[i][qp]*phi_face[j][qp]);
                     Ke_ww_n(i,j) -= JxW_face[qp] * (0.5 * w * qface_normals[qp](2) * phi_face[i][qp]*phi_face[j][qp]);
                   }
               }
               Fe_u(i) -= JxW_face[qp] * (viscosity/density) * (U_bc(0)*(qface_normals[qp]*dphi_face[i][qp]));
               Fe_v(i) -= JxW_face[qp] * (viscosity/density) * (U_bc(1)*(qface_normals[qp]*dphi_face[i][qp]));
               Fe_w(i) -= JxW_face[qp] * (viscosity/density) * (U_bc(2)*(qface_normals[qp]*dphi_face[i][qp]));

               Fe_u(i) += JxW_face[qp] * (viscosity/density) * (penalty/h) * phi_face[i][qp]*U_bc(0);
               Fe_v(i) += JxW_face[qp] * (viscosity/density) * (penalty/h) * phi_face[i][qp]*U_bc(1);
               Fe_w(i) += JxW_face[qp] * (viscosity/density) * (penalty/h) * phi_face[i][qp]*U_bc(2);

               Fe_u(i) -= JxW_face[qp] * 0.5 * u_bc_n * phi_face[i][qp] * U_bc(0);
               Fe_v(i) -= JxW_face[qp] * 0.5 * u_bc_n * phi_face[i][qp] * U_bc(1);
               Fe_w(i) -= JxW_face[qp] * 0.5 * u_bc_n * phi_face[i][qp] * U_bc(2);
               if (step > n_washin_steps)
               {
                 Fe_u(i) -= JxW_face[qp] * 0.5 * u_n * phi_face[i][qp] * u;
                 Fe_v(i) -= JxW_face[qp] * 0.5 * u_n * phi_face[i][qp] * v;
                 Fe_w(i) -= JxW_face[qp] * 0.5 * u_n * phi_face[i][qp] * w;
               }
            }          
          }
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
	  unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
          fe_neighbor_face->side_map (neighbor, elem_side.get(), side_neighbor, qface.get_points(), qface_neighbor_point);
          fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

          const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
          const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
          const double side_order_float = (elem_b_order + neighbor_b_order)/2.;
          Real sideVolume = 0.;
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            sideVolume += JxW_face[qp];
          }
          double h = (2.*elemVolume)/(sideVolume*elem->n_sides()*pow(side_order_float,2.));
          dof_map.dof_indices (neighbor, neighbor_dof_indices_u, u_var);
          dof_map.dof_indices (neighbor, neighbor_dof_indices_v, v_var);
          dof_map.dof_indices (neighbor, neighbor_dof_indices_w, w_var);
          const unsigned int n_uvwp_n_dofs = neighbor_dof_indices_u.size();
   
          Kne_uu.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_uu_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_uv_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_uw_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_vu_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_vv_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_vw_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_wu_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_wv_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Kne_ww_n.resize(n_uvwp_n_dofs, n_uvwp_e_dofs);
          Ken_uu.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_uu_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_uv_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_uw_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_vu_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_vv_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_vw_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_wu_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_wv_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Ken_ww_n.resize(n_uvwp_e_dofs, n_uvwp_n_dofs);
          Kee_uu.resize(n_uvwp_e_dofs, n_uvwp_e_dofs); 
          Kee_uu_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_uv_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_uw_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_vu_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_vv_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_vw_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_wu_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_wv_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Kee_ww_n.resize(n_uvwp_e_dofs, n_uvwp_e_dofs);
          Knn_uu.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_uu_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_uv_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_uw_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_vu_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_vv_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_vw_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_wu_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_wv_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          Knn_ww_n.resize(n_uvwp_n_dofs, n_uvwp_n_dofs);
          
          Fee_u.resize(n_uvwp_e_dofs);
          Fee_v.resize(n_uvwp_e_dofs);
          Fee_w.resize(n_uvwp_e_dofs);
          Fnn_u.resize(n_uvwp_n_dofs);
          Fnn_v.resize(n_uvwp_n_dofs);
          Fnn_w.resize(n_uvwp_n_dofs);
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
             Real u = 0.;
             Real v = 0.;
             Real w = 0.;
             for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
             {
                u += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_u[i]);
                v += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_v[i]);
                w += phi_face[i][qp] * systemAdvDiff.current_solution(dof_indices_w[i]);
             }
            Real u_neighbor = 0.;
            Real v_neighbor = 0.;
            Real w_neighbor = 0.;
            for (unsigned int i=0; i<n_uvwp_e_dofs; i++)
            {
               u_neighbor += phi_neighbor_face[i][qp] * systemAdvDiff.current_solution(neighbor_dof_indices_u[i]);
               v_neighbor += phi_neighbor_face[i][qp] * systemAdvDiff.current_solution(neighbor_dof_indices_v[i]);
               w_neighbor += phi_neighbor_face[i][qp] * systemAdvDiff.current_solution(neighbor_dof_indices_w[i]);
            }
            const Real u_n = u * qface_normals[qp](0) + v * qface_normals[qp](1) + w * qface_normals[qp](2);
            const Real u_n_neighbor = u_neighbor * qface_normals[qp](0) + v_neighbor * qface_normals[qp](1) + w_neighbor * qface_normals[qp](2);
            const Real u_n_mean = 0.5 * (u_n + u_n_neighbor);
            const Real u_n_jump = u_n - u_n_neighbor;

            const Real u_jump = u - u_neighbor;
            const Real v_jump = v - v_neighbor;
            const Real w_jump = w - w_neighbor;
            const Real u_mean = 0.5 * (u + u_neighbor);
            const Real v_mean = 0.5 * (v + v_neighbor);
            const Real w_mean = 0.5 * (w + w_neighbor);

            for (unsigned int i=0; i<n_uvwp_e_dofs; i++)          // elem elem matrix
            {
              for (unsigned int j=0; j<n_uvwp_e_dofs; j++)
              {
                 Kee_uu(i,j) -= 0.5 * JxW_face[qp] * (viscosity/density) * (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]));
                 Kee_uu(i,j) -= 0.5 * JxW_face[qp] * (viscosity/density) * (phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));
                 Kee_uu(i,j) += JxW_face[qp] * (viscosity/density) * (penalty/h) * (phi_face[j][qp]*phi_face[i][qp]);
                 Kee_uu(i,j) -= 0.5 * JxW_face[qp] * u_n * (phi_face[j][qp]*phi_face[i][qp]);
                 if (step > n_washin_steps)
                 {
                   Kee_uu_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](0) * (2*u - u_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_uv_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](1) * (2*u - u_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_uw_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](2) * (2*u - u_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_vu_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](0) * (2*v - v_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_vv_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](1) * (2*v - v_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_vw_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](2) * (2*v - v_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_wu_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](0) * (2*w - w_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_wv_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](1) * (2*w - w_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                   Kee_ww_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](2) * (2*w - w_neighbor) * (phi_face[j][qp]*phi_face[i][qp]);
                 }
              }
              if (step > n_washin_steps)
              { 
                Fee_u(i) -= 0.5 * JxW_face[qp] * (u_n_mean * u_jump + 0.5 * u_n_jump * u_mean + 0.25 * u_n_jump * u_jump) * phi_face[i][qp];
                Fee_v(i) -= 0.5 * JxW_face[qp] * (u_n_mean * v_jump + 0.5 * u_n_jump * v_mean + 0.25 * u_n_jump * v_jump) * phi_face[i][qp];
                Fee_w(i) -= 0.5 * JxW_face[qp] * (u_n_mean * w_jump + 0.5 * u_n_jump * w_mean + 0.25 * u_n_jump * w_jump) * phi_face[i][qp];
              }
            }

            for (unsigned int i=0; i<n_uvwp_n_dofs; i++) // neighbor neighbor matrix
            {
              for (unsigned int j=0; j<n_uvwp_n_dofs; j++)
              {
                 Knn_uu(i,j) += 0.5 * JxW_face[qp] * (viscosity/density) * (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));
                 Knn_uu(i,j) += 0.5 * JxW_face[qp] * (viscosity/density) * (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                 Knn_uu(i,j) += JxW_face[qp] * (viscosity/density) * (penalty/h) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                 Knn_uu(i,j) += 0.5 * JxW_face[qp] * u_n_neighbor * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                 if (step > n_washin_steps)
                 {
                   Knn_uu_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](0) * (2*u_neighbor - u) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_uv_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](1) * (2*u_neighbor - u) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_uw_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](2) * (2*u_neighbor - u) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_vu_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](0) * (2*v_neighbor - v) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_vv_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](1) * (2*v_neighbor - v) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_vw_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](2) * (2*v_neighbor - v) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_wu_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](0) * (2*w_neighbor - w) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_wv_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](1) * (2*w_neighbor - w) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                   Knn_ww_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](2) * (2*w_neighbor - w) * (phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]);
                 }
              } 
              if (step > n_washin_steps)
              {
                Fnn_u(i) -= 0.5 * JxW_face[qp] * (u_n_mean * u_jump + 0.5 * u_n_jump * u_mean - 0.25 * u_n_jump * u_jump) * phi_neighbor_face[i][qp];
                Fnn_v(i) -= 0.5 * JxW_face[qp] * (u_n_mean * v_jump + 0.5 * u_n_jump * v_mean - 0.25 * u_n_jump * v_jump) * phi_neighbor_face[i][qp];
                Fnn_w(i) -= 0.5 * JxW_face[qp] * (u_n_mean * w_jump + 0.5 * u_n_jump * w_mean - 0.25 * u_n_jump * w_jump) * phi_neighbor_face[i][qp];
              }
            }

            for (unsigned int i=0; i<n_uvwp_n_dofs; i++) // neighbor elem matrix
            {
              for (unsigned int j=0; j<n_uvwp_e_dofs; j++)
              {
                 Kne_uu(i,j) -= 0.5 * JxW_face[qp] * (viscosity/density) * (phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));
                 Kne_uu(i,j) += 0.5 * JxW_face[qp] * (viscosity/density) * (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));
                 Kne_uu(i,j) -= JxW_face[qp] * (viscosity/density) * (penalty/h) * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                 Kne_uu(i,j) -= 0.5 * JxW_face[qp] * u_n_mean * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                 if (step > n_washin_steps)
                 {
                   Kne_uu_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](0) * u * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_uv_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](1) * u * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_uw_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](2) * u * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_vu_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](0) * v * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_vv_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](1) * v * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_vw_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](2) * v * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_wu_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](0) * w * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_wv_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](1) * w * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                   Kne_ww_n(i,j) -= 0.25 * JxW_face[qp] * qface_normals[qp](2) * w * (phi_face[j][qp]*phi_neighbor_face[i][qp]);
                 }
              }
            }
            
            for (unsigned int i=0; i<n_uvwp_e_dofs; i++)         // elem neighbor
            {
              for (unsigned int j=0; j<n_uvwp_n_dofs; j++)
              {
                 Ken_uu(i,j) += 0.5 * JxW_face[qp] * (viscosity/density) * (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]));
                 Ken_uu(i,j) -= 0.5 * JxW_face[qp] * (viscosity/density) * (phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                 Ken_uu(i,j) -= JxW_face[qp] * (viscosity/density) * (penalty/h) * (phi_face[i][qp]*phi_neighbor_face[j][qp]);
                 Ken_uu(i,j) += 0.5 * JxW_face[qp] * u_n_mean * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                 if (step > n_washin_steps)
                 {
                   Ken_uu_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](0) * u_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_uv_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](1) * u_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_uw_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](2) * u_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_vu_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](0) * v_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_vv_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](1) * v_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_vw_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](2) * v_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_wu_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](0) * w_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_wv_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](1) * w_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                   Ken_ww_n(i,j) += 0.25 * JxW_face[qp] * qface_normals[qp](2) * w_neighbor * (phi_neighbor_face[j][qp]*phi_face[i][qp]);
                 }
              }
            }
          }
          systemAdvDiff.matrix->add_matrix(Kne_uu,neighbor_dof_indices_u,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Ken_uu,dof_indices_u,neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kee_uu,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Knn_uu,neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kne_uu,neighbor_dof_indices_v,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Ken_uu,dof_indices_v,neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kee_uu,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Knn_uu,neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kne_uu,neighbor_dof_indices_w,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Ken_uu,dof_indices_w,neighbor_dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Kee_uu,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Knn_uu,neighbor_dof_indices_w);
          
          systemAdvDiff.matrix->add_matrix(Kne_uu_n,neighbor_dof_indices_u,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Ken_uu_n,dof_indices_u,neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kee_uu_n,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Knn_uu_n,neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kne_vv_n,neighbor_dof_indices_v,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Ken_vv_n,dof_indices_v,neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kee_vv_n,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Knn_vv_n,neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kne_ww_n,neighbor_dof_indices_w,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Ken_ww_n,dof_indices_w,neighbor_dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Kee_ww_n,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Knn_ww_n,neighbor_dof_indices_w);

          systemAdvDiff.matrix->add_matrix(Kne_uv_n,neighbor_dof_indices_u,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Ken_uv_n,dof_indices_u,neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kee_uv_n,dof_indices_u, dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Knn_uv_n,neighbor_dof_indices_u, neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kne_uw_n,neighbor_dof_indices_u,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Ken_uw_n,dof_indices_u,neighbor_dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Kee_uw_n,dof_indices_u, dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Knn_uw_n,neighbor_dof_indices_u, neighbor_dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Kne_vu_n,neighbor_dof_indices_v,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Ken_vu_n,dof_indices_v,neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kee_vu_n,dof_indices_v,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Knn_vu_n,neighbor_dof_indices_v, neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kne_vw_n,neighbor_dof_indices_v,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Ken_vw_n,dof_indices_v,neighbor_dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Kee_vw_n,dof_indices_v,dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Knn_vw_n,neighbor_dof_indices_v, neighbor_dof_indices_w);
          systemAdvDiff.matrix->add_matrix(Kne_wu_n,neighbor_dof_indices_w,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Ken_wu_n,dof_indices_w,neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kee_wu_n,dof_indices_w,dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Knn_wu_n,neighbor_dof_indices_w, neighbor_dof_indices_u);
          systemAdvDiff.matrix->add_matrix(Kne_wv_n,neighbor_dof_indices_w,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Ken_wv_n,dof_indices_w,neighbor_dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Kee_wv_n,dof_indices_w,dof_indices_v);
          systemAdvDiff.matrix->add_matrix(Knn_wv_n,neighbor_dof_indices_w, neighbor_dof_indices_v);

          systemAdvDiff.rhs->add_vector(Fee_u, dof_indices_u);
          systemAdvDiff.rhs->add_vector(Fee_v, dof_indices_v);
          systemAdvDiff.rhs->add_vector(Fee_w, dof_indices_w);
          systemAdvDiff.rhs->add_vector(Fnn_u, neighbor_dof_indices_u);
          systemAdvDiff.rhs->add_vector(Fnn_v, neighbor_dof_indices_v);
          systemAdvDiff.rhs->add_vector(Fnn_w, neighbor_dof_indices_w);
        }
      }
    }
    systemAdvDiff.matrix->add_matrix(Ke_uu, dof_indices_u);
    systemAdvDiff.matrix->add_matrix(Ke_uu, dof_indices_v);
    systemAdvDiff.matrix->add_matrix(Ke_uu, dof_indices_w);
    systemAdvDiff.matrix->add_matrix(Ke_uu_n, dof_indices_u);
    systemAdvDiff.matrix->add_matrix(Ke_uv_n, dof_indices_u, dof_indices_v);
    systemAdvDiff.matrix->add_matrix(Ke_uw_n, dof_indices_u, dof_indices_w);
    systemAdvDiff.matrix->add_matrix(Ke_vv_n, dof_indices_v);
    systemAdvDiff.matrix->add_matrix(Ke_vu_n, dof_indices_v, dof_indices_u);
    systemAdvDiff.matrix->add_matrix(Ke_vw_n, dof_indices_v, dof_indices_w);
    systemAdvDiff.matrix->add_matrix(Ke_ww_n, dof_indices_w);
    systemAdvDiff.matrix->add_matrix(Ke_wu_n, dof_indices_w, dof_indices_u);
    systemAdvDiff.matrix->add_matrix(Ke_wv_n, dof_indices_w, dof_indices_v);
    systemAdvDiff.rhs->add_vector(Fe_u, dof_indices_u);
    systemAdvDiff.rhs->add_vector(Fe_v, dof_indices_v);
    systemAdvDiff.rhs->add_vector(Fe_w, dof_indices_w);
  }
  std::cout<<"done"<<std::endl;
}

void NS_DG_Solver::assemble_p_proj(EquationSystems& es, const std::string& system_name)
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
  for ( ; el != end_el; ++el)
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
	int lid;
        bcHelper.bc_type_and_lid(es, boundary_id, bc_type, lid);
	int dirichletProfile = bc_type.compare("dirichlet_profile");
	int dirichletWall    = bc_type.compare("dirichlet_wall");
        if (dirichletProfile  == 0 || dirichletWall == 0)
	{
          fe_elem_face->reinit(elem, side);
          fe_elem_face_p->reinit(elem, side);
          RealVectorValue U_bc;
          if (dirichletProfile == 0)
	    bcHelper.init_dirichletprofile_bc(es, lid);
	  else if (dirichletWall == 0)
	    U_bc.zero();
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            if (dirichletProfile == 0)
              bcHelper.compute_dirichletprofile_bc(qface_points[qp],t,U_bc);
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
              {
                Ke_p(i,j) += JxW_face[qp] * penalty_bc * psi_face[i][qp] * psi_face[j][qp];
              }
              Fe_p(i) += JxW_face[qp] * pressure_outlet * penalty_bc * psi_face[i][qp];
            }          
          }
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
            {
              Fee_p(i) += 0.5 * (density/dt) * JxW_face[qp] * u_n_jump * psi_face[i][qp];
            }

            for (unsigned int i=0; i<n_p_n_dofs; i++) // neighbor neighbor matrix
            {
              Fnn_p(i) += 0.5  * (density/dt) * JxW_face[qp] * u_n_jump * psi_neighbor_face[i][qp];
            }
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
