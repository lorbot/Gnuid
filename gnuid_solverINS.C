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
#include "gnuid_vtkWriter.h"
#include "gnuid_bcHelper.h"
#include "gnuid_spHelper.h"

//#define QORDER TENTH 

using namespace libMesh;

void GnuidSolver::init()
{
  EquationSystems& es = this->system();
  if (es.parameters.get<bool>("read es file"))
  {
     TransientLinearImplicitSystem & systemAdvDiff = es.get_system<TransientLinearImplicitSystem> ("AdvDiff");
     systemAdvDiff.update();
     systemAdvDiff.attach_assemble_function (GnuidSolver::assemble_adv_diff);
     DofMap & dof_map = systemAdvDiff.get_dof_map();
     MeshBase& mesh = es.get_mesh();
     AugmentSparsityPatternDefectiveBC augmenter(mesh, dof_map, es);
     dof_map.attach_extra_sparsity_object(augmenter);
     TransientLinearImplicitSystem & systemPProj = es.get_system<TransientLinearImplicitSystem> ("PProj");
     systemPProj.update();
     systemPProj.attach_assemble_function (GnuidSolver::assemble_p_proj);
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
     systemAdvDiff.attach_assemble_function (GnuidSolver::assemble_adv_diff);
     DofMap & dof_map = systemAdvDiff.get_dof_map();
     MeshBase& mesh = es.get_mesh();
     AugmentSparsityPatternDefectiveBC augmenter(mesh, dof_map, es);
     dof_map.attach_extra_sparsity_object(augmenter);
     TransientLinearImplicitSystem & systemPProj = es.add_system<TransientLinearImplicitSystem> ("PProj");
     systemPProj.add_variable ("p", p_libmesh, LAGRANGE);
     systemPProj.attach_assemble_function (GnuidSolver::assemble_p_proj);
     es.init();
  }
}

void GnuidSolver::solve()
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
  GnuidVTKWriter vtk_io;
  
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


