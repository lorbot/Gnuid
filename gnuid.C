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

#include <iostream>
#include <algorithm>
#include <math.h>

#include "libmesh.h"

#include "mesh.h"
#include "mesh_data.h"
#include "point.h"
#include "perf_log.h"
#include "getpot.h"

#include "mesh_refinement.h"
#include "mesh_tools.h"
#include "equation_systems.h"
#include "ns_dg_solver.h"

using namespace libMesh;
int main (int argc, char** argv)
{
  LibMeshInit init (argc, argv);
  {    
    PerfLog perf_log("gnuid");

    std::string working_directory = ".";
    if (argc >= 2)
    {
      working_directory = argv[1];
    }

    GetPot input_file(working_directory + "/" + "gnuid.in");
    std::vector<std::string> identified_variables;

    const unsigned int dim = 3;     
    
    Mesh mesh (dim);

    std::string mesh_file;
    if (input_file.have_variable("mesh_file"))
    {
       mesh_file = input_file("mesh_file", "");
       identified_variables.push_back("mesh_file");
    }
    else std::cerr<<"Missing variable mesh_file. Check your input file."<<std::endl;

    mesh.read (working_directory + "/" + mesh_file);
    
    bool second_order_mesh = false;
    if (input_file.have_variable("first_to_second_order_mesh"))
    {
       second_order_mesh = input_file("first_to_second_order_mesh", false);
       identified_variables.push_back("first_to_second_order_mesh");
    }
    else std::cerr<<"Missing variable first_to_second_order_mesh, using default value. Check your input file"<<std::endl;
    if (second_order_mesh)
    {
       mesh.all_second_order();
    }
    
    mesh.print_info();
    EquationSystems equation_systems (mesh);
    NS_DG_Solver ns_dg_Solver(equation_systems);
    
    if (input_file.have_variable("penalty"))
    {
       equation_systems.parameters.set<Real>("penalty") = input_file("penalty", 4.); // velocity stabilization 
       identified_variables.push_back("penalty");
    }
    else std::cerr<<"Missing variable penalty, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("density"))
    {
       equation_systems.parameters.set<Real>("density") = input_file("density", 1.05); // g/cm^3
       identified_variables.push_back("density");
    }
    else std::cerr<<"Missing variable density, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("viscosity"))
    {
       equation_systems.parameters.set<Real>("viscosity") = input_file("viscosity", 0.035); // poise
       identified_variables.push_back("viscosity");                          
    }                                                                        
    else std::cerr<<"Missing variable viscosity, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("dt"))                                      
    {
       equation_systems.parameters.set<Real>("dt") = input_file("dt", 0.001); // s
       identified_variables.push_back("dt");
    }
    else std::cerr<<"Missing variable dt, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("time"))
    {
       equation_systems.parameters.set<Real>("time") = input_file("time", 0.0); // s
       identified_variables.push_back("time");
    }
    else std::cerr<<"Missing variable time, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("period"))
    {
       equation_systems.parameters.set<Real>("t period") = input_file("period", 0.0); //s
       identified_variables.push_back("period");
    }
    else std::cerr<<"Missing variable period, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("step"))
    {
       equation_systems.parameters.set<unsigned int>("step") = input_file("step", 0); 
       identified_variables.push_back("step");
    }
    else std::cerr<<"Missing variable step, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("Re_number"))
    {
       equation_systems.parameters.set<Real>("Re number") = input_file("Re_number", 0.0);
       identified_variables.push_back("Re_number");
    }
    else std::cerr<<"Missing variable Re_number, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("mass_flow"))
    {
       equation_systems.parameters.set<Real>("mass flow") = input_file("mass_flow" , 0.0);
       identified_variables.push_back("mass_flow");
    }
    else std::cerr<<"Missing variable mass_flow, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("n_uniform_refinement_steps"))
    {
       equation_systems.parameters.set<unsigned int>("n uniform refinement steps") = input_file("n_uniform_refinement_steps", 0); 
       identified_variables.push_back("n_uniform_refinement_steps");
    }
    else std::cerr<<"Missing variable n_uniform_refinement_steps, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("velocity_approximation_order"))
    {
       equation_systems.parameters.set<unsigned int>("velocity approx order") = input_file("velocity_approximation_order", 1); 
       identified_variables.push_back("velocity_approximation_order");
    }
    else std::cerr<<"Missing variable velocity_approximation_order, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("pressure_approximation_order"))
    {
       equation_systems.parameters.set<unsigned int>("pressure approx order") = input_file("pressure_approximation_order", 1); 
       identified_variables.push_back("pressure_approximation_order");
    }
    else std::cerr<<"Missing variable pressure_approximation_order, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("n_timesteps"))
    {
       equation_systems.parameters.set<unsigned int>("n timesteps") = input_file("n_timesteps", 1000);
       identified_variables.push_back("n_timesteps");
    }
    else std::cerr<<"Missing variable n_timesteps, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("second_order_start_step"))
    {
       equation_systems.parameters.set<unsigned int>("second order start step") = input_file("second_order_start_step", 1);
       identified_variables.push_back("second_order_start_step");                          
    }                                               
    else std::cerr<<"Missing variable second_order_start_step, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("n_washin_steps"))
    {
       equation_systems.parameters.set<unsigned int>("n washin steps") = input_file("n_washin_steps", 0);
       identified_variables.push_back("n_washin_steps");                          
    }                                               
    else std::cerr<<"Missing variable n_washin_steps, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("linear_solver_maximum_iterations"))
    {
       equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = input_file("linear_solver_maximum_iterations", 3000);
       identified_variables.push_back("linear_solver_maximum_iterations");
    }
    else std::cerr<<"Missing variable linear_solver_maximum_iterations, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("linear_solver_tolerance"))
    {
       equation_systems.parameters.set<Real>("linear solver tolerance") = input_file("linear_solver_tolerance", 1.e-10);
       identified_variables.push_back("linear_solver_tolerance");
    }
    else std::cerr<<"Missing variable linear_solver_tolerance, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("nonlinear_solver_maximum_iterations"))
    {
       equation_systems.parameters.set<unsigned int>("nonlinear solver maximum iterations") = input_file("nonlinear_solver_maximum_iterations", 10);
       identified_variables.push_back("nonlinear_solver_maximum_iterations");
    }
    else std::cerr<<"Missing variable nonlinear_solver_maximum_iterations, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("nonlinear_solver_tolerance"))
    {
       equation_systems.parameters.set<Real>("nonlinear solver tolerance") = input_file("nonlinear_solver_tolerance", 1.e-1);
       identified_variables.push_back("nonlinear_solver_tolerance");
    }
    else std::cerr<<"Missing variable nonlinear_solver_tolerance, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("n_steps_between_solution_output"))
    {
       equation_systems.parameters.set<unsigned int>("write output interval") = input_file("n_steps_between_solution_output",100);
       identified_variables.push_back("n_steps_between_solution_output");
    }
    else std::cerr<<"Missing variable n_steps_between_solution_output, using default value. Check your input file."<<std::endl;
    if (input_file.have_variable("n_steps_between_equation_system_output"))
    {
       equation_systems.parameters.set<unsigned int>("write es interval") = input_file("n_steps_between_equation_system_output",100);
       identified_variables.push_back("n_steps_between_equation_system_output");
    }
    else std::cerr<<"Missing variable n_steps_between_equation_system_output, using default value. Check your input file."<<std::endl;
    
    if (input_file.have_variable("wall_id"))
    {
        equation_systems.parameters.set<unsigned short>("wall id") = input_file("wall_id", 0);
        identified_variables.push_back("wall_id");
    }
    else std::cerr<<"Missing variable wall_id, using default value. Check your input file"<<std::endl;
    if (input_file.have_variable("n_velocity_boundaries"))
    {
        equation_systems.parameters.set<unsigned int>("n vel bound") = input_file("n_velocity_boundaries", 1);
        identified_variables.push_back("n_velocity_boundaries");
    }
    else std::cerr<<"Missing variable n_velocity_boundaries, using default value. Check your input file"<<std::endl;
    unsigned int n_vel_bound = equation_systems.parameters.get<unsigned int>("n vel bound");
    for (unsigned int i = 0; i < n_vel_bound; i++)
    {
        char vel_bc_type[1024];
        char vel_bc_name[1024];
        char vel_bc_name_output[1024];
        char vel_bc_id[1024];
        char vel_bc_normal[1024];
        char vel_bc_center[1024];
        char vel_bc_radius[1024];
        char vel_bc_scaling[1024];
        char vel_bc_n_fourier_modes[1024];
        sprintf(vel_bc_type,"vel_bc_%02d_type",i);
        sprintf(vel_bc_id,"vel_bc_%02d_id",i);
        sprintf(vel_bc_center,"vel_bc_%02d_center",i);
        sprintf(vel_bc_radius,"vel_bc_%02d_radius",i);
        sprintf(vel_bc_normal,"vel_bc_%02d_normal",i);
        sprintf(vel_bc_scaling,"vel_bc_%02d_scaling",i);
        sprintf(vel_bc_name,"vel_bc_%02d/",i);
        sprintf(vel_bc_name_output,"vel_bc_%02d",i);
        sprintf(vel_bc_n_fourier_modes,"vel_bc_%02d_n_fourier_modes",i);
        
        std::string prefix = std::string(vel_bc_name);
        input_file.set_prefix(vel_bc_name);
        if (input_file.have_variable("inlet_outlet"))
        {
           equation_systems.parameters.set<std::string>(vel_bc_type) = input_file("inlet_outlet", "inlet");
           identified_variables.push_back(prefix + "inlet_outlet");
        }
        else std::cerr<<"Missing variable inlet_outlet in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        if (input_file.have_variable("vel_bc_id"))
        {
           equation_systems.parameters.set<unsigned short>(vel_bc_id) = input_file("vel_bc_id", 0);
           identified_variables.push_back(prefix + "vel_bc_id");
        }
        else std::cerr<<"Missing variable vel_bc_id in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        if (input_file.have_variable("vel_bc_center_x") && input_file.have_variable("vel_bc_center_y") && input_file.have_variable("vel_bc_center_z"))
        {
           equation_systems.parameters.set<Point>(vel_bc_center) =
                                          Point(input_file("vel_bc_center_x",0.0),input_file("vel_bc_center_y",0.0),input_file("vel_bc_center_z",0.0));
           identified_variables.push_back(prefix + "vel_bc_center_x");
           identified_variables.push_back(prefix + "vel_bc_center_y");
           identified_variables.push_back(prefix + "vel_bc_center_z");
        }
        else std::cerr<<"Missing variable vel_bc_center in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        if (input_file.have_variable("vel_bc_radius"))
        {
           equation_systems.parameters.set<Real>(vel_bc_radius) = input_file("vel_bc_radius", 0.0);
           identified_variables.push_back(prefix + "vel_bc_radius");
        }
        else std::cerr<<"Missing variable vel_bc_radius in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        if (input_file.have_variable("vel_bc_normal_i") && input_file.have_variable("vel_bc_normal_j") && input_file.have_variable("vel_bc_normal_k"))
        {
           equation_systems.parameters.set<VectorValue<Real> >(vel_bc_normal) =
                                          VectorValue<Real>(input_file("vel_bc_normal_i",0.0),input_file("vel_bc_normal_j",0.0),input_file("vel_bc_normal_k",0.0)).unit();
           identified_variables.push_back(prefix + "vel_bc_normal_i");
           identified_variables.push_back(prefix + "vel_bc_normal_j");
           identified_variables.push_back(prefix + "vel_bc_normal_k");
        }
        else std::cerr<<"Missing variable vel_bc_normal in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        if (input_file.have_variable("scaling"))
        {
           equation_systems.parameters.set<Real>(vel_bc_scaling) = input_file("scaling", 1.0);
           identified_variables.push_back(prefix + "scaling");
        }
        else std::cerr<<"Missing variable scaling in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        if (input_file.have_variable("n_fourier_modes"))
        {
           equation_systems.parameters.set<unsigned int>(vel_bc_n_fourier_modes) = input_file("n_fourier_modes", 1);
           identified_variables.push_back(prefix + "n_fourier_modes");
        }
        else std::cerr<<"Missing variable n_fourier_modes in section ["<<vel_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        for (unsigned int n = 0; n < equation_systems.parameters.get<unsigned int>(vel_bc_n_fourier_modes); n++)
        {
            char vel_bc_f_mode[1024];
            char f_mode_re[1024];
            char f_mode_im[1024];
            sprintf(vel_bc_f_mode,"vel_bc_%02d_f_mode_%02d",i,n);
            sprintf(f_mode_re,"f_mode_%02d_re",n);
            sprintf(f_mode_im,"f_mode_%02d_im",n);
            std::string f_m_re = std::string(f_mode_re);
            std::string f_m_im = std::string(f_mode_im);
            if (input_file.have_variable(f_mode_re) && input_file.have_variable(f_mode_im))
            {
               equation_systems.parameters.set<Complex>(vel_bc_f_mode) = Complex(input_file(f_mode_re, 0.0),input_file(f_mode_im, 0.0));
               identified_variables.push_back(prefix + f_m_re);
               identified_variables.push_back(prefix + f_m_im);
            }
            else std::cerr<<"Missing variable "<<f_mode_re<<" or "<<f_mode_im<<" in section ["<<vel_bc_name_output
                          <<"], using default value. Check your input file sections and variables."<<std::endl;
        }
        input_file.set_prefix("");
    }
    if (input_file.have_variable("n_pressure_boundaries"))
    {
       equation_systems.parameters.set<unsigned int>("n press bound") = input_file("n_pressure_boundaries", 1);
       identified_variables.push_back("n_pressure_boundaries");
    }
    else std::cerr<<"Missing variable n_pressure_boundaries, using default value. Check your input file"<<std::endl;
    unsigned int n_press_bound = equation_systems.parameters.get<unsigned int>("n press bound");
    for (unsigned int i = 0; i < n_press_bound; i++)
    {
        char press_bc_name[1024];
        char press_bc_name_output[1024];
        char press_bc_id[1024];
        sprintf(press_bc_id,"press_bc_%02d_id",i);
        sprintf(press_bc_name,"press_bc_%02d/",i);
        sprintf(press_bc_name_output,"press_bc_%02d",i);
        
        std::string prefix = std::string(press_bc_name);
        input_file.set_prefix(press_bc_name);
        if (input_file.have_variable("press_bc_id"))
        {
           equation_systems.parameters.set<unsigned short>(press_bc_id) = input_file("press_bc_id", 0);
           identified_variables.push_back(prefix + "press_bc_id");
        }
        else std::cerr<<"Missing variable press_bc_id in section ["<<press_bc_name_output
                      <<"], using default value. Check your input file sections and variables."<<std::endl;
        input_file.set_prefix("");
    }
    
    bool read_es_file = false;
    if (input_file.have_variable("read_es_file"))
    {
       read_es_file = input_file("read_es_file", false);
       equation_systems.parameters.set<bool>("read es file") = read_es_file;
       identified_variables.push_back("read_es_file");
    }
    else std::cerr<<"Missing variable read_es_file, using default value. Check your input file"<<std::endl;
    if (read_es_file == true)
    {
      std::string es_file;
      if (input_file.have_variable("es_file"))
      {
         es_file = input_file("es_file", "");
         identified_variables.push_back("es_file");
      }
      else std::cerr<<"Missing variable es_file. Check your input file."<<std::endl;
      equation_systems.read(working_directory + "/" + es_file, libMeshEnums::DECODE, 1 | 2 | 4);    
    }
    else 
    {
      identified_variables.push_back("es_file");
    }
    std::vector<std::string> unidentified_variables = input_file.unidentified_variables(identified_variables);
    for (unsigned int num = 0; num < unidentified_variables.size(); num++)
    {
       std::cout<<"unidentified variable: "<<unidentified_variables[num]<<" Check your input file."<<std::endl;
    }
    input_file.print();
    
    ns_dg_Solver.init();
    ns_dg_Solver.working_directory = working_directory;

    perf_log.start_event("gnuid");
    ns_dg_Solver.solve();
    perf_log.stop_event("gnuid");
    
    std::cout<<equation_systems.parameters<<std::endl;
  }

  return 0;
}
