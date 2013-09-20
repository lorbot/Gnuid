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

#include "gnuid_bcHelper.h"
#include "libmesh/equation_systems.h"
#include "gnuid_cnHelper.h"

void GnuidBCHelper::init_bcData(const EquationSystems& es, const unsigned int boundary_id, std::string & bc_type)
{
  const unsigned short wall_id = es.parameters.get<unsigned short>("wall id");
  if (boundary_id == wall_id)
  {
    bc_type = "dirichlet_wall";
    return;
  }
  for (unsigned int i = 0; i < es.parameters.get<unsigned int>("n press bound"); i++)
  {
    char press_bc_id[1024];
    sprintf(press_bc_id,"press_bc_%02d_id",i);
    if (boundary_id == es.parameters.get<unsigned short>(press_bc_id))
    {
      bc_type = "neumann";
      return;
    }
  }
  for (unsigned int i = 0; i < es.parameters.get<unsigned int>("n vel bound"); i++)
  {
    char vel_bc_id[1024];
    sprintf(vel_bc_id,"vel_bc_%02d_id",i);
    if (boundary_id == es.parameters.get<unsigned short>(vel_bc_id))
    {
      init_dirichletIOProfile(es,i);
      bc_type = "dirichlet_profile";
      return;
    }
  }
  for (unsigned int i = 0; i < es.parameters.get<unsigned int>("n flow bound"); i++)
  {
    char flow_bc_id[1024];
    sprintf(flow_bc_id,"flow_bc_%02d_id",i);
    if (boundary_id == es.parameters.get<unsigned short>(flow_bc_id))
    {
      init_dirichletIODefective(es,i);
      bc_type = "dirichlet_defective";
      return;
    }
  }
  std::cerr<<"Fatal error: boundary ids do not match"<<std::endl;
  libmesh_error();
}

void GnuidBCHelper::init_dirichletIOProfile(const EquationSystems& es, const unsigned int& lid)
{
  _lid = lid;
  Real Re_number = es.parameters.get<Real>("Re number");
  Real mass_flow = es.parameters.get<Real>("mass flow");
  Real density = es.parameters.get<Real>("density");
  Real viscosity = es.parameters.get<Real>("viscosity");
  _t_period = es.parameters.get<Real>("t period");
  char vel_bc_type[1024];
  char vel_bc_center[1024];
  char vel_bc_radius[1024];
  char vel_bc_normal[1024];
  char vel_bc_scaling[1024];
  char vel_bc_n_fourier_modes[1024];
  sprintf(vel_bc_type,"vel_bc_%02d_type",lid);
  sprintf(vel_bc_center,"vel_bc_%02d_center",lid);
  sprintf(vel_bc_radius,"vel_bc_%02d_radius",lid);
  sprintf(vel_bc_normal,"vel_bc_%02d_normal",lid);
  sprintf(vel_bc_scaling,"vel_bc_%02d_scaling",lid);
  sprintf(vel_bc_n_fourier_modes,"vel_bc_%02d_n_fourier_modes",lid);
  _bc_Radius = es.parameters.get<Real>(vel_bc_radius);
  _bc_Type = es.parameters.get<std::string>(vel_bc_type);
  _bc_Normal = es.parameters.get<VectorValue<Real> >(vel_bc_normal);
  _bc_Center = es.parameters.get<Point>(vel_bc_center);
  _scaling = es.parameters.get<Real>(vel_bc_scaling);
  _alpha = _bc_Radius * sqrt((2.*libMesh::pi/_t_period)/(viscosity/density));
  unsigned int n_f_modes = es.parameters.get<unsigned int>(vel_bc_n_fourier_modes);
  _f_modes.resize(n_f_modes);
  for (unsigned int k =0; k < n_f_modes; k++)
  {
     char vel_bc_f_mode[1024];
     sprintf(vel_bc_f_mode,"vel_bc_%02d_f_mode_%02d",lid,k);
     _f_modes[k] = es.parameters.get<Complex>(vel_bc_f_mode);
  }
  _u_mean = 0.;
  if (Re_number != 0 && mass_flow == 0)
    _u_mean = _scaling * (Re_number * viscosity)/(density * 2.0 * _bc_Radius);
  else if (Re_number == 0 && mass_flow != 0)
    _u_mean = _scaling * mass_flow/(density * libMesh::pi * _bc_Radius * _bc_Radius); 
}
    
void GnuidBCHelper::compute_dirichletIOProfile(const Point& point, const Real& t, RealVectorValue& U_bc)
{
  const Real r_sq = (point - _bc_Center).size_sq();
  const Real R_sq = _bc_Radius * _bc_Radius;
  if (R_sq - r_sq > 0.0)
  {
    Real u = 2. * ((R_sq - r_sq)/R_sq) * _f_modes[0].real();
    for (unsigned int k =1; k < _f_modes.size(); k++)
    {
      Complex c_i = Complex(0,1);

      Complex c_a = C_op::RCmul(_alpha * sqrt(1.0*k), C_op::Cmul(c_i, sqrt(c_i)));  
      Complex c_a_2 = C_op::RCmul(sqrt(r_sq)/_bc_Radius, c_a);
                                                               
      Complex c_t = Complex(0, 2.*libMesh::pi*k*t/_t_period);
      Complex c_f10 = C_op::RCmul(2., C_op::Cdiv(C_op::Cbes(1,c_a), C_op::Cmul(C_op::Cbes(0,c_a),c_a)));
      Complex c_const = C_op::Cdiv(C_op::Cmul(_f_modes[k], C_op::Cexp(c_t)), Complex(1. - c_f10.real(), 0. - c_f10.imag()));
      Complex c_u = C_op::Cmul(c_const, Complex(1.0 - C_op::Cdiv(C_op::Cbes(0,c_a_2), C_op::Cbes(0,c_a)).real(),
                                                0.0 - C_op::Cdiv(C_op::Cbes(0,c_a_2), C_op::Cbes(0,c_a)).imag()));
      u += c_u.real();        
    }                                                            
    if (_bc_Type == "inlet")                               
      U_bc.assign((-1.) * _bc_Normal * u * _u_mean);
    else if (_bc_Type == "outlet")
      U_bc.assign(_bc_Normal * u * _u_mean);
    else
      libmesh_error();
  }
  else
    U_bc.zero();
}

void GnuidBCHelper::init_dirichletIODefective(const EquationSystems& es, const unsigned int& lid)
{
  _lid = lid;
  Real Re_number = es.parameters.get<Real>("Re number");
  Real mass_flow = es.parameters.get<Real>("mass flow");
  Real density = es.parameters.get<Real>("density");
  Real viscosity = es.parameters.get<Real>("viscosity");
  _t_period = es.parameters.get<Real>("t period");
  char flow_bc_type[1024];
  char flow_bc_center[1024];
  char flow_bc_radius[1024];
  char flow_bc_normal[1024];
  char flow_bc_scaling[1024];
  char flow_bc_n_fourier_modes[1024];
  sprintf(flow_bc_type,"flow_bc_%02d_type",lid);
  sprintf(flow_bc_center,"flow_bc_%02d_center",lid);
  sprintf(flow_bc_radius,"flow_bc_%02d_radius",lid);
  sprintf(flow_bc_normal,"flow_bc_%02d_normal",lid);
  sprintf(flow_bc_scaling,"flow_bc_%02d_scaling",lid);
  sprintf(flow_bc_n_fourier_modes,"flow_bc_%02d_n_fourier_modes",lid);
  _bc_Radius = es.parameters.get<Real>(flow_bc_radius);
  _bc_Type = es.parameters.get<std::string>(flow_bc_type);
  _bc_Normal = es.parameters.get<VectorValue<Real> >(flow_bc_normal);
  _bc_Center = es.parameters.get<Point>(flow_bc_center);
  _scaling = es.parameters.get<Real>(flow_bc_scaling);
  _alpha = _bc_Radius * sqrt((2.*libMesh::pi/_t_period)/(viscosity/density));
  unsigned int n_f_modes = es.parameters.get<unsigned int>(flow_bc_n_fourier_modes);
  _f_modes.resize(n_f_modes);
  for (unsigned int k =0; k < n_f_modes; k++)
  {
     char flow_bc_f_mode[1024];
     sprintf(flow_bc_f_mode,"flow_bc_%02d_f_mode_%02d",lid,k);
     _f_modes[k] = es.parameters.get<Complex>(flow_bc_f_mode);
  }
  _u_mean = 0.;
  if (Re_number != 0 && mass_flow == 0)
    _u_mean = _scaling * (Re_number * viscosity)/(density * 2.0 * _bc_Radius);
  else if (Re_number == 0 && mass_flow != 0)
    _u_mean = _scaling * mass_flow/(density * libMesh::pi * _bc_Radius * _bc_Radius); 
}

void GnuidBCHelper::compute_dirichletDefectiveData(const Real& t, RealVectorValue& U_bc, unsigned int& lid)
{
  Real u_mean_t = 0.;
  for (unsigned int k = 0; k < _f_modes.size(); k++)
  {
    Complex c_t = Complex(0, 2.*libMesh::pi*k*t/_t_period);
    u_mean_t += _u_mean * C_op::Cmul(_f_modes[k],C_op::Cexp(c_t)).real();
  }
  if (_bc_Type == "inlet")
    U_bc.assign((-1) * _bc_Normal * u_mean_t);
  else if (_bc_Type == "outlet")
    U_bc.assign(_bc_Normal * u_mean_t);
  lid = _lid;
}
