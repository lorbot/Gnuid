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

#ifndef __complex_operations__
#define __complex_operations__

#include "math.h"
#include "libmesh/libmesh.h"

namespace C_op 
{
  double Cabs(Complex z)
  {
     double x,y,ans,temp;
  
     x=fabs(z.real());
     y=fabs(z.imag());
     if (x == 0.0)
        ans=y;
     else if (y == 0.0)
        ans=x;
     else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
     } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
     }
    return (double) ans;
  }
  
  Complex Cadd(Complex a, Complex b)
  {
     double creal=a.real()+b.real();
     double cimag=a.imag()+b.imag();
     return Complex(creal,cimag);
  }
  
  Complex Cdiv(Complex a, Complex b)
  {
     double creal, cimag;
     double r, den;
  
     if (fabs(b.real()) >= fabs(b.imag())) 
     {
        r=b.imag()/b.real();
        den=b.real()+r*b.imag();
        creal=(a.real()+r*a.imag())/den;
        cimag=(a.imag()-r*a.real())/den;
     }
     else
     {
        r=b.real()/b.imag();
        den=b.imag()+r*b.real();
        creal=(a.real()*r+a.imag())/den;
        cimag=(a.imag()*r-a.real())/den;
     }
   return Complex(creal,cimag);
  }
  
  Complex Cexp(Complex a)
  {
     double im = exp(a.real());
     double re = exp(a.real());
  
     re *= cos(a.imag());
     im *= sin(a.imag());
  
     return Complex(re,im);
  }
  
  Complex Cmul(Complex a, Complex b)
  {
     double creal=a.real()*b.real()-a.imag()*b.imag();
     double cimag=a.imag()*b.real()+a.real()*b.imag();
     return Complex(creal,cimag);
  }
  
  Complex RCmul(double x, Complex a)
  {
     double creal=x*a.real();
     double cimag=x*a.imag();
     return Complex(creal,cimag);
  }

  Complex Csqrt(Complex z)
  {
     double x, y, w, r;
     double creal, cimag;  

     if ((z.real() == 0.0) && (z.imag() == 0.0))
     {
        return Complex(0.0,0.0);
     }
     else
     {
        x=fabs(z.real());
        y=fabs(z.imag());
        if (x >= y) 
        {
           r=y/x;
           w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } 
        else
        {
           r=x/y;
           w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.real() >= 0.0) 
        {
           creal=w;
           cimag=z.imag()/(2.0*w);
        }
        else
        {
           cimag=(z.imag() >= 0) ? w : -w;
           creal=z.imag()/(2.0*cimag);
        }
        return Complex(creal,cimag);
     }
  }

  Complex Cbes(short n, Complex argument)
  {
     Complex z = Complex(1.0, 0.0);
     Complex zarg;
     Complex zproduct = Complex(1.0, 0.0);
  
     zarg = RCmul(-0.25, Cmul(argument, argument));
     Complex zanswer(1.0,0.0);
  
     for (unsigned int i=1; i <= 10000; i++)
     {
       z = RCmul(1./(static_cast<double>(i) * static_cast<double>(i+n)), Cmul(z, zarg));
       if (Cabs(z) <= 1.e-20)
         break;
       zanswer = Cadd(zanswer,z);
     }
  
     for (int i=1; i <= n; zproduct = Cmul(zproduct, RCmul(0.5, argument)), i++);
  
     zanswer = Cmul(zanswer, zproduct);
  
     return zanswer;
  }
}

#endif
