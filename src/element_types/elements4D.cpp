/*****************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and 
to permit others to do so.


This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
    
    1.  Redistributions of source code must retain the above copyright notice, this list of 
        conditions and the following disclaimer.
 
    2.  Redistributions in binary form must reproduce the above copyright notice, this list of 
        conditions and the following disclaimer in the documentation and/or other materials 
        provided with the distribution.
 
    3.  Neither the name of the copyright holder nor the names of its contributors may be used 
        to endorse or promote products derived from this software without specific prior 
        written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************************/

/**************************************************************************************
 *  Elements contains many of the basis functions required to implement a wide range of 
 *  numerical methods for computational physics. Currently 2D and 3D serendipity basis set 
 *  is included, as well as the arbitrary order 2D and 3D tensor product hexahedral elements.
 *  There is also a 4D tesseract element that can be used for space-time methods. Make sure 
 *  include the elements.h file in your code to access the elememnts library.
 *****************************************************************************************/

#include "common/utilities.h"
#include "element_types/elements4D.h"

#include <iostream>  // std::cout etc.
#include <cmath>

#define EPSILON 1.0e-12

using namespace utils;

namespace elements {

  /*
    ==========================
    4D Tesseract element
    ==========================

    The finite element local point numbering for a 16 node Tesseract is
    based on the 3D Hex8 Ensight element
 

    _.15-------------------------------------14
    _.+<    |\                              . >-"/ |
    _ .+>         | \                         .>"" ./    |
    .>""              |  \                     <""    /      |
    12----------------------+------------------13    ./        |
    | )<=               |    \               / | _/""          |
    |     )\+           |     \            / /"|               |
    |         (\=       |   _. 7---------+--6  |               |
    |             \>   .|+<    |       / . "|  |               |
    |               '4--+------+------5'    |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |   _ .3------+-----2_ |               |
    |                |  | "   /       |   /'  '| \= _          |
    |                0--+---+---------1*"      |     ""\       |
    |             ./'   |  "           \       |         "">   |
    |            /      |/              \      |             ".|
    |         /        11----------------+-----+--------------10
    |      ./    .+<""                    )    |           .</
    |    /(   /(                            \  |     _.+</
    | ./  /"                                 \ |  >(
    8------------------------------------------9'

    j
    ^        k
    |      /
    |    / 
    |  /
    |/
    +---------->i

    i = Xi
    j = Eta
    k = Mu
    t = Tau 

  */




  real_t Tess16::ref_vert[ Tess16::num_verts* Tess16::num_dim] = // listed as {Xi, Eta, Mu, Tau}
    {
      // Interior cube bottom
      -1.0, -1.0, -1.0, -1.0,
      +1.0, -1.0, -1.0, -1.0,
      +1.0, -1.0, +1.0, -1.0,
      -1.0, -1.0, +1.0, -1.0,
      // Interior cube top
      -1.0, +1.0, -1.0, -1.0,
      +1.0, +1.0, -1.0, -1.0,
      +1.0, +1.0, +1.0, -1.0,
      -1.0, +1.0, +1.0, -1.0,
      // Exterior cube bottom
      -1.0, -1.0, -1.0, +1.0,
      +1.0, -1.0, -1.0, +1.0,
      +1.0, -1.0, +1.0, +1.0,
      -1.0, -1.0, +1.0, +1.0,
      // Exterior cube top
      -1.0, +1.0, -1.0, +1.0,
      +1.0, +1.0, -1.0, +1.0,
      +1.0, +1.0, +1.0, +1.0,
      -1.0, +1.0, +1.0, +1.0,
    };


  // calculate a physical position in an element for a given xi,eta,mu
  void Tess16::physical_position(
                                 ViewCArray <real_t> &x_point,
                                 const ViewCArray <real_t> &xi_point,
                                 const ViewCArray <real_t> &vertices){

    real_t basis_a[num_verts];
    auto basis = ViewCArray <real_t> (basis_a, num_verts);
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts, num_dim);
   
    // calculate the shape functions from each vertex for (xi,eta,mu, tau)
    for(int this_vert = 0; this_vert < num_verts; this_vert++){
      basis(this_vert) = 1.0/16.0
        * (1.0 + xi_point(0)*ref_verts(this_vert, 0)) 
        * (1.0 + xi_point(1)*ref_verts(this_vert, 1)) 
        * (1.0 + xi_point(2)*ref_verts(this_vert, 2)) 
        * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    } // end for shape functions

    // calculate the position in physical space
    for (int dim = 0; dim < 4; dim++){
      x_point(dim) = 0.0;
    }

    for (int this_vert = 0; this_vert < 16; this_vert++ ){
      for (int dim = 0; dim < 4; dim++){
        x_point(dim) += vertices(this_vert, dim)*basis(this_vert);
      } // end for dim
    } // end for this_vert

  } // End physical position function


  // calculate the value for the basis at each node for a given xi,eta,mu,tau
  void Tess16::basis(
                     ViewCArray <real_t>  &basis,
                     const ViewCArray <real_t>  &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts, num_dim);

    // calculate the basis functions from each vertex for (xi,eta,mu, tau)
    for(int this_vert = 0; this_vert < num_verts; this_vert++){
      basis(this_vert) = 1.0/16.0
        * (1.0 + xi_point(0)*ref_verts(this_vert, 0)) 
        * (1.0 + xi_point(1)*ref_verts(this_vert, 1)) 
        * (1.0 + xi_point(2)*ref_verts(this_vert, 2)) 
        * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    } // end for this_vert

  }


  // Partial derivative of shape functions with respect to Xi at Xi_point
  void Tess16::partial_xi_basis(
                                ViewCArray <real_t> &partial_xi, 
                                const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){
      partial_xi(this_vert) = 1.0/16.0
        * (ref_verts(this_vert, 0))
        * (1.0 + xi_point(1)*ref_verts(this_vert, 1))
        * (1.0 + xi_point(2)*ref_verts(this_vert, 2))
        * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    }   // end for this_vert 

  } // end partial Xi function


  // Partial derivative of shape functions with respect to Eta
  void Tess16::partial_eta_basis(
                                 ViewCArray <real_t> &partial_eta, 
                                 const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){  
      partial_eta(this_vert) = 1.0/16.0
        * (1.0 + xi_point(0)*ref_verts(this_vert, 0))
        * (ref_verts(this_vert, 1))
        * (1.0 + xi_point(2)*ref_verts(this_vert, 2))
        * (1.0 + xi_point(3)*ref_verts(this_vert, 3));               
    }   // end for this_vert 

  }  // End partial eta function


  // Partial derivative of shape functions with respect to Mu
  void Tess16::partial_mu_basis(
                                ViewCArray <real_t> &partial_mu, 
                                const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){  
      partial_mu(this_vert) = 1.0/16.0
        * (1.0 + xi_point(0)*ref_verts(this_vert, 0))
        * (1.0 + xi_point(1)*ref_verts(this_vert, 1))
        * (ref_verts(this_vert, 2))
        * (1.0 + xi_point(3)*ref_verts(this_vert, 3));
    } // end for this_vert 

  } // end partial Mu fuction


  // Partial derivative of shape functions with respect to Tau
  void Tess16::partial_tau_basis(
                                 ViewCArray <real_t> &partial_tau, 
                                 const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts, num_dim);

    for (int this_vert = 0; this_vert < num_verts; this_vert++){  
      partial_tau(this_vert) = 1.0/16.0
        * (1.0 + xi_point(0)*ref_verts(this_vert, 0))
        * (1.0 + xi_point(1)*ref_verts(this_vert, 1))
        * (1.0 + xi_point(2)*ref_verts(this_vert, 2))
        * (ref_verts(this_vert, 3));
    } // end for this_vert   

  } // End partial tau function

} // end namespace elements
