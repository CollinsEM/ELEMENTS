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
#include "element_types/elements3D.h"
#include "element_types/quadrature.h"

#include <iostream>  // std::cout etc.
#include <cmath>

#define EPSILON 1.0e-12

using namespace utils;

namespace elements {



  /* 
     .-------------------------------. 
     | .----------------------------. |
     | |    ______      ________    | |
     | |   / ____ `.   |_   ___ `.  | |
     | |   `'  __) |     | |   `. \ | |
     | |   _  |__ '.     | |    | | | |
     | |  | \____) |    _| |___.' / | |
     | |   \______.'   |________.'  | |
     | |                            | |
     | '----------------------------' |
     '------------------------------' 
  */



  /*
    ==========================
    Hex 8
    ==========================

    The finite element local vertex numbering for a 8 node Hexahedral is
    as follows

    Mu (k)
    |     Eta (j)    
    |    /
    |   /
    6---+----7
    /|   |   /|
    / |   |  / |
    4--------5  |
    |  |    -|--+---> Xi (i)
    |  |     |  |
    |  2-----|--3
    | /      | /       
    |/       |/
    0----*----1
 
  */

  Hex8::Hex8() : Element3D(){
    nsurfaces = 6;
    CArray<size_t> strides(nsurfaces);
    for(int istride=0; istride < nsurfaces; istride++)
      strides(istride) = 4;
    //list of local ids to basis functions needed to interpolate throughout a given element surface
    surface_to_dof_lid = RaggedRightArray<int>(strides);

    //st planes first
    surface_to_dof_lid(0,0) = 0;
    surface_to_dof_lid(0,1) = 1;
    surface_to_dof_lid(0,2) = 2;
    surface_to_dof_lid(0,3) = 3;
  
    surface_to_dof_lid(1,0) = 4;
    surface_to_dof_lid(1,1) = 5;
    surface_to_dof_lid(1,2) = 6;
    surface_to_dof_lid(1,3) = 7;
  
    //sw planes second
    surface_to_dof_lid(2,0) = 0;
    surface_to_dof_lid(2,1) = 1;
    surface_to_dof_lid(2,2) = 4;
    surface_to_dof_lid(2,3) = 5;

    surface_to_dof_lid(3,0) = 2;
    surface_to_dof_lid(3,1) = 3;
    surface_to_dof_lid(3,2) = 6;
    surface_to_dof_lid(3,3) = 7;
  
    //tw planes third
    surface_to_dof_lid(4,0) = 0;
    surface_to_dof_lid(4,1) = 2;
    surface_to_dof_lid(4,2) = 4;
    surface_to_dof_lid(4,3) = 6;

    surface_to_dof_lid(5,0) = 1;
    surface_to_dof_lid(5,1) = 3;
    surface_to_dof_lid(5,2) = 5;
    surface_to_dof_lid(5,3) = 7;

  }

  Hex8::~Hex8(){

  }

  real_t Hex8::ref_vert[Hex8::num_verts_*Hex8::num_dim_] = 
    {// listed as {Xi, Eta, Mu}
      // Bottom Nodes
      -1.0, -1.0, -1.0,// 0
      +1.0, -1.0, -1.0,// 1
      -1.0, +1.0, -1.0,// 2
      +1.0, +1.0, -1.0,// 3
      // Top Nodes
      -1.0, -1.0, +1.0,// 4
      +1.0, -1.0, +1.0,// 5
      -1.0, +1.0, +1.0,// 6
      +1.0, +1.0, +1.0 // 7
    };

  const int Hex8::vert_to_node[Hex8::num_verts_] = 
    {
      0,
      2,
      6,
      8,
      18,
      20,
      24,
      24
    };

  int Hex8::num_verts()
  {
    return Hex8::num_verts_;
  }
  int Hex8::num_nodes()
  {
    return Hex8::num_nodes_;
  }
  int Hex8::num_basis()
  {
    return Hex8::num_basis_;
  }


  // get the physical location for a given xi_point
  void Hex8::physical_position (
                                ViewCArray <real_t>  &x_point, 
                                const ViewCArray <real_t>  &xi_point, 
                                const ViewCArray <real_t>  &vertices){

    real_t basis_a[num_verts_];
    auto basis = ViewCArray <real_t> (basis_a, num_verts_);
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the shape functions from each vertex for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
      basis(vert_lid) = 1.0/8.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim_; dim++){
      x_point(dim) = 0.0;
    }

    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
      for (int dim = 0; dim < num_dim_; dim++){
        x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
      } // end for dim
    } // end for vert_lid

  } // end of function


  // calculate the value for the basis at each node for a given xi,eta, mu
  void Hex8::basis(
                   ViewCArray <real_t>  &basis,
                   const ViewCArray <real_t>  &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the shape functions from each vertex for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
      basis(vert_lid) = 1.0/8.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

  } // end of hex8 basis functions


  // calculate the partials of the shape function 
  // with respect to Xi
  void Hex8::partial_xi_basis(
                              ViewCArray <real_t>  &partial_xi, 
                              const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // std::cout << "Inside partial xi" << std::endl;
    
    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++){

      partial_xi(vert_lid) = (1.0/8.0)
        * (ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));

    } // end for vert_lid
  } // end of partial Xi function


  // with respect to eta
  void Hex8::partial_eta_basis(
                               ViewCArray <real_t> &partial_eta, 
                               const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++){
      partial_eta(vert_lid) = (1.0/8.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

  } //end of partial eta function 


  // with repsect to mu
  void Hex8::partial_mu_basis(
                              ViewCArray <real_t> &partial_mu, 
                              const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++){
      partial_mu(vert_lid) = (1.0/8.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (ref_verts(vert_lid, 2));
    } // end for vert_lid

  } // end of partial mu function

  // Map from vertex to node
  inline int Hex8::vert_node_map( const int vert_lid){
    
    return vert_to_node[vert_lid];

  };


  inline real_t& Hex8::ref_locs(const int vert_lid, const int dim){
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    return ref_verts(vert_lid, dim);
  }




  /*
    ==========================
    Hex 20
    ==========================

    The finite element local point numbering for a 20 node Hexahedral is 
    as follows

    Mu (k)
    |     Eta (j)
    |    /
    |   /

    7----14----6
    /|         /|
    15 |       13 |
    / 19       /  18
    4----12----5   |
    |   |      |   |  --> Xi (i)
    |   |      |   |
    |   3---10-|---2
    16  /      17  /
    | 11       | 9         
    |/         |/
    0-----8----1

  */

  Hex20::Hex20() : Element3D(){
    nsurfaces = 6;
    CArray<size_t> strides(nsurfaces);
    for(int istride=0; istride < nsurfaces; istride++)
      strides(istride) = 8;
    //list of local ids to basis functions needed to interpolate throughout a given element surface
    surface_to_dof_lid = RaggedRightArray<int>(strides);

    //st planes first
    surface_to_dof_lid(0,0) = 0;
    surface_to_dof_lid(0,1) = 8;
    surface_to_dof_lid(0,2) = 1;
    surface_to_dof_lid(0,3) = 11;
    surface_to_dof_lid(0,4) = 9;
    surface_to_dof_lid(0,5) = 3;
    surface_to_dof_lid(0,6) = 10;
    surface_to_dof_lid(0,7) = 2;
  
    surface_to_dof_lid(1,0) = 4;
    surface_to_dof_lid(1,1) = 12;
    surface_to_dof_lid(1,2) = 5;
    surface_to_dof_lid(1,3) = 15;
    surface_to_dof_lid(1,4) = 13;
    surface_to_dof_lid(1,5) = 7;
    surface_to_dof_lid(1,6) = 14;
    surface_to_dof_lid(1,7) = 6;
  
    //sw planes second
    surface_to_dof_lid(3,0) = 0;
    surface_to_dof_lid(3,1) = 8;
    surface_to_dof_lid(3,2) = 1;
    surface_to_dof_lid(3,3) = 16;
    surface_to_dof_lid(3,4) = 17;
    surface_to_dof_lid(3,5) = 4;
    surface_to_dof_lid(3,6) = 12;
    surface_to_dof_lid(3,7) = 5;

    surface_to_dof_lid(3,0) = 3;
    surface_to_dof_lid(3,1) = 10;
    surface_to_dof_lid(3,2) = 2;
    surface_to_dof_lid(3,3) = 19;
    surface_to_dof_lid(3,4) = 18;
    surface_to_dof_lid(3,5) = 7;
    surface_to_dof_lid(3,6) = 14;
    surface_to_dof_lid(3,7) = 6;
  
    //tw planes third
    surface_to_dof_lid(4,0) = 0;
    surface_to_dof_lid(4,1) = 11;
    surface_to_dof_lid(4,2) = 3;
    surface_to_dof_lid(4,3) = 16;
    surface_to_dof_lid(4,4) = 19;
    surface_to_dof_lid(4,5) = 4;
    surface_to_dof_lid(4,6) = 15;
    surface_to_dof_lid(4,7) = 7;

    surface_to_dof_lid(5,0) = 1;
    surface_to_dof_lid(5,1) = 9;
    surface_to_dof_lid(5,2) = 2;
    surface_to_dof_lid(5,3) = 17;
    surface_to_dof_lid(5,4) = 18;
    surface_to_dof_lid(5,5) = 5;
    surface_to_dof_lid(5,6) = 13;
    surface_to_dof_lid(5,7) = 6;

  }

  Hex20::~Hex20(){

  }

  real_t Hex20::ref_vert[Hex20::num_verts_*Hex20::num_dim_] = // listed as {Xi, Eta, Mu}
    // new indices for right hand coordinates
    {
      // bottom corners
      -1.0, -1.0, -1.0,// 0
      +1.0, -1.0, -1.0,// 1
      +1.0, +1.0, -1.0,// 2
      -1.0, +1.0, -1.0,// 3
      // top corners
      -1.0, -1.0, +1.0,// 4
      +1.0, -1.0, +1.0,// 5
      +1.0, +1.0, +1.0,// 6
      -1.0, +1.0, +1.0,// 7
      // bottom edges
      0.0, -1.0, -1.0,// 8
      +1.0,  0.0, -1.0,// 9
      0.0, +1.0, -1.0,// 10
      -1.0,  0.0, -1.0,// 11
      // top edges
      0.0, -1.0, +1.0,// 12
      +1.0,  0.0, +1.0,// 13
      0.0, +1.0, +1.0,// 14
      -1.0,  0.0, +1.0,// 15
      // middle edges
      -1.0, -1.0,  0.0,// 16
      +1.0, -1.0,  0.0,// 17
      +1.0, +1.0,  0.0,// 18
      -1.0, +1.0,  0.0// 19
    };

  const int Hex20::vert_to_node[Hex20::num_verts_] = 
    {
      0,
      4,
      24,
      20,
      100,
      104,
      124,
      120,
      2,
      14,
      22,
      10,
      102,
      114,
      122,
      110,
      50,
      54,
      74,
      70
    };

  int Hex20::num_verts()
  {
    return Hex20::num_verts_;
  }
  int Hex20::num_nodes()
  {
    return Hex20::num_nodes_;
  }
  int Hex20::num_basis()
  {
    return Hex20::num_basis_;
  }

  // get the physical location for a given xi_point
  void Hex20::physical_position (
                                 ViewCArray <real_t>  &x_point, 
                                 const ViewCArray <real_t>  &xi_point, 
                                 const ViewCArray <real_t>  &vertices){

    real_t basis_a[num_verts_];
    auto basis = ViewCArray <real_t> (basis_a, num_verts_);
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid=0; vert_lid<8; vert_lid++){
      basis(vert_lid) = 1.0/8.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
        * (xi_point(0)*ref_verts(vert_lid, 0)
           +  xi_point(1)*ref_verts(vert_lid, 1)
           +  xi_point(2)*ref_verts(vert_lid, 2) - 2.0);
    } // end for vert_lid

    // calculate the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
      basis(vert_lid) = 1.0/4.0
        * (1.0 - xi_point(0)*xi_point(0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the k=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
      basis(vert_lid) = 1.0/4.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 - xi_point(2)*xi_point(2));

    } // end for vert_lid

    // calculate the j=0 edge shape functions pts=[9,11,15,13]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
      basis(vert_lid) = 1.0/4.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim_; dim++) x_point(dim) = 0.0;

    for (int dim = 0; dim < num_dim_; dim++){
      for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
        x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
      }   
    } // end for dim

  } // end of physical position function


  // calculate the value for the basis at each node for a given xi,eta, mu
  void Hex20::basis(
                    ViewCArray <real_t>  &basis,
                    const ViewCArray <real_t>  &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      basis(vert_lid) = 1.0/8.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
        * (xi_point(0)*ref_verts(vert_lid, 0)
           +  xi_point(1)*ref_verts(vert_lid, 1)
           +  xi_point(2)*ref_verts(vert_lid, 2) - 2.0);
    } // end for vert_lid

    // calculate the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
      basis(vert_lid) = 1.0/4.0
        * (1.0 - xi_point(0)*xi_point(0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the k=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
      basis(vert_lid) = 1.0/4.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 - xi_point(2)*xi_point(2)); 
    } // end for vert_lid

    // calculate the j=0 edge shape functions pts=[9,11,15,13]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
      basis(vert_lid) = 1.0/4.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));           
    } // end for vert_lid

  } // end of hex20 basis functions


  // Calculate the partials of the shape functions
  // with respect to Xi
  void  Hex20::partial_xi_basis(
                                ViewCArray <real_t>  &partial_xi, 
                                const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // For 8 Corner shape functions pts=[0,1,2,3,4,5,6,7]
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      partial_xi(vert_lid) = (1.0/8.0) 
        * (ref_verts(vert_lid, 0))
        * (1.0 + (xi_point(1)*ref_verts(vert_lid, 1)))
        * (1.0 + (xi_point(2)*ref_verts(vert_lid, 2)))
        * (2.0 * (xi_point(0)*ref_verts(vert_lid, 0))
           + xi_point(1)*ref_verts(vert_lid, 1)
           + xi_point(2)*ref_verts(vert_lid, 2) - 1.0);  
    } // end for vert_lid

    // for the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
      partial_xi(vert_lid) = (-1.0/2.0)
        * (xi_point(0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the k=0 shape functions pts=[9,11,13,15]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
      partial_xi(vert_lid) = (1.0/4.0)
        * (ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 - xi_point(2)*xi_point(2));
    } // end for vert_lid


    // for the j=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
      partial_xi(vert_lid) = (1.0/4.0)
        * (ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

  } // end of partial Xi function


  // with respect to Eta
  void Hex20::partial_eta_basis(
                                ViewCArray <real_t> &partial_eta, 
                                const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // For 8 Corner shape functions pts=[0,1,2,3,4,5,6,7]
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      partial_eta(vert_lid) = (1.0/8.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
        * (xi_point(0)*ref_verts(vert_lid, 0)
           +  2.0 * xi_point(1)*ref_verts(vert_lid, 1)
           + xi_point(2)*ref_verts(vert_lid, 2) - 1.0);
    } // end for vert_lid

    // for the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
      partial_eta(vert_lid) = (1.0/4.0)
        * (1.0 - (xi_point(0)*xi_point(0)))
        * (ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the j=0 shape functions pts=[9,11,13,15]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
      partial_eta(vert_lid) = (-1.0/2.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (xi_point(1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the k=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
      partial_eta(vert_lid) = (1.0/4.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1))
        * (1.0 - (xi_point(2)*xi_point(2)));
    } // end for vert_lid

  } // end of partial Eta function


  // with repsect to mu
  void Hex20::partial_mu_basis(
                               ViewCArray <real_t> &partial_mu, 
                               const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // For 8 Corner shape functions pts=[0,1,2,3,4,5,6,7]
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      partial_mu(vert_lid) = (1.0/8.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (ref_verts(vert_lid, 2))
        * ((xi_point(0)*ref_verts(vert_lid, 0))
           + (xi_point(1)*ref_verts(vert_lid, 1))
           + (2.0 * xi_point(2)*ref_verts(vert_lid, 2)) - 1.0);
    } // end for vert_lid

    // for the i=0 edge shape functions pts=[8,10,12,14]
    for (int vert_lid = 8; vert_lid <= 14; vert_lid += 2){
      partial_mu(vert_lid) = (1.0/4.0)
        * (1.0 - (xi_point(0)*xi_point(0))) 
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (ref_verts(vert_lid, 2));
    }

    // for the j=0 shape functions pts=[9,11,13,15]
    for (int vert_lid = 9; vert_lid <= 15; vert_lid += 2){
      partial_mu(vert_lid) = (1.0/4.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1))
        * (ref_verts(vert_lid, 2));
    } // end for vert_lid

    // for the j=0 edge shape functions pts=[16,17,18,19]
    for (int vert_lid = 16; vert_lid <= 19; vert_lid++){
      partial_mu(vert_lid) = (-1.0/2.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (xi_point(2));
    } // end for vert_lid

  } // end of partial Mu function



  // Map from vertex to node
  inline int Hex20::vert_node_map( const int vert_lid){

    return vert_to_node[vert_lid];
  };


  inline real_t& Hex20::ref_locs(const int vert_lid, const int dim){
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    return ref_verts(vert_lid, dim);
  }

  /* 
     ==========================
     Hex 32
     ==========================

     The finite element local point numbering for a 32 node Hexahedral is 
     shown below


     Mu (k)
     ^         Eta (j)
     |        /
     |       /
     /
     7----23------22----6
     /|                 /|
     15 |               14 |
     /  |               /  |
     12  31             13   30 
     /    |             /    |
     4-----20-----21----5     |
     |     |            |     |   ----> Xi (i)
     |    27            |     26  
     |     |            |     |
     28     |           29     |
     |     3----19------|18---2
     |    /             |    /
     |  11              |   10
     24  /              25  /
     | 8                | 9         
     |/                 |/
     0----16------17----1
  */


  Hex32::Hex32() : Element3D(){
    nsurfaces = 6;
    CArray<size_t> strides(nsurfaces);
    for(int istride=0; istride < nsurfaces; istride++)
      strides(istride) = 12;
    //list of local ids to basis functions needed to interpolate throughout a given element surface
    surface_to_dof_lid = RaggedRightArray<int>(strides);

    //st planes first
    surface_to_dof_lid(0,0) = 0;
    surface_to_dof_lid(0,1) = 16;
    surface_to_dof_lid(0,2) = 17;
    surface_to_dof_lid(0,3) = 1;
    surface_to_dof_lid(0,4) = 8;
    surface_to_dof_lid(0,5) = 9;
    surface_to_dof_lid(0,6) = 11;
    surface_to_dof_lid(0,7) = 10;
    surface_to_dof_lid(0,8) = 3;
    surface_to_dof_lid(0,9) = 19;
    surface_to_dof_lid(0,10) = 18;
    surface_to_dof_lid(0,11) = 2;
  
    surface_to_dof_lid(1,0) = 4;
    surface_to_dof_lid(1,1) = 20;
    surface_to_dof_lid(1,2) = 21;
    surface_to_dof_lid(1,3) = 5;
    surface_to_dof_lid(1,4) = 12;
    surface_to_dof_lid(1,5) = 13;
    surface_to_dof_lid(1,6) = 15;
    surface_to_dof_lid(1,7) = 14;
    surface_to_dof_lid(1,8) = 7;
    surface_to_dof_lid(1,9) = 23;
    surface_to_dof_lid(1,10) = 22;
    surface_to_dof_lid(1,11) = 6;
  
    //sw planes second
    surface_to_dof_lid(2,0) = 0;
    surface_to_dof_lid(2,1) = 16;
    surface_to_dof_lid(2,2) = 17;
    surface_to_dof_lid(2,3) = 1;
    surface_to_dof_lid(2,4) = 24;
    surface_to_dof_lid(2,5) = 25;
    surface_to_dof_lid(2,6) = 28;
    surface_to_dof_lid(2,7) = 29;
    surface_to_dof_lid(2,8) = 4;
    surface_to_dof_lid(2,9) = 20;
    surface_to_dof_lid(2,10) = 21;
    surface_to_dof_lid(2,11) = 25;

    surface_to_dof_lid(3,0) = 3;
    surface_to_dof_lid(3,1) = 19;
    surface_to_dof_lid(3,2) = 18;
    surface_to_dof_lid(3,3) = 2;
    surface_to_dof_lid(3,4) = 27;
    surface_to_dof_lid(3,5) = 26;
    surface_to_dof_lid(3,6) = 31;
    surface_to_dof_lid(3,7) = 30;
    surface_to_dof_lid(3,8) = 7;
    surface_to_dof_lid(3,9) = 23;
    surface_to_dof_lid(3,10) = 22;
    surface_to_dof_lid(3,11) = 6;
  
    //tw planes third
    surface_to_dof_lid(4,0) = 0;
    surface_to_dof_lid(4,1) = 8;
    surface_to_dof_lid(4,2) = 11;
    surface_to_dof_lid(4,3) = 3;
    surface_to_dof_lid(4,4) = 24;
    surface_to_dof_lid(4,5) = 27;
    surface_to_dof_lid(4,6) = 28;
    surface_to_dof_lid(4,7) = 31;
    surface_to_dof_lid(4,8) = 4;
    surface_to_dof_lid(4,9) = 12;
    surface_to_dof_lid(4,10) = 15;
    surface_to_dof_lid(4,11) = 7;

    surface_to_dof_lid(5,0) = 1;
    surface_to_dof_lid(5,1) = 9;
    surface_to_dof_lid(5,2) = 10;
    surface_to_dof_lid(5,3) = 2;
    surface_to_dof_lid(5,4) = 25;
    surface_to_dof_lid(5,5) = 26;
    surface_to_dof_lid(5,6) = 29;
    surface_to_dof_lid(5,7) = 30;
    surface_to_dof_lid(5,8) = 5;
    surface_to_dof_lid(5,9) = 13;
    surface_to_dof_lid(5,10) = 14;
    surface_to_dof_lid(5,11) = 6;

  }

  Hex32::~Hex32(){

  }

  real_t Hex32::ref_vert[Hex32::num_verts_*Hex32::num_dim_] = // listed as {Xi, Eta, Mu}
    {
      -1.0, -1.0, -1.0,// 0
      +1.0, -1.0, -1.0,// 1
      +1.0, +1.0, -1.0,// 2
      -1.0, +1.0, -1.0,// 3
      -1.0, -1.0, +1.0,// 4
      +1.0, -1.0, +1.0,// 5
      +1.0, +1.0, +1.0,// 6
      -1.0, +1.0, +1.0,// 7
      // Xi/Mu = +- 1/3
      -1.0, -1./3., -1.0,// 8
      1.0, -1./3., -1.0,// 9
      1.0, +1./3., -1.0,// 10
      -1.0, +1./3., -1.0,// 11
      -1.0, -1./3., +1.0,// 12
      1.0, -1./3., +1.0,// 13
      1.0, +1./3., +1.0,// 14
      -1.0, +1./3., +1.0,// 15
      // Eta/Mu = +- 1/3
      -1./3., -1.0, -1.0,// 16
      1./3., -1.0, -1.0,// 17
      1./3., +1.0, -1.0,// 18
      -1./3., +1.0, -1.0,// 19
      -1./3., -1.0,  1.0,// 20
      1./3., -1.0,  1.0,// 21
      1./3., +1.0,  1.0,// 22
      -1./3., +1.0,  1.0,// 23
      // Xi/Eta = +- 1/3
      -1.0, -1.0, -1./3.,// 24
      1.0, -1.0, -1./3.,// 25
      1.0,  1.0, -1./3.,// 26
      -1.0,  1.0, -1./3.,// 27
      -1.0, -1.0,  1./3.,// 28
      1.0, -1.0,  1./3.,// 29
      1.0,  1.0,  1./3.,// 30
      -1.0,  1.0,  1./3.,// 31
    };

  const int Hex32::vert_to_node[Hex32::num_verts_] = 
    {
      0,
      6,
      48,
      42,
      294,
      300,
      342,
      336,
      14,
      20,
      32,
      28,
      308,
      314,
      328,
      322,
      2,
      4,
      46,
      44,
      296,
      298,
      340,
      338,
      98,
      104,
      146,
      140,
      196,
      202,
      244,
      298
    };

  int Hex32::num_verts()
  {
    return Hex32::num_verts_;
  }
  int Hex32::num_nodes()
  {
    return Hex32::num_nodes_;
  }
  int Hex32::num_basis()
  {
    return Hex32::num_basis_;
  }


  // get the physical location for a given xi_point
  void Hex32::physical_position (
                                 ViewCArray <real_t>  &x_point, 
                                 const ViewCArray <real_t>  &xi_point, 
                                 const ViewCArray <real_t>  &vertices){

    real_t basis_a[num_verts_];
    auto basis = ViewCArray <real_t> (basis_a, num_verts_);
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      basis(vert_lid) = (1.0/64.0) 
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
        * (9.0 * xi_point(0)*xi_point(0)
           +  9.0 * xi_point(1)*xi_point(1)
           +  9.0 * xi_point(2)*xi_point(2) - 19.0);  
    } // end for vert_lid

    // calculate the edge shape functions for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
      basis(vert_lid) = (9.0/64.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1))
        * (1.0 + 9.0 * xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid
    
    // calculate the edge shape functions for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
      basis(vert_lid) = (9.0/64.0)
        * (1.0 - xi_point(0)*xi_point(0))
        * (1.0 + 9.0 * xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2)); 
    } // end for vert_lid

    // calculate the edge shape functions for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
      basis(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + 9.0 * xi_point(2)*ref_verts(vert_lid, 2))
        * (1.0 - xi_point(2)*xi_point(2)); 
    } // end for vert_lid


    // calculate the position in physical space
    for (int dim = 0; dim < num_dim_; dim++){
      x_point(dim) = 0.0;
    }

    for (int vert_lid = 0; vert_lid <= num_verts_; vert_lid++ ){
      //std::cout << "Vert :" << vert_lid << std::endl;
      for (int dim = 0; dim < num_dim_; dim++){
        x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
      }
    }

  } // end of physical position function


  void Hex32::basis(
                    ViewCArray <real_t>  &basis,
                    const ViewCArray <real_t>  &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner shape functions for (xi,eta,mu)
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      basis(vert_lid) = (1.0/64.0) 
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2))
        * (9.0 * xi_point(0)*xi_point(0)
           +  9.0 * xi_point(1)*xi_point(1)
           +  9.0 * xi_point(2)*xi_point(2) - 19.0);  
    } // end for vert_lid

    // calculate the edge shape functions for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
      basis(vert_lid) = (9.0/64.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1))
        * (1.0 + 9.0 * xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the edge shape functions for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
      basis(vert_lid) = (9.0/64.0)
        * (1.0 - xi_point(0)*xi_point(0))
        * (1.0 + 9.0 * xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2)*ref_verts(vert_lid, 2)); 
    } // end for vert_lid

    // calculate the edge shape functions for pts=[24-31]
    for (int vert_lid = 24; vert_lid < num_verts_; vert_lid++){
      basis(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (1.0 + 9.0 * xi_point(2)*ref_verts(vert_lid, 2))
        * (1.0 - xi_point(2)*xi_point(2)); 
    } // end for vert_lid

  } // end of hex20 basis functions

  // Calculate the partials of the shape functions
  // with respect to Xi
  void  Hex32::partial_xi_basis(
                                ViewCArray <real_t>  &partial_xi, 
                                const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner partial wrt Xi 
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      partial_xi(vert_lid) = (1.0/64.0) 
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
        *((9.0 * (ref_verts(vert_lid, 0))
           * (xi_point(0)*xi_point(0) +  xi_point(1)*xi_point(1) + xi_point(2)*xi_point(2)))
          + (18.0 * xi_point(0) * (1.0 + xi_point(0)*ref_verts(vert_lid, 0)))
          - (19.0 * ref_verts(vert_lid, 0)));
    } // end for vert_lid

    // calculate the edge partial wrt Xi for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
      partial_xi(vert_lid) = (9.0/64.0) 
        * (ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        * (1.0 + 9.0 * xi_point(2) * ref_verts(vert_lid, 2))
        * (1.0 - xi_point(2) * xi_point(2));
    } // end for vert_lid

    // calculate the edge partial wrt Xi for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
      partial_xi(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
        * (9.0 * ref_verts(vert_lid, 0) * (1.0 - 3.0 * xi_point(0) * xi_point(0))
           - (2 * xi_point(0)));
    } // end for vert_lid

    // calculate the edge partial wrt Xi for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
      partial_xi(vert_lid) = (9.0/64.0) 
        * (ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1) * xi_point(1))
        * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2) * ref_verts(vert_lid, 2));
    } // end for vert_lid

  } // end of partial Xi function


  // with respect to Eta
  // functions for [18-15] and [24-31] were switched 
  void Hex32::partial_eta_basis(
                                ViewCArray <real_t> &partial_eta, 
                                const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner partial wrt Eta 
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      partial_eta(vert_lid) = (1.0/64.0) 
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
        *((9.0 * ref_verts(vert_lid, 1)
           * (xi_point(0)*xi_point(0) 
              +  xi_point(1)*xi_point(1)
              +  xi_point(2)*xi_point(2)))
          + (18.0 * xi_point(1) * (1.0 + xi_point(1)*ref_verts(vert_lid, 1)))
          - (19.0 * ref_verts(vert_lid, 1)));
    } // end for vert_lid

    // calculate the edge partial wrt Eta for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
      partial_eta(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 + xi_point(2) * ref_verts(vert_lid, 2))
        *((9.0 * ref_verts(vert_lid, 1) * (1.0 - 3.0 * xi_point(1) * xi_point(1)))
          - (2.0 * xi_point(1)));
    } // end for vert_lid

    // calculate the edge partial wrt Eta for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
      partial_eta(vert_lid) = (9.0/64.0) 
        * (1.0 - xi_point(0) * xi_point(0))
        * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1))
        * (1.0 + xi_point(2) * ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the edge partial wrt Eta for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
      partial_eta(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1))
        * (1.0 + 9.0 * xi_point(2) * ref_verts(vert_lid, 2))
        * (1.0 - xi_point(2) * xi_point(2));
    } // end for vert_lid

  } // end of partial Eta function


  // with repsect to mu
  // functions for [18-15] and [24-31] were switched 
  void Hex32::partial_mu_basis(
                               ViewCArray <real_t> &partial_mu, 
                               const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the 8 corner partial wrt Mu 
    for (int vert_lid = 0; vert_lid < 8; vert_lid++){
      partial_mu(vert_lid) = (1.0/64.0)
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        *((9.0 * (ref_verts(vert_lid, 2))
           * (xi_point(0)*xi_point(0) 
              +  xi_point(1)*xi_point(1)
              +  xi_point(2)*xi_point(2))) 
          + (18.0 * xi_point(2) * (1 + xi_point(2)*ref_verts(vert_lid, 2)))
          - (19.0 * ref_verts(vert_lid, 2)));
    } // end for vert_lid

    // calculate the edge partial wrt Mu for pts=[8-15]
    for (int vert_lid = 8; vert_lid <= 15; vert_lid++){
      partial_mu(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1) * xi_point(1))
        * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1))
        * (ref_verts(vert_lid, 2));

    } // end for vert_lid

    // calculate the edge partial wrt Mu for pts=[16-23]
    for (int vert_lid = 16; vert_lid <= 23; vert_lid++){
      partial_mu(vert_lid) = (9.0/64.0) 
        * (1.0 - xi_point(0) * xi_point(0))
        * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        * (ref_verts(vert_lid, 2));
    } // end for vert_lid

    // calculate the edge partial wrt Mu for pts=[24-31]
    for (int vert_lid = 24; vert_lid <= 31; vert_lid++){
      partial_mu(vert_lid) = (9.0/64.0) 
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        *((9.0 * ref_verts(vert_lid, 2) 
           * (1.0 -  3.0 * xi_point(2) * xi_point(2)))
          - (2.0 * xi_point(2)));
    } // end for vert_lid

  } // end of partial Mu function


  // Map from vertex to node
  int Hex32::vert_node_map( const int vert_lid){

    return vert_to_node[vert_lid];
  };


  inline real_t& Hex32::ref_locs(const int vert_lid, const int dim){
    
    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    return ref_verts(vert_lid, dim);
  }



  /*
    _   _           _   _ 
    | | | | _____  _| \ | |
    | |_| |/ _ \ \/ /  \| |
    |  _  |  __/>  <| |\  |
    |_| |_|\___/_/\_\_| \_|
                       
    representative linear element for visualization
   
    j
    |     k    
    |    /
    |   /
    6---+----7
    /|   |   /|
    / |   |  / |
    4--------5  |
    |  |    -|--+---> i
    |  |     |  |
    |  2-----|--3
    | /      | /       
    |/       |/
    0--------1
    

  */


  void HexN::setup_HexN(int elem_order){

        

    if(elem_order == 0){

      num_nodes_1d_ = 2;
      num_nodes_ = pow(num_nodes_1d_, 3);

      HexN_Nodes_1d_ = CArray <real_t> (num_nodes_);
      HexN_Nodes_ = CArray <real_t> (num_nodes_, 3);


      // Vertices
      num_verts_1d_ = 2;
      num_verts_ = pow(num_verts_1d_, 3);
      num_basis_ = pow(num_verts_1d_, 3);
            
      HexN_Verts_1d_ = CArray <real_t> (num_verts_);
      HexN_Verts_ = CArray <real_t> (num_verts_, 3);


      Vert_Node_map_ = CArray <size_t> (num_verts_);

      order_ = elem_order+1;


    }


    else{
            
      // Nodes
      num_nodes_1d_ = 2 * elem_order + 1;
      num_nodes_ = pow(num_nodes_1d_, 3);

      HexN_Nodes_1d_ = CArray <real_t> (num_nodes_);
      HexN_Nodes_ = CArray <real_t> (num_nodes_, 3);


      // Vertices
      num_verts_1d_ = elem_order + 1;
      num_verts_ = pow(num_verts_1d_, 3);
      num_basis_ = pow(num_verts_1d_, 3);
            
      HexN_Verts_1d_ = CArray <real_t> (num_verts_);
      HexN_Verts_ = CArray <real_t> (num_verts_, 3);


      Vert_Node_map_ = CArray <size_t> (num_verts_);

      order_ = elem_order;

    }
        
    create_lobatto_nodes(elem_order);


    // Set the vertex to node map (every other node)
    if(elem_order == 0){

      int vert_rid = 0;
      for(int k = 0; k < num_nodes_1d_; k++){
        for(int j = 0; j < num_nodes_1d_; j++){
          for(int i = 0; i < num_nodes_1d_; i++){

            int node_id = node_rid(i, j, k);
                        
            Vert_Node_map_(vert_rid) = node_id;

            vert_rid++;                        
          }   
        }
      }


    }


    if (elem_order >= 1){

      int vert_rid = 0;
      for(int k = 0; k < num_nodes_1d_; k=k+2){
        for(int j = 0; j < num_nodes_1d_; j=j+2){
          for(int i = 0; i < num_nodes_1d_; i=i+2){

            int node_id = node_rid(i, j, k);
            Vert_Node_map_(vert_rid) = node_id;
            vert_rid++;
          }   
        }
      }
    }
    
  }

  int HexN::num_verts()
  {
    return HexN::num_verts_;
  };
  int HexN::num_nodes()
  {
    return HexN::num_nodes_;
  };
  int HexN::num_basis()
  {
    return HexN::num_basis_;
  };

  real_t& HexN::node_coords(int node_rlid, int this_dim)
  {
    return HexN_Nodes_(node_rlid, this_dim);
  };


  int HexN::vert_node_map(int vert_rid) const
  {
    return Vert_Node_map_(vert_rid);
  }


  int HexN::node_rid(int i, int j, int k) const 
  {
    return i + j*num_nodes_1d_ + k*num_nodes_1d_*num_nodes_1d_;
  };

  int HexN::vert_rid(int i, int j, int k) const 
  {
    return i + j*num_verts_1d_ + k*num_verts_1d_*num_verts_1d_;
  };

  void HexN::basis(CArray <real_t> &basis, CArray <real_t> &point)
  {

    auto val_1d = CArray <real_t> (num_verts_1d_);
    auto val_3d = CArray <real_t> (num_verts_1d_, 3);

    // Calculate 1D basis for the X coordinate of the point
    lagrange_basis_1D(val_1d, point(0));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
      val_3d(i, 0) = val_1d(i);
      val_1d(i) = 0.0;
    }

    // Calculate 1D basis for the Y coordinate of the point
    lagrange_basis_1D(val_1d, point(1));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
      val_3d(i, 1) = val_1d(i);
      val_1d(i) = 0.0;
    }

    // Calculate 1D basis for the Z coordinate of the point
    lagrange_basis_1D(val_1d, point(2));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
      val_3d(i, 2) = val_1d(i);
      val_1d(i) = 0.0;
    }
        
    // Multiply the i, j, k components of the basis from each node
    // to get the tensor product basis for the node
    for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
        for(int i = 0; i < num_verts_1d_; i++){

          int vert_rlid = vert_rid(i,j,k);
          basis(vert_rlid) = val_3d(i, 0)*val_3d(j, 1)*val_3d(k, 2);
        }
      }
    }
  };


  void HexN::partial_xi_basis(CArray <real_t> &partial_xi, CArray <real_t> &point)
  {


    auto val_1d = CArray <real_t> (num_verts_1d_);
    auto val_3d = CArray <real_t> (num_verts_1d_, 3);

    auto Dval_1d = CArray <real_t> (num_verts_1d_);
    auto Dval_3d = CArray <real_t> (num_verts_1d_, 3);

    // Calculate 1D partial w.r.t. xi for the X coordinate of the point
    lagrange_derivative_1D(Dval_1d, point(0));


    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      Dval_3d(i,0) = Dval_1d(i);
      Dval_1d(i) = 0.0;
    }


    // Calculate 1D basis for the Y coordinate of the point
    lagrange_basis_1D(val_1d, point(1));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      val_3d(i,1) = val_1d(i);
      val_1d(i) = 0.0;
    }


    // Calculate 1D basis for the Z coordinate of the point
    lagrange_basis_1D(val_1d, point(2));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      val_3d(i,2) = val_1d(i);
      val_1d(i) = 0.0;
    }

    // Multiply the i, j, k components of the basis and partial_xi from each node
    // to get the tensor product partial derivatives of the basis at each node
    for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
        for(int i = 0; i < num_verts_1d_; i++){
                    
          int vert_rlid = vert_rid(i,j,k);

          // Partial w.r.t xi
          partial_xi(vert_rlid) = Dval_3d(i, 0)*val_3d(j, 1)*val_3d(k, 2);

        }
      }
    }
  };


  void HexN::partial_eta_basis(CArray <real_t> &partial_eta, CArray <real_t> &point)
  {   

    auto val_1d = CArray <real_t> (num_verts_1d_);
    auto val_3d = CArray <real_t> (num_verts_1d_, 3);

    auto Dval_1d = CArray <real_t> (num_verts_1d_);
    auto Dval_3d = CArray <real_t> (num_verts_1d_, 3);

    // Calculate 1D basis for the Y coordinate of the point
    lagrange_basis_1D(val_1d, point(0));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      val_3d(i,0) = val_1d(i);
      val_1d(i) = 0.0;
    }


    // Calculate 1D partial w.r.t. eta for the Y coordinate of the point
    lagrange_derivative_1D(Dval_1d, point(1));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      Dval_3d(i,1) = Dval_1d(i);
      Dval_1d(i) = 0.0;
    }


    // Calculate 1D basis for the Z coordinate of the point
    lagrange_basis_1D(val_1d, point(2));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      val_3d(i,2) = val_1d(i);
      val_1d(i) = 0.0;
    }

    // Multiply the i, j, k components of the basis and partial_eta from each node
    // to get the tensor product partial derivatives of the basis at each node
    for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
        for(int i = 0; i < num_verts_1d_; i++){
                    
          int vert_rlid = vert_rid(i,j,k);

          // Partial w.r.t xi
          partial_eta(vert_rlid) = val_3d(i, 0)*Dval_3d(j, 1)*val_3d(k, 2);

        }
      }
    }
  };


  void HexN::partial_mu_basis(CArray <real_t> &partial_mu, CArray <real_t> &point){


    auto val_1d = CArray <real_t> (num_verts_1d_);
    auto val_3d = CArray <real_t> (num_verts_1d_, 3);

    auto Dval_1d = CArray <real_t> (num_verts_1d_);
    auto Dval_3d = CArray <real_t> (num_verts_1d_, 3);

    // Calculate 1D basis for the X coordinate of the point
    lagrange_basis_1D(val_1d, point(0));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      val_3d(i,0) = val_1d(i);
      val_1d(i) = 0.0;
    }


    // Calculate 1D basis for the Y coordinate of the point
    lagrange_basis_1D(val_1d, point(1));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      val_3d(i,1) = val_1d(i);
      val_1d(i) = 0.0;
    }


    // Calculate 1D partial w.r.t. mu for the Z coordinate of the point
    lagrange_derivative_1D(Dval_1d, point(2));
        
    // Save the basis value at the point to a temp array and zero out the temp array
    for(int i = 0; i < num_verts_1d_; i++){
            
      Dval_3d(i,2) = Dval_1d(i);
      val_1d(i) = 0.0;
    }

    // Multiply the i, j, k components of the basis and partial_xi from each node
    // to get the tensor product partial derivatives of the basis at each node
    for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
        for(int i = 0; i < num_verts_1d_; i++){
                    
          int vert_rlid = vert_rid(i,j,k);

          // Partial w.r.t mu
          partial_mu(vert_rlid) = val_3d(i, 0)*val_3d(j, 1)*Dval_3d(k, 2);

        }
      }
    }
  };



  void HexN::lagrange_basis_1D(
                               CArray <real_t> &interp,    // interpolant from each basis
                               const real_t &x_point){     // point of interest in element
                         
        
    // calculate the basis value associated with each node_i
    for(int vert_i = 0; vert_i < num_verts_1d_; vert_i++){ 
            
      real_t numerator = 1.0;         // placeholder numerator
      real_t denominator = 1.0;       // placeholder denominator
      real_t interpolant = 1.0;       // placeholder value of numerator/denominator
            

      for(int vert_j = 0; vert_j < num_verts_1d_; vert_j++){  // looping over the verts !=vert_i
        if (vert_j != vert_i ){
                    
          // Calculate the numerator
          numerator = numerator*(x_point - HexN_Verts_1d_(vert_j));
                    
          // Calculate the denominator 
          denominator = denominator*(HexN_Verts_1d_(vert_i) - HexN_Verts_1d_(vert_j));
                
        }//end if
                
        interpolant = numerator/denominator; // storing a single value for interpolation for node vert_i
                
      } // end looping over nodes != vert_i

      // writing value to vectors for later use
      interp(vert_i)   = interpolant;           // Interpolant value at given point

    } // end loop over all nodes
  } // end of Legrange_1D function

  void HexN::lagrange_derivative_1D(
                                    CArray <real_t> &derivative,    // derivative
                                    const real_t &x_point){         // point of interest in element

    for(int vert_i = 0; vert_i < num_verts_1d_; vert_i++){ // looping over the nodes
        

      real_t denominator = 1.0;       // placeholder denominator
      real_t num_gradient = 0.0;      // placeholder for numerator of the gradient
      real_t gradient = 0.0;

      for(int vert_j = 0; vert_j < num_verts_1d_; vert_j++){

        // std::cout<<"HexN 1D Vert "<<vert_j<< " = "<<HexN_Verts_1d_(vert_j)<<std::endl;
      }

      for(int vert_j = 0; vert_j < num_verts_1d_; vert_j++){  // looping over the nodes !=vert_i
        if (vert_j != vert_i ){

          // Calculate the denominator that is the same for 
          // both the basis and the gradient of the basis
          denominator = denominator*(HexN_Verts_1d_(vert_i) - HexN_Verts_1d_(vert_j));
                    
          real_t product_gradient = 1.0;
                    
          // Calculate the numerator of the gradient
          for(int N = 0; N < num_verts_1d_; N++){  // looping over the nodes !=vert_i
                        
            if (N != vert_j && N != vert_i ){
              product_gradient = product_gradient * (x_point - HexN_Verts_1d_(N));
            }// end if
          }//end for
                    
          // Sum over the product of the numerator 
          // contributions from each node
          num_gradient += product_gradient; 
                
        }//end if
                
        gradient = (num_gradient/denominator); // storing the derivative of the interpolating function
            
      } // end looping over nodes != vert_i

      // writing value to vectors for later use
      derivative(vert_i)  = gradient;    // derivative of each function

    } // end loop over all nodes
  } // end of Legrange_1D function

  void HexN::create_lobatto_nodes(int element_order){

    int num_nodes_1d = 0;
    int num_nodes_3d = 0;

    if( element_order == 0){

      num_nodes_1d = 2;
      num_nodes_3d = pow(num_nodes_1d, 3);

    }

    else{
      num_nodes_1d = 2.0 * element_order + 1;
      num_nodes_3d = pow(num_nodes_1d, 3);
    }

    // --- build gauss nodal positions and weights ---
    // std::cout<< "Num Nodes passed to create lobatto nodes "<<std::endl;
    //elements::lobatto_nodes_1D(HexN_Nodes_1d_, num_nodes_1d);
    lobatto_nodes_1D(HexN_Nodes_1d_, num_nodes_1d);
        
    for(int num_k = 0; num_k < num_nodes_1d; num_k++){
      for(int num_j = 0; num_j < num_nodes_1d; num_j++){
        for(int num_i = 0; num_i < num_nodes_1d; num_i++){
        
          int node_rlid = node_rid(num_i, num_j, num_k);

          HexN_Nodes_(node_rlid, 0) = HexN_Nodes_1d_(num_i);
          HexN_Nodes_(node_rlid, 1) = HexN_Nodes_1d_(num_j);
          HexN_Nodes_(node_rlid, 2) = HexN_Nodes_1d_(num_k);
        }
      }
    }

    // Saving vertex positions in 1D
    if( element_order == 0){

      int vert_id = 0;
      for(int i = 0; i < num_nodes_1d; i++){

        HexN_Verts_1d_(vert_id) = HexN_Nodes_1d_(i);

        vert_id++;
      }  

    }

    else{
            
      int vert_id = 0;
      for(int i = 0; i < num_nodes_1d; i=i+2){

        HexN_Verts_1d_(vert_id) = HexN_Nodes_1d_(i);

        vert_id++;
      }  
    }
  }

} // end namespace elements
