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
#include "element_types/elements2D.h"

#include <iostream>  // std::cout etc.
#include <cmath>

#define EPSILON 1.0e-12

using namespace utils;

namespace elements {

  /*
    .-------------------------------. 
    | .----------------------------. |
    | |    _____       ________    | |
    | |   / ___ `.    |_   ___ `.  | |
    | |  |_/___) |      | |   `. \ | |
    | |   .'____.'      | |    | | | |
    | |  / /____       _| |___.' / | |
    | |  |_______|    |________.'  | |
    | |                            | |
    | '----------------------------' |
    '-------------------------------' 
  */

  /*
    ===========================
    2D Quad 4 Elements
    ===========================


    The finite element local point numbering for a 4 node Quadrilateral is
    as follows

    Eta
    ^
    |
    3------+-----2
    |      |     |
    |      |     |
    |      |     |
    |      ------+------> Xi   
    |            |
    |            |
    0------------1

  */

  Quad4::Quad4() : Element2D(){
    nsurfaces = 4;
    CArray<size_t> strides(nsurfaces);
    for(int istride=0; istride < nsurfaces; istride++)
      strides(istride) = 2;
    //list of local ids to basis functions needed to interpolate throughout a given element surface
    surface_to_dof_lid = RaggedRightArray<int>(strides);
    surface_to_dof_lid(0,0) = 0;
    surface_to_dof_lid(0,1) = 1;

    surface_to_dof_lid(1,0) = 3;
    surface_to_dof_lid(1,1) = 2;

    surface_to_dof_lid(2,0) = 0;
    surface_to_dof_lid(2,1) = 3;

    surface_to_dof_lid(3,0) = 1;
    surface_to_dof_lid(3,1) = 2;

  }

  Quad4::~Quad4(){

  }

  real_t Quad4::ref_vert[Quad4::num_verts_*Quad4::num_dim_] = 
    
    {// listed as {Xi, Eta}
      -1.0, -1.0,// 0
      1.0, -1.0,// 1
      1.0,  1.0,// 2
      -1.0,  1.0,// 3
    };

  const int Quad4::vert_to_node[Quad4::num_verts_] = 
    {
      0,
      2,
      6,
      8
    };


  int Quad4::num_verts()
  {
    return Quad4::num_verts_;
  }
  int Quad4::num_nodes()
  {
    return Quad4::num_nodes_;
  }
  int Quad4::num_basis()
  {
    return Quad4::num_basis_;
  }


  // calculate a physical position in an element for a given xi,eta
  void Quad4::physical_position(
                                ViewCArray <real_t> &x_point, 
                                const ViewCArray <real_t> &xi_point, 
                                const ViewCArray <real_t> &vertices){

    real_t basis_a[num_verts_];
    auto basis = ViewCArray <real_t> (basis_a, num_verts_);

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the shape functions from each vertex for 0 through num_verts_(xi,eta)
    for( int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
      basis(vert_lid) = 1.0/4.0
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    }// end for

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim_; dim++) x_point(dim) = 0.0;
    
    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
      for (int dim = 0; dim < num_dim_; dim++){
        x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
      }// end for dim
    } // end for vert_lid

  } // end of physical position functionfunction


  // calculate the value for the basis at each node for a given xi,eta
  void Quad4::basis(
                    ViewCArray <real_t> &basis,
                    const ViewCArray <real_t> &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);
    
    // calculate the shape functions from each vertex for 0 through num_verts_(xi,eta)
    for( int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){

      basis(vert_lid) = (1.0/4.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));

    }// end for

  }// end of quad4 basis functions 


  // Partial derivative of shape functions with respect to Xi
  void  Quad4::partial_xi_basis(
                                ViewCArray <real_t>  &partial_xi, 
                                const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);
    
    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++){
      partial_xi(vert_lid) = (1.0/4.0)
        * (ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    }

  }// end of partial xi funciton


  // Partial derivative of shape functions with respect to Eta
  void  Quad4::partial_eta_basis(
                                 ViewCArray <real_t> &partial_eta, 
                                 const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++){
      partial_eta(vert_lid) = (1.0/4.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1));
    }

  }// end of partial eta function

  int Quad4::vert_node_map(const int vert_lid){

    return vert_to_node[vert_lid];

  }



  /*
    ===========================
    2D Quad 8 Elements
    ===========================


    The finite element local point numbering for a 8 node Quadrilateral is
    as follows

    Eta
    ^
    |
    3-------6------2
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    7       +------5-----> Xi   
    |              |
    |              |
    |              |
    0------4-------1

  */

  Quad8::Quad8() : Element2D(){
    nsurfaces = 4;
    CArray<size_t> strides(nsurfaces);
    for(int istride=0; istride < nsurfaces; istride++)
      strides(istride) = 3;
    //list of local ids to basis functions needed to interpolate throughout a given element surface
    surface_to_dof_lid = RaggedRightArray<int>(strides);

    surface_to_dof_lid(0,0) = 0;
    surface_to_dof_lid(0,1) = 4;
    surface_to_dof_lid(0,2) = 1;

    surface_to_dof_lid(1,0) = 3;
    surface_to_dof_lid(1,1) = 6;
    surface_to_dof_lid(1,2) = 2;

    surface_to_dof_lid(2,0) = 0;
    surface_to_dof_lid(2,1) = 7;
    surface_to_dof_lid(2,2) = 3;

    surface_to_dof_lid(3,0) = 1;
    surface_to_dof_lid(3,1) = 5;
    surface_to_dof_lid(3,2) = 2;

  }

  Quad8::~Quad8(){

  }

  real_t Quad8::ref_vert[Quad8::num_verts_*Quad8::num_dim_] = // listed as {Xi, Eta}
    {// listed as {Xi, Eta}
      -1.0, -1.0, // 0  
      1.0, -1.0, // 1
      1.0,  1.0, // 2
      -1.0,  1.0, // 3
      // midline nodes
      0.0, -1.0, // 4
      1.0,  0.0, // 5
      0.0,  1.0, // 6
      -1.0,  0.0, // 7
    };

  const int Quad8::vert_to_node[Quad8::num_verts_] = 
    {
      0,
      4,
      24,
      20,
      2,
      14,
      23,
      10
    };

  int Quad8::num_verts()
  {
    return Quad8::num_verts_;
  }
  int Quad8::num_nodes()
  {
    return Quad8::num_nodes_;
  }
  int Quad8::num_basis()
  {
    return Quad8::num_basis_;
  }

  // calculate a physical position in an element for a given xi,eta,
  void Quad8::physical_position(
                                ViewCArray <real_t> &x_point, 
                                const ViewCArray <real_t> &xi_point, 
                                const ViewCArray <real_t> &vertices){

    real_t basis_a[num_verts_];
    auto basis = ViewCArray <real_t> (basis_a, num_verts_);

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      basis(vert_lid) = (1.0/4.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (xi_point(0)*ref_verts(vert_lid, 0) 
           +  xi_point(1)*ref_verts(vert_lid, 1) - 1.0);
    } // end for vert_lid

    // calculate the shape functions for node 4,6(xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
      basis(vert_lid) = (1.0/2.0)
        * (1.0 - xi_point(0)*xi_point(0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the shape functions for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
      basis(vert_lid) = (1.0/2.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1));
    } // end for vert_lid

    // calculate the position in physical space
    for (int dim = 0; dim < num_dim_; dim++) x_point(dim) = 0.0;

    for (int vert_lid = 0; vert_lid < num_verts_; vert_lid++ ){
      for (int dim = 0; dim < num_dim_; dim++){
        x_point(dim) += vertices(vert_lid, dim)*basis(vert_lid);
      } // end for dim
    } // end for vert_lid

  } // end of function


  // calculate the value for the basis at each node for a given xi,eta
  void Quad8::basis(
                    ViewCArray <real_t> &basis,
                    const ViewCArray <real_t> &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);
    
    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      basis(vert_lid) = (1.0/4.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1))
        * (xi_point(0)*ref_verts(vert_lid, 0) 
           +  xi_point(1)*ref_verts(vert_lid, 1) - 1.0);
    } // end for vert_lid


    // calculate the shape functions for node 4,6(xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
      basis(vert_lid) = (1.0/2.0)
        * (1.0 - xi_point(0)*xi_point(0))
        * (1.0 + xi_point(1)*ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the shape functions for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
      basis(vert_lid) = (1.0/2.0)
        * (1.0 + xi_point(0)*ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1));

    } // end for vert_lid

  }// end of quad8 basis functions


  // Partial derivative of shape functions with respect to Xi
  void Quad8::partial_xi_basis(
                               ViewCArray <real_t>  &partial_xi, 
                               const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the Xi partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      partial_xi(vert_lid) = 1.0/4.0
        * (ref_verts(vert_lid, 0))
        * (1.0 + ref_verts(vert_lid, 1)*xi_point(1))
        *((2.0 * ref_verts(vert_lid, 0)*xi_point(0)) 
          + (ref_verts(vert_lid, 1)*xi_point(1)));
    } // end for vert_lid


    // calculate the Xi partials for node 4,6 (xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
      partial_xi(vert_lid) = -1.0
        * (xi_point(0))
        * (1 + ref_verts(vert_lid, 1)*xi_point(1));
    } // end for vert_lid

    // calculate the Xi partials for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
      partial_xi(vert_lid) = 1.0/2.0
        * (ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1)*xi_point(1));

    } // end for vert_lid

  } // end partial Xi function


  // Partial derivative of shape functions with respect to Eta
  void Quad8::partial_eta_basis(
                                ViewCArray <real_t>  &partial_eta, 
                                const ViewCArray <real_t> &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the Eta partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      partial_eta(vert_lid) = (1.0/4.0)
        * (1.0 + ref_verts(vert_lid, 0)*xi_point(0))
        * (ref_verts(vert_lid, 1))
        *((ref_verts(vert_lid, 0)*xi_point(0))
          + (2.0 * ref_verts(vert_lid, 1)*xi_point(1))); 
    } // end for vert_lid

    // calculate the Eta partials for node 4,6 (xi,eta)
    for( int vert_lid = 4; vert_lid <= 6; vert_lid += 2 ){
      partial_eta(vert_lid) = (1.0/2.0)
        * (1.0 - xi_point(0)*xi_point(0))
        * (ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the Eta partials for node 5,7 (xi,eta)
    for( int vert_lid = 5; vert_lid <= 7; vert_lid += 2 ){
      partial_eta(vert_lid) = (-1.0)
        * (1.0 + ref_verts(vert_lid, 0)*xi_point(0))
        * (xi_point(1));
    } // end for vert_lid

  } // end partial Eta function

  int Quad8::vert_node_map(const int vert_lid){

    return vert_to_node[vert_lid];

  }


  /*
    ===========================
    2D Quad 12 Elements
    ===========================


    The finite element local point numbering for a 8 node Quadrilateral is
    as follows

    Eta
    ^
    |
    3---7------6---2
    |       |      |
    |       |      |
    11       |      10
    |       |      |
    |       +------|-----> Xi   
    |              |
    8              9
    |              |
    0----4-----5---1

  */

  Quad12::Quad12() : Element2D(){
    nsurfaces = 4;
    CArray<size_t> strides(nsurfaces);
    for(int istride=0; istride < nsurfaces; istride++)
      strides(istride) = 4;
    //list of local ids to basis functions needed to interpolate throughout a given element surface
    surface_to_dof_lid = RaggedRightArray<int>(strides);
    surface_to_dof_lid(0,0) = 0;
    surface_to_dof_lid(0,1) = 4;
    surface_to_dof_lid(0,2) = 5;
    surface_to_dof_lid(0,3) = 1;

    surface_to_dof_lid(1,0) = 3;
    surface_to_dof_lid(1,1) = 7;
    surface_to_dof_lid(1,2) = 6;
    surface_to_dof_lid(1,3) = 2;

    surface_to_dof_lid(2,0) = 0;
    surface_to_dof_lid(2,1) = 8;
    surface_to_dof_lid(2,2) = 11;
    surface_to_dof_lid(2,3) = 3;

    surface_to_dof_lid(3,0) = 1;
    surface_to_dof_lid(3,1) = 9;
    surface_to_dof_lid(3,2) = 10;
    surface_to_dof_lid(3,3) = 2;

  }

  Quad12::~Quad12(){

  }
      
  real_t Quad12::ref_vert[Quad12::num_verts_*Quad12::num_dim_] = 
    {// listed as {Xi, Eta}
      //corner nodes
      -1.0, -1.0 ,// 0
      1.0, -1.0 ,// 1
      1.0,  1.0 ,// 2
      -1.0,  1.0 ,// 3
      // Eta +- 1./3.
      -1./3., -1.0 ,// 4
      1./3., -1.0 ,// 5
      1./3.,  1.0 ,// 6
      -1./3.,  1.0 ,// 7
      // Xi +- 1./3.
      -1.0, -1./3. ,// 8
      1.0, -1./3. ,// 9
      1.0,  1./3. ,// 10
      -1.0,  1./3. ,// 11
    };

  const int Quad12::vert_to_node[Quad12::num_verts_] = 
    {
      0,
      6,
      48,
      42,
      2,
      4,
      46,
      44,
      14,
      20,
      34,
      28
    };


  int Quad12::num_verts()
  {
    return Quad12::num_verts_;
  }
  int Quad12::num_nodes()
  {
    return Quad12::num_nodes_;
  }
  int Quad12::num_basis()
  {
    return Quad12::num_basis_;
  }



  // calculate a physical position in an element for a given xi,eta,
  void Quad12::physical_position(
                                 ViewCArray <real_t>  &x_point, 
                                 const ViewCArray <real_t>  &xi_point, 
                                 const ViewCArray <real_t>  &vertices){

    real_t basis_a[num_verts_];
    auto basis = ViewCArray <real_t> (basis_a, num_verts_);


    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      basis(vert_lid) = 1.0/32.0
        * (1.0 + ref_verts(vert_lid, 0) * xi_point(0))
        * (1.0 + ref_verts(vert_lid, 1) * xi_point(1))
        * (9.0 * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1))
           - (10.0));

    } // end for vert_lid

    // calculate the shape functions for node 4-7(xi,eta)
    for( int vert_lid = 4; vert_lid <= 7; vert_lid++ ){
      basis(vert_lid) = 9.0/32.0
        * (1.0 - xi_point(0) * xi_point(0))
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0));
    } // end for vert_lid

    // calculate the shape functions for node 8-11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++ ){
      basis(vert_lid) = 9.0/32.0
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1) * xi_point(1))
        * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1));
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


  // calculate the value for the basis at each node for a given xi,eta
  void Quad12::basis(
                     ViewCArray <real_t>  &basis,
                     const ViewCArray <real_t>  &xi_point){

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the shape functions for node 0,1,2,3(xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      basis(vert_lid) = 1.0/32.0
        * (1.0 + ref_verts(vert_lid, 0) * xi_point(0))
        * (1.0 + ref_verts(vert_lid, 1) * xi_point(1))
        * (9.0 * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1))
           - (10.0));

    } // end for vert_lid

    // calculate the shape functions for node 4-7(xi,eta)
    for( int vert_lid = 4; vert_lid <= 7; vert_lid++ ){
      basis(vert_lid) = 9.0/32.0
        * (1.0 - xi_point(0) * xi_point(0))
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0));
    } // end for vert_lid

    // calculate the shape functions for node 8-11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++ ){
      basis(vert_lid) = 9.0/32.0
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1) * xi_point(1))
        * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1));
    } // end for vert_lid

  }// end of quad12 basis functions


  // Partial derivative of shape functions with respect to Xi
  void Quad12::partial_xi_basis(
                                ViewCArray <real_t>  &partial_xi, 
                                const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);

    // calculate the Xi partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      partial_xi(vert_lid) = 1.0/32.0
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        *((9.0 * ref_verts(vert_lid, 0) 
           * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1)))
          + (18.0 * xi_point(0) * (1.0 + xi_point(0) * ref_verts(vert_lid, 0)))
          - (10.0 * ref_verts(vert_lid, 0)));
    } // end for vert_lid

    // calculate the Xi partials for node 4,5,6,7 (xi,eta)
    for( int vert_lid = 4; vert_lid < 8; vert_lid++ ){
      partial_xi(vert_lid) = (9.0/32.0) 
        * (1.0 + xi_point(1) * ref_verts(vert_lid, 1))
        *((9.0 * ref_verts(vert_lid, 0) 
           * (1.0 - 3.0 * xi_point(0)*xi_point(0)))
          - (2.0 * xi_point(0)));
    } // end for vert_lid

    // calculate the Xi partials for node 8,9,10,11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++){
      partial_xi(vert_lid) = 9.0/32.0
        * (ref_verts(vert_lid, 0))
        * (1.0 - xi_point(1) * xi_point(1))
        * (1.0 + 9.0 * xi_point(1) * ref_verts(vert_lid, 1));
    } // end for vert_lid

  } // end partial Xi function


  // Partial derivative of shape functions with respect to Eta
  void Quad12::partial_eta_basis(
                                 ViewCArray <real_t> &partial_eta, 
                                 const ViewCArray <real_t>  &xi_point) {

    auto ref_verts = ViewCArray<real_t> (ref_vert, num_verts_, num_dim_);
    // calculate the Eta partials for node 0,1,2,3 (xi,eta)
    for( int vert_lid = 0; vert_lid < 4; vert_lid++ ){
      partial_eta(vert_lid) = 1.0/32.0
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        *((9.0 * ref_verts(vert_lid, 1) 
           * (xi_point(0) * xi_point(0) + xi_point(1) * xi_point(1)))
          + (18.0 * xi_point(1) * (1.0 + xi_point(1) * ref_verts(vert_lid, 1)))
          - (10.0 * ref_verts(vert_lid, 1)));
    } // end for vert_lid

    // calculate the Eta partials for node 4,5,6,7 (xi,eta)
    for( int vert_lid = 4; vert_lid <= 7; vert_lid++ ){
      partial_eta(vert_lid) = 9.0/32.0
        * (1.0 - xi_point(0) * xi_point(0))
        * (1.0 + 9.0 * xi_point(0) * ref_verts(vert_lid, 0))
        * (ref_verts(vert_lid, 1));
    } // end for vert_lid

    // calculate the Eta partials for node 8,9,10,11 (xi,eta)
    for( int vert_lid = 8; vert_lid <= 11; vert_lid++){
      partial_eta(vert_lid) = 9.0/32.0
        * (1.0 + xi_point(0) * ref_verts(vert_lid, 0))
        *((9.0 * ref_verts(vert_lid, 1) * (1.0 - 3.0 * xi_point(1)*xi_point(1)))
          - (2.0 * xi_point(1)));

    } // end for vert_lid

  } // end partial Eta function

  int Quad12::vert_node_map(const int vert_lid){

    return vert_to_node[vert_lid];
  }



  /*
    ==========================
    Arbitrary Order Elements
    ==========================


    / __ \                | | \ | |
    | |  | |_   _  __ _  __| |  \| |
    | |  | | | | |/ _` |/ _` | . ` |
    | |__| | |_| | (_| | (_| | |\  |
    \___\_\\__,_|\__,_|\__,_|_| \_| 

    Representative linear element for visualization
   
    Eta (j)
    ^
    |
    3--------------2
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       +------|-----> Xi (i) 
    |              |
    |              |
    |              |
    0--------------1
  */


  // Lagrange Interp in 1D, returns interpolants and derivative
  // works with any nodal spacing
  void QuadN::lagrange_1D(
                          ViewCArray <real_t> &interp,          // interpolant
                          ViewCArray <real_t> &Dinterp,         // derivative of function
                          const real_t &x_point,                  // point of interest in element
                          const ViewCArray <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
                          const int &orderN){                     // order of element

    real_t num_a[orderN+1];
    auto num = ViewCArray <real_t> (num_a, orderN+1); // numerator of interpolant

    real_t denom_a[orderN+1];
    auto denom = ViewCArray <real_t> (denom_a, orderN+1); // denomenator of interpolant
  
    real_t q = 0.0;
   
    for(int i = 0; i < orderN + 1; i++){ // looping over the nodes
      real_t n = 1.0;         // placeholder numerator
      real_t d = 1.0;         // placeholder denominator
      real_t c = 1.0;         // placeholder value of n/d
      real_t p = 0.0;         // placeholder for derivative values
      real_t s = 1.0;
        
      for(int j = 0; j < orderN + 1; j++){  // looping over the nodes !=i
        if (j != i ){
          n = n*(x_point - xi_point(j));
          d = d*(xi_point(i) - xi_point(j));
          real_t s = 1.0;
                
          for(int N = 0; N < orderN + 1; N++){  // looping over the nodes !=i
            if (N != j && N != i ){
              s = s * (x_point - xi_point(N));
            }// end if
          }//end for
                
          p += s; 
        }//end if
         
        c = n/d; // storing a single value for interpolation for node i
        q = (p/d); // storing the derivative of the interpolating function
      } // end looping over nodes != i

        // writing value to vectors for later use
      interp(i)   = c;     // Interpolant value at given point
      Dinterp(i)  = q;     // derivative of each function
    } // end loop over all nodes

  } // end of Legrange_1D function


  // Corners of Lagrange element for mapping
  void QuadN::corners (
                       ViewCArray <real_t> &lag_nodes,   // Nodes of Lagrange elements 
                       ViewCArray <real_t> &lag_corner,  // corner nodes of QuadN element
                       const int &orderN){                 // Element order

    /*
      This image represents the corner mapping notation of an arbitrary ordered
      Lagrange element. The corner function takes in the element order and nodal positions and
      returns a vector containing the indices of the corner in alphabetical order.

      Eta
      ^
      |
      C------+-----D
      |      |     |
      |      |     |
      |      |     |
      |      ------+------> Xi   
      |            |
      |            |
      A------------B

    */

    int num_corners = 4;
    int N = orderN + 1;      //number of nodes in each direction
    int corner_ids[num_corners];

    corner_ids[0] = 0;                      
    corner_ids[1] = N - 1.0;                 
    corner_ids[2] = (N*N) - N;         
    corner_ids[3] = (N*N)-1.0;               

    for(int corner = 0; corner < num_corners; corner++){
      for(int dim = 0; dim < num_dim; dim++){
        lag_corner(corner, dim) = lag_nodes(corner_ids[corner], dim);
      }
    }


  }// end of corner mapping function


  // Functions for mapping reference position to physical position for any 
  // point in an arbitrary order 3D lagrange element
  void QuadN::physical_position (
                                 ViewCArray <real_t> &x_point,             // location in real space
                                 const ViewCArray <real_t> &lag_nodes,     // Nodes of Lagrange elements 
                                 const ViewCArray <real_t> &lag_basis_2d,  // 2D basis values 
                                 const int &orderN){                         // order of the element

    int nodes = orderN + 1;
    int Nnodes_2d = nodes * nodes;

    for (int this_vert = 0; this_vert < Nnodes_2d; this_vert++ ){
      for (int dim = 0; dim < num_dim; dim++){
        x_point(dim) += lag_nodes(this_vert, dim)*lag_basis_2d(this_vert);
      } // end for dim
    } // end for this_vert

  }// end physical position function


  void QuadN::basis_partials (
                              ViewCArray <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
                              ViewCArray <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
                              ViewCArray <real_t> &val_1d,          // Interpolant Value in 1D
                              ViewCArray <real_t> &DVal_1d,         // Derivateive of basis in 1D
                              ViewCArray <real_t> &val_2d,          // for holding the interpolant in each direction
                              ViewCArray <real_t> &DVal_2d,         // for holding the derivatives in each direction
                              ViewCArray <real_t> &lag_basis_2d,    // 3D basis values 
                              ViewCArray <real_t> &lag_partial,     // Partial of basis 
                              const ViewCArray <real_t> &xi_point,  // point of interest
                              const int &orderN){                     // Element order

    /*

      representative linear element for visualization

      Eta
      ^
      |
      2------+-----3
      |      |     |
      |      |     |
      |      |     |
      |      ------+------> Xi   
      |            |
      |            |
      0------------1


    */

    int N = orderN + 1;      //number of nodes in each direction
    int tot_pts = (N*N);     // total nodes in 2D

    real_t sumi = 0.0;
    real_t sumj = 0.0;

    //Setting nodes for Lagrange Elements
    for (int m = 0; m < tot_pts; m++) {

      int i, j;
      // sets up the i and j indices for the nodes of an 
      j = floor(m/N)+1; 
      i = (m+1) - N*(j-1);

      // xi direction
      lag_nodes(m, 0) = nodes_1d(i-1); 

      // eta direction
      lag_nodes(m, 1) = nodes_1d(j-1); 


      // calling function to assign nodal values for basis and derivative
        
      //evaluating Lagrange interpolants for each function at xi_point
      lagrange_1D(val_1d, DVal_1d, xi_point(0), nodes_1d, orderN);
      val_2d(m, 0)  = val_1d(i-1); 
      DVal_2d(m, 0) = DVal_1d(i-1);

      // resetting to zero
      for(int i = 0.0; i < N; i++){
        val_1d(i)  = 0.0;
        DVal_1d(i) = 0.0;
      }


      //evaluating Legrange interpolants for each function at xi_point
      lagrange_1D(val_1d, DVal_1d, xi_point(1), nodes_1d, orderN);
      val_2d(m, 1)  = val_1d(j-1); 
      DVal_2d(m, 1) = DVal_1d(j-1);

      // resetting to zero
      for(int i = 0.0; i < N; i++){
        val_1d(i)  = 0.0;
        DVal_1d(i) = 0.0;
      }
      // Assigning and storing the Basis
      lag_basis_2d(m) = val_2d(m, 0) * val_2d(m, 1);

      // Assigning and storing the partials
      lag_partial(m, 0)  = DVal_2d(m, 0) * val_2d(m, 1);
      lag_partial(m, 1)  = val_2d(m, 0) * DVal_2d(m, 1);
    } // end for 

  }// end basis_partials function

} // end namespace elements
