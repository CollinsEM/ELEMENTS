#pragma once

#include "element_types/jacobi_polynomials.h"

#include "common/matar_blas_lapack_interface.h"
#include "common/error.h"
#include "common/common.h"

namespace elements {

  /**
   * The Gauss-Jacobi (also Jacobi or Mehler) quadrature rule is a 1D
   * quadrature rule over the interval: $x\in[-1,1]$ with the weighting
   * function:
   *
   * \[
   *   W(x) = (1-x)^\alpha * (1+x)^\beta
   * \]
   *
   * for $\alpha,\beta>-1$. For of a quadrature rule with
   * order-of-accuracy $2n$, the $n$ abscissas are given by the roots
   * of the Jacobi polynomials $P_n^{alpha,beta}(x)$.
   *
   * Guass-Jacobi quadrature can be used to approximate integrals with
   * singularities at the end points of the interval.
   *
   * The Gauss-Legendre quadrature rule can be recovered from:
   * $\alpha=0, \beta=0$.
   *
   * Gauss-Gegenbauer quadrature rules can be obtained by setting
   * $\alpha=\beta$, in which case the Jacobi polynomials become
   * Gegenbauer polynomials.
   *
   * Gauss-Chebyshev quadrature of the first kind can be obtained from
   * $\alpha=-0.5, \beta=-0.5$, and of the second kind from:
   * $\alpha= 0.5, \beta= 0.5$.
   */
  struct GaussJacobiQuadrature {
    SizeType N;
    CArray<Real> points;
    CArray<Real> weights;
    /**
     * @param n     Number of quadrature points generated
     * @param alpha Jacobi polynomial parameter
     * @param beta  Jacobi polynomial parameter
     */
    GaussJacobiQuadrature( SizeType n, Real alpha, Real beta );
  };

} // end namespace elements
