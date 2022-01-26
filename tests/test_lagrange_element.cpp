#include "element_types/lagrange_element.h"
#include "element_types/point_distributions.h"

#include <cstdlib>   // rand, RAND_MAX
#include <ctime>     // time
#include <iostream>

// There are 3 tests that cover the functionality of the Lagrange element
// routines. These are described below

// Note: in the spirit of intellectual honesty, it should be disclaimed that
// the tolerances chosen in these tests were chosen so that the tests would
// pass when the results are "pretty damn close"; that is, the discrepancies
// between the results and the expectation or truth are presumed to come from
// finite precision errors. That may or may not be good enough. If in the
// course of using these routines elsewhere their numerical stabililty is
// called into question, these discrepancies and their sources may need to be
// investigated.

/* Test parameters */
template <typename NumType>
struct TestParams {
  LagrangeElement<NumType> *elem;

  NumType *c;

  NumType X[3];
  NumType Xv[3];

  TestParams(SizeType order) {
    // Test parameters
    SizeType Np = order;
    SizeType N = order + 1;

    // Generate a set of points between -1 and 1
    NumType Zl = -1.0;
    NumType Zr = 1.0;
    NumType *Z = new NumType[N];
    equispaced_points(N, Zl, Zr, Z);

    // Create Lagrange element
    elem = new LagrangeElement<NumType>(Np, Z);

    // Generate random array of coefficients between 0 and 1
    c = new NumType[elem->Ne];
    for (SizeType i = 0; i < elem->Ne; i++) {
      c[i] = Real(rand())/RAND_MAX;
    }

    // Select coordinates between -1 and 1
    X[0] = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);
    X[1] = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);
    X[2] = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);

    // Select a random vertex in the element
    SizeType I = std::round((N - 1)*Real(rand())/Real(RAND_MAX));
    SizeType J = std::round((N - 1)*Real(rand())/Real(RAND_MAX));
    SizeType K = std::round((N - 1)*Real(rand())/Real(RAND_MAX));

    Xv[0] = Z[I];
    Xv[1] = Z[J];
    Xv[2] = Z[K];
  };

  ~TestParams() { delete elem; }
};

/*
 * Test 1 
 * ------
 * This test checks the consistency of Lagrange tensor-product basis functions
 * with the Lagrange tensor-product interpolant. There are two ways to evaluate
 * the sum of the products of Lagrange tensor-product basis functions and
 * coefficients. One is to evaluate the basis functions individually, scale
 * them by the corresponding coefficients, and sum. The other is to use the
 * interpolation routine. Either way you do it, you should get the same answer;
 * that is, they should be consistent.
 *
 * Additionally, there are two cases of interest when using the Lagrange
 * polynomials in barycentric form, as we are here. The first case is when the
 * input coordinate is not coincident with any of the nodes. The second case is
 * when the input coordinate is coincident, which causes singularity to arise
 * in the barycentric form and must be handled differently.
 */
bool test1(TestParams<Real> &p) {
  Real tol = 1e-14;

  // Case 1: random coordinates, not coincident with any vertices
  bool case1 = false;
  Real sum1 = 0.0;

  for (SizeType i = 0; i < p.elem->Ne; i++) {
    Real phi_i = p.elem->eval_basis(i, p.X);
    sum1 += p.c[i]*phi_i;
  }

  Real sum2 = p.elem->eval_approx(p.c, p.X);

  Real rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    case1 = true;
  } else {
    std::cout << "Test 1, Case 1, error: " << rel_error << std::endl; 
  }

  // Case 2: random vertex, coincident case
  bool case2 = false;
  sum1 = 0.0;

  for (SizeType i = 0; i < p.elem->Ne; i++) {
    Real phi_i = p.elem->eval_basis(i, p.Xv);
    sum1 += p.c[i]*phi_i;
  }

  sum2 = p.elem->eval_approx(p.c, p.Xv);

  rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    case2 = true;
  } else {
    std::cout << "Test 1, Case 2, error: " << rel_error << std::endl; 
  }

  return case1 && case2;
}

/*
 * Test 2
 * ------
 * This test checks the consistency of the (first) derivatives of the Lagrange
 * polynomials with the derivative of the Lagrange interpolant. As in Test 1,
 * there are two ways to obtain the derivative of the sum of the products of
 * the Lagrange polynomials and coefficients. The results from both should be
 * consistent.
 *
 * The test addresses the same two cases as in Test 1.
 */
bool test2(TestParams<Real> &p) {
  Real tol = 1e-10;

  // Case 1: random coordinate, not coincident with any vertices
  bool case1 = false;

  Real d0f = 0.0;
  Real d1f = 0.0;
  Real d2f = 0.0;

  for (SizeType i = 0; i < p.elem->Ne; i++) {
    Real grad_phi[3];
    p.elem->eval_grad_basis(i, p.X, grad_phi);
    d0f += p.c[i]*grad_phi[0];
    d1f += p.c[i]*grad_phi[1];
    d2f += p.c[i]*grad_phi[2];
  }

  Real grad_f[3];
  p.elem->eval_grad_approx(p.c, p.X, grad_f);

  Real rel_error0 = common::abs((d0f - grad_f[0])/d0f);
  Real rel_error1 = common::abs((d1f - grad_f[1])/d1f);
  Real rel_error2 = common::abs((d2f - grad_f[2])/d2f);

  if (rel_error0 < tol && rel_error1 < tol && rel_error2 < tol) {
    case1 = true;
  } else {
    std::cout << "Test 2, Case 1, errors: " 
              << rel_error0 << " " 
              << rel_error1 << " " 
              << rel_error2 << " " 
              << std::endl; 
  }

  // Case 2: random vertex, coincident case
  bool case2 = false;

  d0f = 0.0;
  d1f = 0.0;
  d2f = 0.0;

  for (SizeType i = 0; i < p.elem->Ne; i++) {
    Real grad_phi[3];
    p.elem->eval_grad_basis(i, p.Xv, grad_phi);
    d0f += p.c[i]*grad_phi[0];
    d1f += p.c[i]*grad_phi[1];
    d2f += p.c[i]*grad_phi[2];
  }

  p.elem->eval_grad_approx(p.c, p.Xv, grad_f);

  rel_error0 = common::abs((d0f - grad_f[0])/d0f);
  rel_error1 = common::abs((d1f - grad_f[1])/d1f);
  rel_error2 = common::abs((d2f - grad_f[2])/d2f);

  if (rel_error0 < tol && rel_error1 < tol && rel_error2 < tol) {
    case2 = true;
  } else {
    std::cout << "Test 2, Case 2, errors: " 
              << rel_error0 << " " 
              << rel_error1 << " " 
              << rel_error2 << " " 
              << std::endl; 
  }

  return case1 && case2;
}

/*
 * Test 3
 * ------
 * This test checks the correctness of the basis gradient routine against an
 * approximation of the derivative based on the complex step method.
 *
 * This test does not check the pathological coincident vertex case. There are
 * two reasons for this. First, it would be impractical to do so. Second,
 * assuming the Lagrange polynomial tests are passing, it should be sufficient
 * to check the non-coincident case here, since the way the Lagrange polynomial
 * routines are called from the Lagrange element routine are the same in either
 * case.
 */
bool test3(TestParams<Complex> &p) {
  Real tol = 1e-10;
  Real h = 1e-30;

  // Case 1: random coordinate, not coincident with any vertices
  bool pass = false;

  Complex Xc[3];
  Complex phi;

  std::copy(p.X, p.X+3, Xc);
  Xc[0] += Complex(0, h);
  phi = p.elem->eval_basis(0, Xc);
  Real d0phi = common::imag(phi)/h;

  std::copy(p.X, p.X+3, Xc);
  Xc[1] += Complex(0, h);
  phi = p.elem->eval_basis(0, Xc);
  Real d1phi = common::imag(phi)/h;

  std::copy(p.X, p.X+3, Xc);
  Xc[2] += Complex(0, h);
  phi = p.elem->eval_basis(0, Xc);
  Real d2phi = common::imag(phi)/h;

  Complex grad_phi[3];
  p.elem->eval_grad_basis(0, p.X, grad_phi);

  Real rel_error0 = common::abs((grad_phi[0] - d0phi)/d0phi);
  Real rel_error1 = common::abs((grad_phi[1] - d1phi)/d1phi);
  Real rel_error2 = common::abs((grad_phi[2] - d2phi)/d2phi);

  if (rel_error0 < tol && rel_error1 < tol && rel_error2 < tol) {
    pass = true;
  } else {
    std::cout << "Test 3, errors: " 
              << rel_error0 << " " 
              << rel_error1 << " " 
              << rel_error2 << " " 
              << std::endl; 
  }

  return pass;
}

int main() {
  std::cout.precision(15);

  std::cout << "TEST LAGRANGE ELEMENT\n" 
            << "---------------------" 
            << std::endl;

  // Initialize random seed
  srand(time(NULL));

  // Generate real-valued test parameters
  TestParams<Real> rp(8);

  // Run test 1
  bool pass1 = test1(rp);
  std::cout << "TEST 1 " << (pass1 ? "PASSED" : "FAILED") << std::endl;

  // Run test 2
  bool pass2 = test2(rp);
  std::cout << "TEST 2 " << (pass2 ? "PASSED" : "FAILED") << std::endl;

  // Generate complex-valued test parameters
  TestParams<Complex> ip(8);

  // Run test 3
  bool pass3 = test3(ip);
  std::cout << "TEST 3 " << (pass3 ? "PASSED" : "FAILED") << std::endl;

  std::cout << std::endl;
  std::cout << "PASSED " 
            << int(pass1) + int(pass2) + int(pass3) 
            << "/3" << std::endl;
  
  return 0;
}
