#include "bernstein_polynomials.h"
#include "point_distributions.h"

#include <cstdlib>
#include <iostream>

// Test Parameters
template <typename NumType>
struct TestParams {
  SizeType Np;
  SizeType N;

  NumType X;
  NumType *c;

  TestParams(SizeType order){
    Np = order;
    N = Np + 1;

  // Generating coefficients for bernstein polynomial expansion
    c = new NumType[N];
    for (SizeType i = 0; i<N; i++){
      c[i] = Real(rand())/RAND_MAX;
    }

 // Selecting a random number between -1 and 1
    NumType Zl = -1.0;
    NumType Zr = 1.0;

    X = Zl + (Zr-Zl)*Real(rand())/Real(RAND_MAX);
    };

    ~TestParams() {delete[] c;}
};

// Test 1: Comparing Bernstein Polynomials with the Bernstein approximation

bool test1(TestParams<Real> &p){
    Real tol = 1e-15;

    bool pass = false;

    Real sum1 = 0.0;
    for (SizeType i = 0; i < p.N; i++){
        Real Bi = bernstein::eval(p.N, i, p.X);
        sum1 += p.c[i]*Bi;
    };

    Real sum2 = bernstein::eval_approx(p.N, p.c, p.X);

    Real rel_error = common::abs((sum1-sum2)/sum1);

    if (rel_error < tol ){
        pass = true;
    }
    else{
        std::cout << "Test 1, error: " << rel_error << std::endl;
    }

    return pass;
};

// Test 2: Checking consistency of the first derivative of the polynomials with derivative of the Bernstein Approximation

bool test2(TestParams<Real> &p) {
    Real tol = 1e-15;

    bool pass = false;
    Real sum1 = 0.0;

    for (SizeType i = 0; i < p.N; i++) {
        Real dBi = bernstein::eval_der(p.N, i, p.X);
        sum1 += p.c[i]*dBi;
    }

    Real sum2 = bernstein::eval_der_approx(p.N, p.c, p.X);

    Real rel_error = common::abs((sum1 - sum2)/sum1);

    if (rel_error < tol) {
        pass = true;
    } 
    else {
        std::cout << "Test 2, error: " << rel_error << std::endl; 
    }

  return pass;
}

// Test 3: This test ensures the a property of the bernstein polynomials where the 
// the sum of the nth degree polynomials from v equals 0 to n-1 is equal to 1.

bool test3(TestParams<Real> &p){
    Real tol = 1e-15;
    bool pass = false;

    // Establishes random number between -1 and 1 to use for each Bernstein test
    Real Zl = -1.0;
    Real Zr = 1.0;

    Real rand_X = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);

    // Sums bernstein polynomials at value rand_X
    Real sum = 0.0;
    for (SizeType i = 0; i < p.N; i++) {
        sum += bernstein::eval(p.N, i, rand_X);
    }
    Real error = common::abs( sum - 1 );
    if (error< tol) {
        pass = true;
    }
    else {
        std::cout << "Test 3, error: " << error << std::endl;
    }

    return pass;
}

// Test 4: Checks for the positivity of the Bernstein polynomials generated by the code

bool test4(TestParams<Real> &p){

  bool pass = false;
  SizeType NumTests = 1000;
  for ( SizeType i = 0; i < NumTests + 1; i ++){

    Real Zl = -1.0;
    Real Zr = 1.0;

    Real rand_X = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);
    for (SizeType j=0; j < p.N; j++){
      Real val = bernstein::eval(p.N, j, rand_X);
      if (0 <= val){
        pass = true;
      }
      else {
        std::cout << "Test 4 Error: At least one of the values tested is negative" << std::endl;
      }
    }
  }  
  return pass;
}

bool test5(TestParams<Real> &p){
    Real tol = 1e-15;
    bool pass = false;

    Real Zl = -1.0;
    Real Zr = 1.0;
    
    Real rand_X = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);
    
    //Sums bernstein polynomial derivatives at value rand_X
    Real sum = 0.0;
    for (SizeType i = 0; i < p.N; i++) {
    sum += bernstein::eval_der(p.N, i, rand_X);
    }
    Real error = common::abs( sum );
    if (error< tol) {
      pass = true;
    }
    else {
      std::cout << "Test 3, error: " << error  << std::endl;
    }
    return pass;
}

 int main(){
    std::cout.precision(15);

    std::cout << "TEST BERNSTEIN POLYNOMIALS\n"
        <<"----------------------------"
        << std::endl;

    // Initialize Random Seed
    srand(time(NULL));

    // Generate the Parameters
    TestParams<Real> rp(8);

    // Run Test 1
    bool pass1 = test1(rp);
    std::cout << "Test 1 " << (pass1 ? "PASSED" : "FAILED") << std::endl;

    // Run Test 2
    bool pass2 = test2(rp);
    std::cout << "Test 2 " << (pass2 ? "PASSED" : "FAILED") << std::endl;

    // Run Test 3
    bool pass3 = test3(rp);
    std::cout << "Test 3 " << (pass3 ? "PASSED" : "FAILED") << std::endl;

    // Run Test 4
    bool pass4 = test4(rp);
    std::cout << "Test 4 " << (pass4 ? "PASSED" : "FAILED") << std::endl;
    bool pass5 = test5(rp);
    std::cout << "Test 5 " << (pass5 ? "PASSED" : "FAILED") << std::endl;

    std::cout << std::endl;
    std::cout << "PASSED "
              << int(pass1) + int(pass2) + int(pass3) + int(pass4) + int(pass5)
              << "/5" << std::endl;

    return 0;
    
 }


 