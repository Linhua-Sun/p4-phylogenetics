#include <iostream>
#include <cmath>

// This just shows that it has a fairly low upper limit in 'n', the
// size of the input tree.  Depends on beta, but seems to be perhaps
// only 150 or so.  Can this be loggified?

// Build using g++, then run the a.out

double bForN(int n)
{
  double prod = 1.;
  if (n > 3) {
    for (int k=4; k<=n; k++) {
      prod *= (double)((2 * k) - 5);
    }
  }
  return prod;
}

double BS2009_Eqn30_ZTApprox(int n, double beta, int cT) 
{
  double myLambda;
  double tester;
  double epsilon;
  double bigANEpsilon;
  double termA;
  double termB;
  double dn;
  double dcT;

  dn = (double)n;
  dcT = (double)cT;
  myLambda = dcT/(2.0 * dn);
  tester = 0.5 * std::log((dn - 3.0)/myLambda);

  epsilon = std::exp(-2.0 * beta);
  bigANEpsilon = 1. + (((2. * dn) - 3.) * epsilon) + (2. * ((dn * dn) - (4. * dn) - 6.) * epsilon * epsilon);
  termA = bigANEpsilon + 6. * dcT * epsilon * epsilon;

  if (beta < tester) {
    termB = std::exp(-(2. * beta) * (dn - 3.) + (myLambda * (std::exp(2. * beta) - 1.)));
    termB *= bForN(n);
    if (termA > termB) {
      return termA;
    } else {
      return termB;
    }
  } else {
    return termA;
  }
}


int main ()
{
  int i;
  double ret;

  // adjust these limits to test
  for (i=10; i<150; i++) {
    // adjust beta and fake a cT to test
    ret = BS2009_Eqn30_ZTApprox(i, 0.5, i/4);
    std::cout << i << "    " << ret << std::endl;
  }
  return 0;
}

