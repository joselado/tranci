#include"fermi_sign.h"
#include "create.h"  // function that creates a particle
#include "destroy.h"  // function that destroys a particle

// overload different operators

// multiply two many body state, yields 0 or nonzero
bool operator*(Basis_vector a, Basis_vector b) { 
  int i; // counter
  if (a.is_zero) return false; // if a is zero
  if (b.is_zero) return false; // if b is zero
  for (i=0;i<a.size;i++) if (!(a.basevec[i]==b.basevec[i])) return false;
  return true;
}

// multiply a many body state by a boolean
Basis_vector operator*(Basis_vector vin, bool ll) {
  Basis_vector vo; // output vector
  vo = vin.copy(); // create a copy of the vector
  if (!ll)  {vo.is_zero = true; return vo; } // if multiply by zero
  else return vo; // return a copy if boolean is true
} 


Basis_vector operator*(bool ll,Basis_vector vin) {
  return vin*ll; // call the previous function
}







// equal operator for basis vectors
bool operator==(Basis_vector a, Basis_vector b) { 
  int i; // counter
  if (a.is_zero or b.is_zero) return false; // if a is zero
  for (i=0;i<a.size;i++) if (!(a.basevec[i]==b.basevec[i])) return false;
  return true;
}






