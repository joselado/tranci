
// applies the creation operator to the vector
Basis_vector destroy(Basis_vector vin,int level) {
  Basis_vector vout; // output vector
  int i; // counter
  bool nonex; // non existent level
  vout.alloc(vin.size); // allocate output vector
  if (vin.is_zero) {  // if input vector is the zero vector
    vout.is_zero = true;
    return vout;}
  nonex = ((level+1)>vin.size) or (level<0);
  if (nonex) {  // if create in a level which does not exist
    vout.is_zero = true;
    return vout;}
  else {  // if input is nonzero
    if (!vin.basevec[level]) { // if occupied level, return zero
      vout.is_zero = true; // is the zero vector
      return vout; // return the vector
    }
    else {  // if unoccupied, copy the vector and occupy
      for (i=0;i<vin.size;i++)  {
       if (vin.basevec[i]) vout.basevec[i] = true; // copy
       else vout.basevec[i] = false; // copy
      }
      vout.basevec[level] = false; // unoccupy level
      vout.sign = vin.sign*fermi_sign(vin,level); // multiply by the fermionic sign 
      return vout;
    }
  } // close nonzero input
}
