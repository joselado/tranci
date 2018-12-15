// gets the sign associated to getting into that level by permutation
int fermi_sign(Basis_vector v, int level) { 
  int i; // counter
  int s; // sign
  s = 1; // initialice
  for (i=0;i<level;i++) if(v.basevec[i]) s*=-1; // for each creator minus sign
  return s; // return the value
}
