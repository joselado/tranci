// library to calculate matrix elements involving the spin operator
// Sz operator




SpMat one_particle_sz(int nm)  { // calculates the one particle Sz matrix
  int n=2*nm; // total dimension
  int i; // counter
  SpMat mat(n,n); // output matrix
  for (i=0;i<nm;i++) {  // loop over orbitals up orbitals
    mat.insert(i,i) = (complex<double>)0.5;  // add contribution (up channel)
    mat.insert(i+nm,i+nm) = (complex<double>)-0.5;  // add contribution (down channel)
  }
  return mat;
}





SpMat one_particle_sp(int nm)  { // calculates the one particle Lp matrix
  int n=2*nm; // total dimension
  int i; // counter
  SpMat mat(n,n); // output matrix
  for (i=0;i<nm;i++) {  // loop over orbitals up orbitals
    mat.insert(i+nm,i) = (complex<double>)1.0;  // add contribution (up channel)
  }
  return mat;
}


// one particle S- matrix
SpMat one_particle_sm(int nm)  { // calculates the one particle Lp matrix
  SpMat mat(nm*2,nm*2); // output matrix
  mat = one_particle_sp(nm);
  mat = mat.adjoint();
  return mat;
}


SpMat one_particle_sx(int nm)  { // calculates the one particle Lx matrix
  int n=nm*2; // dimension of the matrix 
  SpMat mat(n,n); // output matrix
  SpMat sp(n,n); // output matrix
  SpMat sm(n,n); // output matrix
  sp = one_particle_sp(nm); // get Lp
  sm = one_particle_sm(nm); // get Lp
  mat = (sp + sm)/2.0; // adjoint the matrix
  return mat;
}




SpMat one_particle_sy(int nm)  { // calculates the one particle Lx matrix
  SpMat mat; // output matrix
  SpMat sp; // output matrix
  SpMat sm; // output matrix
  sp = one_particle_sp(nm); // get Lp
  sm = one_particle_sm(nm); // get Lp
  mat = -im*(sp - sm)/2.0; // adjoint the matrix
  return mat;
}






SpMat get_sz_matrix(Basis b)  {
  int snorb; // number of orbitals (times spin)
  int norb; // number of orbitals (without spin)
  SpMat mat(b.size,b.size); // output matrix
  snorb = b.basis[0].size; // number of orbitals
  if (!(snorb%2==0)) throw 20; // check that there are even number of orbitals
  norb = snorb/2; // orbitals without spin
  SpMat opsz; // one particle Sz
  opsz = one_particle_sz(norb); // get the one particle operator
  mat = one2many(b,opsz); // one particle to many body
  return mat; // return the matrix
}




SpMat get_sx_matrix(Basis b)  {
  int snorb; // number of orbitals (times spin)
  int norb; // number of orbitals (without spin)
  SpMat mat(b.size,b.size); // output matrix
  snorb = b.basis[0].size; // number of orbitals
  if (!(snorb%2==0)) throw 20; // check that there are even number of orbitals
  norb = snorb/2; // orbitals without spin
  SpMat opsz; // one particle Sz
  opsz = one_particle_sx(norb); // get the one particle operator
  mat = one2many(b,opsz); // one particle to many body
  return mat; // return the matrix
}





SpMat get_sy_matrix(Basis b)  {
  int snorb; // number of orbitals (times spin)
  int norb; // number of orbitals (without spin)
  SpMat mat(b.size,b.size); // output matrix
  snorb = b.basis[0].size; // number of orbitals
  if (!(snorb%2==0)) throw 20; // check that there are even number of orbitals
  norb = snorb/2; // orbitals without spin
  SpMat opsz; // one particle Sz
  opsz = one_particle_sy(norb); // get the one particle operator
  mat = one2many(b,opsz); // one particle to many body
  return mat; // return the matrix
}










