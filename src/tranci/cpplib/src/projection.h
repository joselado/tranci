// get the projection operators





// projection operator
SpMat get_projection_matrix(Basis b, int ii )  {
  SpMat mat(b.size,b.size); // output matrix
  int norb;
  norb = b.basis[0].size; // get the number of orbitals
  if (ii>=norb) throw 20; // check that there are 2*nm
  SpMat om(norb,norb); // one particle 
  om.insert(ii,ii) = 1.0; // insert one in that orbital
  mat = one2many(b,om); // one particle to many body
  return mat; // return the matrix
}



