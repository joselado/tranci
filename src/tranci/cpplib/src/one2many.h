// transform a one body operator into a many body operator

SpMat one2many(Basis b, SpMat op)  {
  int i,j; // counter
  int snorb=b.size; // number of orbitals times spin
  complex<double> c; // coefficient to multiply
  SpMat mat(snorb,snorb); // output matrix
  SpMat tmpm(snorb,snorb); // temporal matrix
  for (int k=0; k<op.outerSize(); ++k)
  for (SparseMatrix<complex<double> >::InnerIterator it(op,k); it; ++it) {
    i = it.row();  // get the row
    j = it.col();  // get the column
    c = it.value();  // get the data
    if ((i>=snorb)or(j>=snorb)) throw 20; // check that within bounds
    tmpm = get_Nexchange_matrix(i,j,b); // move upwards the level (up)
    mat += c*tmpm;  // add contribution (up channel)
  }
  return mat; // return the matrix
}





