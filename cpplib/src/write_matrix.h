void write_matrix(string filename, SpMat mat) { // write the matrix in a file
  const char * c = filename.c_str();  // variable for filename
  ofstream myfile; // variable for file
  myfile.open (c); // open file
  myfile << "# SIZE = ";
  myfile << mat.rows() << endl;
  myfile << "# i   j    real    imag\n";
  myfile.precision(14); // write with a certain precision
  for (int k=0; k<mat.outerSize(); ++k)
    for (SparseMatrix<complex<double> >::InnerIterator it(mat,k); it; ++it)
  {
    myfile << it.row() << "     ";   // row index
    myfile << it.col() << "     ";   // col index (here it is equal to k)
    myfile << it.value().real() << "    ";  // real part
    myfile << it.value().imag();  // imaginary part
    myfile << endl; // end the line
    it.index(); // inner index, here it is equal to it.row()
  }
  myfile.close(); // close file
}

