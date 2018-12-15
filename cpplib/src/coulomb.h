

// re4ad the coulomb matrix element from file
double ****read_coulomb_d() { // read the coulomb matrix
  int i,j,k,l; // indexes
  double z; // coulomb matrix elements (read from file)
  double ****v; // coulomb matrix elements (read from file)
  v = new double*** [5];
  for (i=0;i<5;i++) {
    v[i] = new double** [5];
    for (j=0;j<5;j++) {
      v[i][j] = new double* [5];
      for (k=0;k<5;k++) {
        v[i][j][k] = new double [5];
        for (l=0;l<5;l++) v[i][j][k][l] = 0.0; // set to zero
      }
    }
  }
  double d; // holder for 
  int nl,il; //number of lines
  ifstream myfile; // file to read
  myfile.open("Vc.in"); // open the file
  myfile >> nl; // read number of lines  
  for (il=0;il<nl;il++) {
    myfile >> i >> j >> k >> l >> d; // read the data
    z = d; // convert to complex
//    cout << i << "    " << j<< "   " << k<<"    " << l<<"   " << d << endl;
    i -=1;j-=1;k-=1;l-=1; // counters start in 0
    v[i][j][k][l] = z; // store in array
  }
  myfile.close(); // close the file
  return v; // return the pointer
}












// gets the density matrix
SpMat get_coulomb_matrix(Basis b) {
  SpMat m(b.size,b.size); // output matrix
  int i,j,k,l; // counters
  int norb=b.basis[0].size; // number of orbitals
  int sln = norb/2; // spinless number of orbitals
  double **** v; // coulomb matrix elements (read from file)
  double vtmp; // coulomb matrix elements (read from file)
  v = read_coulomb_d(); // read coulomb matrix elements
  for (i=0;i<norb/2;i++) { // first loop
    for (j=0;j<norb/2;j++) { // second loop
      for (k=0;k<norb/2;k++) { // third loop
        for (l=0;l<norb/2;l++) { // fourth loop
        vtmp = v[i][j][k][l]; // get the matrix
        if (vtmp!=0.0) 
        m += (complex<double>)vtmp*get_4field_matrix(i,j,k,l,b); // get the matrix
        m += (complex<double>)vtmp*get_4field_matrix(i+sln,j,k,l+sln,b); // get the matrix
//        m += (complex<double>)vtmp*get_4field_matrix(i+sln,j,k+sln,l,b); // get the matrix
        m += (complex<double>)vtmp*get_4field_matrix(i,j+sln,k+sln,l,b); // get the matrix
//        m += (complex<double>)vtmp*get_4field_matrix(i,j+sln,k,l+sln,b); // get the matrix
        m += (complex<double>)vtmp*get_4field_matrix(i+sln,j+sln,k+sln,l+sln,b); // get the matrix
//        cout << i<<j<<k<<l<<m<<norb<<endl;
//        m += (complex<double>)vtmp*get_4field_matrix(i,j,k,l,b); // get the matrix
        }
      }
    }
  }
  m.makeCompressed(); // remove zeros
  return m;
}

