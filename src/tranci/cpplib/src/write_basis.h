void write_basis(Basis b)  { // write the basis in a file
  ofstream myfile; // variable for file
  int i,j; // counters
  myfile.open("basis.out"); // open the file
  myfile << "# SIZE = "<<b.size<<endl; // size of the basis
  for (i=0;i<b.size;i++) {
    for (j=0;j<b.basis[i].size;j++) myfile << b.basis[i].basevec[j] << " ";
    myfile << "      #" << i; // new line
    myfile << endl; // new line
  }
  myfile.close(); // close file
}
