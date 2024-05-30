



// gets a particular creation matrix
SpMat get_Creation_matrix(int level, Basis b) {
  SpMat m(b.size,b.size); // output matrix
  int i,j; // counters
  bool nz; // if it is nonzero
  cout<<"###In CreaDes###\n";
  for (i=0;i<b.size;i++) { // first loop
    for (j=0;j<b.size;j++) { // second loop
      nz = b.basis[i] * create(b.basis[j],level);
      if (nz) m.insert(i,j) = 1.0;  // put one if there is overlap
    }
  }
  return m;
}




// gets a particular destruction matrix
SpMat get_Destruction_matrix(int level, Basis b) {
  SpMat m(b.size,b.size); // output matrix
  int i,j; // counters
  bool nz; // if it is nonzero
  for (i=0;i<b.size;i++) { // first loop
    for (j=0;j<b.size;j++) { // second loop
      nz = b.basis[i] * destroy(b.basis[j],level);
      if (nz) m.insert(i,j) = 1.0;  // put one if there is overlap
    }
  }
  return m;
}






// overload the CreaDes function, so that can be called
// by number and true/false
SpMat get_CreaDes_matrix(int level,bool creates, Basis b) {
  SpMat m(b.size,b.size); 
  if (creates) m = get_Creation_matrix(level,b);
  else m = get_Destruction_matrix(level,b);
  return m;
}




// gets a particular creation matrix
SpMat get_Nexchange_matrix(int level1, int level2, Basis b) {
  SpMat m(b.size,b.size); // output matrix
  Basis_vector v; // basis vector object
  v.alloc(b.basis[0].size); // allocate the vector
  int i,j; // counters
  for (i=0;i<b.size;i++) { // first loop
      v = create(destroy(b.basis[i],level1),level2); // matrix
      j = get_number(v); // matrix
      if (j!=-1)  { // if it is not the vacuum
      j = b.indexes[j]; // get the index in the reduced basis
      if (j!=-1) m.insert(j,i) = 1.0*v.sign; } // multiply by the fermionic sign
  }
  return m;
}




// gets a 4 field matrix, first destroys then creates
SpMat get_4field_matrix(int c1, int c2, int d1, int d2, Basis b) {
  SpMat m1(b.size,b.size); // output matrix
  SpMat m2(b.size,b.size); // output matrix
  SpMat mo(b.size,b.size); // output matrix
  Basis_vector v; // basis vector object
  v.alloc(b.basis[0].size); // allocate the vector
  int i,j; // counters
//  cout << c1 << "   "<< c2<<"    "<< d1<<"   "<<d2<< endl;
  for (i=0;i<b.size;i++) { // first loop
      // first combination
//      v = create(destroy(create(destroy(b.basis[i],d2),c2),d1),c1);
//      j = get_number(v);
//      if (j!=-1)  { // if it is not the vacuum
//      j = b.indexes[j]; // get the index in the reduced basis
//      if (j!=-1) {m1.insert(i,j) = 1.0*v.sign; } }
//      // second combination
//      v = create(destroy(create(destroy(b.basis[i],d1),c2),d2),c1);
//      j = get_number(v);
//      if (j!=-1)  { // if it is not the vacuum
//      j = b.indexes[j]; // get the index in the reduced basis
//      if (j!=-1) {m2.insert(i,j) = -1.0*v.sign; } }
      // third combination
      v = create(create(destroy(destroy(b.basis[i],d1),d2),c2),c1);
      j = get_number(v);
      if (j!=-1)  { // if it is not the vacuum
        j = b.indexes[j]; // get the index in the reduced basis
        if (j!=-1) m2.insert(i,j) = -1.0*v.sign;
      }
  }
  mo = m2; // sum both contributions
  return mo;
}

//// Gets an operator diagonal in the orbital
//SpMat get_diagonal_operator(Basis b,int norb, double * vals)  {
//  int i; // counter
//  SpMat m(b.size,b.size); // output matrix
//  SpMat tmpm(b.size,b.size); // temporal matrix
//  if (!(b.basis[0].size==norb)) throw 20; // check that there are incorrect
//  for (i=0;i<norb;i++) {  // loop over orbitals
//    tmpm = get_Nexchange_matrix(i,i,b); // get the orbital density















//// Gets an operator diagonal in the orbital
//SpMat get_diagonal_operator(Basis b,int norb, double * vals)  {
//  int i; // counter
//  SpMat m(b.size,b.size); // output matrix
//  SpMat tmpm(b.size,b.size); // temporal matrix
//  if (!(b.basis[0].size==norb)) throw 20; // check that there are incorrect
//  for (i=0;i<norb;i++) {  // loop over orbitals
//    tmpm = get_Nexchange_matrix(i,i,b); // get the orbital density
//    m += vals[]*tmpm; } // add contribution (up channel)
//  for (i=5;i<10;i++) {  // loop over down orbitals
//    tmpm = get_Nexchange_matrix(i,i,b); // get the orbital density
//    m += (i-7)*tmpm; } // add contribution (down channel)
//  return m; // return the matrix
//}
//
//


