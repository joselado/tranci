// dependencies
// objects.h


SpMat one_particle_lp(int l)  { // calculates the one particle Lp matrix
  int nm=2*l+1; // number of m values
  int n=2*nm; // total dimension
  double c; // factor to multiply the exchange
  int i; // counter
  double m; // m number
  SpMat mat(n,n); // output matrix
  for (i=0;i<(nm-1);i++) {  // loop over orbitals up orbitals
    m = (double)(i-l); // m number 
    c = sqrt(((double)l-m)*((double)l+m+1)); // coefficient
    mat.insert(i,i+1) = c;  // add contribution (up channel)
    mat.insert(i+nm,i+1+nm) = c;  // add contribution (down channel)
  }
  return mat;
}


SpMat one_particle_lm(int l)  { // calculates the one particle Lp matrix
  SpMat mat; // output matrix
  mat = one_particle_lp(l); // get Lp
  mat = mat.adjoint(); // adjoint the matrix
  return mat; 
}




SpMat one_particle_lx(int l)  { // calculates the one particle Lx matrix
  int n=2*(2*l+1);
  SpMat mat(n,n); // output matrix
  SpMat lp(n,n); // output matrix
  SpMat lm(n,n); // output matrix
  lp = one_particle_lp(l); // get Lp
  lm = one_particle_lm(l); // get Lp
  mat = (lp + lm)/2.0; // adjoint the matrix
  return mat; 
}




SpMat one_particle_ly(int l)  { // calculates the one particle Lx matrix
  int n=2*(2*l+1);
  SpMat mat(n,n); // output matrix
  SpMat lp(n,n); // output matrix
  SpMat lm(n,n); // output matrix
  lp = one_particle_lp(l); // get Lp
  lm = one_particle_lm(l); // get Lp
  mat = -im*(lp - lm)/2.0; // adjoint the matrix
  return mat; 
}




SpMat one_particle_lz(int l)  { // calculates the one particle Lp matrix
  int nm=2*l+1; // number of m values
  int n=2*nm; // total dimension
  double c; // factor to multiply the exchange
  int i; // counter
  double m; // m number
  SpMat mat(n,n); // output matrix
  for (i=0;i<nm;i++) {  // loop over orbitals
    m = (double)(i-l); // m number 
    c = m; // coefficient
    mat.insert(i,i) = c;  // add contribution (up channel)
    mat.insert(i+nm,i+nm) = c;  // add contribution (down channel)
  }
  return mat;
}




//
//
//// L+ operator
//SpMat get_lz_matrix(Basis b)  {
//  int i; // counter
//  complex<double> c; // coefficient to multiply
//  double l,m,nm; // angular number
//  double guessl; // guess the angular number
//  SpMat mat(b.size,b.size); // output matrix
//  SpMat tmpm(b.size,b.size); // temporal matrix
//  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
//  l = guessl; // for dlevels
//  nm = 2*l+1; // number of m levels
//  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
//  for (i=0;i<(int)nm;i++) {  // loop over orbitals up orbitals
//    m = (double)i-l; // m number 
//    tmpm = get_Nexchange_matrix(i,i,b); // move upwards the level (up)
//    mat += m*tmpm;  // add contribution (up channel)
//    tmpm = get_Nexchange_matrix(i+(int)nm,i+(int)nm,b); // move upwards the level (down)
//    mat += m*tmpm; } // add contribution (down channel)
//  return mat; // return the matrix
//}
//
//








// L+ operator
SpMat get_lz_matrix(Basis b)  {
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oplz; // one particle Lz
  oplz = one_particle_lz((int)l); // get the one particle operator
  mat = one2many(b,oplz); // one particle to many body
  return mat; // return the matrix
}


// Lx operator
SpMat get_lx_matrix(Basis b)  {
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oplz; // one particle Lz
  oplz = one_particle_lx((int)l); // get the one particle operator
  mat = one2many(b,oplz); // one particle to many body
  return mat; // return the matrix
}




// L+ operator
SpMat get_ly_matrix(Basis b)  {
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oplz; // one particle Lz
  oplz = one_particle_ly((int)l); // get the one particle operator
  mat = one2many(b,oplz); // one particle to many body
  return mat; // return the matrix
}





// Lz to a certain power
SpMat get_zn_matrix(Basis b, int nexp)  {
  int i; // counter
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oplz,oplzn; // one particle Lz
  oplz = one_particle_lz((int)l); // get the one particle operator
  oplzn = one_particle_lz((int)l); // get the one particle operator
  for (i=0;i<(nexp-1);i++) oplzn = oplzn * oplz; // perform the power
  mat = one2many(b,oplzn); // one particle to many body
  return mat; // return the matrix
}





// Lx to a certain power
SpMat get_xn_matrix(Basis b, int nexp)  {
  int i; // counter
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oplx,oplxn; // one particle Lz
  oplx = one_particle_lx((int)l); // get the one particle operator
  oplxn = one_particle_lx((int)l); // get the one particle operator
  for (i=0;i<(nexp-1);i++) oplxn = oplxn * oplx; // perform the power
  mat = one2many(b,oplxn); // one particle to many body
  return mat; // return the matrix
}





// Lx to a certain power
SpMat get_xnyn_matrix(Basis b, int nexp)  {
  int i; // counter
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oplx,oplxn; // one particle Lz
  oplx = one_particle_lx((int)l)*one_particle_ly((int)l); // get the one particle operator
  oplxn = one_particle_lx((int)l)*one_particle_ly((int)l); // get the one particle operator
  for (i=0;i<(nexp-1);i++) oplxn = oplxn * oplx; // perform the power
  SpMat m2; // one particle Lz
  m2 = oplxn.adjoint();    // add the hermitian
  m2 += oplxn ;
  mat = one2many(b,m2); // one particle to many body
  return mat; // return the matrix
}














// Ly to a certain power
SpMat get_yn_matrix(Basis b, int nexp)  {
  int i; // counter
  double l,nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  l = guessl; // for dlevels
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==2*nm)) throw 20; // check that there are 2*nm
  SpMat oply,oplyn; // one particle Lz
  oply = one_particle_ly((int)l); // get the one particle operator
  oplyn = one_particle_ly((int)l); // get the one particle operator
  for (i=0;i<(nexp-1);i++) oplyn = oplyn * oply; // perform the power
  mat = one2many(b,oplyn); // one particle to many body
  return mat; // return the matrix
}












//
//// L+ operator
//SpMat get_lp_matrix(Basis b)  {
//  int i; // counter
//  complex<double> c; // coefficient to multiply
//  double l,m,nm; // angular number
//  double guessl; // guess the angular number
//  SpMat mat(b.size,b.size); // output matrix
//  SpMat tmpm(b.size,b.size); // temporal matrix
//  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
//  l = guessl; // for dlevels
//  nm = 2*l+1; // number of m levels
//  if (!(b.basis[0].size==(int)2*nm)) throw 20; // check that there are 2*nm
//  for (i=0;i<((int)nm-1);i++) {  // loop over orbitals up orbitals
//    m = (double)(i-l); // m number 
//    c = sqrt((l-m)*(l+m+1)); // coefficient
//    tmpm = get_Nexchange_matrix(i,i+1,b); // move upwards the level (up)
//    mat += c*tmpm;  // add contribution (up channel)
//    tmpm = get_Nexchange_matrix(i+(int)nm,i+(int)nm+1,b); // move upwards the level (down)
//    mat += c*tmpm; } // add contribution (down channel)
//  return mat; // return the matrix
//}
//
//




// L+ operator
SpMat get_lp_matrix(Basis b)  {
  complex<double> c; // coefficient to multiply
  double nm; // angular number
  double guessl; // guess the angular number
  SpMat mat(b.size,b.size); // output matrix
  SpMat tmpm(b.size,b.size); // temporal matrix
  guessl = (double)(int)sqrt((double)(b.basis[0].size)/2); // square root
  int l; // l number
  l = (int)guessl; // 
  SpMat oplp; // one particle L+
  oplp = one_particle_lp(l); // get the one particle L+ operator
  nm = 2*l+1; // number of m levels
  if (!(b.basis[0].size==(int)2*nm)) throw 20; // check that there are 2*nm
  mat = one2many(b,oplp); // get the many body operator
  return mat; // return the matrix
}


// S- matrix
SpMat get_lm_matrix(Basis b)  {
  return  get_lp_matrix(b).adjoint(); // S+ matrix is real, no need to conj
}


// get operator L**2
SpMat get_l2_matrix(Basis b)  {
  SpMat mat(b.size,b.size); // output matrix
  SpMat tmpmat(b.size,b.size); // temporal matrix
  tmpmat = get_lp_matrix(b); // get L+
  mat = tmpmat*tmpmat.transpose(); // L+L- 
  mat += tmpmat.transpose()*tmpmat; // L-L+ 
  tmpmat = get_lz_matrix(b); // get Lz
  mat /= 2; // divide by 1/2
  mat += tmpmat*tmpmat; // L-L+ 
  return mat;
}

//
//
//// Sx matrix
//SpMat get_lx_matrix(Basis b)  {
//  SpMat mat2(b.size,b.size); // output matrix
//  SpMat mat(b.size,b.size); // output matrix
//  mat = get_lp_matrix(b); // S+ matrix is real, no need to conj
//  mat2 = mat.transpose();
//  mat = (mat + mat2)/2.0;
//  return mat; // return matrix
//}
//
//
//
//
//
//
//// Sy matrix
//SpMat get_ly_matrix(Basis b)  {
//   // complex unit
//  complex<double> im;
//  im = -1.0;
//  im = sqrt(im);
//  SpMat mat2(b.size,b.size); // output matrix
//  SpMat mat(b.size,b.size); // output matrix
//  mat = get_lp_matrix(b); // S+ matrix is real, no need to conj
//  mat2 = mat.transpose();
//  mat = -im*(mat - mat2)/2.0;
//  return mat; // return matrix
//}
//





