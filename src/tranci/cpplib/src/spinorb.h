// generate the spin orbit coupling operator

SpMat get_ls_matrix(Basis b) {
  int n=b.basis[0].size; // number of one aprticle orbitals
  SpMat sx,sy,sz,lx,ly,lz,ls; // one particle operators
  SpMat mat(b.size,b.size); // output matrix
  int l=(int)sqrt(((double)n)/2.0); // l number
  lx = one_particle_lx(l);  // Lx
  ly = one_particle_ly(l);  // Ly
  lz = one_particle_lz(l);  // Lz
  sx = one_particle_sx(n/2);  // Sx
  sy = one_particle_sy(n/2);  // Sy
  sz = one_particle_sz(n/2);  // Sz
  ls = lx*sx + ly*sy + lz*sz; // one particle LS
  mat = one2many(b,ls); // many body LS
  return mat; // return LS matrix
}



