
// generate the basis object for a certain number of levels
Basis generate_basis(int nlvl) { // number of electrons and number of levels
  int nwf; // number of wavefunctions
  int i,j; // counters
  bool lval;
  Basis obasis; // output basis

  
  nwf = (int)pow(2,nlvl) ; // total number of wavefunctions
  obasis.alloc(nwf);  // allocate basis
  
  for (j=0;j<nwf;j++) { // loop over wavefunctions
    obasis.basis[j].alloc(nlvl);  // allocate wavefunction
    for (i = 0; i < nlvl; ++i) {  // assuming a 32 bit int
        lval =  (j & (1 << i) ? 1 : 0); // logical value
        obasis.basis[j].basevec[i] = (lval!=0); // store in the vector
    };
  };

  return obasis;  
}



