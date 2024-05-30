// this subroutine has as input a basis and
// outputs another basis with a certain number of electrons
//#include"objects.h"  // objects for the CI calculation
Basis get_basis_nelectrons(Basis bin, int nel) {
Basis bout; // output basis
//Basis bout; // output basis
int i,j; //counters
int dim; //length of wavevectors
int nout; //number of output wavevectors
int tmpne; // temporal number of electrons
int desired[bin.size];// array for the desired wavectors 

dim = bin.basis[0].size; // length of each basis element
if (nel>dim) throw 20; // throw error if incorrect number of electrons
nout = 0;
for (i=0;i<bin.size;i++) { // loop over input basis
  tmpne = 0; // initialice the number of electrons
  for (j=0;j<dim;j++) {
    if (bin.basis[i].basevec[j]) tmpne += 1; // increase number of electrons
  };
  if (tmpne==nel) {
    desired[nout] = i; // store which vector is the good one
    nout +=1; // increase number of output vectors
  };
};
// now save in the vector which will be the true output
bout.alloc(nout); // allocate true dimension of the full basis
for (i=0;i<nout;i++) { // loop over input basis
  bout.basis[i].alloc(dim); // allocate each element of the basis
  for (j=0;j<dim;j++) {  // loop over components
    bout.basis[i].basevec[j] = bin.basis[desired[i]].basevec[j]; // store
  };
};

return bout; // return the cropped basis

}

