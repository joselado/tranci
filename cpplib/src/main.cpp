#include "ci.h" // ci library




int main() {
//Basis bin; // basis variable
Basis bin,MB_basis; // basis variable
int norb; // number of obitals
norb = 10;
bin = generate_basis(norb); // generate all-electron basis
cout << "Begin print basis\n";


//bin.print() ;
//cout << "End print basis\n";
ifstream myfile; // file to read
myfile.open("nelectrons.in"); // open the file
int ne; // number of electrons
myfile >> ne;
myfile.close();
MB_basis = get_basis_nelectrons(bin,ne);
MB_basis.get_indexes(); // generate the indexes
//for (i=0;i<MB_basis.size;i++) MB_basis.basis[i] = create(MB_basis.basis[i],0);
//MB_basis.print() ;
//cout << get_CreaDes_matrix(0,true,MB_basis);
//cout << get_Nexchange_matrix(0,1,bin);
//////////////////////////
// define several matrices
//////////////////////////
int d; // dimension of the basis
d = MB_basis.size; // dimension
SpMat lx(d,d); // Lx
SpMat ly(d,d); // Ly
SpMat lz(d,d); // Lz
SpMat x4(d,d); // Lx2
SpMat y4(d,d); // Ly2
SpMat z4(d,d); // Lz2
SpMat x2y2(d,d); // Lx2
SpMat x2(d,d); // Lx2
SpMat y2(d,d); // Ly2
SpMat z2(d,d); // Lz2
SpMat sx(d,d); // Sx
SpMat sy(d,d); // Sy
SpMat sz(d,d); // Sz
SpMat ls(d,d); // LS
SpMat jx(d,d); // Jx
SpMat jy(d,d); // Jy
SpMat jz(d,d); // Jz
SpMat j2(d,d); // Jz
SpMat vc(d,d); // Coulomb matrix
sx = get_sx_matrix(MB_basis); // get Sx term
sy = get_sy_matrix(MB_basis); // get Sy term
sz = get_sz_matrix(MB_basis); // get Sz term
lx = get_lx_matrix(MB_basis); // get Lx term
ly = get_ly_matrix(MB_basis); // get Ly term
lz = get_lz_matrix(MB_basis); // get Lz term
x2 = get_xn_matrix(MB_basis,2); // get Lx2 term
y2 = get_yn_matrix(MB_basis,2); // get Ly2 term
z2 = get_zn_matrix(MB_basis,2); // get Lz2 term
x2y2 = get_xnyn_matrix(MB_basis,2); // get Lx2 term
x4 = get_xn_matrix(MB_basis,4); // get Lx2 term
y4 = get_yn_matrix(MB_basis,4); // get Ly2 term
z4 = get_zn_matrix(MB_basis,4); // get Lz2 term
ls = get_ls_matrix(MB_basis); // get Lz term
jx = sx + lx; // Jx term
jy = sy + ly; // Jy term
jz = sz + lz; // Jz term
j2 = jx*jx + jy*jy + jz*jz; // J^2 term
vc = get_coulomb_matrix(MB_basis);
//////////////////////////
// projection operators
SpMat * proji; // projections
proji = new SpMat [norb]; // allocate matrices
int iproj; // index
 // generate all the projections 
for (iproj=0;iproj<norb;iproj++) {proji[iproj] = get_projection_matrix(MB_basis,iproj);};

//////////////////////////
write_basis(MB_basis); // write basis in a file
cout << "Coulomb\n";
#include"write_all.h"  // write all the operators in a file
SpMat ham(d,d); // hamiltonian
ham = -im*(lx*ly - ly*lx) ; // get LS term
ham = ls ; // get LS term
//ham = ly ; // get LS term
//ham = jz ; // get LS term
//ham = sx*sx+sy*sy+sz*sz ; // get LS term
//ham = j2 ; // get LS term
DMat dham(MB_basis.size,MB_basis.size); // dense hamiltonian
dham = DMat(ham);  // convert to dense matrix
//cout << dham << endl;  // print in screen
SelfAdjointEigenSolver<MatrixXcf> eigensolver(dham);
//cout << eigensolver.eigenvalues() << endl;
//SelfAdjointEigenSolver<complex<double> > eigensolver;
//cout << eigensolver.eigenvalues();
//cout << get_ls_matrix(MB_basis);
cout << "DONE!"<<endl;
return 0;
}
