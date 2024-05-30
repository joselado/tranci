// common libraries
//#include </usr/include/python2.7/Python.h>
#include<vector>
#include <fstream>
#include <math.h>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>



using namespace std;
using namespace Eigen;
using Eigen::SparseMatrix;  // Complex sparse matrix
typedef SparseMatrix<complex<double> > SpMat;  // Complex sparse matrix
typedef MatrixXcf DMat;  // Complex sparse matrix
typedef vector< vector< vector< vector<complex<double> >  > > >  Tensor4;  // Complex sparse matrix
//typedef vector< vector<int> > FourTen;  // Complex sparse matrix
complex<double> im(0,1); //complex unit


// This files includes all the files needed by the ci calculation

// objects
#include "objects.h" // global objects


// special operators
#include "overload.h" // specialized operators




// functions
#include "generate_basis.h" // get full basis
#include "get_basis_nelectrons.h" // crop the basis
//#include "create.h"  // function that creates a particle
//#include "destroy.h"  // function that destroys a particle
#include "build_matrix.h" // build creation/anhilation operator



// one body to many body converter
#include "one2many.h"



// special level routines
#include "projection.h" // routines for angular momentum stuff
#include "angular.h" // routines for angular momentum stuff
#include "spin.h" // routines for spin operators
#include "spinorb.h" // routines for spin orbit coupling
#include "coulomb.h" // coulomb matrix elements


// writing functions
#include "write_matrix.h" // write matrix in a file
#include "write_basis.h" // write matrix in a file
