//#include<stdlib.h>

//using namespace std;


// vector object, contains a single vector
class Basis_vector {  // create basis object
public:
  bool * basevec; // variable to store one wavefunction
  int size; // dimension of the basis
  bool is_zero; // if it is the zero vector
  int sign; // sign of the basis vector
  Basis_vector() { // constructor
    size=1; // default 1 size
    is_zero = false; // assume non zero state
    basevec = new bool [size]; // allocate size
    sign = 1; // sign of the vector
  }
  void alloc (int n) { // allocate 
    basevec = new bool [n]; // allocate size
    size = n;
    };
  void print() { // prints the vector in screen
    int i; // counter
    if (is_zero) cout << "ZERO";
    else for (i=0;i<size;i++) cout << basevec[i] << "";
    cout << "\n";
  };
  void set(int i,bool occ) { 
    basevec[i] = occ; };
  Basis_vector copy() { // copy method
    Basis_vector vo; // output vector
    int i; // counter
    vo.alloc(size); // set size
    vo.is_zero = is_zero; // set zero
    for (i=0;i<size;i++) vo.basevec[i] = basevec[i]; // copy components
    return vo; // return the copy
  };
}; // end basis object


int get_number(Basis_vector bvec) {  // get the number of this vector
  int ind; // index to return
  int acc; // accumulator
  int i; // counters
  ind = 0;
  acc = 1;
  if (bvec.is_zero) return -1; // if the vector is null
  for (i=0;i<bvec.size;i++)  {if (bvec.basevec[i]) ind += acc; acc*=2;};
  return ind;
}


//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////


// basis objects, contains all the vectors
class Basis  {  
public:
  Basis_vector * basis; // basis of the system
  int size; // size of the basis
  int * indexes; // indexes of the wavefunctions
  bool ** basevec; // different vectors
//  map<int, int> wf2ind;  // map for getting the index
  Basis() { // constructor
    size=1;
    basis = new Basis_vector[size]; } // allocate size
  void alloc (int n) { // allocate the basis
    basis = new Basis_vector[n];  // allocate size
    basevec = new bool*[n]; // allocate the different vectors
    if (basis==NULL) cout << "SHIT" ;
    size = n; // size of the basis
  };
  void print () {
    int i; // counter
    cout << size << "  basis vectors\n";
    for (i=0;i<size;i++) {
      cout <<get_number(basis[i]) <<  "      ";
      basis[i].print() ;
    };
  };
  void get_indexes()  {// creates the map of the wavevectors
    int td; // total dimension
    td = (int)pow(2,basis[0].size); // total dimension
    indexes = new int[td]; // allocate huge array
    int i; // counter
    for (i=0;i<td;i++) indexes[i] = -1; // default value
    for (i=0;i<size;i++) indexes[get_number(basis[i])] = i; // get number
  };
} ;


// class for creation/destruction operators
class CreaDes_operator {
  public:
  int level; // level in which creates
};

// creation operator
class Creation_operator: public CreaDes_operator {
  public:
  int level; // level in which creates
};



// destruction operator
class Destruction_operator: public CreaDes_operator {
  public:
  int level; // level in which destroys
};


