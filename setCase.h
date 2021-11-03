

// -----------------------------------
// Default command-line option values
// -----------------------------------

string caseName = "";
double alpha = 0.;

// ----------------------------
// Process command-line options
// ----------------------------

for (int count =  0 ; count < argc; ++count)
  {
    if ( !strcmp(argv[count],"-case"  ) ) caseName  =argv[count+1];
    if ( !strcmp(argv[count],"-alpha" ) ) alpha     =atof(argv[count+1]);
  }

if ( caseName == ""      ) { cout << "Choose a code with -case <ID>\n"; exit(0); }
if ( fabs(alpha) < 1.e-10) { cout << "Set initial guess, alpha using -alpha <value>\n"; exit(0); }

// --------------------------
// Initialize the array
// --------------------------

// p287

int m1 = 4;
int n1 = 4;
double A1[m1][n1] = {  {  5 ,   0  ,  0 ,  0 } ,
		       {  0 ,   5  ,  0 ,  0 } ,
		       {  1 ,   4  , -3 ,  0 } ,
		       { -1 ,  -2  ,  0 , -3 } } ;
double b1[m1]    =     {  1.  ,  3. ,  4. ,  2.     };


// p285

int   m2 = 3;
int   n2 = 3;
double A2[m2][n2] = { {  1 , 3  ,  3  } ,
		      { -3 ,-5  , -3  } ,
		      {  3 , 3  ,  1  }   };
double b2[m2]    = { 1 , 4 , 4 };

// p270

int   m3 = 3;
int   n3 = 3;
double A3[m3][n3] = { {  4 , -1 , 6   },
		      {  2 ,  1 , 6   },
		      {  2 , -1 , 8   } };
  double b3[m3]    =  { 1  ,  4 ,  4  };

// --------------------------
// Set the case
// --------------------------

m = 0;
n = 0;

if ( caseName == "p287"      )  {      m = m1;      n = n1;    }
if ( caseName == "p285"      )  {      m = m2;      n = n2;    }
if ( caseName == "p270"      )  {      m = m3;      n = n3;    }

if ( m == 0 ) { cout << "Bad case selection\n"; exit(0);}

double b[m];
double **A = Array2D_double(m,n);  // Note new method of acquiring 2D array, see iPow.h

if ( caseName == "p287"      )   iLOOP  { b[i] = b1[i];  jLOOP A[i][j] = A1[i][j];   }
if ( caseName == "p285"      )   iLOOP  { b[i] = b2[i];  jLOOP A[i][j] = A2[i][j];   }  
if ( caseName == "p270"      )   iLOOP  { b[i] = b3[i];  jLOOP A[i][j] = A3[i][j];   }

double **L  = Array2D_double(m,n);   // Note new method of acquiring 2D array, see iPow.h
double **LU = Array2D_double(m,n);   //   "   "    "    "      "     "    "     "    "

for ( int i = 0 ; i < m ; ++i ) for ( int j = 0 ; j < m ; ++j ) L[i][j] = 0.;

cout << "\n";
cout << "\n";
cout << "\n";
cout << "========================================================\n";
cout << "========================================================\n";
cout << "Running case " + caseName << " with alpha = " << alpha << "\n";
cout << "========================================================\n";
cout << "========================================================\n";
cout << "\n";
cout << "\n";
cout << "\n";


  // --------------------------
  // Save Original Matrix
  // --------------------------

double ** Asave              = Array2D_double(m,n);
double ** Aoriginal          = Array2D_double(m,n);

iLOOP jLOOP Aoriginal[i][j]  = A[i][j];

double bsave[n]   ;
int    pivotCol;
  
printSystem("Original Matrix" , m , n , A , &b[0] );


