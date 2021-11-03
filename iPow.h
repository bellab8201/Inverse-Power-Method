#include <iostream>
#include <string>
#include "stdio.h"
#include "math.h"
#include "string.h"

using std :: cout;
using std :: endl;
using std :: string;

#define iLOOP   for ( int i = 0 ; i < m ; ++i )
#define jLOOP   for ( int j = 0 ; j < n ; ++j )


// =========================================================================
// Allocate 2D Array
// =========================================================================

double  * * Array2D_double(int nRows,int nCols)
{
  double *myArray;
  myArray = new double [ nRows * nCols ];
    
  // Create a pointer that points to the beginning of each new row

  double * * myArray_ptr;
  myArray_ptr = new double * [nRows];

  int count = 0;

  for ( int row = 0 ; row < nRows ; ++ row )
    {
      myArray_ptr[row] = &myArray[ count*nCols ];
      ++count;
    }

  // Return that pointer
  
  return myArray_ptr;

}


// =========================================================================
// Small Number Filters
// =========================================================================

double filter(double val)
{
  double smallVal = 1.e-10;

  if ( fabs(val) < smallVal ) return 0.;
  return val;
}


bool isZero(double val)
{
  double smallVal = 1.e-10;
  if ( fabs(val) < smallVal ) return true;
  return false;
  
}

bool isOne(double val)
{
  double smallVal = 1.e-10;
  if ( fabs(val-1.0) < smallVal ) return true;
  return false;
  
}

// =========================================================================
// Compute the ratio of numer/denom
// =========================================================================

string ratio( double numer , double denom )
{
  if ( ! isZero(denom) )
    {
      char str[20];
      sprintf(str, "%f", numer/denom);
      return str;
    }

  if ( isZero(numer) ) return "----"; 

  return "inf";
}


// =========================================================================
// printSystem and results
// =========================================================================

void printSystem(string message, int m , int n , double **A , double *b)
{
  
  cout << endl;
  cout << "--------------------------------------" << endl;
  cout << message << endl;
  cout << "--------------------------------------" << endl;
  cout << endl;

  for ( int i = 0 ; i < m ; ++i )
    {
      for ( int j = 0 ; j < n ; ++j )
	{
	  double out = A[i][j];
	  printf("%5.2f ", filter(out));
	}
      printf("      |    %5.2f ", filter(b[i]) );
      cout << endl;
    }
    
}

// =========================================================================
// findPivot
// =========================================================================

int findPivot( int m , int n , double **A , double *b , int row , int col )
{

  for ( int colTry = col ; colTry < n ; ++colTry )
    {
      // Find largest pivot in this column, not counting rows above "row"

      int largest = row;
  
      for ( int i = row ; i < m ; ++i )
	  if ( fabs( A[i][colTry] ) > fabs( A[ largest][ colTry ] ) ) largest = i;

      if ( !isZero( A[row][ colTry] ) ) return colTry;

    }

  return -1;

}


// =========================================================================
// mat-mat multiply
//
//    A ( mA x nAB ) * B ( nAB x nB ) = C( mA x nB )
//
// =========================================================================

void MatMat( double **A , int mA , double **B, int nB , int nAB ,  double **C )
{
  for ( int i = 0 ; i < mA ; ++i )
    for ( int j = 0 ; j < nB ; ++j )
      {
	C[i][j ] = 0.;

	for ( int k = 0 ; k < nAB; ++k ) C[ i ] [ j ] +=  A[i][k] * B [k][j];

      }
}


