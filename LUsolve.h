
// ---------------------------------------------
// Solve (A - labmda I) xhat = x
// ---------------------------------------------

// A = LU, so solve the above system in a 2-step process

void LUsolve(double **L ,
	     double **U ,
	     double *Sol ,
	     double *RHS , int m , int n , double **A)
{
double **temp = Array2D_double(m,n);
for(int i = 0; i < n; i++){
	for(int j =0; j <n; j++)
	{
		temp[i][j] = L[i][j];
	}
}
double **temp2 = Array2D_double(m,n);
for(int i = 0; i < n; i++){
	for(int j =0; j <n; j++)
	{
		temp2[i][j] = U[i][j];
	}
}
  double y   [n];
for(int p = 0; p<n; p++){y[p]=RHS[p];}
//iLOOP cout << RHS[i] << endl;
  // Solve L*y = x
	for ( int pivotRow = 0 ; pivotRow < m ; ++pivotRow ) 
  {
    int pivotCol = findPivot(m,n, temp, y, pivotRow,pivotRow);
   
    
    for ( int elimRow = pivotRow + 1 ; elimRow < m ; ++ elimRow )
      {
	double fac = -temp[elimRow][pivotCol]/temp[pivotRow][pivotCol];
	
	for ( int col = pivotCol ; col < n ; ++col ) temp[elimRow][col] += temp[pivotRow][col] * fac;
	
	y[elimRow] += y[pivotRow] * fac;
      }
  }
//printSystem("** post check **" , m , n , temp  , y );
  
  // Solve U*Sol = RHS
for ( int pivotRow = m-1 ; pivotRow >= 0 ; --pivotRow )
    {
		int pivotCol = findPivot(m,n, temp2, y, pivotRow,pivotRow);
      if ( ! isZero(temp2[pivotRow][pivotCol] ) )
      for ( int elimRow = pivotRow - 1 ; elimRow >= 0 ; -- elimRow )
	{
	  
	  double fac = -temp2[elimRow][pivotCol]/temp2[pivotRow][pivotCol];
		  
	  for ( int col = 0 ; col < n ; ++col ) temp2[elimRow][col] += temp2[pivotRow][col] * fac;
		  
	  y[elimRow] += y[pivotRow] * fac;
	}
    }
	
for ( int i = 0 ; i < m ; ++i )
    {
		int pivotCol = findPivot(m,n, U, y, i,i);
      double denom = temp2[i][pivotCol];
      if ( ! isZero(denom) )
      {
	y[i] /= denom;
	for ( int j = 0 ; j < n ; ++j ) temp2[i][j] /= denom;
      }
    }
for(int p = 0; p<n; p++){Sol[p]=y[p];}
// Check solution LU * Sol = RHS, or Asave*Sol = RHS

  if ( n == m ) // Must be true for there to be a unique solution
    
    for ( int i = 0 ; i < m ; ++i )
      {
	double error = RHS[i];
	for ( int j = 0 ; j < n ; ++j )  error -= A[i][j]*Sol[j];

	if ( fabs(error) > 1.e-10 )
	  {
	    cout << "LU Solve failed\n";
	    exit(0);
	  }
      }
}

