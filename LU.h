

// --------------------------
// Forward Phase
// --------------------------



for ( int pivotRow = 0 ; pivotRow < m ; ++pivotRow ) 
  {
    pivotCol = findPivot(m,n, A, b, pivotRow,pivotRow);
    for(int j = 0; j < m; j++){
      if(j < pivotRow) {L[j][pivotRow] = 0;}
    else{
      L[j][pivotRow]=A[j][pivotRow]/A[pivotCol][pivotRow];
      for(int k=0; k<pivotRow; k++){
        L[j][pivotRow] = A[j][pivotRow]/A[pivotRow][pivotCol];
      }
    }
  }
    // Populate L during elimination
    L[pivotRow][pivotCol] = 1.;
    
    
    // Now, eliminate
    
    for ( int elimRow = pivotRow + 1 ; elimRow < m ; ++ elimRow )
      {
	double fac = -A[elimRow][pivotCol]/A[pivotRow][pivotCol];
	
	for ( int col = pivotCol ; col < n ; ++col ) A[elimRow][col] += A[pivotRow][col] * fac;
	
	b[elimRow] += b[pivotRow] * fac;
      }
  }

// --------------------------
// Print Echelon Form
// --------------------------

printSystem("** Echelon Form **" , m , n , A  , b );
printSystem("** L **  "          , m , m , L  , b );


// --------------------------
// Check LU
// --------------------------

MatMat( L , m , A , n , m , LU );
printSystem("** LU Check **  "  , m , n , LU  , b );
