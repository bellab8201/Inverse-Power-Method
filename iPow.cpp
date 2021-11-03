#include "iPow.h"
#include "LUsolve.h"




// =========================================================================
// M A I N   P R O G R A M
// =========================================================================

int main(int argc, char *argv[])
{
  int m; // number of rows
  int n; // number of columns

  #include "setCase.h"

  // ---------------------------------------------
  // Set up Inverse Power Method
  // ---------------------------------------------

  int    iter      = 0;
  int    maxIter   = 100;
  int    converged = 0;
  double tol       = 0.00001;
  double mu = 0;
  double mu_old;
  double lambda;

  double x[n];     // Used during Inverse Power Method Iteration, see class notes
  double xhat[n];  //   "   "       "       "      "      "        "    "     "
  double temp[n];
  // Compute initial guess at x
  
  for ( int i = 0; i < n; i++ ){
    x[i] = alpha;
    xhat[i] =alpha;
  } 

  // Compute (A - alpha*I) and store in A.

  for (int i = 0; i<n; i++) A[i][i] = A[i][i] - alpha;

  // Save new A = (A - alpha*I) in Asave:

  for (int i = 0; i<n; i++) 
  {
    for(int j=0; j < n; j++)
    {
      Asave[i][j] = A[i][j];
    }
  }
  // Print (A - alpha*I)
  printSystem("A-alpha*I",m,n,Asave,&b[0]  );
  
  // ---------------------------------------------
  // Perform LU factorization of new A
  // ---------------------------------------------
  #include "LU.h"

  // A now stores what.?  Point to the upper triangular A with  U:

    double **U = A;

  // ---------------------------------------------
  // Inverse Power Method Iterations
  // ---------------------------------------------
  while (  converged != 1 )
    { 
      //Max iteration check
      
      if (iter >= maxIter){
        
        break;
        }

      // Step 2a: Compute xhat (Solve  (A - alpha I) * xhat = x)
      LUsolve(L , U , xhat , x , m , n , Asave );   // Asave is passed to check solution
      // Step 2b: Scale xhat so that the largest value = 1
        mu_old = mu;
        mu = 0;
        
      for(int i = 0; i <n; i++){
        if(fabs(xhat[i])>mu){
          mu = xhat[i];
        }
      }
      for(int j = 0; j < n; j++){
        xhat[j] /= mu;
      }
      for(int j = 0; j < n; j++){
        x[j] = xhat[j];
      }
      // Step 2c: Check for convergence
      
      if(fabs(mu-mu_old) < tol){
        converged = 1;
        }
        iter++;
      cout << "iter = " << iter << " mu = " << mu << endl;
      
    }

  // Compute final estimate for eigenvalue of original A, Aoriginal:

    lambda = alpha + 1/mu;
  
    cout << caseName << ": Inverse Power Method Converged in " << iter << " iterations." << endl;
    cout << caseName << ": =========================== Results " << endl;
    cout << caseName << ": mu = " << mu << endl;
    cout << caseName << ": lambda = " << lambda << endl;
    cout << caseName << ": ===========================         " << endl;
  
  
  // ---------------------------------------------
  // Test Eigenvalue/Eigenvector
  // ---------------------------------------------

    // Print eigenvector and Aoriginal * eigenvector, computing the ratio
    // to see how consistent that ratio is.

    double Atest[n];

    for ( int i = 0 ; i < m ; ++i )
    //for ( int j = 0 ; j < 1 ; ++j )
      {
	  Atest[i] = 0.;
    
	  for ( int k = 0 ; k < n; ++k ) {
        Atest[ i ] +=  Aoriginal[i][k] * xhat [k];
        //cout << "a " << Aoriginal[i][k] << " xhat " << xhat[k] << endl;
    }
      }
      
    //MatMat(Aoriginal,m,xhat2,n,n,Atest);
    
    iLOOP printf("%s: E-vector Check:     x[%d] = %5.2f     Ax[%d] = %5.2f     ratio = %s\n",
		caseName.c_str(),i,xhat[i] ,i,Atest[i],ratio(Atest[i], xhat[i]).c_str()   );

    
    
  
}

