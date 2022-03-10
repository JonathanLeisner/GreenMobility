//------------------------------------------------------------------------------
// CUtils.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Utilities for use in std. error module.
//------------------------------------------------------------------------------

// Function to compute emax given coefficients from model solution
double Emax_Hat(double *PI_Coef, double *rsk, int exper, double lastwage)
{
	double Emax_Out;
    
	Emax_Out =   1                   * PI_Coef[0]
               + rsk[0]              * PI_Coef[1]
               + rsk[1]              * PI_Coef[2]
               + rsk[2]              * PI_Coef[3]
               + rsk[3]              * PI_Coef[4]
               + rsk[4]              * PI_Coef[5]
               + lastwage            * PI_Coef[6]
               + exper               * PI_Coef[7]
               + rsk[0] * rsk[0]     * PI_Coef[8]
               + rsk[1] * rsk[1]     * PI_Coef[9]
               + rsk[2] * rsk[2]     * PI_Coef[10]
               + rsk[3] * rsk[3]     * PI_Coef[11]
               + rsk[4] * rsk[4]     * PI_Coef[12]
               + lastwage * lastwage * PI_Coef[13]
               + exper  * exper      * PI_Coef[14]
               + rsk[0] * rsk[1]     * PI_Coef[15]
               + rsk[0] * rsk[2]     * PI_Coef[16]
               + rsk[0] * rsk[3]     * PI_Coef[17]
               + rsk[0] * rsk[4]     * PI_Coef[18]
               + rsk[0] * lastwage   * PI_Coef[19]
               + rsk[0] * exper      * PI_Coef[20]
               + rsk[1] * rsk[2]     * PI_Coef[21]
               + rsk[1] * rsk[3]     * PI_Coef[22]
               + rsk[1] * rsk[4]     * PI_Coef[23]
               + rsk[1] * lastwage   * PI_Coef[24]
               + rsk[1] * exper      * PI_Coef[25]
               + rsk[2] * rsk[3]     * PI_Coef[26]
               + rsk[2] * rsk[4]     * PI_Coef[27]
               + rsk[2] * lastwage   * PI_Coef[28]
               + rsk[2] * exper      * PI_Coef[29]
               + rsk[3] * rsk[4]     * PI_Coef[30]
               + rsk[3] * lastwage   * PI_Coef[31]
               + rsk[3] * exper      * PI_Coef[32]
               + rsk[4] * lastwage   * PI_Coef[33]
               + rsk[4] * exper      * PI_Coef[34]
               + lastwage * exper    * PI_Coef[35];
    
	return Emax_Out;
}


// Function to find max over array
double Find_Max(double *V, int dim)
{
	double max;
	int i;
	max = V[0];

	for (i = 1; i < dim; i++)
	{
		if (V[i] > max)
		{
			max = V[i];
		}
	}

	return max;
}

// Function to find min over array
double Find_Min(double *V, int dim)
{
	double min;
	int i;
	min = V[0];

	for (i = 1; i < dim; i++)
	{
		if (V[i] < min && V[i] > 0)
		{
			min = V[i];
		}
	}

	return min;
}

// Function to compare sign of two floats
int SameSign(double a, double b) {
	int y;
    if (a*b >=0)
    {
    	return y=1;
    }
    else
    {
    	return y=0;
    }
}


//--------------------------------------------------------------------------------------
// Function to compute the inverse of a positive definite square matrix X
//--------------------------------------------------------------------------------------
int fcn_invX(double *X, int n, double *invX)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------

	// Counters and info
	int i, j, info;

	// Upper triangular placeholder, U, and X'X
	double *U;

	// Indicator for upper triangular
	char s;


	//---------------------------------------------------------------
	// Memory Allocation
	//---------------------------------------------------------------
	U  = calloc(n*n,sizeof(double));


	//---------------------------------------------------------------
	// Compute inverse of X
	//---------------------------------------------------------------
	for (i = 0; i < n*n; i++)
    {
    	U[i] = X[i];
    }

    // Upper Cholesky factorization
	s = 'U'; // upper
	dpotrf_(&s, &n, U, &n, &info);
	// info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, s, n, U, n);

	for (i = 0; i < n*n; i++)
    {
    	invX[i] = U[i];
    }

    // Get inverse
    dpotri_(&s, &n, invX, &n, &info);
    // info = LAPACKE_dpotri(LAPACK_COL_MAJOR, s, n, invX, n);

    for (i = 0; i < n; i++)
    {
    	for (j = 0; j < n; j++)
    	{
    		if (i > j)
    		{
    			invX[i + j * n] = invX[j + i * n];
    		}
    	}
    }


	//---------------------------------------------------------------
	// Free Memory
	//---------------------------------------------------------------
	free(U);

	return info;

}



//--------------------------------------------------------------------------------------
// Function to compute the inverse of X'X
//--------------------------------------------------------------------------------------
int fcn_invXX(double *X, int rows, int cols, double *invXX)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------

	// Counters and info
	int i, j, info;

	// Upper triangular placeholder, U, and X'X
	double *UU, *XX;

	// Indicator for upper triangular
	char s;

	// Alpha and beta for matrix multiplication - see dgemm documentation
	double alpha, beta;


	//---------------------------------------------------------------
	// Memory Allocation
	//---------------------------------------------------------------
	UU = calloc(cols*cols,sizeof(double));
	XX = calloc(cols*cols,sizeof(double));


	//---------------------------------------------------------------
	// Compute inverse of X'X
	//---------------------------------------------------------------
	alpha = 1.0;
	beta  = 0.0;

	// Get X'X
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, cols, cols, rows, alpha, X, rows, X, rows, beta, XX, cols);

	for (i = 0; i < cols*cols; i++)
    {
    	UU[i] = XX[i];
    }

    // Upper Cholesky factorization
	s = 'U'; // upper
	dpotrf_(&s, &cols, UU, &cols, &info);
	// info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, s, cols, UU, cols);

	for (i = 0; i < cols*cols; i++)
    {
    	invXX[i] = UU[i];
    }

    // Get inverse
    dpotri_(&s, &cols, invXX, &cols, &info);
    // info = LAPACKE_dpotri(LAPACK_COL_MAJOR, s, cols, invXX, cols);

    for (i = 0; i < cols; i++)
    {
    	for (j = 0; j < cols; j++)
    	{
    		if (i > j)
    		{
    			invXX[i + j * cols] = invXX[j + i * cols];
    		}
    	}
    }


	//---------------------------------------------------------------
	// Free Memory
	//---------------------------------------------------------------
	free(UU);
	free(XX);

	return info;

}


//--------------------------------------------------------------------------------------
// Function to compute X'Y
//--------------------------------------------------------------------------------------
void fcn_XY(double *X, double *Y, int rows, int colsX, double *XY)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------

	// Cols of Y
	int colsY;

	// Alpha and beta for matrix multiplication - see dgemm documentation
	double alpha, beta;

	//---------------------------------------------------------------
	// Compute inverse of X'Y
	//---------------------------------------------------------------
	alpha = 1.0;
	beta  = 0.0;
	colsY = 1;

	// Get X'X
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, colsX, colsY, rows, alpha, X, rows, Y, rows, beta, XY, colsX);

}


//--------------------------------------------------------------------------------------
// Function to compute A*B
//--------------------------------------------------------------------------------------
void fcn_matmul(double *A, double *B, int m, int n, int k, double *AB)
{
	// A[m,n],		B[n,k],		AB[m,k]

	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------

	// Alpha and beta for matrix multiplication - see dgemm documentation
	double alpha, beta;

	//---------------------------------------------------------------
	// Compute A*B
	//---------------------------------------------------------------
	alpha = 1.0;
	beta  = 0.0;

	// Get X'X
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, A, m, B, n, beta, AB, m);

}

//--------------------------------------------------------------------------------------
// Function to do linear regression
//--------------------------------------------------------------------------------------
int fcn_linreg(double *Y, double *X, int m, int n, double *coef, double *SST, double *SSE)
{
	//--------------------------------------------------------------
	// X[m,n],		Y[m,1],		coef[n,1]
	//--------------------------------------------------------------

	//--------------------------------------------------------------
	// Declarations
	//--------------------------------------------------------------
	int i, info;
	int one = 1;
	double *invXX, *XY, *Y_hat;
	double Y_bar;

	*SSE  = 0.0;
	*SST  = 0.0;
	Y_bar = 0.0;


	//--------------------------------------------------------------
	// Allocate memory
	//--------------------------------------------------------------
	invXX = calloc(n*n,sizeof(double));
	XY    = calloc(m,sizeof(double));
	Y_hat = calloc(m,sizeof(double));


	//--------------------------------------------------------------
	// Do OLS
	//--------------------------------------------------------------
	info = fcn_invXX(X, m, n, invXX);
    fcn_XY(X, Y, m, n, XY);
    fcn_matmul(invXX, XY, n, n, one, coef);


    //--------------------------------------------------------------
    // Basic statistics
    //--------------------------------------------------------------
    fcn_matmul(X, coef, m, n, one, Y_hat);

    for (i = 0; i < m; i++)
    {
    	*SSE  += (Y[i] - Y_hat[i]) * (Y[i] - Y_hat[i]);
    	Y_bar += Y[i] / m;
    }

    for (i = 0; i < m; i++)
    {
    	*SST += (Y[i] - Y_bar) * (Y[i] - Y_bar);
    }


    //--------------------------------------------------------------
    // Free allocated memory
    //--------------------------------------------------------------
    free(invXX);
    free(XY);
    free(Y_hat);

    return info;
	
}

//--------------------------------------------------------------------------------------
// Function to unpack parameters
//--------------------------------------------------------------------------------------
int unpackPars(double *param, double *beta, double *sigma, double *gamma, double *xi, double *kappa, double *nu, double *lambda, double *phi)
{
    int iParam, iSec, lagsec;

    for (iParam = 0; iParam < nParamBeta*nSector; iParam++)
    {
        beta[iParam] = param[iParam];
    }

    for (iSec = 0; iSec < nSector; iSec++)
    {
        sigma[iSec] = param[iSec + nParamBeta*nSector];
    }

    for (iSec = 0; iSec < nSector + 1; iSec++)
    {
        gamma[iSec] = param[iSec + (nParamBeta+1)*nSector];
    }

    xi[0] = 0.0;
    for (lagsec = 1; lagsec < nSector; lagsec++)
    {
    	xi[lagsec + 0 * nSector] = param[lagsec + (nParamBeta+2)*nSector];
    }

    for (iSec = 0; iSec < nSector; iSec++)
    {
        xi[iSec + 1 * nSector] = param[iSec + (nParamBeta+3)*nSector];
    }

    for (iParam = (nParamBeta+4)*nSector; iParam < (nParamBeta+5)*nSector - 1; iParam++)
    {
        kappa[iParam - (nParamBeta+4)*nSector] = param[iParam];
    }

    nu[0] = param[(nParamBeta+5)*nSector - 1];

    for (iParam = (nParamBeta+5)*nSector; iParam < (nParamBeta+6)*nSector; iParam++)
    {
    	lambda[iParam - (nParamBeta+5)*nSector] = param[iParam];
    }

    for (iParam = (nParamBeta+6)*nSector; iParam < (nParamBeta+7)*nSector - 1; iParam++)
    {
    	phi[iParam - (nParamBeta+6)*nSector] = param[iParam];
    }

    return 0;
}
