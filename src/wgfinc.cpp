#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include "mex.h"
#include "FLT.h"
#include "PoissonSolver.h"
#include "Helper_U.h"
#include "BFM.h"
#include "BFMSlow.h"
#include "BFMIncomp.h"
using namespace std;


void create_dat_file(const double* A, int size, const string& filename, const int n, const int saveData){

    if(saveData == 0) return;

    char filename_tmp[100];
    sprintf(filename_tmp, "%s/rho-%04d.dat", filename.c_str(), n);

    // save as bin files
    ofstream out(filename_tmp, ios::out | ios::binary);
    if(!out) {
        cout << "Cannot open file.";
        return;
    }

    out.write((char *) A, size*sizeof(double));
    out.close();
}

void check_saveData_input(char* input_char, const mxArray *prhs[], int charSize, int idx){

    mxGetString(prhs[idx], input_char, charSize);
    
    const size_t len = strlen(input_char);
    if( input_char[len-1] == '\\' || input_char[len-1] == '/' )
    {
        // Terminate the string earlier
        input_char[len-1] = 0;
    }
}


void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    double *initmu   = mxGetPr(prhs[0]);
    double *V        = mxGetPr(prhs[1]);
    double *obstacle = mxGetPr(prhs[2]);
    int max_iteration= (int) mxGetScalar(prhs[3]);
    double tolerance = (double) mxGetScalar(prhs[4]);
    int nt           = (int) mxGetScalar(prhs[5]);
    double tau       = (double) mxGetScalar(prhs[6]);

    int    charSize   = 100;
    char*  input_char = new char[charSize];

    double m = 0;
    double gamma      = 0;
    int    verbose    = 0;

    int saveData = 1;
    check_saveData_input(input_char, prhs, charSize, 7);

    string filename = input_char;

    verbose = (int) mxGetScalar(prhs[8]);

    int n1=mxGetM(prhs[0]);
    int n2=mxGetN(prhs[0]);
    
    int pcount=n1*n2;
    
    plhs[0] = mxCreateDoubleMatrix(n1,n2,mxREAL);
    
    double *result = mxGetPr(plhs[0]); // result mu after nt outer iterations

    double sum = 0;

    double* mu = new double[n1*n2];
    memcpy(mu, initmu, n1*n2*sizeof(double));

    for(int i=0;i<pcount;i++){
        if (mu[i]<0) mexErrMsgTxt("Initial density contains negative values");
        sum += mu[i];
    }

    if(verbose > 0){
        cout << "XXX Incompressible Flow XXX\n\n";
        
        cout << "n    : " << n1 <<"\n";
        cout << "nt   : " << nt <<"\n";
        cout << "tau  : " << tau << "\n";

        cout << "Max Iteration : " << max_iteration <<"\n";
        cout << "Tolerance     : " << tolerance <<"\n";
    }

    BackAndForth* bf  = nullptr;
    Helper_U* helper_f= nullptr;

    // Initialize the Back-and-Forth method
    bf       = new BackAndForthIncompressible(n1,n2,max_iteration,tolerance,tau);
    // Initialize the internal energy functional U
    helper_f = new Helper_U_incompressible(n1,n2,tau,mu);   

    // Set a potential V
    helper_f->set_V(V);
    // Set an obstacle
    helper_f->set_obstacle(obstacle);

    create_dat_file(mu,n1*n2,filename,0,saveData);

    clock_t time;
    time=clock();

    for(int n=0;n<nt;++n){
        bf->start_OT(helper_f, mu, n, verbose);
        helper_f->calculate_DUstar_normalized(bf->phi_);
        memcpy(mu,helper_f->DUstar_,n1*n2*sizeof(double));
        create_dat_file(mu,n1*n2,filename,n+1,saveData);
    }

    memcpy(result, mu, n1*n2*sizeof(double));

    time=clock()-time;
    if(verbose>0) printf ("\nCPU time for GF: %.1f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    delete[] input_char;
}