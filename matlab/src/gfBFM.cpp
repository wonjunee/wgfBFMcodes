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

void create_dat_file(const double* A, int size, string filename, const int n, const int saveData){

    if(saveData == 0) return;

    string filename_tmp="./data/" + filename + "-" +to_string(n) + ".dat";

    cout << "filename_tmp : " << filename_tmp << "\n";

    if(saveData == 1){
        // save as csv files
        ofstream out(filename_tmp);
        if(!out) {
            cout << "Cannot open file.";
            return;
        }
        for(int i=0;i<size;++i){
            out << A[i];
            if(i<size-1) out << ",";
        }
        out.close();    
    }else if(saveData == 2){
        // save as bin files
        ofstream out(filename_tmp, ios::out | ios::binary);
        if(!out) {
            cout << "Cannot open file.";
            return;
        }

        out.write((char *) A, size*sizeof(double));
        out.close();
    }
    
}

bool compare_char_string(const char c[], const string& s){
    for(int i=0;i<s.length();++i){
        if(c[i] != s[i]) return false;
    }
    return true;
}

void check_saveData_input(int& saveData, string& filename, char* input_char, const mxArray *prhs[], int charSize, int idx){
    if(compare_char_string(input_char,"csv")){
        saveData = 1;
    }else if(compare_char_string(input_char,"bin")){
        saveData = 2;
    }else if(compare_char_string(input_char,"0")){
        saveData = 0;
    }else{
        mexErrMsgTxt("wrong input for saveData.\nsaveData should be one of the following: 'csv', 'bin', '0'");
    }

    if(saveData != 0) mxGetString(prhs[idx], input_char, charSize);
    
    filename = input_char;
    // filename = "./data/" + filename;
}

int check_flow_type(char* input_char, const mxArray *prhs[], int charSize, int idx){
    mxGetString(prhs[idx], input_char, charSize);

    if(compare_char_string(input_char,"slow")){
        return 0;
    }else if(compare_char_string(input_char,"incompressible")){
        return 1;
    }else{
        mexErrMsgTxt("wrong input for flowType.\nflowType should be 'slow' or 'incompressible'");
    }
    
    return 0;
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
    string filename   = "";

    int flowType = check_flow_type(input_char, prhs, charSize, 7);

    double m = 0;
    double gamma      = 0;
    int    verbose    = 0;
    int    saveData   = 0;

    if(flowType == 0){ // slow diffusion
        m = (double) mxGetScalar(prhs[8]);

        assert( m > 1 );

        gamma   = (double) mxGetScalar(prhs[9]);
        verbose = (int)    mxGetScalar(prhs[10]);
        mxGetString(prhs[11], input_char, charSize);
        check_saveData_input(saveData, filename, input_char, prhs, charSize, 12);
    }else{ // incompressible
        verbose = (int)    mxGetScalar(prhs[8]);
        mxGetString(prhs[9], input_char, charSize);
        check_saveData_input(saveData, filename, input_char, prhs, charSize, 10); 
    }

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
        if(flowType == 0){
            cout << "XXX Slow Diffusion XXX\n\n";    
        }else{
            cout << "XXX Incompressible Flow XXX\n\n";
        }
        
        cout << "n    : " << n1 <<"\n";
        cout << "nt   : " << nt <<"\n";
        cout << "tau  : " << tau << "\n";

        if(flowType == 0){
            cout << "m    : " << m << "\n";
            cout << "gamma: " << gamma << "\n";    
        }

        cout << "Max Iteration : " << max_iteration <<"\n";
        cout << "Tolerance     : " << tolerance <<"\n";
    }

    BackAndForth* bf  = nullptr;
    Helper_U* helper_f= nullptr;

    if(flowType == 0){
        // Initialize the Back-and-Forth method
        bf       = new BackAndForthSlow(n1,n2,max_iteration,tolerance,gamma,tau,m);
        // Initialize the internal energy functional U
        helper_f = new Helper_U_slow(n1,n2,gamma,tau,m,mu);     
    }else{
        // Initialize the Back-and-Forth method
        bf       = new BackAndForthIncompressible(n1,n2,max_iteration,tolerance,tau);
        // Initialize the internal energy functional U
        helper_f = new Helper_U_incompressible(n1,n2,tau,mu);   
    }
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
    if(verbose>0) printf ("\nCPU time for GF: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    delete[] input_char;
}