/**
 * 
 * Poisson solver using fftw.
 *
 */

#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include <iostream>
#include <fftw3.h>

using namespace std;

class poisson_solver{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *workspace; // This will contain the coefficients for complex numbers
    double *kernel;

    int n1;
    int n2;
    double dx;
    double dy;

    double eps;

    double coeff;

    poisson_solver(){
        workspace=NULL;
        kernel=NULL;
    }

    poisson_solver(int n1, int n2) {
        // initialize(n1,n2);
        this->n1=n1;
        this->n2=n2;

        // eps = sqrt(200.0/n1);
        eps = 0.05;

        workspace =(double*) fftw_malloc(n1*n2*sizeof(double));

        planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
        planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);
        
        kernel=new double[n1*n2];
        create_negative_laplacian_kernel_2d();
    }


    ~poisson_solver(){
        fftw_cleanup();
        delete[] kernel;
        fftw_free(workspace);
        fftw_destroy_plan(planIn);
        fftw_destroy_plan(planOut);
    }

    void create_negative_laplacian_kernel_2d(){
        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
                double negativeLaplacian=2.0*n1*n1*(1-cos(M_PI*(j)/n1)) + 2*n2*n2*(1-cos(M_PI*(i)/n2));
                kernel[i*n1+j]=negativeLaplacian;
            }
        }
    }

    void get_fourier_coefficients(const double* push_mu, const double* DFstar){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                workspace[i*n1+j]=-push_mu[i*n1+j]+DFstar[i*n1+j];
            }
        }

        fftw_execute(planIn);
    }

    void perform_inverse_laplacian(const double* push_mu, const double* DFstar, const double c1, const double c2, const double sigma){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                workspace[i*n1+j]=-push_mu[i*n1+j]+DFstar[i*n1+j];
            }
        }

        fftw_execute(planIn);
        // workspace[0]=0;

        for(int i=0;i<n1*n2;++i){
            workspace[i]/=4*(n1)*(n2)*(c1+c2*kernel[i])/sigma;
        }

        fftw_execute(planOut);
    }

    void perform_convolution(const double* u){
        memcpy(workspace, u, n1*n2*sizeof(double));

        fftw_execute(planIn);

        for(int i=0;i<n1*n2;++i){
            workspace[i] *= exp(-kernel[i]*eps*eps) / (4*n1*n2);
        }

        fftw_execute(planOut);
    }

    void solve_heat_equation(double* rho, double tau){
        memcpy(workspace,rho,n1*n2*sizeof(double));

        fftw_execute(planIn);
        
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                workspace[i*n1+j]*=exp(-2*tau*kernel[i*n1+j]);
            }
        }

        for(int i=0;i<n1*n2;++i){
            workspace[i]/=4.0*(n1)*(n2);
        }

        fftw_execute(planOut);

        memcpy(rho,workspace,n1*n2*sizeof(double));
    }


}; // Poisson Solver

#endif