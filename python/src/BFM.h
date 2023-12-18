#ifndef BFM_H
#define BFM_H

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <vector>

namespace py = pybind11;
using namespace std;

class BFM{
public:
    int n1;
    int n2;

    double beta_1;
    double beta_2;
    double alpha_1;
    double alpha_2;

    double tau_;


    // Initialize phi and psi
    double* gradx;
    double* grady;

    double* vxx_;
    double* vyy_;
    double* vxy_;

    Pushforward_mapping pushforward;

    BFM(int n1, int n2, double tau): n1(n1), n2(n2), tau_(tau){

        gradx=new double[n1*n2];
        grady=new double[n1*n2];

        vxx_=new double[n1*n2];
        vxy_=new double[n1*n2];
        vyy_=new double[n1*n2];

        beta_1=0.05;
        beta_2=0.95;
        alpha_1=1.2;
        alpha_2=0.8;

        // pushforward.initialize(n1,n2);
        // flt2d.initialize(n1,n2);

        pushforward = Pushforward_mapping(n1,n2);
    }

    ~BFM(){
        delete[] gradx;
        delete[] grady;
    }


    // This function will provide S1(x,) where x and y are n [0,1] double values
    double interpolate_function(double x,double y,const double* func){
        double indj=fmin(n1-1,fmax(0,x*n1-0.5));
        double indi=fmin(n2-1,fmax(0,y*n2-0.5));

        double lambda1=indj-(int)indj;
        double lambda2=indi-(int)indi;

        double x00 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj))];
        double x01 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj+1))];
        double x10 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj))];
        double x11 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj+1))];

        double interpolated_value = (1-lambda1)*(1-lambda2)*x00+(lambda1)*(1-lambda2)*x01
                                   +(1-lambda1)*(lambda2)  *x10+(lambda1)*(lambda2)  *x11;
        return interpolated_value;  
    }




    void calculate_push_pull_rho(double* push_rho, const double* rho, const double* vxx,const double* vyy,const double* vxy,const double* phi,const double det_threshold=0.9){

        double eps = pow(1.0/n1, 0.5);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
 
                // centered difference
                double vxval=0.5*n1*(phi[i*n1+(int)fmin(n1-1,j+1)]-phi[i*n1+(int)fmax(0,j-1)]);
                double vyval=0.5*n2*(phi[(int)fmin(n2-1,i+1)*n1+j]-phi[(int)fmax(0,i-1)*n1+j]);                

                double x=(j+0.5)/(1.0*n1)-tau_*vxval;
                double y=(i+0.5)/(1.0*n2)-tau_*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){

                    double vxx_val = interpolate_function(x,y,vxx);
                    double vyy_val = interpolate_function(x,y,vyy);
                    double vxy_val = interpolate_function(x,y,vxy);

                    double det = fabs((1.0-tau_*vxx_val) * (1.0-tau_*vyy_val) - tau_*tau_ * vxy_val * vxy_val);

                    if(det > det_threshold){
                        push_rho[i*n1+j] = rhovalue/det;    
                    }else{
                        int jpp = fmin(n1-1,j+2);
                        int jp  = fmin(n1-1,j+1);
                        int jm  = fmax(0,j-1);
                        int jmm = fmax(0,j-2);

                        int ipp = fmin(n2-1,i+2);
                        int ip  = fmin(n2-1,i+1);
                        int im  = fmax(0,i-1);
                        int imm = fmax(0,i-2);

                        vxx_val = 0.25*n1*n1* (phi[i*n1+jpp] - 2.*phi[i*n1+j] + phi[i*n1+jmm]);
                        vyy_val = 0.25*n2*n2* (phi[ipp*n1+j] - 2.*phi[i*n1+j] + phi[imm*n1+j]);
                        vxy_val = 0.25*n1*n2* (phi[ip*n1+jp] - phi[ip*n1+jm] - phi[im*n1+jp] + phi[im*n1+jm]);

                        det = fabs((1.0-tau_*vxx_val) * (1.0-tau_*vyy_val) - tau_*tau_ * vxy_val * vxy_val);
                        det = fmin(1.0/eps, det);
                        push_rho[i*n1+j] = rhovalue*det;    
                    }
                }else{
                    push_rho[i*n1+j]=0;
                }

            }
        }
    }


    void calculate_gradient_vxx_vyy_vxy(double* vxx, double* vyy, double* vxy, const double* phi){
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                int jpp = fmin(n1-1,j+2);
                int jp  = fmin(n1-1,j+1);
                int jm  = fmax(0,j-1);
                int jmm = fmax(0,j-2);

                int ipp = fmin(n2-1,i+2);
                int ip  = fmin(n2-1,i+1);
                int im  = fmax(0,i-1);
                int imm  = fmax(0,i-2);

                vxx[i*n1+j] = 0.25*n1*n1* (phi[i*n1+jpp] - 2.*phi[i*n1+j] + phi[i*n1+jmm]);
                vyy[i*n1+j] = 0.25*n2*n2* (phi[ipp*n1+j] - 2.*phi[i*n1+j] + phi[imm*n1+j]);
                vxy[i*n1+j] = 0.25*n1*n2* (phi[ip*n1+jp] - phi[ip*n1+jm] - phi[im*n1+jp] + phi[im*n1+jm]);
            }
        }
    }

    void calculate_gradient(const double* phi_c, double* gradx, double* grady){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1-1;++j){
                gradx[i*n1+j]=1.0*n1*(phi_c[i*n1+j+1]-phi_c[i*n1+j]);
            }
        }
        for(int j=0;j<n1;++j){
            for(int i=0;i<n2-1;++i){
                grady[i*n1+j]=1.0*n2*(phi_c[(i+1)*n1+j]-phi_c[i*n1+j]);
            }
        }
    }

    double calculate_h_minus_1(poisson_solver& fftps, const double* push_mu, const double* DFstar){
        double error=0;
        for(int i=0;i<n1*n2;++i){
            double value=-push_mu[i]+DFstar[i];
            error+=value*fftps.workspace[i];
        }
        return error/(1.0*n1*n2);
    }

    double initialize_sigma(const double* mu){
        double sigma=1.0;
        for(int i=0;i<n1*n2;++i){
            sigma=fmax(sigma,mu[i]);
        }
        return tau_/sigma;
    }

    double update_sigma(double sigma, const double W2_value, const double W2_value_previous, const double error)
    {
        // Update sigma
        if(W2_value_previous-W2_value>-sigma*beta_1*error){
            sigma*=alpha_2;
        }else if(W2_value_previous-W2_value<-sigma*beta_2*error){
            sigma*=alpha_1;
        }

        return sigma;
    }

    void compute_pull_back(
        py::array_t<double, py::array::c_style | py::array::forcecast> push_mu_np,
        py::array_t<double, py::array::c_style | py::array::forcecast> phi_np,
        py::array_t<double, py::array::c_style | py::array::forcecast> psi_np,
        py::array_t<double, py::array::c_style | py::array::forcecast> DUstar_np
    ){
        py::buffer_info push_mu_buf = push_mu_np.request();
        py::buffer_info phi_buf     = phi_np.request();
        py::buffer_info psi_buf     = psi_np.request();
        py::buffer_info DUstar_buf  = DUstar_np.request();
        double *push_mu             = static_cast<double *>(push_mu_buf.ptr);
        double *phi                 = static_cast<double *>(phi_buf.ptr);
        double *psi                 = static_cast<double *>(psi_buf.ptr);
        double *DUstar              = static_cast<double *>(DUstar_buf.ptr);
        // calculate_gradient(phi, gradx, grady);
        // pushforward.run_pushforward_on_density(DUstar,push_mu,gradx,grady);

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, phi);
        calculate_push_pull_rho(push_mu, DUstar, vxx_, vyy_, vxy_, psi, 0.99);
    }

    void compute_push_forth(
        py::array_t<double, py::array::c_style | py::array::forcecast> push_mu_np,
        py::array_t<double, py::array::c_style | py::array::forcecast> phi_np,
        py::array_t<double, py::array::c_style | py::array::forcecast> psi_np,
        py::array_t<double, py::array::c_style | py::array::forcecast> mu_np
    ){
        py::buffer_info push_mu_buf = push_mu_np.request();
        py::buffer_info phi_buf = phi_np.request();
        py::buffer_info psi_buf = psi_np.request();
        py::buffer_info mu_buf = mu_np.request();
        double *push_mu             = static_cast<double *>(push_mu_buf.ptr);
        double *phi                 = static_cast<double *>(phi_buf.ptr);
        double *psi                 = static_cast<double *>(psi_buf.ptr);
        double *mu                  = static_cast<double *>(mu_buf.ptr);
        // calculate_gradient(psi, gradx, grady);
        // pushforward.run_pushforward_on_density(mu,push_mu,gradx,grady);

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, psi);
        calculate_push_pull_rho(push_mu, mu, vxx_, vyy_, vxy_, phi, 0.2);
    }

    double calculate_F(const double lambda, const double* phi, const double* V, const double constant, const double M, const double mprime){
        double sum=0;

        for(int i=0;i<n1*n2;++i){
            double eval= - phi[i] - V[i] + lambda;
            if(eval>0){
                sum+=exp((mprime-1)*log(eval)); // IMPORTANT    
            }    
        }
        sum/=constant;

        return sum /(1.0*n1*n2)- M;
    }

    void calculate_DUstar(
            py::array_t<double, py::array::c_style | py::array::forcecast> DUstar_np,
            py::array_t<double, py::array::c_style | py::array::forcecast> phi_np,
            py::array_t<double, py::array::c_style | py::array::forcecast> V_np,
            double m, double M){
        
        py::buffer_info DUstar_buf = DUstar_np.request();
        py::buffer_info phi_buf = phi_np.request();
        py::buffer_info V_buf   = V_np.request();
        double *DUstar          = static_cast<double *>(DUstar_buf.ptr);
        double *phi             = static_cast<double *>(phi_buf.ptr);
        double *V               = static_cast<double *>(V_buf.ptr);

        const double tolerance=1e-4;
        int max_iteration=20;
        bool do_bisection = true;

        const double gamma = 1.0;

        double mprime = m/(m-1);

        double lambda_a=-phi[0]-V[0];
        double lambda_b=-phi[0]-V[0];

        double exponent=1.0/(m-1);

        double gammaprime = pow(gamma * mprime, mprime - 1);

        double lambda = 0;
        double val_at_0 = calculate_F(lambda,phi,V,gammaprime, M, mprime);

        if(fabs(val_at_0) < M * tolerance){
            do_bisection = false;
        }else if(val_at_0 > 0){
            lambda_b = 0;

            double t = -0.05*M;
            while(calculate_F(t,phi,V,gammaprime, M, mprime) > 0){
                t *= 2;
            }
            lambda_a = t;
        }else{
            lambda_a = 0;

            double t =  0.05*M;
            while(calculate_F(t,phi,V,gammaprime, M, mprime) < 0){
                t *= 2;
            }
            lambda_b = t;
        }
        
        if(do_bisection){
            lambda = 0.5 * (lambda_a + lambda_b);

            for(int iter=0;iter<max_iteration;++iter){
                double val = calculate_F(lambda,phi,V,gammaprime, M, mprime);

                if(fabs(val)<tolerance){
                    break;
                }

                if(val==0){
                    break;
                }else if(val<0){
                    lambda_a = lambda;
                }else{
                    lambda_b = lambda;
                }

                lambda = 0.5 * (lambda_a + lambda_b);
            }
        }

        for(int i=0;i<n1*n2;++i) phi[i] -= lambda;

        for(int i=0;i<n1*n2;++i){
            DUstar[i] = 1.0/gammaprime * pow(fmax(0, - phi[i] - V[i]), exponent);
        }
    }
}; // Back and Forth



#endif