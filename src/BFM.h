/**
 * 
 * One outer loop of JKO scheme using the Back-and-Forth (BFM) method.
 *
 */

#ifndef BFM_H
#define BFM_H

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include "FLT.h"
#include "PoissonSolver.h"
#include "Helper_U.h"
#include "constants.h"

using namespace std;


class BackAndForth{
public:
    int n1;
    int n2;

    double tau_;

    int max_iteration_;
    double tolerance_;
    double W2_value_;

    double C_tr1_;
    double C_tr2_;

    double phi_c1_;
    double phi_c2_;
    double psi_c1_;
    double psi_c2_;

    double beta_1_;
    double beta_2_;
    double alpha_1_;
    double alpha_2_;

    double SCALER_forth_;
    double SCALER_back_ ;
    
    double C_phi_;
    double C_psi_;

    double* vxx_;
    double* vyy_;
    double* vxy_;

    double* phi_;
    double* psi_;

    double* push_mu_;

    double* R_arr_;

    poisson_solver* fftps_;
    FLT2D*          flt2d_;

    BackAndForth(int n1, int n2, double tau, int max_iteration, double tolerance){

        this->n1=n1;
        this->n2=n2;
        tau_ = tau;
        max_iteration_ =max_iteration;
        tolerance_     =tolerance;

        vxx_=new double[n1*n2];
        vyy_=new double[n1*n2];
        vxy_=new double[n1*n2];

        phi_=new double[n1*n2];
        psi_=new double[n1*n2];

        push_mu_= new double[n1*n2];

        flt2d_  = new FLT2D(n1,n2);

        R_arr_ = new double[n1];

        // set a constant for the trace theorem
        C_phi_ = 1;
        C_psi_ = 1;

        clock_t time;
        time=clock();
        fftps_ = new poisson_solver(n1,n2);
        time=clock()-time;
        printf ("\nCPU time for FFT: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);
    }

    virtual ~BackAndForth(){
        delete[] vxx_;
        delete[] vyy_;
        delete[] vxy_;
        delete[] phi_;
        delete[] psi_;
        delete[] push_mu_;
        delete[] R_arr_;

        delete flt2d_;
        delete fftps_;
    }

    void calculate_gradient(double* vx, double* vy, const double* phi_c){
        /* centered difference*/
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                vx[i*n1+j]=0.5*n1*(phi_c[i*n1+(int)fmin(n1-1,j+1)]-phi_c[i*n1+(int)fmax(0,j-1)]);
                vy[i*n1+j]=0.5*n2*(phi_c[(int)fmin(n2-1,i+1)*n1+j]-phi_c[(int)fmax(0,i-1)*n1+j]);
            }
        }
    }

    double calculate_dual_value(Helper_U* helper_f, const double* phi, const double* psi, const double* mu){

         // Here psi is assumed to correspond to the c-transform of phi
        int pcount=n1*n2;
        
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term1 += helper_f->calculate_Ustar(phi,i,j);
            }
        }

        for(int i=0;i<n1*n2;++i){
            term2+=psi[i]*mu[i];
        }
        
        return (term2 - term1) / pcount;
    }

    double calculate_h_minus_1(const double* push_mu, const double* DUstar_){
        double error=0;
        for(int i=0;i<n1*n2;++i){
            double value=-push_mu[i]+DUstar_[i];
            error+=value*fftps_->workspace[i];
        }
        return error/(1.0*n1*n2);
    }

    double calculate_L1_error(const double* push_mu, const double* mu){
        double error = 0;
        for(int i=0;i<n1*n2;++i) error += fabs(push_mu[i] - mu[i]);
        return error/(1.0*n1*n2);
    }

    void calculate_push_rho(double* push_rho, const double* rho, const double* vxx,const double* vyy,const double* vxy, const double* phi){

        double eps = pow(1.0/n1, 0.7);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

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
                    det = fmax(eps,det);
                    push_rho[i*n1+j] = rhovalue/det;
                    
                }else{
                    push_rho[i*n1+j]=0;
                }

            }
        }
    }

    void calculate_pull_rho(double* push_rho, const double* rho, const double* vxx, const double* vyy, const double* vxy, const double* phi){

        double eps = pow(1.0/n1, 0.7);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=0.5*n1*(phi[i*n1+(int)fmin(n1-1,j+1)]-phi[i*n1+(int)fmax(0,j-1)]);
                double vyval=0.5*n2*(phi[(int)fmin(n2-1,i+1)*n1+j]-phi[(int)fmax(0,i-1)*n1+j]);

                double x=(j+0.5)/(1.0*n1)-tau_*vxval;
                double y=(i+0.5)/(1.0*n2)-tau_*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){
                    double vxx_val = vxx[i*n1+j];
                    double vyy_val = vyy[i*n1+j];
                    double vxy_val = vxy[i*n1+j];

                    double det = fabs((1.0-tau_*vxx_val) * (1.0-tau_*vyy_val) - tau_*tau_ * vxy_val * vxy_val);
                    det = fmin(1.0/eps, det);
                    push_rho[i*n1+j] = rhovalue*det;
                    
                }else{
                    push_rho[i*n1+j]=0;
                }

            }
        }
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

    // This function will provide S1(x,) where x and y are n [0,1] double values
    double interpolate_function_v(double x,double y,const double* func){

        if(x>1 || x<0 || y>1 || y<0) return 0;

        double indj=fmin(n1-1,fmax(0,x*n1-0.5));
        double indi=fmin(n2-1,fmax(0,y*n2-0.5));

        double lambda1=indj-(int)indj;
        double lambda2=indi-(int)indi;

        double x00 = func[(int)(indi)*n1+(int)(indj)];
        double x01 = 0;
        double x10 = 0;
        double x11 = 0;
        if(indj+1 <= n1-1 && indi+1 <= n2-1){
            x11 = func[(int)(indi+1)*n1+(int)(indj+1)];
            x01 = func[(int)(indi)*n1+(int)(indj+1)];
            x10 = func[(int)(indi+1)*n1+(int)(indj)];
        }else if(indj+1 <= n1-1){
            x01 = func[(int)(indi)*n1+(int)(indj+1)];
        }else if(indi+1 <= n2-1){
            x10 = func[(int)(indi+1)*n1+(int)(indj)];
        }

        double interpolated_value = (1-lambda1)*(1-lambda2)*x00+(lambda1)*(1-lambda2)*x01
                                   +(1-lambda1)*(lambda2)  *x10+(lambda1)*(lambda2)  *x11;
        return interpolated_value;  
    }

    /* Update sigma based on Goldstein scheme */
    double update_sigma(const double sigma, const double W2_value, const double W2_value_previous, const double error){

        if(W2_value_previous-W2_value>-beta_1_*error){
            return alpha_2_;
        }else if(W2_value_previous-W2_value<-beta_2_*error){
            return alpha_1_;
        }

        if(W2_value - W2_value_previous < beta_1_*error){
            return alpha_2_;
        }else if(W2_value - W2_value_previous > beta_2_*error){
            return alpha_1_;
        }
        return 1;
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

    double perform_OT_iteration_back_det(Helper_U* helper_f,double& sigma,double& W2_value,const double* mu, const int iter){
        // ------------------------------------------------------------
        double W2_value_previous = 0;
        double error_nu = 0;

        flt2d_->find_c_concave(psi_,phi_,tau_);
        flt2d_->find_c_concave(phi_,psi_,tau_);
    
        helper_f->calculate_DUstar_normalized(phi_);    

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, phi_);

        // pushforward helper_f->DUstar_ -> push_mu_
        calculate_push_pull_rho(push_mu_, helper_f->DUstar_, vxx_, vyy_, vxy_, psi_, 0.99);

        fftps_->perform_inverse_laplacian(push_mu_, mu, psi_c1_, psi_c2_, sigma);

        W2_value_previous=calculate_dual_value(helper_f,phi_,psi_,mu);
        double error_nu_h=calculate_h_minus_1(push_mu_,mu);

        error_nu=calculate_L1_error(push_mu_,mu);

        for(int i=0;i<n1*n2;++i){
            psi_[i] += fftps_->workspace[i];
        }

        flt2d_->find_c_concave(psi_,phi_,tau_);

        W2_value=calculate_dual_value(helper_f,phi_,psi_,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu_h);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(Helper_U* helper_f,double& sigma,double& W2_value,const double* mu, const int iter){
        // ------------------------------------------------------------
        double W2_value_previous = 0;
        double error_mu = 0;

        flt2d_->find_c_concave(phi_,psi_,tau_);
        flt2d_->find_c_concave(psi_,phi_,tau_);

        helper_f->calculate_DUstar_normalized(phi_);    
            
        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, psi_);

        // pushforward mu -> push_mu_
        calculate_push_pull_rho(push_mu_, mu, vxx_, vyy_, vxy_, phi_, 0.2);

        fftps_->perform_inverse_laplacian(push_mu_, helper_f->DUstar_, phi_c1_, phi_c2_, sigma);

        W2_value_previous=calculate_dual_value(helper_f,phi_,psi_,mu);
        double error_mu_h=calculate_h_minus_1(push_mu_,helper_f->DUstar_);

        error_mu=calculate_L1_error(push_mu_, helper_f->DUstar_);

        for(int i=0;i<n1*n2;++i){
            phi_[i] += fftps_->workspace[i];
        }


        flt2d_->find_c_concave(phi_,psi_,tau_);

        W2_value=calculate_dual_value(helper_f,phi_,psi_,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu_h);

        return error_mu;
    }

    void display_iteration(const int iter,const double W2_value,const double error_mu,const double error_nu, const double C_phi, const double C_psi, const double infgradphi) const{
        cout << setprecision(6);
        cout << fixed;
        cout <<setw(5)<<iter+1 << " C : " << C_phi << "  infgradphi : " << infgradphi << "  c1 : " << phi_c1_ << " " << psi_c1_ << " coeff : " << setw(10) << phi_c2_/phi_c1_ << " " << setw(6) << psi_c2_/psi_c1_ << "  W2 : " << scientific << setw(13) << W2_value << "  L1 error : "<<scientific<<setw(13) << error_mu << " " << error_nu<<"\n";
    }

    /* Calculate a = 0.1 * max(-phi) */
    virtual double calculate_lambda(Helper_U* helper_f) const{

        double phimax = 1e-5;
        
        for(int i=0;i<n1*n2;++i){
            if(helper_f->obstacle_[i] == 0){
                phimax = fmax(phimax,-phi_[i]-helper_f->V_[i]);
            }
        }
        return fmax(0.1, phimax * 0.01);
    }

    /* Calculate infgradphi = inf(|nabla phi|) */
    virtual double calculate_infgradphi_on_level_set(const double lambda, const double* V, const double* obstacle, const double* mu, const double thres=20){
        double infgradphi= 10000;
        int count = 0;
        for(int i=0;i<n2-1;++i){
            for(int j=0;j<n1-1;++j){
                if(obstacle[i*n1+j] == 0 && obstacle[i*n1+j+1] == 0 && obstacle[(i+1)*n1+j] == 0){
                    if(-phi_[i*n1+j]-V[i*n1+j] > 0 && -phi_[i*n1+j]-V[i*n1+j] < lambda){
                        double gradxphi = n1*(-phi_[i*n1+(j+1)]+phi_[i*n1+j]-V[i*n1+(j+1)]+V[i*n1+j]);
                        double gradyphi = n2*(-phi_[(i+1)*n1+j]+phi_[i*n1+j]-V[(i+1)*n1+j]+V[i*n1+j]);
                        double eval = gradxphi*gradxphi + gradyphi*gradyphi;

                        if(count == 0) {
                            infgradphi = eval;
                            count = 1;
                        } else {
                            infgradphi = fmin(infgradphi, eval);
                        }
                        
                    }
                }
            }
        }

        return fmax(thres,sqrt(infgradphi));
    }

    virtual void calculate_c_transform_constant(double& C, const double* phi, const double* mu){

        C = 0;

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, phi);

        for(int i=1;i<n2-1;++i){
            for(int j=1;j<n1-1;++j){

                double vxval=0.5*n1*(phi[i*n1+(int)fmin(n1-1,j+1)]-phi[i*n1+(int)fmax(0,j-1)]);
                double vyval=0.5*n2*(phi[(int)fmin(n2-1,i+1)*n1+j]-phi[(int)fmax(0,i-1)*n1+j]);

                double x = (j+0.5)/n1 - tau_ * vxval;
                double y = (i+0.5)/n2 - tau_ * vyval;

                double mu_val = interpolate_function(x,y,mu);
                // double mu_val = 1;

                if(mu_val > 0){
                    /* calculate eigen values */
                        /* calculate trace */
                    double trace = fabs(2 - tau_*vxx_[i*n1+j] - tau_*vyy_[i*n1+j]);
                        /* calculate det */
                    double det   = fabs((1.0-tau_*vxx_[i*n1+j]) * (1.0-tau_*vyy_[i*n1+j]) - tau_*tau_ * vxy_[i*n1+j] * vxy_[i*n1+j]);

                    double t1 = 0.5 * fabs(trace + sqrt(fabs(trace*trace - 4 * det)));
                    double t2 = 0.5 * fabs(trace - sqrt(fabs(trace*trace - 4 * det)));
                    
                    C = fmax(C, mu_val*fmax(t1,t2));
                }
            }
        }
    }

    void sub_function_calculate_trace_constant_matt(double& C_tr1, double& C_tr2, Helper_U* helper_f, const double* push_mu, const double lambda){

        memset(R_arr_, 0, n1*sizeof(double));

        // calculate the laplacian
        for(int i=1;i<n2-1;++i){
            for(int j=1;j<n1-1;++j){
                if(push_mu_[i*n1+j]<0){
                    double laplacian = n1*n1*(push_mu[i*n1+(j+1)]-2.0*push_mu[i*n1+j]+push_mu[i*n1+(j-1)])
                                      +n2*n2*(push_mu[(i+1)*n1+j]-2.0*push_mu[i*n1+j]+push_mu[(i-1)*n1+j]);
                    int ind = - push_mu[i*n1+j]*n1;
                    if(ind >= 0 && ind < n1){
                        R_arr_[ind] = fmax(R_arr_[ind], fmax(0,laplacian));
                    }
                }
            }
        }
        
        for(int ri=4;ri<n1;++ri) R_arr_[ri] = fmax(R_arr_[ri], R_arr_[ri-1]);    

        double max_distance = 0; for(int i=0;i<n1*n2;++i) max_distance = fmax(max_distance, -push_mu_[i]);
        double R = LARGE_VALUE;
        for(int ri=3;ri<n1;++ri){
            if(1.0*ri/n1 > max_distance) break;
            R = fmin(R, R_arr_[ri] + 0.1*n1/(ri+1));
        }

        double C = fmax(R, 1);
        double eps = n1/2;
        C_tr1 = C / eps + C;
        C_tr2 = eps / C;
    }

    void calculate_trace_theorem_constant_matt(Helper_U* helper_f, double lambda, int iter, const double* mu){

        // outside 
        {
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double add = 0;
                    int indim = i;
                    if(helper_f->obstacle_[indim*n1+j] > 0) {add = 2*LARGE_VALUE;}
                    double evalim = - phi_[indim*n1+j] - helper_f->V_[indim*n1+j] - add;

                    add = 0;
                    int indjm = j;
                    if(helper_f->obstacle_[i*n1+indjm] > 0) {add = 2*LARGE_VALUE;}
                    double evaljm = - phi_[i*n1+indjm] - helper_f->V_[i*n1+indjm] - add;

                    add = 0;
                    int indi = fmin(n2-1,i+1);
                    if(helper_f->obstacle_[indi*n1+j] > 0) {add = 2*LARGE_VALUE;}
                    double evali = - phi_[indi*n1+j] - helper_f->V_[indi*n1+j] - add;

                    add = 0;
                    int indj = fmin(n1-1,j+1);
                    if(helper_f->obstacle_[i*n1+indj] > 0) {add = 2*LARGE_VALUE;}
                    double evalj = - phi_[i*n1+indj] - helper_f->V_[i*n1+indj] - add;

                    if((evalim>lambda && evali<-LARGE_VALUE) || (evali>lambda && evalim<-LARGE_VALUE) || (evaljm>lambda && evalj<-LARGE_VALUE) || (evalj>lambda && evaljm<-LARGE_VALUE))
                        push_mu_[i*n1+j] = -LARGE_VALUE;
                    else  if((evalim-lambda) * (evali-lambda) < 0 || (evaljm-lambda) * (evalj-lambda) < 0 ||
                        (evalim) * (evali) < 0 || (evaljm) * (evalj) < 0){
                        push_mu_[i*n1+j] = 0;
                    }else{
                        push_mu_[i*n1+j] = -LARGE_VALUE;
                    }

                }
            }

            flt2d_->find_c_concave(push_mu_, push_mu_, 0.5);

            for(int i=0;i<n1*n2;++i){
                double add = 0;
                if(helper_f->obstacle_[i] > 0) add = LARGE_VALUE;
                double eval = - phi_[i] - helper_f->V_[i] - add;
                if(eval < lambda && eval > 0) push_mu_[i] =  sqrt(push_mu_[i]);
                else push_mu_[i] = -sqrt(push_mu_[i]);
            }

            sub_function_calculate_trace_constant_matt(C_tr1_, C_tr2_, helper_f, push_mu_, lambda);
        }

        

        // inside
        {
            for(int i=0;i<n1*n2;++i) push_mu_[i] *= -1;

            double C_tr1=0;
            double C_tr2=0;

            sub_function_calculate_trace_constant_matt(C_tr1, C_tr2, helper_f, push_mu_, lambda);

            if(C_tr1 < C_tr1_) C_tr1_ = C_tr1;
            if(C_tr2 < C_tr2_) C_tr2_ = C_tr2;
        }
    }

    void calculate_trace_theorem_constant_flavien(){
        double K_gamma = 0.01;

        C_tr1_ = 12 * K_gamma;
        C_tr2_ = 6 * K_gamma;
    }

    virtual void calculate_trace_theorem_constant(Helper_U* helper_f, double lambda, int iter, const double* mu){
        // calculate_trace_theorem_constant_flavien();
        calculate_trace_theorem_constant_matt(helper_f,lambda, iter, mu);
    }

    virtual void start_OT(Helper_U* helper_f, const double* mu, const int outer_iter, const bool verbose) = 0;

}; // Back and Forth

#endif