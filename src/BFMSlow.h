/**
 * 
 * One outer loop of JKO scheme using the Back-and-Forth (BFM) method.
 * Slow Diffusion
 *
**/

#ifndef BFMSLOW_H
#define BFMSLOW_H

#include "BFM.h"

using namespace std;

class BackAndForthSlow:public BackAndForth{
public:
    
    double gamma_;
    double tau_;
    double m_;
    double mprime_;

    BackAndForthSlow(int n1, int n2, int max_iteration, double tolerance, double gamma, double tau, double m):BackAndForth(n1, n2, tau, max_iteration, tolerance){

        gamma_  = gamma;
        tau_    = tau;
        m_      = m;
        mprime_ = m/(m-1);

    }

    virtual ~BackAndForthSlow(){
    }

    void display_iteration(const int iter,const double W2_value,const double error_mu,const double error_nu, const double C_phi, const double C_psi, const double infgradphi) const{
        printf("%5d \tdual value: %5.6f\tL1 error: %5.4f\n", iter+1, W2_value,  error_mu);
        // printf("%5d C: %5.4f\tc1: %8.4f %8.4f\tc2: %8.4f %8.4f\tdual: %8.4f\tL1 error: %8.4f %8.4f\n", iter+1, C_phi, phi_c1_, psi_c1_, phi_c2_, psi_c2_, W2_value, error_mu, error_nu);
    }

    void set_coeff_m_2(double& c1, double& c2, const double mu_max, const double C_c_transform){
        c1 = 1/gamma_;
        c2 = C_c_transform * tau_;
    }

    void set_coeff(double& c1, double& c2, const double C, const double mu_max, const double d1, const double d2, const double d3, const bool verbose, const double C_c_transform){
        c1 = C_tr1_ * d1 + d2;
        c2 = C_tr2_ * d1 + C_c_transform * tau_;
    }

    void initialize_phi(Helper_U* helper_f,const double* mu, const int outer_iter){

        for(int i=0;i<n1*n2;++i){
            if(mu[i] > 0) push_mu_[i] = 0;
            else          push_mu_[i] = -LARGE_VALUE;
        }

        flt2d_->find_c_concave(push_mu_, push_mu_, 1);

        for(int i=0;i<n1*n2;++i){
            phi_[i] = - gamma_ * mprime_ * pow(mu[i],m_-1) - helper_f->V_[i];
            phi_[i] +=  pow(push_mu_[i], 0.5);
        }
    }

    void calculate_d1_d2_d3(double& d1, double& d2, double& d3, const double lambda, const double infgradphi, const double* nu){
        double eval = pow(gamma_ * mprime_, 1 - mprime_);
        d1 = eval * pow(lambda, mprime_-1) / infgradphi;
        d2 = eval * pow(lambda, mprime_-2) * (mprime_ - 1);
    }

    virtual void start_OT(Helper_U* helper_f, const double* mu, const int outer_iter, const bool verbose){

        const int skip = 50; // frequency of printout

        double error_mu = 1.0;
        double error_nu = 1.0;
        double error=1.0;

        /* Initialize coefficients for siga update */

        beta_1_ =0.1;
        beta_2_ =0.9;
        alpha_1_=1.05;
        alpha_2_=0.95;

        /* Calculate sup(mu) */

        double mu_max = 1;
        for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu[i]);

        cout << "Iter : " << outer_iter + 1 << "\n";

        /* Initialize the constants */

        double lambda = 1;
        double infgradphi  = 1;
        double sigma_forth = 1;
        double sigma_back  = 1;

        double d1 = 1;
        double d2 = 1;
        double d3 = 1;

        double C_c_transform = 1;

        initialize_phi(helper_f,mu,outer_iter);
        
        /* Starting the loop */
        
        for(int iter = 0; iter < max_iteration_; ++iter){


            /* Determinant version pushforward */

            sigma_forth = 1;
            sigma_back  = 1;


            if((iter) % 20 == 0 && iter >= 0){
                if(m_ == 2){
                    calculate_c_transform_constant(C_c_transform, phi_, mu);
                    set_coeff_m_2(phi_c1_, phi_c2_, mu_max, C_c_transform);
                }else{
                    lambda = calculate_lambda(helper_f);
                    infgradphi = calculate_infgradphi_on_level_set(lambda,helper_f->V_,helper_f->obstacle_,mu);
                    calculate_d1_d2_d3(d1, d2, d3, lambda, infgradphi, helper_f->V_);
                    calculate_trace_theorem_constant(helper_f, lambda, iter, mu);
                    calculate_c_transform_constant(C_c_transform, phi_, mu);
                    set_coeff(phi_c1_, phi_c2_, C_phi_, mu_max, d1, d2, d3, false, C_c_transform);
                }
            }
            
            error_mu = perform_OT_iteration_forth_det(helper_f,sigma_forth,W2_value_,mu,iter);

            if((iter) % 20 == 0 && iter >= 0){
                if(m_ == 2){
                    calculate_c_transform_constant(C_c_transform, psi_, helper_f->DUstar_);
                    set_coeff_m_2(psi_c1_, psi_c2_, mu_max, C_c_transform);    
                }else{
                    lambda = calculate_lambda(helper_f);
                    infgradphi = calculate_infgradphi_on_level_set(lambda,helper_f->V_,helper_f->obstacle_,mu);
                    calculate_d1_d2_d3(d1, d2, d3, lambda, infgradphi, helper_f->V_);
                    calculate_trace_theorem_constant(helper_f, lambda, iter, mu);
                    set_coeff(psi_c1_, psi_c2_, C_psi_, mu_max, d1, d2, d3, false, C_c_transform);    
                }
            }

            error_nu = perform_OT_iteration_back_det (helper_f,sigma_back, W2_value_,mu,iter);
            
            error=fmax(error_mu,error_nu);

            /*  Stopping Condition */

            if(((abs(error)<tolerance_ && abs(error)>0 && iter>=0) || iter==max_iteration_-1) ){
                display_iteration(iter,W2_value_,error_mu,error_nu,C_phi_,C_psi_,infgradphi);
                break;
            }

            /* Display the result per iterations */

            if(iter%skip==skip-1){
                display_iteration(iter,W2_value_,error_mu,error_nu,C_phi_,C_psi_,infgradphi);
                cout << flush;
            }
        }
    }


}; // Back and Forth Slow Diffusion

#endif