/**
 * 
 * One outer loop of JKO scheme using the Back-and-Forth (BFM) method.
 * Incompressible Flow
 *
**/

#ifndef BFMINCOMP_H
#define BFMINCOMP_H

#include "BFM.h"

using namespace std;

class BackAndForthIncompressible:public BackAndForth{
public:
    
    double tau_;

    // Pushforward_mapping* pushforward;

    BackAndForthIncompressible(int n1, int n2, int max_iteration, double tolerance, double tau):BackAndForth(n1,n2,tau,max_iteration,tolerance){
        tau_=tau;
    }

    virtual ~BackAndForthIncompressible() {}

    double calculate_WENO(const double a, const double b, const double c, const double d){

        /* setting constants for WENO. The constants need to satisfy w0 + w1 + w2 = 1 */

        double eps = 1e-6;
        double IS0 = 13.0 * (a-b)*(a-b) + 3.0 * (a-3*b)*(a-3*b);
        double IS1 = 13.0 * (b-c)*(b-c) + 3.0 * (b+c)*(b+c);
        double IS2 = 13.0 * (c-d)*(c-d) + 3.0 * (3*c-d)*(3*c-d);
        double alpha0 = 1.0/((eps + IS0)*(eps + IS0));
        double alpha1 = 6.0/((eps + IS1)*(eps + IS1));
        double alpha2 = 3.0/((eps + IS2)*(eps + IS2));
        double w0 = alpha0 / (alpha0 + alpha1 + alpha2);
        double w2 = alpha2 / (alpha0 + alpha1 + alpha2);

        /* return the value */
        return w0/3.0 * (a - 2.0 * b + c) + (w2 - 0.5)/6.0 * (b - 2.0 * c + d);
    }

    double calculate_gradx_WENO(const double* phi, const int i, const int j){
        double phi_j_pre_2  = phi[i*n1+ (int)fmax(0,j-2)];
        double phi_j_pre_1  = phi[i*n1+ (int)fmax(0,j-1)];
        double phi_j        = phi[i*n1+j];
        double phi_j_post_1 = phi[i*n1+ (int)fmin(n1-1,j+1)];
        double phi_j_post_2 = phi[i*n1+ (int)fmin(n1-1,j+2)];
        double phi_j_post_3 = phi[i*n1+ (int)fmin(n1-1,j+3)];
        double eval = n1/12.0 * (
                                 -         (phi_j_pre_1  - phi_j_pre_2)
                                 + 7.0  *  (phi_j        - phi_j_pre_1)
                                 + 7.0  *  (phi_j_post_1 - phi_j)
                                 -         (phi_j_post_2 - phi_j_post_1)
                                );
        double a = n1 * (phi_j_post_3 - 2.0 * phi_j_post_2 + phi_j_post_1);
        double b = n1 * (phi_j_post_2 - 2.0 * phi_j_post_1 + phi_j);
        double c = n1 * (phi_j_post_1 - 2.0 * phi_j + phi_j_pre_1);
        double d = n1 * (phi_j - 2.0 * phi_j_pre_1  + phi_j_pre_2);
        eval += calculate_WENO(a,b,c,d);
        return eval;
    }

    double calculate_grady_WENO(const double* phi, const int i, const int j){
        double phi_i_pre_2  = phi[(int)fmax(0,i-2)*n1+j];
        double phi_i_pre_1  = phi[(int)fmax(0,i-1)*n1+j];
        double phi_i        = phi[i*n1+j];
        double phi_i_post_1 = phi[(int)fmin(n2-1,i+1)*n1+j];
        double phi_i_post_2 = phi[(int)fmin(n2-1,i+2)*n1+j];
        double phi_i_post_3 = phi[(int)fmin(n2-1,i+3)*n1+j];
        double eval = n2/12.0 * (
                                 -         (phi_i_pre_1  - phi_i_pre_2)
                                 + 7.0  *  (phi_i        - phi_i_pre_1)
                                 + 7.0  *  (phi_i_post_1 - phi_i)
                                 -         (phi_i_post_2 - phi_i_post_1)
                                );
        double a = n2 * (phi_i_post_3 - 2.0 * phi_i_post_2 + phi_i_post_1);
        double b = n2 * (phi_i_post_2 - 2.0 * phi_i_post_1 + phi_i);
        double c = n2 * (phi_i_post_1 - 2.0 * phi_i + phi_i_pre_1);
        double d = n2 * (phi_i - 2.0 * phi_i_pre_1  + phi_i_pre_2);
        eval += calculate_WENO(a,b,c,d);
        return eval;
    }

    double calculate_gradx_WENO_(const double* phi, const int i, const int j){
        double phi_j_pre_3  = phi[i*n1+ (int)fmax(0,j-3)];
        double phi_j_pre_2  = phi[i*n1+ (int)fmax(0,j-2)];
        double phi_j_pre_1  = phi[i*n1+ (int)fmax(0,j-1)];
        double phi_j        = phi[i*n1+j];
        double phi_j_post_1 = phi[i*n1+ (int)fmin(n1-1,j+1)];
        double phi_j_post_2 = phi[i*n1+ (int)fmin(n1-1,j+2)];
        double eval = n1/12.0 * (
                                 -         (phi_j_pre_1  - phi_j_pre_2)
                                 + 7.0  *  (phi_j        - phi_j_pre_1)
                                 + 7.0  *  (phi_j_post_1 - phi_j)
                                 -         (phi_j_post_2 - phi_j_post_1)
                                );
        double a = n1 * (phi_j_pre_1 - 2.0 * phi_j_pre_2 + phi_j_pre_3);
        double b = n1 * (phi_j - 2.0 * phi_j_pre_1 + phi_j_pre_2);
        double c = n1 * (phi_j_post_1 - 2.0 * phi_j + phi_j_pre_1);
        double d = n1 * (phi_j_post_2 - 2.0 * phi_j_post_1 + phi_j);
        eval -= calculate_WENO(a,b,c,d);
        return eval;

    }

    void display_iteration(const int iter,const double W2_value_,const double error_mu,const double error_nu,const double solution_error, const double C_phi_, const double C_psi_) const{
        printf("%5d \tdual: %5.4f\tL1 error: %5.4f\n", iter+1, W2_value_,  error_mu);
    }

    /**
        Calculate infgradphi = inf(|nabla phi|)
    */
    double calculate_infgradphi_on_level_set(const double lambda,const double* nu, Helper_U* helper_f){

        flt2d_->find_c_concave(phi_,psi_,tau_);
        flt2d_->find_c_concave(psi_,phi_,tau_);

        double infgradphi= 1;
        int count = 0;
        for(int i=1;i<n2-1;++i){
            for(int j=1;j<n1-1;++j){
                if(helper_f->obstacle_[(i-1)*n1+j] == 0 && helper_f->obstacle_[i*n1+j-1] == 0 && helper_f->obstacle_[i*n1+j+1] ==0 && helper_f->obstacle_[(i+1)*n1+j] == 0){
                    if( (-phi_[(i-1)*n1+j]-nu[(i-1)*n1+j]) * (-phi_[(i+1)*n1+j]-nu[(i+1)*n1+j]) < 0 ||  (-phi_[i*n1+j-1]-nu[i*n1+j-1]) * (-phi_[i*n1+j+1]-nu[i*n1+j+1]) < 0){
                        double gradxphi = 0.5 * n1 * (-phi_[i*n1+j+1]+phi_[i*n1+j-1]-nu[i*n1+j+1]+nu[i*n1+j-1]);
                        double gradyphi = 0.5 * n2 * (-phi_[(i+1)*n1+j]+phi_[(i-1)*n1+j]-nu[(i+1)*n1+j]+nu[(i-1)*n1+j]);
                        double eval = gradxphi*gradxphi + gradyphi*gradyphi;
                        if(count == 0){
                            infgradphi = eval;
                            count = 1;
                        }else{
                            infgradphi = fmin(infgradphi, eval);    
                        }
                    }    
                }
            }
        }

        // return infgradphi;
        return fmax(0.1,sqrt(infgradphi));
    }


    void set_coeff(double& c1, double& c2, const double C, const double mu_max, const double d1, const double d2, const bool verbose, const double C_tr){
        c1 = 10 * C * d1;
        c2 = C * d1 + C_tr * tau_;
    } 

    void initialize_phi(Helper_U* helper_f,const double* mu){
        for(int i=0;i<n1*n2;++i){
            if(mu[i] > 0) push_mu_[i] = 0;
            else          push_mu_[i] = -LARGE_VALUE;
        }

        flt2d_->find_c_concave(push_mu_, push_mu_, 0.5);

        for(int i=0;i<n1*n2;++i){
            phi_[i] = - mu[i] - helper_f->V_[i];
            phi_[i] +=  pow(push_mu_[i],0.5);
        }
    }

    void calculate_d1_d2(double& d1, double& d2, const double lambda, const double infgradphi){
        d1 = 1.0 / infgradphi;
        d2 = 10 / infgradphi;
    }

    void start_OT(Helper_U* helper_f, const double* mu, const int outer_iter, const bool verbose){

        int skip = 10; // frequency of printout

        double error_mu,error_nu,error=1.0;

        /* Initialize coefficients for siga update */

        beta_1_ =0.2;
        beta_2_ =0.8;
        alpha_1_=1.05;
        alpha_2_=0.95;

        /* Initialize the tolerance_ based on tau_^2 */

        double mu_max = 1;

        double tol_modified = tolerance_;

        double mid_tolerance_ = 5e-2;

        cout << "Iter : " << outer_iter + 1 << " Tolerance : " << tol_modified << "\n";

        /* Initialize the coefficients for fftps */

        double lambda = 1;
        double infgradphi  = 1;
        double sigma_forth = 1;
        double sigma_back  = 1;

        double d1 = 1;
        double d2 = 1;

        double C_tr = 1;

        C_phi_ = 1;
        C_psi_ = 1;

        if(outer_iter == 0) initialize_phi(helper_f,mu); // intiailize phi in the first outer iteration

        double solution_error = 1;


        if(outer_iter == 0){
            phi_c1_ = 20;
            psi_c1_ = 20;

            phi_c2_ = 1;
            psi_c2_ = 1;
        }
        
        /* Starting the loop */

        for(int iter=0;iter<max_iteration_;++iter){

            double high_thres = 0.1;
            double low_thres  = 0.05;
            C_phi_ = fmin(high_thres, fmax(low_thres, C_phi_ / sigma_forth));
            C_psi_ = fmin(high_thres, fmax(low_thres, C_psi_ / sigma_back));

            if(iter % 10 == 0){
                infgradphi = calculate_infgradphi_on_level_set(0,helper_f->V_,helper_f);
                calculate_d1_d2(d1, d2, 0, infgradphi);
                calculate_c_transform_constant(C_tr, phi_, mu);

                set_coeff(phi_c1_, phi_c2_, C_phi_, mu_max, d1, d2, false, C_tr);
                set_coeff(psi_c1_, psi_c2_, C_psi_, mu_max, d1, d2, false, C_tr);
            }

            /* Determinant version pushforward */

            sigma_forth = 1;
            sigma_back  = 1;
            
            error_mu = perform_OT_iteration_forth_det(helper_f,sigma_forth,W2_value_,mu,iter);
            error_nu = perform_OT_iteration_back_det(helper_f,sigma_back,W2_value_,mu,iter);

            /* Calculating the relative error */
            
            error_mu /= helper_f->M_;
            error_nu /= helper_f->M_;            
            error=fmax(error_mu,error_nu);
            
            /* Stopping Condition */

            if(((fabs(error ) <tol_modified && abs(error)>0 && iter>=0) || iter==max_iteration_-1) ){
                display_iteration(iter,W2_value_,error_mu,error_nu,solution_error,C_phi_,C_psi_);
                cout<<"Tolerance met!"<<endl;
                cout << flush;
                break;
            }

            /* Display the result per iterations */
            
            if(iter%skip==skip-1){
                /* Compare with actual solution */
                display_iteration(iter,W2_value_,error_mu,error_nu,solution_error,C_phi_,C_psi_);
                cout << flush;
            }
        }
    }


}; // Back and Forth Incompressible Flow

#endif
