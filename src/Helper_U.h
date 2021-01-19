/**
 * 
 * A class that calculates the energy functional U and delta U.
 *
 */

#ifndef HELPER_U_H
#define HELPER_U_H

#include <iostream>

class Helper_U{
public:
    int n1;
    int n2;

    double M_; // the initial mass of given density.

    double* obstacle_;
    double* V_;
    double* DUstar_;

    Helper_U(){
        V_        =NULL;
        DUstar_   =NULL;
        obstacle_ =NULL;
    }

    Helper_U(int n1,int n2,const double* mu){
        this->n1 = n1;
        this->n2 = n2;

        V_        =NULL;
        DUstar_   =new double[n1*n2];
        obstacle_ =NULL;

        set_M_(mu);

        cout << "Total mass : " << M_ << "\n";
    }

    virtual ~Helper_U(){
        delete[] DUstar_;
    }

    // Keeping the mass of the initial density
    void set_M_(const double* mu){
        M_ = 0;
        for(int i=0;i<n1*n2;++i) M_ += mu[i];
        M_ /= n1*n2;
    }

    void set_obstacle(double* obstacle){
        obstacle_ = obstacle;
    }
    void set_V(double* V){
        V_ = V;
    }

    virtual double calculate_Ustar(const double* phi, const int i, const int j) const = 0;
    virtual double calculate_F(const double lambda, const double* phi, const double constant) = 0;
    virtual void calculate_DUstar_(const double* phi) = 0;
    virtual void calculate_DUstar_normalized(double* phi, const double tolerance=1e-6) = 0;

}; // Helper_U

class Helper_U_slow: public Helper_U{
public:
    double gamma_;
    double tau_;
    double m_;
    double mprime_;

    Helper_U_slow(): Helper_U() {}

    Helper_U_slow(int n1,int n2,double gamma,double tau,double m,const double* mu): Helper_U(n1, n2, mu){
        gamma_ = gamma;
        tau_   = tau;
        m_     = m;
        mprime_= m_/(m_-1);
    }

    virtual ~Helper_U_slow() {}

    double calculate_Ustar(const double* phi, const int i, const int j) const
    {
        if(V_[i*n1+j] < 0) return 0;

        double c_phi = - phi[i*n1+j] - V_[i*n1+j];

        if(c_phi > 0 && obstacle_[i*n1+j] == 0){
            double constant = 1.0 / (pow(gamma_ * mprime_, mprime_ - 1) * mprime_);
            return constant * exp(mprime_*log(c_phi));    
        }
        return 0;
    }

    void calculate_DUstar_(const double* phi){
        
        double exponent = 1.0/(m_-1);
        double gammaprime = pow(gamma_ * mprime_, exponent);

        for(int i=0;i<n1*n2;++i){
            // DUstar_[i] = 1.0/gamma_ * pow(fmax(0, - phi[i]/tau - V_[i]), 1.0/(m-1));
            if(obstacle_[i] == 0){
                DUstar_[i] = 1.0/gammaprime * exp(exponent * log(fmax(0, - phi[i] - V_[i])));    
            }else{
                DUstar_[i] = 0;
            }
        }
    }

    double calculate_F(const double lambda, const double* phi, const double constant){
        double sum=0;

        for(int i=0;i<n1*n2;++i){
            if(obstacle_[i] == 0){
                double eval=- phi[i] - V_[i] + lambda;
                if(eval>0){
                    sum+=exp((mprime_-1)*log(eval)); // IMPORTANT    
                }    
            }
        }
        sum/=constant;

        return sum /(1.0*n1*n2)- M_;
    }

    void calculate_DUstar_normalized(double* phi, const double tolerance=1e-6){
        
        int max_iteration=100;
        bool do_bisection = true;

        double lambda_a=-phi[0]-V_[0];
        double lambda_b=-phi[0]-V_[0];

        double exponent=1.0/(m_-1);

        double gammaprime = pow(gamma_ * mprime_, mprime_ - 1);

        double lambda = 0;
        double val_at_0 = calculate_F(lambda,phi,gammaprime);

        if(fabs(val_at_0) < M_ * tolerance){
            do_bisection = false;
        }else if(val_at_0 > 0){
            lambda_b = 0;

            double t = -0.05*M_;
            while(calculate_F(t,phi,gammaprime) > 0){
                t *= 2;
            }
            lambda_a = t;
        }else{
            lambda_a = 0;

            double t =  0.05*M_;
            while(calculate_F(t,phi,gammaprime) < 0){
                t *= 2;
            }
            lambda_b = t;
        }
        
        if(do_bisection){
            lambda = 0.5 * (lambda_a + lambda_b);

            for(int iter=0;iter<max_iteration;++iter){
                double val = calculate_F(lambda,phi,gammaprime);

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
            if(obstacle_[i] == 0){
                DUstar_[i] = 1.0/gammaprime * pow(fmax(0, - phi[i] - V_[i]), exponent);
            }else{
                DUstar_[i] = 0;
            }
        }
    }

}; // Helper_U




class Helper_U_incompressible: public Helper_U{
public:
    double tau_;

    Helper_U_incompressible(): Helper_U() {}
    Helper_U_incompressible(int n1,int n2,double tau,double* mu): Helper_U(n1,n2,mu){
        tau_=tau;
    }

    virtual ~Helper_U_incompressible(){
    }

    virtual double calculate_Ustar(const double* phi, const int i, const int j) const
    {
        double c_phi = - phi[i*n1+j] - V_[i*n1+j];
        if(obstacle_[i*n1+j] == 0){
            if(c_phi >= 0 ){
                return c_phi;
            }    
        }
        return 0;
    }

    virtual void calculate_DUstar_(const double* phi){

        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - V_[i];
            if(obstacle_[i] == 0){
                if(eval>0){
                    DUstar_[i] = 1;
                }else{
                    DUstar_[i] = 0;
                }    
            }else{
                DUstar_[i] = 0;
            }
        }
    }

    virtual double calculate_F(const double lambda, const double* phi, const double constant = 1){
        double sum=0;

        for(int i=0;i<n1*n2;++i){
            double eval =- phi[i] - V_[i] + lambda;
            if(obstacle_[i] == 0){
                if(eval>0){
                    sum += 1.0; // IMPORTANT    
                }    
            }
        }

        return sum /(1.0*n1*n2) - M_;
    }

    virtual void calculate_DUstar_normalized(double* phi, const double tolerance=1e-3){
        
        double TOL = 1e-3;

        double lambda_a=-phi[0]-V_[0];
        double lambda_b=-phi[0]-V_[0];

        bool do_bisection = true;

        double lambda = 0;

        double val_at_0 = calculate_F(lambda,phi);
        if(fabs(val_at_0) < M_ * TOL){
            do_bisection = false;
        }else if(val_at_0 > 0){
            lambda_b = 0;

            double t = -TOL*M_;
            while(calculate_F(t,phi) > 0){
                t *= 2;
            }
            lambda_a = t;
        }else{
            lambda_a = 0;

            double t =  TOL*M_;
            while(calculate_F(t,phi) < 0){
                t *= 2;
            }
            lambda_b = t;
        }

        if(do_bisection){
            int max_iteration=50;

            lambda = 0.5 * (lambda_a + lambda_b);

            for(int iter=0;iter<max_iteration;++iter){
                double val = calculate_F(lambda, phi);

                if(fabs(val) < M_ * TOL || iter == max_iteration-1){
                    break;
                }

                if(val<0){
                    lambda_a = lambda;
                }else{
                    lambda_b = lambda;
                }

                lambda = 0.5 * (lambda_a + lambda_b);
            }
        }

        // update phi
        for(int i=0;i<n1*n2;++i){
            phi[i] -= lambda;
        } 

        // update DUstar_
        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - V_[i];
            if(obstacle_[i]==0){
                if(eval>0){
                    DUstar_[i] = 1.0;
                }else{
                    DUstar_[i] = 0.0;
                }
            }else{
                DUstar_[i] = 0;
            }
        }
    }

    void calculate_DUstar_normalized_bisection(double* phi, const double tolerance=1e-11){
        
        double lambda_a=-phi[0]-V_[0];
        double lambda_b=-phi[0]-V_[0];

        for(int i=0;i<n1*n2;++i){
            lambda_a = fmax(lambda_a, - phi[i] - V_[i]);
            lambda_b = fmin(lambda_b, - phi[i] - V_[i]);
        }

        lambda_a = - lambda_a - 10;
        lambda_b = - lambda_b + 10;

        int max_iteration=1000;

        double lambda = 0.5 * (lambda_a + lambda_b);

        for(int iter=0;iter<max_iteration;++iter){
            double sum=0;

            for(int i=0;i<n1*n2;++i){
                double eval =- phi[i] - V_[i] + lambda;
                if(V_[i] >= 0){
                    if(eval>0){
                        sum += 1.0; // IMPORTANT    
                    }    
                }
            }

            double val =sum /(1.0*n1*n2) - M_;

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

        // update phi
        for(int i=0;i<n1*n2;++i) phi[i] -= lambda;

        // update DUstar_
        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - V_[i];
            if(V_[i]>=0){
                if(eval>0){
                    DUstar_[i] = 1.0;
                }else{
                    DUstar_[i] = 0.0;
                }
            }else{
                DUstar_[i] = 0;
            }
        }
    }

}; // Helper_U_incompressible



#endif