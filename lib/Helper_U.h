#ifndef HELPER_U_H
#define HELPER_U_H

#include <iostream>

class Helper_U{
public:
    int n1_;
    int n2_;

    double gamma_;
    double tau_;
    double m_;
    double mprime_;

    double M_; // the initial mass of given density.

    unsigned char* obstacle_;
    double* V_;
    double* DUstar_;

    Helper_U(){
        V_=NULL;
        DUstar_=NULL;
        obstacle_ = NULL;
    }

    Helper_U(int n1,int n2,double gamma,double tau,double m,const double* mu){
        n1_    = n1;
        n2_    = n2;
        gamma_ = gamma;
        tau_   = tau;
        m_     = m;
        mprime_= m_/(m_-1);

        DUstar_ =new double[n1*n2];
        V_      =new double[n1*n2];

        set_M_(mu);
    }

    ~Helper_U(){
        delete[] DUstar_;
        delete[] V_;
        delete[] obstacle_;
    }

    // Keeping the mass of the initial density
    void set_M_(const double* mu){
    	M_ = 0;
    	for(int i=0;i<n1_*n2_;++i) M_ += mu[i];
    	M_ /= n1_*n2_;
    }

    void set_obstacle(unsigned char* obstacle){
        obstacle_ = obstacle;
    }

    double calculate_U(const double* phi, const int i, const int j) const
    {
        if(V_[i*n1_+j] < 0) return 0;

        double c_phi = - phi[i*n1_+j] - V_[i*n1_+j];

        if(c_phi > 0 && obstacle_[i*n1_+j] == 0){
        	double constant = 1.0 / (pow(gamma_ * mprime_, mprime_ - 1) * mprime_);
            return constant * exp(mprime_*log(c_phi));    
        }
        return 0;
    }

    void calculate_DUstar_(const double* phi){
        
        double exponent = 1.0/(m_-1);
        double gammaprime = pow(gamma_ * mprime_, exponent);

        for(int i=0;i<n1_*n2_;++i){
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

        for(int i=0;i<n1_*n2_;++i){
            if(obstacle_[i] == 0){
                double eval=- phi[i] - V_[i] + lambda;
                if(eval>0){
                    sum+=exp((mprime_-1)*log(eval)); // IMPORTANT    
                }    
            }
        }
        sum/=constant;

        return sum /(1.0*n1_*n2_)- M_;
    }

    void calculate_DUstar_normalized(double* phi, const double tolerance=1e-7){
        
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

        for(int i=0;i<n1_*n2_;++i) phi[i] -= lambda;

        for(int i=0;i<n1_*n2_;++i){
            if(obstacle_[i] == 0){
                DUstar_[i] = 1.0/gammaprime * pow(fmax(0, - phi[i] - V_[i]), exponent);
            }else{
                DUstar_[i] = 0;
            }
        }
    }

}; // Helper_U


#endif