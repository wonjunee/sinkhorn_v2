#ifndef BFM_H
#define BFM_H

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include "Accessory.h"
#include "Plotting.h"

using namespace std;

const double LARGE_VALUE = 999999;

class BackAndForth{
public:
    int n1;
    int n2;

    int n_mu;
    int n_nu;

    int DIM_;

    int max_iteration_;
    double tolerance_;
    double W2_value_;
    double W2_value_back_;

    double dual_forth_;
    double dual_back_;

    double sigma_; // Coefficients for trace theorem
    double sigma_forth_;
    double sigma_back_;

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

    double* phi_g_;
    double* psi_f_;

    Points*  push_mu_;
    Points*  push_nu_;

    int* push_mu_idx_;
    int* push_nu_idx_;

    double* C_mat_; // cost matrix

    double* C_mu_;  // cost matrix for mu
    double* C_nu_;  // cost matrix for mu

    double* C_mu_sum_;
    double* C_nu_sum_;

    double epsilon_;

    poisson_solver* fftps_;

    BackAndForth(int n1, int n2, int n_mu, int n_nu, int max_iteration, double tolerance, double sigma){

        epsilon_ = 1e-8;

        this->n1=n1;
        this->n2=n2;

        this->n_mu=n_mu;
        this->n_nu=n_nu;

        DIM_ = 2;
        max_iteration_ =max_iteration;
        tolerance_     =tolerance;
        sigma_ = sigma;

        printf("n_mu : %d n_mu : %d", n_mu, n_nu);

        phi_g_=new double[n_nu];
        psi_f_=new double[n_mu];

        memset(phi_g_, 0, n_nu*sizeof(double));
        memset(psi_f_, 0, n_mu*sizeof(double));

        C_mat_ = new double[n_mu * n_nu];
        C_mu_  = new double[n_mu * n_mu];
        C_nu_  = new double[n_nu * n_nu];

        C_mu_sum_  = new double[n_mu];
        C_nu_sum_  = new double[n_nu];

        push_mu_idx_ = new int[n_mu];
        push_nu_idx_ = new int[n_nu];

        // set a constant for the trace theorem
        C_tr1_ = 1;
        C_tr2_ = 1;

        phi_c1_ = 0;
        psi_c1_ = 0;

        phi_c2_ = 200;
        psi_c2_ = 200;

        clock_t time;
        time=clock();
        fftps_ = new poisson_solver(n1,n2);
        time=clock()-time;
        printf ("\nCPU time for FFT: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);
    }

    ~BackAndForth(){

        printf("bfm deconstruction\n");
        delete[] C_mat_;

        printf("there\n");

        delete[] C_mu_sum_;
        delete[] C_nu_sum_;
        
        delete[] push_mu_idx_;
        delete[] push_nu_idx_;

        delete fftps_;

        printf("bfm deconstruction done\n");
    }

    void Initialize_C_mat_(Points* mu, Points* nu){
        /** 
         *       Create a cost function 
         * C_mat_ is a matrix C (num_points_nu x num_points_mu)
         * C_{ij} is a distance between mu_j and nu_i
        **/

        for(int i=0;i<n_nu;++i){
            for(int j=0;j<n_mu;++j){
                double mux = (*mu)(j,0);
                double muy = (*mu)(j,1);
                double nux = (*nu)(i,0);
                double nuy = (*nu)(i,1);
                double px  = mux - nux;
                double py  = muy - nuy;
                C_mat_[i*n_mu+j] = 0.5 * (px*px + py*py);
            }
        }

        for(int i=0;i<n_nu;++i){
            double nux = (*nu)(i,0);
            double nuy = (*nu)(i,1);
            for(int j=i;j<n_nu;++j){
                double mux = (*nu)(j,0);
                double muy = (*nu)(j,1);
                double px  = mux - nux;
                double py  = muy - nuy;
                C_nu_[i*n_nu+j] = 0.5 * (px*px + py*py);
                C_nu_[j*n_nu+i] = C_nu_[i*n_nu+j];
            }
        }

        for(int i=0;i<n_mu;++i){
            double nux = (*mu)(i,0);
            double nuy = (*mu)(i,1);
            for(int j=i;j<n_mu;++j){
                double mux = (*mu)(j,0);
                double muy = (*mu)(j,1);
                double px  = mux - nux;
                double py  = muy - nuy;
                C_mu_[i*n_mu+j] = 0.5 * (px*px + py*py);
                C_mu_[j*n_mu+i] = C_mu_[i*n_mu+j];
            }
        }

        for(int j=0;j<n_mu;++j){
            double val = 0;
            for(int j1=0;j1<n_mu;++j1){
                double eval = sqrt(C_mu_[j*n_mu+j1]*2);
                val += - 1.0/(2.0*M_PI) * log(epsilon_ + eval);
            }
            C_mu_sum_[j] = val;
        }

        for(int i=0;i<n_nu;++i){
            double val = 0;
            for(int i1=0;i1<n_nu;++i1){
                double eval = sqrt(C_nu_[i*n_nu+i1]*2);
                val += - 1.0/(2.0*M_PI) * log(epsilon_ + eval);
            }
            C_nu_sum_[i] = val;
        }
    }

    void compute_psi_c_transform(double* phi_g_, int* push_nu_idx_, const double* psi_f_, Points* mu, Points* nu){
        for(int i=0;i<n_nu;++i){
            double eval = C_mat_[i*n_mu+0] - psi_f_[0];
            double min_val = eval;
            int    min_idx = 0;
            for(int j=1;j<n_mu;++j){
                eval = C_mat_[i*n_mu+j] - psi_f_[j];
                if(eval < min_val){
                    min_val = eval;
                    min_idx = j;
                }
            }
            phi_g_[i] = min_val;

            push_nu_idx_[i] = min_idx;
        }
    }

    void compute_phi_c_transform(double* psi_f_, int* push_mu_idx_, double* phi_g_, Points* mu, Points* nu){
        for(int j=0;j<n_mu;++j){
            double eval = C_mat_[0*n_mu+j] - phi_g_[0];
            double min_val = eval;
            int    min_idx = 0;
            for(int i=1;i<n_nu;++i){
                eval = C_mat_[i*n_mu+j] - phi_g_[i];
                if(eval < min_val){
                    min_val = eval;
                    min_idx = i;
                }
            }
            psi_f_[j] = min_val;

            push_mu_idx_[j] = min_idx;
        }
    }

    double calculate_dual_value(const double* phi, const double* psi, const double* mu, const double* nu){

         // Here psi is assumed to correspond to the c-transform of phi
        
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term1 += phi[i*n1+j]*nu[i*n1+j];
                term2 += psi[i*n1+j]*mu[i*n1+j];
            }
        }
        
        
        return (term2 + term1)/(1.0*n1*n2);
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


    /**
        Update sigma based on Goldstein scheme
    */
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

    double calculate_dual_value(const double* phi_g_, const double* psi_f_, Points* mu, Points* nu){
        double sum = 0;
        for(int i=0;i<n_nu;++i){
            sum += phi_g_[i];
        }
        for(int j=0;j<n_mu;++j){
            sum += psi_f_[j];
        }
        return sum;
    }

    double calculate_W2_value(int* push_mu_idx_){
        double sum = 0;
        for(int j=0;j<n_mu;++j){
            sum += C_mat_[push_mu_idx_[j]*n_mu+j];
        }
        return sum;
    }

    double perform_OT_iteration_back_det(double& sigma, double& W2_value, double& dual_value, Points* mu, Points* nu, const int iter){
        double W2_value_previous = 0;
        double error_nu = 0;

        compute_phi_c_transform(psi_f_, push_mu_idx_, phi_g_, mu, nu);
        compute_psi_c_transform(phi_g_, push_nu_idx_, psi_f_, mu, nu);


        for(int j=0;j<n_mu;++j){

            double update_val = sigma*C_mu_sum_[j];

            for(int i=0;i<n_nu;++i){
                double eval = sqrt(C_mu_[j*n_mu+push_nu_idx_[i]]*2);
                update_val += sigma/(2.0*M_PI) * log(epsilon_ + eval);
            }
            psi_f_[j] += update_val;
        }

        compute_psi_c_transform(phi_g_, push_nu_idx_, psi_f_, mu, nu);
        compute_phi_c_transform(psi_f_, push_mu_idx_, phi_g_, mu, nu);

        dual_value = calculate_dual_value(phi_g_, psi_f_, mu, nu);
        W2_value   = calculate_W2_value(push_mu_idx_);

        // sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu_h);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(double& sigma, double& W2_value, double& dual_value, Points* mu, Points* nu, const int iter){
        double W2_value_previous = 0;
        double error_mu = 1;

        compute_psi_c_transform(phi_g_, push_nu_idx_, psi_f_, mu, nu);
        compute_phi_c_transform(psi_f_, push_mu_idx_, phi_g_, mu, nu);


        for(int i=0;i<n_nu;++i){

            double update_val = sigma*C_nu_sum_[i];

            for(int j=0;j<n_mu;++j){
                double eval = sqrt(C_nu_[i*n_nu+push_mu_idx_[j]]*2);
                update_val += sigma/(2.0*M_PI) * log(epsilon_ + eval);
            }
            phi_g_[i] += update_val;
        }

        compute_phi_c_transform(psi_f_, push_mu_idx_, phi_g_, mu, nu);
        compute_psi_c_transform(phi_g_, push_nu_idx_, psi_f_, mu, nu);
        

        dual_value = calculate_dual_value(phi_g_, psi_f_, mu, nu);
        W2_value   = calculate_W2_value(push_mu_idx_);

        // sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu_h);

        return error_mu;
    }

// display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);

    void display_iteration(const int iter,const double W2_value,const double W2_value_back,const double dual_forth,const double dual_back,const double error_mu,const double error_nu) const{
        printf("iter: %5d    dual: %8.4f %8.4f    W2: %8.4f %8.4f\n", iter+1, dual_forth, dual_back, W2_value, W2_value_back);
    }

    void set_coeff(double& c1, double& c2, const double mu_max, const double C_c_transform){
        c1 = 0;
        c2 = C_c_transform;
    }

    bool check_collides(int* push_mu_idx_){
        for(int j=0;j<n_mu;++j){
            int j_idx = push_mu_idx_[j];

            for(int j1=j+1;j1<n_mu;++j1){
                int j1_idx = push_mu_idx_[j1];

                if(j_idx == j1_idx){
                    return false;
                }
            }
        }
        return true;
    }

    void start_OT(Points* mu, Points* nu, Plotting& plt){

        Initialize_C_mat_(mu, nu);

        const int skip = 1; // frequency of printout

        double error_mu = 1.0;
        double error_nu = 1.0;
        double error=1.0;

        /* Initialize coefficients for siga update */

        beta_1_ =0.1;
        beta_2_ =0.9;
        alpha_1_=1.05;
        alpha_2_=0.95;

        /* Calculate sup(mu) */

        cout << " Tolerance : " << tolerance_ << "\n";

        /* Initialize the constants */

        // double sigma_forth = sigma_;
        // double sigma_back  = sigma_;
        
        /* Starting the loop */
        
        for(int iter = 0; iter < max_iteration_; ++iter){


            /* Determinant version pushforward */
            
            error_mu = perform_OT_iteration_forth_det(sigma_forth_,W2_value_,dual_forth_,mu,nu,iter);
            error_nu = perform_OT_iteration_back_det (sigma_back_,W2_value_back_,dual_back_,mu,nu,iter);
            
            // error=fmax(error_mu,error_nu);

            /* Stopping Condition */

            if(check_collides(push_mu_idx_)){
                display_iteration(iter,W2_value_,W2_value_back_,dual_forth_,dual_back_,error_mu,error_nu);
                printf("Tolerance Met!\n");
                break;
            }

            /* Display the result per iterations */

            if(iter%skip==skip-1){
                display_iteration(iter,W2_value_,W2_value_back_,dual_forth_,dual_back_,error_mu,error_nu);
                cout << flush;

                // for(int p=0;p<nu->num_points_;++p){
                //     (*nu)(p,0) += 1e-2*(1.0*rand()/RAND_MAX-0.5);
                //     (*nu)(p,1) += 1e-2*(1.0*rand()/RAND_MAX-0.5);
                // }
                // Initialize_C_mat_(mu, nu);
                
            }
        }
    }

    void create_interpolate_video(Points* mu, Points* nu, int nt, Plotting& plt){
        string figurename = "video";

        plt.save_image_opencv(push_mu_idx_, mu, nu, figurename,0, nt);

        for(int n=1;n<nt+1;++n){
            plt.save_image_opencv(push_mu_idx_, mu, nu, figurename, n, nt);
        }
    }


}; // Back and Forth


#endif