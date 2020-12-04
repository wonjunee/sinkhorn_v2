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
// #include "FLT_bf.h"
#include "FLT.h"
#include "PoissonSolver.h"
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
    double* push_mu_grid_;

    double* mu_grid_;
    double* nu_grid_;

    double* C_mat_; // cost matrix

    double* C_mu_;  // cost matrix for mu
    double* C_nu_;  // cost matrix for mu

    double epsilon_;

    poisson_solver* fftps_;
    FLT2D*          flt2d_;


    BackAndForth(int n1, int n2, int n_mu, int n_nu, int max_iteration, double tolerance, double sigma){

        epsilon_ = 1e-9;

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

        mu_grid_ = new double[n1*n2];
        nu_grid_ = new double[n1*n2];

        push_mu_grid_ = new double[n1*n2];

        flt2d_  = new FLT2D(n1,n2);

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

        printf("here\n");
        delete push_mu_;
        delete push_nu_;

        printf("there\n");

        delete[] mu_grid_;
        delete[] nu_grid_;

        delete flt2d_;
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
            for(int j=0;j<n_nu;++j){
                double mux = (*nu)(j,0);
                double muy = (*nu)(j,1);
                double nux = (*nu)(i,0);
                double nuy = (*nu)(i,1);
                double px  = mux - nux;
                double py  = muy - nuy;
                C_nu_[i*n_nu+j] = 0.5 * (px*px + py*py);
            }
        }

        for(int i=0;i<n_mu;++i){
            for(int j=0;j<n_mu;++j){
                double mux = (*mu)(j,0);
                double muy = (*mu)(j,1);
                double nux = (*mu)(i,0);
                double nuy = (*mu)(i,1);
                double px  = mux - nux;
                double py  = muy - nuy;
                C_mu_[i*n_mu+j] = 0.5 * (px*px + py*py);
            }
        }
    }

    void compute_psi_c_transform(double* phi_g_, Points* push_nu_, const double* psi_f_, Points* mu, Points* nu){
        for(int i=0;i<n_nu;++i){
            double eval = C_mat_[i*n_mu+0] - psi_f_[0];
            double min_val = eval;
            int    min_idx = 0;
            for(int j=1;j<n_mu;++j){
                eval = C_mat_[i*n_mu+j] - psi_f_[j];
                if(eval <= min_val){
                    min_val = eval;
                    min_idx = j;
                }
            }
            phi_g_[i] = min_val;
            (*push_nu_)(i,0) = (*mu)(min_idx,0);
            (*push_nu_)(i,1) = (*mu)(min_idx,1);
        }
    }

    void compute_phi_c_transform(double* psi_f_, Points* push_mu_, double* phi_g_, Points* mu, Points* nu){
        for(int j=0;j<n_mu;++j){
            double eval = C_mat_[0*n_mu+j] - phi_g_[0];
            double min_val = eval;
            int    min_idx = 0;
            for(int i=1;i<n_nu;++i){
                eval = C_mat_[i*n_mu+j] - phi_g_[i];
                if(eval <= min_val){
                    min_val = eval;
                    min_idx = i;
                }
            }
            psi_f_[j] = min_val;
            (*push_mu_)(j,0) = (*nu)(min_idx,0);
            (*push_mu_)(j,1) = (*nu)(min_idx,1);
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

    double perform_OT_iteration_back_det(double& sigma, double& W2_value, Points* mu, Points* nu, const int iter){
        double W2_value_previous = 0;
        double error_nu = 0;

        compute_phi_c_transform(psi_f_, push_mu_, phi_g_, mu, nu);
        compute_psi_c_transform(phi_g_, push_nu_, psi_f_, mu, nu);


        for(int j=0;j<n_mu;++j){

            double update_val = 0;
            // update_val += - sigma/(2.0*M_PI) * log(epsilon_);
            // update_val += - sigma/(2.0*M_PI) * (-LARGE_VALUE);
            double mux = (*mu)(j,0);
            double muy = (*mu)(j,1);

            for(int j1=0;j1<n_mu;++j1){
                double nux = (*mu)(j1,0);
                double nuy = (*mu)(j1,1);
                double px  = mux - nux;
                double py  = muy - nuy;

                double eval= sqrt(px*px + py*py);

                update_val += - sigma/(2.0*M_PI) * log(epsilon_ + eval);
                // if(eval == 0) update_val += LARGE_VALUE;
            }

            for(int i=0;i<n_nu;++i){
                double nux = (*push_nu_)(i,0);
                double nuy = (*push_nu_)(i,1);
                double px  = mux - nux;
                double py  = muy - nuy;

                double eval= sqrt(px*px + py*py);

                update_val += sigma/(2.0*M_PI) * log(epsilon_ + eval);
                // if(eval == 0) update_val -= LARGE_VALUE;

                // if(eval > 0) update_val += sigma/(2.0*M_PI) * log(eval);
                // else         update_val += sigma/(2.0*M_PI) * (-LARGE_VALUE);

            }
            psi_f_[j] += update_val;
        }

        compute_psi_c_transform(phi_g_, push_nu_, psi_f_, mu, nu);
        compute_phi_c_transform(psi_f_, push_mu_, phi_g_, mu, nu);

        // flt2d_->find_c_concave(psi_,phi_,1);
        // flt2d_->find_c_concave(phi_,psi_,1);

        // calculate_gradient(vx_, vy_, phi_);

        // // pushforward helper_f.DUstar_ -> push_mu_
        // calculate_push_rho(push_mu_grid_, push_nu_, nu, vx_, vy_, phi_);

        // fftps_->perform_inverse_laplacian(push_mu_grid_, mu_grid_, psi_c1_, psi_c2_, sigma);

        // W2_value_previous=calculate_dual_value(phi_,psi_,mu_grid_,nu_grid_);
        // double error_nu_h=calculate_h_minus_1(push_mu_grid_,mu_grid_);

        // error_nu=calculate_L1_error(push_mu_grid_,mu_grid_);

        // for(int i=0;i<n1*n2;++i) psi_[i] += fftps_->workspace[i];

        // flt2d_->find_c_concave(psi_,phi_,1);

        W2_value=calculate_dual_value(phi_g_, psi_f_, mu, nu);

        // sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu_h);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(double& sigma, double& W2_value, Points* mu, Points* nu, const int iter){
        double W2_value_previous = 0;
        double error_mu = 1;

        compute_psi_c_transform(phi_g_, push_nu_, psi_f_, mu, nu);
        compute_phi_c_transform(psi_f_, push_mu_, phi_g_, mu, nu);


        for(int i=0;i<n_nu;++i){

            double update_val = 0;
            // update_val += - sigma/(2.0*M_PI) * log(epsilon_);
            // phi_g_[i] += - sigma/(2.0*M_PI) * (-LARGE_VALUE);
            double nux = (*nu)(i,0);
            double nuy = (*nu)(i,1);

            for(int i1=0;i1<n_nu;++i1){
                double mux = (*nu)(i1,0);
                double muy = (*nu)(i1,1);
                double px  = mux - nux;
                double py  = muy - nuy;

                double eval= sqrt(px*px + py*py);

                update_val += - sigma/(2.0*M_PI) * log(epsilon_ + eval);
                // if(eval == 0) update_val += LARGE_VALUE;
            }

            for(int j=0;j<n_mu;++j){
                double mux = (*push_mu_)(j,0);
                double muy = (*push_mu_)(j,1);
                double px  = mux - nux;
                double py  = muy - nuy;

                double eval= sqrt(px*px + py*py);

                update_val += sigma/(2.0*M_PI) * log(epsilon_ + eval);
                // if(eval == 0) update_val -= LARGE_VALUE;

                // if(eval > 0){
                //     phi_g_[i] += sigma/(2.0*M_PI) * log(eval);
                // }else{
                //     phi_g_[i] += sigma/(2.0*M_PI) * (-LARGE_VALUE);
                // }
            }

            phi_g_[i] += update_val;
        }



        compute_phi_c_transform(psi_f_, push_mu_, phi_g_, mu, nu);
        compute_psi_c_transform(phi_g_, push_nu_, psi_f_, mu, nu);

        


        // flt2d_->find_c_concave(phi_,psi_,1);
        // flt2d_->find_c_concave(psi_,phi_,1);


        // calculate_gradient(vx_, vy_, psi_);

        // // pushforward mu -> push_mu_
        // calculate_push_rho(push_mu_grid_, push_mu_, mu, vx_, vy_, psi_);
            
        // fftps_->perform_inverse_laplacian(push_mu_grid_, nu_grid_, phi_c1_, phi_c2_, sigma);

        // W2_value_previous=calculate_dual_value(phi_,psi_,mu_grid_,nu_grid_);
        // double error_mu_h=calculate_h_minus_1(push_mu_grid_,nu_grid_);

        // error_mu=calculate_L1_error(push_mu_grid_, nu_grid_);

        // for(int i=0;i<n1*n2;++i) phi_[i] += fftps_->workspace[i];

        // flt2d_->find_c_concave(phi_,psi_,1);

        W2_value=calculate_dual_value(phi_g_, psi_f_, mu, nu);

        // sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu_h);

        return error_mu;
    }

// display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);

    void display_iteration(const int iter,const double W2_value,const double W2_value_back,const double error_mu,const double error_nu) const{
        printf("%*d c2: %*.2f %*.2f\tdual: %*f %*f\tL1 error: %*f %*f\n",5,iter+1, 8, phi_c2_, 8, psi_c2_, 8, W2_value, 8, W2_value_back, 8, error_mu, 8, error_nu);
    }

    void set_coeff(double& c1, double& c2, const double mu_max, const double C_c_transform){
        c1 = 0;
        c2 = C_c_transform;
    }

    bool check_collides(Points* push_mu_){
        for(int j=0;j<n_mu;++j){
            double x0 = (*push_mu_)(j,0);
            double y0 = (*push_mu_)(j,1);

            for(int j1=j+1;j1<n_mu;++j1){
                double x1 = (*push_mu_)(j1,0);
                double y1 = (*push_mu_)(j1,1);

                double eval = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);

                if(eval == 0){
                    return false;
                }
            }
        }
        return true;
    }

    void start_OT(Points* mu, Points* nu, Plotting& plt){

        push_mu_ = new Points(mu->DIM(), mu->num_points());
        push_nu_ = new Points(nu->DIM(), nu->num_points());

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

        double mu_max = 1;
        for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu_grid_[i]);

        cout << " Tolerance : " << tolerance_ << "\n";

        /* Initialize the constants */

        // double sigma_forth = sigma_;
        // double sigma_back  = sigma_;

        double C_c_transform = 1;

        {
            string figurename = "mu";
            plt.save_image_opencv(mu_grid_, figurename, 0);
        }

        {
            string figurename = "nu";
            plt.save_image_opencv(nu_grid_, figurename, 0);
        }
        
        /* Starting the loop */
        
        for(int iter = 0; iter < max_iteration_; ++iter){


            /* Determinant version pushforward */
            
            error_mu = perform_OT_iteration_forth_det(sigma_forth_,W2_value_,mu,nu,iter);

            // if(iter%skip==skip-1){
            //     string figurename = "output";
            //     plt.save_image_opencv(push_mu_ ,nu, figurename,(iter+1)/skip);
            // }

            error_nu = perform_OT_iteration_back_det (sigma_back_,W2_value_back_,mu,nu,iter);
            
            error=fmax(error_mu,error_nu);


            /* Stopping Condition */

            if(check_collides(push_mu_)){
                // cout<<"iteration: " << iter << " dual: " << W2_value_ << "" << " Tolerance met!"<<"\n";
                display_iteration(iter,W2_value_,W2_value_back_,error_mu,error_nu);
                printf("Tolerance Met!\n");
                break;
            }

            /* Display the result per iterations */

            if(iter%skip==skip-1){
                display_iteration(iter,W2_value_,W2_value_back_,error_mu,error_nu);
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

        plt.save_image_opencv(push_mu_, mu, nu, figurename,0, nt);

        for(int n=1;n<nt+1;++n){
            plt.save_image_opencv(push_mu_, mu, nu, figurename, n, nt);
        }
    }


}; // Back and Forth


#endif