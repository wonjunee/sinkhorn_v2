#ifndef Sinkhorn_H
#define Sinkhorn_H

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include "Accessory.h"

using namespace std;

const double LARGE_VALUE = 999999;

double dist2_points(const double* p1, const double* p2, int DIM){
    double dist = 0;

    for(int p=0;p<DIM;++p){
        double diff = p1[p] - p2[p];
        dist += diff*diff;
    }

    return dist;
}

class Sinkhorn{
public:
    int n_mu;
    int n_nu;

    int DIM_;

    double vol_unit_ball_;

    int max_iteration_;
    double tolerance_;
    double W2_value_;
    double W2_value_back_;

    double dual_value_;

    double sigma_; // Coefficients for trace theorem
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

    double lambda_;

    double* phi_;
    double* psi_;
    double* rho_;

    int* push_mu_idx_;

    // double* C_mat_; // cost matrix
    // double* G_mat_; // cost matrix

    double* dual_value_list_;

    double epsilon_;

    Sinkhorn(const int DIM, const int n_mu, const int n_nu, const int max_iteration, const double tolerance, const double lambda, const double sigma){

        epsilon_ = 1e-20;

        this->n_mu=n_mu;
        this->n_nu=n_nu;

        DIM_ = DIM;
        max_iteration_ = max_iteration;
        tolerance_     = tolerance;
        lambda_        = lambda;
        sigma_         = sigma;

        printf("n_mu : %d n_mu : %d DIM: %d\n", n_mu, n_nu, DIM_);

        phi_=new double[n_nu];
        psi_=new double[n_mu];
        rho_=new double[n_mu]; // if n_mu != n_nu then I need to create two different rho_; or make rho_[fmax(n_mu,n_mu)]

        dual_value_list_=new double[max_iteration_];

        for(int i=0;i<n_nu;++i) phi_[i] = 0;
        for(int j=0;j<n_mu;++j) psi_[j] = 0;

        // C_mat_ = new double[n_mu * n_nu];
        // G_mat_ = new double[n_nu * n_nu];

        push_mu_idx_ = new int[n_mu];

        // set a constant for the trace theorem
        C_tr1_ = 1;
        C_tr2_ = 1;

        phi_c1_ = 0;
        psi_c1_ = 0;

        phi_c2_ = 200;
        psi_c2_ = 200;

        vol_unit_ball_ = 1;

        cout << "lambda: " << lambda_ << " sigma: " << sigma_ << endl;
    }

    ~Sinkhorn(){

        delete[] phi_;
        delete[] psi_;
        delete[] rho_;
        // delete[] C_mat_;        
        // delete[] G_mat_;
        delete[] push_mu_idx_;
        delete[] dual_value_list_;

        printf("bfm deconstruction done\n");
    }

    double calculate_log(double eval){
        return - 1.0/(4*M_PI) * log(1e-12 + eval);
    }

    double calculate_log_cutoff(double eval){
        double cutoff = 10000;
        if(eval == 0) return cutoff;
        return fmin(cutoff, - 1.0/(4*M_PI) * log(eval));
    }

    double no_precondition(double eval){
        if(eval == 0) return 1;
        return 0;
    }

    double calculate_guassian(double eval){
        double sigma = 0.5;
        return exp(-eval/(sigma*sigma));
    }

    // eval = |x - y|^2
    double G_(double eval){

        // return no_precondition(eval);
        // return calcaulte_power(eval);
        return calculate_log(eval);
        // return calculate_guassian(eval);
    }

        /** 
         *       Create a cost function 
         *
         * This function is used only once in the beginning.
         *
         * C_mat_ is a matrix C (num_points_nu x num_points_mu)
         * C_{ij} is a distance between mu_j and nu_i
         *
         * C_mu_ is a symmetric (n_mu x n_mu) matrix.
         * C_mu_{ij} is a distance between mu_i and mu_j
         *
         * C_nu_ is a symmetric (n_nu x n_nu) matrix.
         * C_nu_{ij} is a distance between nu_i and nu_j
         *
         * C_mu_sum_ is a (n_mu) vector, each index is a sum of row of C_mu_.
         * C_mu_sum_i = \sum_{j=1}^{n_mu} (-\Delta)^{-1} \delta_{x_j} (x_i)
         *
         * C_nu_sum_ is a (n_nu) vector, each index is a sum of row of C_nu_.
         * C_nu_sum_i = \sum_{j=1}^{n_nu} (-\Delta)^{-1} \delta_{y_j} (y_i)
         *
        **/
/*
    void Initialize_C_mat_(const double* pos_mu, const double* pos_nu){
        for(int i=0;i<n_nu;++i){
            for(int j=0;j<n_mu;++j){
                // double dist2 = dist2_points(&(mu->data[j*DIM_]),&(nu->data[i*DIM_]),DIM_);
                double dist2 = dist2_points(&pos_mu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                C_mat_[i*n_mu+j] = 0.5 * dist2;
            }
        }


        for(int i=0;i<n_nu;++i){
            for(int j=i;j<n_nu;++j){
                double dist2 = dist2_points(&pos_nu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                G_mat_[i*n_nu+j] = G_(dist2);
                G_mat_[j*n_nu+i] = G_mat_[i*n_nu+j];
            }
        }

    }
*/


    /**
        Update sigma based on Goldstein scheme
    */
    double update_sigma(const double sigma, const double W2_value, const double W2_value_previous, const double error){

        if(W2_value_previous-W2_value>-sigma*beta_1_*error){
            return sigma*alpha_2_;
        }else if(W2_value_previous-W2_value<-sigma*beta_2_*error){
            return sigma*alpha_1_;
        }
        return sigma;
    }

    double calculate_dual_value(const double* phi_, const double* psi_, const double* mu, const double* nu) const{
        double sum = 0;
        for(int i=0;i<n_nu;++i){
            sum += phi_[i] * nu[i];
        }
        for(int j=0;j<n_mu;++j){
            sum -= psi_[j] * mu[j];
        }
        return sum/n_mu;
    }

    double calc_L(double eval){
        return - log(1e-12+eval); // log (original)
        // double sig = 0.01; return exp(-(eval*eval)/(sig*sig)); // log (original)
        // return 1.0/(fabs(eval));
    }
    // modify psi_
    void compute_psi(const double* pos_mu, const double* pos_nu){
        double eval = 0;

        // compute psi
        for(int j=0;j<n_mu;++j){
            // eval    = phi_[0] - C_mat_[0*n_mu+j];
            double dist2 = dist2_points(&pos_mu[j*DIM_],&pos_nu[0*DIM_],DIM_);
            eval    = phi_[0] - dist2/2.0;
            for(int i=1;i<n_nu;++i){ 
                // eval = fmax(eval, phi_[i] - C_mat_[i*n_mu+j]); 
                double dist2 = dist2_points(&pos_mu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                eval = fmax(eval, phi_[i] - dist2/2.0); 
            }
            psi_[j] = eval;
        }
        // update psi
        for(int j=0;j<n_mu;++j){
            eval = 0;
            for(int i=0;i<n_nu;++i){
                double dist2 = dist2_points(&pos_mu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                // eval += exp( - (C_mat_[i*n_mu+j] + psi_[j] - phi_[i])/ lambda_ );
                eval += exp( - (dist2/2.0 + psi_[j] - phi_[i])/ lambda_ );
            }
            psi_[j] = lambda_ * log(eval) + psi_[j];
        }
    }
    // modify phi_
    void compute_phi_sinkhorn(const double* pos_mu, const double* pos_nu){
        double eval = 0;
        // update phi
        for(int i=0;i<n_nu;++i){
            eval = 0;
            for(int j=0;j<n_mu;++j){
                double dist2 = dist2_points(&pos_mu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                // eval += exp( - (C_mat_[i*n_mu+j] + psi_[j] - phi_[i])/ lambda_ );
                eval += exp( - (dist2/2.0 + psi_[j] - phi_[i])/ lambda_ );

            } 
            phi_[i] -= sigma_ * (calc_L(1) - calc_L(eval));
        }
    }
    // modify rho_
    void compute_rho(const double* pos_mu, const double* pos_nu){
        // define rho
        double eval = 0;
        for(int i=0;i<n_nu;++i){
            eval = 0;
            for(int j=0;j<n_mu;++j) {
                double dist2 = dist2_points(&pos_mu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                // eval += exp(- (C_mat_[i*n_mu+j] + psi_[j] - phi_[i])/ lambda_ );
                eval += exp(- (dist2/2.0 + psi_[j] - phi_[i])/ lambda_ );
            }
            rho_[i] = eval;
        }
    }
    // modify phi_ and return error
    double compute_phi_inverse_lap(const double* pos_mu, const double* pos_nu){
        double error = 0;
        double eval  = 0;
        // inverse laplacian
        for(int i=0;i<n_nu;++i){
            eval = 0;
            for(int j=0;j<n_nu;++j) {
                double dist2 = dist2_points(&pos_nu[j*DIM_],&pos_nu[i*DIM_],DIM_);
                // eval += (1-rho_[j]) * G_mat_[i*n_nu+j];
                eval += (1-rho_[j]) * G_(dist2);
            }
            phi_[i] += sigma_ * eval;
            error += eval;
        }
        return error / n_nu;
    }
    // original
    double perform_sinkhorn_iteration(const double* mu, const double* pos_mu, const double* nu, const double* pos_nu, const int iter){
        compute_psi(pos_mu, pos_nu);

        /* --- uncomment this to run sinkhorn --- */
        
        // compute_phi_sinkhorn(pos_mu,pos_nu); // sinkhorn
        // dual_value_ = calculate_dual_value(phi_, psi_, mu, nu);

        /* --- uncomment this to run laplacian version --- */

        double dual_value_previous = dual_value_; // update the previous dual value
        compute_rho(pos_mu,pos_nu); 
        double error = compute_phi_inverse_lap(pos_mu,pos_nu); // laplacian
        dual_value_ = calculate_dual_value(phi_, psi_, mu, nu);
        sigma_ = update_sigma(sigma_, dual_value_, dual_value_previous, error);

        /* --- this is the end of the laplacian version --- */

        return dual_value_;
    }

    void display_iteration(const int iter,const double dual_forth,const double rel_error) const{
        printf("iter: %5d dual: %8.4f rel error: %8.4e sigma: %8.4e\n", iter+1, dual_forth, rel_error, sigma_);
    }

    bool check_collides(int* push_mu_idx_){
        for(int j=0;j<n_mu;++j){
            int j_idx = push_mu_idx_[j];

            for(int j1=j+1;j1<n_mu;++j1){
                int j1_idx = push_mu_idx_[j1];

                if(j_idx == j1_idx) return false;
            }
        }
        return true;
    }

    void start_OT(const double* mu, const double* pos_mu, const double* nu, const double* pos_nu){

        beta_1_ =0.21;
        beta_2_ =0.8;
        alpha_1_=1.01;
        alpha_2_=0.9;

        // Initialize_C_mat_(pos_mu, pos_nu);

        const int skip = 100; // frequency of printout
        cout << " Tolerance : " << tolerance_ << "\n";

        /* Initialize the constants */
        double previous_dual = 1;
        
        /* Starting the loop */
        for(int iter = 0; iter < max_iteration_; ++iter){
            /* Determinant version pushforward */
            dual_value_ = perform_sinkhorn_iteration(mu,pos_mu,nu,pos_nu,iter);

            double rel_error = fabs((dual_value_ - previous_dual) / previous_dual);
            previous_dual =  dual_value_;

            /* Stopping Condition */
            if(rel_error < tolerance_){
                display_iteration(iter,dual_value_,rel_error);
                printf("Tolerance Met!\n");
                break;
            }

            /* Display the result per iterations */
            if(iter%skip==skip-1){
                display_iteration(iter,dual_value_,rel_error);
                cout << flush;
            }

            dual_value_list_[iter] = dual_value_;

            sigma_ = fmin(0.01,fmax(1e-5, sigma_));
        }
    }

}; // Back and Forth


#endif