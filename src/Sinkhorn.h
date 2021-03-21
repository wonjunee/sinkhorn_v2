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
#include "Plotting.h"

using namespace std;

const double LARGE_VALUE = 999999;

double dist2_points(double* p1, double* p2, int DIM){
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

    int* push_mu_idx_;

    double* C_mat_; // cost matrix

    double epsilon_;

    Sinkhorn(int DIM, int n_mu, int n_nu, int max_iteration, double tolerance, double lambda){

        epsilon_ = 1e-20;

        this->n_mu=n_mu;
        this->n_nu=n_nu;

        DIM_ = DIM;
        max_iteration_ =max_iteration;
        tolerance_     =tolerance;
        lambda_ = lambda;

        printf("n_mu : %d n_mu : %d DIM: %d\n", n_mu, n_nu, DIM_);

        phi_=new double[n_nu];
        psi_=new double[n_mu];

        for(int i=0;i<n_nu;++i) phi_[i] = 1;
        for(int j=0;j<n_mu;++j) psi_[j] = 1;

        C_mat_ = new double[n_mu * n_nu];

        push_mu_idx_ = new int[n_mu];

        // set a constant for the trace theorem
        C_tr1_ = 1;
        C_tr2_ = 1;

        phi_c1_ = 0;
        psi_c1_ = 0;

        phi_c2_ = 200;
        psi_c2_ = 200;

        vol_unit_ball_ = 1;

        cout << "lamdba : " << lambda_ << endl;

        // if(DIM_ > 2){
        //     int k = DIM_/2;
        //     if(DIM_ % 2 == 0){
        //         for(int i=1;i<k+1;++i){
        //             vol_unit_ball_ *= M_PI/i;
        //         }
        //     }else{
        //         vol_unit_ball_ = 2.0 / (2*k+1);
        //         for(int i=2*k-1;i>0;i-=2){
        //             vol_unit_ball_ *= (2.0*M_PI)/i;
        //             cout << "\nvol_unit_ball_: " << vol_unit_ball_ << endl;                    
        //         }
        //     }
        // }
    }

    ~Sinkhorn(){

        delete[] phi_;
        delete[] psi_;
        delete[] C_mat_;        
        delete[] push_mu_idx_;

        printf("bfm deconstruction done\n");
    }

    void Initialize_C_mat_(Points* mu, Points* nu){
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


        for(int i=0;i<n_nu;++i){
            for(int j=0;j<n_mu;++j){
                double dist2 = dist2_points(&(mu->data[j*DIM_]),&(nu->data[i*DIM_]),DIM_);
                C_mat_[i*n_mu+j] = 0.5 * dist2;
            }
        }
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

    double calculate_dual_value(const double* phi_, const double* psi_, Points* mu, Points* nu){
        double sum = 0;
        for(int i=0;i<n_nu;++i){
            sum += phi_[i];
        }
        for(int j=0;j<n_mu;++j){
            sum += psi_[j];
        }
        return sum/n_mu;
    }

    double calculate_W2_value(int* push_mu_idx_){
        double sum = 0;
        for(int j=0;j<n_mu;++j){
            sum += C_mat_[push_mu_idx_[j]*n_mu+j];
        }
        return sum;
    }

    double perform_sinkhorn_iteration(Points* mu, Points* nu, const int iter){
        // update psi
        for(int j=0;j<n_mu;++j){
            double eval = 0;
            for(int i=0;i<n_nu;++i) eval += exp(-1.0/lambda_ * (C_mat_[i*n_mu+j] - phi_[i]));
            psi_[j] = - lambda_ * log(eval);
        }
        // update phi
        for(int i=0;i<n_nu;++i){
            double eval = 0;
            for(int j=0;j<n_mu;++j) eval += exp(-1.0/lambda_ * (C_mat_[i*n_mu+j] - psi_[j] - phi_[i]));
            phi_[i] += lambda_ * ( - log(eval));
        }
        // calculate dual value
        double dual_value = 0;
        for(int i=0;i<n_nu;++i){
            for(int j=0;j<n_mu;++j){
                double eval = phi_[i] - C_mat_[i*n_mu+j] + psi_[j];
                dual_value += C_mat_[i*n_mu+j] * exp(eval/lambda_);
            }
        }
        return dual_value / n_mu;
    }

    void display_iteration(const int iter,const double dual_forth,const double rel_error) const{
        printf("iter: %5d  dual: %8.4f  rel error: %8.4e\n", iter+1, dual_forth, rel_error);
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

    void start_OT(Points* mu, Points* nu, Plotting& plt){

        Initialize_C_mat_(mu, nu);

        const int skip = 1; // frequency of printout
        cout << " Tolerance : " << tolerance_ << "\n";

        /* Initialize the constants */
        double previous_dual = 1;
        
        /* Starting the loop */
        for(int iter = 0; iter < max_iteration_; ++iter){
            /* Determinant version pushforward */
            dual_value_ = perform_sinkhorn_iteration(mu,nu,iter);

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
        }
    }

    void create_interpolate_video(Points* mu, Points* nu, int nt, Plotting& plt){
        string figurename = "video";

        plt.save_image_opencv(push_mu_idx_, mu, nu, figurename, 0, nt);

        for(int n=1;n<nt+1;++n){
            plt.save_image_opencv(push_mu_idx_, mu, nu, figurename, n, nt);
        }
    }


}; // Back and Forth


#endif