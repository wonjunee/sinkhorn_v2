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
#include "FLT_bf.h"
#include "PoissonSolver.h"
#include "Accessory.h"
#include "Plotting.h"

using namespace std;

const double LARGE_VALUE = 999999;

class BackAndForth{
public:
    int n1;
    int n2;
    int DIM_;

    int max_iteration_;
    double tolerance_;
    double W2_value_;

    double smooth_; // Coefficients for trace theorem

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

    double* phi_;
    double* psi_;

    Points*  push_mu_;
    Points*  push_nu_;
    double* push_mu_grid_;

    double* mu_grid_;
    double* nu_grid_;

    double* vxx_;
    double* vyy_;
    double* vxy_;

    double* vx_;
    double* vy_;

    poisson_solver* fftps_;
    FLT2D*          flt2d_;


    BackAndForth(int n1, int n2, int max_iteration, double tolerance, double smooth){

        this->n1=n1;
        this->n2=n2;
        DIM_ = 2;
        max_iteration_ =max_iteration;
        tolerance_     =tolerance;
        smooth_ = smooth;

        vx_=new double[n1*n2];
        vy_=new double[n1*n2];

        vxx_=new double[n1*n2];
        vyy_=new double[n1*n2];
        vxy_=new double[n1*n2];

        phi_=new double[n1*n2];
        psi_=new double[n1*n2];

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
        delete[] vx_;
        delete[] vy_;
        delete[] vxx_;
        delete[] vyy_;
        delete[] vxy_;
        delete[] phi_;
        delete[] psi_;
        delete[] push_mu_grid_;

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

    void calculate_gradient(double* vx, double* vy, const double* phi_c){
        /* centered difference*/
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                vx[i*n1+j]=0.5*n1*(phi_c[i*n1+(int)fmin(n1-1,j+1)]-phi_c[i*n1+(int)fmax(0,j-1)]);
                vy[i*n1+j]=0.5*n2*(phi_c[(int)fmin(n2-1,i+1)*n1+j]-phi_c[(int)fmax(0,i-1)*n1+j]);
            }
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

    void calculate_push_rho(double* push_rho_grid, Points* push_rho, Points* rho, const double* vx, const double* vy, const double* phi_c, const double tau=1){

        for(int p=0;p<rho->num_points();++p){
            double px = (*rho)(p,0);
            double py = (*rho)(p,1);

            double dx_phi_c = interpolate_function(px,py,vx);
            double dy_phi_c = interpolate_function(px,py,vy);

            (*push_rho)(p,0) = px - tau * dx_phi_c;
            (*push_rho)(p,1) = py - tau * dy_phi_c;
        }

        setup_grid(push_rho_grid, push_rho);
    }

    void calculate_pull_rho(double* push_rho, const double* rho, const double* vxx, const double* vyy, const double* vxy, const double* phi){

        double eps = pow(1.0/n1, 0.7);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=0.5*n1*(phi[i*n1+(int)fmin(n1-1,j+1)]-phi[i*n1+(int)fmax(0,j-1)]);
                double vyval=0.5*n2*(phi[(int)fmin(n2-1,i+1)*n1+j]-phi[(int)fmax(0,i-1)*n1+j]);

                double x=(j+0.5)/(1.0*n1)-vxval;
                double y=(i+0.5)/(1.0*n2)-vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){
                    double vxx_val = vxx[i*n1+j];
                    double vyy_val = vyy[i*n1+j];
                    double vxy_val = vxy[i*n1+j];

                    double det = fabs((1.0-vxx_val) * (1.0-vyy_val) - vxy_val * vxy_val);
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

                double x=(j+0.5)/(1.0*n1)-vxval;
                double y=(i+0.5)/(1.0*n2)-vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){

                    double vxx_val = interpolate_function(x,y,vxx);
                    double vyy_val = interpolate_function(x,y,vyy);
                    double vxy_val = interpolate_function(x,y,vxy);

                    double det = fabs((1.0-vxx_val) * (1.0-vyy_val) - vxy_val * vxy_val);

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

                        det = fabs((1.0-vxx_val) * (1.0-vyy_val) - vxy_val * vxy_val);
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

    double perform_OT_iteration_back_det(double& sigma,double& W2_value,Points* mu, Points* nu, const int iter){
        // ------------------------------------------------------------
        double W2_value_previous = 0;
        double error_nu = 0;

        flt2d_->find_c_concave(psi_,phi_,1);
        flt2d_->find_c_concave(phi_,psi_,1);

        calculate_gradient(vx_, vy_, phi_);

        // pushforward helper_f.DUstar_ -> push_mu_
        calculate_push_rho(push_mu_grid_, push_nu_, nu, vx_, vy_, phi_);

        fftps_->perform_inverse_laplacian(push_mu_grid_, mu_grid_, psi_c1_, psi_c2_, sigma);

        W2_value_previous=calculate_dual_value(phi_,psi_,mu_grid_,nu_grid_);
        double error_nu_h=calculate_h_minus_1(push_mu_grid_,mu_grid_);

        error_nu=calculate_L1_error(push_mu_grid_,mu_grid_);

        for(int i=0;i<n1*n2;++i) psi_[i] += fftps_->workspace[i];

        flt2d_->find_c_concave(psi_,phi_,1);

        W2_value=calculate_dual_value(phi_,psi_,mu_grid_,nu_grid_);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu_h);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(double& sigma, double& W2_value, Points* mu, Points* nu, const int iter){
        // ------------------------------------------------------------
        double W2_value_previous = 0;
        double error_mu = 0;


        flt2d_->find_c_concave(phi_,psi_,1);
        flt2d_->find_c_concave(psi_,phi_,1);


        calculate_gradient(vx_, vy_, psi_);

        // pushforward mu -> push_mu_
        calculate_push_rho(push_mu_grid_, push_mu_, mu, vx_, vy_, psi_);
            
        fftps_->perform_inverse_laplacian(push_mu_grid_, nu_grid_, phi_c1_, phi_c2_, sigma);

        W2_value_previous=calculate_dual_value(phi_,psi_,mu_grid_,nu_grid_);
        double error_mu_h=calculate_h_minus_1(push_mu_grid_,nu_grid_);

        error_mu=calculate_L1_error(push_mu_grid_, nu_grid_);

        for(int i=0;i<n1*n2;++i) phi_[i] += fftps_->workspace[i];

        flt2d_->find_c_concave(phi_,psi_,1);

        W2_value=calculate_dual_value(phi_,psi_,mu_grid_,nu_grid_);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu_h);

        return error_mu;
    }

// display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);

    void display_iteration(const int iter,const double W2_value,const double error_mu,const double error_nu) const{
        printf("%*d c2: %*.2f %*.2f\tdual: %*f\tL1 error: %*f %*f\n",5,iter+1, 8, phi_c2_, 8, psi_c2_, 8, W2_value, 8, error_mu, 8, error_nu);
    }

    void set_coeff(double& c1, double& c2, const double mu_max, const double C_c_transform){
        c1 = 0;
        c2 = C_c_transform;
    }

    void calculate_c_transform_constant(double& C, const double* phi, const double* mu){

        C = 0;

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, phi);

        for(int i=1;i<n2-1;++i){
            for(int j=1;j<n1-1;++j){

                double vxval=0.5*n1*(phi[i*n1+(int)fmin(n1-1,j+1)]-phi[i*n1+(int)fmax(0,j-1)]);
                double vyval=0.5*n2*(phi[(int)fmin(n2-1,i+1)*n1+j]-phi[(int)fmax(0,i-1)*n1+j]);

                double x = (j+0.5)/n1 - vxval;
                double y = (i+0.5)/n2 - vyval;

                double mu_val = interpolate_function(x,y,mu);
                // double mu_val = 1;

                if(mu_val > 0){
                    /* calculate eigen values */
                        /* calculate trace */
                    double trace = fabs(2 - vxx_[i*n1+j] - vyy_[i*n1+j]);
                        /* calculate det */
                    double det   = fabs((1.0-vxx_[i*n1+j]) * (1.0-vyy_[i*n1+j]) - vxy_[i*n1+j] * vxy_[i*n1+j]);

                    double t1 = 0.5 * fabs(trace + sqrt(fabs(trace*trace - 4 * det)));
                    double t2 = 0.5 * fabs(trace - sqrt(fabs(trace*trace - 4 * det)));
                    
                    C = fmax(C, mu_val*fmax(t1,t2));
                }
            }
        }
        C = fmax(C, 1);
        // double mu_max = 0; for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu[i]); C *= mu_max;
    }

    void setup_grid(double* mu_grid, Points* mu){
        memset(mu_grid, 0, n1*n2*sizeof(double));

        for(int p=0;p<mu->num_points();++p){
            double px = (*mu)(p,0);
            double py = (*mu)(p,1);
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double x = (j+.5)/n1;
                    double y = (i+.5)/n2;
                    mu_grid[i*n1+j] += 1.0/(smooth_*smooth_*2*M_PI) * exp(-1.0/(2.0*smooth_*smooth_) * (pow(x-px,2) + pow(y-py,2)));
                }
            }
        }
    }

    void start_OT(Points* mu, Points* nu, Plotting& plt){

        push_mu_ = new Points(mu->DIM(), mu->num_points());
        push_nu_ = new Points(nu->DIM(), nu->num_points());

        setup_grid(mu_grid_,mu);
        setup_grid(nu_grid_,nu);

        const int skip = 10; // frequency of printout

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

        double sigma_forth = 1;
        double sigma_back  = 1;

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

            sigma_forth = 1;
            sigma_back  = 1;


            if((iter) % 10 == 0 && iter >= 0){
                calculate_c_transform_constant(C_c_transform, phi_, mu_grid_);
                set_coeff(phi_c1_, phi_c2_, mu_max, C_c_transform);
            }

            
            error_mu = perform_OT_iteration_forth_det(sigma_forth,W2_value_,mu,nu,iter);

            if(iter%skip==skip-1){
                string figurename = "output";
                calculate_gradient(vx_, vy_, psi_);
                calculate_push_rho(push_mu_grid_, push_mu_, mu, vx_, vy_, psi_, 1);
                plt.save_image_opencv(push_mu_,nu, figurename,(iter+1)/skip);
            }

            // string figurename = "output";
            // calculate_gradient(vx_, vy_, psi_);
            // calculate_push_rho(push_mu_grid_, push_mu_, mu, vx_, vy_, psi_);
            // plt.save_image_opencv(push_mu_, nu, figurename,(iter+1)/1);

            // {
            //     string figurename = "push_mu_";
            //     plt.save_image_opencv(push_mu_grid_, figurename, iter);
            // }

            if((iter) % 10 == 0 && iter >= 0){
                calculate_c_transform_constant(C_c_transform, phi_, nu_grid_);
                set_coeff(psi_c1_, psi_c2_, mu_max, C_c_transform);
            }



            error_nu = perform_OT_iteration_back_det (sigma_back,W2_value_,mu,nu,iter);
            
            error=fmax(error_mu,error_nu);

            /* Stopping Condition */

            if(((abs(error)<tolerance_ && abs(error)>0 && iter>=50) || iter==max_iteration_-1) ){
                cout<<"Tolerance met!"<<"\n";
                display_iteration(iter,W2_value_,error_mu,error_nu);
                break;
            }

            /* Display the result per iterations */

            cout<<"-"<<flush;
            if(iter%skip==skip-1){
                cout<<"|"; display_iteration(iter,W2_value_,error_mu,error_nu);
                cout << flush;
            }
        }

        delete push_mu_;
        delete push_nu_;
        push_mu_ = NULL;
        push_nu_ = NULL;
    }

    void create_interpolate_video(Points* mu, Points* nu, int nt, Plotting& plt){
        push_mu_ = new Points(mu->DIM(), mu->num_points());
        push_nu_ = new Points(nu->DIM(), nu->num_points());

        string figurename = "video";

        plt.save_image_opencv(mu, nu, figurename,0);

        for(int n=1;n<nt;++n){
            calculate_gradient(vx_, vy_, psi_);
            calculate_push_rho(push_mu_grid_, push_mu_, mu, vx_, vy_, psi_, 1.0*(n+1)/nt);
            plt.save_image_opencv(push_mu_, nu, figurename, n);
        }
    }


}; // Back and Forth


#endif