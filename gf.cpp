#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include "Pushforward.h"
// #include "FLT.h"
#include "FLT_bf.h"
#include "PoissonSolver.h"
#include "Accessory.h"
#include "Barenblatt.h"
#include "Initializer.h"
#include "Helper_U.h"
using namespace std;



class BackAndForth{
public:
    int n1;
    int n2;

    int max_iteration_;
    double tolerance_;
    double W2_value_;

    double C; // Coefficients for trance theorem

    double phi_c1_;
    double phi_c2_;
    double psi_c1_;
    double psi_c2_;

    double beta_1_;
    double beta_2_;
    double alpha_1_;
    double alpha_2_;

    double gamma_;
    double tau_;
    double m_;
    double mprime_;

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

    poisson_solver* fftps_;
    FLT2D*          flt2d_;

    BackAndForth(int n1, int n2, int max_iteration, double tolerance, double gamma, double tau, double m, double C){

        this->n1=n1;
        this->n2=n2;
        max_iteration_ =max_iteration;
        tolerance_     =tolerance;

        gamma_  = gamma;
        tau_    = tau;
        m_      = m;
        mprime_ = m/(m-1);

        vxx_=new double[n1*n2];
        vyy_=new double[n1*n2];
        vxy_=new double[n1*n2];

        phi_=new double[n1*n2];
        psi_=new double[n1*n2];

        push_mu_= new double[n1*n2];

        flt2d_  = new FLT2D(n1,n2);

        // set a constant for the trace theorem
        C_phi_ = C;
        C_psi_ = C;

        clock_t time;
        time=clock();
        fftps_ = new poisson_solver(n1,n2);
        time=clock()-time;
        printf ("\nCPU time for FFT: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);
    }

    ~BackAndForth(){
        delete[] vxx_;
        delete[] vyy_;
        delete[] vxy_;
        delete[] phi_;
        delete[] psi_;
        delete[] push_mu_;

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

    double calculate_dual_value(Helper_U& helper_f, const double* phi, const double* psi, const double* mu){

         // Here psi is assumed to correspond to the c-transform of phi
        int pcount=n1*n2;
        
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term1 += helper_f.calculate_U(phi,i,j);
            }
        }

        term1/=pcount;

        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){                
                term2+=psi[i*n1+j]*mu[i*n1+j];
            }
        }
        
        term2/=pcount;
        
        return term2 - term1;
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

    double perform_OT_iteration_back_det(Helper_U& helper_f,double& sigma,double& W2_value,const double* mu, const int iter){
        // ------------------------------------------------------------
        double W2_value_previous = 0;
        double error_nu = 0;

        flt2d_->find_c_concave(psi_,phi_,tau_);
        flt2d_->find_c_concave(phi_,psi_,tau_);
    
        helper_f.calculate_DUstar_normalized(phi_);    

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, phi_);

        // pushforward helper_f.DUstar_ -> push_mu_
        calculate_push_pull_rho(push_mu_, helper_f.DUstar_, vxx_, vyy_, vxy_, psi_, 0.99);


        // calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, psi_);

        // // pushforward helper_f.DUstar_ -> push_mu_
        // calculate_pull_rho(push_mu_, helper_f.DUstar_, vxx_, vyy_, vxy_, psi_);


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


    double perform_OT_iteration_forth_det(Helper_U& helper_f,double& sigma,double& W2_value,const double* mu, const int iter){
        // ------------------------------------------------------------
        double W2_value_previous = 0;
        double error_mu = 0;

        flt2d_->find_c_concave(phi_,psi_,tau_);
        flt2d_->find_c_concave(psi_,phi_,tau_);

        helper_f.calculate_DUstar_normalized(phi_);    
            
        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, psi_);

        // pushforward mu -> push_mu_
        calculate_push_pull_rho(push_mu_, mu, vxx_, vyy_, vxy_, phi_, 0.2);


        // calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, psi_);

        // // pushforward mu -> push_mu_
        // calculate_push_rho(push_mu_, mu, vxx_, vyy_, vxy_, phi_);

        // calculate_gradient_vxx_vyy_vxy(phi, vxx_, vyy_, vxy_);
        // calculate_pull_rho(mu, push_mu,vx,vy,vxx_,vyy_,vxy_);    
            
        fftps_->perform_inverse_laplacian(push_mu_, helper_f.DUstar_, phi_c1_, phi_c2_, sigma);


        W2_value_previous=calculate_dual_value(helper_f,phi_,psi_,mu);
        double error_mu_h=calculate_h_minus_1(push_mu_,helper_f.DUstar_);

        error_mu=calculate_L1_error(push_mu_, helper_f.DUstar_);

        for(int i=0;i<n1*n2;++i){
            phi_[i] += fftps_->workspace[i];
        }


        flt2d_->find_c_concave(phi_,psi_,tau_);

        W2_value=calculate_dual_value(helper_f,phi_,psi_,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu_h);

        return error_mu;
    }

// display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);

    void display_iteration(const int iter,const double W2_value,const double error_mu,const double error_nu, const double C_phi, const double C_psi, const double infgradphi) const{
        cout << setprecision(6);
        cout << fixed;
        cout <<setw(5)<<iter+1 << " C : " << C_phi << "  infgradphi : " << infgradphi << "  c1 : " << phi_c1_ << " " << psi_c1_ << " coeff : " << setw(10) << phi_c2_/phi_c1_ << " " << setw(6) << psi_c2_/psi_c1_ << "  W2 : " << scientific << setw(13) << W2_value << "  L1 error : "<<scientific<<setw(13) << error_mu << " " << error_nu<<"\n";
    }

    /**
        Calculate a = 0.1 * max(-phi)
    */
    double calculate_lambda(Helper_U& helper_f) const{

        double phimax = 1;
        
        for(int i=0;i<n1*n2;++i){
            if(helper_f.obstacle_[i] == 0){
                phimax = fmax(phimax,-phi_[i]-helper_f.V_[i]);
            }
        }

        return fmin(0.1, phimax * 0.1);
        // return fmax(0.1, phimax * 0.05);
    }

    /**
        Calculate infgradphi = inf(|nabla phi|)
    */
    double calculate_infgradphi_on_level_set(const double lambda, const double* V, const unsigned char* obstacle){

        double infgradphi= 10000;
        int count = 0;
        for(int i=1;i<n2-1;++i){
            for(int j=1;j<n1-1;++j){
                if(obstacle[i*n1+j] == 0 && obstacle[i*n1+j+1] == 0 && obstacle[(i+1)*n1+j] == 0  && obstacle[i*n1+j-1] == 0 && obstacle[(i-1)*n1+j] == 0){
                    if(-phi_[i*n1+j]-V[i*n1+j] > 0 && -phi_[i*n1+j]-V[i*n1+j] < lambda){
                        double gradxphi = 0.5*n1*(phi_[i*n1+j+1]-phi_[i*n1+j-1]-V[i*n1+j+1]+V[i*n1+j-1]);
                        double gradyphi = 0.5*n2*(phi_[(i+1)*n1+j]-phi_[(i-1)*n1+j]-V[(i+1)*n1+j]+V[(i-1)*n1+j]);
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

        return fmax(5.0,sqrt(infgradphi));
    }

    void set_coeff_m_2(double& c1, double& c2, const double mu_max, const double C_c_transform){
        c1 = 1/gamma_;
        c2 = C_c_transform * tau_;
    }

    void set_coeff(double& c1, double& c2, const double C, const double mu_max, const double d1, const double d2, const double d3, const bool verbose, const double C_c_transform){
        c1 = C * d1 + d2;
        c2 = C * d1 + C_c_transform * tau_;
    }

    void initialize_phi(Helper_U& helper_f,const double* mu, const int outer_iter){

        for(int i=0;i<n1*n2;++i){
            if(mu[i] > 0) push_mu_[i] = 0;
            else          push_mu_[i] = -1.0/tau_*100;
        }

        flt2d_->find_c_concave(push_mu_, push_mu_, tau_);

        for(int i=0;i<n1*n2;++i){
            phi_[i] = - gamma_ * mprime_ * pow(mu[i],m_-1) + push_mu_[i] - helper_f.V_[i];
        }
    }

    void calculate_d1_d2_d3(double& d1, double& d2, double& d3, const double lambda, const double infgradphi, const double* nu){
        double eval = pow(gamma_ * mprime_, 1 - mprime_);

        d1 = eval * pow(lambda, mprime_-1) / infgradphi;
        d2 = eval * pow(lambda, mprime_-2) * (mprime_ - 1);
    }

    void calculate_c_transform_constant(double& C, const double* phi, const double* mu){

        C = 0;

        calculate_gradient_vxx_vyy_vxy(vxx_, vyy_, vxy_, phi);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=0.5*n1*(phi[i*n1+(int)fmin(n1-1,j+1)]-phi[i*n1+(int)fmax(0,j-1)]);
                double vyval=0.5*n2*(phi[(int)fmin(n2-1,i+1)*n1+j]-phi[(int)fmax(0,i-1)*n1+j]);

                double x = (j+0.5)/n1 - tau_ * vxval;
                double y = (i+0.5)/n2 - tau_ * vyval;

                double mu_val = interpolate_function(x,y,mu);

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

        // double mu_max = 0; for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu[i]); C *= mu_max;
    }

    void start_OT(Helper_U& helper_f, const double* mu, const int outer_iter, Initializer& init){

        const int skip = 50; // frequency of printout

        double error_mu = 1.0;
        double error_nu = 1.0;
        double error=1.0;

        /*
            Initialize coefficients for siga update
        */

        beta_1_ =0.1;
        beta_2_ =0.9;
        alpha_1_=1.05;
        alpha_2_=0.95;

        /*
            Initialize the tolerance based on tau^2
        */

        double mu_max = 1;
        for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu[i]);

        cout << "Iter : " << outer_iter + 1 << " Tolerance : " << tolerance_ << "\n";

        /*
            Initialize the coefficients for fftps
        */

        double lambda = 1;
        double infgradphi  = 1;
        double sigma_forth = 1;
        double sigma_back  = 1;

        double d1 = 1;
        double d2 = 1;
        double d3 = 1;

        double C_c_transform = 1;

        initialize_phi(helper_f,mu,outer_iter); // intiailize phi in the first outer iteration

        {
            string figurename = "output-phi";
            init.save_image_opencv(phi_, figurename,outer_iter, -1);
        }
        
        /*
            Starting the loop
        */
        
        for(int iter = 0; iter < max_iteration_; ++iter){
            
            /* 
                Calculating the relative error
            */

            if(iter % 50 == 0 && iter >= 0){
                if(m_ == 2){
                    calculate_c_transform_constant(C_c_transform, phi_, mu);
                    set_coeff_m_2(phi_c1_, phi_c2_, mu_max, C_c_transform);
                    calculate_c_transform_constant(C_c_transform, psi_, helper_f.DUstar_);
                    set_coeff_m_2(psi_c1_, psi_c2_, mu_max, C_c_transform);    
                }else{
                    lambda = calculate_lambda(helper_f);
                    infgradphi = calculate_infgradphi_on_level_set(lambda,helper_f.V_,helper_f.obstacle_);
                    calculate_d1_d2_d3(d1, d2, d3, lambda, infgradphi, helper_f.V_);

                    calculate_c_transform_constant(C_c_transform, phi_, mu);
                    set_coeff(phi_c1_, phi_c2_, C_phi_, mu_max, d1, d2, d3, false, C_c_transform);
                    // calculate_c_transform_constant(C_c_transform, psi_, helper_f.DUstar_);
                    set_coeff(psi_c1_, psi_c2_, C_psi_, mu_max, d1, d2, d3, false, C_c_transform);    
                }
            }

            /*
                Determinant version pushforward
            */

            sigma_forth = 1;
            sigma_back  = 1;
            
            error_mu = perform_OT_iteration_forth_det(helper_f,sigma_forth,W2_value_,mu,iter);
            error_nu = perform_OT_iteration_back_det (helper_f,sigma_back, W2_value_,mu,iter);
            
            error=fmin(error_mu,error_nu);

            /* 
                Stopping Condition 
            */

            if(((abs(error)<tolerance_ && abs(error)>0 && iter>=0) || iter==max_iteration_-1) ){
                // cout << "error : " << error << " iter : " << iter << "\n";
                // cout<<"Tolerance met!"<<"\n";
                display_iteration(iter,W2_value_,error_mu,error_nu,C_phi_,C_psi_,infgradphi);
                break;
            }

            /*
                Display the result per iterations
            */

            // cout<<"-"<<flush;

            if(iter%skip==skip-1){
                // cout<<"|";
                display_iteration(iter,W2_value_,error_mu,error_nu,C_phi_,C_psi_,infgradphi);

                string figurename = "output";
                for(int i=0;i<n1*n2;++i) push_mu_[i] = fabs(push_mu_[i] - mu[i]);
                init.save_image_opencv(push_mu_, figurename,(iter+1)/skip, mu_max);
                // init.save_image_opencv(helper_f.DUstar_,figurename,(iter+1)/skip, mu_max);
            }
        }
    }


}; // Back and Forth

int main(int argc, char** argv){
    if(argc!=11){
        cout<<"Do the following:"<<"\n";
        cout<<"./gf [n1] [n2] [max_iteration] [tolerance] [nt] [tau] [gamma] [m] [C trace] [opencv 0 or 1]"<<"\n";
        return 0;
    }

    int n1=stoi(argv[1]);
    int n2=stoi(argv[2]);
    int max_iteration=stoi(argv[3]);
    double tolerance=stod(argv[4]);
    int nt=stoi(argv[5]);
    double tau=stod(argv[6]);
    double gamma=stod(argv[7]);
    double m=stod(argv[8]);
    double C=stod(argv[9]);
    int plot=stoi(argv[10]);

    cout << "n1   : " << n1 <<"\n";
    cout << "n2   : " << n2 <<"\n";
    cout << "nt   : " << nt <<"\n";
    cout << "tau  : " << tau << "\n";
    cout << "gamma: " << gamma << "\n";
    cout << "m    : " << m << "\n";
    cout << "Max Iteration : " << max_iteration <<"\n";
    cout << "tolerance     : " << tolerance <<"\n";

    /* Initialize Initializer */
    Initializer init(n1,n2);

    create_csv_parameters(n1,n2,nt,tau,gamma,m);

    // Initialize mu and obstacle
    double* mu=new double[n1*n2];
    unsigned char* obstacle = new unsigned char[n1*n2];

    create_mu_square(mu,0.2,0.2,0.1,n1,n2);
    // create_mu_from_image(mu,n1,n2);
    // create_mu_from_image2(mu,n1,n2);

    for(int i=0;i<n1*n2;++i) mu[i] *= 0.1;

    cout << "XXX Starting Gradient Flow XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);
  
    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));
    cout << "\n";

    

    string data_folder = "data";

    Helper_U helper_f(n1,n2,gamma,tau,m,mu);

    double a = 5;
    // init_entropy_sine(helper_f.V_, 5, 3, a, n1, n2);
    init_entropy_quadratic(helper_f.V_, 0.5, 0.5, a, n1, n2);

    // init_obstacle_from_image(obstacle, n1, n2);
    // init_obstacle_two_moons(obstacle, n1, n2);
    // init_obstacle_pac_man(obstacle, n1, n2);
    init_obstacle_circle(obstacle, n1, n2);

    helper_f.set_obstacle(obstacle);

    cout << setprecision(6);

    BackAndForth bf(n1,n2,max_iteration,tolerance,gamma,tau,m,C);

    string filename="./data/mu-"+to_string(0)+".csv";
    create_bin_file(mu,n1*n2,filename);

    string figurename = "barenblatt";

    if(plot > 0) init.save_image_opencv(mu,obstacle,figurename,0);
    // if(plot > 0) init.save_image_opencv(mu,figurename,0);

    clock_t time;
    time=clock();

    double sum = 0;

    for(int n=0;n<nt;++n){
        clock_t time_outer = clock();
        bf.start_OT(helper_f, mu, n, init);
        helper_f.calculate_DUstar_normalized(bf.phi_);
        memcpy(mu,helper_f.DUstar_,n1*n2*sizeof(double));

        /* print out sum */
        sum=0; for(int i=0;i<n1*n2;++i) sum+=mu[i]; cout << "sum : " << sum/(n1*n2) << "\n";

        filename="./data/mu-"+to_string(n+1)+".csv";
        create_bin_file(mu,n1*n2,filename);
        time_outer=clock()-time_outer;
        printf ("\nCPU time for outer iteration: %f seconds.\n\n",((float)time_outer)/CLOCKS_PER_SEC);
        if(plot > 0) init.save_image_opencv(mu,obstacle,figurename,n+1);
        // if(plot > 0) init.save_image_opencv(mu,figurename,n+1);

    }

    time=clock()-time;
    printf ("\nCPU time for GF: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    delete[] mu;
}