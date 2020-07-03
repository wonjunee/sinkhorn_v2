#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include "Pushforward.h"
#include "FLT.h"
#include "PoissonSolver.h"
#include "Accessory.h"
#include "Barenblatt.h"

using namespace std;

class Helper_E{
public:
    int n1;
    int n2;

    double tau;

    double original_sum;

    double* nu;
    double* DEstar;

    Helper_E(){
        nu=NULL;
        DEstar=NULL;
    }
    Helper_E(int n1,int n2,double tau){
        initialize(n1,n2,tau);
    }

    void initialize(int n1,int n2,double tau){
        this->n1=n1;
        this->n2=n2;
        this->tau=tau;

        DEstar=new double[n1*n2];
        nu    =new double[n1*n2];
    }

    ~Helper_E(){
        delete[] DEstar;
        delete[] nu;
    }

    double calculate_E(const double* phi, const int i, const int j) const
    {
        double c_phi = - phi[i*n1+j] - nu[i*n1+j];
        if(nu[i*n1+j] >= 0){
            if(c_phi > 0 ){
                return c_phi;
            }    
        }
        return 0;
    }


    void setup_original_sum(const double* mu){
        double s=0;
        for(int i=0;i<n1*n2;++i){
            s += mu[i];
        }
        original_sum = s/(1.0*n1*n2);
    }

    void calculate_DEstar(const double* phi){

        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - nu[i];
            if(nu[i] >= 0){
                if(eval>0){
                    DEstar[i] = 1;    
                }else{
                    DEstar[i] = 0;
                }    
            }else{
                DEstar[i] = 0;
            }
        }
    }

    void calculate_DEstar(const double* phi,const double* mu){

        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - nu[i];
            if(nu[i] >= 0){
                if(eval>0){
                    DEstar[i] = 1;    
                }else{
                    DEstar[i] = fmin(1,mu[i]);
                }
            }else{
                DEstar[i] = 0;
            }
        }
    }

    void calculate_DEstar_normalized(const double* phi, const double tolerance=1e-10){
        
        double lambda_a=-phi[0]-nu[0];
        double lambda_b=-phi[0]-nu[0];

        for(int i=0;i<n1*n2;++i){
            lambda_a = fmax(lambda_a, - phi[i] - nu[i]);
            lambda_b = fmin(lambda_b, - phi[i] - nu[i]);
        }

        lambda_a = - lambda_a - 10;
        lambda_b = - lambda_b + 10;

        int max_iteration=1000;

        double lambda = 0.5 * (lambda_a + lambda_b);

        for(int iter=0;iter<max_iteration;++iter){
            double sum=0;

            for(int i=0;i<n1*n2;++i){
                double eval =- phi[i] - nu[i] + lambda;
                if(nu[i] >= 0){
                    if(eval>0){
                        sum += 1.0; // IMPORTANT    
                    }    
                }
            }

            double val =sum /(1.0*n1*n2) - original_sum;

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

        

        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - nu[i] + lambda;
            if(nu[i]>=0){
                if(eval>0){
                    DEstar[i] = 1.0;
                }else{
                    DEstar[i] = 0.0;
                }
            }else{
                DEstar[i] = 0;
            }
            
        }
    }


    void calculate_DEstar_normalized(const double* phi, const double* mu, const double tolerance=1e-7){
        
        double lambda_a=-phi[0]-nu[0];
        double lambda_b=-phi[0]-nu[0];

        for(int i=0;i<n1*n2;++i){
            lambda_a = fmax(lambda_a, - phi[i] - nu[i]);
            lambda_b = fmin(lambda_b, - phi[i] - nu[i]);
        }

        lambda_a = - lambda_a - 10;
        lambda_b = - lambda_b + 10;

        int max_iteration=1000;

        double lambda = 0.5 * (lambda_a + lambda_b);

        for(int iter=0;iter<max_iteration;++iter){
            double sum=0;

            for(int i=0;i<n1*n2;++i){
                double eval =- phi[i] - nu[i] + lambda;
                if(nu[i] >= 0){
                    if(eval>0){
                        sum += 1.0; // IMPORTANT    
                    }else{
                        sum += fmin(mu[i],1);
                    }
                }
            }

            double val =sum /(1.0*n1*n2) - original_sum;

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

        

        for(int i=0;i<n1*n2;++i){
            double eval = - phi[i] - nu[i] + lambda;
            if(nu[i]>=0){
                if(eval>0){
                    DEstar[i] = 1.0;
                }else{
                    DEstar[i] = fmin(mu[i],1);
                }
            }else{
                DEstar[i] = 0;
            }
            
        }
    }

}; // Helper_E


class BackAndForth{
public:
    int n1;
    int n2;

    int max_iteration;
    double tolerance_scale;
    double W2_value;

    double C; // Coefficients for trance theorem
    double c1;
    double c2;

    double beta_1;
    double beta_2;
    double alpha_1;
    double alpha_2;

    double c1_phi;
    double c2_phi;

    double c1_psi;
    double c2_psi;

    double tau;

    double* vx;
    double* vy;

    // Initialize phi and psi
    double* phi;
    double* psi;

    double* gradx;

    double* push_mu;

    poisson_solver* fftps;
    FLT2D* flt2d;

    Pushforward_mapping* pushforward;

    BackAndForth(int n1, int n2, int max_iteration, double tolerance_scale, double tau){

        this->n1=n1;
        this->n2=n2;
        this->max_iteration=max_iteration;
        this->tolerance_scale=tolerance_scale;

        this->tau=tau;

        vx=new double[n1*n2];
        vy=new double[n1*n2];

        phi=new double[n1*n2];
        psi=new double[n1*n2];

        gradx=new double[n1*n2];

        push_mu=new double[n1*n2];

        beta_1=0.05;
        beta_2=0.95;
        alpha_1=1.0/0.99;
        alpha_2=1.0/alpha_1;
        alpha_1=1.2;
        alpha_2=0.8;

        // pushforward->initialize(n1,n2);
        // flt2d->initialize(n1,n2);

        pushforward = new Pushforward_mapping(n1,n2);
        flt2d = new FLT2D(n1,n2);

        clock_t time;
        time=clock();
        fftps = new poisson_solver(n1,n2);
        time=clock()-time;
        printf ("\nCPU time for FFT: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);
    }

    ~BackAndForth(){
        delete[] vx;
        delete[] vy;
        delete[] phi;
        delete[] psi;
        delete[] gradx;
        delete[] push_mu;

        delete pushforward;
        delete flt2d;
        delete fftps;
    }

    void calculate_gradient(const double* phi_c, double* gradx, double* grady){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1-1;++j){
                gradx[i*n1+j]=1.0*n1*(phi_c[i*n1+j+1]-phi_c[i*n1+j]);
            }
        }
        for(int j=0;j<n1;++j){
            for(int i=0;i<n2-1;++i){
                grady[i*n1+j]=1.0*n2*(phi_c[(i+1)*n1+j]-phi_c[i*n1+j]);
            }
        }
    }

    double calculate_dual_value(Helper_E& helper_e, const double* phi, const double* psi, const double* mu){

         // Here psi is assumed to correspond to the c-transform of phi
        int pcount=n1*n2;
        
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term1 += helper_e.calculate_E(phi,i,j);
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

    double calculate_h_minus_1(poisson_solver* fftps, const double* push_mu, const double* DEstar){
        double error=0;
        for(int i=0;i<n1*n2;++i){
            double value=-push_mu[i]+DEstar[i];
            error+=value*fftps->workspace[i];
        }
        return error/(1.0*n1*n2);
    }

    double initialize_sigma(const double* mu){
        double sigma=1.0;
        for(int i=0;i<n1*n2;++i){
            sigma=fmax(sigma,mu[i]);
        }
        return tau/sigma;
    }

    void calculate_push_rho(const double* rho, double* push_rho,const double* vx,const double* vy){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=0;
                double vyval=0;

                if(j==0){
                    vxval=0.5*(vx[i*n1+j]);
                }else{
                    vxval=0.5*(vx[i*n1+j]+vx[i*n1+j-1]);
                }

                if(i==0){
                    vyval=0.5*(vy[i*n1+j]);
                }else{
                    vyval=0.5*(vy[i*n1+j]+vy[(i-1)*n1+j]);
                }

                double x,y;

                x=(j+0.5)/(1.0*n1)-tau*vxval;
                y=(i+0.5)/(1.0*n2)-tau*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){

                    x=(j+0.5)/(1.0*n1);
                    y=(i+0.5)/(1.0*n2);

                    double xpost=x+1.0/n1;
                    double xpre =x;
                    double gradx_vx_post=interpolate_function(xpost,y,vx);
                    double gradx_vx_pre =interpolate_function(xpre,y,vx);
                    
                    double ypost=y+1.0/n2;
                    double ypre =y;
                    double grady_vy_post=interpolate_function(x,ypost,vy);
                    double grady_vy_pre =interpolate_function(x,ypre,vy);

                    double gradx_vx=1.0*n1*(gradx_vx_post-gradx_vx_pre);
                    double grady_vy=1.0*n2*(grady_vy_post-grady_vy_pre);

                    push_rho[i*n1+j]=rhovalue*fabs((1-tau*gradx_vx)*(1-tau*grady_vy)); 
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
                                   +(1-lambda1)*(lambda2)*x10+(lambda1)*(lambda2)*x11;

        return interpolated_value;  
    }

    double interpolate_function_gradient_vx(double x,double y,const double* func){
        double indj=x*n1-1.0;
        double indi=y*n2-0.5;

        double lambda1=indj-(int)indj;
        double lambda2=indi-(int)indi;

        double x00 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj))];
        double x01 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj+1))];
        double x10 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj))];
        double x11 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj+1))];

        if(indj>n1-1 || indj<0){
            x00=0;
            x10=0;
        }
        if(indj+1>n1-1 || indj+1<0){
            x01=0;
            x11=0;
        }
        if(indi>n2-1 || indi<0){
            x00=0;
            x01=0;
        }
        if(indi+1>n2-1 || indi+1<0){
            x10=0;
            x11=0;
        }

        double interpolated_value = (1-lambda1)*(1-lambda2)*x00+(lambda1)*(1-lambda2)*x01
                                   +(1-lambda1)*(lambda2)*x10+(lambda1)*(lambda2)*x11;

        return interpolated_value; 
    }


    double interpolate_function_gradient_vy(double x,double y,const double* func){
        double indj=x*n1-0.5;
        double indi=y*n2-1.0;

        double lambda1=indj-(int)indj;
        double lambda2=indi-(int)indi;

        double x00 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj))];
        double x01 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj+1))];
        double x10 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj))];
        double x11 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj+1))];

        if(indj>n1-1 || indj<0){
            x00=0;
            x10=0;
        }
        if(indj+1>n2-1 || indj+1<0){
            x01=0;
            x11=0;
        }
        if(indi>n2-1 || indi<0){
            x00=0;
            x01=0;
        }
        if(indi+1>n2-1 || indi+1<0){
            x10=0;
            x11=0;
        }

        double interpolated_value = (1-lambda1)*(1-lambda2)*x00+(lambda1)*(1-lambda2)*x01
                                   +(1-lambda1)*(lambda2)*x10+(lambda1)*(lambda2)*x11;

        return interpolated_value; 
    }

    /**
        Update sigma based on Goldstein scheme
    */
    double update_sigma(double sigma, const double W2_value, const double W2_value_previous, const double error)
    {
        // if(W2_value_previous-W2_value>-sigma*beta_1*error){
        //     return alpha_2;
        // }else if(W2_value_previous-W2_value<-sigma*beta_2*error){
        //     return alpha_1;
        // }
        // return 1;

        if(W2_value_previous-W2_value>-sigma*beta_1*error){
            sigma*=alpha_2;
        }else if(W2_value_previous-W2_value<-sigma*beta_2*error){
            sigma*=alpha_1;
        }
        return sigma;
    }


    double perform_OT_iteration_back_det(Helper_E& helper_e,double& sigma,double& W2_value,const double* mu){
        // ------------------------------------------------------------
        double W2_value_previous, error_nu;

        flt2d->find_c_concave(psi,phi,tau);
        flt2d->find_c_concave(phi,psi,tau);

        helper_e.calculate_DEstar(phi,mu);
        // helper_e.calculate_DEstar_normalized(phi,mu);

        calculate_gradient(psi, vx, vy);
        calculate_push_rho(helper_e.DEstar, push_mu,vx,vy);

        fftps->perform_inverse_laplacian(push_mu,mu,c1_psi,c2_psi);

        W2_value_previous=calculate_dual_value(helper_e,phi,psi,mu);
        error_nu=calculate_h_minus_1(fftps,push_mu,mu);

        for(int i=0;i<n1*n2;++i){
            psi[i]+=sigma*fftps->workspace[i];
        }

        flt2d->find_c_concave(psi,phi,tau);

        W2_value=calculate_dual_value(helper_e,phi,psi,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(Helper_E& helper_e,double& sigma,double& W2_value,const double* mu){
        // ------------------------------------------------------------
        double W2_value_previous, error_mu;

        flt2d->find_c_concave(phi,psi,tau);
        flt2d->find_c_concave(psi,phi,tau);
            
        calculate_gradient(phi, vx, vy);
        calculate_push_rho(mu, push_mu,vx,vy);
        
        helper_e.calculate_DEstar(phi,mu);
        // helper_e.calculate_DEstar_normalized(phi,mu);


        fftps->perform_inverse_laplacian(push_mu,helper_e.DEstar,c1_phi,c2_phi);


        W2_value_previous=calculate_dual_value(helper_e,phi,psi,mu);
        error_mu=calculate_h_minus_1(fftps,push_mu,helper_e.DEstar);

        for(int i=0;i<n1*n2;++i){
            phi[i]+=sigma*fftps->workspace[i];
        }


        flt2d->find_c_concave(phi,psi,tau);

        W2_value=calculate_dual_value(helper_e,phi,psi,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu);

        return error_mu;
    }

    void display_iteration(const int iter,const double sigma_forth,const double sigma_back,const double W2_value,const double error_mu,const double error_nu,const double C_phi,const double C_psi,const double SCALER_forth,const double SCALER_back) const{
        cout <<setw(5)<<iter+1<<"  sigma : "<<scientific << setw(13) << sigma_forth  << setw(13) << sigma_back  << " C : " << C_phi << " " << C_psi << "  fftps->coeff : "<<scientific<<setw(13) << c2_phi/c1_phi << " " << c2_psi/c1_psi <<"  W2 : " <<scientific<<setw(13) << W2_value << "  h-1 : "<<scientific<<setw(13) << error_mu << " "<< error_nu << endl;
        cout << fixed;
    }

    /**
        Calculate a = 0.1 * max(-phi)
    */
    double calculate_a() const{
        // double MIN = 0.05;
        double MIN = 1;
        double phimax = -phi[0];
        double phimin = -phi[0];

        for(int i=1;i<n1*n2;++i){
            phimax = fmax(phimax, -phi[i]);
            phimin = fmin(phimin, -phi[i]);
        }

        return 0.02*phimax + 0.98*phimin;
    }

    /**
        calculate the area of a level set given an indicator function u
    */
    double calculate_area_of_level_set(const double* u){
        double area = 0;
        for(int i=0;i<n1*n2;++i){
            if(u[i] < FLT_MAX) area += u[i];
        }
        return area/(n1*n2);
    }
    /**
        calculate the perimenter of a level set given an indicator function u
    */
    double calculate_perimeter_of_level_set(const double* u){

        fftps->perform_convolution(u);
        fftps->eps = 1.0/pow(n1,0.3);
        // fftps->eps = 0.01;

        double per = 0;
        for(int i=0;i<n1*n2;++i){
            per += (1-u[i]) * fmax(0,fftps->workspace[i]);
        }
        return per/(n1*n2*fftps->eps);
    }

    /**
        calculate F for Newton's method for bound
    */
    void calculate_F_Newton_bound(double& F, double& DF, const double w, const double area, const double per, const double* g){
        F  = 0;
        DF = 0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                if(i==0 && j==0) continue;
                double negativeLaplacian=2.0*n1*n1*(1-cos(M_PI*(j)/n1)) + 2*n2*n2*(1-cos(M_PI*(i)/n2));
                double eval  = 1/w + per/area + w * negativeLaplacian;
                double eval2 = -1/(w*w) + negativeLaplacian;
                F  += - g[i*n1+j]*g[i*n1+j] * exp(-2*log(eval)) * eval2;
                DF += g[i*n1+j]*g[i*n1+j] * (
                        2 * exp(-3*log(eval)) * eval2 * eval2 - exp(-2*log(eval)) * 2 / (w * w * w)
                    );
            }
        }

        F  /= n1*n2;
        DF /= n1*n2;
    }

    double calculate_F_bisection(const double w, const double area, const double per, const double* g, const double gradunorm){
        double F  = 0;

        for(int i=1; i<n1*n2; ++i){
            double negativeLaplacian=fftps->kernel[i];
            double eval  = gradunorm/w + per/area + (w * gradunorm + tau) * negativeLaplacian;
            double eval2 = -gradunorm/(w*w) + gradunorm * negativeLaplacian;
            F  += - g[i]*g[i] / (eval * eval) * eval2;
        }

        return F/(n1*n2);
    }

    /**
        optimize w using newton's method / bijection method
    */
    double optimize_w(const double area, const double per, const double* mu, const double* g, const double gradunorm){

        int MAX_ITERATION = 50;
        double TOLERANCE  = 1e-8;

        double min_bound = 1e-3;
        double max_bound = 0.3;
        
        double min_val = calculate_F_bisection(min_bound,area,per,g,gradunorm);
        double max_val = calculate_F_bisection(max_bound,area,per,g,gradunorm);

        // min_val needs to be positive
        while(min_val < 0 && min_bound > 1e-5){
            min_bound *= 0.7;
            min_val = calculate_F_bisection(min_bound,area,per,g,gradunorm);
        }

        // max_val needs to be negative
        while(max_val > 0 && max_bound < 100){
            max_bound *= 1.5;
            max_val = calculate_F_bisection(max_bound,area,per,g,gradunorm);
        }

        // cout << "min : " << min_bound << " " << min_val << " max : " << max_bound << " " << max_val << endl;

        if(max_val * min_val > 0){
            return 0.01;  
        } 

        // int N = 2000;
        // for(int i=0;i<N;++i){
        //     double t = 1e-3 + 1e-3*i;
        //     cout << calculate_F_bisection(t,area,per,g,gradunorm);
        //     if(i<N-1){
        //         cout << ",";
        //     }
        // }
        // cout << endl;

        double mid_bound = 0.5 * (min_bound + max_bound);
        double mid_val = calculate_F_bisection(mid_bound,area,per,g,gradunorm);

        for(int iter=0; iter<MAX_ITERATION; ++iter){

            if(fabs(mid_val) < TOLERANCE) break;

            if(mid_val < 0){
                max_bound = mid_bound;
            }else{
                min_bound = mid_bound;
            }

            mid_bound = 0.5 * (min_bound + max_bound);

            mid_val = calculate_F_bisection(mid_bound,area,per,g,gradunorm);


        }

        return mid_bound;
    }

    /** 
        calculate area and perimenter
    */
    void calculate_area_perimeter(const double* nu, double& area, double& per){

        // double a = 1e-4; // for incompressible
        // double infgradphi = calculate_infgradphi_on_level_set_a(a);

        // u is an indicator function (will use gradx to avoid creating additional double array)
        for(int i=0;i<n1*n2;++i){
            gradx[i] = 0;
            if(nu[i] >= 0){
                double eval = -phi[i] - nu[i];
                if(eval > 0){
                    gradx[i] = 1;
                }
            }
        }

        // calculate the area
        area = calculate_area_of_level_set(gradx);

        // calculate the perimeter
        per  = calculate_perimeter_of_level_set(gradx);        
    }

    /** 
        calculate c1 and c2 for coefficients for the operator c1 Id + c2 Laplacian
        C is constant for the trace theorem
    */
    void calculate_c1_c2_coeff_psi(const double area, const double per, const double C, const bool verbose, const double* mu, Helper_E& helper_e, const double infgradphi){

        double w = 0.1;
        c1_psi = C/w/infgradphi;

        // double mumax = 1; // incompressible case

        if(area > 0){

            // optimize w using newton's method
            helper_e.calculate_DEstar(phi);
            calculate_gradient(psi, vx, vy);
            calculate_push_rho(helper_e.DEstar, push_mu,vx,vy);
            fftps->get_fourier_coefficients(push_mu,mu);
            
            w = optimize_w(area,per,mu,fftps->workspace,C);
            // w = 0.05;
            c1_psi = (C/w + per/area)/infgradphi;
        }

        c2_psi = C*w/infgradphi;

        if(verbose){
            double phimax = -psi[0]-helper_e.nu[0];
            double phimin = -psi[0]-helper_e.nu[0];
            for(int i=1;i<n1*n2;++i){
                if(helper_e.nu[i] >= 0){
                    phimax=fmax(phimax,-psi[i]-helper_e.nu[i]);
                    phimin=fmin(phimin,-psi[i]-helper_e.nu[i]);
                }
            }

            cout << "PSI w : " << w << " area : " << area << " per : " << per << " c1 : " << c1_psi << " c2 : " << c2_psi  << " c2/c1 : " << c2_psi/c1_psi << endl;    
            cout << endl;
        }
        
    }

    /** 
    calculate c1 and c2 for coefficients for the operator c1 Id + c2 Laplacian
    C is constant for the trace theorem
    */
    void calculate_c1_c2_coeff_phi(const double area, const double per, const double C, const bool verbose, const double* mu, Helper_E& helper_e, const double infgradphi){

        double w = 0.1;
        c1_phi = C/w;

        if(area > 0){
            // optimize w using newton's method
            calculate_gradient(phi, vx, vy);
            calculate_push_rho(mu, push_mu,vx,vy);
            helper_e.calculate_DEstar(phi);
            fftps->get_fourier_coefficients(push_mu,helper_e.DEstar);
            
            w = optimize_w(area,per,mu,fftps->workspace,C);
            // w = 0.05;
            c1_phi = C/w + per/area;
        }
        
        c1_phi /= infgradphi;

        c2_phi = C*w/infgradphi;

        if(verbose){
            double phimax = -phi[0]-helper_e.nu[0];
            double phimin = -phi[0]-helper_e.nu[0];
            for(int i=1;i<n1*n2;++i){
                if(helper_e.nu[i] >= 0){
                    phimax=fmax(phimax,-phi[i]-helper_e.nu[i]);
                    phimin=fmin(phimin,-phi[i]-helper_e.nu[i]);
                }
            }

            cout << "PHI w : " << w << " area : " << area << " per : " << per << " c1 : " << c1_phi << " c2 : " << c2_phi  << " c2/c1 : " << c2_phi/c1_phi << " phimin : " << phimin << " phimax : " << phimax << endl;    
        }
        
    }

    /**
        caclulate the coefficient C for the operator I + C Laplacian
    */
    void initialize_coeff_for_fftps(poisson_solver* fftps, const double* mu, const double c1, const double c2, const double C){

        double maxrho=1.0;
        for(int i=0;i<n1*n2;++i){
            maxrho=fmax(maxrho,mu[i]);
        }
        fftps->coeff = (c2 + tau * maxrho) / c1;    
    }

    void initialize_phi(Helper_E& helper_e,const double* mu){
        // TODO
        for(int i=0;i<n1*n2;++i){
            phi[i] = -(mu[i] + helper_e.nu[i]);
        }
    }

    /** 
        Calculate infgradphi = inf(|nabla phi|)
    */
    double calculate_infgradphi_on_level_set_a(const double a) const{

        double MIN = 1e-2;

        double infgradphi=1000;

        int count = 0;

        for(int i=1;i<n2;++i){
            for(int j=1;j<n1;++j){
                // if( 
                //       (
                //          (   (-phi[i*n1+j]-a) > 0 && (-phi[i*n1+j-1]-a) < 0  ) || (   (-phi[i*n1+j]-a) < 0 && (-phi[i*n1+j-1]-a) > 0  ) 
                //       )
                //       || 
                //       (
                //          (   (-phi[i*n1+j]-a) > 0 && (-phi[(i-1)*n1+j]-a) < 0  ) || (   (-phi[i*n1+j]-a) < 0 && (-phi[(i-1)*n1+j]-a) > 0  ) 
                //       )
                // )

                if(fabs(phi[i*n1+j]) < 1e-4)
                {
                    double gradxphi = 1.0*n1*(phi[i*n1+j]-phi[i*n1+(int) fmax(0,j-1)]);
                    double gradyphi = 1.0*n2*(phi[i*n1+j]-phi[(int) fmax(0,i-1)*n1+j]);
                    double eval = sqrt(gradxphi*gradxphi + gradyphi*gradyphi);

                    if(eval > 0){
                        infgradphi = fmin(infgradphi, eval);
                    }   
                }
            }
        }

        return fmax(0.001,infgradphi);
    }

    void start_OT(Helper_E& helper_e, const double* mu, const int outer_iter){

        int skip = 10; // frequency of printout

        double error_mu,error_nu,error=1.0;

        /*
            Initialize coefficients for siga update
        */

        beta_1=0.1;
        beta_2=0.9;

        alpha_1=1.1;
        alpha_2=0.9;

        /*
            Initialize the tolerance based on tau^2
        */

        double tol_modified = 0;
        for(int i=0;i<n1*n2;++i){
            tol_modified = fmax(tol_modified, mu[i]);
        }
        tol_modified = tolerance_scale * tol_modified*tau*tau;

        cout << "iter = " << outer_iter << "\ttolerance = " << scientific << tol_modified  << endl;
        cout << fixed;

        /*
            Initialize phi for the first iteration only
        */

        double SCALER_forth = 0.1;
        double SCALER_back  = 0.1;
        
        double C_phi = 10;
        double C_psi = 10;

        if(outer_iter==0){
            // initialize_phi(helper_e,mu); // intiailize phi in the first outer iteration            
            c1_phi = 10;
            c2_phi = 1e-3;

            c1_psi = 10;
            c2_psi = 1e-3;
            SCALER_forth = 0.1;
            SCALER_back  = 0.1;
        }

        /*
            Initialize the coefficients for fftps
        */
        double area=0;
        double per =0;

        double infgradphi = calculate_infgradphi_on_level_set_a(0);
        calculate_area_perimeter(helper_e.nu,area,per);
        calculate_c1_c2_coeff_phi(area,per,C_phi,true,mu,helper_e,infgradphi);
        calculate_c1_c2_coeff_psi(area,per,C_psi,true,mu,helper_e,infgradphi);

        /*
            Initialize the sigma
        */

        double sigma_forth = SCALER_forth;
        double sigma_back  = SCALER_back;

        /*
            Starting the loop
        */
    
        for(int iter=0;iter<max_iteration;++iter){

            if(iter<20 && outer_iter == 0){
                fftps->coeff = 0.01;
                sigma_forth=0.001;
                sigma_back=0.001;
            }

                // Perform OT iteration back and forth         
        
            error_mu=perform_OT_iteration_forth_det(helper_e,sigma_forth,W2_value,mu);
            error_nu=perform_OT_iteration_back_det(helper_e,sigma_back,W2_value,mu);
            
                // Calculating the relative error
                        
            error=fmin(error_mu,error_nu);

            if(sigma_forth/SCALER_forth>1){
                C_phi = fmax(0.1,C_phi*0.99);
            }else if(sigma_forth/SCALER_forth<1){
                C_phi = fmin(10,C_phi*1.01);
            }

            if(sigma_back/SCALER_back>1){
                C_psi = fmax(0.1,C_psi*0.99);
            }else if(sigma_back/SCALER_back<1){
                C_psi = fmin(10,C_psi*1.01);
            }

            if(iter%1 == 0 && iter>0){
                // cout << "C_phi : " << C_phi << " C_psi : " << C_psi << endl;
                infgradphi = calculate_infgradphi_on_level_set_a(0);
                calculate_area_perimeter(helper_e.nu,area,per);
                calculate_c1_c2_coeff_phi(area,per,C_phi,false,mu,helper_e,infgradphi);
                calculate_c1_c2_coeff_psi(area,per,C_psi,false,mu,helper_e,infgradphi); 
            }            

            sigma_forth = fmin(20,fmax(10, sigma_forth));
            sigma_back  = fmin(20,fmax(10, sigma_back ));

            SCALER_forth = sigma_forth;
            SCALER_back  = sigma_back ;

            /* Display the result per iterations */

            cout<<"-"<<flush;   

            if(iter%skip==skip-1){
                cout<<"|";
                display_iteration(iter,sigma_forth,sigma_back,W2_value,error_mu,error_nu,C_phi,C_psi,SCALER_forth,SCALER_back);

                /* calculate infgradphi */
                infgradphi = calculate_infgradphi_on_level_set_a(0);
                cout << "infgradphi : " << infgradphi << "\n";
            }

            /* Stopping Condition */

            // if(W2_value>0 && ((abs(error)<tolerance && iter>=0) || iter==max_iteration-1 || sigma_forth <1e-9) ){
            if(((abs(error)<tol_modified && abs(error)>0 && iter>=5) || iter==max_iteration-1) ){




                cout<<"Tolerance met!"<<endl;
                display_iteration(iter,sigma_forth,sigma_back,W2_value,error_mu,error_nu,C_phi,C_psi,SCALER_forth,SCALER_back);


                /* calculate infgradphi */
                double infgradphi = calculate_infgradphi_on_level_set_a(0);
                cout << "infgradphi : " << infgradphi << "\n";

                break;
            }
        }
    }


}; // Back and Forth

int main(int argc, char** argv){
    if(argc!=7){
        cout<<"Do the following:"<<endl;
        cout<<"./main.exe [n1] [n2] [max_iteration] [tolerance_scale] [nt] [tau]"<<endl;
        return 0;
    }

    int n1=stoi(argv[1]);
    int n2=stoi(argv[2]);
    int max_iteration=stoi(argv[3]);
    double tolerance_scale=stod(argv[4]);
    int nt=stoi(argv[5]);
    double tau=stod(argv[6]);

    string data_folder = "./data";

    create_csv_parameters(n1,n2,nt,tau,data_folder);

    // Initialize mu and nu
    double* mu=new double[n1*n2]; 

    cout << "XXX Starting Incompressible Gradient Flow XXX" << endl;

    cout << endl;

    cout << "The files will be saved in " << data_folder << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);
  
    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));
    cout << endl;


    /*
        Print out parameters
    */

    cout << "n1   : " << n1 <<endl;
    cout << "n2   : " << n2 <<endl;
    cout << "nt   : " << nt <<endl;
    cout << "tau  : " << tau << endl;
    cout << "Max Iteration   : " << max_iteration <<endl;
    cout << "tolerance_scale : " << tolerance_scale <<endl;

    Helper_E helper_e(n1,n2,tau);

    // Initialize mu
    // create_mu_gaussian(mu,0.25,0.25,0.1,n1,n2);
    // create_mu_square(mu,0.5,0.5,0.1,n1,n2);
    // create_mu_one_square(mu,n1,n2);
    // create_mu_two_square(mu,n1,n2);
    // create_mu_four_square(mu,n1,n2);
    // create_mu_one_circle(mu,n1,n2);
    // create_mu_two_circles(mu,n1,n2);
    // create_mu_four_circles(mu,n1,n2);
    // create_mu_from_image(mu,n1,n2);
    // create_mu_flavien_test(mu,n1,n2);

    cout << setprecision(6);

    unsigned char* obstacle = new unsigned char[n1*n2];

    init_obstacle_from_image(obstacle, n1, n2);
    // init_obstacle_two_moons(obstacle, n1, n2);
    // init_obstacle_pac_man(obstacle, n1, n2);
    init_entropy_image_obstacle_opencv(helper_e.nu, obstacle, n1, n2, data_folder);

    // for(int i=0;i<n1*n2;++i){
    //     if(obstacle[i] > 0){
    //         mu[i] = 0;
    //     }
    // }

    helper_e.setup_original_sum(mu);

    // Creating Video

    int pcount=n1*n2;
    unsigned char* pixels=new unsigned char[pcount];

    double max_video=0;
    for(int i=0;i<n1*n2;++i){
        max_video=fmax(max_video,mu[i]);
    }

    cout << "max : " << max_video << endl;

    char command[1000]; // Make it large enough.
    sprintf(command, "ffmpeg -hide_banner -loglevel panic -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s %dx%d -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 interpolate.mp4", n2,n1);
    FILE *pipeout=popen(command, "w");

    prepare_pixels(mu,pixels,pcount,max_video);
    fwrite(pixels, sizeof(unsigned char), pcount, pipeout);

    BackAndForth bf(n1,n2,max_iteration,tolerance_scale,tau);

    string filename = data_folder + "/mu-"+to_string(0)+".csv";
    create_bin_file(mu,n1*n2,filename);

    string figurename = "circle-potential-center";
    save_image_opencv(mu,obstacle,figurename,0,n1,n2);

    clock_t time;
    time=clock();

    double sum=0;

    for(int n=0;n<nt;++n){

        bf.start_OT(helper_e, mu, n);
        // helper_e.calculate_DEstar(bf.phi);
        // helper_e.calculate_DEstar_normalized(bf.phi);
        helper_e.calculate_DEstar_normalized(bf.phi,mu);
        memcpy(mu,helper_e.DEstar,n1*n2*sizeof(double));

        sum=0;
        for(int i=0;i<n1*n2;++i) sum += mu[i];

        cout << "sum : " << sum/(n1*n2) << endl;

        prepare_pixels(mu,pixels,pcount,max_video); fwrite(pixels, sizeof(unsigned char), pcount, pipeout);

        filename = data_folder + "/mu-"+to_string(n+1)+".csv"; create_bin_file(mu,n1*n2,filename);
        create_bin_file(mu,n1*n2,filename);

        // save_image_opencv(mu,obstacle,figurename,n+1,n1,n2);
        save_image_opencv(mu,obstacle,figurename,n+1,n1,n2);
    }



    time=clock()-time;
    printf ("\nCPU time for GF: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    // filename = data_folder + "/nu.csv";
    // create_csv_file(helper_e.nu,n1*n2,filename);
    // filename = data_folder + "/phi.csv";
    // create_csv_file(bf.phi,n1*n2,filename);
    // filename = data_folder + "/psi.csv";
    // create_csv_file(bf.psi,n1*n2,filename);
    // filename = data_folder + "/push_mu.csv";
    // create_csv_file(mu,n1*n2,filename);

    fclose (pipeout);

    delete[] mu;
    delete[] pixels;
}