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

using namespace std;

class Helper_E{
public:
    int n1;
    int n2;

    double gamma;
    double tau;
    double m;

    double M;

    double* nu;
    double* DEstar;

    Helper_E(){
        nu=NULL;
        DEstar=NULL;
    }

    Helper_E(int n1,int n2,double gamma,double tau,double m,double M){
        this->n1=n1;
        this->n2=n2;
        this->gamma=gamma;
        this->tau=tau;
        this->m=m;
        this->M=M;

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
        double exponent=m/(m-1);
        if(c_phi > 0 ){
            return 1.0/(exponent*pow(gamma,1/(m-1))) * exp(exponent*log(c_phi));    
        }
        return 0;
    }

    void calculate_DEstar(const double* phi){
        
        double exponent = 1.0/(m-1);
        double gammaprime = pow(gamma,1/(m-1));

        for(int i=0;i<n1*n2;++i){
            // DEstar[i] = 1.0/gamma * pow(fmax(0, - phi[i]/tau - nu[i]), 1.0/(m-1));
            DEstar[i] = 1.0/gammaprime * exp(exponent * log(fmax(0, - phi[i] - nu[i])));
        }
    }

    void calculate_DEstar_normalized(const double* phi, const double tolerance=1e-7){
        
        double lambda_a=-phi[0]-nu[0];
        double lambda_b=-phi[0]-nu[0];
        double gammaprime = pow(gamma,1/(m-1));

        for(int i=0;i<n1*n2;++i){
            lambda_a = fmax(lambda_a, - phi[i] - nu[i]);
            lambda_b = fmin(lambda_b, - phi[i] - nu[i]);
        }

        lambda_a = - lambda_a - 1;
        lambda_b = - lambda_b + gammaprime + 1;

        int max_iteration=100;

        double lambda = 0.5 * (lambda_a + lambda_b);

        double exponent=1.0/(m-1);

        for(int iter=0;iter<max_iteration;++iter){
            double sum=0;

            for(int i=0;i<n1*n2;++i){
                double eval=fmax(0, - phi[i] - nu[i] + lambda);
                if(eval>0){
                    sum+=exp(exponent*log(eval)); // IMPORTANT    
                }
            }
            sum/=gammaprime;

            double val =sum /(1.0*n1*n2)- M;

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
            DEstar[i] = 1.0/gammaprime * pow(fmax(0, - phi[i] - nu[i] + lambda), 1.0/(m-1));
        }
    }

}; // Helper_E


class BackAndForth{
public:
    int n1;
    int n2;

    int max_iteration;
    double tolerance;
    double W2_value;

    double C; // Coefficients for trance theorem
    double c1;
    double c2;

    double phi_c1;
    double phi_c2;
    double psi_c1;
    double psi_c2;


    double beta_1;
    double beta_2;
    double alpha_1;
    double alpha_2;

    double gamma;
    double tau;
    double m;
    double mprime;

    double SCALER_forth;
    double SCALER_back ;
    
    double C_phi;
    double C_psi;

    double* vx;
    double* vy;

    double* vxtmp;
    double* vytmp;

    // Initialize phi and psi
    double* phi;
    double* psi;

    double* push_mu;

    poisson_solver* fftps;
    FLT2D* flt2d;

    // Pushforward_mapping* pushforward;

    BackAndForth(int n1, int n2, int max_iteration, double tolerance, double gamma, double tau, double m){

        this->n1=n1;
        this->n2=n2;
        this->max_iteration=max_iteration;
        this->tolerance=tolerance;

        this->gamma=gamma;
        this->tau=tau;
        this->m=m;

        mprime = m/(m-1);

        vx=new double[n1*n2];
        vy=new double[n1*n2];

        vxtmp=new double[n1*n2];
        vytmp=new double[n1*n2];

        phi=new double[n1*n2];
        psi=new double[n1*n2];

        push_mu=new double[n1*n2];

        // pushforward = new Pushforward_mapping(n1,n2);

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
        delete[] vxtmp;
        delete[] vytmp;
        delete[] phi;
        delete[] psi;
        delete[] push_mu;

        // delete pushforward;
        delete flt2d;
        delete fftps;
    }

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

    double calculate_grady_WENO_(const double* phi, const int i, const int j){
        double phi_i_pre_3  = phi[(int)fmax(0,i-3)*n1+j];
        double phi_i_pre_2  = phi[(int)fmax(0,i-2)*n1+j];
        double phi_i_pre_1  = phi[(int)fmax(0,i-1)*n1+j];
        double phi_i        = phi[i*n1+j];
        double phi_i_post_1 = phi[(int)fmin(n2-1,i+1)*n1+j];
        double phi_i_post_2 = phi[(int)fmin(n2-1,i+2)*n1+j];
        double eval = n2/12.0 * (
                                 -         (phi_i_pre_1  - phi_i_pre_2)
                                 + 7.0  *  (phi_i        - phi_i_pre_1)
                                 + 7.0  *  (phi_i_post_1 - phi_i)
                                 -         (phi_i_post_2 - phi_i_post_1)
                                );
        double a = n2 * (phi_i_pre_1 - 2.0 * phi_i_pre_2 + phi_i_pre_3);
        double b = n2 * (phi_i - 2.0 * phi_i_pre_1 + phi_i_pre_2);
        double c = n2 * (phi_i_post_1 - 2.0 * phi_i + phi_i_pre_1);
        double d = n2 * (phi_i_post_2 - 2.0 * phi_i_post_1 + phi_i);
        eval -= calculate_WENO(a,b,c,d);

        return eval;

        
    }

    double calculate_gradx_vx(const double* phi, const int i, const int j){
        double phi_j_pre_3  = phi[i*n1+ (int)fmax(0,j-3)];
        double phi_j_pre_2  = phi[i*n1+ (int)fmax(0,j-2)];
        double phi_j_pre_1  = phi[i*n1+ (int)fmax(0,j-1)];
        double phi_j        = phi[i*n1+j];
        double phi_j_post_1 = phi[i*n1+ (int)fmin(n1-1,j+1)];
        double phi_j_post_2 = phi[i*n1+ (int)fmin(n1-1,j+2)];
        double phi_j_post_3 = phi[i*n1+ (int)fmin(n1-1,j+3)];

        double eval = 3.0/4.0 * (phi_j_post_1 - phi_j_pre_1) - 3.0/20.0 * (phi_j_post_2 - phi_j_pre_2) + 1.0/60.0 * (phi_j_post_3 - phi_j_pre_3);
        // double eval = 1.0/4.0 * ((phi_j_post_1 - phi_j_pre_1) + (phi_j_post_2 - phi_j_pre_2) + (phi_j_post_3 - phi_j_pre_3));

        return 1.0*n1*eval;
    }

    double calculate_grady_vy(const double* phi, const int i, const int j){
        double phi_i_pre_3  = phi[(int)fmax(0,i-3)*n1+j];
        double phi_i_pre_2  = phi[(int)fmax(0,i-2)*n1+j];
        double phi_i_pre_1  = phi[(int)fmax(0,i-1)*n1+j];
        double phi_i        = phi[i*n1+j];
        double phi_i_post_1 = phi[(int)fmin(n2-1,i+1)*n1+j];
        double phi_i_post_2 = phi[(int)fmin(n2-1,i+2)*n1+j];
        double phi_i_post_3 = phi[(int)fmin(n2-1,i+3)*n1+j];

        double eval = 3.0/4.0 * (phi_i_post_1 - phi_i_pre_1) - 3.0/20.0 * (phi_i_post_2 - phi_i_pre_2) + 1.0/60.0 * (phi_i_post_3 - phi_i_pre_3);
        // double eval = 1.0/4.0 * ((phi_i_post_1 - phi_i_pre_1) + (phi_i_post_2 - phi_i_pre_2) + (phi_i_post_3 - phi_i_pre_3));

        return 1.0*n2*eval;
        
    }

    void calculate_gradient2(const double* phi, double* vx, double* vy){
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                // vx[i*n1+j] = calculate_gradx_WENO(phi, i, j);
                // vy[i*n1+j] = calculate_grady_WENO(phi, i, j);
                vx[i*n1+j] = calculate_gradx_vx(phi, i, j);
                vy[i*n1+j] = calculate_grady_vy(phi, i, j);
            }
        }
    }

    void calculate_gradient(const double* phi_c, double* vx, double* vy){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1-1;++j){
                vx[i*n1+j]=1.0*n1*(phi_c[i*n1+j+1]-phi_c[i*n1+j]);
            }
        }
        for(int j=0;j<n1;++j){
            for(int i=0;i<n2-1;++i){
                vy[i*n1+j]=1.0*n2*(phi_c[(i+1)*n1+j]-phi_c[i*n1+j]);
            }
        }
    }

    double calculate_dual_value(Helper_E& helper_f, const double* phi, const double* psi, const double* mu){

         // Here psi is assumed to correspond to the c-transform of phi
        int pcount=n1*n2;
        
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term1 += helper_f.calculate_E(phi,i,j);
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

    void calculate_push_rho(const double* rho, double* push_rho,const double* vx,const double* vy,const double* vxtmp,const double* vytmp){

        double eps = pow(1.0/n1, 0.25);

        double xpost,ypost,xpre,ypre;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=vx[i*n1+j];
                double vyval=vy[i*n1+j];

                double x=(j+0.5)/(1.0*n1)-tau*vxval;
                double y=(i+0.5)/(1.0*n2)-tau*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){
                    xpre = x - 1.0/n1;
                    ypre = y - 1.0/n2;
                    double vxval_pre_1 = interpolate_function_v(xpre,y,vxtmp);
                    double vyval_pre_1 = interpolate_function_v(x,ypre,vytmp);
                    // xpre = x - 2.0/n1;
                    // ypre = y - 2.0/n2;
                    // double vxval_pre_2 = interpolate_function_v(xpre,y,vxtmp);
                    // double vyval_pre_2 = interpolate_function_v(x,ypre,vytmp);
                    // xpre = x - 3.0/n1;
                    // ypre = y - 3.0/n2;
                    // double vxval_pre_3 = interpolate_function_v(xpre,y,vxtmp);
                    // double vyval_pre_3 = interpolate_function_v(x,ypre,vytmp);

                    xpost = x + 1.0/n1;
                    ypost = y + 1.0/n2;
                    double vxval_post_1 = interpolate_function_v(xpost,y,vxtmp);
                    double vyval_post_1 = interpolate_function_v(x,ypost,vytmp);
                    // xpost = x + 2.0/n1;
                    // ypost = y + 2.0/n2;
                    // double vxval_post_2 = interpolate_function_v(xpost,y,vxtmp);
                    // double vyval_post_2 = interpolate_function_v(x,ypost,vytmp);
                    // xpost = x + 3.0/n1;
                    // ypost = y + 3.0/n2;
                    // double vxval_post_3 = interpolate_function_v(xpost,y,vxtmp);
                    // double vyval_post_3 = interpolate_function_v(x,ypost,vytmp);

                    // double gradx_vx = 1.0*n1* (3.0/4.0 * (vxval_post_1 - vxval_pre_1) - 3.0/20.0 * (vxval_post_2 - vxval_pre_2) + 1.0/60.0 * (vxval_post_3 - vxval_pre_3));
                    // double grady_vy = 1.0*n2* (3.0/4.0 * (vyval_post_1 - vyval_pre_1) - 3.0/20.0 * (vyval_post_2 - vyval_pre_2) + 1.0/60.0 * (vyval_post_3 - vyval_pre_3));

                    // double gradx_vx = 1.0*n1* (3.0/4.0 * (vxval_post_1 - vxval_pre_1) - 3.0/20.0 * (vxval_post_2 - vxval_pre_2) );
                    // double grady_vy = 1.0*n2* (3.0/4.0 * (vyval_post_1 - vyval_pre_1) - 3.0/20.0 * (vyval_post_2 - vyval_pre_2) );

                    double gradx_vx = 1.0*n1* (0.5 * (vxval_post_1 - vxval_pre_1));
                    double grady_vy = 1.0*n2* (0.5 * (vyval_post_1 - vyval_pre_1));

                    // push_rho[i*n1+j]=rhovalue/fabs((1.0-tau*gradx_vx) * (1.0-tau*grady_vy)); 
                    double eval = fabs((1.0-tau*gradx_vx) * (1.0-tau*grady_vy));
                    eval = fmax(eps,eval);
                    push_rho[i*n1+j] = rhovalue/eval;


                ////////////////////////////////////
                    // double vxval_post=vx[i*n1+(int)fmin(n1-1,j+1)];
                    // double vyval_post=vy[(int)fmin(n2-1,i+1)*n1+j];

                    // double gradx_vx=1.0*n1*(vxval_post-vxval);
                    // double grady_vy=1.0*n2*(vyval_post-vyval);

                    // double gradx_vx=calculate_gradx_WENO_(vx, i, j);
                    // double grady_vy=calculate_grady_WENO_(vy, i, j);
                ////////////////////////////////////
                    // double gradx_vx=calculate_gradx_vx(vx, i, j);
                    // double grady_vy=calculate_grady_vy(vy, i, j);

                    // push_rho[i*n1+j]=rhovalue*fabs((1.0-tau*gradx_vx) * (1.0-tau*grady_vy)); 
                ////////////////////////////////////

                    
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

        // double x02 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj+2))];
        // double x12 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj+2))];
        // double x20 = func[(int)fmin(n2-1,fmax(0,indi+2))*n1+(int)fmin(n1-1,fmax(0,indj))];
        // double x21 = func[(int)fmin(n2-1,fmax(0,indi+2))*n1+(int)fmin(n1-1,fmax(0,indj+1))];
        // double x22 = func[(int)fmin(n2-1,fmax(0,indi+2))*n1+(int)fmin(n1-1,fmax(0,indj+2))];
        
        // double x0 = ((int)indj + 0.5)/n1;
        // double x1 = x0 + 1.0/n1;
        // double x2 = x0 + 2.0/n1;

        // double y0 = ((int)indi + 0.5)/n2;
        // double y1 = y0 + 1.0/n2;
        // double y2 = y0 + 2.0/n2;

        // double interpolated_value = ((x1-x)*(x2-x))/((x1-x0)*(x2-x0)) * ((y1-y)*(y2-y))/((y1-y0)*(y2-y0)) * x00
        //                            +((x0-x)*(x2-x))/((x0-x1)*(x2-x1)) * ((y1-y)*(y2-y))/((y1-y0)*(y2-y0)) * x01
        //                            +((x1-x)*(x0-x))/((x1-x2)*(x0-x2)) * ((y1-y)*(y2-y))/((y1-y0)*(y2-y0)) * x02

        //                            +((x1-x)*(x2-x))/((x1-x0)*(x2-x0)) * ((y0-y)*(y2-y))/((y0-y1)*(y2-y1)) * x10
        //                            +((x0-x)*(x2-x))/((x0-x1)*(x2-x1)) * ((y0-y)*(y2-y))/((y0-y1)*(y2-y1)) * x11
        //                            +((x1-x)*(x0-x))/((x1-x2)*(x0-x2)) * ((y0-y)*(y2-y))/((y0-y1)*(y2-y1)) * x12

        //                            +((x1-x)*(x2-x))/((x1-x0)*(x2-x0)) * ((y1-y)*(y0-y))/((y1-y2)*(y0-y2)) * x20
        //                            +((x0-x)*(x2-x))/((x0-x1)*(x2-x1)) * ((y1-y)*(y0-y))/((y1-y2)*(y0-y2)) * x21
        //                            +((x1-x)*(x0-x))/((x1-x2)*(x0-x2)) * ((y1-y)*(y0-y))/((y1-y2)*(y0-y2)) * x22;
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
            // x11 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj+1))];
            // x01 = func[(int)fmin(n2-1,fmax(0,indi))*n1+(int)fmin(n1-1,fmax(0,indj+1))];
            // x10 = func[(int)fmin(n2-1,fmax(0,indi+1))*n1+(int)fmin(n1-1,fmax(0,indj))];

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
    double update_sigma(double sigma, const double W2_value, const double W2_value_previous, const double error){

        if(W2_value_previous-W2_value>-beta_1*error){
            return alpha_2;
        }else if(W2_value_previous-W2_value<-beta_2*error){
            return alpha_1;
        }
        return 1;
    }


    double perform_OT_iteration_back_det(Helper_E& helper_f,double& sigma,double& W2_value,const double* mu){
        // ------------------------------------------------------------
        double W2_value_previous, error_nu;

        flt2d->find_c_concave(psi,phi,tau);
        flt2d->find_c_concave(phi,psi,tau);

        helper_f.calculate_DEstar(phi);
        // helper_f.calculate_DEstar_normalized(phi);

        calculate_gradient(psi, vx, vy);
        calculate_gradient(phi, vxtmp, vytmp);
        calculate_push_rho(helper_f.DEstar, push_mu,vx,vy,vxtmp,vytmp);

        fftps->perform_inverse_laplacian(push_mu,mu,psi_c1,psi_c2,sigma);

        W2_value_previous=calculate_dual_value(helper_f,phi,psi,mu);
        error_nu=calculate_h_minus_1(fftps,push_mu,mu);

        for(int i=0;i<n1*n2;++i){
            psi[i] += fftps->workspace[i];
        }

        flt2d->find_c_concave(psi,phi,tau);

        W2_value=calculate_dual_value(helper_f,phi,psi,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(Helper_E& helper_f,double& sigma,double& W2_value,const double* mu){
        // ------------------------------------------------------------
        double W2_value_previous, error_mu;

        flt2d->find_c_concave(phi,psi,tau);
        flt2d->find_c_concave(psi,phi,tau);
            
        calculate_gradient(phi, vx, vy);
        calculate_gradient(psi, vxtmp, vytmp);
        calculate_push_rho(mu, push_mu,vx,vy,vxtmp,vytmp);
        
        helper_f.calculate_DEstar(phi);
        // helper_f.calculate_DEstar_normalized(phi);

        fftps->perform_inverse_laplacian(push_mu,helper_f.DEstar,phi_c1,phi_c2,sigma);


        W2_value_previous=calculate_dual_value(helper_f,phi,psi,mu);
        error_mu=calculate_h_minus_1(fftps,push_mu,helper_f.DEstar);

        for(int i=0;i<n1*n2;++i){
            phi[i] += fftps->workspace[i];
        }


        flt2d->find_c_concave(phi,psi,tau);

        W2_value=calculate_dual_value(helper_f,phi,psi,mu);

        sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu);


        return error_mu;
    }

// display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);

    void display_iteration(const int iter,const double W2_value,const double error_mu,const double error_nu,const double solution_error, const double C_phi, const double C_psi) const{
        cout << setprecision(6);
        cout << fixed;
        cout <<setw(5)<<iter+1 << " C : " << C_phi << " " << C_psi << " coeff : " << setw(10) << phi_c2/phi_c1 << " " << setw(6) << psi_c2/psi_c1 << "  W2 : " << scientific << setw(13) << W2_value << "  h-1 : "<<scientific<<setw(13) << error_mu << " " << error_nu << " solution_error : " << solution_error <<endl;
    }

    /**
        Calculate a = 0.1 * max(-phi)
    */
    double calculate_lambda() const{
        double phimax = 0;
        
        for(int i=0;i<n1*n2;++i){
            phimax = fmax(phimax,-phi[i]);
        }
        if(phimax > 5) return 2;
        return phimax * 0.2;
    }

    /**
        Calculate infgradphi = inf(|nabla phi|)
    */
    double calculate_infgradphi_on_level_set(const double lambda){

        double infgradphi= 1000;
        int count = 0;
        for(int i=0;i<n2-1;++i){
            for(int j=0;j<n1-1;++j){
                // if(fabs(fabs(phi[i*n1+j]/lambda) - 1)  < 1e-2){
                if(-phi[i] > 0 && -phi[i] < lambda){
                    double gradxphi = 1.0*n1*(phi[i*n1+j+1]-phi[i*n1+j]);
                    double gradyphi = 1.0*n2*(phi[(i+1)*n1+j]-phi[i*n1+j]);
                    // double gradxphi = calculate_gradx_vx(phi,i,j);
                    // double gradyphi = calculate_grady_vy(phi,i,j);
                    double eval = gradxphi*gradxphi + gradyphi*gradyphi;
                    infgradphi = fmin(infgradphi, eval);
                }
            }
        }

        return fmax(1,sqrt(infgradphi));
    }

    /** 
        calculate c1 and c2 for coefficients for the operator c1 Id + c2 Laplacian
        C is constant for the trace theorem
    */

    void calculate_D_coeff_psi(Helper_E& helper_f, double& c1, double& c2, const double C, const bool verbose, const double* mu, const double lambda, const double infgradphi, const double mu_max){

        double mprime = m/(m-1);
        double lambda_mprime = pow(lambda, mprime-1);
        double D2 = lambda_mprime/(mprime*infgradphi);
        double D3 = lambda_mprime/lambda;
        c1 = D2 + D3/C_phi;
        c2 = D2 + tau*mu_max;

        if(verbose){
            cout << "D2 : " << D2 << " D3/C : " << D3/C << "\n";
            cout << fixed << "PSI - lambda : " << lambda << " C : " << C << " infgradphi : " << infgradphi << " c1 : " << c1 << " c2 : " << c2 << " c2/c1 : " << c2/c1 << endl;   
        }
    }

    void set_coeff(double& c1, double& c2, const double C, const double mu_max, const double d1, const double d2, const bool verbose){
        c1 = C * d1 + d2;
        c2 = C * d1 + tau * mu_max;
    }

    void initialize_phi(Helper_E& helper_f,const double* mu){
        for(int i=0;i<n1*n2;++i){
            phi[i] = - (gamma * pow(mu[i],m-1) + helper_f.nu[i]);    
        }
    }

    double compute_barenblatt_solution_error(Helper_E& helper_f, Barenblatt& solution, const double* phi, const int outer_iter){
        // Compare with actual solution
        solution.calc_solution_at_n(outer_iter+1);

        helper_f.calculate_DEstar_normalized(phi);

        double sum=0;
        for(int i=0;i<n1*n2;++i){
            sum += fabs(solution.sol[i] - helper_f.DEstar[i]);
        }

        return sum/(1.0*n1*n2);
    }

    void calculate_d1_d2(double& d1, double& d2, const double lambda, const double infgradphi){
        double area = 0;
        for(int i=0;i<n1*n2;++i){
            if(-phi[i] > lambda){
                area += 1;
            }
        }
        area /= n1*n2;
        double eval = pow(lambda / gamma, mprime - 1);
        d1 = eval / infgradphi;
        d2 = eval / lambda * (mprime - 1);
    }

    double start_OT(Helper_E& helper_f, const double* mu, Barenblatt& solution, const int outer_iter, Initializer& init){

        int skip = 10; // frequency of printout

        double error_mu,error_nu,error=1.0;

        /*
            Initialize coefficients for siga update
        */

        beta_1 =0.1;
        beta_2 =0.9;
        alpha_1=1.1;
        alpha_2=0.9;

        /*
            Initialize the tolerance based on tau^2
        */

        double mu_max = 1;
        for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu[i]);
        double tol_modified = tolerance * mu_max *tau*tau;

        cout << "Iter : " << outer_iter + 1 << " Tolerance : " << tol_modified << "\n";

        /*
            Initialize the coefficients for fftps
        */

        double lambda = 1;
        double infgradphi  = 1;
        double sigma_forth = 1;
        double sigma_back  = 1;

        double d1 = 1;
        double d2 = 1;

        
        

        if(outer_iter==0){
            initialize_phi(helper_f,mu); // intiailize phi in the first outer iteration
            phi_c1 = 100;
            psi_c1 = 100;

            phi_c2 = 1;
            psi_c2 = 1;
        }

        double solution_error = 1;
        
        /*
            Starting the loop
        */
    

        for(int iter=0;iter<max_iteration;++iter){            
            
            /* 
                Calculating the relative error
            */
                        
            error=fmin(error_mu,error_nu);

            if(iter % 10 == 0){
                lambda = calculate_lambda();
                infgradphi = calculate_infgradphi_on_level_set(lambda);
                calculate_d1_d2(d1, d2, lambda, infgradphi);

                // C_phi = fmax(0.15,fmin(0.3,C_phi/sigma_forth));
                // C_psi = fmax(0.15,fmin(0.3,C_psi/sigma_back));
                set_coeff(phi_c1, phi_c2, C_phi, mu_max, d1, d2, false);
                set_coeff(psi_c1, psi_c2, C_psi, mu_max, d1, d2, false);
            }

            /*
                Determinant version pushforward
            */

            sigma_forth = 1;
            sigma_back  = 1;
            error_mu = perform_OT_iteration_forth_det(helper_f,sigma_forth,W2_value,mu);
            error_nu = perform_OT_iteration_back_det(helper_f,sigma_back,W2_value,mu);

            /*
                Display the result per iterations
            */

            cout<<"-"<<flush;

            if(iter%skip==skip-1){
                cout<<"|";
                /* Compare with actual solution */
                // solution_error = compute_barenblatt_solution_error(helper_f, solution, phi, outer_iter);
                display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);
                cout << "infgradphi : " << infgradphi << " c1 : " << phi_c1 << " " << psi_c1 << "\n";

                // string figurename = "output";
                // for(int i=0;i<n1*n2;++i) push_mu[i] = fabs(push_mu[i] - mu[i]);
                // init.save_image_opencv(push_mu,figurename,(iter+1)/skip, mu_max);
            }

            /* 
                Stopping Condition 
            */

            // if(W2_value>0 && ((abs(error)<tolerance && iter>=0) || iter==max_iteration-1 || sigma_forth <1e-9) ){
            if(((abs(error)<tol_modified && abs(error)>0 && iter>=1) || iter==max_iteration-1) ){
                cout<<"Tolerance met!"<<endl;
                /* Compare with actual solution */
                solution_error = compute_barenblatt_solution_error(helper_f, solution, phi, outer_iter);
                display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);
                break;
            }
        }

        return solution_error;
    }


}; // Back and Forth

int main(int argc, char** argv){
    if(argc!=9){
        cout<<"Do the following:"<<endl;
        cout<<"./gf [n1] [n2] [max_iteration] [tolerance] [nt] [tau] [m] [C trace]"<<endl;
        return 0;
    }

    int n1=stoi(argv[1]);
    int n2=stoi(argv[2]);
    int max_iteration=stoi(argv[3]);
    double tolerance=stod(argv[4]);
    int nt=stoi(argv[5]);
    double tau=stod(argv[6]);
    double m=stod(argv[7]);
    double C=stod(argv[8]);

    double M = 1.0; // initial mass

    Barenblatt solution(n1,n2,tau,m,M);
    solution.calc_solution_at_n(0);

    // Setting gamma value to correspond Barenblatt solution
    double gamma = m/(m-1);

    /* Initialize Initializer */
    Initializer init(n1,n2);

    create_csv_parameters(n1,n2,nt,tau,gamma,m,M);

    // Initialize mu
    double* mu=new double[n1*n2];
    memcpy(mu, solution.sol, n1*n2*sizeof(double));

    cout << "XXX Starting Gradient Flow XXX" << endl;

    cout << endl;
    // declaring argument of time() 
    time_t my_time = time(NULL);
  
    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));
    cout << endl;

    cout << "n1   : " << n1 <<endl;
    cout << "n2   : " << n2 <<endl;
    cout << "nt   : " << nt <<endl;
    cout << "tau  : " << tau << endl;
    cout << "s    : " << solution.s << endl; // Barenblatt s
    cout << "gamma: " << gamma << endl;
    cout << "m    : " << m << endl;
    cout << "Max Iteration : " << max_iteration <<endl;
    cout << "tolerance     : " << tolerance <<endl;

    Helper_E helper_f(n1,n2,gamma,tau,m,M);

    cout << setprecision(6);

    BackAndForth bf(n1,n2,max_iteration,tolerance,gamma,tau,m);

    bf.C_phi = C;
    bf.C_psi = C;

    string filename="./data/mu-"+to_string(0)+".csv";
    create_bin_file(mu,n1*n2,filename);

    string figurename = "barenblatt";
    // init.save_image_opencv(mu,figurename,0);

    clock_t time;
    time=clock();

    double sum, sum_solution=0;

    for(int n=0;n<nt;++n){

        double solution_error = bf.start_OT(helper_f, mu, solution, n, init);
        // helper_f.calculate_DEstar(bf.phi);
        helper_f.calculate_DEstar_normalized(bf.phi);
        memcpy(mu,helper_f.DEstar,n1*n2*sizeof(double));

        sum=0;
        for(int i=0;i<n1*n2;++i){
            sum+=mu[i];
        }
        cout << "sum : " << sum/(n1*n2) << endl;

        filename="./data/mu-"+to_string(n+1)+".csv";
        create_bin_file(mu,n1*n2,filename);

        sum_solution += solution_error;

        // init.save_image_opencv(mu,figurename,n+1);
    }

    time=clock()-time;
    printf ("\nCPU time for GF: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    cout << fixed;
    cout << "L1 Error : " << setprecision(10) << sum_solution / nt  <<endl;

    delete[] mu;
}