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
        if(nu[i*n1+j] < 0) return 0;

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
            if(nu[i] >= 0){
                DEstar[i] = 1.0/gammaprime * exp(exponent * log(fmax(0, - phi[i] - nu[i])));    
            }else{
                DEstar[i] = 0;
            }
        }
    }

    void calculate_DEstar_normalized(double* phi, const double tolerance=1e-7){
        
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
                if(nu[i] >= 0){
                    double eval=- phi[i] - nu[i] + lambda;
                    if(eval>0){
                        sum+=exp(exponent*log(eval)); // IMPORTANT    
                    }    
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
            if(nu[i] >= 0){
                DEstar[i] = 1.0/gammaprime * pow(fmax(0, - phi[i] - nu[i] + lambda), 1.0/(m-1));    
            }else{
                DEstar[i] = 0;
            }
        }
        for(int i=0;i<n1*n2;++i) phi[i] -= lambda;
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

    double* vxx;
    double* vyy;
    double* vxy;

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

        vxx=new double[n1*n2];
        vyy=new double[n1*n2];
        vxy=new double[n1*n2];

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
        delete[] vxx;
        delete[] vyy;
        delete[] vxy;
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

        /* forward difference */
        // for(int i=0;i<n2;++i){
        //     for(int j=0;j<n1-1;++j){
        //         vx[i*n1+j]=1.0*n1*(phi_c[i*n1+j+1]-phi_c[i*n1+j]);
        //     }
        // }
        // for(int j=0;j<n1;++j){
        //     for(int i=0;i<n2-1;++i){
        //         vy[i*n1+j]=1.0*n2*(phi_c[(i+1)*n1+j]-phi_c[i*n1+j]);
        //     }
        // }

        /* centered difference*/
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                vx[i*n1+j]=0.5*n1*(phi_c[i*n1+(int)fmin(n1-1,j+1)]-phi_c[i*n1+(int)fmax(0,j-1)]);
                vy[i*n1+j]=0.5*n2*(phi_c[(int)fmin(n2-1,i+1)*n1+j]-phi_c[(int)fmax(0,i-1)*n1+j]);
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

    double calculate_L1_error(const double* push_mu, const double* mu){
        double error = 0;
        for(int i=0;i<n1*n2;++i) error += fabs(push_mu[i] - mu[i]);
        return error/(1.0*n1*n2);
    }

    void calculate_push_rho(const double* rho, double* push_rho,const double* vx,const double* vy,const double* vxx,const double* vyy,const double* vxy){

        double eps = pow(1.0/n1, 0.7);

        double xpost,ypost,xpre,ypre;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=vx[i*n1+j];
                double vyval=vy[i*n1+j];

                double x=(j+0.5)/(1.0*n1)-tau*vxval;
                double y=(i+0.5)/(1.0*n2)-tau*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){

                    double vxx_val = interpolate_function(x,y,vxx);
                    double vyy_val = interpolate_function(x,y,vyy);
                    double vxy_val = interpolate_function(x,y,vxy);

                    // push_rho[i*n1+j]=rhovalue/fabs((1.0-tau*gradx_vx) * (1.0-tau*grady_vy)); 
                    double det = fabs((1.0-tau*vxx_val) * (1.0-tau*vyy_val) - tau*tau * vxy_val * vxy_val);
                    // double det = fabs((1.0-tau*gradx_vx) * (1.0-tau*grady_vy));
                    det = fmax(eps,det);
                    push_rho[i*n1+j] = rhovalue/det;


                ////////////////////////////////////
                    // double vxval_post=vx[i*n1+(int)fmin(n1-1,j+1)];
                    // double vyval_post=vy[(int)fmin(n2-1,i+1)*n1+j];

                    // double gradx_vx=1.0*n1*(vxval_post-vxval);
                    // double grady_vy=1.0*n2*(vyval_post-vyval);

                    // double gradx_vx=calculate_gradx_WENO_(vx, i, j);
                    // double grady_vy=calculate_grady_WENO_(vy, i, j);
                ////////////////////////////////////

                    
                }else{
                    push_rho[i*n1+j]=0;
                }

            }
        }
    }

    void calculate_pull_rho(const double* rho, double* push_rho,const double* vx,const double* vy,const double* vxx,const double* vyy,const double* vxy){

        double eps = pow(1.0/n1, 0.7);

        double xpost,ypost,xpre,ypre;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=vx[i*n1+j];
                double vyval=vy[i*n1+j];

                double x=(j+0.5)/(1.0*n1)-tau*vxval;
                double y=(i+0.5)/(1.0*n2)-tau*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){
                    double vxx_val = vxx[i*n1+j];
                    double vyy_val = vyy[i*n1+j];
                    double vxy_val = vxy[i*n1+j];

                    double det = fabs((1.0-tau*vxx_val) * (1.0-tau*vyy_val) - tau*tau * vxy_val * vxy_val);
                    det = fmin(1.0/eps, det);
                    push_rho[i*n1+j] = rhovalue*det;
                    
                }else{
                    push_rho[i*n1+j]=0;
                }

            }
        }
    }


    void calculate_push_rho(const double* rho, double* push_rho,const double* vx,const double* vy,const double* vxx,const double* vyy,const double* vxy,const double* phi){

        double eps = pow(1.0/n1, 0.7);
        double det_threshold = 0.99;

        double xpost,ypost,xpre,ypre;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double vxval=vx[i*n1+j];
                double vyval=vy[i*n1+j];

                double x=(j+0.5)/(1.0*n1)-tau*vxval;
                double y=(i+0.5)/(1.0*n2)-tau*vyval;

                double rhovalue=interpolate_function(x,y,rho);

                if(rhovalue>0){

                    double vxx_val = interpolate_function(x,y,vxx);
                    double vyy_val = interpolate_function(x,y,vyy);
                    double vxy_val = interpolate_function(x,y,vxy);

                    double det = fabs((1.0-tau*vxx_val) * (1.0-tau*vyy_val) - tau*tau * vxy_val * vxy_val);

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

                        det = fabs((1.0-tau*vxx_val) * (1.0-tau*vyy_val) - tau*tau * vxy_val * vxy_val);
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
    double update_sigma(double sigma, const double W2_value, const double W2_value_previous, const double error){

        if(W2_value_previous-W2_value>-beta_1*error){
            return alpha_2;
        }else if(W2_value_previous-W2_value<-beta_2*error){
            return alpha_1;
        }
        return 1;
    }

    void calculate_gradient_vxx_vyy_vxy(const double* phi, double* vxx, double* vyy, double* vxy){
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

    double perform_OT_iteration_back_det(Helper_E& helper_f,double& sigma,double& W2_value,const double* mu, const int iter){
        // ------------------------------------------------------------
        // double W2_value_previous = 0;
        double error_nu = 0;

        flt2d->find_c_concave(psi,phi,tau);
        flt2d->find_c_concave(phi,psi,tau);
    
        helper_f.calculate_DEstar_normalized(phi);    

        calculate_gradient(psi, vx, vy);

        calculate_gradient_vxx_vyy_vxy(phi, vxx, vyy, vxy);
        calculate_push_rho(helper_f.DEstar, push_mu,vx,vy,vxx,vyy,vxy,psi);

        // calculate_gradient_vxx_vyy_vxy(psi, vxx, vyy, vxy);
        // calculate_pull_rho(helper_f.DEstar, push_mu,vx,vy,vxx,vyy,vxy);       

        fftps->perform_inverse_laplacian(push_mu,mu,psi_c1,psi_c2,sigma);

        // W2_value_previous=calculate_dual_value(helper_f,phi,psi,mu);
        // error_nu=calculate_h_minus_1(fftps,push_mu,mu);
        error_nu=calculate_L1_error(push_mu,mu);

        for(int i=0;i<n1*n2;++i){
            psi[i] += fftps->workspace[i];
        }

        flt2d->find_c_concave(psi,phi,tau);

        W2_value=calculate_dual_value(helper_f,phi,psi,mu);

        // sigma = update_sigma(sigma, W2_value, W2_value_previous, error_nu);

        return error_nu;
    }


    double perform_OT_iteration_forth_det(Helper_E& helper_f,double& sigma,double& W2_value,const double* mu, const int iter){
        // ------------------------------------------------------------
        // double W2_value_previous = 0;
        double error_mu = 0;

        flt2d->find_c_concave(phi,psi,tau);
        flt2d->find_c_concave(psi,phi,tau);
            
        calculate_gradient(phi, vx, vy);

        calculate_gradient_vxx_vyy_vxy(psi, vxx, vyy, vxy);
        calculate_push_rho(mu, push_mu,vx,vy,vxx,vyy,vxy,phi);

        // calculate_gradient_vxx_vyy_vxy(phi, vxx, vyy, vxy);
        // calculate_pull_rho(mu, push_mu,vx,vy,vxx,vyy,vxy);    
            
        helper_f.calculate_DEstar_normalized(phi);    

        fftps->perform_inverse_laplacian(push_mu,helper_f.DEstar,phi_c1,phi_c2,sigma);


        // W2_value_previous=calculate_dual_value(helper_f,phi,psi,mu);
        // error_mu=calculate_h_minus_1(fftps,push_mu,helper_f.DEstar);
        error_mu=calculate_L1_error(push_mu,helper_f.DEstar);

        for(int i=0;i<n1*n2;++i){
            phi[i] += fftps->workspace[i];
        }


        flt2d->find_c_concave(phi,psi,tau);

        W2_value=calculate_dual_value(helper_f,phi,psi,mu);

        // sigma = update_sigma(sigma, W2_value, W2_value_previous, error_mu);

        return error_mu;
    }

// display_iteration(iter,W2_value,error_mu,error_nu,solution_error,C_phi,C_psi);

    void display_iteration(const int iter,const double W2_value,const double error_mu,const double error_nu, const double C_phi, const double C_psi, const double infgradphi) const{
        cout << setprecision(6);
        cout << fixed;
        cout <<setw(5)<<iter+1 << " C : " << C_phi << "  infgradphi : " << infgradphi << "  c1 : " << phi_c1 << " " << psi_c1 << " coeff : " << setw(10) << phi_c2/phi_c1 << " " << setw(6) << psi_c2/psi_c1 << "  W2 : " << scientific << setw(13) << W2_value << "  L1 error : "<<scientific<<setw(13) << error_mu << " " << error_nu<<"\n";
    }

    /**
        Calculate a = 0.1 * max(-phi)
    */
    double calculate_lambda(const double* nu) const{
        double phimax = 1;
        
        for(int i=0;i<n1*n2;++i){
            if(nu[i] >= 0){
                phimax = fmax(phimax,-phi[i]-nu[i]);    
            }
        }
        // if(phimax > 2) return 0.1;
        // return fmin(2,phimax * 0.5);

        return fmax(0.01,phimax * 0.001);
    }

    /**
        Calculate infgradphi = inf(|nabla phi|)
    */
    double calculate_infgradphi_on_level_set(const double lambda, const double* nu){

        double infgradphi= 10000;
        int count = 0;
        for(int i=1;i<n2-1;++i){
            for(int j=1;j<n1-1;++j){
                // if(fabs(fabs(phi[i*n1+j]/lambda) - 1)  < 1e-2){
                if(nu[i*n1+j] >= 0 && nu[i*n1+j+1] >= 0 && nu[(i+1)*n1+j] >= 0  && nu[i*n1+j-1] >= 0 && nu[(i-1)*n1+j] >= 0){
                    if(-phi[i*n1+j]-nu[i*n1+j] > 0 && -phi[i*n1+j]-nu[i*n1+j] < lambda){
                        double gradxphi = 0.5*n1*(phi[i*n1+j+1]-phi[i*n1+j-1]-nu[i*n1+j+1]+nu[i*n1+j-1]);
                        double gradyphi = 0.5*n2*(phi[(i+1)*n1+j]-phi[(i-1)*n1+j]-nu[(i+1)*n1+j]+nu[(i-1)*n1+j]);
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

        return fmax(0.1,sqrt(infgradphi));
    }

    void set_coeff_m_2(double& c1, double& c2, const double mu_max, const double C_c_transform){
        c1 = 1/gamma;
        c2 = C_c_transform * tau;
    }

    void set_coeff(double& c1, double& c2, const double C, const double mu_max, const double d1, const double d2, const double d3, const bool verbose, const double C_c_transform){
        // c1 = C * d1 + d2;
        // c2 = C * d1 + C_c_transform * tau;

        c1 = d3 * C * d1 + d2;
        c2 = d3 * C * d1 + C_c_transform * tau;
    }

    void initialize_phi(Helper_E& helper_f,const double* mu){
        for(int i=0;i<n1*n2;++i){
            if(helper_f.nu[i] >= 0){
                phi[i] = - (gamma * pow(mu[i],m-1) + helper_f.nu[i]);
            } else {
                phi[i] = 0;
            }
        }
    }

    void calculate_d1_d2_d3(double& d1, double& d2, double& d3, const double lambda, const double infgradphi, const double* nu){
        // double area = 0;
        // for(int i=0;i<n1*n2;++i){
        //     if(nu[i] >= 0){
        //         if(-phi[i]-nu[i] > lambda){
        //             area += 1;
        //         }    
        //     }
        // }
        // area /= n1*n2;
        
        double eval = 1.0 / pow(gamma, mprime - 1);

        d1 = eval * pow(lambda, mprime-1) / infgradphi;
        d2 = eval * pow(lambda, mprime-2) * (mprime - 1);

        // cout << "d1 : " << d1 << " d2 : " << d2 << "\n";
    }

    void calculate_c_transform_constant(double& C, const double* phi, const double* mu){

        C = 0;

        // calculate_gradient(phi, vx, vy);

        double mu_max = 0; for(int i=0;i<n1*n2;++i) mu_max = fmax(mu_max, mu[i]);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double x = (j+0.5)/n1 - tau * vx[i*n1+j];
                double y = (i+0.5)/n2 - tau * vy[i*n1+j];
                // double mu_val = interpolate_function(x,y,mu);
                double mu_val = 1;

                if(mu_val > 0){
                    /* calculate eigen values */

                    /* calculate trace */
                    double trace = fabs(2 - tau*vxx[i*n1+j] - tau*vyy[i*n1+j]);
                    /* calculate det */
                    double det   = fabs((1.0-tau*vxx[i*n1+j]) * (1.0-tau*vyy[i*n1+j]) - tau*tau * vxy[i*n1+j] * vxy[i*n1+j]);


                    double t1 = 0.5 * fabs(tau + sqrt(fabs(tau*tau - 4 * det)));
                    double t2 = 0.5 * fabs(tau - sqrt(fabs(tau*tau - 4 * det)));

                    
                    C = fmax(C, mu_val*fmax(t1,t2));
                }
            }
        }

        C *= mu_max;
    }

    void start_OT(Helper_E& helper_f, const double* mu, const int outer_iter, Initializer& init){

        const int skip = 50; // frequency of printout

        double error_mu = 1.0;
        double error_nu = 1.0;
        double error=1.0;

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
        // double tol_modified = tolerance * mu_max *tau*tau;
        double tol_modified = tolerance;

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
        double d3 = 1;

        double C_c_transform = 1;

        
        
        double max_iteration_tmp = max_iteration;

        initialize_phi(helper_f,mu); // intiailize phi in the first outer iteration
        
        /*
            Starting the loop
        */
        
        for(int iter = 0; iter < max_iteration_tmp; ++iter){
            
            /* 
                Calculating the relative error
            */

            if(iter % 1 == 0 && iter >= 0){
                if(m == 2){
                    calculate_c_transform_constant(C_c_transform, phi, mu);
                    set_coeff_m_2(phi_c1, phi_c2, mu_max, C_c_transform);
                    calculate_c_transform_constant(C_c_transform, psi, helper_f.DEstar);
                    set_coeff_m_2(psi_c1, psi_c2, mu_max, C_c_transform);    
                }else{
                    lambda = calculate_lambda(helper_f.nu);
                    infgradphi = calculate_infgradphi_on_level_set(lambda,helper_f.nu);
                    calculate_d1_d2_d3(d1, d2, d3, lambda, infgradphi, helper_f.nu);
                    

                    // C_phi = fmax(0.15,fmin(0.3,C_phi/sigma_forth));
                    // C_psi = fmax(0.15,fmin(0.3,C_psi/sigma_back));

                    d3 = pow(1.1,(iter+1)/100);
                    calculate_c_transform_constant(C_c_transform, phi, mu);
                    set_coeff(phi_c1, phi_c2, C_phi, mu_max, d1, d2, d3, false, C_c_transform);
                    calculate_c_transform_constant(C_c_transform, psi, helper_f.DEstar);
                    set_coeff(psi_c1, psi_c2, C_psi, mu_max, d1, d2, d3, false, C_c_transform);    
                }
            }

            if(outer_iter==0 && iter < 10){
            
                phi_c1 = 0.1;
                psi_c1 = 0.1;

                phi_c2 = 0.1;
                psi_c2 = 0.1;            
            }

            /*
                Determinant version pushforward
            */

            sigma_forth = 1;
            sigma_back  = 1;
            
            error_mu = perform_OT_iteration_forth_det(helper_f,sigma_forth,W2_value,mu,iter);
            error_nu = perform_OT_iteration_back_det(helper_f,sigma_back,W2_value,mu,iter);
            
            error=fmin(error_mu,error_nu);

            /* 
                Stopping Condition 
            */

            // if(W2_value>0 && ((abs(error)<tolerance && iter>=0) || iter==max_iteration-1 || sigma_forth <1e-9) ){
            if(((abs(error)<tol_modified && abs(error)>0 && iter>=0) || iter==max_iteration_tmp-1) ){
                // cout << "error : " << error << " iter : " << iter << "\n";
                // cout<<"Tolerance met!"<<"\n";
                display_iteration(iter,W2_value,error_mu,error_nu,C_phi,C_psi,infgradphi);
                break;
            }

            /*
                Display the result per iterations
            */

            // cout<<"-"<<flush;

            if(iter%skip==skip-1){
                // cout<<"|";
                display_iteration(iter,W2_value,error_mu,error_nu,C_phi,C_psi,infgradphi);

                string figurename = "output";
                for(int i=0;i<n1*n2;++i) push_mu[i] = fabs(push_mu[i] - mu[i]);
                init.save_image_opencv(push_mu,figurename,(iter+1)/skip, mu_max);
                // init.save_image_opencv(helper_f.DEstar,figurename,(iter+1)/skip, mu_max);
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

    double M = 1;


    cout << "n1   : " << n1 <<"\n";
    cout << "n2   : " << n2 <<"\n";
    cout << "nt   : " << nt <<"\n";
    cout << "tau  : " << tau << "\n";
    cout << "gamma: " << gamma << "\n";
    cout << "m    : " << m << "\n";
    cout << "Max Iteration : " << max_iteration <<"\n";
    cout << "tolerance     : " << tolerance <<"\n";


    /* modify gamma */
    gamma = gamma * m /(m-1);

    cout << "Modified gamma : " << gamma << "\n";

    /* Initialize Initializer */
    Initializer init(n1,n2);

    create_csv_parameters(n1,n2,nt,tau,gamma,m,M);

    // Initialize mu
    double* mu=new double[n1*n2];
    // create_mu_square(mu,0.2,0.2,0.1,n1,n2);
    create_mu_from_image(mu,n1,n2);

    cout << "XXX Starting Gradient Flow XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);
  
    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));
    cout << "\n";

    // unsigned char* obstacle = new unsigned char[n1*n2];

    string data_folder = "data";

    Helper_E helper_f(n1,n2,gamma,tau,m,M);

    double a = 0.1;
    init_entropy_sine(helper_f.nu, 5, 3, a, n1, n2);

    cout << "a for entropy sine : " << a << "\n";


    // init_obstacle_from_image(obstacle, n1, n2);
    // init_obstacle_two_moons(obstacle, n1, n2);
    // init_obstacle_pac_man(obstacle, n1, n2);
    // init_obstacle_circle(obstacle, n1, n2);

    // init.init_entropy_image_obstacle_opencv(helper_f.nu, obstacle, data_folder);

    cout << setprecision(6);

    BackAndForth bf(n1,n2,max_iteration,tolerance,gamma,tau,m);

    bf.C_phi = C;
    bf.C_psi = C;

    string filename="./data/mu-"+to_string(0)+".csv";
    create_bin_file(mu,n1*n2,filename);

    string figurename = "barenblatt";

    // if(plot > 0) init.save_image_opencv(mu,obstacle,figurename,0);
    if(plot > 0) init.save_image_opencv(mu,figurename,0);

    clock_t time;
    time=clock();

    double sum = 0;

    for(int n=0;n<nt;++n){

        bf.start_OT(helper_f, mu, n, init);
        // helper_f.calculate_DEstar(bf.phi);
        helper_f.calculate_DEstar_normalized(bf.phi);
        memcpy(mu,helper_f.DEstar,n1*n2*sizeof(double));

        /* print out sum */
        sum=0; for(int i=0;i<n1*n2;++i) sum+=mu[i]; cout << "sum : " << sum/(n1*n2) << "\n";

        filename="./data/mu-"+to_string(n+1)+".csv";
        create_bin_file(mu,n1*n2,filename);

        if(plot > 0) init.save_image_opencv(mu,figurename,n+1);
    }

    time=clock()-time;
    printf ("\nCPU time for GF: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    delete[] mu;
}