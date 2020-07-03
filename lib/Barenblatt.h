#ifndef BARENBLATT_H
#define BARENBLATT_H

#include <iostream>

class Barenblatt{
public:
    int n1;
    int n2;
    double tau;
    double m;
    double M;
    double mprime;
    double d;
    double alpha;
    double k;
    double xi0;
    double C;
    double s;
    double h0;

    double* sol;

    Barenblatt(){
        sol = NULL;
    }
    Barenblatt(int n1, int n2, double tau, double m, double M){

        sol = new double[n1*n2];

        this->n1=n1;
        this->n2=n2;
        this->tau=tau;
        this->m=m;
        this->M=M;

        mprime = m/(m-1); // conjugate
        d = 2; // 2 dimision

        alpha = d/(d*(m-1) + 2);
        k     = (m-1) * alpha/(2*m*d);
        C     = pow(M * mprime * k / M_PI, 1/mprime);


        // xi0 as a free parameter

        // xi0   = 0.2;
        // s     = pow(xi0*xi0*k/C,1.0/alpha);

        // s as a free parameter

        // s     = 1e-4;
        // xi0   = pow((C * pow(s,alpha)/k),0.5);

        // initial height as a free parameter

        h0 = 15;
        s  = pow(h0 * pow(C,-1.0/(m-1)), -1.0/alpha);
    }

    ~Barenblatt(){
        delete[] sol;
    }

    double F(double x){
        return pow(fmax(0, C - k*x*x), 1.0/(m-1));
    }

    void calc_solution_at_n(double n){
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double x= (j+0.5)/n1 - 0.5;
                double y= (i+0.5)/n2 - 0.5;

                double r = sqrt(x*x+y*y);

                double tval = n*tau + s;
                double eval = r * pow(tval, -alpha/d);
                sol[i*n1+j] = pow(tval,-alpha) * F(eval);
            }
        }
    }

}; // Barenblatt

#endif