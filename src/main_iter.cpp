#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include <random>
#include "Accessory.h"
#include "Sinkhorn.h"

using namespace std;

double calculate_actual_answer(int DIM, double muA, double muB, double sigmaA, double sigmaB, double lambda){
    double cons = sqrt(4*sigmaA*sigmaA*sigmaB*sigmaB + lambda*lambda*lambda*lambda);
    double eval = sigmaA*sigmaA + sigmaB*sigmaB - cons + lambda*lambda*(1-log(2*lambda*lambda)) + lambda*lambda*log(cons+lambda*lambda);
    return 0.5*DIM*((muA-muB)*(muA-muB) + eval);
}

int main(int argc, char** argv){
    if(argc!=6){
        cout<<"Do the following:"<<"\n";
        cout<< argv[0] << " [DIM] [max_iteration] [tolerance] [lambda] [sigma]"<<"\n";
        return 0;
    }

    int DIM=stoi(argv[1]);
    int max_iteration=stoi(argv[2]);
    double tolerance=stod(argv[3]);
    double lambda=stod(argv[4]);
    double sigma=stod(argv[5]);

    if(DIM < 2){
        cout << "DIM needs to be greater or equal to 2\n";
        return -1;
    }

    int n1=64;
    int n2=64;

    string data_folder = "data";

    cout << "n1            : " << n1 <<"\n";
    cout << "n2            : " << n2 <<"\n";
    cout << "DIM           : " << DIM <<"\n";
    cout << "Max Iteration : " << max_iteration <<"\n";
    cout << "tolerance     : " << tolerance <<"\n";

    create_csv_parameters(n1,n2);

    /* Initialize mu and nu */
    int n_mu = 0;
    int n_nu = 0;

    srand(time(NULL));

    DIM   = 2;

    n_mu  = n1*n2;
    n_nu  = n1*n2;

    string filename = "data/original_error_DIM_8.dat";
    ofstream outfile;
    outfile.open(filename);

    double* pos_mu = new double[n_mu*DIM];
    double* pos_nu = new double[n_nu*DIM];
    double* mu = new double[n_mu];
    double* nu = new double[n_nu];

    double sigmaA = 0.1;
    double sigmaB = 0.05;

    double muA = 0.3;
    double muB = 0.7;

    // doing 2d example
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            // setup position of mu and nu
            int p = i*n1 + j;
            pos_mu[p*DIM+1] = (i+0.5)/n2;
            pos_mu[p*DIM+0] = (j+0.5)/n1;
            // setup mu and nu
            double x = (j+0.5)/n1;
            double y = (i+0.5)/n2;
            mu[i*n1+j] = 1.0/(2*M_PI*sigmaA*sigmaA) * exp(-0.5/(sigmaA*sigmaA) * (pow(x-muA,2)+pow(y-muA,2)));
            nu[i*n1+j] = 1.0/(2*M_PI*sigmaB*sigmaB) * exp(-0.5/(sigmaB*sigmaB) * (pow(x-muB,2)+pow(y-muB,2)));
        }
    }

    

    double actual_answer = calculate_actual_answer(DIM, muA, muB, sigmaA, sigmaB, lambda);

    cout << "XXX Starting particle BFM XXX" << "\n\n";

    // declaring argument of time()
    time_t my_time = time(NULL);

    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));

    // Initialize BFM
    Sinkhorn bf(DIM, n_mu, n_nu, max_iteration, tolerance, lambda, sigma);

    clock_t time;
    time=clock();

    bf.start_OT(mu, pos_mu, nu, pos_nu);

    time=clock()-time;
    double elapsed_time = ((float)time)/CLOCKS_PER_SEC;
    printf ("\nCPU time for particle BFM: %f seconds.\n\n", elapsed_time);

    double calculated_W2_value = bf.dual_value_;
    double L1_error = fabs(calculated_W2_value - actual_answer);

    printf ("L1 error: %f BFM: %f Actual: %f\n\n", L1_error, calculated_W2_value, actual_answer);

    // write dualvalue from iter=0,...,max_iter into outfile
    for(int i_=0;i_<max_iteration;++i_){
        // outfile << fabs(bf.dual_value_list_[i_] - actual_answer);
        outfile << bf.dual_value_list_[i_];
        if(i_<max_iteration-1){ outfile << ','; }
    }
    outfile << '\n';

    // delete mu;
    // delete nu;

    delete[] pos_mu;
    delete[] pos_nu;
    delete[] mu;
    delete[] nu;


    outfile.close();
} // end of main