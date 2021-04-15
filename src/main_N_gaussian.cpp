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
#include "Plotting.h"
#include "Points.h"
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

    int n1=1024;
    int n2=1024;

    if(DIM < 2){
        cout << "DIM needs to be greater or equal to 2\n";
        return -1;
    }

    

    string data_folder = "data";

    cout << "n1            : " << n1 <<"\n";
    cout << "n2            : " << n2 <<"\n";
    cout << "DIM           : " << DIM <<"\n";
    cout << "Max Iteration : " << max_iteration <<"\n";
    cout << "tolerance     : " << tolerance <<"\n";

    /* Initialize plotting tool */
    Plotting plt(n1,n2);

    create_csv_parameters(n1,n2);

    /* Initialize mu and nu */
    int n_mu = 0;
    int n_nu = 0;

    srand(time(NULL));

DIM = 32;
string filename = "data/gaussian_0.5_error_DIM_32.dat";
ofstream outfile_error;
outfile_error.open(filename);

filename = "data/gaussian_0.5_time_check_DIM_32.dat";
ofstream outfile_time;
outfile_time.open(filename);

int num_x = 5; // number of items in N_list
int N_list[] = {100, 200, 400, 800, 1600}; // Dimensions list
// double sigma_list[] = {0.01, 0.01, 0.01, 0.01, 0.005}; // log
double sigma_list[] = {0.005, 0.005, 0.005, 0.005, 0.001}; // gaussian 0.05


for(int idx_DIM=0;idx_DIM<num_x;++idx_DIM){
    sigma = sigma_list[idx_DIM];
    n_mu = N_list[idx_DIM];
    n_nu = n_mu;
    // run the simulation 50 times per dimension
    for(int idx_multiple=0;idx_multiple<40;++idx_multiple){

        if(idx_multiple != 0){
            outfile_error << ",";
            outfile_time  << ",";
        }

    Points* mu = new Points(DIM, n_mu);
    Points* nu = new Points(DIM, n_nu);

    double sigmaA = 0.05;
    double sigmaB = 0.05;

    double muA = 0.3;
    double muB = 0.7;

    { // define mu
        std::random_device rd;
        std::default_random_engine generator;
        generator.seed( rd() );
        std::normal_distribution<double> distribution(muA,sigmaA);
        for(int p=0;p<n_mu;++p){
            for(int i=0;i<DIM;++i){
                (*mu)(p,i) = distribution(generator);
            }
        }
    }

    // for(int i=0;i<n_mu;++i){
    //     printf("mu[%4d] = ", i);
    //     for(int j=0;j<DIM;++j) printf("%4.2f ", mu->data[i*DIM+j]);
    //     printf("\n");
    // }

    { // define nu
        std::random_device rd;
        std::default_random_engine generator;
        generator.seed( rd() );
        std::normal_distribution<double> distribution(muB,sigmaB);
        for(int p=0;p<n_mu;++p){
            for(int i=0;i<DIM;++i){
                (*nu)(p,i) = distribution(generator);
            }
        }
    }

    double actual_answer = calculate_actual_answer(DIM, muA, muB, sigmaA, sigmaB, lambda);

    cout << "XXX Starting particle BFM run: " << idx_multiple << " XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);

    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));

    /* Initialize BFM */
    Sinkhorn bf(DIM, n_mu, n_nu, max_iteration, tolerance, lambda, sigma);

    string figurename = "barenblatt";

    clock_t time;
    time=clock();

    bf.start_OT(mu, nu, plt);

    time=clock()-time;
    double elapsed_time = ((float)time)/CLOCKS_PER_SEC;
    printf ("\nCPU time for particle BFM: %f seconds.\n\n", elapsed_time);

    double calculated_W2_value = bf.dual_value_;
    double L1_error = fabs(calculated_W2_value - actual_answer);

    printf ("L1 error: %f BFM: %f Actual: %f\n\n", L1_error, calculated_W2_value, actual_answer);

    outfile_error << L1_error;
    outfile_time  << elapsed_time;

    delete mu;
    delete nu;
} outfile_error << "\n";
  outfile_time  << "\n";}

cout << "Herere\n";

outfile_error.close();
outfile_time .close();

    

} // end of main