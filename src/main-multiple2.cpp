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

int main(int argc, char** argv){
    if(argc!=6){
        cout<<"Do the following:"<<"\n";
        cout<< argv[0] << " [DIM] [max_iteration] [tolerance] [sigma] [nt]"<<"\n";
        return 0;
    }

    int DIM=stoi(argv[1]);
    int max_iteration=stoi(argv[2]);
    double tolerance=stod(argv[3]);
    double sigma=stod(argv[4]);
    int nt=stoi(argv[5]);

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
    int n_mu = 100;
    int n_nu = 100;

    // srand(1);

bool bool_for_estimator = false; // false: plugin true: H
string filename = "time_check_" + to_string(sigma) + ".dat";
ofstream outfile;
outfile.open(filename);

int num_N_list = 9;
double N_list[] = {1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.3, 3.7};
for(int idx_N=0;idx_N<num_N_list;++idx_N){
    double power_N = N_list[idx_N];
    n_mu = pow(10, power_N);
    n_nu = n_mu;

    for(int idx_multiple=0;idx_multiple<30;++idx_multiple){

        if(idx_multiple != 0){
            outfile << ",";
        }

    Points* mu = new Points(DIM, n_mu);
    Points* nu = new Points(DIM, n_nu);

    for(int p=0;p<n_mu;++p){
        for(int i=0;i<DIM;++i){
            // (*mu)(p,i) = 0.2 * (1.0*rand()/RAND_MAX-0.5) + 0.5;    
            (*mu)(p,i) = 1.0/8.0 * (2.0*rand()/RAND_MAX-1.0) + 1.0/4.0; 
            (*nu)(p,i) = 1.0/4.0 * (2.0*rand()/RAND_MAX-1.0) + 3.0/4.0; 
        }
    }

    if(bool_for_estimator){
        // perturb the points to get \hat{H}
        std::default_random_engine generator;
        double mean_normal_distribution  = 0;
        double stdev_normal_distribution = 0.1;
        std::normal_distribution<double> distribution(mean_normal_distribution, stdev_normal_distribution);



        for(int p=0;p<n_mu;++p){
            for(int i=0;i<DIM;++i){
                double number = distribution(generator);
                (*mu)(p,i) += number;
                (*nu)(p,i) += number;
            }
        }
    }


    cout << "XXX Starting Sinkhorn XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);

    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));

    /* Initialize BFM */
    Sinkhorn bf(DIM, n_mu, n_nu, max_iteration, tolerance, sigma);
    bf.sigma_forth_ = sigma;
    bf.sigma_back_  = sigma;

    string figurename = "barenblatt";

    clock_t time;
    time=clock();

    bf.start_OT(mu, nu, plt);

    time=clock()-time;
    double elapsed_time = ((float)time)/CLOCKS_PER_SEC;
    printf ("\nCPU time for particle BFM: %f seconds.\n\n", elapsed_time);

    double calculated_W2_value = bf.dual_forth_;
    double L1_error = fabs(calculated_W2_value - 0.5* DIM * 49./192.);

    printf ("L1 error: %f\n", L1_error);    

    outfile << L1_error;

    // bf.create_interpolate_video(mu, nu, nt, plt);

    delete mu;
    delete nu;


    // TODO: interpolate video
} outfile << "\n";}

cout << "Herere\n";

outfile.close();

    

} // end of main