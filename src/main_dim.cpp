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
#include "Points.h"
#include "Sinkhorn.h"

using namespace std;

int main(int argc, char** argv){
    if(argc!=5){
        cout<<"Do the following:"<<"\n";
        cout<< argv[0] << " [DIM] [max_iteration] [tolerance] [lambda]"<<"\n";
        return 0;
    }

    int DIM=stoi(argv[1]);
    int max_iteration=stoi(argv[2]);
    double tolerance=stod(argv[3]);
    double lambda=stod(argv[4]);

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
    int n_mu = 2000;
    int n_nu = 2000;

    srand(1);

string filename = "error_N_2000.dat";
ofstream outfile_error;
outfile_error.open(filename);

filename = "time_check_N_2000.dat";
ofstream outfile_time;
outfile_time.open(filename);

int num_x = 10; // number of items in DIM_list
int DIM_list[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024}; // Dimensions list
for(int idx_DIM=0;idx_DIM<num_x;++idx_DIM){
    DIM = DIM_list[idx_DIM];
    // run the simulation 50 times per dimension
    for(int idx_multiple=0;idx_multiple<30;++idx_multiple){

        if(idx_multiple != 0){
            outfile_error << ",";
            outfile_time  << ",";
        }

    Points* mu = new Points(DIM, n_mu);
    Points* nu = new Points(DIM, n_nu);

    for(int p=0;p<n_mu;++p){
        for(int i=0;i<DIM;++i){
            (*mu)(p,i) = 1.0/8.0 * (2.0*rand()/RAND_MAX-1.0) + 1.0/4.0; 
            (*nu)(p,i) = 1.0/4.0 * (2.0*rand()/RAND_MAX-1.0) + 3.0/4.0; 
        }
    }

    cout << "XXX Starting particle BFM XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);

    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));

    /* Initialize BFM */
    Sinkhorn bf(DIM, n_mu, n_nu, max_iteration, tolerance, lambda);

    string figurename = "barenblatt";

    clock_t time;
    time=clock();

    bf.start_OT(mu, nu, plt);

    time=clock()-time;
    double elapsed_time = ((float)time)/CLOCKS_PER_SEC;
    printf ("\nCPU time for particle BFM: %f seconds.\n\n", elapsed_time);

    double calculated_W2_value = fmax(bf.dual_back_, bf.dual_forth_);
    double eval = 0.5* DIM * 49./192.; // actual value
    double L1_error = fabs(calculated_W2_value - eval);

    printf ("L1 error: %f BFM: %f Actual: %f\n\n", L1_error, calculated_W2_value, eval);

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