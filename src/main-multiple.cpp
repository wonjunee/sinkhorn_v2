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
#include "BFM.h"

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

    srand(1);

string filename = "time_check.dat";
ofstream outfile;
outfile.open(filename);

int num_x = 5;
int DIM_list[] = {25, 50, 100, 250, 500};
for(int idx_DIM=0;idx_DIM<num_x;++idx_DIM){
    DIM = DIM_list[idx_DIM];
    for(int idx_multiple=0;idx_multiple<50;++idx_multiple){

        if(idx_multiple != 0){
            outfile << ",";
        }


    Points* mu = new Points(DIM, n_mu);
    Points* nu = new Points(DIM, n_nu);

    for(int p=0;p<n_mu;++p){
        for(int i=0;i<DIM;++i){
            // (*mu)(p,i) = 0.2 * (1.0*rand()/RAND_MAX-0.5) + 0.5;    
            (*mu)(p,i) = 0.3 * (1.0*rand()/RAND_MAX-0.5) + 0.3;     
            (*nu)(p,i) = 0.3 * (1.0*rand()/RAND_MAX-0.5) + 0.6;     
        }
    }

    for(int p=0;p<n_mu;++p){
        for(int i=0;i<DIM;++i){
            // (*mu)(p,i) = 0.2 * (1.0*rand()/RAND_MAX-0.5) + 0.5;    
            (*mu)(p,i) = 0.3 * (1.0*rand()/RAND_MAX-0.5) + 0.3;     
            (*nu)(p,i) = 0.3 * (1.0*rand()/RAND_MAX-0.5) + 0.6;     
        }
    }

    // for(int p=0;p<n_nu;++p){
    //     double r = 0.2*1.0*rand()/RAND_MAX + 0.3;
    //     double theta = 2.0*M_PI*rand()/RAND_MAX;

    //     (*nu)(p,0) = r * cos(theta) + 0.5;
    //     (*nu)(p,1) = r * sin(theta) + 0.5;
    //     for(int d=2;d<DIM;++d){
    //         (*nu)(p,d) = 0.2 * (1.0*rand()/RAND_MAX-0.5) + 0.5;
    //     }
    // }

    if(false){
        (*mu)(0,0) = 0.1;
        (*mu)(0,1) = 0.1;

        (*mu)(1,0) = 0.1;
        (*mu)(1,1) = 0.4;

        (*mu)(2,0) = 0.1;
        (*mu)(2,1) = 0.5;

        (*nu)(0,0) = 0.7;
        (*nu)(0,1) = 0.2;

        (*nu)(1,0) = 0.7;
        (*nu)(1,1) = 0.200001;

        (*nu)(2,0) = 0.7;
        (*nu)(2,1) = 0.200002;
    }

    double sigma_back = sigma;

    if(false){
        double minval = 1;
        for(int j=0;j<n_mu;++j){
            double mux = (*mu)(j,0);
            double muy = (*mu)(j,1);
            for(int j1=j+1;j1<n_mu;++j1){
                double mu1x = (*mu)(j1,0);
                double mu1y = (*mu)(j1,1);

                minval = fmin(minval, pow(mux-mu1x,2) + pow(muy-mu1y,2));
            }
        }    
        printf("mu min: %f  ", minval);

        double minval2 = 1;
        for(int i=0;i<n_nu;++i){
            double nux = (*nu)(i,0);
            double nuy = (*nu)(i,1);
            for(int i1=i+1;i1<n_nu;++i1){
                double nu1x = (*nu)(i1,0);
                double nu1y = (*nu)(i1,1);

                minval2 = fmin(minval2, pow(nux-nu1x,2) + pow(nuy-nu1y,2));
            }
        }    

        printf("nu min: %f\n", minval2);

        // double minval3 = 1;
        // for(int i=0;i<n_nu;++i){
        //     double nux = (*nu)(i,0);
        //     double nuy = (*nu)(i,1);
        //     for(int j=i+1;j<n_mu;++j){
        //         double nu1x = (*mu)(j,0);
        //         double nu1y = (*mu)(j,1);

        //         minval3 = fmin(minval3, pow(nux-nu1x,2) + pow(nuy-nu1y,2));
        //     }
        // }    

        // printf("nu min: %f\n", minval3);

        sigma = minval;
        sigma_back = minval2;

        // sigma = 0.3*fmin(minval3, sigma);
    }    

    cout << "XXX Starting particle BFM XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);

    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));

    /* Initialize BFM */
    BackAndForth bf(DIM, n_mu, n_nu, max_iteration, tolerance, sigma);
    bf.sigma_forth_ = sigma;
    bf.sigma_back_  = sigma_back;

    string figurename = "barenblatt";

    clock_t time;
    time=clock();

    bf.start_OT(mu, nu, plt);

    time=clock()-time;
    double elapsed_time = ((float)time)/CLOCKS_PER_SEC;
    printf ("\nCPU time for particle BFM: %f seconds.\n\n", elapsed_time);

    outfile << elapsed_time;

    // bf.create_interpolate_video(mu, nu, nt, plt);

    delete mu;
    delete nu;


    // TODO: interpolate video
} outfile << "\n";}

cout << "Herere\n";

outfile.close();

    

} // end of main