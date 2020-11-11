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
#include "Plotting.h"
#include "Points.h"
#include "BFM.h"

using namespace std;

const int DIM = 2;

int main(int argc, char** argv){
    if(argc!=7){
        cout<<"Do the following:"<<"\n";
        cout<< argv[0] << " [n1] [n2] [max_iteration] [tolerance] [smooth] [nt]"<<"\n";
        return 0;
    }

    int n1=stoi(argv[1]);
    int n2=stoi(argv[2]);
    int max_iteration=stoi(argv[3]);
    double tolerance=stod(argv[4]);
    double smooth=stod(argv[5]);
    int nt=stoi(argv[6]);

    

    string data_folder = "data";

    cout << "n1            : " << n1 <<"\n";
    cout << "n2            : " << n2 <<"\n";
    cout << "Max Iteration : " << max_iteration <<"\n";
    cout << "tolerance     : " << tolerance <<"\n";

    /* Initialize plotting tool */
    Plotting plt(n1,n2);

    create_csv_parameters(n1,n2);

    /* Initialize mu and nu */
    int num_points_mu = 4;
    int num_points_nu = 4;

    Points* mu = new Points(DIM, num_points_mu);
    Points* nu = new Points(DIM, num_points_nu);

    for(int p=0;p<mu->num_points_;++p){
        (*mu)(p,0) = 0.9*rand()/RAND_MAX + 0.05;
        (*mu)(p,1) = 0.9*rand()/RAND_MAX + 0.05;
    }

    for(int p=0;p<nu->num_points_;++p){
        (*nu)(p,0) = 0.9*rand()/RAND_MAX + 0.05;
        (*nu)(p,1) = 0.9*rand()/RAND_MAX + 0.05;
    }

    cout << "XXX Starting particle BFM XXX" << "\n";

    cout << "\n";
    // declaring argument of time() 
    time_t my_time = time(NULL);

    // ctime() used to give the present time 
    printf("%s", ctime(&my_time));

    /* Initialize BFM */
    BackAndForth bf(n1,n2,max_iteration,tolerance,smooth);

    string figurename = "barenblatt";

    plt.save_image_opencv(mu,nu,figurename,0);

    clock_t time;
    time=clock();

    bf.start_OT(mu, nu, plt);

    time=clock()-time;
    printf ("\nCPU time for particle BFM: %f seconds.\n\n",((float)time)/CLOCKS_PER_SEC);

    // TODO: interpolate video

    bf.create_interpolate_video(mu, nu, nt, plt);

}