#ifndef ACCESSORY_H
#define ACCESSORY_H

#include <iostream>
#include <fstream>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace std;
using namespace cv;

void create_csv_parameters(int n1,int n2){
    ofstream outfile;
    string filename = "./data//parameters.csv";
    outfile.open(filename);
    outfile<<n1<<"\n"<<n2;
    outfile.close();
}

void create_bin_file(const double* A, int size, string filename){
    ofstream out(filename, ios::out | ios::binary);
    if(!out) {
        cout << "Cannot open file.";
        return;
    }

    out.write((char *) A, size*sizeof(double));
    out.close();
}


void prepare_pixels(double* rho, unsigned char* pixels, int pcount, double max){

    max=0;
    for(int i=0;i<pcount;++i){
        max=fmax(max,rho[i]);
    }
    for(int i=0;i<pcount;i++){
        double val=255*rho[i]/max;
        pixels[i]=fmin(255,val);
    }
}



void create_mu_gaussian(double* mu,double px,double py,double var,int n1,int n2){
    double sum=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            mu[i*n1+j]=exp(-(pow(x-px,2)+pow(y-py,2))/pow(var,2));    
            sum+=mu[i*n1+j];
        }
    }
    for(int i=0;i<n1*n2;++i){
        mu[i]/=sum/(1.0*n1*n2);
    }
}

void create_mu_square(double* mu,double px,double py,double width,int n1,int n2){
    double sum=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(fabs(x-px)<width && fabs(y-py)<width){
                mu[i*n1+j] = 1;    
            }
            sum+=mu[i*n1+j];
        }
    }
    for(int i=0;i<n1*n2;++i){
        mu[i]/=sum/(1.0*n1*n2);
    }
}


void create_mu_PME(double* mu,const double kappa, const int n1, const int n2, const double xi0){
    double sum=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            double xx = x-0.5;
            double yy = y-0.5;

            mu[i*n1+j] = pow(fmax(0, xi0*xi0 - xx*xx - yy*yy),1.0/(kappa-1));
            sum+=mu[i*n1+j];
        }
    }

    sum = 1.0*n1*n2/sum;
    for(int i=0;i<n1*n2;++i){
        mu[i]*=sum;
    }

    std::cout  << "C = " <<setprecision(100)<< sum << std::endl;
}

void create_mu_one_circle(double* mu,int n1,int n2){

    double r = 0.15;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(pow(x-0.2,2)+pow(y-0.2,2)<pow(r,2)){
                mu[i*n1+j] = 1;    
            }
        }
    }
}

void create_mu_flavien_test(double* mu,int n1,int n2){

    double r = 0.2;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            // mu[i*n1+j]  = 1.6*exp( - (pow(x-0.5,2)+pow(y-0.5,2))/(2*0.15*0.15));
            mu[i*n1+j]  = 1.5*exp( - (pow(x-0.25,2)+pow(y-0.25,2))/(2*0.1*0.1));
            mu[i*n1+j]  += 1.5*exp( - (pow(x-0.75,2)+pow(y-0.75,2))/(2*0.1*0.1));
            mu[i*n1+j]  += 1.5*exp( - (pow(x-0.80,2)+pow(y-0.38,2))/(2*0.06*0.06));
            mu[i*n1+j]  += 1.5*exp( - (pow(x-0.66,2)+pow(y-0.18,2))/(2*0.06*0.06));
            mu[i*n1+j]  += 1.5*exp( - (pow(x-0.24,2)+pow(y-0.69,2))/(2*0.082*0.082));
            // mu[i*n1+j]  = 1.6*exp( - (pow(x-0.2,2)+pow(y-0.2,2))/(2*0.15*0.15));
            // mu[i*n1+j] += 1.6*exp( - (pow(x-0.8,2)+pow(y-0.8,2))/(2*0.15*0.15));
            // mu[i*n1+j] += 1.6*exp( - (pow(x-0.2,2)+pow(y-0.8,2))/(2*0.1*0.1));
            // mu[i*n1+j] += 1.6*exp( - (pow(x-0.8,2)+pow(y-0.2,2))/(2*0.1*0.1));
            // mu[i*n1+j] += 1.6*exp( - (pow(x-0.8,2)+pow(y-0.8,2))/(2*0.1*0.1));
        }
    }
}

void create_mu_two_circles(double* mu,int n1,int n2){

    // double r = 0.15;
    double r = 0.125;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(pow(x-0.175,2)+pow(y-0.175,2)<pow(r,2)){
                mu[i*n1+j] = 1;
            }
            if(pow(x-0.825,2)+pow(y-0.825,2)<pow(r,2)){
                mu[i*n1+j] = 1;
            }

            // if(pow(x-0.9,2)+pow(y-0.1,2)<pow(r,2)){
            //     mu[i*n1+j] = 1;
            // }
        }
    }
}

void create_mu_four_circles(double* mu,int n1,int n2){

    // double r = 0.15;
    double r = 0.1;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(pow(x-0.1,2)+pow(y-0.1,2)<pow(r,2)){
                mu[i*n1+j] = 1;
            }
            if(pow(x-0.9,2)+pow(y-0.9,2)<pow(r,2)){
                mu[i*n1+j] = 1;
            }
            if(pow(x-0.1,2)+pow(y-0.9,2)<pow(r,2)){
                mu[i*n1+j] = 1;
            }
            if(pow(x-0.9,2)+pow(y-0.1,2)<pow(r,2)){
                mu[i*n1+j] = 1;
            }

            // if(pow(x-0.9,2)+pow(y-0.1,2)<pow(r,2)){
            //     mu[i*n1+j] = 1;
            // }
        }
    }
}
void create_mu_one_square(double* mu,int n1,int n2){
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(fabs(x-0.2)<0.1 && fabs(y-0.2)<0.1){
                mu[i*n1+j] = 1;    
            }
        }
    }
}

void create_mu_two_square(double* mu,int n1,int n2){
    double sum=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(fabs(x-0.4)<0.1 && fabs(y-0.4)<0.1){
                mu[i*n1+j] = 1;    
            }
            if(fabs(x-0.6)<0.1 && fabs(y-0.6)<0.1){
                mu[i*n1+j] = 1;    
            }
            sum+=mu[i*n1+j];
        }
    }
    for(int i=0;i<n1*n2;++i){
        mu[i]/=sum/(1.0*n1*n2);
    }
}

void create_mu_four_square(double* mu,int n1,int n2){
    double sum=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;
            if(fabs(x-0.45)<0.1 && fabs(y-0.45)<0.1){
                mu[i*n1+j] += 1;    
            }
            if(fabs(x-0.55)<0.1 && fabs(y-0.55)<0.1){
                mu[i*n1+j] += 1;    
            }

            if(fabs(x-0.3)<0.1 && fabs(y-0.7)<0.1){
                mu[i*n1+j] += 1;    
            }

            if(fabs(x-0.7)<0.1 && fabs(y-0.3)<0.1){
                mu[i*n1+j] += 1;    
            }
            sum+=mu[i*n1+j];
        }
    }
    for(int i=0;i<n1*n2;++i){
        mu[i]/=sum/(1.0*n1*n2);
    }
}

void create_mu_from_image(double* mu,int n1,int n2){

    double nuSum=0;
    int pcount=n1*n2;

    ifstream infile;
    infile.open("../external-images/starobstacle_gray_512.dat");

    double x;
    int ind=0;
    while(infile>>x){
        mu[ind]=x;
        nuSum+=x;
        ind++;
    }
    infile.close();

    nuSum/=pcount;    
    
    for(int k=0;k<pcount;k++){
        mu[k]/=nuSum;
    } 
}

void create_mu_from_image2(double* mu,int n1,int n2){

    double nuSum=0;
    int pcount=n1*n2;

    ifstream infile;
    infile.open("../external-images/starobstacle_gray_512.dat");

    double x;
    int ind=0;
    while(infile>>x){
        mu[ind]=x;
        nuSum+=x;
        ind++;
    }
    infile.close();

    double val = mu[n2/2*n1 + n1/2];

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x = (j+0.5)/n1;
            double y = (i+0.5)/n2;

            if(fabs(x-0.2)<0.1 && fabs(y-0.2)<0.1){
                mu[i*n1+j] = val;
            }

            if(pow(x-0.2,2) + pow(y-0.8,2) < pow(0.1,2)){
                mu[i*n1+j] = val;   
            }

            if(pow(x-0.8,2) + pow(y-0.2,2) < pow(0.1,2)){
                mu[i*n1+j] = val;   
            }

            nuSum+=mu[i*n1+j];
        }
    }

    nuSum/=pcount;
    
    for(int k=0;k<pcount;k++){
        mu[k]/=nuSum;
    } 
}

// Quadratic potential V(x,y) = 1/2 a ( (x-xc)^2 + (y-yc)^2 )
void init_entropy_quadratic(double* nu,double xc, double yc, double a, int n1, int n2){    
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);

            nu[i*n1+j]=a*((x-xc)*(x-xc)+(y-yc)*(y-yc));
        }
    }
}


void init_entropy_sine(double* nu,double s1, double s2, double a, int n1, int n2){        
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);

            nu[i*n1+j]=a*(-sin(s1*M_PI*x)*sin(s2*M_PI*y)+1);
        }
    }
}

// obstacle (disk/square) with radius/side r centered at the middle
void init_entropy_one_circles_obstacle(double* nu, int n1, int n2){

    // Quadratic potential V = 1/2 a ( (x-xc)^2 + (y-yc)^2 )
    double xc=0.8;
    double yc=0.8;

    double r = 0.2;
    double a=  0.5;

    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);

            // bool a1=(x-0.3)*(x-0.3)+(y-0.7)*(y-0.7)<r*r;
            // bool a2=(x-0.7)*(x-0.7)+(y-0.3)*(y-0.3)<r*r;

            bool a1=(x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<r*r;

            if(!a1){
                nu[i*n1+j]=a*((x-xc)*(x-xc)+(y-yc)*(y-yc));
            }else{
                nu[i*n1+j]=-1;
            }                                  
        }
    }
}

// obstacle (disk/square) with radius/side r centered at the middle
void init_entropy_two_circles_obstacle(double* nu, double r, double a, int n1, int n2){

    // Quadratic potential V = 1/2 a ( (x-xc)^2 + (y-yc)^2 )
    double xc=0.8;
    double yc=0.8;

    r = 0.25;

    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);

            bool a1=(x-0.3)*(x-0.3)+(y-0.7)*(y-0.7)<r*r;
            bool a2=(x-0.7)*(x-0.7)+(y-0.3)*(y-0.3)<r*r;

            // bool a1=(x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<r*r;
            // bool a2=false;

            if(!a1 && !a2){
                nu[i*n1+j]=a*((x-xc)*(x-xc)+(y-yc)*(y-yc));
            }else{
                nu[i*n1+j] = -1;
            }
        }
    }
}

// obstacle (disk/square) with radius/side r centered at the middle
void init_entropy_moonlight_obstacle(double* nu, int n1, int n2){

    // Quadratic potential V = 1/2 a ( (x-xc)^2 + (y-yc)^2 )
    double xc=0.8;
    double yc=0.8;

    double r1 = 0.2;
    double r2 = 0.15;

    double a = 0.5;

    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);

            bool a1=(x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<r1*r1;
            bool a2= x < 0.5;
            bool a3= y < 0.5;

            if(!(a1 && !(a2 && a3))) {
                nu[i*n1+j]=a*((x-xc)*(x-xc)+(y-yc)*(y-yc));
            }else{
                nu[i*n1+j] = -1;
            }
        }
    }
}

void init_entropy_image_obstacle(double* nu, double gamma, double r, double a, int n1, int n2){

    // Quadratic potential V = 1/2 a ( (x-xc)^2 + (y-yc)^2 )

    double* image=new double[n1*n2];

    double xc=0.9;
    double yc=0.9;

    double nuSum=0;
    int pcount=n1*n2;

    ifstream infile;
    infile.open("maze4_gray.dat");

    double x;
    int ind=0;
    while(infile>>x){
        int i=ind/n1;
        int j=ind%n1;
        image[i*n1+j]=x;
        ind++;
    }
    infile.close();

    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);

            if(image[i*n1+j]>0){
                nu[i*n1+j]=exp(-a*((x-xc)*(x-xc)+(y-yc)*(y-yc))/(2*gamma));
                nuSum+=nu[i*n1+j];
            }                                        
        }
    }
    nuSum/=pcount;    
    
    for(int k=0;k<pcount;k++){
        nu[k]/=nuSum;
    } 

    delete[] image;
}


/* create a circle obstacle */
void init_obstacle_circle(unsigned char* obstacle, double px, double py, double r, int n1, int n2){

    // double px = 0.5;
    // double py = 0.5;
    // double r  = 0.2;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x = (j+.5)/n1;
            double y = (i+.5)/n2;

            if((x-px)*(x-px) + (y-py)*(y-py) < r*r) obstacle[i*n1+j] = 100;
            else obstacle[i*n1+j] = 0;
        }
    }
}


/* create a circle obstacle */
void init_obstacle_two_circle(unsigned char* obstacle, double px, double py, double r, double px2, double py2, double r2, int n1, int n2){
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x = (j+.5)/n1;
            double y = (i+.5)/n2;

            if((x-px)*(x-px) + (y-py)*(y-py) < r*r) obstacle[i*n1+j] = 100;
            else if((x-px2)*(x-px2) + (y-py2)*(y-py2) < r2*r2) obstacle[i*n1+j] = 100;
            else obstacle[i*n1+j] = 0;
        }
    }
}


/* create a circle obstacle */
void init_obstacle_four_circle(unsigned char* obstacle, double px, double py, double r, double px2, double py2, double r2, double px3, double py3, double r3, double px4, double py4, double r4, int n1, int n2){
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x = (j+.5)/n1;
            double y = (i+.5)/n2;

            if((x-px)*(x-px) + (y-py)*(y-py) < r*r) obstacle[i*n1+j] = 100;
            else if((x-px2)*(x-px2) + (y-py2)*(y-py2) < r2*r2) obstacle[i*n1+j] = 100;
            else if((x-px3)*(x-px3) + (y-py3)*(y-py3) < r3*r3) obstacle[i*n1+j] = 100;
            else if((x-px4)*(x-px4) + (y-py4)*(y-py4) < r4*r4) obstacle[i*n1+j] = 100;
            else obstacle[i*n1+j] = 0;
        }
    }
}

/* create a circle obstacle */
void init_obstacle_two_moons(unsigned char* obstacle, int n1, int n2){

    double px = 0.5;
    double py = 0.5;
    double r  = 0.2;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x = (j+.5)/n1;
            double y = (i+.5)/n2;

            bool a1 = pow(x-0.35,2)+pow(y-0.35,2) < pow(0.1,2);
            bool a2 = pow(x-0.30,2)+pow(y-0.30,2) < pow(0.08,2);
            bool b1 = pow(x-0.65,2)+pow(y-0.65,2) < pow(0.1,2);
            bool b2 = pow(x-0.70,2)+pow(y-0.70,2) < pow(0.08,2);
            // bool c1 = pow(x-0.75,2)+pow(y-0.25,2) < pow(0.1,2);
            // bool c2 = pow(x-0.80,2)+pow(y-0.20,2) < pow(0.08,2);
            // bool d1 = pow(x-0.25,2)+pow(y-0.75,2) < pow(0.1,2);
            // bool d2 = pow(x-0.20,2)+pow(y-0.80,2) < pow(0.08,2);

            if( (a1 && !a2) || (b1 && !b2) ) obstacle[i*n1+j] = (unsigned char) 200;
            // if( a1 ) obstacle[i*n1+j] = (unsigned char) 100;
            else obstacle[i*n1+j] = (unsigned char)  0;
        }
    }
}

/* create a pac man obstacle */
void init_obstacle_pac_man(unsigned char* obstacle, int n1, int n2){

    double px = 0.5;
    double py = 0.5;
    double r  = 0.2;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x = (j+.5)/n1;
            double y = (i+.5)/n2;

            bool a1 = pow(x-0.5,2)+pow(y-0.5,2) < pow(0.2,2);
            // bool c1 = pow(x-0.75,2)+pow(y-0.25,2) < pow(0.1,2);
            // bool c2 = pow(x-0.80,2)+pow(y-0.20,2) < pow(0.08,2);
            // bool d1 = pow(x-0.25,2)+pow(y-0.75,2) < pow(0.1,2);
            // bool d2 = pow(x-0.20,2)+pow(y-0.80,2) < pow(0.08,2);

            if( a1 && !(x<0.5 && y<0.5) ) obstacle[i*n1+j] = (unsigned char) 200;
            // if( a1 ) obstacle[i*n1+j] = (unsigned char) 100;
            else obstacle[i*n1+j] = (unsigned char)  0;
        }
    }
}


/* create an obstacle from an external image */
void init_entropy_from_image(double* V, int n1, int n2){

    // Quadratic potential V = 1/2 a ( (x-xc)^2 + (y-yc)^2 )

    // Mat src = imread("../external-images/paw.png");
    // Mat src = imread("../external-images/two-moons.jpeg");
    Mat src = imread("../external-images/fruits.png");

    Mat src_gray;

    resize(src,src_gray,Size(n2,n1));

    cvtColor(src_gray, src_gray, COLOR_BGR2GRAY );

    src_gray.convertTo(src_gray, CV_8U);

    Mat dst;

    int DELAY_BLUR = 100;
    int MAX_KERNEL_LENGTH = 5;

    for ( int i = 1; i < MAX_KERNEL_LENGTH; i = i + 2 )
    {
        GaussianBlur( src_gray, dst, Size( i, i ), 0, 0 );
    }

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            int a = (int) dst.ptr<unsigned char>(i)[j];

            V[i*n1+j] = a/255.0;
        }
    }

    if(true){  // show image
        const char* source_window = "Potential";
        namedWindow( source_window );
        imshow( source_window, dst );

        waitKey(0);
    }
}



/* create an obstacle from an external image */
void init_obstacle_from_image(unsigned char* obstacle, int n1, int n2){

    // Quadratic potential V = 1/2 a ( (x-xc)^2 + (y-yc)^2 )

    // Mat src = imread("../external-images/paw.png");
    // Mat src = imread("../external-images/two-moons.jpeg");
    Mat src = imread("../external-images/starobstacle.png");

    Mat src_gray;

    resize(src,src_gray,Size(n2,n1));

    cvtColor(src_gray, src_gray, COLOR_BGR2GRAY );

    src_gray.convertTo(src_gray, CV_8U);

    /// Erode

    int erosion_size = 1;

    int erosion_type = MORPH_RECT;

    //![kernel]
    Mat element = getStructuringElement( erosion_type,
                         Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                         Point( erosion_size, erosion_size ) );
    //![kernel]

    Mat src_eroded;
    /// Apply the erosion operation
    erode(src_gray, src_eroded, element);

    

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            int a = (int) src_eroded.ptr<unsigned char>(i)[j];

            if(a > 100){
                obstacle[i*n1+j] = 255;
            }else{
                obstacle[i*n1+j] = 0;
            }
        }
    }

    if(false){  // show image
        const char* source_window = "Potential";
        namedWindow( source_window );
        imshow( source_window, src_eroded );

        waitKey(0);
    }
}


#endif