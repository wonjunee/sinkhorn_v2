#ifndef PLOTTING_H
#define PLOTTING_H

#include <iostream>
#include <fstream>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "Points.h"

using namespace std;
using namespace cv;

class Plotting{
public:
    int n1;
    int n2;
    int DIM_;

    Mat img_in;
    Mat img_color;
    unsigned char* ar;

    Plotting(int n1, int n2){
        this->n1 = n1;
        this->n2 = n2;

        DIM_ = 2;

        img_in = Mat(n2,n1,CV_8U);
        memset(img_in.data, 0, n1*n2*sizeof(unsigned char));
        img_color = Mat(n2,n1,CV_8UC3, Scalar(10, 100, 150));
        ar = new unsigned char[n1*n2];
    }

    ~Plotting(){
        printf("deconst plotting\n");
        delete[] ar;
        printf("deconst plotting done\n");
    }

    void save_image_opencv(const double* A, const string& filename, const int iter){
        string filename_iter = "./figures/" + filename;
        filename_iter = filename_iter + "-" + to_string(iter) + ".png";

        double min_A = A[0];
        double max_A = A[0];
        max_A = A[0];
        for(int i=1;i<n1*n2;++i){
            max_A = fmax(max_A, A[i]);
            min_A = fmin(min_A, A[i]);
        }
        
        for(int i=0;i<n1*n2;++i) ar[i] = (unsigned char) ((A[i]-min_A)/(max_A-min_A) * 255);

        memcpy(img_in.data, ar, n1*n2*sizeof(unsigned char));

        // Apply the colormap:
        applyColorMap(img_in, img_color, COLORMAP_INFERNO);

        imwrite(filename_iter, img_color);  
    }

    void save_image_opencv(Points* push_mu_, Points* mu, Points* nu, const string& filename, const int iter, const int max_iter){
        string filename_iter = "./figures/" + filename;
        filename_iter = filename_iter + "-" + to_string(iter) + ".png";

        // memcpy(img_in.data, ar, n1*n2*sizeof(unsigned char));
        // applyColorMap(img_in, img_color, COLORMAP_INFERNO);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double x = (j+.5)/n1;
                double y = (i+.5)/n2;

                Vec3b & color = img_color.at<Vec3b>(i,j);

                color[0] = 0;
                color[1] = 0;
                color[2] = 0;  

                for(int p=0;p<mu->num_points();++p){
                    double mux = mu->get(p,0);
                    double muy = mu->get(p,1);

                    double pushx = push_mu_->get(p,0);
                    double pushy = push_mu_->get(p,1);

                    double px = 1.0*iter/max_iter*pushx + (1-1.0*iter/max_iter)*mux;
                    double py = 1.0*iter/max_iter*pushy + (1-1.0*iter/max_iter)*muy;

                    if(pow(x-px,2) + pow(y-py,2) < pow(0.01,2)){
                        color[0] = 255;
                        color[1] = 0;
                        color[2] = 0;  
                    }
                }

                for(int p=0;p<nu->num_points();++p){
                    double px = nu->get(p,0);
                    double py = nu->get(p,1);

                    if(pow(x-px,2) + pow(y-py,2) < pow(0.005,2)){
                        color[0] = 0;
                        color[1] = 0;
                        color[2] = 255;  
                    }
                }
            }
        }

        imwrite(filename_iter, img_color);
    }

    void save_image_opencv(Points* A, Points* B, const string& filename, const int iter){
        string filename_iter = "./figures/" + filename;
        filename_iter = filename_iter + "-" + to_string(iter) + ".png";

        // memcpy(img_in.data, ar, n1*n2*sizeof(unsigned char));
        // applyColorMap(img_in, img_color, COLORMAP_INFERNO);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){

                double x = (j+.5)/n1;
                double y = (i+.5)/n2;

                Vec3b & color = img_color.at<Vec3b>(i,j);

                color[0] = 0;
                color[1] = 0;
                color[2] = 0;  

                for(int p=0;p<A->num_points();++p){
                    double px = A->get(p,0);
                    double py = A->get(p,1);

                    if(pow(x-px,2) + pow(y-py,2) < pow(0.01,2)){
                        color[0] = 255;
                        color[1] = 0;
                        color[2] = 0;  
                    }
                }

                for(int p=0;p<B->num_points();++p){
                    double px = B->get(p,0);
                    double py = B->get(p,1);

                    if(pow(x-px,2) + pow(y-py,2) < pow(0.005,2)){
                        color[0] = 0;
                        color[1] = 0;
                        color[2] = 255;  
                    }
                }
            }
        }

        imwrite(filename_iter, img_color);
    }


    void init_entropy_image_obstacle_opencv(double* nu, const unsigned char* obstacle, const string& data_folder){

        double xc=0.9;
        double yc=0.9;

        double a = 10;

        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
                
                double x=(j+.5)/(n1*1.0);
                double y=(i+.5)/(n2*1.0);

                int cc = (int) obstacle[i*n1+j];

                if(cc>10){
                    nu[i*n1+j]=-1;
                }else{
                    nu[i*n1+j]=a*((x-xc)*(x-xc)+(y-yc)*(y-yc));
                }                  
            }
        }
    }

};

#endif