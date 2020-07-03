#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <iostream>
#include <fstream>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace std;
using namespace cv;

class Initializer{
public:
    int n1;
    int n2;

    Mat img_in;
    Mat img_color;
    unsigned char* ar;

    Initializer(int n1, int n2){
        this->n1 = n1;
        this->n2 = n2;

        img_in = Mat(n2,n1,CV_8U);
        ar = new unsigned char[n1*n2];
    }

    ~Initializer(){
        delete[] ar;
    }

    void save_image_opencv(const double* A, const string& filename, const int iter){
        string filename_iter = "./figures/" + filename;
        filename_iter = filename_iter + "-" + to_string(iter);
        filename_iter = filename_iter + ".png";

        double max_A = 0;
        for(int i=0;i<n1*n2;++i) max_A = fmax(max_A, A[i]);
        for(int i=0;i<n1*n2;++i) ar[i] = (unsigned char) (A[i]/max_A * 255);

        memcpy(img_in.data, ar, n1*n2*sizeof(unsigned char));

        // Apply the colormap:
        applyColorMap(img_in, img_color, COLORMAP_INFERNO);

        imwrite(filename_iter, img_color);  
    }

    void save_image_opencv(const double* A,const unsigned char* obstacle,const string& filename, const int iter){
        string filename_iter = "./figures/" + filename;
        filename_iter = filename_iter + "-" + to_string(iter);
        filename_iter = filename_iter + ".png";

        double max_A = 0;
        for(int i=0;i<n1*n2;++i) max_A = fmax(max_A, A[i]);
        for(int i=0;i<n1*n2;++i) ar[i] = (unsigned char) (A[i]/max_A * 255);

        memcpy(img_in.data, ar, n1*n2*sizeof(unsigned char));

        // Apply the colormap:
        applyColorMap(img_in, img_color, COLORMAP_INFERNO);

        for(int i=0;i<n1;++i){
            for(int j=0;j<n2;++j){
                if(obstacle[i*n1+j] != 0){
                    Vec3b & color = img_color.at<Vec3b>(i,j);
                    color[0] = 255;
                    color[1] = 255;
                    color[2] = 255;    
                }
            }
        }

        imwrite(filename_iter, img_color);  
    }


    void init_entropy_image_obstacle_opencv(double* nu, const unsigned char* obstacle, const string& data_folder){

        double xc=0.8;
        double yc=0.8;

        double a = 0;

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