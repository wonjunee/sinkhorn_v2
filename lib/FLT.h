#ifndef FLT_H
#define FLT_H

#include <iostream>
#include <vector>
#include "vector2d.h"

using namespace std;

class FLT{
public:
    int n1;

    double* points_x;
    double* points_y;
    double* convex_hull_vertices_x;
    double* convex_hull_vertices_y;

    vector2d* lower;

    int convex_hull_size;

    FLT(){
        points_x=NULL;
        points_y=NULL;
        convex_hull_vertices_x=NULL;
        convex_hull_vertices_y=NULL;

        lower=NULL;
    }

    FLT(int n1){
        this->n1=n1;

        points_x=new double[n1];
        points_y=new double[n1];
        convex_hull_vertices_x=new double[n1];
        convex_hull_vertices_y=new double[n1];
        lower = new vector2d[n1];
    }

    virtual ~FLT(){
        delete[] points_x;
        delete[] points_y;
        delete[] convex_hull_vertices_x;
        delete[] convex_hull_vertices_y;
        delete[] lower;
    }

    void setup_points(double* y_array){
        for(int i=0;i<n1;++i){
            points_x[i]=(i+0.5)/(1.0*n1);
            points_y[i]=y_array[i];
        }
    }

    double cross(vector2d& o, vector2d& a, vector2d& b){
        return 1.0*(a.x-o.x)*(b.y-o.y) - 1.0*(a.y-o.y)*(b.x-o.x);
    }

    void convex_hull(){

        if(n1<=1){
            return;
        }
        // build lower hull
        // vector<vector2d> lower;

        int counter=0;
        for(int i=0;i<n1;++i){
            vector2d p(points_x[i],points_y[i]);
            while(counter>=2){
                if(cross(lower[counter-2],lower[counter-1], p)>0){
                    break;
                }
                counter--;
            }
            lower[counter] = p;
            counter++;
        }
        convex_hull_size=0;
        for(int i=0;i<counter;++i){
            convex_hull_vertices_x[i]=lower[i].x;
            convex_hull_vertices_y[i]=lower[i].y;
            convex_hull_size++;
        }
    }   

    int find_index(double m, int start_ind){
        double first_slope=(1.0*convex_hull_vertices_y[1]-convex_hull_vertices_y[0])/(convex_hull_vertices_x[1]-convex_hull_vertices_x[0]);
        if(m<=first_slope){
            return 0;
        }
        for(int i=fmax(1,start_ind);i<convex_hull_size-1;++i){
            double slope_current=(1.0*convex_hull_vertices_y[i+1]-convex_hull_vertices_y[i])/(convex_hull_vertices_x[i+1]-convex_hull_vertices_x[i]);
            double slope_previous=(1.0*convex_hull_vertices_y[i]-convex_hull_vertices_y[i-1])/(convex_hull_vertices_x[i]-convex_hull_vertices_x[i-1]);
            if(slope_previous<m && m<=slope_current){
                return i;
            }
        }
        return convex_hull_size-1;
    }
    
    // Expect a size of n1
    void run_flt(double* y_star){
        convex_hull();
        int ind=1;
        for(int i=0;i<n1;++i){
            double m=(i+0.5)/(1.0*n1);
            ind=find_index(m,ind);
            double x=convex_hull_vertices_x[ind];
            double y=convex_hull_vertices_y[ind];
            y_star[i]=x*m-y;
        }
    }
}; // FLT

class FLT2D{
public:
    int n1;
    int n2;

    FLT* flt1;
    FLT* flt2;

    double* y_array;
    double* v_s;
    double* ustartmp;

    double* PHI;

    FLT2D(){
        flt1=NULL;
        flt2=NULL;
        y_array=NULL;
        v_s=NULL;
        ustartmp=NULL;
        PHI=NULL;
    }

    FLT2D(int n1, int n2){
        this->n1=n1;
        this->n2=n2;

        flt1 = new FLT(n1);
        flt2 = new FLT(n2);

        y_array=new double[n1];
        v_s=new double[n1*n2];

        ustartmp=new double[n1*n2];

        PHI=new double[n1*n2];
    }

    virtual ~FLT2D(){
        delete[] y_array;
        delete[] v_s;
        delete[] ustartmp;
        delete[] PHI;

        delete flt1;
        delete flt2;
    }

    /** 
        Find ustar from u 2d matrix
    */
    void run_flt2d(const double* u, double* ustar){

        for(int x2_index=0;x2_index<n2;++x2_index){
            // create y array for x1
            for(int j=0;j<n1;++j){
                y_array[j]=u[x2_index*n1+j];
            }

            // Step 3: setup points
            flt1->setup_points(y_array);

            // Step 4: Perform Legendre transform i.e. y_array*
            flt1->run_flt(&v_s[x2_index*n1]);    
        }
        

        for(int x1_index=0;x1_index<n1;++x1_index){
            // create y array for x1
            for(int i=0;i<n2;++i){
                y_array[i]=-v_s[i*n1+x1_index];
            }

            // Step 3: setup points
            flt2->setup_points(y_array);

            // Step 4: Perform Legendre transform i.e. y_array*
            flt2->run_flt(&ustartmp[x1_index*n2]);   
        }
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                ustar[i*n1+j]=ustartmp[j*n2+i];
            }
        }   
    }


    void find_c_concave(const double* phi, double* phi_c, const double tau){
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double x=(j+0.5)/(1.0*n1);
                double y=(i+0.5)/(1.0*n2);
                PHI[i*n1+j]=0.5*(x*x+y*y)-tau*phi[i*n1+j];
            }
        }

        run_flt2d(PHI,phi_c);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double x=(j+0.5)/(1.0*n1);
                double y=(i+0.5)/(1.0*n2);
                phi_c[i*n1+j]=1.0/tau*(0.5*(x*x+y*y)-phi_c[i*n1+j]);
            }
        }
    }
}; // FLT2d


#endif