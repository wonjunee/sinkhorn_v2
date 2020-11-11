#ifndef POINTS_H
#define POINTS_H

#include <iostream>

class Points{
public:
    int num_points_;
    int DIM_;
    double* mu;
    
    Points(){
        mu = NULL;
    }

    Points(int DIM, int num_points){

        num_points_ = num_points;
        DIM_ = DIM;

        mu = new double[DIM_ * num_points_];
    }

    ~Points(){
        delete[] mu;
        mu = NULL;
    }

    int DIM() const{
        return num_points_;
    }

    int num_points() const{
        return num_points_;
    }

    double get(int p, int n) const{
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < num_points_);
     
        return mu[p*DIM_+n];
    }
    
    double& operator()(int p, int n)
    {
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < num_points_);
     
        return mu[p*DIM_+n];
    }
     
    double operator()(int p, int n) const
    {
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < num_points_);
     
        return mu[p*DIM_+n];
    }

}; // Points

#endif