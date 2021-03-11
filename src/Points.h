#ifndef POINTS_H
#define POINTS_H

#include <iostream>

class Points{
public:
    int n_points_;
    int DIM_;
    double* data;
    
    Points(){
        data = NULL;
    }

    Points(int DIM, int n_points){

        n_points_ = n_points;
        DIM_ = DIM;

        // initialize points
        data = new double[n_points_*DIM_]; // data[point_idx * DIM_ + dimension_idx]
    }

    ~Points(){
        delete[] data;
        data = NULL;
    }

    int DIM() const{
        return DIM_;
    }

    int num_points() const{
        return n_points_;
    }

    double get(int p, int n) const{
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < n_points_);
     
        // return (*data[p])(n);
        return data[p*DIM_+n];
    }
    
    double& operator()(int p, int n)
    {
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < n_points_);
     
        // return (*data[p])(n);
        return data[p*DIM_+n];
    }

    double* operator()(int p)
    {
        assert(p >= 0 && p < n_points_);
        return &data[p*DIM_];
    }
     
    double operator()(int p, int n) const
    {
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < n_points_);
        return data[p*DIM_+n];
    }

}; // Points

#endif