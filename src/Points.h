#ifndef POINTS_H
#define POINTS_H

#include <iostream>

class PointND{
public:
    int DIM_;
    double* data_;

    PointND(){
        data_ = NULL;
    }

    PointND(int DIM){
        DIM_ = DIM;
        data_ = new double[DIM_];
        memset(data_, 0, DIM_*sizeof(double));
    }

    ~PointND(){
        delete[] data_;
    }

    double get(int n) const{
        assert(n >= 0 && n < DIM_);
        return data_[n];
    }

    double& operator()(int n)
    {
        assert(n >= 0 && n < DIM_);
        return data_[n];
    }
     
    double operator()(int n) const
    {
        assert(n >= 0 && n < DIM_);
        return data_[n];
    }
};

class Points{
public:
    int n_points_;
    int DIM_;
    // PointND** mu;
    double* mu;
    
    Points(){
        mu = NULL;
    }

    Points(int DIM, int n_points){

        n_points_ = n_points;
        DIM_ = DIM;

        // initialize points
        mu = new double[n_points_*DIM_]; // mu[point_idx * DIM_ + dimension_idx]
        // mu = new PointND*[n_points_];
        // for(int p=0;p<n_points_;++p){
        //     mu[p] = new PointND(DIM_);
        // }
    }

    ~Points(){
        delete[] mu;

        // for(int p=0;p<n_points_;++p){
        //     delete mu[p];
        // }
        // delete[] mu;
        mu = NULL;
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
     
        // return (*mu[p])(n);
        return mu[p*DIM_+n];
    }

    // PointND* get(int p) const{
    //     assert(p >= 0 && p < n_points_);
     
    //     return mu[p];
    // }
    
    double& operator()(int p, int n)
    {
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < n_points_);
     
        // return (*mu[p])(n);
        return mu[p*DIM_+n];
    }

    double* operator()(int p)
    {
        assert(p >= 0 && p < n_points_);
        return &mu[p*DIM_];
    }
     
    double operator()(int p, int n) const
    {
        assert(n >= 0 && n < DIM_);
        assert(p >= 0 && p < n_points_);
     
        // return (*mu[p])(n);
        return mu[p*DIM_+n];
    }

    // PointND* operator()(int p) const
    // {
    //     assert(p >= 0 && p < n_points_);
     
    //     return mu[p];
    // }

}; // Points

#endif