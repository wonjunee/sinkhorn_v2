#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <iostream>

class vector2d{
public:
    double x;
    double y;
    vector2d():x(0),y(0){};
    vector2d(double x_, double y_):x(x_),y(y_){};
}; // vector 2d

#endif