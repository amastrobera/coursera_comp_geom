#ifndef POINT_H
#define POINT_H

#include <sstream>

namespace crs {

struct Point {
    double x,y;
    
    Point() {} 
    
    Point(double c1, double c2) : x(c1), y(c2) { }
    
    Point& operator=(Point const& cp) {
        if (this != &cp) {
            x = cp.x;
            y = cp.y;
        }
        return *this;
    }
    
    Point(Point const& cp) {
        *this = cp;
    }
};


bool operator==(Point const& p1, Point const& p2) {
    return p1.x == p2.x && p1.y == p2.y;
}

bool operator!=(Point const& p1, Point const& p2) {
    return !(p1 == p2);
}

std::ostream& operator<<(std::ostream& os, Point p) {
    os << p.x << " " << p.y;
    return os;
}


} 


#endif