#ifndef POINT_H
#define POINT_H

#include <sstream>

namespace crs {
struct Point {
    double x,y;
    
    Point() : x(0), y(0) {} 
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
    
    struct hash{
        size_t operator()(const Point &p) const {
            size_t h1 = std::hash<double>()(p.x);
            size_t h2 = std::hash<double>()(p.y);
            return h1 ^ (h2 << 1);
        }
    };
};


bool operator==(Point const& p1, Point const& p2) {
    return p1.x == p2.x && p1.y == p2.y;
}

bool operator!=(Point const& p1, Point const& p2) {
    return !(p1 == p2);
}

bool operator<(Point const& p1, Point const& p2) {
    return p1.y < p2.y || (p1.y == p2.y && p1.x < p2.x);
}


std::ostream& operator<<(std::ostream& os, Point p) {
    os << p.x << " " << p.y;
    return os;
}


} 


#endif