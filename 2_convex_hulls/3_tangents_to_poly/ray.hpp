#ifndef RAY_H
#define RAY_H

#include "point.hpp"
#include "segment.hpp"

#include <sstream>
#include <vector>

namespace crs {

struct Ray {
    // a (horizontal, left) ray passing by point p
    Point p;
    
    bool intersects(Segment const& s) {

        double ymax = (s.b.y >= s.a.y) ? s.b.y : s.a.y;
        double ymin = (s.b.y < s.a.y) ? s.b.y : s.a.y;
        
        if (p.y >= ymin && p.y <= ymax) {

            // horizontal segment
            if (s.a.y == s.b.y) {
                double xmin = (s.a.x < s.b.x) ? s.a.x : s.b.x;
                if (p.x >= xmin)
                    return true;
                return false;
            }
            
            // vertical segment
            if (s.a.x == s.b.x) {
                if (p.x >= s.a.x)
                    return true;
                return false;
            }
            
            // diagonal segment. make sure it points up, decide whether
            // it point of the ray is on the right (or on the segment)
            Segment scmp;            
            if (s.a.y < s.b.y) {
                scmp.a = s.a;
                scmp.b = s.b;
            } else {
                scmp.a = s.b;
                scmp.b = s.a;
            }
            std::string location = scmp.point_location(p);
            
            if (location == "RIGHT" || location == "ON_SEGMENT" )
                return true;
            return false;      
        }

        return false;
    }
    
    std::vector<Point> intersections(Segment const& s) {
        // return 1 or 2 intersecting points of the segment with the Ray
        std::vector<Point> int_points;

        double ymax = (s.b.y >= s.a.y) ? s.b.y : s.a.y;
        double ymin = (s.b.y < s.a.y) ? s.b.y : s.a.y;
        
        if (p.y >= ymin && p.y <= ymax) {

            // horizontal segment
            if (s.a.y == s.b.y) {
                double xmin = (s.a.x < s.b.x) ? s.a.x : s.b.x;
                if (p.x >= s.a.x)
                    int_points.push_back(s.a);
                if (p.x >= s.b.x)
                    int_points.push_back(s.b);
            } else {
                // vertical segment
                if (s.a.x == s.b.x) {
                    if (p.x >= s.a.x)
                        int_points.push_back(Point(s.a.x, p.y));
                } else {            
                    // diagonal segment. make sure it points up, decide whether
                    // it point of the ray is on the right (or on the segment)
                    Segment scmp;            
                    if (s.a.y < s.b.y) {
                        scmp.a = s.a;
                        scmp.b = s.b;
                    } else {
                        scmp.a = s.b;
                        scmp.b = s.a;
                    }
                    std::string location = scmp.point_location(p);
                    
                    if (location == "RIGHT" || location == "ON_SEGMENT" ) {
                        double m = (scmp.b.y - scmp.a.y) / (scmp.b.x - scmp.a.x);
                        double x = scmp.a.x + (p.y - scmp.a.y) / m;
                        int_points.push_back(Point(x, p.y));
                    }
                }
            
            }
        }

        return std::move(int_points);
    }
};


std::ostream& operator<<(std::ostream& os, Ray r) {
    os << "ray (" << r.p << ")";
    return os;
}

    
}

#endif