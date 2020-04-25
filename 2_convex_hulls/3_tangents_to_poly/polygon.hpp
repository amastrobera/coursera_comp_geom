#ifndef POLYGON_H
#define POLYGON_H

#include "point.hpp"
#include "segment.hpp"
#include "ray.hpp"

 
#include <iostream> // error msg
#include <sstream> // string operators
#include <vector> // main container 
#include <unordered_set> // non cvx function
#include <utility> 

namespace crs {

struct Polygon {

    std::vector<Point> vertices;

    std::string type() {
        if (is_convex())
            return "CONVEX";
        return "NOT_CONVEX";
    }

    bool is_convex() {
        // checks that the polygon is convex by using right/left turn definition
        // assumes vertices are in counter clockwise order
        bool is_ccw = true;
        size_t n = vertices.size();
        for (size_t i1 = 0; i1 < n; ++i1) {
            size_t i2 = (i1 + 1) % n; // circular list! 
            size_t i3 = (i1 + 2) % n;
            std::string loc = Segment(vertices[i1], vertices[i2]
                                    ).point_location(vertices[i3]);
            if (loc == "RIGHT" || loc == "ON_SEGMENT") {
                is_ccw = false;
                break;
            }
        }
        return is_ccw;
    }

    bool is_valid() {
        // checks that the polygon has at least 3 edges
        //      and those are in counter-clockwise
        //      and no 3 consecutive vertices are collinear
        //  note: I am not checking if the polygon is 
        //        not self touching / self intersecting
        //        it seems hard and we can look at it later.  
        if (vertices.size() < 3)
            return false;   
        double sum_ccw = 0;
        // ccw is achieved if sum((x2-x1)(y2+y1)) < 0 (if > 0 cw)
        // https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        bool is_not_collinear = true;
        size_t n = vertices.size();
        for (size_t i1 = 0; i1 < n; ++i1) {
            size_t i2 = (i1 + 1) % n; // circular list! 
            size_t i3 = (i1 + 2) % n;
            if (Segment(vertices[i1], vertices[i2]).point_location(vertices[i3]) 
                == "ON_LINE") {
                is_not_collinear = false;
                break;
            }
            sum_ccw += (vertices[i2].x - vertices[i1].x)*(vertices[i2].y + vertices[i1].y);
        }
        if (sum_ccw >= 0) {
            std::cerr << "error: no counter clockwise vertices" << std::endl;
            return false;
        }
        if (!is_not_collinear) {
            std::cerr << "error: collinear points" << std::endl;
            return false;
        }
        return true;
    }


    std::string point_location(Point const& p) {
        if (is_convex())
            return point_location_cvx(p); // nicer looking, used only for convex poly
        return point_location_not_cvx(p); // used in any case
    }
    
    std::pair<Point, Point> tangent_points_to_point(Point const& p) {
        // pick up two tangents in linear time
        //   tangets are points t1, t2 such that 
        //      t1->p followed by p->t2 for two new edges of the
        //      new convex_hull that contains p
        std::pair<Point, Point> tangents;
        size_t n = vertices.size();
        
        if (n) {
            if (n > 1) {
                if (n > 2) {
                    Point first, second;
                    // first point
                    std::vector<Point>::const_iterator it1 = vertices.begin();
                    for (; it1 != vertices.end(); ++it1) {
                        auto it2 = (it1 != vertices.end()-1) ? it1+1 : vertices.begin(); // circular list
                        std::string loc = Segment(*it1, p).point_location(*it2);
                        if ( loc != "RIGHT") {
                            first = *it1;
                            break;
                        }
                    }
                    // second point 
                    if (it1 != vertices.end()) {
                        std::vector<Point>::const_iterator it2 = (it1 != vertices.end()-1) ? it1+1 : vertices.begin(); 
                        for (; it2 != it1+n; ++it2) {
                            auto it3 = (it2 != vertices.end()-1) ? it2+1 : vertices.begin(); // circular list
                            std::string loc = Segment(p, *it2).point_location(*it3);
                            if ( loc != "RIGHT") {
                                second = *it2;
                                break;
                            }
                        }
                        if (it2 != it1+n) {
                            tangents.first =  first;
                            tangents.second = second;
                        } else 
                            std::cerr << "second tanget point not found. skipping calc" << std::endl;
                    } else 
                        std::cerr << "first tanget point not found. skipping calc" << std::endl;

                } else 
                    tangents = {vertices[0], vertices[1]};
            } else 
                tangents = {vertices[0], vertices[0]};
        }
        
        return std::move(tangents);
    }

private:

    std::string point_location_cvx(Point const& p) {
        
        bool is_inside = true;
        bool is_border = false;
        size_t n = vertices.size();
        for (size_t i1 = 0; i1 < n; ++i1) {
            size_t i2 = (i1 + 1) % n; // circular list! 
            std::string loc = Segment(vertices[i1], vertices[i2]).point_location(p);
            if (loc == "LEFT")
                is_inside &= true;
            else if (loc == "ON_SEGMENT") {
                is_border = true;
                break;
            } else {
                is_inside = false;
                break;
            }
        }
        if (is_border)
            return "BORDER";
        if (is_inside)
            return "INSIDE";
        return "OUTSIDE";
    }
    
    std::string point_location_not_cvx(Point const& p) {
    
        // logic: 
        //      1. check if point is a a border --> BORDER
        //      2. trace a ray and count the intersections
        //          2.1 UNIQUE intersections (handles the case of ray intersecting 
        //              a vertex (2 edges)
        //          2.2 if num(intersections) is odd --> INSIDE
        //              also handles the case of ray intersecting a horizontal edge
        //              (2 x 2 unique points --> 2 points is even)
    
        bool is_border = false;
        std::unordered_set<Point, Point::hash> unique_intersections;        
                
        Ray ray;
        ray.p = p;
        
        size_t n = vertices.size();
        
        for (size_t i1 = 0; i1 < n; ++i1) {
        
            size_t i2 = (i1 + 1) % n; // circular list! 
            
            Segment seg(vertices[i1], vertices[i2]);
            
            std::string loc = seg.point_location(p);

            // handles degenerate cases of border points
            if (loc == "ON_SEGMENT") {
                is_border = true;
                break;
            }
            
            std::vector<Point> inters = ray.intersections(seg);
            for (Point const& p : inters) {
                //debug
                // std::cout << ray << " intersects " << seg << std::endl;
                unique_intersections.emplace(p);
            }
        }
    
        if (is_border)
            return "BORDER";
        
        // N(int) is odd -> inside, if N(int) is even -> outside
        if (unique_intersections.size() % 2 != 0) 
            return "INSIDE";
            
        return "OUTSIDE";
    }


};

std::ostream& operator<<(std::ostream& os, Polygon p) {
    os << "polygon (";
    for (size_t i = 0; i < p.vertices.size(); ++i) {
        os << p.vertices[i];
        if (i < p.vertices.size() - 1)
            os << ",";
    }
    os << ")";
    return os;
}



} 


#endif