#ifndef ALGOS_H
#define ALGOS_H


#include "point.hpp"
#include "segment.hpp"
#include "polygon.hpp"

#include <algorithm> //sort etc
#include <cmath> // atan, PI etc
#include <iostream> // error msg 
#include <iterator> // next, prev
#include <vector> // main container 


namespace crs {

#define _USE_MATH_DEFINES

double point_angle_2d(double x, double y) {
    /*
        input:     x, y:
        output:    angle (radians) of the line (make of two points)
                        (the first point is the origin, and the second is the end of the radius)
        example:
                      *         |        * (q)
               135o      .      |     .         45o
                            .   |  .
                   -------- (p) * ----------->
               225o         .   |  .
                         .      |     .       315o
                      *         |       *
    */
    if (x == 0)
        return (y > 0) ? M_PI/2 : ((y < 0) ? M_PI*3/2 : 0); 
    else if (y == 0)
        return (x < 0) ? M_PI : 0;

    if (x > 0)
        if (y > 0)
            return std::atan(y/x);
        else
            return 2*M_PI - std::atan(-y/x);
    else
        if (y > 0)
            return M_PI/2. + std::atan(-x/y);
        else
            return M_PI + std::atan(y/x);
    return 0; 
}



Polygon convex_hull(std::vector<Point> points) {
    // Graham's scan
    Polygon cvpoly;
    
    size_t n = points.size();    
    if (n > 2) {

        // create z = mid point 2d
        double xavg = 0;
        double yavg = 0;
        for (auto it = points.begin(); it != points.end(); ++it) {
            xavg += it->x;
            yavg += it->y;
        }
        xavg /= n;
        yavg /= n;
        Point z(xavg, yavg);
    
        // sort ASC by angle against z
        std::sort(points.begin(), points.end(), 
                    [z](Point const& p1, Point const& p2) {
                        return point_angle_2d(p1.x-z.x, p1.y-z.y) < 
                               point_angle_2d(p2.x-z.x, p2.y-z.y); 
                     });


        std::vector<Point> cvpoints;
        cvpoints.push_back(*(points.end()-1));  // point with max angle
        cvpoints.push_back(*(points.begin()));  // followed by point with min angle
        
        for (std::vector<Point>::const_iterator it = points.begin()+1;
                it != points.end(); ++it) {
            
            cvpoints.push_back(*it);

            auto it3 = std::prev(cvpoints.end(), 1);
            auto it2 = std::prev(it3, 1);
            auto it1 = std::prev(it3, 2);
            
            Segment seg(*it1, *it2);
            std::string loc = seg.point_location(*it3);
            
            while (it1 != std::prev(cvpoints.begin(), 1) && 
                    (loc == "RIGHT" || loc == "ON_LINE")) {
                
                *it2 = *it3;
                cvpoints.erase(std::next(it2, 1));
                
                --it1;
                --it2;
                --it3;

                seg.a = *it1;
                seg.b = *it2;
                loc = seg.point_location(*it3);
            }            
        }

        if (*(cvpoints.end()-1) == *(cvpoints.begin()))
            cvpoints.erase(cvpoints.end()-1);

        cvpoly.vertices = std::move(cvpoints);
        
    } else
        for (auto it = points.begin(); it!= points.end(); ++it)
            cvpoly.vertices.push_back(*it);
    
    
    return std::move(cvpoly);
}


} 



#endif