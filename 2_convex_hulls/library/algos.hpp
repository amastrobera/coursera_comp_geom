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
#include <unordered_set>


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


Polygon merge_cvx_poly(Polygon const& poly1, Polygon const& poly2) {
    // takes 2 convex polygons and merges them into another convex polygon

    // two ways of doing this:
    //  a. find the tangents, remove inner points
    //  b. transform all into points and make again convex hull (easier but slower ?)
    
    // TODO: (a)
    
    // temporary solution (b)
    std::vector<Point> points;
    std::unordered_set<Point, Point::hash> unique_points;
    for (Point const& p : poly1.vertices) 
        unique_points.insert(p);
    for (Point const& p : poly2.vertices) 
        unique_points.insert(p);
    for (Point const& p : unique_points) 
        points.push_back(p);
    Polygon merged = convex_hull(points);
    
    
    return std::move(merged);
}


Polygon merge_cvx_poly(std::vector<Polygon> const& polys) {
    // takes 1+ convex polygons and merges them into another convex polygon
    Polygon merged;
    
    size_t npoly = polys.size();
    if (npoly) {
        if (npoly > 1) {
            Polygon semimerged = merge_cvx_poly(polys[0], polys[1]);
            for (size_t i = 2; i < npoly; ++i) {
                Polygon tmp = merge_cvx_poly(semimerged, polys[i]);
                semimerged = std::move(tmp);
            }
            merged = std::move(semimerged);            
        } else
            merged = polys[0];
    }
    
    return std::move(merged);
}


} 



#endif