/*
c++ std: 11
compiler: g++ v 7.4.0

compile with 
g++ -I.  main.cpp

launch with
cat data.txt | ./a.out 
*/

#include <algorithm> // sort
#include <cmath>
#include <cstdio> // exit
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


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



struct Segment {

    Point a, b;
    
    Segment(Point const& p1, Point const& p2) : a(p1), b(p2) {}

    std::string point_location(Point c) {
        
        double d = det(a, b, c);
        
        if (d > 0)
            return "LEFT";
        if (d < 0)
            return "RIGHT";
        
        // determinant is zero. check if the point is on the line of on segment
        //      according to the location of AB
        // I want to simply check xmin <= c.x <= xmax AND ymin <= c.y <= ymax
        //      so I determine the min/max for AB, then apply the check
        double xmin, xmax;
        if (b.x >= a.x) {
            xmin = a.x;
            xmax = b.x;
        } else {
            xmin = b.x;
            xmax = a.x;
        }
        double ymin, ymax;
        if (b.y >= a.y) {
            ymin = a.y;
            ymax = b.y;
        } else {
            ymin = b.y;
            ymax = a.y;
        }
        
        if (c.x >= xmin && c.x <= xmax && 
            c.y >= ymin && c.y <= ymax)
            return "ON_SEGMENT";
        return "ON_LINE";
    }
    
    double length() {
        return std::sqrt(std::pow(b.x-a.x, 2) + std::pow(b.y-a.y, 2));
    }

private:
    double det(Point a, Point b, Point c) {
        return (b.x - a.x)*(c.y- a.y) - (c.x - a.x)*(b.y - a.y);
    }
};

bool operator==(Segment const& s1, Segment const& s2) {
    // if (a,b) are coincident the segments are equal
    //  note: I am not considering that (a,b) == (b,a) 
    //        equal size but different direction
    return s1.a.x == s2.a.x && s1.a.y == s2.a.y &&  
           s1.b.x == s2.b.x && s1.b.y == s2.b.y;
}

bool operator!=(Segment const& s1, Segment const& s2) {
    // if (a,b) are coincident the segments are equal
    //  note: I am not considering that (a,b) == (b,a) 
    //        equal size but different direction
    return !(s1 == s2);
}


std::ostream& operator<<(std::ostream& os, Segment s) {
    os << "line (" << s.a << "," << s.b  << ")";
    return os;
}




struct Polygon {

    std::vector<Point> vertices;

    std::string type() {
        if (is_convex())
            return "CONVEX";
        return "NOT_CONVEX";
    }

    bool is_convex() {
        // checks that the polygon has at least 3 edges
        //      and those are in counter-clockwise

        bool is_ccw = true;
        size_t n = vertices.size();
        
        for (size_t i1 = 0; i1 < n; ++i1) {
        
            size_t i2 = (i1 + 1) % n; // circular list! 
            size_t i3 = (i1 + 2) % n;
            
            Segment seg(vertices[i1], vertices[i2]);
            
            if (seg.point_location(vertices[i3]) == "LEFT")
                is_ccw &= true;
            else {
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
    
        bool is_inside = true;
        bool is_border = false;
        
        size_t n = vertices.size();
        
        for (size_t i1 = 0; i1 < n; ++i1) {
        
            size_t i2 = (i1 + 1) % n; // circular list! 
            
            Segment seg(vertices[i1], vertices[i2]);
            
            std::string loc = seg.point_location(p);
            
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


Polygon convex_hull(std::vector<Point> points) {


    Polygon cvpoly;
    
    // sorted in DESC order (max .... min)
    std::sort(points.begin(), points.end(), 
            [](Point const& p1, Point const& p2) {
                return p1.x > p2.x;
            });
    
    std::vector<Point>& cvpoints = cvpoly.vertices;
    
    size_t n = points.size();
    if (n) 
        cvpoints.push_back(points[0]);

    if (n > 1) 
        cvpoints.push_back(points[1]);
        
    if (n > 2) {
        // UpperHull(P)
        for (size_t i = 2; i < n; ++i) {
            size_t i1 = i - 2;
            size_t i2 = i - 1;
            size_t i3 = i;
            cvpoints.push_back(points[i3]);
            
            Segment seg(cvpoints[i1], cvpoints[i2]);
            std::string loc = seg.point_location(cvpoints[i3]);
            while (i1 > 0 && (loc == "RIGHT" || loc == "ON_LINE")) {
                
                cvpoints[i2] = cvpoints[i2+1];
                cvpoints.erase(cvpoints.end()-1);
                
                --i1;
                --i2;
                --i3;

                seg.a = cvpoints[i1];
                seg.b = cvpoints[i2];
                loc = seg.point_location(cvpoints[i3]);

            }
            
            // end condition of UpperHull(P): reached xmin
            if (*(cvpoints.end()-1) == *(points.end()-1))
                break;
        }
        
        // LowerHull(P)
        cvpoints.push_back(points[points.size()-2]);
        for (size_t i = n-3; i >= 0; --i) {
            size_t i1 = i + 2;
            size_t i2 = i + 1;
            size_t i3 = i;
            cvpoints.push_back(points[i3]);
            
            Segment seg(cvpoints[i1], cvpoints[i2]);
            std::string loc = seg.point_location(cvpoints[i3]);
            while (i1 < n && (loc == "RIGHT" || loc == "ON_LINE")) {
                
                cvpoints[i2] = cvpoints[i2+1];
                cvpoints.erase(cvpoints.end()-1);
                
                ++i1;
                ++i2;
                ++i3;

                seg.a = cvpoints[i1];
                seg.b = cvpoints[i2];
                loc = seg.point_location(cvpoints[i3]);

            }
            
            // end condition of LowerHull(P): reached xmin
            if (*(cvpoints.end()-1) == *(points.begin()))
                break;
        }
    }
    
    return std::move(cvpoly);
}



} 


int main() {

    using namespace std;


    // find how many vertices there are
    int numverts;
    string inputNumVerts;
    getline(cin, inputNumVerts);
    istringstream(inputNumVerts) >> numverts;

    if (!(numverts >= 3 && numverts <= 5000)) {
        cerr << "bad input n points " 
             << numverts << " not in [3, 5000]" << endl;
        exit(1);
    }

    // get vertices  
    string strPoly;
    getline(cin, strPoly);
    istringstream streamPoly(strPoly);

    std::vector<crs::Point> points;

    bool input_error = false;
    double input_limit = 10^6;

    while (numverts--) {

        crs::Point point;
        string strCoord;
        
        getline(streamPoly, strCoord, ' ');
        istringstream(strCoord) >> point.x;
        
        getline(streamPoly, strCoord, ' ');
        istringstream(strCoord) >> point.y;
        
        if (abs(point.x) > input_limit || abs(point.y) > input_limit) {
            input_error = true;
            break;
        }             
        
        points.push_back(point);
        
    }
    
    // debug
    cout << "input" << endl;
    size_t m = points.size();
    cout << m << endl;
    for (size_t i = 0; i < m; ++i) {
        cout << points[i].x << " " << points[i].y;
        if (i < m -1)
            cout << " ";
        else
            cout << endl;
    }
    
    
    if (input_error) {    
        cerr << "bad input, invalid point exceeds limit " << input_limit << endl;
        exit(1);
    }
    

    crs::Polygon cvpoly = crs::convex_hull(points);
    
    size_t n = cvpoly.vertices.size();
    
    // debug
    cout << "output" << endl;
    
    cout << n << endl;
    for (size_t i = 0; i < n; ++i) {
        cout << cvpoly.vertices[i].x << " " << cvpoly.vertices[i].y;
        if (i < n -1)
            cout << " ";
        else
            cout << endl;
    }

    

    return 0;
}
