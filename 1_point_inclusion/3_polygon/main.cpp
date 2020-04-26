/*
c++ std: 11
compiler: g++ v 7.4.0

compile with 
g++ -I.  main.cpp

launch with
cat data.txt | ./a.out 
*/

#include <cmath>
#include <cstdio> // exit
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>


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



struct Segment {

    Point a, b;
    
    Segment() {}
    Segment(Point const& p1, Point const& p2) : a(p1), b(p2) {}
    Segment& operator=(Segment const& cp) {
        if (this != &cp) {
            a = cp.a;
            b = cp.b;
        }
        return *this;
    }
    Segment(Segment const& cp) {
        *this = cp;
    }


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





struct Polygon {

    std::vector<Point> vertices;

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


int main() {

    using namespace std;


    // find how many vertices there are
    int numverts;
    string inputNumVerts;
    getline(cin, inputNumVerts);
    istringstream(inputNumVerts) >> numverts;

    if (!(numverts >= 3 && numverts <= 200)) {
        cerr << "bad input n(vertices) " 
             << numverts << " not in [3, 200]" << endl;
        exit(1);
    }

    // get vertices  
    string strPoly;
    getline(cin, strPoly);
    istringstream streamPoly(strPoly);

    crs::Polygon poly;
    
    while (numverts--) {

        crs::Point vert;
        string strCoord;
        
        getline(streamPoly, strCoord, ' ');
        istringstream(strCoord) >> vert.x;
        
        getline(streamPoly, strCoord, ' ');
        istringstream(strCoord) >> vert.y;

        poly.vertices.push_back(vert);
        
    }


    if (!poly.is_valid()) {
        cerr << "bad input, invalid " << poly << endl;
        exit(1);
    }

    int numpoints;
    string inputNumPoints;
    getline(cin, inputNumPoints);
    istringstream(inputNumPoints) >> numpoints;

    if (!(numpoints >= 1 && numpoints <= 200)) {
        cerr << "bad input n(points) " 
             << numpoints << " not in [1, 200]" << endl;
         exit(1);                
    }
    
    
    double input_limit = 10000;
    
    while (numpoints--) {

        string inputpoint;
        getline(cin, inputpoint);
        istringstream ptStream(inputpoint);

        // put segment data into Segment class from input stream
        string inputxy;
        crs::Point point;
        getline(ptStream, inputxy, ' ');
        istringstream(inputxy) >> point.x;
        getline(ptStream, inputxy, ' ');
        istringstream(inputxy) >> point.y;


        if (abs(point.x) > input_limit || abs(point.y) > input_limit) {
            cerr << "bad input, invalid entry, " 
                 << point << " exceeds limit " << input_limit << endl;
            continue;
        }

        string location = poly.point_location(point);

        // debug
        // cout << point << " is " << location << endl;

        // result
        cout << location << endl;
    
    }
    

    return 0;
}
