/*
c++ std: 11
compiler: g++ v 7.4.0

compile with 
g++ -I.  main.cpp

launch with
cat data.txt | ./a.out 
*/


#include <algorithm> //sort etc
#include <cmath> // atan, PI etc
#include <cstdio> // exit
#include <iostream> // error msg 
#include <iterator> // next, prev
#include <sstream>
#include <vector> // main container 
#include <unordered_set> // non cvx function
#include <utility> 


namespace crs {

// POINT
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



// SEGMENT

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


// RAY

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


// POLYGON
struct Polygon {

    std::vector<Point> vertices;

    Polygon()  {} 
    Polygon& operator=(Polygon const& cp) {
        if (this != &cp) {
            vertices = cp.vertices;
        }
        return *this;
    }
    Polygon(Polygon const& cp) {
        *this = cp;
    }

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



// ALGORITHMS 


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




int main() {

    using namespace std;

    int numpoly;
    string inputNumPolys;
    getline(cin, inputNumPolys);
    istringstream(inputNumPolys) >> numpoly;

    if (!(numpoly >= 2)) {
        cerr << "bad input n polygons " 
             << numpoly << " not in [2, inf[" << endl;
        exit(1);
    }

    std::vector<crs::Polygon> polys;
    
    while (numpoly--) {
        
        int numverts;
        string inputNumVerts;
        getline(cin, inputNumVerts);
        istringstream(inputNumVerts) >> numverts;

        if (!(numverts >= 3 && numverts <= 100000)) {
            cerr << "bad input (" << numverts << ") n vertices " 
                 << numverts << " not in [3, 100,000]" << endl;
            exit(1);
        }

        // get vertices  
        string strPoly;
        getline(cin, strPoly);
        istringstream streamPoly(strPoly);

        crs::Polygon poly;

        double input_limit = 1000000;

        while (numverts--) {

            crs::Point point;
            string strCoord;
            
            getline(streamPoly, strCoord, ' ');
            istringstream(strCoord) >> point.x;

            getline(streamPoly, strCoord, ' ');
            istringstream(strCoord) >> point.y;

            if (abs(point.x) > input_limit || abs(point.y) > input_limit) {
                cerr << "bad input (" << numpoly << ") coordinate ("
                     << numverts << ") point exceeds limit " << input_limit << endl;
                exit(1);
            }

            poly.vertices.push_back(point);
        }
        
        if (!poly.is_convex()) {
            cerr << "bad input (" << numpoly 
                 << ") not convex " << poly << endl;
            exit(1);
        }
        
        polys.push_back(poly);
    }

    crs::Polygon mergepoly = crs::merge_cvx_poly(polys);
    
    if (!mergepoly.is_convex()) {    
        cerr << "bad output, not convex\n " << mergepoly << endl;
        exit(1);
    }

    size_t n = mergepoly.vertices.size();    
    cout << n << endl;
    for (size_t i = 0; i < n; ++i) {
        cout << mergepoly.vertices[i].x << " " << mergepoly.vertices[i].y;
        if (i < n -1)
            cout << " ";
        else
            cout << endl;
    }

  
    return 0;
}
