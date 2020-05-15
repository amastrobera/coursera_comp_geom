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
#include <sstream>
#include <vector> // main container 
#include <unordered_set> // non cvx function
#include <utility> 



namespace crs {

// SIMPLE ALGORITHMS

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
    
    double distance(Point const& p) const {
        return std::sqrt(std::pow(p.x - x, 2) + std::pow(p.y - y, 2)); 
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
        
     
                    bool found = false;
        
                    for (size_t i = 0; i < n + 1; ++i) {
                        
                        size_t idx0 = (i - 1) % n;
                        size_t idx1 = (i) % n;
                        size_t idx2 = (i + 1) % n;
                        
                        if (Segment(vertices[idx1], p).point_location(vertices[idx0]) != "RIGHT" 
                            && 
                            Segment(vertices[idx1], p).point_location(vertices[idx2]) != "RIGHT") {
                            // found first (right-most tangent v[i]-p)
                            
                            // find the farmost tangent : search for points on line                            
                            size_t di = 1;
                            while (Segment(vertices[(idx1 - di) % n], vertices[idx1])
                                    .point_location(p) == "ON_LINE")
                                ++di;
                            
                            tangents.first = vertices[(idx1 - di + 1) % n];
                            
                            for (size_t j = 1; j < n + 1; ++j) {
                                idx0 = (i + j - 1) % n;
                                idx1 = (i + j) % n;
                                idx2 = (i + j + 1) % n;

                                if (Segment(p, vertices[idx1]).point_location(vertices[idx0]) != "RIGHT" 
                                    && 
                                    Segment(p, vertices[idx1]).point_location(vertices[idx2]) == "LEFT") {
                                    // found second (right-most tangent p-v[j])
                                    
                                    // find the farmost tangent : search for points on line                            
                                    size_t di = 1;
                                    while (Segment(p, vertices[idx1]).point_location(
                                            vertices[(idx1 + di) % n]) == "ON_LINE")
                                        ++di;
                                    found = true;
                                    tangents.second = vertices[(idx1 + di -1) % n];
                                    break;
                                }
                            }
                            if (found)
                                break;
                        }
                    }
 
                } else {
                    // try to create a counter clockwise order triangle
                    //  from a segment ab and the point p
                    std::string loc = Segment(
                                    vertices[0], vertices[1]).point_location(p);
                    if (loc == "LEFT") {
                        tangents = {vertices[1], vertices[0]};
                    } else if (loc == "RIGHT" || loc == "ON_SEGMENT") {
                        tangents = {vertices[0], vertices[1]};
                    } else {
                        if (p.distance(vertices[0]) < p.distance(vertices[1])) {
                            tangents = {vertices[1], vertices[0]};
                        } else {
                            tangents = {vertices[0], vertices[1]};
                        }
                    }
                    tangents = {vertices[0], vertices[1]};
                }
                    
            } else {
                // create a segment from a point p and the only vertex 
                tangents = {vertices[0], vertices[0]};
            }
                
        }
        
        // std::cout <<"   tangents = " << tangents.first << " " << tangents.second << std::endl;
        
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
        // pick the lowest point
        size_t i0 = 0;
        double ymin = points[0].y;
        for (size_t i = 1; i < n; ++i)
            if (points[i].y < ymin) {
                ymin = points[i].y;
                i0 = i;
            }

        // initialize with the smallest point and the point after
        std::vector<Point> cvpoints;
        cvpoints.push_back(points[i0 % n]);
        cvpoints.push_back(points[(i0+1) % n]);
        
        for (size_t i = 2; i < n + 1 ; ++i) {            
            cvpoints.push_back(points[(i0+i) % n]);
            
            size_t m = cvpoints.size();
            size_t i3 = m-1;
            size_t i2 = i3-1;
            size_t i1 = i3-2;            

            while (i1 >= 0 && 
                   Segment(cvpoints[i1], cvpoints[i2])
                        .point_location(cvpoints[i3]) != "LEFT") {
                cvpoints[i2] = cvpoints[i3];
                cvpoints.pop_back();
                --i1;
                --i2;
                --i3;
            }            
        }
        
        // end condition
        if (*cvpoints.begin() == *(cvpoints.end()-1))
            cvpoints.erase(cvpoints.end()-1);        

        cvpoly.vertices = std::move(cvpoints);
        
    } else
        for (auto it = points.begin(); it!= points.end(); ++it)
            cvpoly.vertices.push_back(*it);
    
    
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

    if (!(numverts >= 3 && numverts <= 1000)) {
        cerr << "bad input n vertices " 
             << numverts << " not in [3, 1000]" << endl;
        exit(1);
    }

    // get vertices  
    string strPoly;
    getline(cin, strPoly);
    istringstream streamPoly(strPoly);

    crs::Polygon poly;

    while (numverts--) {

        crs::Point point;
        string strCoord;
        
        getline(streamPoly, strCoord, ' ');
        istringstream(strCoord) >> point.x;
        
        getline(streamPoly, strCoord, ' ');
        istringstream(strCoord) >> point.y;

        poly.vertices.push_back(point);
        
    }
    
    if (!poly.is_convex()) {    
        cerr << "bad input, not convex\n " << poly << endl;
        exit(1);
    }
    
    // find how many vertices there are
    int numpoints;
    string inputNumPoints;
    getline(cin, inputNumPoints);
    istringstream(inputNumPoints) >> numpoints;

    if (!(numpoints >= 1 && numpoints <= 1000)) {
        cerr << "bad input n points " 
             << numpoints << " not in [1, 1000]" << endl;
        exit(1);
    }

    while (numpoints--) {

        crs::Point point;

        string strLineCoords;
        getline(cin, strLineCoords);
        istringstream streamCoords(strLineCoords);
        string strCoord;
       
        getline(streamCoords, strCoord, ' ');
        istringstream(strCoord) >> point.x;
        
        getline(streamCoords, strCoord, ' ');
        istringstream(strCoord) >> point.y;
        
        if (poly.point_location(point) == "INSIDE") {
            cerr << "invalid input (" << point << ") is inside the polygon" << endl;
            continue;
        }
        
        auto tangents = poly.tangent_points_to_point(point);

        cout << tangents.first.x << " " << tangents.first.y << " "
             << tangents.second.x << " " << tangents.second.y << endl;
 
    }
  
    return 0;
}
