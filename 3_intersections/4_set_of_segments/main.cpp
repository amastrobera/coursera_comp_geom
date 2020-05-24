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
#include <unordered_set>
#include <utility> 


namespace crs {

// POINT
struct Point {
    int x,y;
    
    Point() : x(0), y(0) {} 
    Point(int c1, int c2) : x(c1), y(c2) { }
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
    ~Point() = default;
    Point(Point && mv) = default;
    
    double distance(Point const& p) const {
        return std::sqrt(std::pow(p.x - x, 2) + std::pow(p.y - y, 2)); 
    }
    
    
    struct hash{
        size_t operator()(const Point &p) const {
            size_t h1 = std::hash<int>()(p.x);
            size_t h2 = std::hash<int>()(p.y);
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


// RANGE
struct Range {
    int xmin, xmax, ymin, ymax;
};

std::ostream& operator<<(std::ostream& os, Range rng) {
    os << "[" << rng.xmin << ", " << rng.xmax << "] x "
       << "[" << rng.ymin << ", " << rng.ymax << "]";
    return os;
}



// SEGMENT

struct Segment {

    Point p1, p2;
    
    Segment() {}
    Segment(Point const& cp1, Point const& cp2) : p1(cp1), p2(cp2) {}
    Segment& operator=(Segment const& cp) {
        if (this != &cp) {
            p1 = cp.p1;
            p2 = cp.p2;
        }
        return *this;
    }
    Segment(Segment const& cp) {
        *this = cp;
    }


    std::string point_location(Point const& c) const {
        
        double d = det(p1, p2, c);
        
        if (d > 0)
            return "LEFT";
        if (d < 0)
            return "RIGHT";
        
        // determinant is zero. check if the point is on the line of on segment
        //      according to the location of AB
        // I want to simply check xmin <= c.x <= xmax AND ymin <= c.y <= ymax
        //      so I determine the min/max for AB, then apply the check
        double xmin, xmax;
        if (p2.x >= p1.x) {
            xmin = p1.x;
            xmax = p2.x;
        } else {
            xmin = p2.x;
            xmax = p1.x;
        }
        double ymin, ymax;
        if (p2.y >= p1.y) {
            ymin = p1.y;
            ymax = p2.y;
        } else {
            ymin = p2.y;
            ymax = p1.y;
        }
        
        if (c.x >= xmin && c.x <= xmax && 
            c.y >= ymin && c.y <= ymax)
            return "ON_SEGMENT";
        return "ON_LINE";
    }
    
    double length() const {
        return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
    }
    
    
    std::string intersect(Segment const& seg2) const {
        
        // step 1: make system of 2 (implicit) equations + 2 unknowns
        //      implicit a*x + b*y = c, for each segment
        //           a = y1 - y2 (1 = first point, 2 = second point)
        //           b = x2 - x1
        //           c = y1*x2 - x1*y2

        // debug only
        // std::cout << "s1:" << p1 << " -> " << p2
                  // << ", s2: " << seg2.p1 << " -> " << seg2.p2 << "; ";

        double a1 = p1.y - p2.y;
        double b1 = p2.x - p1.x;
        double c1 = p1.y * p2.x - p1.x * p2.y;

        double a2 = seg2.p1.y - seg2.p2.y;
        double b2 = seg2.p2.x - seg2.p1.x;
        double c2 = seg2.p1.y * seg2.p2.x - seg2.p1.x * seg2.p2.y;

        // step 2: check for parallelism
        //       d = | a1  b1 |
        //           | a2  b2 |
         
        double d = a1*b2 - a2*b1;
        if (d == 0) {
            // lines are parallel

            // step 2.3: check if a segment is inside another
            if (point_location(seg2.p1) == "ON_SEGMENT" || 
                point_location(seg2.p2) == "ON_SEGMENT" ) {
                return "A common segment of non-zero length.";
            }
            
            return "No common points.";
            
        } else {
            
            // step 3: find the intersection of the lines
            //         Cramer law, here 
            //       x = | c1  b1 | / d    y = | a1  c1 | / d
            //           | c2  b2 |            | a2  c2 |

            double xi = (c1*b2 - c2*b1) / d;
            double yi = (a1*c2 - a2*c1) / d;
            Point inter(xi, yi);
            // debug only
            // std::cout << "xi= " << inter.x << " yi=" << inter.y << "  ";
            
            // step 3.1: check that the intersection is inside the domain
            if (     point_location(inter) == "ON_SEGMENT" && 
                seg2.point_location(inter) == "ON_SEGMENT") {
                std::ostringstream os;
                os << "The intersection point is (" 
                   << inter.x << ", " << inter.y << ").";
                return os.str();                
            }
            
            //           ... if it's not, no intersection 
            return "No common points.";
        }
        
        return "";
    }
    

private:
    double det(Point const& a, Point const& b, Point const& c) const {
        return (b.x - a.x)*(c.y- a.y) - (c.x - a.x)*(b.y - a.y);
    }

    Range domain() {
        // returns the domain of the segment
        Range r;
        r.xmin = (p1.x < p2.x) ? p1.x : p2.x;
        r.xmax = (p1.x < p2.x) ? p2.x : p1.x;
        r.ymin = (p2.y < p1.y) ? p2.y: p1.y;
        r.ymax = (p2.y < p1.y) ? p1.y: p2.y;
        return std::move(r);
    }

};

bool operator==(Segment const& s1, Segment const& s2) {
    // if (a,b) are coincident the segments are equal
    //  note: I am not considering that (a,b) == (b,a) 
    //        equal size but different direction
    return s1.p1.x == s2.p1.x && s1.p1.y == s2.p1.y &&  
           s1.p2.x == s2.p2.x && s1.p2.y == s2.p2.y;
}

bool operator!=(Segment const& s1, Segment const& s2) {
    // if (a,b) are coincident the segments are equal
    //  note: I am not considering that (a,b) == (b,a) 
    //        equal size but different direction
    return !(s1 == s2);
}


std::ostream& operator<<(std::ostream& os, Segment s) {
    os << "line (" << s.p1 << "," << s.p2  << ")";
    return os;
}


// RAY
struct Ray {
    // a (horizontal, left) ray passing by point p
    Point p;
    
    bool intersects(Segment const& s) const {

        int ymax = (s.p2.y >= s.p1.y) ? s.p2.y : s.p1.y;
        int ymin = (s.p2.y < s.p1.y) ? s.p2.y : s.p1.y;
        
        if (p.y >= ymin && p.y <= ymax) {

            // horizontal segment
            if (s.p1.y == s.p2.y) {
                int xmin = (s.p1.x < s.p2.x) ? s.p1.x : s.p2.x;
                if (p.x >= xmin)
                    return true;
                return false;
            }
            
            // vertical segment
            if (s.p1.x == s.p2.x) {
                if (p.x >= s.p1.x)
                    return true;
                return false;
            }
            
            // diagonal segment. make sure it points up, decide whether
            // it point of the ray is on the right (or on the segment)
            Segment scmp;            
            if (s.p1.y < s.p2.y) {
                scmp.p1 = s.p1;
                scmp.p2 = s.p2;
            } else {
                scmp.p1 = s.p2;
                scmp.p2 = s.p1;
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

        int ymax = (s.p2.y >= s.p1.y) ? s.p2.y : s.p1.y;
        int ymin = (s.p2.y < s.p1.y) ? s.p2.y : s.p1.y;
        
        if (p.y >= ymin && p.y <= ymax) {

            // horizontal segment
            if (s.p1.y == s.p2.y) {
                int xmin = (s.p1.x < s.p2.x) ? s.p1.x : s.p2.x;
                if (p.x >= s.p1.x)
                    int_points.push_back(s.p1);
                if (p.x >= s.p2.x)
                    int_points.push_back(s.p2);
            } else {
                // vertical segment
                if (s.p1.x == s.p2.x) {
                    if (p.x >= s.p1.x)
                        int_points.push_back(Point(s.p1.x, p.y));
                } else {            
                    // diagonal segment. make sure it points up, decide whether
                    // it point of the ray is on the right (or on the segment)
                    Segment scmp;            
                    if (s.p1.y < s.p2.y) {
                        scmp.p1 = s.p1;
                        scmp.p2 = s.p2;
                    } else {
                        scmp.p1 = s.p2;
                        scmp.p2 = s.p1;
                    }
                    std::string location = scmp.point_location(p);
                    
                    if (location == "RIGHT" || location == "ON_SEGMENT" ) {
                        double m = (scmp.p2.y - scmp.p1.y) / (scmp.p2.x - scmp.p1.x);
                        double x = scmp.p1.x + (p.y - scmp.p1.y) / m;
                        int_points.push_back(Point(x, p.y)); // will be casted into int
                        // it's so ugly to use "ints" instead of doubles for x,y
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

    Polygon() = default;
    ~Polygon() = default;
    Polygon& operator=(Polygon const& cp) {
        if (this != &cp) {
            vertices = cp.vertices;
        }
        return *this;
    }
    Polygon(Polygon const& cp) {
        *this = cp;
    }
    
    Polygon (Polygon && mv) = default;

    std::string type() {
        if (is_convex())
            return "CONVEX";
        return "NOT_CONVEX";
    }

    bool is_convex() const {
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

    bool is_valid() const {
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


    std::string point_location(Point const& p) const {
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

    std::string point_location_cvx(Point const& p) const {
        
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
    
    std::string point_location_not_cvx(Point const& p) const {
    
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

        // sort ASC by angle against Z mean point
        //      precision to 9th decimal 
        //      (sort ints, not doubles)
        double xavg = 0, yavg = 0;
        for (Point const& p: points) {
            xavg += p.x;
            yavg += p.y;
        }
        xavg /= n;
        yavg /= n;
        
        std::sort(points.begin(), points.end(), 
                    [xavg, yavg](Point const& p1, Point const& p2) {
                        return (int)1000000000*point_angle_2d(p1.x - xavg, p1.y - yavg) < 
                               (int)1000000000*point_angle_2d(p2.x - xavg, p2.y - yavg); 
                     });

        // pick the lowest point
        size_t i0 = 0;
        int ymin = points[0].y;
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
            // int used to check for non-negativity in the loop
            int i3 = m - 1;
            int i2 = i3 - 1;
            int i1 = i3 - 2;            
            
            // std::cout << "cvpoints(" << i << "/" << n << ") = [" 
                      // << cvpoints[i1] << ", " 
                      // << cvpoints[i2] << ", " 
                      // << cvpoints[i3] << "]"
                      // << std::endl;

            while (i1 >= 0 && 
                   Segment(cvpoints[i1], cvpoints[i2])
                        .point_location(cvpoints[i3]) != "LEFT") {

                // std::cout << "... pop(" << i2 << ") = " << cvpoints[i2] << std::endl;
                
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
        
    } else {
        for (auto it = points.begin(); it!= points.end(); ++it)
            cvpoly.vertices.push_back(*it);
    }

    return std::move(cvpoly);
}


Polygon merge_cvx_poly(Polygon const& poly1, Polygon const& poly2) {
    // takes 2 convex polygons and merges them into another convex polygon

    // two ways of doing this:
    //  a. find the tangents, remove inner points (how ??? )
    //  b. transform all into points and make again convex hull (easier but slower ?)
    
    // TODO: can we try (a) ? 
    
    // applied solution (b)
    std::vector<Point> points;
    for (Point const& p : poly1.vertices) 
        points.push_back(p);
    for (Point const& p : poly2.vertices) 
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
            
            std::vector<Point> points;
            for (Polygon const& poly : polys) {
                for (Point const& p : poly.vertices) {
                    points.push_back(p);
                }
            }
            merged = convex_hull(points);

        } else
            for (Point const& p : polys[0].vertices)
                merged.vertices.push_back(p);
    }
    return std::move(merged);
}


size_t sweep_intersections_number(std::vector<Segment>) {
    
    
    return 0;
}


std::vector<Segment> sweep_intersection_lines(std::vector<Segment>) {
    
    std::vector<Segment> inters;
    
    
    return std::move(inters);
}

} 





int main() {

    using namespace std;

    double limit_coord = 1000000;

    vector<crs::Segment> segments;
    
    
    size_t n;
    string numString;
    getline(cin, numString);
    istringstream(numString) >> n;
    
    
    if (n < 1 || n > 50000)  {
        cerr << "bad input coordinates: n (" << n 
             << ") not in range [1, 50,000]" << endl;
        exit(1);
    }
    

    while (n--) {
        
        crs::Segment s;
        
        string segString;
        getline(cin, segString);
        istringstream segStream(segString);
        
        string coordString;
        getline(segStream, coordString, ' ');
        istringstream(coordString) >> s.p1.x;
        getline(segStream, coordString, ' ');
        istringstream(coordString) >> s.p1.y;
        
        getline(segStream, coordString, ' ');
        istringstream(coordString) >> s.p2.x;
        getline(segStream, coordString, ' ');
        istringstream(coordString) >> s.p2.y;
        
        if (abs(s.p1.x) > limit_coord || abs(s.p1.y) > limit_coord ||
            abs(s.p2.x) > limit_coord || abs(s.p2.y) > limit_coord ) {
            cerr << "bad input coordinates: " << coordString 
                 << " (|coordinate| > " << limit_coord << ")" << endl;
            exit(1);
        }
        
        segments.push_back(s);
    }


    vector<crs::Segment> inters = crs::sweep_intersection_lines(segments);


    size_t m = inters.size();
    cout << m << endl;
    
    for (size_t i = 0; i < m; ++i) {
        
        cout << inters[i].p1.x << " "
             << inters[i].p1.y << " "
             << inters[i].p2.x << " "
             << inters[i].p2.y << endl;

    }
  
    return 0;
}
