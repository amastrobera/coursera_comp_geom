/*
c++ std: any
compiler: g++ v 7.4.0

compile with 
g++ -I.  main.cpp

launch with
cat data.txt | ./a.out 
*/

#include <iostream>
#include <sstream>
#include <string>


namespace crs {

struct Point {
    double x,y;
    
    Point& operator=(Point const& cp) {
        if (this != &cp) {
            x = cp.x;
            y = cp.y;
        }
        return *this;
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
    // a (horizontal, left) ray passing by point a
    Point p;
    
    bool intersects(Segment s) {

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
};


std::ostream& operator<<(std::ostream& os, Ray r) {
    os << "ray (" << r.p << ")";
    return os;
}




struct Triangle {

    Segment s1, s2, s3;


    bool is_valid() {
        // checks that the triangle is a closed set of points
        //      (does not check if segments are in counter-clockwise)
        return    s1 != s2 && s2 != s3 && s3 != s1 
               && s1.b == s2.a && s2.b == s3.a && s3.b == s1.a;
    }

    
    std::string point_location(Point const& p) {
    
        // I want to use the definition of LEFT / ON_SEGMENT for each side
        std::string ss1 = s1.point_location(p);
        std::string ss2 = s2.point_location(p);
        std::string ss3 = s3.point_location(p);
        
        // handles degenerate cases of border points
        int numborder = 0;
        if (ss1 == "ON_SEGMENT")
            ++numborder;
        if (ss2 == "ON_SEGMENT")
            ++numborder;
        if (ss3 == "ON_SEGMENT")
            ++numborder;
                        
        if (numborder > 0)
            return "BORDER";

        // handle plain cases: if the ray intersects the triangle an odd number
        //              of times, the point is inside
        int num = 0;
        Ray ray;
        ray.p = p;
        if (ray.intersects(s1)) {
            // debug
            //std::cout << ray << " intersects " << s1 << std::endl;
            ++num;
        }
        if (ray.intersects(s2)) {
            // debug
            //std::cout << ray << " intersects " << s2 << std::endl;
            ++num;
        }
        if (ray.intersects(s3)) {
            // debug
            //std::cout << ray << " intersects " << s3 << std::endl;
            ++num;
        }        
        if (num == 1) return "INSIDE";
        return "OUTSIDE";
    }

};

std::ostream& operator<<(std::ostream& os, Triangle t) {
    os << "triangle (" << t.s1 << "," << t.s2 << "," << t.s3 << ")";
    return os;
}



} 


int main() {

    using namespace std;

    // load the segment from input
    string strTriangle;
    getline(cin, strTriangle);
    istringstream streamTriangle(strTriangle);

    crs::Point p1, p2, p3;
    string strTrianglePoint;

    getline(streamTriangle, strTrianglePoint, ' ');
    istringstream(strTrianglePoint) >> p1.x;
    getline(streamTriangle, strTrianglePoint, ' ');
    istringstream(strTrianglePoint) >> p1.y;

    getline(streamTriangle, strTrianglePoint, ' ');
    istringstream(strTrianglePoint) >> p2.x;
    getline(streamTriangle, strTrianglePoint, ' ');
    istringstream(strTrianglePoint) >> p2.y;

    getline(streamTriangle, strTrianglePoint, ' ');
    istringstream(strTrianglePoint) >> p3.x;
    getline(streamTriangle, strTrianglePoint, ' ');
    istringstream(strTrianglePoint) >> p3.y;

    crs::Triangle triangle;
    triangle.s1.a = p1;
    triangle.s1.b = p2;
    triangle.s2.a = p2;
    triangle.s2.b = p3;
    triangle.s3.a = p3;
    triangle.s3.b = p1;


    if (triangle.is_valid()) {   

        // debug
        //cout << "got " << triangle << endl;

        // find how many input there are 
        int numdata;
        string inputnum;
        getline(cin, inputnum);
        istringstream(inputnum) >> numdata;

        // get points and evaluate each one
        while (numdata--) {
            
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
     
            string location = triangle.point_location(point);

            // debug
            // cout << point << " is " << location << endl;

            // result
            cout << location << endl;

        }
        
    } else {
    
        cerr << triangle << " is invalid" << endl;
    }


    return 0;
}
