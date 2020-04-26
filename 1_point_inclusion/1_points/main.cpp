/*
c++ std: 11
compiler: g++ v 7.4.0

compile with 
g++ -I. main.cpp

launch with
cat data.txt | ./a.out 
*/

#include <cmath> // abs
#include <cstdio> // exit
#include <iostream>
#include <sstream>
#include <string>


namespace crs {

struct Point {
    double x,y;
};

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
    
    double length() {
        return std::sqrt(std::pow(b.x-a.x, 2) + std::pow(b.y-a.y, 2));
    }

private:
    double det(Point a, Point b, Point c) {
        return (b.x - a.x)*(c.y- a.y) - (c.x - a.x)*(b.y - a.y);
    }
};

} 


int main() {

    using namespace std;

    // load the segment from input
    string inputseg;
    getline(cin, inputseg);
    istringstream segStream(inputseg);

    // put segment data into Segment class from input stream
    crs::Segment segment;
    string inputdata;
    getline(segStream, inputdata, ' ');
    istringstream(inputdata) >> segment.a.x;
    getline(segStream, inputdata, ' ');
    istringstream(inputdata) >> segment.a.y;
    getline(segStream, inputdata, ' ');
    istringstream(inputdata) >> segment.b.x;
    getline(segStream, inputdata, ' ');
    istringstream(inputdata) >> segment.b.y;
    
    if (!segment.length()) {
        cerr << "invalid segment: dengenerate, (size = 0)" << endl; 
        exit(1);
    }


    // debug
    //cout << "got segment: (" << segment.a.x << "," << segment.a.y 
    //                        << ") --> (" 
    //                        << segment.b.x << "," << segment.b.y << ")" << endl;

    // find how many input there are 
    int numdata;
    string inputnum;
    getline(cin, inputnum);
    istringstream(inputnum) >> numdata;


    if (!(numdata >= 1 && numdata <= 10000)) {
        cerr << "invalid input n(points) > 1000" << endl; 
        exit(1);
    }

    double input_limit = 10000;

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
 
        if (abs(point.x) > input_limit || abs(point.y) > input_limit) {
            cerr << "bad input, invalid entry, " 
                 << point << " exceeds limit " << input_limit << endl;
            continue;
        }

        string location = segment.point_location(point);

        // debug
        //cout << "point (" << point.x << "," << point.y << ") is " << location << endl;

        // result
        cout << location << endl;

    }
    

    return 0;
}
