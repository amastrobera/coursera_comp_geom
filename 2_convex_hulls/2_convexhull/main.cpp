/*
c++ std: 11
compiler: g++ v 7.4.0

compile with 
g++ -I.  main.cpp

launch with
cat data.txt | ./a.out 
*/

#include "point.hpp"
#include "segment.hpp"
#include "polygon.hpp"
#include "algos.hpp"

#include <cstdio> // exit
#include <iostream>


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
    /*
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
    */
        
    if (input_error) {    
        cerr << "bad input, invalid point exceeds limit " << input_limit << endl;
        exit(1);
    }
    
    crs::Polygon cvpoly = crs::convex_hull(points);
    
    if (!cvpoly.is_convex()) {    
        cerr << "bad output, not convex\n " << cvpoly << endl;
        exit(1);
    }
    
    // debug
    //cout << "output" << endl;

    size_t n = cvpoly.vertices.size();    
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
