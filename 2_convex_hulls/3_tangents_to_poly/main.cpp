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
#include <sstream>


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
