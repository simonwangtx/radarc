#include <iostream>
#include "coord.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    RadarConv radarConv;
    GEO geo = GEO{31.20740278,121.34019167,34.85};
    Radar radar = radarConv.Geo2Radar(geo);
    cout<<"radar.x: "<<radar.x<<endl;
    cout<<"radar.y: "<<radar.y<<endl;
    cout<<"radar.z: "<<radar.z<<endl;
    cout<<"radar.theta: "<<radar.theta<<endl;
    cout<<"radar.range: "<<radar.range<<endl;
    cout<<"radar.rangeXY: "<<radar.rangeXY<<endl;
    cout<<"radar.level: "<<radar.level<<endl;
    return 0;
}