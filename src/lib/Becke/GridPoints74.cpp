#include<iostream>
#include <Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints74()
{
    _Npts = 74;
    _Lmax = 13;
    _L2max = 6;
    _LebedevGridPoints = {
        // (   theta,            phi,            weight) 
        {  1.570796326795,  0.000000000000,  0.000513067180},
        {  1.570796326795,  3.141592653590,  0.000513067180},
        {  1.570796326795,  1.570796326795,  0.000513067180},
        {  1.570796326795,  4.712388980385,  0.000513067180},
        {  0.000000000000,  3.141592653590,  0.000513067180},
        {  3.141592653590,  3.141592653590,  0.000513067180},
        {  0.785398163397,  1.570796326795,  0.016604069566},
        {  0.785398163397,  4.712388980385,  0.016604069566},
        {  2.356194490192,  1.570796326795,  0.016604069566},
        {  2.356194490192,  4.712388980385,  0.016604069566},
        {  0.785398163397,  0.000000000000,  0.016604069566},
        {  0.785398163397,  3.141592653590,  0.016604069566},
        {  2.356194490192,  0.000000000000,  0.016604069566},
        {  2.356194490192,  3.141592653590,  0.016604069566},
        {  1.570796326795,  0.785398163397,  0.016604069566},
        {  1.570796326795,  2.356194490192,  0.016604069566},
        {  1.570796326795,  5.497787143782,  0.016604069566},
        {  1.570796326795,  3.926990816987,  0.016604069566},
        {  0.955316618125,  0.785398163397, -0.029586038961},
        {  0.955316618125,  2.356194490192, -0.029586038961},
        {  0.955316618125,  5.497787143782, -0.029586038961},
        {  0.955316618125,  3.926990816987, -0.029586038961},
        {  2.186276035465,  0.785398163397, -0.029586038961},
        {  2.186276035465,  2.356194490192, -0.029586038961},
        {  2.186276035465,  5.497787143782, -0.029586038961},
        {  2.186276035465,  3.926990816987, -0.029586038961},
        {  0.746898593069,  0.785398163397,  0.026576207082},
        {  0.746898593069,  2.356194490192,  0.026576207082},
        {  0.746898593069,  5.497787143782,  0.026576207082},
        {  0.746898593069,  3.926990816987,  0.026576207082},
        {  2.394694060521,  0.785398163397,  0.026576207082},
        {  2.394694060521,  2.356194490192,  0.026576207082},
        {  2.394694060521,  5.497787143782,  0.026576207082},
        {  2.394694060521,  3.926990816987,  0.026576207082},
        {  1.069703313530,  0.991156586431,  0.026576207082},
        {  1.069703313530,  2.150436067159,  0.026576207082},
        {  1.069703313530,  5.292028720748,  0.026576207082},
        {  1.069703313530,  4.132749240021,  0.026576207082},
        {  2.071889340060,  0.991156586431,  0.026576207082},
        {  2.071889340060,  2.150436067159,  0.026576207082},
        {  2.071889340060,  5.292028720748,  0.026576207082},
        {  2.071889340060,  4.132749240021,  0.026576207082},
        {  1.069703313530,  0.579639740364,  0.026576207082},
        {  1.069703313530,  2.561952913226,  0.026576207082},
        {  1.069703313530,  5.703545566816,  0.026576207082},
        {  1.069703313530,  3.721232393953,  0.026576207082},
        {  2.071889340060,  0.579639740364,  0.026576207082},
        {  2.071889340060,  2.561952913226,  0.026576207082},
        {  2.071889340060,  5.703545566816,  0.026576207082},
        {  2.071889340060,  3.721232393953,  0.026576207082},
        {  1.570796326795,  1.244251195420,  0.016522170994},
        {  1.570796326795,  1.897341458170,  0.016522170994},
        {  1.570796326795,  5.038934111760,  0.016522170994},
        {  1.570796326795,  4.385843849009,  0.016522170994},
        {  1.570796326795,  0.326545131375,  0.016522170994},
        {  1.570796326795,  2.815047522215,  0.016522170994},
        {  1.570796326795,  5.956640175804,  0.016522170994},
        {  1.570796326795,  3.468137784965,  0.016522170994},
        {  0.326545131375,  0.000000000000,  0.016522170994},
        {  0.326545131375,  3.141592653590,  0.016522170994},
        {  2.815047522215,  0.000000000000,  0.016522170994},
        {  2.815047522215,  3.141592653590,  0.016522170994},
        {  1.244251195420,  0.000000000000,  0.016522170994},
        {  1.244251195420,  3.141592653590,  0.016522170994},
        {  1.897341458170,  0.000000000000,  0.016522170994},
        {  1.897341458170,  3.141592653590,  0.016522170994},
        {  0.326545131375,  1.570796326795,  0.016522170994},
        {  0.326545131375,  4.712388980385,  0.016522170994},
        {  2.815047522215,  1.570796326795,  0.016522170994},
        {  2.815047522215,  4.712388980385,  0.016522170994},
        {  1.244251195420,  1.570796326795,  0.016522170994},
        {  1.244251195420,  4.712388980385,  0.016522170994},
        {  1.897341458170,  1.570796326795,  0.016522170994},
        {  1.897341458170,  4.712388980385,  0.016522170994}
    };
}