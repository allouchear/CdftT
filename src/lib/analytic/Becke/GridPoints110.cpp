#include<iostream>
#include<analytic/Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints110()
{
    _Npts = 110;
    _Lmax = 17;
    _L2max = 8;
    _LebedevGridPoints = {
        // (   theta,            phi,            weight) 
        {  1.570796326795,  0.000000000000,  0.003828270495},
        {  1.570796326795,  3.141592653590,  0.003828270495},
        {  1.570796326795,  1.570796326795,  0.003828270495},
        {  1.570796326795,  4.712388980385,  0.003828270495},
        {  0.000000000000,  3.141592653590,  0.003828270495},
        {  3.141592653590,  3.141592653590,  0.003828270495},
        {  0.955316618125,  0.785398163397,  0.009793737512},
        {  0.955316618125,  2.356194490192,  0.009793737512},
        {  0.955316618125,  5.497787143782,  0.009793737512},
        {  0.955316618125,  3.926990816987,  0.009793737512},
        {  2.186276035465,  0.785398163397,  0.009793737512},
        {  2.186276035465,  2.356194490192,  0.009793737512},
        {  2.186276035465,  5.497787143782,  0.009793737512},
        {  2.186276035465,  3.926990816987,  0.009793737512},
        {  0.264879572022,  0.785398163397,  0.008211737283},
        {  0.264879572022,  2.356194490192,  0.008211737283},
        {  0.264879572022,  5.497787143782,  0.008211737283},
        {  0.264879572022,  3.926990816987,  0.008211737283},
        {  2.876713081568,  0.785398163397,  0.008211737283},
        {  2.876713081568,  2.356194490192,  0.008211737283},
        {  2.876713081568,  5.497787143782,  0.008211737283},
        {  2.876713081568,  3.926990816987,  0.008211737283},
        {  1.384606796719,  1.381292828874,  0.008211737283},
        {  1.384606796719,  1.760299824716,  0.008211737283},
        {  1.384606796719,  4.901892478306,  0.008211737283},
        {  1.384606796719,  4.522885482463,  0.008211737283},
        {  1.756985856871,  1.381292828874,  0.008211737283},
        {  1.756985856871,  1.760299824716,  0.008211737283},
        {  1.756985856871,  4.901892478306,  0.008211737283},
        {  1.756985856871,  4.522885482463,  0.008211737283},
        {  1.384606796719,  0.189503497921,  0.008211737283},
        {  1.384606796719,  2.952089155668,  0.008211737283},
        {  1.384606796719,  6.093681809258,  0.008211737283},
        {  1.384606796719,  3.331096151511,  0.008211737283},
        {  1.756985856871,  0.189503497921,  0.008211737283},
        {  1.756985856871,  2.952089155668,  0.008211737283},
        {  1.756985856871,  6.093681809258,  0.008211737283},
        {  1.756985856871,  3.331096151511,  0.008211737283},
        {  1.353124175903,  0.785398163397,  0.009942814891},
        {  1.353124175903,  2.356194490192,  0.009942814891},
        {  1.353124175903,  5.497787143782,  0.009942814891},
        {  1.353124175903,  3.926990816987,  0.009942814891},
        {  1.788468477687,  0.785398163397,  0.009942814891},
        {  1.788468477687,  2.356194490192,  0.009942814891},
        {  1.788468477687,  5.497787143782,  0.009942814891},
        {  1.788468477687,  3.926990816987,  0.009942814891},
        {  0.808725400927,  0.303149695109,  0.009942814891},
        {  0.808725400927,  2.838442958480,  0.009942814891},
        {  0.808725400927,  5.980035612070,  0.009942814891},
        {  0.808725400927,  3.444742348699,  0.009942814891},
        {  2.332867252663,  0.303149695109,  0.009942814891},
        {  2.332867252663,  2.838442958480,  0.009942814891},
        {  2.332867252663,  5.980035612070,  0.009942814891},
        {  2.332867252663,  3.444742348699,  0.009942814891},
        {  0.808725400927,  1.267646631686,  0.009942814891},
        {  0.808725400927,  1.873946021904,  0.009942814891},
        {  0.808725400927,  5.015538675494,  0.009942814891},
        {  0.808725400927,  4.409239285275,  0.009942814891},
        {  2.332867252663,  1.267646631686,  0.009942814891},
        {  2.332867252663,  1.873946021904,  0.009942814891},
        {  2.332867252663,  5.015538675494,  0.009942814891},
        {  2.332867252663,  4.409239285275,  0.009942814891},
        {  0.593890307361,  0.785398163397,  0.009595471336},
        {  0.593890307361,  2.356194490192,  0.009595471336},
        {  0.593890307361,  5.497787143782,  0.009595471336},
        {  0.593890307361,  3.926990816987,  0.009595471336},
        {  2.547702346229,  0.785398163397,  0.009595471336},
        {  2.547702346229,  2.356194490192,  0.009595471336},
        {  2.547702346229,  5.497787143782,  0.009595471336},
        {  2.547702346229,  3.926990816987,  0.009595471336},
        {  1.163977851408,  1.125357546251,  0.009595471336},
        {  1.163977851408,  2.016235107339,  0.009595471336},
        {  1.163977851408,  5.157827760929,  0.009595471336},
        {  1.163977851408,  4.266950199841,  0.009595471336},
        {  1.977614802182,  1.125357546251,  0.009595471336},
        {  1.977614802182,  2.016235107339,  0.009595471336},
        {  1.977614802182,  5.157827760929,  0.009595471336},
        {  1.977614802182,  4.266950199841,  0.009595471336},
        {  1.163977851408,  0.445438780544,  0.009595471336},
        {  1.163977851408,  2.696153873046,  0.009595471336},
        {  1.163977851408,  5.837746526635,  0.009595471336},
        {  1.163977851408,  3.587031434134,  0.009595471336},
        {  1.977614802182,  0.445438780544,  0.009595471336},
        {  1.977614802182,  2.696153873046,  0.009595471336},
        {  1.977614802182,  5.837746526635,  0.009595471336},
        {  1.977614802182,  3.587031434134,  0.009595471336},
        {  1.570796326795,  1.071999817948,  0.009694996362},
        {  1.570796326795,  2.069592835642,  0.009694996362},
        {  1.570796326795,  5.211185489232,  0.009694996362},
        {  1.570796326795,  4.213592471538,  0.009694996362},
        {  1.570796326795,  0.498796508847,  0.009694996362},
        {  1.570796326795,  2.642796144743,  0.009694996362},
        {  1.570796326795,  5.784388798333,  0.009694996362},
        {  1.570796326795,  3.640389162437,  0.009694996362},
        {  0.498796508847,  0.000000000000,  0.009694996362},
        {  0.498796508847,  3.141592653590,  0.009694996362},
        {  2.642796144743,  0.000000000000,  0.009694996362},
        {  2.642796144743,  3.141592653590,  0.009694996362},
        {  1.071999817948,  0.000000000000,  0.009694996362},
        {  1.071999817948,  3.141592653590,  0.009694996362},
        {  2.069592835642,  0.000000000000,  0.009694996362},
        {  2.069592835642,  3.141592653590,  0.009694996362},
        {  0.498796508847,  1.570796326795,  0.009694996362},
        {  0.498796508847,  4.712388980385,  0.009694996362},
        {  2.642796144743,  1.570796326795,  0.009694996362},
        {  2.642796144743,  4.712388980385,  0.009694996362},
        {  1.071999817948,  1.570796326795,  0.009694996362},
        {  1.071999817948,  4.712388980385,  0.009694996362},
        {  2.069592835642,  1.570796326795,  0.009694996362},
        {  2.069592835642,  4.712388980385,  0.009694996362}    
    };
}