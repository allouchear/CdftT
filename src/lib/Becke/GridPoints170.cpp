#include<iostream>
#include <Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints170()
{
    _Npts = 170;
    _Lmax = 21;
    _L2max = 10;
    _LebedevGridPoints = {
        // (   theta,            phi,            weight) 
        {  1.570796326795,  0.000000000000,  0.005544842902},
        {  1.570796326795,  3.141592653590,  0.005544842902},
        {  1.570796326795,  1.570796326795,  0.005544842902},
        {  1.570796326795,  4.712388980385,  0.005544842902},
        {  0.000000000000,  3.141592653590,  0.005544842902},
        {  3.141592653590,  3.141592653590,  0.005544842902},
        {  0.785398163397,  1.570796326795,  0.006071332771},
        {  0.785398163397,  4.712388980385,  0.006071332771},
        {  2.356194490192,  1.570796326795,  0.006071332771},
        {  2.356194490192,  4.712388980385,  0.006071332771},
        {  0.785398163397,  0.000000000000,  0.006071332771},
        {  0.785398163397,  3.141592653590,  0.006071332771},
        {  2.356194490192,  0.000000000000,  0.006071332771},
        {  2.356194490192,  3.141592653590,  0.006071332771},
        {  1.570796326795,  0.785398163397,  0.006071332771},
        {  1.570796326795,  2.356194490192,  0.006071332771},
        {  1.570796326795,  5.497787143782,  0.006071332771},
        {  1.570796326795,  3.926990816987,  0.006071332771},
        {  0.955316618125,  0.785398163397,  0.006383674774},
        {  0.955316618125,  2.356194490192,  0.006383674774},
        {  0.955316618125,  5.497787143782,  0.006383674774},
        {  0.955316618125,  3.926990816987,  0.006383674774},
        {  2.186276035465,  0.785398163397,  0.006383674774},
        {  2.186276035465,  2.356194490192,  0.006383674774},
        {  2.186276035465,  5.497787143782,  0.006383674774},
        {  2.186276035465,  3.926990816987,  0.006383674774},
        {  0.369127250133,  0.785398163397,  0.005183387588},
        {  0.369127250133,  2.356194490192,  0.005183387588},
        {  0.369127250133,  5.497787143782,  0.005183387588},
        {  0.369127250133,  3.926990816987,  0.005183387588},
        {  2.772465403456,  0.785398163397,  0.005183387588},
        {  2.772465403456,  2.356194490192,  0.005183387588},
        {  2.772465403456,  5.497787143782,  0.005183387588},
        {  2.772465403456,  3.926990816987,  0.005183387588},
        {  1.312819076651,  1.303777789065,  0.005183387588},
        {  1.312819076651,  1.837814864525,  0.005183387588},
        {  1.312819076651,  4.979407518115,  0.005183387588},
        {  1.312819076651,  4.445370442655,  0.005183387588},
        {  1.828773576939,  1.303777789065,  0.005183387588},
        {  1.828773576939,  1.837814864525,  0.005183387588},
        {  1.828773576939,  4.979407518115,  0.005183387588},
        {  1.828773576939,  4.445370442655,  0.005183387588},
        {  1.312819076651,  0.267018537730,  0.005183387588},
        {  1.312819076651,  2.874574115860,  0.005183387588},
        {  1.312819076651,  6.016166769449,  0.005183387588},
        {  1.312819076651,  3.408611191320,  0.005183387588},
        {  1.828773576939,  0.267018537730,  0.005183387588},
        {  1.828773576939,  2.874574115860,  0.005183387588},
        {  1.828773576939,  6.016166769449,  0.005183387588},
        {  1.828773576939,  3.408611191320,  0.005183387588},
        {  1.265271650081,  0.785398163397,  0.006317929010},
        {  1.265271650081,  2.356194490192,  0.006317929010},
        {  1.265271650081,  5.497787143782,  0.006317929010},
        {  1.265271650081,  3.926990816987,  0.006317929010},
        {  1.876321003509,  0.785398163397,  0.006317929010},
        {  1.876321003509,  2.356194490192,  0.006317929010},
        {  1.876321003509,  5.497787143782,  0.006317929010},
        {  1.876321003509,  3.926990816987,  0.006317929010},
        {  0.830698505928,  0.419558379604,  0.006317929010},
        {  0.830698505928,  2.722034273986,  0.006317929010},
        {  0.830698505928,  5.863626927575,  0.006317929010},
        {  0.830698505928,  3.561151033194,  0.006317929010},
        {  2.310894147661,  0.419558379604,  0.006317929010},
        {  2.310894147661,  2.722034273986,  0.006317929010},
        {  2.310894147661,  5.863626927575,  0.006317929010},
        {  2.310894147661,  3.561151033194,  0.006317929010},
        {  0.830698505928,  1.151237947191,  0.006317929010},
        {  0.830698505928,  1.990354706399,  0.006317929010},
        {  0.830698505928,  5.131947359989,  0.006317929010},
        {  0.830698505928,  4.292830600781,  0.006317929010},
        {  2.310894147661,  1.151237947191,  0.006317929010},
        {  2.310894147661,  1.990354706399,  0.006317929010},
        {  2.310894147661,  5.131947359989,  0.006317929010},
        {  2.310894147661,  4.292830600781,  0.006317929010},
        {  0.657053154535,  0.785398163397,  0.006201670007},
        {  0.657053154535,  2.356194490192,  0.006201670007},
        {  0.657053154535,  5.497787143782,  0.006201670007},
        {  0.657053154535,  3.926990816987,  0.006201670007},
        {  2.484539499055,  0.785398163397,  0.006201670007},
        {  2.484539499055,  2.356194490192,  0.006201670007},
        {  2.484539499055,  5.497787143782,  0.006201670007},
        {  2.484539499055,  3.926990816987,  0.006201670007},
        {  1.124207897725,  1.071447091197,  0.006201670007},
        {  1.124207897725,  2.070145562393,  0.006201670007},
        {  1.124207897725,  5.211738215983,  0.006201670007},
        {  1.124207897725,  4.213039744787,  0.006201670007},
        {  2.017384755865,  1.071447091197,  0.006201670007},
        {  2.017384755865,  2.070145562393,  0.006201670007},
        {  2.017384755865,  5.211738215983,  0.006201670007},
        {  2.017384755865,  4.213039744787,  0.006201670007},
        {  1.124207897725,  0.499349235598,  0.006201670007},
        {  1.124207897725,  2.642243417992,  0.006201670007},
        {  1.124207897725,  5.783836071582,  0.006201670007},
        {  1.124207897725,  3.640941889188,  0.006201670007},
        {  2.017384755865,  0.499349235598,  0.006201670007},
        {  2.017384755865,  2.642243417992,  0.006201670007},
        {  2.017384755865,  5.783836071582,  0.006201670007},
        {  2.017384755865,  3.640941889188,  0.006201670007},
        {  1.570796326795,  1.306331088687,  0.005477143385},
        {  1.570796326795,  1.835261564903,  0.005477143385},
        {  1.570796326795,  4.976854218493,  0.005477143385},
        {  1.570796326795,  4.447923742276,  0.005477143385},
        {  1.570796326795,  0.264465238108,  0.005477143385},
        {  1.570796326795,  2.877127415481,  0.005477143385},
        {  1.570796326795,  6.018720069071,  0.005477143385},
        {  1.570796326795,  3.406057891698,  0.005477143385},
        {  0.264465238108,  0.000000000000,  0.005477143385},
        {  0.264465238108,  3.141592653590,  0.005477143385},
        {  2.877127415481,  0.000000000000,  0.005477143385},
        {  2.877127415481,  3.141592653590,  0.005477143385},
        {  1.306331088687,  0.000000000000,  0.005477143385},
        {  1.306331088687,  3.141592653590,  0.005477143385},
        {  1.835261564903,  0.000000000000,  0.005477143385},
        {  1.835261564903,  3.141592653590,  0.005477143385},
        {  0.264465238108,  1.570796326795,  0.005477143385},
        {  0.264465238108,  4.712388980385,  0.005477143385},
        {  2.877127415481,  1.570796326795,  0.005477143385},
        {  2.877127415481,  4.712388980385,  0.005477143385},
        {  1.306331088687,  1.570796326795,  0.005477143385},
        {  1.306331088687,  4.712388980385,  0.005477143385},
        {  1.835261564903,  1.570796326795,  0.005477143385},
        {  1.835261564903,  4.712388980385,  0.005477143385},
        {  0.546370868070,  0.282146391436,  0.005968383988},
        {  0.546370868070,  2.859446262153,  0.005968383988},
        {  0.546370868070,  6.001038915743,  0.005968383988},
        {  0.546370868070,  3.423739045026,  0.005968383988},
        {  2.595221785520,  0.282146391436,  0.005968383988},
        {  2.595221785520,  2.859446262153,  0.005968383988},
        {  2.595221785520,  6.001038915743,  0.005968383988},
        {  2.595221785520,  3.423739045026,  0.005968383988},
        {  1.425623870147,  1.042166590272,  0.005968383988},
        {  1.425623870147,  2.099426063318,  0.005968383988},
        {  1.425623870147,  5.241018716908,  0.005968383988},
        {  1.425623870147,  4.183759243862,  0.005968383988},
        {  1.715968783443,  1.042166590272,  0.005968383988},
        {  1.715968783443,  2.099426063318,  0.005968383988},
        {  1.715968783443,  5.241018716908,  0.005968383988},
        {  1.715968783443,  4.183759243862,  0.005968383988},
        {  0.546370868070,  1.288649935359,  0.005968383988},
        {  0.546370868070,  1.852942718231,  0.005968383988},
        {  0.546370868070,  4.994535371821,  0.005968383988},
        {  0.546370868070,  4.430242588948,  0.005968383988},
        {  2.595221785520,  1.288649935359,  0.005968383988},
        {  2.595221785520,  1.852942718231,  0.005968383988},
        {  2.595221785520,  4.994535371821,  0.005968383988},
        {  2.595221785520,  4.430242588948,  0.005968383988},
        {  1.048299574758,  1.403074664045,  0.005968383988},
        {  1.048299574758,  1.738517989544,  0.005968383988},
        {  1.048299574758,  4.880110643134,  0.005968383988},
        {  1.048299574758,  4.544667317635,  0.005968383988},
        {  2.093293078832,  1.403074664045,  0.005968383988},
        {  2.093293078832,  1.738517989544,  0.005968383988},
        {  2.093293078832,  4.880110643134,  0.005968383988},
        {  2.093293078832,  4.544667317635,  0.005968383988},
        {  1.425623870147,  0.528629736523,  0.005968383988},
        {  1.425623870147,  2.612962917067,  0.005968383988},
        {  1.425623870147,  5.754555570657,  0.005968383988},
        {  1.425623870147,  3.670222390113,  0.005968383988},
        {  1.715968783443,  0.528629736523,  0.005968383988},
        {  1.715968783443,  2.612962917067,  0.005968383988},
        {  1.715968783443,  5.754555570657,  0.005968383988},
        {  1.715968783443,  3.670222390113,  0.005968383988},
        {  1.048299574758,  0.167721662749,  0.005968383988},
        {  1.048299574758,  2.973870990840,  0.005968383988},
        {  1.048299574758,  6.115463644430,  0.005968383988},
        {  1.048299574758,  3.309314316339,  0.005968383988},
        {  2.093293078832,  0.167721662749,  0.005968383988},
        {  2.093293078832,  2.973870990840,  0.005968383988},
        {  2.093293078832,  6.115463644430,  0.005968383988},
        {  2.093293078832,  3.309314316339,  0.005968383988}    
    };
}