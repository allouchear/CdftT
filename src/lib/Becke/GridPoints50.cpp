#include<iostream>
#include <Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints50()
{
    _Npts = 50;
    _Lmax = 11;
    _L2max = 5;
    _LebedevGridPoints = {
        // (   theta,            phi,            weight) 
        {  1.570796326795,  0.000000000000,  0.012698412698},
        {  1.570796326795,  3.141592653590,  0.012698412698},
        {  1.570796326795,  1.570796326795,  0.012698412698},
        {  1.570796326795,  4.712388980385,  0.012698412698},
        {  0.000000000000,  3.141592653590,  0.012698412698},
        {  3.141592653590,  3.141592653590,  0.012698412698},
        {  0.785398163397,  1.570796326795,  0.022574955908},
        {  0.785398163397,  4.712388980385,  0.022574955908},
        {  2.356194490192,  1.570796326795,  0.022574955908},
        {  2.356194490192,  4.712388980385,  0.022574955908},
        {  0.785398163397,  0.000000000000,  0.022574955908},
        {  0.785398163397,  3.141592653590,  0.022574955908},
        {  2.356194490192,  0.000000000000,  0.022574955908},
        {  2.356194490192,  3.141592653590,  0.022574955908},
        {  1.570796326795,  0.785398163397,  0.022574955908},
        {  1.570796326795,  2.356194490192,  0.022574955908},
        {  1.570796326795,  5.497787143782,  0.022574955908},
        {  1.570796326795,  3.926990816987,  0.022574955908},
        {  0.955316618125,  0.785398163397,  0.021093750000},
        {  0.955316618125,  2.356194490192,  0.021093750000},
        {  0.955316618125,  5.497787143782,  0.021093750000},
        {  0.955316618125,  3.926990816987,  0.021093750000},
        {  2.186276035465,  0.785398163397,  0.021093750000},
        {  2.186276035465,  2.356194490192,  0.021093750000},
        {  2.186276035465,  5.497787143782,  0.021093750000},
        {  2.186276035465,  3.926990816987,  0.021093750000},
        {  0.440510663005,  0.785398163397,  0.020173335538},
        {  0.440510663005,  2.356194490192,  0.020173335538},
        {  0.440510663005,  5.497787143782,  0.020173335538},
        {  0.440510663005,  3.926990816987,  0.020173335538},
        {  2.701081990585,  0.785398163397,  0.020173335538},
        {  2.701081990585,  2.356194490192,  0.020173335538},
        {  2.701081990585,  5.497787143782,  0.020173335538},
        {  2.701081990585,  3.926990816987,  0.020173335538},
        {  1.264518957625,  1.249045772398,  0.020173335538},
        {  1.264518957625,  1.892546881192,  0.020173335538},
        {  1.264518957625,  5.034139534781,  0.020173335538},
        {  1.264518957625,  4.390638425988,  0.020173335538},
        {  1.877073695965,  1.249045772398,  0.020173335538},
        {  1.877073695965,  1.892546881192,  0.020173335538},
        {  1.877073695965,  5.034139534781,  0.020173335538},
        {  1.877073695965,  4.390638425988,  0.020173335538},
        {  1.264518957625,  0.321750554397,  0.020173335538},
        {  1.264518957625,  2.819842099193,  0.020173335538},
        {  1.264518957625,  5.961434752783,  0.020173335538},
        {  1.264518957625,  3.463343207986,  0.020173335538},
        {  1.877073695965,  0.321750554397,  0.020173335538},
        {  1.877073695965,  2.819842099193,  0.020173335538},
        {  1.877073695965,  5.961434752783,  0.020173335538},
        {  1.877073695965,  3.463343207986,  0.020173335538}
    };
}