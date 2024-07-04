#include<iostream>
#include <Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints230()
{
    _Npts = 230;
    _Lmax = 25;
    _L2max = 12;
    _LebedevGridPoints = {
        // (   theta,            phi,            weight) 
        {  1.570796326795,  0.000000000000, -0.055226399197},
        {  1.570796326795,  3.141592653590, -0.055226399197},
        {  1.570796326795,  1.570796326795, -0.055226399197},
        {  1.570796326795,  4.712388980385, -0.055226399197},
        {  0.000000000000,  3.141592653590, -0.055226399197},
        {  3.141592653590,  3.141592653590, -0.055226399197},
        {  0.955316618125,  0.785398163397,  0.004450274607},
        {  0.955316618125,  2.356194490192,  0.004450274607},
        {  0.955316618125,  5.497787143782,  0.004450274607},
        {  0.955316618125,  3.926990816987,  0.004450274607},
        {  2.186276035465,  0.785398163397,  0.004450274607},
        {  2.186276035465,  2.356194490192,  0.004450274607},
        {  2.186276035465,  5.497787143782,  0.004450274607},
        {  2.186276035465,  3.926990816987,  0.004450274607},
        {  0.688359441480,  0.785398163397,  0.004496841068},
        {  0.688359441480,  2.356194490192,  0.004496841068},
        {  0.688359441480,  5.497787143782,  0.004496841068},
        {  0.688359441480,  3.926990816987,  0.004496841068},
        {  2.453233212110,  0.785398163397,  0.004496841068},
        {  2.453233212110,  2.356194490192,  0.004496841068},
        {  2.453233212110,  5.497787143782,  0.004496841068},
        {  2.453233212110,  3.926990816987,  0.004496841068},
        {  1.104921612004,  1.043976411203,  0.004496841068},
        {  1.104921612004,  2.097616242387,  0.004496841068},
        {  1.104921612004,  5.239208895977,  0.004496841068},
        {  1.104921612004,  4.185569064792,  0.004496841068},
        {  2.036671041586,  1.043976411203,  0.004496841068},
        {  2.036671041586,  2.097616242387,  0.004496841068},
        {  2.036671041586,  5.239208895977,  0.004496841068},
        {  2.036671041586,  4.185569064792,  0.004496841068},
        {  1.104921612004,  0.526819915592,  0.004496841068},
        {  1.104921612004,  2.614772737998,  0.004496841068},
        {  1.104921612004,  5.756365391587,  0.004496841068},
        {  1.104921612004,  3.668412569182,  0.004496841068},
        {  2.036671041586,  0.526819915592,  0.004496841068},
        {  2.036671041586,  2.614772737998,  0.004496841068},
        {  2.036671041586,  5.756365391587,  0.004496841068},
        {  2.036671041586,  3.668412569182,  0.004496841068},
        {  0.364456068749,  0.785398163397,  0.005049153450},
        {  0.364456068749,  2.356194490192,  0.005049153450},
        {  0.364456068749,  5.497787143782,  0.005049153450},
        {  0.364456068749,  3.926990816987,  0.005049153450},
        {  2.777136584841,  0.785398163397,  0.005049153450},
        {  2.777136584841,  2.356194490192,  0.005049153450},
        {  2.777136584841,  5.497787143782,  0.005049153450},
        {  2.777136584841,  3.926990816987,  0.005049153450},
        {  1.316006579721,  1.307307814523,  0.005049153450},
        {  1.316006579721,  1.834284839067,  0.005049153450},
        {  1.316006579721,  4.975877492657,  0.005049153450},
        {  1.316006579721,  4.448900468112,  0.005049153450},
        {  1.825586073869,  1.307307814523,  0.005049153450},
        {  1.825586073869,  1.834284839067,  0.005049153450},
        {  1.825586073869,  4.975877492657,  0.005049153450},
        {  1.825586073869,  4.448900468112,  0.005049153450},
        {  1.316006579721,  0.263488512272,  0.005049153450},
        {  1.316006579721,  2.878104141317,  0.005049153450},
        {  1.316006579721,  6.019696794907,  0.005049153450},
        {  1.316006579721,  3.405081165862,  0.005049153450},
        {  1.825586073869,  0.263488512272,  0.005049153450},
        {  1.825586073869,  2.878104141317,  0.005049153450},
        {  1.825586073869,  6.019696794907,  0.005049153450},
        {  1.825586073869,  3.405081165862,  0.005049153450},
        {  1.411825387672,  0.785398163397,  0.003976408018},
        {  1.411825387672,  2.356194490192,  0.003976408018},
        {  1.411825387672,  5.497787143782,  0.003976408018},
        {  1.411825387672,  3.926990816987,  0.003976408018},
        {  1.729767265918,  0.785398163397,  0.003976408018},
        {  1.729767265918,  2.356194490192,  0.003976408018},
        {  1.729767265918,  5.497787143782,  0.003976408018},
        {  1.729767265918,  3.926990816987,  0.003976408018},
        {  0.797929269312,  0.222962425829,  0.003976408018},
        {  0.797929269312,  2.918630227761,  0.003976408018},
        {  0.797929269312,  6.060222881350,  0.003976408018},
        {  0.797929269312,  3.364555079419,  0.003976408018},
        {  2.343663384278,  0.222962425829,  0.003976408018},
        {  2.343663384278,  2.918630227761,  0.003976408018},
        {  2.343663384278,  6.060222881350,  0.003976408018},
        {  2.343663384278,  3.364555079419,  0.003976408018},
        {  0.797929269312,  1.347833900966,  0.003976408018},
        {  0.797929269312,  1.793758752624,  0.003976408018},
        {  0.797929269312,  4.935351406214,  0.003976408018},
        {  0.797929269312,  4.489426554555,  0.003976408018},
        {  2.343663384278,  1.347833900966,  0.003976408018},
        {  2.343663384278,  1.793758752624,  0.003976408018},
        {  2.343663384278,  4.935351406214,  0.003976408018},
        {  2.343663384278,  4.489426554555,  0.003976408018},
        {  1.198789539797,  0.785398163397,  0.004401400650},
        {  1.198789539797,  2.356194490192,  0.004401400650},
        {  1.198789539797,  5.497787143782,  0.004401400650},
        {  1.198789539797,  3.926990816987,  0.004401400650},
        {  1.942803113793,  0.785398163397,  0.004401400650},
        {  1.942803113793,  2.356194490192,  0.004401400650},
        {  1.942803113793,  5.497787143782,  0.004401400650},
        {  1.942803113793,  3.926990816987,  0.004401400650},
        {  0.851652805844,  0.504215580725,  0.004401400650},
        {  0.851652805844,  2.637377072865,  0.004401400650},
        {  0.851652805844,  5.778969726454,  0.004401400650},
        {  0.851652805844,  3.645808234315,  0.004401400650},
        {  2.289939847746,  0.504215580725,  0.004401400650},
        {  2.289939847746,  2.637377072865,  0.004401400650},
        {  2.289939847746,  5.778969726454,  0.004401400650},
        {  2.289939847746,  3.645808234315,  0.004401400650},
        {  0.851652805844,  1.066580746070,  0.004401400650},
        {  0.851652805844,  2.075011907520,  0.004401400650},
        {  0.851652805844,  5.216604561110,  0.004401400650},
        {  0.851652805844,  4.208173399659,  0.004401400650},
        {  2.289939847746,  1.066580746070,  0.004401400650},
        {  2.289939847746,  2.075011907520,  0.004401400650},
        {  2.289939847746,  5.216604561110,  0.004401400650},
        {  2.289939847746,  4.208173399659,  0.004401400650},
        {  0.057144733819,  0.785398163397,  0.017245443505},
        {  0.057144733819,  2.356194490192,  0.017245443505},
        {  0.057144733819,  5.497787143782,  0.017245443505},
        {  0.057144733819,  3.926990816987,  0.017245443505},
        {  3.084447919771,  0.785398163397,  0.017245443505},
        {  3.084447919771,  2.356194490192,  0.017245443505},
        {  3.084447919771,  5.497787143782,  0.017245443505},
        {  3.084447919771,  3.926990816987,  0.017245443505},
        {  1.530399900229,  1.530366898943,  0.017245443505},
        {  1.530399900229,  1.611225754647,  0.017245443505},
        {  1.530399900229,  4.752818408237,  0.017245443505},
        {  1.530399900229,  4.671959552533,  0.017245443505},
        {  1.611192753361,  1.530366898943,  0.017245443505},
        {  1.611192753361,  1.611225754647,  0.017245443505},
        {  1.611192753361,  4.752818408237,  0.017245443505},
        {  1.611192753361,  4.671959552533,  0.017245443505},
        {  1.530399900229,  0.040429427852,  0.017245443505},
        {  1.530399900229,  3.101163225738,  0.017245443505},
        {  1.530399900229,  6.242755879328,  0.017245443505},
        {  1.530399900229,  3.182022081442,  0.017245443505},
        {  1.611192753361,  0.040429427852,  0.017245443505},
        {  1.611192753361,  3.101163225738,  0.017245443505},
        {  1.611192753361,  6.242755879328,  0.017245443505},
        {  1.611192753361,  3.182022081442,  0.017245443505},
        {  1.570796326795,  0.949137761921,  0.004231083095},
        {  1.570796326795,  2.192454891668,  0.004231083095},
        {  1.570796326795,  5.334047545258,  0.004231083095},
        {  1.570796326795,  4.090730415511,  0.004231083095},
        {  1.570796326795,  0.621658564874,  0.004231083095},
        {  1.570796326795,  2.519934088716,  0.004231083095},
        {  1.570796326795,  5.661526742306,  0.004231083095},
        {  1.570796326795,  3.763251218463,  0.004231083095},
        {  0.621658564874,  0.000000000000,  0.004231083095},
        {  0.621658564874,  3.141592653590,  0.004231083095},
        {  2.519934088716,  0.000000000000,  0.004231083095},
        {  2.519934088716,  3.141592653590,  0.004231083095},
        {  0.949137761921,  0.000000000000,  0.004231083095},
        {  0.949137761921,  3.141592653590,  0.004231083095},
        {  2.192454891668,  0.000000000000,  0.004231083095},
        {  2.192454891668,  3.141592653590,  0.004231083095},
        {  0.621658564874,  1.570796326795,  0.004231083095},
        {  0.621658564874,  4.712388980385,  0.004231083095},
        {  2.519934088716,  1.570796326795,  0.004231083095},
        {  2.519934088716,  4.712388980385,  0.004231083095},
        {  0.949137761921,  1.570796326795,  0.004231083095},
        {  0.949137761921,  4.712388980385,  0.004231083095},
        {  2.192454891668,  1.570796326795,  0.004231083095},
        {  2.192454891668,  4.712388980385,  0.004231083095},
        {  1.570796326795,  1.208323206480,  0.005198069864},
        {  1.570796326795,  1.933269447110,  0.005198069864},
        {  1.570796326795,  5.074862100700,  0.005198069864},
        {  1.570796326795,  4.349915860070,  0.005198069864},
        {  1.570796326795,  0.362473120315,  0.005198069864},
        {  1.570796326795,  2.779119533275,  0.005198069864},
        {  1.570796326795,  5.920712186865,  0.005198069864},
        {  1.570796326795,  3.504065773905,  0.005198069864},
        {  0.362473120315,  0.000000000000,  0.005198069864},
        {  0.362473120315,  3.141592653590,  0.005198069864},
        {  2.779119533275,  0.000000000000,  0.005198069864},
        {  2.779119533275,  3.141592653590,  0.005198069864},
        {  1.208323206480,  0.000000000000,  0.005198069864},
        {  1.208323206480,  3.141592653590,  0.005198069864},
        {  1.933269447110,  0.000000000000,  0.005198069864},
        {  1.933269447110,  3.141592653590,  0.005198069864},
        {  0.362473120315,  1.570796326795,  0.005198069864},
        {  0.362473120315,  4.712388980385,  0.005198069864},
        {  2.779119533275,  1.570796326795,  0.005198069864},
        {  2.779119533275,  4.712388980385,  0.005198069864},
        {  1.208323206480,  1.570796326795,  0.005198069864},
        {  1.208323206480,  4.712388980385,  0.005198069864},
        {  1.933269447110,  1.570796326795,  0.005198069864},
        {  1.933269447110,  4.712388980385,  0.005198069864},
        {  0.566775719462,  1.133830530511,  0.004695720973},
        {  0.566775719462,  2.007762123079,  0.004695720973},
        {  0.566775719462,  5.149354776669,  0.004695720973},
        {  0.566775719462,  4.275423184100,  0.004695720973},
        {  2.574816934128,  1.133830530511,  0.004695720973},
        {  2.574816934128,  2.007762123079,  0.004695720973},
        {  2.574816934128,  5.149354776669,  0.004695720973},
        {  2.574816934128,  4.275423184100,  0.004695720973},
        {  1.062755843178,  1.307707352875,  0.004695720973},
        {  1.062755843178,  1.833885300715,  0.004695720973},
        {  1.062755843178,  4.975477954304,  0.004695720973},
        {  1.062755843178,  4.449300006465,  0.004695720973},
        {  2.078836810412,  1.307707352875,  0.004695720973},
        {  2.078836810412,  1.833885300715,  0.004695720973},
        {  2.078836810412,  4.975477954304,  0.004695720973},
        {  2.078836810412,  4.449300006465,  0.004695720973},
        {  0.566775719462,  0.436965796284,  0.004695720973},
        {  0.566775719462,  2.704626857305,  0.004695720973},
        {  0.566775719462,  5.846219510895,  0.004695720973},
        {  0.566775719462,  3.578558449874,  0.004695720973},
        {  2.574816934128,  0.436965796284,  0.004695720973},
        {  2.574816934128,  2.704626857305,  0.004695720973},
        {  2.574816934128,  5.846219510895,  0.004695720973},
        {  2.574816934128,  3.578558449874,  0.004695720973},
        {  1.341576135361,  1.047737898046,  0.004695720973},
        {  1.341576135361,  2.093854755544,  0.004695720973},
        {  1.341576135361,  5.235447409134,  0.004695720973},
        {  1.341576135361,  4.189330551636,  0.004695720973},
        {  1.800016518229,  1.047737898046,  0.004695720973},
        {  1.800016518229,  2.093854755544,  0.004695720973},
        {  1.800016518229,  5.235447409134,  0.004695720973},
        {  1.800016518229,  4.189330551636,  0.004695720973},
        {  1.062755843178,  0.263088973920,  0.004695720973},
        {  1.062755843178,  2.878503679670,  0.004695720973},
        {  1.062755843178,  6.020096333260,  0.004695720973},
        {  1.062755843178,  3.404681627509,  0.004695720973},
        {  2.078836810412,  0.263088973920,  0.004695720973},
        {  2.078836810412,  2.878503679670,  0.004695720973},
        {  2.078836810412,  6.020096333260,  0.004695720973},
        {  2.078836810412,  3.404681627509,  0.004695720973},
        {  1.341576135361,  0.523058428749,  0.004695720973},
        {  1.341576135361,  2.618534224841,  0.004695720973},
        {  1.341576135361,  5.760126878430,  0.004695720973},
        {  1.341576135361,  3.664651082339,  0.004695720973},
        {  1.800016518229,  0.523058428749,  0.004695720973},
        {  1.800016518229,  2.618534224841,  0.004695720973},
        {  1.800016518229,  5.760126878430,  0.004695720973},
        {  1.800016518229,  3.664651082339,  0.004695720973}    
    };
}