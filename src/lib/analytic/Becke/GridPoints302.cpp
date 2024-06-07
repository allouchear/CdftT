#include<iostream>
#include<analytic/Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints302()
{
    _Npts = 302;
    _Lmax = 29;
    _L2max = 14;
    _LebedevGridPoints = {
        // (   theta,            phi,            weight) 
        {  1.570796326795,  0.000000000000,  0.000854591173},
        {  1.570796326795,  3.141592653590,  0.000854591173},
        {  1.570796326795,  1.570796326795,  0.000854591173},
        {  1.570796326795,  4.712388980385,  0.000854591173},
        {  0.000000000000,  3.141592653590,  0.000854591173},
        {  3.141592653590,  3.141592653590,  0.000854591173},
        {  0.955316618125,  0.785398163397,  0.003599119285},
        {  0.955316618125,  2.356194490192,  0.003599119285},
        {  0.955316618125,  5.497787143782,  0.003599119285},
        {  0.955316618125,  3.926990816987,  0.003599119285},
        {  2.186276035465,  0.785398163397,  0.003599119285},
        {  2.186276035465,  2.356194490192,  0.003599119285},
        {  2.186276035465,  5.497787143782,  0.003599119285},
        {  2.186276035465,  3.926990816987,  0.003599119285},
        {  0.520353205918,  0.785398163397,  0.003449788424},
        {  0.520353205918,  2.356194490192,  0.003449788424},
        {  0.520353205918,  5.497787143782,  0.003449788424},
        {  0.520353205918,  3.926990816987,  0.003449788424},
        {  2.621239447672,  0.785398163397,  0.003449788424},
        {  2.621239447672,  2.356194490192,  0.003449788424},
        {  2.621239447672,  5.497787143782,  0.003449788424},
        {  2.621239447672,  3.926990816987,  0.003449788424},
        {  1.211555061487,  1.185820355096,  0.003449788424},
        {  1.211555061487,  1.955772298494,  0.003449788424},
        {  1.211555061487,  5.097364952084,  0.003449788424},
        {  1.211555061487,  4.327413008686,  0.003449788424},
        {  1.930037592103,  1.185820355096,  0.003449788424},
        {  1.930037592103,  1.955772298494,  0.003449788424},
        {  1.930037592103,  5.097364952084,  0.003449788424},
        {  1.930037592103,  4.327413008686,  0.003449788424},
        {  1.211555061487,  0.384975971699,  0.003449788424},
        {  1.211555061487,  2.756616681891,  0.003449788424},
        {  1.211555061487,  5.898209335481,  0.003449788424},
        {  1.211555061487,  3.526568625289,  0.003449788424},
        {  1.930037592103,  0.384975971699,  0.003449788424},
        {  1.930037592103,  2.756616681891,  0.003449788424},
        {  1.930037592103,  5.898209335481,  0.003449788424},
        {  1.930037592103,  3.526568625289,  0.003449788424},
        {  1.190673880275,  0.785398163397,  0.003604822601},
        {  1.190673880275,  2.356194490192,  0.003604822601},
        {  1.190673880275,  5.497787143782,  0.003604822601},
        {  1.190673880275,  3.926990816987,  0.003604822601},
        {  1.950918773315,  0.785398163397,  0.003604822601},
        {  1.950918773315,  2.356194490192,  0.003604822601},
        {  1.950918773315,  5.497787143782,  0.003604822601},
        {  1.950918773315,  3.926990816987,  0.003604822601},
        {  0.854450640997,  0.514328753078,  0.003604822601},
        {  0.854450640997,  2.627263900512,  0.003604822601},
        {  0.854450640997,  5.768856554102,  0.003604822601},
        {  0.854450640997,  3.655921406668,  0.003604822601},
        {  2.287142012592,  0.514328753078,  0.003604822601},
        {  2.287142012592,  2.627263900512,  0.003604822601},
        {  2.287142012592,  5.768856554102,  0.003604822601},
        {  2.287142012592,  3.655921406668,  0.003604822601},
        {  0.854450640997,  1.056467573717,  0.003604822601},
        {  0.854450640997,  2.085125079873,  0.003604822601},
        {  0.854450640997,  5.226717733463,  0.003604822601},
        {  0.854450640997,  4.198060227307,  0.003604822601},
        {  2.287142012592,  1.056467573717,  0.003604822601},
        {  2.287142012592,  2.085125079873,  0.003604822601},
        {  2.287142012592,  5.226717733463,  0.003604822601},
        {  2.287142012592,  4.198060227307,  0.003604822601},
        {  0.732579039339,  0.785398163397,  0.003576729662},
        {  0.732579039339,  2.356194490192,  0.003576729662},
        {  0.732579039339,  5.497787143782,  0.003576729662},
        {  0.732579039339,  3.926990816987,  0.003576729662},
        {  2.409013614251,  0.785398163397,  0.003576729662},
        {  2.409013614251,  2.356194490192,  0.003576729662},
        {  2.409013614251,  5.497787143782,  0.003576729662},
        {  2.409013614251,  3.926990816987,  0.003576729662},
        {  1.078211020817,  1.004259022179,  0.003576729662},
        {  1.078211020817,  2.137333631411,  0.003576729662},
        {  1.078211020817,  5.278926285001,  0.003576729662},
        {  1.078211020817,  4.145851675769,  0.003576729662},
        {  2.063381632772,  1.004259022179,  0.003576729662},
        {  2.063381632772,  2.137333631411,  0.003576729662},
        {  2.063381632772,  5.278926285001,  0.003576729662},
        {  2.063381632772,  4.145851675769,  0.003576729662},
        {  1.078211020817,  0.566537304616,  0.003576729662},
        {  1.078211020817,  2.575055348974,  0.003576729662},
        {  1.078211020817,  5.716648002564,  0.003576729662},
        {  1.078211020817,  3.708129958206,  0.003576729662},
        {  2.063381632772,  0.566537304616,  0.003576729662},
        {  2.063381632772,  2.575055348974,  0.003576729662},
        {  2.063381632772,  5.716648002564,  0.003576729662},
        {  2.063381632772,  3.708129958206,  0.003576729662},
        {  0.136446414324,  0.785398163397,  0.002352101414},
        {  0.136446414324,  2.356194490192,  0.002352101414},
        {  0.136446414324,  5.497787143782,  0.002352101414},
        {  0.136446414324,  3.926990816987,  0.002352101414},
        {  3.005146239266,  0.785398163397,  0.002352101414},
        {  3.005146239266,  2.356194490192,  0.002352101414},
        {  3.005146239266,  5.497787143782,  0.002352101414},
        {  3.005146239266,  3.926990816987,  0.002352101414},
        {  1.474464319498,  1.474014208161,  0.002352101414},
        {  1.474464319498,  1.667578445428,  0.002352101414},
        {  1.474464319498,  4.809171099018,  0.002352101414},
        {  1.474464319498,  4.615606861751,  0.002352101414},
        {  1.667128334092,  1.474014208161,  0.002352101414},
        {  1.667128334092,  1.667578445428,  0.002352101414},
        {  1.667128334092,  4.809171099018,  0.002352101414},
        {  1.667128334092,  4.615606861751,  0.002352101414},
        {  1.474464319498,  0.096782118634,  0.002352101414},
        {  1.474464319498,  3.044810534956,  0.002352101414},
        {  1.474464319498,  6.186403188546,  0.002352101414},
        {  1.474464319498,  3.238374772223,  0.002352101414},
        {  1.667128334092,  0.096782118634,  0.002352101414},
        {  1.667128334092,  3.044810534956,  0.002352101414},
        {  1.667128334092,  6.186403188546,  0.002352101414},
        {  1.667128334092,  3.238374772223,  0.002352101414},
        {  0.319303392340,  0.785398163397,  0.003108953122},
        {  0.319303392340,  2.356194490192,  0.003108953122},
        {  0.319303392340,  5.497787143782,  0.003108953122},
        {  0.319303392340,  3.926990816987,  0.003108953122},
        {  2.822289261250,  0.785398163397,  0.003108953122},
        {  2.822289261250,  2.356194490192,  0.003108953122},
        {  2.822289261250,  5.497787143782,  0.003108953122},
        {  2.822289261250,  3.926990816987,  0.003108953122},
        {  1.346967533858,  1.341139739736,  0.003108953122},
        {  1.346967533858,  1.800452913854,  0.003108953122},
        {  1.346967533858,  4.942045567444,  0.003108953122},
        {  1.346967533858,  4.482732393325,  0.003108953122},
        {  1.794625119732,  1.341139739736,  0.003108953122},
        {  1.794625119732,  1.800452913854,  0.003108953122},
        {  1.794625119732,  4.942045567444,  0.003108953122},
        {  1.794625119732,  4.482732393325,  0.003108953122},
        {  1.346967533858,  0.229656587059,  0.003108953122},
        {  1.346967533858,  2.911936066530,  0.003108953122},
        {  1.346967533858,  6.053528720120,  0.003108953122},
        {  1.346967533858,  3.371249240649,  0.003108953122},
        {  1.794625119732,  0.229656587059,  0.003108953122},
        {  1.794625119732,  2.911936066530,  0.003108953122},
        {  1.794625119732,  6.053528720120,  0.003108953122},
        {  1.794625119732,  3.371249240649,  0.003108953122},
        {  1.441195151732,  0.785398163397,  0.003650045808},
        {  1.441195151732,  2.356194490192,  0.003650045808},
        {  1.441195151732,  5.497787143782,  0.003650045808},
        {  1.441195151732,  3.926990816987,  0.003650045808},
        {  1.700397501858,  0.785398163397,  0.003650045808},
        {  1.700397501858,  2.356194490192,  0.003650045808},
        {  1.700397501858,  5.497787143782,  0.003650045808},
        {  1.700397501858,  3.926990816987,  0.003650045808},
        {  0.793749869014,  0.182271146595,  0.003650045808},
        {  0.793749869014,  2.959321506995,  0.003650045808},
        {  0.793749869014,  6.100914160585,  0.003650045808},
        {  0.793749869014,  3.323863800185,  0.003650045808},
        {  2.347842784576,  0.182271146595,  0.003650045808},
        {  2.347842784576,  2.959321506995,  0.003650045808},
        {  2.347842784576,  6.100914160585,  0.003650045808},
        {  2.347842784576,  3.323863800185,  0.003650045808},
        {  0.793749869014,  1.388525180200,  0.003650045808},
        {  0.793749869014,  1.753067473390,  0.003650045808},
        {  0.793749869014,  4.894660126980,  0.003650045808},
        {  0.793749869014,  4.530117833790,  0.003650045808},
        {  2.347842784576,  1.388525180200,  0.003650045808},
        {  2.347842784576,  1.753067473390,  0.003650045808},
        {  2.347842784576,  4.894660126980,  0.003650045808},
        {  2.347842784576,  4.530117833790,  0.003650045808},
        {  1.570796326795,  1.303198744718,  0.002982344963},
        {  1.570796326795,  1.838393908872,  0.002982344963},
        {  1.570796326795,  4.979986562462,  0.002982344963},
        {  1.570796326795,  4.444791398308,  0.002982344963},
        {  1.570796326795,  0.267597582077,  0.002982344963},
        {  1.570796326795,  2.873995071513,  0.002982344963},
        {  1.570796326795,  6.015587725103,  0.002982344963},
        {  1.570796326795,  3.409190235667,  0.002982344963},
        {  0.267597582077,  0.000000000000,  0.002982344963},
        {  0.267597582077,  3.141592653590,  0.002982344963},
        {  2.873995071513,  0.000000000000,  0.002982344963},
        {  2.873995071513,  3.141592653590,  0.002982344963},
        {  1.303198744718,  0.000000000000,  0.002982344963},
        {  1.303198744718,  3.141592653590,  0.002982344963},
        {  1.838393908872,  0.000000000000,  0.002982344963},
        {  1.838393908872,  3.141592653590,  0.002982344963},
        {  0.267597582077,  1.570796326795,  0.002982344963},
        {  0.267597582077,  4.712388980385,  0.002982344963},
        {  2.873995071513,  1.570796326795,  0.002982344963},
        {  2.873995071513,  4.712388980385,  0.002982344963},
        {  1.303198744718,  1.570796326795,  0.002982344963},
        {  1.303198744718,  4.712388980385,  0.002982344963},
        {  1.838393908872,  1.570796326795,  0.002982344963},
        {  1.838393908872,  4.712388980385,  0.002982344963},
        {  1.570796326795,  0.961981553560,  0.003600820932},
        {  1.570796326795,  2.179611100030,  0.003600820932},
        {  1.570796326795,  5.321203753620,  0.003600820932},
        {  1.570796326795,  4.103574207150,  0.003600820932},
        {  1.570796326795,  0.608814773235,  0.003600820932},
        {  1.570796326795,  2.532777880355,  0.003600820932},
        {  1.570796326795,  5.674370533945,  0.003600820932},
        {  1.570796326795,  3.750407426825,  0.003600820932},
        {  0.608814773235,  0.000000000000,  0.003600820932},
        {  0.608814773235,  3.141592653590,  0.003600820932},
        {  2.532777880355,  0.000000000000,  0.003600820932},
        {  2.532777880355,  3.141592653590,  0.003600820932},
        {  0.961981553560,  0.000000000000,  0.003600820932},
        {  0.961981553560,  3.141592653590,  0.003600820932},
        {  2.179611100030,  0.000000000000,  0.003600820932},
        {  2.179611100030,  3.141592653590,  0.003600820932},
        {  0.608814773235,  1.570796326795,  0.003600820932},
        {  0.608814773235,  4.712388980385,  0.003600820932},
        {  2.532777880355,  1.570796326795,  0.003600820932},
        {  2.532777880355,  4.712388980385,  0.003600820932},
        {  0.961981553560,  1.570796326795,  0.003600820932},
        {  0.961981553560,  4.712388980385,  0.003600820932},
        {  2.179611100030,  1.570796326795,  0.003600820932},
        {  2.179611100030,  4.712388980385,  0.003600820932},
        {  0.994564953323,  1.266795091971,  0.003571540554},
        {  0.994564953323,  1.874797561619,  0.003571540554},
        {  0.994564953323,  5.016390215209,  0.003571540554},
        {  0.994564953323,  4.408387745561,  0.003571540554},
        {  2.147027700266,  1.266795091971,  0.003571540554},
        {  2.147027700266,  1.874797561619,  0.003571540554},
        {  2.147027700266,  5.016390215209,  0.003571540554},
        {  2.147027700266,  4.408387745561,  0.003571540554},
        {  0.643379849978,  1.139105867439,  0.003571540554},
        {  0.643379849978,  2.002486786151,  0.003571540554},
        {  0.643379849978,  5.144079439741,  0.003571540554},
        {  0.643379849978,  4.280698521029,  0.003571540554},
        {  2.498212803612,  1.139105867439,  0.003571540554},
        {  2.498212803612,  2.002486786151,  0.003571540554},
        {  2.498212803612,  5.144079439741,  0.003571540554},
        {  2.498212803612,  4.280698521029,  0.003571540554},
        {  0.994564953323,  0.304001234824,  0.003571540554},
        {  0.994564953323,  2.837591418766,  0.003571540554},
        {  0.994564953323,  5.979184072356,  0.003571540554},
        {  0.994564953323,  3.445593888414,  0.003571540554},
        {  2.147027700266,  0.304001234824,  0.003571540554},
        {  2.147027700266,  2.837591418766,  0.003571540554},
        {  2.147027700266,  5.979184072356,  0.003571540554},
        {  2.147027700266,  3.445593888414,  0.003571540554},
        {  1.317079548059,  0.597875683843,  0.003571540554},
        {  1.317079548059,  2.543716969747,  0.003571540554},
        {  1.317079548059,  5.685309623336,  0.003571540554},
        {  1.317079548059,  3.739468337433,  0.003571540554},
        {  1.824513105531,  0.597875683843,  0.003571540554},
        {  1.824513105531,  2.543716969747,  0.003571540554},
        {  1.824513105531,  5.685309623336,  0.003571540554},
        {  1.824513105531,  3.739468337433,  0.003571540554},
        {  0.643379849978,  0.431690459356,  0.003571540554},
        {  0.643379849978,  2.709902194234,  0.003571540554},
        {  0.643379849978,  5.851494847823,  0.003571540554},
        {  0.643379849978,  3.573283112946,  0.003571540554},
        {  2.498212803612,  0.431690459356,  0.003571540554},
        {  2.498212803612,  2.709902194234,  0.003571540554},
        {  2.498212803612,  5.851494847823,  0.003571540554},
        {  2.498212803612,  3.573283112946,  0.003571540554},
        {  1.317079548059,  0.972920642952,  0.003571540554},
        {  1.317079548059,  2.168672010638,  0.003571540554},
        {  1.317079548059,  5.310264664228,  0.003571540554},
        {  1.317079548059,  4.114513296542,  0.003571540554},
        {  1.824513105531,  0.972920642952,  0.003571540554},
        {  1.824513105531,  2.168672010638,  0.003571540554},
        {  1.824513105531,  5.310264664228,  0.003571540554},
        {  1.824513105531,  4.114513296542,  0.003571540554},
        {  0.445390437897,  1.280399762345,  0.003392312205},
        {  0.445390437897,  1.861192891244,  0.003392312205},
        {  0.445390437897,  5.002785544834,  0.003392312205},
        {  0.445390437897,  4.421992415935,  0.003392312205},
        {  2.696202215693,  1.280399762345,  0.003392312205},
        {  2.696202215693,  1.861192891244,  0.003392312205},
        {  2.696202215693,  5.002785544834,  0.003392312205},
        {  2.696202215693,  4.421992415935,  0.003392312205},
        {  1.145300544069,  1.434948238656,  0.003392312205},
        {  1.145300544069,  1.706644414933,  0.003392312205},
        {  1.145300544069,  4.848237068523,  0.003392312205},
        {  1.145300544069,  4.576540892246,  0.003392312205},
        {  1.996292109521,  1.434948238656,  0.003392312205},
        {  1.996292109521,  1.706644414933,  0.003392312205},
        {  1.996292109521,  4.848237068523,  0.003392312205},
        {  1.996292109521,  4.576540892246,  0.003392312205},
        {  0.445390437897,  0.290396564450,  0.003392312205},
        {  0.445390437897,  2.851196089140,  0.003392312205},
        {  0.445390437897,  5.992788742730,  0.003392312205},
        {  0.445390437897,  3.431989218039,  0.003392312205},
        {  2.696202215693,  0.290396564450,  0.003392312205},
        {  2.696202215693,  2.851196089140,  0.003392312205},
        {  2.696202215693,  5.992788742730,  0.003392312205},
        {  2.696202215693,  3.431989218039,  0.003392312205},
        {  1.447126475279,  1.141810028851,  0.003392312205},
        {  1.447126475279,  1.999782624738,  0.003392312205},
        {  1.447126475279,  5.141375278328,  0.003392312205},
        {  1.447126475279,  4.283402682441,  0.003392312205},
        {  1.694466178311,  1.141810028851,  0.003392312205},
        {  1.694466178311,  1.999782624738,  0.003392312205},
        {  1.694466178311,  5.141375278328,  0.003392312205},
        {  1.694466178311,  4.283402682441,  0.003392312205},
        {  1.145300544069,  0.135848088139,  0.003392312205},
        {  1.145300544069,  3.005744565451,  0.003392312205},
        {  1.145300544069,  6.147337219041,  0.003392312205},
        {  1.145300544069,  3.277440741728,  0.003392312205},
        {  1.996292109521,  0.135848088139,  0.003392312205},
        {  1.996292109521,  3.005744565451,  0.003392312205},
        {  1.996292109521,  6.147337219041,  0.003392312205},
        {  1.996292109521,  3.277440741728,  0.003392312205},
        {  1.447126475279,  0.428986297943,  0.003392312205},
        {  1.447126475279,  2.712606355646,  0.003392312205},
        {  1.447126475279,  5.854199009236,  0.003392312205},
        {  1.447126475279,  3.570578951533,  0.003392312205},
        {  1.694466178311,  0.428986297943,  0.003392312205},
        {  1.694466178311,  2.712606355646,  0.003392312205},
        {  1.694466178311,  5.854199009236,  0.003392312205},
        {  1.694466178311,  3.570578951533,  0.003392312205}    
    };
}