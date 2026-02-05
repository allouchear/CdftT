#include<iostream>
#include<cmath>
#include "../Utils/LM.h"

using namespace std;

LXYZ::LXYZ()
{
    _coef = sqrt(1.0/(4*M_PI));
    _l = vector<int> (3, 0);
}

LXYZ::LXYZ(vector<int> L, double C)
{
    _l=L;
    _coef=C;
}

Zlm::Zlm()
{
    setZlm0();
}

Zlm::Zlm(int L, int M, Binomial& Bin)
{
    _bino=Bin;
    if(L==0 && M == 0)
        setZlm0();
    if(L<0)
        setZlm0();
    if(abs(M)>L)
        setZlm0();

    _l = L;
    _m = M;
    setCoefZlm();
}

void Zlm::setZlm0()
{
    _l = 0;
    _m = 0;
    _numberOfCoefficients = 1;
    _lxyz = vector<LXYZ> (_numberOfCoefficients);
    _lxyz[0] = LXYZ();
}

void Zlm::setCoefZlm()
{
    int Nc = 0;
    double Norm;
    double tmp;
    int l = _l;
    int m = _m;
    int absm = abs(m);
    int t;
    int u;
    int v2;

    /*Norm = sqrt((2*l+1)/(4*PI))*sqrt(factorial(l+absm)/factorial(l-absm))*factorial(absm)/doubleFactorial(2*absm); */
    Norm = sqrt((2*l+1)/(4*M_PI))*sqrt(_bino.fact().factorial(l+absm)*_bino.fact().factorial(l-absm))/_bino.fact().factorial(l)/pow(2.0,absm);

    if(m != 0)
        Norm *= sqrt(2.0);

    for (t=0; t <= (l - absm)/2; t++)
        for (u=0; u<=t; u++)
        {
            int v2m;
            if(m >= 0)
                v2m = 0;
            else
                v2m = 1;
            for(v2 = v2m; v2 <= absm; v2+=2)
                Nc++;
        }

    _numberOfCoefficients = Nc;
    _lxyz = vector<LXYZ> (Nc);
    Nc= -1;
    for (t=0; t <= (l - absm)/2; t++)
    {
        for (u=0; u<=t; u++)
        {
            int v2m;
            if(m >= 0)
                v2m = 0;
            else
                v2m = 1;
            for(v2 = v2m; v2 <= absm; v2+=2)
            {
                Nc++;
                int sign;
                if((t + (v2-v2m)/2)%2)
                    sign = -1;
                else
                    sign = 1;
                tmp = _bino.binomial(l,t)*_bino.binomial(l-t,absm+t)*_bino.binomial(t,u)*_bino.binomial(absm,v2);
                tmp /= pow(4.0,t);

                double C = Norm*tmp*sign;
                vector<int> L (3);
                L[0] = 2*t + absm - 2*u - v2;
                L[1] = 2*u + v2;
                L[2] = l - L[0] - L[1];
                _lxyz[Nc] = LXYZ (L, C);
            }
        }
    }
}

void getlTable(int L, vector<int>& nCoefs, vector<vector<double>>& coefs, vector<vector<vector<int>>>& l, Binomial& Bin)
{
    int m;
    int l1,l2,l3;
    int n;
    /* Spehrical D,F, ...*/
    if(L<-1)
    {
        int klbeg=0, klend=-L, klinc=1;
        int inc;
        int kl;
        int M;
        int c;
        Zlm Stemp;
        int m = -1;
        L = abs(L);
        if(L==1)
        {
            klbeg = L;
            klend = 0;
            klinc = -1;
        }
        else
        {
            klbeg = 0;
             klend = L;
            klinc = +1;
        }
        for(kl = klbeg;(klbeg == 0 && kl<=klend) || (klbeg == L && kl>=klend);kl +=klinc)
        {
            if(kl!=0)
                inc = 2*kl;    
            else
                inc = 1;
            for(M=kl;M>=-kl;M -=inc)
            {
                if(L==1)
                    m = M+abs(L);
                else
                    m++;
                Stemp = Zlm(L,M, Bin);
                nCoefs[m] = Stemp.numberOfCoefficients();
                for(n=0;n<Stemp.numberOfCoefficients();n++)
                {
                    coefs[m][n] = Stemp.lxyz()[n].coef();
                    for(c=0;c<3;c++) 
                        l[c][m][n] = Stemp.lxyz()[n].l()[c];
                }
                if(L==0)
                    break;
            }
            if(L==0)
                break;
        }
        return;
    }
    /* Cartesian S,P,D,F,..*/
    L = abs(L); /* for P */
    for(m=0;m<(L+1)*(L+2)/2;m++) 
    {
        nCoefs[m] = 1;
        coefs[m][0] = 1.0;
    }
    switch(L)
    {
        case 0 :
            m=0;
            l[0][m][0] = 0;l[1][m][0] = 0;l[2][m][0] = 0;
            break;
        case 1 :
            m=0;
            l[0][m][0] = 1;l[1][m][0] = 0;l[2][m][0] = 0; /* X */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 1;l[2][m][0] = 0; /* Y */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 0;l[2][m][0] = 1; /* Z */
            break;
        case 2 :
            m=0;
            l[0][m][0] = 2;l[1][m][0] = 0;l[2][m][0] = 0; /* XX */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 2;l[2][m][0] = 0; /* YY */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 0;l[2][m][0] = 2; /* ZZ */
            m++;
            l[0][m][0] = 1;l[1][m][0] = 1;l[2][m][0] = 0; /* XY */
            m++;
            l[0][m][0] = 1;l[1][m][0] = 0;l[2][m][0] = 1; /* XZ */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 1;l[2][m][0] = 1; /* YZ */
            break;
        case 3 :
            m=0;
            l[0][m][0] = 3;l[1][m][0] = 0;l[2][m][0] = 0; /* XXX */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 3;l[2][m][0] = 0; /* YYY */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 0;l[2][m][0] = 3; /* ZZZ */
            m++;
            l[0][m][0] = 1;l[1][m][0] = 2;l[2][m][0] = 0; /* XYY */
            m++;
            l[0][m][0] = 2;l[1][m][0] = 1;l[2][m][0] = 0; /* XXY */
            m++;
            l[0][m][0] = 2;l[1][m][0] = 0;l[2][m][0] = 1; /* XXZ */
            m++;
            l[0][m][0] = 1;l[1][m][0] = 0;l[2][m][0] = 2; /* XZZ */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 1;l[2][m][0] = 2; /* YZZ */
            m++;
            l[0][m][0] = 0;l[1][m][0] = 2;l[2][m][0] = 1; /* YYZ */
            m++;
            l[0][m][0] = 1;l[1][m][0] = 1;l[2][m][0] = 1; /* XYZ */
            break;
        default :
            m=0;
            for(l3=abs(L);l3>=0;l3--)
                for(l2=abs(L)-l3;l2>=0;l2--)
                {
                     l1 = abs(L)-l2-l3;
                    l[0][m][0] = l1;
                    l[1][m][0] = l2;
                    l[2][m][0] = l3;
                    m++;
                }
    }
}
