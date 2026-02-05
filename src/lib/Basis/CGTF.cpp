#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../Basis/CGTF.h"


CGTF::CGTF()
{
    _gtf.resize(0);
    _numberOfFunctions=0;
    _coefficients=std::vector<double> ();
    _bino = Binomial();
}

CGTF::CGTF(std::vector<GTF> gtfs) : _gtf(gtfs)
{
    _numberOfFunctions=_gtf.size();
    _coefficients=std::vector<double> (_numberOfFunctions, 1);
    _bino = gtfs[0].get_bino();
}

void CGTF::setCoef(double c)
{
    _coefficients.push_back(c);
}

double CGTF::ERICGTF(CGTF &q, CGTF &r, CGTF &s)
{
    int np,nq;
    int nr,ns;
    double sum = 0.0;

    for(np=0;np<_numberOfFunctions;np++)
        for(nq=0;nq<q.numberOfFunctions();nq++)
            for(nr=0;nr<r.numberOfFunctions();nr++)
                for(ns=0;ns<s.numberOfFunctions();ns++)
                    sum += _gtf[np].ERIGTF(q.gtf()[nq],r.gtf()[nr],s.gtf()[ns]); 

    return sum;
}

void CGTF::normaliseCGTF()
{
    int n,np;
    double sum=0.0;

    for(n=0 ; n<_numberOfFunctions ; n++)
        _gtf[n].normaliseRadialGTF();

    for(n=0 ; n<_numberOfFunctions ; n++)
        sum += _coefficients[n]*_coefficients[n]* _gtf[n].overlapGTF(_gtf[n]);

    for(n=0;n<_numberOfFunctions-1 ; n++)
        for(np=n+1; np<_numberOfFunctions; np++)
            sum += 2*_coefficients[n]*_coefficients[np]*_gtf[n].overlapGTF(_gtf[np]);

    if(sum>1.e-20)
    {
        sum = std::sqrt(sum);
        for(n=0 ; n<_numberOfFunctions ; n++)
            _coefficients[n]/=sum;
    }
    else
    {
        std::cout<<"A Contacted Gaussian Type function is nul"<<std::endl;
        exit(1);
    }
}

void CGTF::denormaliseCGTF()
{
    int n;

    for(n=0 ; n<_numberOfFunctions ; n++)
        _gtf[n].denormaliseRadialGTF();
}

double CGTF::overlapCGTF(CGTF &right)
{
    double sum = 0.0;

    //std::cout<<"Number of Function GTF in CGTF = "<<_numberOfFunctions<<std::endl;

    for(int n = 0; n < _numberOfFunctions; ++n)
    {
        for(int np = 0; np < right._numberOfFunctions; ++np)
        {
            sum += _coefficients[n] * right._coefficients[np] * _gtf[n].overlapGTF(right._gtf[np]);
        }
    }

    //std::cout<<"Test overlap in CGTF "<<std::endl;

    return sum;
}

double CGTF::overlap3CGTF(CGTF &midle, CGTF &right)
{
    double sum=0.0;
    int n;
    int np;
    int ns;

    for(n=0;n<_numberOfFunctions;n++)
        for(np=0;np<midle.numberOfFunctions();np++)
            for(ns=0;ns<right.numberOfFunctions();ns++)
                sum += _gtf[n].overlap3GTF(midle.gtf()[np],right.gtf()[ns]);

    return sum;
}

double CGTF::overlap4CGTF(CGTF &middleLeft, CGTF &middleRight, CGTF &right)
{
    double sum=0.0;
    int np;
    int nq;
    int nr;
    int ns;

    for (np = 0; np < _numberOfFunctions; np++)
        for (nq = 0; nq < middleLeft.numberOfFunctions(); nq++)
            for (nr = 0; nr < middleRight.numberOfFunctions(); nr++)
                for (ns = 0; ns < right.numberOfFunctions(); ns++)
                    sum += _gtf[np].overlap4GTF(middleLeft.gtf()[nq], middleRight.gtf()[nr], right.gtf()[ns]);

    return sum;
}

double CGTF::CGTFxyzCGTF(CGTF& right, int ix, int iy, int iz)
{
    double sum=0.0;
    int n;
    int ns;
    std::array<double, 3> C({ 0.0, 0.0, 0.0 });
    std::vector<int> l {ix, iy, iz};
    GTF m(0.0, 1.0, C, l, _bino);

    for(n=0;n<_numberOfFunctions;n++)
        for(ns=0;ns<right.numberOfFunctions();ns++)
            sum += _gtf[n].overlap3GTF(m,right.gtf()[ns]);

    return sum;
}

double CGTF::kineticCGTF(const CGTF& otherCGTF)
{
    double sum = 0.0;

    for(int n = 0; n < _numberOfFunctions; n++)
    {
        for(int np = 0; np < otherCGTF.numberOfFunctions(); ++np)
        {
            sum += _coefficients[n] * otherCGTF._coefficients[np] * _gtf[n].kineticGTF(otherCGTF.gtf()[np]);
        }
    }

    return sum;
}

double CGTF::ionicPotentialCGTF(const CGTF& otherCGTF, const std::array<double, 3>& position, double charge, bool debug)
{
    int n;
    int np;
    double sum = 0.0;

    if (debug)
    {
        std::cout << "Left CGTF:" << std::endl;
        std::cout << *this << std::endl;
        std::cout << "Right CGTF:" << std::endl;
        std::cout << otherCGTF << std::endl;
    }

    for(n = 0; n < _numberOfFunctions; ++n)
    {
        for(np = 0; np < otherCGTF.numberOfFunctions(); ++np)
        {
            sum += _coefficients[n] * otherCGTF._coefficients[np] * _gtf[n].ionicPotentialGTF(otherCGTF.gtf()[np], position, charge, debug);
        }
    }

    return sum;
}

double CGTF::CGTFstarCGTF(CGTF& right)
{
    int n;
    int np;
    double sum=0.0;

    for(n=0;n<_numberOfFunctions;n++)
        for(np=0;np<right._numberOfFunctions;np++)
            sum += _gtf[n].GTFstarGTF(right.gtf()[np]);

    return sum;
}
/*
bool CGTF::CGTFEqCGTF(CGTF& t2)
{
    int i;
    int c;
    if(_numberOfFunctions != t2.numberOfFunctions()) return false;
    for(i=0;i<3;i++)
        if(_gtf[0].l()[i] != t2.gtf()[0].l()[i]) return false;
    for(i=0;i<_numberOfFunctions;i++)
    {
        if(fabs(_gtf[i].exposant()-t2.gtf()[i].exposant())>1e-10) return false;
        if(fabs(_gtf[i].coefficient()-t2.gtf()[i].coefficient())>1e-10) return false;
        for(c=0;c<3;c++)
            if(fabs(_gtf[i].coord()[c]-t2.gtf()[i].coord()[c])>1e-10) return false;
    }
    return true;
}
*/

void CGTF::push_back(GTF& gtf)
{
    _gtf.push_back(gtf);
    _numberOfFunctions++;
}

void CGTF::setNumCenter(int nC)
{
    _num_center=nC;
}

void CGTF::setLtype(std::string s)
{
    _l_type=s;
}

void CGTF::setFactorCoef(double d)
{
    _factor_coef=d;
}

void CGTF::setFormat(std::string format)
{
    _l_format=format;
}

double CGTF::func(double x, double y, double z) const
{
    double r = 0.0;

    for(int i = 0; i < _numberOfFunctions; ++i)
    {
        r += _coefficients[i] * _gtf[i].func(x, y, z);
    }

    return r;
}

bool operator==(const CGTF& left, const CGTF& right)
{
    bool equal = false;

    if (left.gtf().size() == right.gtf().size())
    {
        size_t c = 0;

        for (size_t i = 0; i < left.gtf().size(); ++i)
        {
            for (size_t j = 0; j < left.gtf().size(); ++j)
            {
                if (left.gtf()[i] == right.gtf()[j])
                {
                    ++c;
                }
            }
        }

        if (c == left.gtf().size())
        {
            equal = true;
        }
    }

    return equal;
}

std::ostream& operator<<(std::ostream &stream, const CGTF &cgtf)
{
    for (int i = 0; i < cgtf.numberOfFunctions(); i++)
        stream << std::setw(20) << cgtf.coefficients()[i] << cgtf.gtf()[i] << std::endl;

    return stream;
}

double operator*(const std::vector<CGTF> &cgtfs, const std::vector<double> &coords)
{
    double r=1.0;

    for (size_t i = 0; i < cgtfs.size(); i++)
        r *= cgtfs[i].func(coords[0], coords[1], coords[2]);

    return r;
}
double CGTF::grad_CGTF(const double& x, const double& y, const double& z, int i)
{
    double v=0;
    for(size_t j=0; j<_gtf.size(); j++) 
    {
        v += _coefficients[j]*_gtf[j].grad_GTF(x,y,z,i);
    }
    return v;
}
