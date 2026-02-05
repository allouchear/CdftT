#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../Basis/GTF.h"


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

GTF::GTF():
    _exponent(0.0),
    _coefficient(0.0),
    _coord(),
    _l(),
    _bino()
{ }

GTF::GTF(const double exponent, const double coefficient, const std::array<double, 3>& coord, const std::vector<int>& l, const Binomial& binomial):
    _exponent(exponent),
    _coefficient(coefficient),
    _coord(coord),
    _l(l),
    _bino(binomial)
{ }


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

double GTF::get_exponent() const
{
    return _exponent;
}

double GTF::get_coefficient() const
{
    return _coefficient;
}

const std::array<double, 3>& GTF::get_coord() const
{
    return _coord;
}

const std::vector<int>& GTF::get_l() const
{
    return _l;
}

Binomial& GTF::get_bino()
{
    return _bino;
}






double GTF::GTFstarGTF (GTF& right)
{
    int i,j;
    std::vector<double> sum(3,0.0);
    double t;
    std::vector<double> PA(3);
    std::vector<double> PB(3);
    double gama=_exponent+right._exponent;
    double R2=0.0;
    double c=0;

    for(j=0; j<3; j++)
    {
        t=(_exponent*_coord[j] + right._exponent*right._coord[j])/gama;
        PA[j]=_coord[j]-t;
        PB[j]=right._coord[j]-t;
        R2+=(_coord[j]-right._coord[j])*(_coord[j]-right._coord[j]);
    }

    c = (M_PI/gama)*std::sqrt(M_PI/gama)*exp(-_exponent*right._exponent/gama*R2);

    for(j=0; j<3; j++)
    {
        for(i=0; i<=(_l[j]+right._l[j])/2; i++)
        {
            sum[j] += f(2 * i, _l[j], right._l[j], PA[j], PB[j], _bino) * _bino.fact().double_factorial(2 * i - 1) / (power(2.0, i) * power(gama, i));
        }
    }

    return c*sum[0]*sum[1]*sum[2];
}

double GTF::GTFstarGTFstarGTF (GTF& mid, GTF& right)
{
    std::vector<double> sum(3);
    double t;
    std::vector<double> PA(3);
    std::vector<double> PB(3);
    std::vector<double> P(3);
    std::vector<double> QP(3);
    std::vector<double> QC(3);
    double gama1=_exponent+right.get_exponent();
    double gama=gama1+mid.get_exponent();
    double R2AB=0.0;
    double R2PC=0.0;
    double c = 0;
    int iAB;
    int i,j;

    for(j=0;j<3;j++)
    {
        t=(_exponent*_coord[j]+right.get_exponent()*right.get_coord()[j])/gama1;
        P[j]=t;
        PA[j]=_coord[j]-t;
        PB[j]=right.get_coord()[j]-t;
        R2AB += (_coord[j]-right.get_coord()[j])*(_coord[j]-right.get_coord()[j]);
    }
    for(j=0;j<3;j++)
    {
        t=(gama1*P[j]+mid.get_exponent()*mid.get_coord()[j])/gama;
        QP[j]=P[j]-t;
        QC[j]=mid.get_coord()[j]-t;
        R2PC += (P[j]-mid.get_coord()[j])*(P[j]-mid.get_coord()[j]);
    }
    c = (M_PI/gama)*std::sqrt(M_PI/gama)*exp(-_exponent*right.get_exponent()/gama1*R2AB)*exp(-gama1*mid.get_exponent()/gama*R2PC);

    for(j=0;j<3;j++)
    {
        sum[j]=0.0;
        for(iAB=0;iAB<=(_l[j]+right.get_l()[j]);iAB++)
        {
            double fiAB = f(iAB,_l[j],right.get_l()[j],PA[j],PB[j], _bino);
            for(i=0;i<=(iAB+mid.get_l()[j])/2;i++)
            {
                sum[j] +=
                fiAB*
                f(2*i,iAB,mid.get_l()[j],QP[j],QC[j], _bino)*
                _bino.fact().double_factorial(2*i-1)/(power(2.0,i)*power(gama,i));
             }
        }
    }
    return  c*sum[0]*sum[1]*sum[2];
}

double GTF::GTFstarGTFstarGTFstarGTF(GTF& B, GTF& C, GTF& D)
{
    std::vector<double> sum(3);
    double t;
    std::vector<double> PA(3);
    std::vector<double> PB(3);
    std::vector<double> QC(3);
    std::vector<double> QD(3);
    std::vector<double> P(3);
    std::vector<double> Q(3);
    std::vector<double> GP(3);
    std::vector<double> GQ(3);
    double gama1=_exponent+B.get_exponent();
    double gama2=C.get_exponent()+D.get_exponent();
    double gama=gama1+gama2;
    double R2AB=0.0;
    double R2CD=0.0;
    double R2PQ=0.0;
    double c = 0;
    int iAB;
    int iCD;
    int i,j;

    for(j=0;j<3;j++)
    {
        t=(_exponent*_coord[j]+B.get_exponent()*B.get_coord()[j])/gama1;
        P[j]=t;
        PA[j]=_coord[j]-t;
        PB[j]=B.get_coord()[j]-t;
        R2AB += (_coord[j]-B.get_coord()[j])*(_coord[j]-B.get_coord()[j]);
    }
    for(j=0;j<3;j++)
    {
        t=(C.get_exponent()*C.get_coord()[j]+D.get_exponent()*D.get_coord()[j])/gama2;
        Q[j]=t;
        QC[j]=C.get_coord()[j]-t;
        QD[j]=D.get_coord()[j]-t;
        R2CD += (C.get_coord()[j]-D.get_coord()[j])*(C.get_coord()[j]-D.get_coord()[j]);
    }
    for(j=0;j<3;j++)
    {
        t=(gama1*P[j]+gama2*Q[j])/gama;
        GP[j]=P[j]-t;
        GQ[j]=Q[j]-t;
        R2PQ += (P[j]-Q[j])*(P[j]-Q[j]);
    }
    c = (M_PI/gama)*std::sqrt(M_PI/gama)
        *exp(-_exponent*B.get_exponent()/gama1*R2AB)
        *exp(-C.get_exponent()*D.get_exponent()/gama2*R2CD)
        *exp(-gama1*gama2/gama*R2PQ);


    for(j=0;j<3;j++)
    {
        sum[j]=0.0;
        for(iAB=0;iAB<=(_l[j]+B.get_l()[j]);iAB++)
        {
            double fiAB = f(iAB,_l[j],B.get_l()[j],PA[j],PB[j],_bino);
            for(iCD=0;iCD<=(C.get_l()[j]+D.get_l()[j]);iCD++)
            {
                double fiCD = f(iCD,C.get_l()[j],D.get_l()[j],QC[j],QD[j],_bino);
                for(i=0;i<=(iAB+iCD)/2;i++)
                {
                    sum[j] +=
                    fiAB*
                    fiCD*
                    f(2*i,iAB,iCD,GP[j],GQ[j],_bino)*
                    _bino.fact().double_factorial(2*i-1)/(power(2.0,i)*power(gama,i));
                 }
            }
        }
    }
    return  c*sum[0]*sum[1]*sum[2];
}

double GTF::normeGTF()
{
    return std::sqrt(2 * _exponent / M_PI * std::sqrt(2 * _exponent / M_PI)
                     * power(4 * _exponent, _l[0] + _l[1] + _l[2])
                     / (_bino.fact().double_factorial(_l[0])
                        * _bino.fact().double_factorial(_l[1])
                        * _bino.fact().double_factorial(_l[2])));

    // Note Ludo : error in the formula? See. MOTEC page 420, eq. (A4)
    // Should be result below instead?

    // return std::sqrt(2 * _exponent / M_PI * std::sqrt(2 * _exponent / M_PI)
    //                  * power(4 * _exponent, _l[0] + _l[1] + _l[2])
    //                  / (_bino.fact().double_factorial(2 * _l[0] - 1)
    //                     * _bino.fact().double_factorial(2 * _l[1] - 1)
    //                     * _bino.fact().double_factorial(2 * _l[2] - 1)));
}

double GTF::normeGTF(GTF& q)
{
    return std::sqrt(2 * q.get_exponent() / M_PI * std::sqrt(2 * q.get_exponent() / M_PI)
                     * power(4 * q.get_exponent(), q.get_l()[0] + q.get_l()[1] + q.get_l()[2])
                     / (q.get_bino().fact().double_factorial(q.get_l()[0])
                        * q.get_bino().fact().double_factorial(q.get_l()[1])
                        * q.get_bino().fact().double_factorial(q.get_l()[2])));

    // Note Ludo : error in the formula? See. MOTEC page 420, eq. (A4)
    // Should be result below instead?

    // return std::sqrt(2 * q.get_exponent() / M_PI * std::sqrt(2 * q.get_exponent() / M_PI)
    //                  * power(4 * q.get_exponent(), q.get_l()[0] + q.get_l()[1] + q.get_l()[2])
    //                  / (q.get_bino().fact().double_factorial(2 * q.get_l()[0] - 1)
    //                     * q.get_bino().fact().double_factorial(2 * q.get_l()[1] - 1)
    //                     * q.get_bino().fact().double_factorial(2 * q.get_l()[2] - 1)));
}

void GTF::normaliseRadialGTF()
{
    int l_bis=_l[0]+_l[1]+_l[2];
    std::vector<int> l (3);
    l[0]=l_bis;
    l[1]=0;
    l[2]=0;
    GTF q(_exponent, _coefficient, _coord, l, _bino);
    _coefficient*=normeGTF(q);
}

void GTF::denormaliseRadialGTF()
{
    int l_bis=_l[0]+_l[1]+_l[2];
    std::vector<int> l (3);
    l[0]=l_bis;
    l[1]=0;
    l[2]=0;
    GTF q(_exponent, _coefficient, _coord, l, _bino);
    _coefficient/=normeGTF(q);
}

void GTF::normaliseGTF()
{
    _coefficient*=normeGTF();
}

double GTF::overlapGTF(GTF& right)
{
    //cout<<"Test overlap in GTF"<<endl;
    return _coefficient * right._coefficient * GTFstarGTF(right);
}

double GTF::overlap3GTF(GTF& mid, GTF& right)
{
    return _coefficient*mid.get_coefficient()*right.get_coefficient()*GTFstarGTFstarGTF(mid, right);
}

double GTF::overlap4GTF(GTF& B, GTF& C, GTF& D)
{
    return _coefficient*B.get_coefficient()*C.get_coefficient()*D.get_coefficient()*GTFstarGTFstarGTFstarGTF(B,C,D);
}

double GTF::GTFxyzGTF(GTF& q, int ix, int iy, int iz)
{
    std::array<double, 3> C({ 0.0, 0.0, 0.0 });
    std::vector<int> l{ix, iy, iz};
    GTF m(0.0, 1.0, C, l, _bino);
    return overlap3GTF(m, q);
}

double GTF::kineticGTF(const GTF& right)
{
    int j;
    GTF b = right;
    std::vector<double> Ti(7);
    std::vector<double> sum(3);
    double T=0.0;

    for(j=0; j<7; j++)
    {
        Ti[j]=0.0;
    }

    Ti[0] = GTFstarGTF(b);

    b.setL(0, right.get_l()[0] + 2);
    b.setL(1, right.get_l()[1]);
    b.setL(2, right.get_l()[2]);

    Ti[1] = GTFstarGTF(b);

    b.setL(0, right.get_l()[0]);
    b.setL(1, right.get_l()[1] + 2);
    b.setL(2, right.get_l()[2]);

    Ti[2] = GTFstarGTF(b);

    b.setL(0, right.get_l()[0]);
    b.setL(1, right.get_l()[1]);
    b.setL(2, right.get_l()[2] + 2);

    Ti[3] = GTFstarGTF(b);

    b.setL(0, right.get_l()[0] - 2);
    b.setL(1, right.get_l()[1]);
    b.setL(2, right.get_l()[2]);

    Ti[4] = GTFstarGTF(b);

    b.setL(0, right.get_l()[0]);
    b.setL(1, right.get_l()[1] - 2);
    b.setL(2, right.get_l()[2]);

    Ti[5] = GTFstarGTF(b);

    b.setL(0, right.get_l()[0]);
    b.setL(1, right.get_l()[1]);
    b.setL(2, right.get_l()[2] - 2);

    Ti[6] = GTFstarGTF(b);

    sum[0] = right.get_exponent() * (2 * (right.get_l()[0] + right.get_l()[1] + right.get_l()[2]) + 3) * Ti[0];
    sum[1] = 0.0;

    for (j = 1; j <= 3; ++j)
    {
        sum[1] += Ti[j];
    }

    sum[1] = -2 * right.get_exponent() * right.get_exponent() * sum[1];
    sum[2] = 0.0;

    for(j = 4; j <= 6; ++j)
    {
        sum[2] += right.get_l()[j - 4] * (right.get_l()[j - 4] - 1) * Ti[j];
    }

    sum[2] *= -0.5;

    for(j = 0; j < 3; ++j)
    {
        T += sum[j];
    }

    return T * _coefficient * right.get_coefficient();
}

double GTF::ionicPotentialGTF(const GTF& right, const std::array<double, 3>& C, double Z, bool debug)
{
    int i, r, u;
    int j, s, v;
    int k, t, w;

    double Sx, Sy, Sz;

    std::array<double, 3> PA;
    std::array<double, 3> PB;
    std::array<double, 3> CP;
    double P;

    double gamma = _exponent + right.get_exponent();
    double AB2 = 0.0;
    double CP2 = 0.0;

    for(j = 0; j < 3; ++j)
    {
        P = (_exponent * _coord[j] + right.get_exponent() * right.get_coord()[j]) / gamma;

        PA[j] = _coord[j] - P;
        PB[j] = right.get_coord()[j] - P;
        CP[j] = P - C[j];

        AB2 += (_coord[j] - right.get_coord()[j]) * (_coord[j] - right.get_coord()[j]);
        CP2 += CP[j] * CP[j];
    }

    std::vector<double> FTable = getFTable(_l[0] + right.get_l()[0] + _l[1] + right.get_l()[1] + _l[2] + right.get_l()[2], gamma * CP2, debug);

    if (debug)
    {
        std::cout << "FTable(" << _l[0] << " + " << right.get_l()[0] << " + " << _l[1] << " + " << right.get_l()[1] << " + " << _l[2] << " + " << right.get_l()[2] << ", " << gamma << " * " << CP2 << ") values:" << std::endl;
        for (size_t bhj = 0; bhj < FTable.size(); ++bhj)
        {
            std::cout << "FTable[" << bhj << "] = " << FTable[bhj] << std::endl;
        }
        std::cout << std::endl;
    }

    double sum = 0.0;
    for(i = 0; i <= _l[0] + right.get_l()[0]; ++i)
    {
        for(r = 0; r <= i / 2; ++r)
        {
            for(u = 0; u <= (i - 2 * r) / 2; ++u)
            {
                Sx = A(i, r, u, _l[0], right.get_l()[0], PA[0], PB[0], CP[0], gamma, _bino);

                for(j = 0; j <= _l[1] + right.get_l()[1]; ++j)
                {
                    for(s = 0; s <= j / 2; ++s)
                    {
                        for(v = 0; v <= (j - 2 * s) / 2; ++v)
                        {
                            Sy = A(j, s, v, _l[1], right.get_l()[1], PA[1], PB[1], CP[1], gamma, _bino);

                            for(k = 0; k <= _l[2] + right.get_l()[2]; ++k)
                            {
                                for(t = 0; t <= k / 2; ++t)
                                {
                                    for(w = 0; w <= (k - 2 * t) / 2; ++w)
                                    {
                                        Sz = A(k, t, w, _l[2], right.get_l()[2], PA[2], PB[2], CP[2], gamma, _bino);

                                        sum += Sx * Sy * Sz * FTable[i + j + k - 2 * (r + s + t) - u - v - w];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    sum *= 2 * M_PI / gamma * std::exp(- _exponent * right.get_exponent() * AB2 / gamma) * _coefficient * right.get_coefficient();

    return -Z * sum;
}

double GTF::ERIGTF(GTF& q, GTF& r, GTF& s)    
{
    int I,Ip,R,Rp,U;
    int J,Jp,S,Sp,N;
    int K,Kp,T,Tp,W;
    double Sx,Sy,Sz;
    std::vector<std::vector<double>> Te;
    Te.resize(2, std::vector<double>(3));
    double temp1,temp2;
    std::vector<double> PA(3);
    std::vector<double> PB(3);

    std::vector<double> QC(3);
    std::vector<double> QD(3);
    std::vector<double> PQ(3);

    
    double g1=_exponent+q.get_exponent();
    double g2=r.get_exponent()+s.get_exponent();
    double d=(1.0/g1+1.0/g2)/4;
    double RAB2=0.0;
    double RCD2=0.0;
    double RPQ2=0.0;
    int j;
    
    double sum;

    for(j=0; j<3; j++)
    {
        temp1=(_exponent*_coord[j]+q.get_exponent()*q.get_coord()[j])/g1;
        PA[j]=_coord[j]-temp1;
        PB[j]=q.get_coord()[j]-temp1;
        
        temp2=(r.get_exponent()*r.get_coord()[j]+s.get_exponent()*s.get_coord()[j])/g2;
        QC[j]=r.get_coord()[j]-temp2;
        QD[j]=s.get_coord()[j]-temp2;

        PQ[j]=temp2-temp1;

        RAB2+=(_coord[j]-q.get_coord()[j])*(_coord[j]-q.get_coord()[j]);
        RCD2+=(r.get_coord()[j]-s.get_coord()[j])*(r.get_coord()[j]-s.get_coord()[j]);
        RPQ2+=PQ[j]*PQ[j];
    }

    sum=0.0;
    for(I=0; I<=_l[0]+q.get_l()[0]; I++)
        for(R=0; R<=I/2; R++)
        {
            Te[0][0]=Theta(I,R,_l[0],q.get_l()[0],PA[0],PB[0],g1,_bino);
            for(Ip=0; Ip<=r.get_l()[0]+s.get_l()[0]; Ip++)
                for(Rp=0; Rp<=Ip/2; Rp++)
                {
                    Te[1][0]=Theta(Ip,Rp,r.get_l()[0],s.get_l()[0],QC[0],QD[0],g2,_bino);
                    for(U=0; U<=(I+2*Ip)/2-R-Rp; U++)
                    {
                        Sx=B(I,Ip,R,Rp,U,PQ[0],d,Te[0][0],Te[1][0],_bino.fact());
                        for(J=0; J<=_l[1]+q.get_l()[1]; J++)
                            for(S=0; S<=J/2; S++)
                            {
                                Te[0][1]=Theta(J,S,_l[1],q.get_l()[1],PA[1],PB[1],g1,_bino);
                                for(Jp=0; Jp<=r.get_l()[1]+s.get_l()[1]; Jp++)
                                    for(Sp=0; Sp<=Jp/2; Sp++)
                                    {
                                        Te[1][1]=Theta(Jp,Sp,r.get_l()[1],s.get_l()[1],QC[1],QD[1],g2,_bino);
                                        for(N=0; N<=(J+Jp)/2-S-Sp; N++)
                                        {
                                            Sy=B(J,Jp,S,Sp,N,PQ[1],d,Te[0][1],Te[1][1],_bino.fact());
                                            for(K=0; K<=_l[2]+q.get_l()[2]; K++)
                                                for(T=0; T<=K/2; T++)
                                                {
                                                    Te[0][2]=Theta(K,T,_l[2],q.get_l()[2],PA[2],PB[2],g1,_bino);
                                                    for(Kp=0; Kp<=r.get_l()[2]+s.get_l()[2]; Kp++)
                                                        for(Tp=0; Tp<=Kp/2; Tp++)
                                                        {
                                                            Te[1][2]=Theta(Kp,Tp,r.get_l()[2],s.get_l()[2],QC[2],QD[2],g2,_bino);
                                                            for(W=0; W<=(K+Kp)/2-T-Tp; W++)
                                                            {
                                                                Sz = B(K, Kp, T, Tp, W, PQ[2], d, Te[0][2], Te[1][2], _bino.fact());
                                                                sum += Sx * Sy * Sz * F(I + Ip + J + Jp + K + Kp - 2 * (R + Rp + S + Sp + T + Tp) - U - N - W, RPQ2 / 4.0 / d);
                                                            }
                                                        }

                                                }
                                        }
                                    }
                            }
                    }
                }
        }

    sum*=2*M_PI*M_PI*std::sqrt(M_PI)/g1/g2/std::sqrt(g1+g2)*exp(-_exponent*q.get_exponent()*RAB2/g1)*exp(-r.get_exponent()*s.get_exponent()*RCD2/g2)*
            _coefficient*q.get_coefficient()*r.get_coefficient()*s.get_coefficient();
    return sum;
}

double GTF::func(double x, double y, double z) const
{
    double xA = x - _coord[0];
    double yA = y - _coord[1];
    double zA = z - _coord[2];

    return _coefficient * power(xA, _l[0]) * power(yA, _l[1]) * power(zA, _l[2]) * exp(- _exponent * (xA * xA + yA * yA + zA * zA));
}

void GTF::operator/=(double c)
{
    _coefficient /= c;
}

void GTF::operator*=(double c)
{
    _coefficient *= c;
}

void GTF::push_back(const double& a, const double& c, const std::array<double, 3>& coord, const std::vector<int>& l, Binomial& B)
{
    _exponent=a;
    _coefficient=c;
    _coord=coord;
    _l=l;
    _bino=B;
}

bool operator==(GTF a, GTF b)
{
    size_t i;
    bool n=true;
    if(abs(a.get_exponent()-b.get_exponent())>10e-10)
        n=false;
    if(abs(a.get_coefficient()-b.get_coefficient())>10e-10)
        n=false;
    for(i=0; i<3; i++)
        if(abs(a.get_l()[i]-b.get_l()[i])>10e10-10 || abs(a.get_coord()[i]-b.get_coord()[i])>10e-10)
            n=false;
    return n;
}

std::ostream& operator<<(std::ostream& flux, const GTF& gtf)
{
    flux << std::left << std::setw(20) << gtf.get_coefficient() << std::setw(20) << gtf.get_exponent() << std::setw(5) << gtf.get_l()[0] << std::setw(5) << gtf.get_l()[1]
    << std::setw(5) << gtf.get_l()[2] << std::setw(20) << gtf.get_coord()[0] << std::setw(20) << gtf.get_coord()[1] << std::setw(20) << gtf.get_coord()[2];
    return flux;
}

double operator*(const std::vector<GTF>& a, const std::vector<double>& b)
{
    double r=1.0;

    for(size_t i=0; i<a.size(); i++)
        r*=a[i].func(b[0],b[1],b[2]);
    
    return r;
}
double GTF::grad_GTF(const double& x, const double& y, const double& z, int i)
{
    double xA[3] = {x-_coord[0], y-_coord[1], z-_coord[2]};
    int k=(i+1)%3;
    int j=(i+2)%3;
    double expo =(xA[i]*xA[i]+xA[j]*xA[j]+xA[k]*xA[k])*_exponent;
    if(expo > 40) return 1e-14;
    if(_l[i]==0)
        return -2*_exponent*xA[i]*_coefficient*power(xA[k], _l[k])*power(xA[j], _l[j])*exp(-expo);
    else
        return _coefficient*exp(-expo)*power(xA[i], _l[i]-1)*power(xA[k], _l[k])*power(xA[j], _l[j])*(_l[i]-2*_exponent*xA[i]*xA[i]);
}
