#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#ifdef ENABLE_OMP
#include <omp.h>
#endif

#include <Cube/Grid.h>
#include <Cube/GridCP.h>
#include <Common/Constants.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Grid::Grid() :
    _domain(),
    _structure(),
    _values()
{
    reset();
}

Grid::Grid(const Domain& d) :
    _domain(d),
    _structure(),
    _values()
{
    reset();
}

Grid::Grid(std::ifstream& nameFile, const PeriodicTable& table)
{
    readFromCube(nameFile, table);    
}


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

Domain Grid::get_domain() const
{
    return _domain;
}

Structure Grid::get_structure() const
{
    return _structure;
}

std::vector<std::vector<std::vector<std::vector<double>>>> Grid::get_values() const
{
    return _values;
}


void Grid::reset()
{
    if (_domain.get_Nval() < 1 || _domain.get_N1() < 1 || _domain.get_N2() < 1 || _domain.get_N3() < 1)
    {
        _values = std::vector<std::vector<std::vector<std::vector<double>>>>();
    }
    else 
    {
        std::vector<double> values(_domain.get_Nval(), 0);
        std::vector<std::vector<double>> z(_domain.get_N3(), values);
        std::vector<std::vector<std::vector<double>>> y(_domain.get_N2(), z);
        _values = std::vector<std::vector<std::vector<std::vector<double>>>>(_domain.get_N1(), y);
    }
}

void Grid::readFromCube(std::ifstream& nameFile, const PeriodicTable& table)
{
    std::string line;
    getline(nameFile, line);
    getline(nameFile, line);

    int nbAtoms;
    nameFile >> nbAtoms;

    _domain = Domain(nameFile);
    _structure = Structure(nameFile, nbAtoms, table);

    reset();

    for(int i = 0; i < _domain.get_N1(); ++i)
    {
        for(int j = 0; j < _domain.get_N2(); ++j)
        {
            for(int k = 0; k < _domain.get_N3(); ++k)
            {
                for(int l = 0; l < _domain.get_Nval(); ++l)
                {
                    nameFile >> _values[i][j][k][l];
                }
            }
        }
    }

    std::cout << "File has been read... Proceeding" << std::endl;
}



void Grid::set_domain(const Domain& d)
{
    _domain = d;
    reset();
    
}

void Grid::set_structure(const Structure& S)
{
    _structure = S;
}

void Grid::set_values(const std::vector<std::vector<std::vector<std::vector<double>>>>& U)
{
    _values = U;
}
void Grid::set_Vijkl(double rho, int i, int j, int k, int l)
{
    _values[i][j][k][l] = rho;
}

Grid Grid::operator+(const Grid& g)
{
    Grid sum(_domain);

    try
    {
        if(g._domain == _domain)
        {
            sum.set_structure(_structure + g._structure);

            #ifdef ENABLE_OMP
            #pragma omp parallel for shared(sum,g,_values)
            #endif
            #ifdef ENABLE_ACC
            #pragma acc kernels loop
            #endif
            for(int i = 0; i < g._domain.get_N1(); ++i)
            {
                for(int j = 0; j < g._domain.get_N2(); ++j)
                {
                    for(int k = 0; k < g._domain.get_N3(); ++k)
                    {
                        for(int l = 0; l < g._domain.get_Nval(); ++l)
                        {
                            sum._values[i][j][k][l] = _values[i][j][k][l] + g._values[i][j][k][l];
                        }
                    }
                }
            }
        }
        else if(g.get_domain().get_Nval() != get_domain().get_Nval())
        {
            throw std::string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for addition");
        }
        else if (g.get_domain().get_N1() != get_domain().get_N1())
        {
            throw std::string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for addition");
        }
        else if (g.get_domain().get_N2() != get_domain().get_N2())
        {
            throw std::string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for addition");
        }
        else if(g.get_domain().get_N3() != get_domain().get_N3())
        {
            throw std::string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for addition");
        }
        else
        {
            throw std::string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for addition.");
        }
    }
    catch (std::string Error)
    {
        std::cerr << Error << std::endl;

        exit(1);
    }

    return sum;
}

Grid Grid::add(const Grid& g)
{
    try
    {
        if(g._domain == _domain)
        {
            set_domain(_domain);
            _structure.add(g._structure);
            _values.resize(g._domain.get_N1());
            
            #ifdef ENABLE_OMP
            #pragma omp parallel for
            #endif
            #ifdef ENABLE_ACC
            #pragma acc kernels loop
            #endif
            for(int i = 0; i < g._domain.get_N1(); ++i)
            {    
                _values[i].resize(g._domain.get_N2());
                for(int j = 0; j < g._domain.get_N2(); ++j)
                {
                    _values[i][j].resize(g._domain.get_N3());
                    for(int k = 0; k < g._domain.get_N3(); ++k)
                    {
                        _values[i][j][k].resize(g._domain.get_Nval());
                        for(int l = 0; l < g._domain.get_Nval(); ++l)
                        {
                            _values[i][j][k][l] = _values[i][j][k][l] + g._values[i][j][k][l];
                        }
                    }
                }
            }
            return *this;
        }
        else if(g.get_domain().get_Nval() != get_domain().get_Nval())
        {
            throw std::string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for addition");
        }
        else if (g.get_domain().get_N1() != get_domain().get_N1())
        {
            throw std::string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for addition");
        }
        else if (g.get_domain().get_N2() != get_domain().get_N2())
        {
            throw std::string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for addition");
        }
        else if(g.get_domain().get_N3() != get_domain().get_N3())
        {
            throw std::string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for addition");
        }
        else
        {
            throw std::string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for addition.");
        }
    }
    catch (std::string Error)
    {
        std::cerr << Error << std::endl;

        exit(1);
    }
}

Grid Grid::operator*(const Grid& g)
{
    Grid product(_domain);

    try
    {
        if(g._domain == _domain)
        {
            product.set_structure(_structure);

            #ifdef ENABLE_OMP
            #pragma omp parallel for
            #endif
            #ifdef ENABLE_ACC
            #pragma acc kernels loop
            #endif
            for(int i = 0; i < g._domain.get_N1(); ++i)
            {
                for(int j = 0; j < g._domain.get_N2(); ++j)
                {
                    for(int k = 0; k < g._domain.get_N3(); ++k)
                    {
                        for(int l = 0; l < g._domain.get_Nval(); ++l)
                        {
                            product._values[i][j][k][l] = _values[i][j][k][l] * g._values[i][j][k][l];
                        }
                    }
                }
            }
        }
        else if(g.get_domain().get_Nval() != get_domain().get_Nval())
        {
            throw std::string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for product");
        }
        else if (g.get_domain().get_N1() != get_domain().get_N1())
        {
            throw std::string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for product");
        }
        else if (g.get_domain().get_N2() != get_domain().get_N2())
        {
            throw std::string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for product");
        }
        else if(g.get_domain().get_N3() != get_domain().get_N3())
        {
            throw std::string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for product");
        }
        else
        {
            throw std::string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for product");
        }
    }
    catch (std::string Error)
    {
        std::cerr << Error << std::endl;

        exit(1);
    }

    return product;
}

Grid Grid::operator-(const Grid& g)
{
    Grid diff(_domain);

    try
    {
        if(g._domain == _domain)
        {
            diff.set_structure(_structure);

            #ifdef ENABLE_OMP
            #pragma omp parallel for
            #endif
            #ifdef ENABLE_ACC
            #pragma acc kernels loop
            #endif
            for(int i = 0; i < g._domain.get_N1(); ++i)
            {
                for(int j = 0; j < g._domain.get_N2(); ++j)
                {
                    for(int k = 0; k < g._domain.get_N3(); ++k)
                    {
                        for(int l = 0; l < g._domain.get_Nval(); ++l)
                        {
                            diff._values[i][j][k][l] = _values[i][j][k][l] - g._values[i][j][k][l];
                        }
                    }
                }
            }
        }
        else if(g.get_domain().get_Nval() != get_domain().get_Nval())
        {
            throw std::string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for difference");
        }
        else if (g.get_domain().get_N1() != get_domain().get_N1())
        {
            throw std::string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for difference");
        }
        else if (g.get_domain().get_N2() != get_domain().get_N2())
        {
            throw std::string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for difference");
        }
        else if(g.get_domain().get_N3() != get_domain().get_N3())
        {
            throw std::string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for difference");
        }
        else
        {
            throw std::string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for difference");
        }
    }
    catch (std::string Error)
    {
        std::cerr << Error << std::endl;

        exit(1);
    }

    return diff;
}

double Grid::sum()
{
    double sum = 0;

    #ifdef ENABLE_OMP
    #pragma omp parallel for reduction (+:sum)
    #endif
    #ifdef ENABLE_ACC
    #pragma acc kernels loop reduction(+:sum)
    #endif
    for(int i = 0; i < _domain.get_N1(); ++i)
    {    
        for(int j = 0; j < _domain.get_N2(); ++j)
        {    
            for(int k = 0; k < _domain.get_N3(); ++k)
            {
                for(int l = 0; l < _domain.get_Nval(); ++l)
                {
                    sum += _values[i][j][k][l];
                }
            }
        }
    }

    return sum;    
}

double Grid::integrateOverDomain()
{
    return sum() * _domain.get_dv();
}

void Grid::reset_Boundary(int nBound)
{
    int NV = (_domain.get_Nval() > 1) ? 1 : 0;

    // Left
    for(int i = 0; i < nBound; ++i)
    {
        for(int j = 0; j < _domain.get_N2() ; ++j)
        {
            for(int k = 0; k < _domain.get_N3(); ++k)
            {    
                for(int l = NV; l < _domain.get_Nval(); ++l)
                {
                    _values[i][j][k][l] = _values[nBound][j][k][l];
                }
            }
        }
    }

    std::cout << "done left" << std::endl;


    // Right
    for(int i = _domain.get_N1() - nBound; i < _domain.get_N1(); ++i)
    {    
        for(int j = 0; j < _domain.get_N2(); ++j)
        {
            for(int k = 0; k < _domain.get_N3(); ++k)
            {    
                for(int l = NV; l < _domain.get_Nval(); ++l)
                {
                    _values[i][j][k][l] = _values[_domain.get_N1() - nBound - 1][j][k][l];
                }
            }
        }
    }
    
    std::cout << "done right" << std::endl;


    // Front
    for(int j = 0; j < nBound; ++j)
    {
        for(int i = 0; i < _domain.get_N1(); ++i)
        {    
            for(int k = 0; k < _domain.get_N3(); ++k)
            {    
                for(int l = NV; l < _domain.get_Nval(); ++l)
                {
                    _values[i][j][k][l] = _values[i][nBound][k][l];
                }
            }
        }
    }

    std::cout << "done front" << std::endl;
    
    
    // Back
    for(int j = _domain.get_N2() - nBound; j < _domain.get_N2(); ++j)
    {
        for(int i = 0; i < _domain.get_N1(); ++i)
        {    
            for(int k = 0; k < _domain.get_N3(); ++k)
            {
                for(int l = NV; l < _domain.get_Nval(); ++l)
                {
                    _values[i][j][k][l] = _values[i][_domain.get_N2() - nBound - 1][k][l];
                }
            }
        }
    }

    std::cout << "done back" << std::endl;
    
    
    // Top
    for(int k = 0; k < nBound; ++k)
    {    
        for(int j = 0; j < _domain.get_N2(); ++j)
        {
            for(int i = 0; i < _domain.get_N1(); ++i)
            {
                for(int l = NV; l < _domain.get_Nval(); ++l)
                {    
                    _values[i][j][k][l] = _values[i][j][nBound][l];
                }
            }
        }
    }

    std::cout << "done top" << std::endl;

    
    // Bottom
    for(int k = _domain.get_N3() - nBound; k < _domain.get_N3(); ++k)
    {    
        for(int j = 0; j < _domain.get_N2(); ++j)
        {
            for(int i = 0; i < _domain.get_N1(); ++i)
            {
                for(int l = NV; l < _domain.get_Nval(); ++l)
                {
                    _values[i][j][k][l] = _values[i][j][_domain.get_N3() - nBound - 1][l];
                }
            }
        }
    }

    std::cout << "done bottom" << std::endl;
}


Grid Grid::coulomb_Grid(double q, std::vector<double> R)
{
    Grid coulombGrid(_domain);

    double x = 0;
    double y = 0;
    double z = 0;
    double value = 0;

    #ifdef ENABLE_OMP
    #pragma omp parallel for
    #endif
    #ifdef ENABLE_ACC
    #pragma acc kernels loop
    #endif
    for(int i = 0; i < _domain.get_N1(); ++i)
    {    
        for(int j = 0; j < _domain.get_N2(); ++j)
        {        
            for(int k = 0; k < _domain.get_N3(); ++k)
            {    
                for(int l = 0; l < _domain.get_Nval(); ++l)
                {
                    x = _domain.get_origin()[0] + i * _domain.get_T()[0][0] + j * _domain.get_T()[0][1] +  k * _domain.get_T()[0][2]; 
                    y = _domain.get_origin()[1] + i * _domain.get_T()[1][0] + j * _domain.get_T()[1][1] +  k * _domain.get_T()[1][2]; 
                    z = _domain.get_origin()[2] + i * _domain.get_T()[2][0] + j * _domain.get_T()[2][1] +  k * _domain.get_T()[2][2];
                    
                    value = q / sqrt( (x - R[0]) * (x - R[0])
                                      + (y - R[1]) * (y - R[1])
                                      + (z - R[2]) * (z - R[2]) );
                    
                    value = value < 1e-10 ? 0 : value;

                    coulombGrid._values[i][j][k][l] = value;
                }
            }
        }
    }

    return coulombGrid;
}

void Grid::coefs_Laplacian(int nBound, std::vector<double>& fcx, std::vector<double>& fcy, std::vector<double>& fcz, double& cc) const
{
    std::vector<double> coefs(nBound+1);

    if(nBound == 1)
    {
        std::vector<double> c = {-2.0, 1.0};
        for(int i = 0; i <= nBound; ++i)
        {
            coefs[i] = c[i];
        }
    }
    else if(nBound == 2)
    {
        double denom = 12.0;

        std::vector<double> c = {-30.0, 16.0, -1.0};
        for(int i = 0; i <= 2; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 3)
    {
        double denom = 180.0;

        std::vector<double> c = {-490.0, 270.0,-27.0, 2.0};
        for(int i = 0; i <= 3; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 4)
    {
        double denom = 5040.0;

        std::vector<double> c = {-14350.0, 8064.0, -1008.0, 128.0, -9.0};
        for(int i = 0; i <= 4; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 5)
    {
        double denom = 25200.0;

        std::vector<double> c = {-73766.0, 42000.0, -6000.0, 1000.0, -125.0, 8.0};
        for(int i = 0; i <= 5; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 6)
    {
        double denom = 831600.0;

        std::vector<double> c = {-2480478.0,1425600.0,-222750.0,44000.0,-7425.0,864.0,-50.0};
        for(int i = 0; i <= 6; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 7)
    {
        double denom = 75675600.0;

        std::vector<double> c = {-228812298.0,132432300.0,-22072050.0,4904900.0,-1003275.0, 160524.0,-17150.0,900.0};
        for(int i = 0; i <= 7; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 8)
    {
        double denom = 302702400.0;

        std::vector<double> c = {-924708642.0,538137600.0,-94174080.0,22830080.0,-5350800.0,1053696.0,-156800.0,15360.0,-735.0};
        for(int i = 0; i <= 8; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 9)
    {
        double denom = 15437822400.0;

        std::vector<double> c = {-47541321542.0,+27788080320.0, -5052378240.0,+1309875840.0,-340063920.0,+77728896.0,-14394240.0,+1982880.0,-178605.0,+7840.0};
        for(int i = 0; i <= 9; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 10)
    {
        double denom = 293318625600.0;

        std::vector<double> c = {-909151481810.0,+533306592000.0, -99994986000.0,+27349056000.0,-7691922000.0,+1969132032.0,-427329000.0,+73872000.0,-9426375.0,+784000.0,-31752.0};
        for(int i = 0; i <= 10; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 11)
    {
        double denom = 3226504881600.0;

        std::vector<double> c = {-10053996959110.0,+5915258949600.0,-1137549798000.0,+325014228000.0,-97504268400.0,+27301195152.0,-6691469400.0,+1365606000.0,-220114125.0,+26087600.0,-2012472.0,+75600.0};
        for(int i = 0; i <= 11; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else if (nBound == 12)
    {
        double denom = 74209612276800.0;

        std::vector<double> c = {-232272619118930.0,+137002361126400.0,-26911178078400.0,+7973682393600.0,-2522922944850.0,+759845028096.0,-205205061600.0,+47609337600.0,-9112724775.0,+1371462400.0,-151484256.0,+10886400.0,-381150.0};
        for(int i = 0; i <= 12; ++i)
        {
            coefs[i] = c[i] / denom;
        }
    }
    else
    {
        exit(1);
    }

    cc = 1 / (_domain.get_dx() *_domain.get_dx())
         + 1 / (_domain.get_dy()*_domain.get_dy())
         + 1 / (_domain.get_dz()*_domain.get_dz());

    cc *= coefs[0];

    for(int i = 0; i <= nBound; ++i)
    {
        fcx[i] = coefs[i] / (_domain.get_dx() * _domain.get_dx());
        fcy[i] = coefs[i] / (_domain.get_dy() * _domain.get_dy());
        fcz[i] = coefs[i] / (_domain.get_dz() * _domain.get_dz());
    }
}

Grid Grid::laplacian(int nBound) const
{
    try
    {
        if(nBound<1 or nBound>12)
        {
            throw std::string("nBound oustide acceptable precision values");
        }
        else
        {
            std::vector<double> fcx(nBound);
            std::vector<double> fcy(nBound); 
            std::vector<double> fcz(nBound);
            double cc=0;
            Grid g(_domain);
            g.set_structure(_structure);
            coefs_Laplacian(nBound, fcx, fcy, fcz, cc);
#ifdef ENABLE_OMP
//std::cout<<"Number of cores = "<<omp_get_num_procs()<<std::endl;
#pragma omp parallel for shared(g,_values,fcx,fcy,fcz,cc,nBound)
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
            //for(int i=nBound;i<g._dom.N1()-nBound;i++) // DEBUG
            for(int ii=0;ii<g._domain.get_N1()-2*nBound;ii++)
            {
                int i=ii+nBound;
                for(int j=nBound;j<g._domain.get_N2()-nBound;j++)
                {    
                    for(int k=nBound;k<g._domain.get_N3()-nBound;k++)
                    {
                        double v = cc*_values[i][j][k][0];
                        for(int n=1; n<=nBound ; n++)
                        {
                            v += fcx[n] *(_values[i-n][j][k][0]+_values[i+n][j][k][0]);
                            v += fcy[n] *(_values[i][j-n][k][0]+_values[i][j+n][k][0]);
                            v += fcz[n] *(_values[i][j][k-n][0]+_values[i][j][k+n][0]);
                        }
                        g._values[i][j][k][0]=v;
                    }
                }
            }
            g.reset_Boundary(nBound);
            return g;
        }
    }
    catch(std::string error)
    {
        std::cout<<error<<std::endl;
        exit(1);
    }                
}

void Grid::coefs_Gradient(int nBound, std::vector<double>& fcx, std::vector<double>& fcy, std::vector<double>& fcz) const
{
    std::vector<double> coefs(nBound+1);    
        if(nBound==1)
        {
            double denom = 2.0;
            std::vector<double> c = {-1.0};
            for(int i=0;i<nBound;i++)
                    coefs[i] = c[i]/denom;
        }
        else if(nBound==2)
        {
            double denom = 12.0;
            std::vector<double> c = { 1.0, -8.0};
            for(int i=0;i<2;i++)
                coefs[i] = c[i]/denom;
        }
        else if (nBound==3)
        {
            double denom = 60.0;
            std::vector<double> c ={ -1.0, +9.0, -45.0};
            for(int i=0;i<3;i++)
                coefs[i] = c[i]/denom;
        }
        else if (nBound==4)
        {
            double denom = 840.0;
            std::vector<double> c = { 3.0, -32.0, +168.0, -672.0};
            for(int i=0;i<4;i++)
                coefs[i] = c[i]/denom;
        }
        else if (nBound==5)
        {
            double denom = 2520.0;
            std::vector<double> c = { -2.0, +25.0, -150.0,+600.0, -2100.0};
            for(int i=0;i<5;i++)
                coefs[i] = c[i]/denom;
        }
        else if (nBound==6)
        {
            double denom = 27720.0;
            std::vector<double> c = { 5.0, -72.0, +495.0, -2200.0, +7425.0, -23760.0};
            for(int i=0;i<6;i++)
                coefs[i] = c[i]/denom;
        }
        else if (nBound==7)
        {
            double denom = 360360.0;
            std::vector<double> c = { -15.0, +245.0, -1911.0, +9555.0, -35035.0, +105105.0, -315315.0};
            for(int i=0;i<7;i++)
                coefs[i] = c[i]/denom;
        }
        else if (nBound==8)
        {
            double denom = 720720.0;
            std::vector<double> c = { 7.0, -128.0, +1120.0, -6272.0, +25480.0, -81536.0, +224224.0, -640640.0};
            for(int i=0;i<8;i++)
                coefs[i] = c[i]/denom;
        }
        else
        {
            exit(1);
        }

        for(int i=0;i<nBound;i++)
        {
            fcx[i] =  coefs[i]/_domain.get_dx();
            fcy[i] =  coefs[i]/_domain.get_dy();
            fcz[i] =  coefs[i]/_domain.get_dz();
        }
}


Grid Grid::gradient(int nBound) const
{
    try
    {
        if(nBound<1 or nBound>8)
        {
            throw std::string("nBound oustide bounds");
        }
        else
        {
            std::vector<double> fcx(nBound);
            std::vector<double> fcy(nBound); 
            std::vector<double> fcz(nBound);
            Domain dg=_domain;
            dg.set_Nval(4);
            Grid g(dg);
            g._structure=_structure;
            
            coefs_Gradient(nBound, fcx, fcy, fcz);
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
            for(int i=nBound;i<g._domain.get_N1()-nBound;i++) // DEBUG
            //for(int ii=0;ii<g._dom.N1()-2*nBound;ii++)
            {
                //int i=ii+nBound;
                for(int j=nBound;j<g._domain.get_N2()-nBound;j++)
                {    
                    for(int k=nBound;k<g._domain.get_N3()-nBound;k++)
                    {
                        g._values[i][j][k][0]=_values[i][j][k][0];
                        double gx=0;
                        double gy=0;
                        double gz=0;
                        for(int n=-nBound, kn=0 ; kn<nBound ; n++, kn++)
                        {
                            gx += fcx[kn] * (_values[i+n][j][k][0]-_values[i-n][j][k][0]);
                            gy += fcy[kn] * (_values[i][j+n][k][0]-_values[i][j-n][k][0]);
                            gz += fcz[kn] * (_values[i][j][k+n][0]-_values[i][j][k-n][0]);
                        }
                        g._values[i][j][k][1]=gx;
                        g._values[i][j][k][2]=gy;
                        g._values[i][j][k][3]=gz;
                    }
                }
            }
            g.reset_Boundary(nBound);
            return g;
        }
    }
    catch(std::string error)
    {
        std::cout<<error<<std::endl;
        exit(1);
    }                
}



/*
void bicub_Coef(std::vector<double> z, std::vector<double> dzdx, std::vector<double> dzdy, std::vector<double> d2zdxdy, double dx, double dy, std::vector<std::vector<double>> c)
{
    std::vector<std::vector<int>> Ainv[16][16]={
        {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},{-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0}{2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},{0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},{0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},{-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},{9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},{-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},{2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},{-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},{4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
    int i,j,k,l;
    double s, dxdy;
    std::vector<double> c1D(16),x(16);
    dxdy= dx*dy;
    for (i=0;i<4;i++)
    {
        x[i]=z[i];
        x[i+4]=dzdx[i]*dx;
        x[i+8]=dzdy[i]*dy;
        x[i+12]=d2zdxdy[i]*dxdy;
    }
    for (i=0;i<16;i++)
    {
        s=0.0;
        for (k=0;k<16;k++) 
        {
            s += Ainv[i][k]*x[k];
        }
        c1D[i]=s;
    }
    l=0;
    for (i=0;i<4;i++)
    {
        for (j=0;j<4;j++)
        {
        c[i][j]=c1D[l++];
        }
    }
}
*/



Grid Grid::finer_Grid()
{
    Grid g(Domain(_domain.get_Nval() ,2*_domain.get_N1()-1,2*_domain.get_N2()-1,2*_domain.get_N3()-1, _domain.get_origin()));
    g._structure=_structure;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            g._domain.set_Tij(_domain.get_Tij(i,j)/2, i, j);
        }
    }
    //coarse grid points to fine grid

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=0;i<_domain.get_N1();i++)
    {    
        for(int j=0;j<_domain.get_N2();j++)
        {    
            for(int k=0;k<_domain.get_N3();k++)
            {
                g._values[2*i][2*j][2*k]=_values[i][j][k];
            }    
        }
    }
        
    //center cube points
#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=1;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=1;j<g._domain.get_N2()-1;j+=2)
        {    
            for(int k=1;k<g._domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval();l++)
                {    
                    g._values[i][j][k][l] = 0.125*g._values[i-1][j-1][k-1][l] + 0.125*g._values[i-1][j-1][k+1][l] + 0.125*g._values[i-1][j+1][k-1][l] + 0.125*g._values[i-1][j+1][k+1][l] + 0.125*g._values[i+1][j-1][k-1][l] + 0.125*g._values[i+1][j-1][k+1][l] + 0.125*g._values[i+1][j+1][k-1][l] + 0.125*g._values[i+1][j+1][k+1][l];
                }    
            }
        }
    }
    
                
    //cubes face points
#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=0;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=0;j<g._domain.get_N2()-1;j+=2)
        {    
            for(int k=1;k<g._domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval();l++)
                {
                    g._values[i][j][k][l]= 0.5*g._values[i][j][k-1][l]+0.5*g._values[i][j][k+1][l];
                }
            }    
        }
    }

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=0;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=1;j<_domain.get_N2()-1;j+=2)
        {    
            for(int k=0;k<_domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval();l++)
                {
                    g._values[i][j][k][l] = 0.5*g._values[i][j-1][k][l]+0.5*g._values[i][j+1][k][l];
                }
            }    
        }
    }

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=1;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=0;j<g._domain.get_N2()-1;j+=2)
        {    
            for(int k=0;k<g._domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval();l++)
                {
                    g._values[i][j][k][l] = 0.5*g._values[i-1][j][k][l]+0.5*g._values[i+1][j][k][l];
                }
            }    
        }
    }

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=0;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=1;j<g._domain.get_N2()-1;j+=2)
        {    
            for(int k=1;k<g._domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval()-1;l++)
                {
                    g._values[i][j][k][l] = 0.25*g._values[i][j-1][k-1][l] + 0.25*g._values[i+1][j-1][k+1][l] + 0.25*g._values[i][j+1][k-1][l] + 0.25*g._values[i][j+1][k+1][l];
                }
            }    
        }
    }
    
#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=1;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=0;j<g._domain.get_N2()-1;j+=2)
        {    
            for(int k=1;k<g._domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval();l++)
                {
                    g._values[i][j][k][l] = 0.25*g._values[i-1][j][k-1][l] + 0.25*g._values[i-1][j][k+1][l] + 0.25*g._values[i+1][j][k-1][l] + 0.25*g._values[i+1][j][k+1][l];
                }
            }    
        }
    }

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i=1;i<g._domain.get_N1()-1;i+=2)
    {    
        for(int j=1;j<g._domain.get_N2()-1;j+=2)
        {    
            for(int k=0;k<g._domain.get_N3()-1;k+=2)
            {
                for(int l=0; l<g._domain.get_Nval();l++)
                {
                    g._values[i][j][k][l] = 0.25*g._values[i-1][j-1][k][l] + 0.25*g._values[i+1][j-1][k][l] + 0.25*g._values[i-1][j+1][k][l] + 0.25*g._values[i+1][j+1][k][l];
                }
            }    
        }
    }
    return g;
}

Grid Grid::coarser_Grid()
{
    int N[3]={_domain.get_N1(),_domain.get_N2(), _domain.get_N3()};
    for(int i=0; i<3; i++)
    {
        if(N[i]%2==1)
        {
            N[i]=N[i]-1;
        }
        N[i] =N[i]/2;
    }
    Grid g(Domain(_domain.get_Nval(), N[0], N[1], N[2], _domain.get_origin()));
    g._structure=_structure;
    double scale = 1.0 / 64.0;
    

    /*printf("Begin restriction\n");*/

    int iXBegin = 1;
    int iXEnd = g._domain.get_N1();
    int iYBegin = 1;
    int iYEnd = g._domain.get_N2();
    int iZBegin = 1;
    int iZEnd = g._domain.get_N3();

#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
    for(int i = iXBegin ; i <= iXEnd-1 ; i++)
    {
        int x0, xp, xm;
        x0 = 2 * i;
        xp = x0 + 1;
        xm = x0 - 1;
        for(int j = iYBegin ; j <= iYEnd-1 ; j++)
        {
            int y0, yp, ym;
            y0 = 2 * j;
            yp = y0 + 1;
            ym = y0 - 1;
            for(int k= iZBegin ; k <= iZEnd-1 ; k++)
            {
                int z0, zp, zm;
                z0 = 2 * k;
                zp = z0 + 1;
                zm = z0 - 1;
                for(int l=0;l<_domain.get_Nval(); l++)
                {     
                    double face, corner, edge;
                    
                    face = _values [xm][y0][z0][l] +_values [xp][y0][z0][l] +_values [x0][ym][z0][l] +_values [x0][yp][z0][l] + _values [x0][y0][zm][l] +_values [x0][y0][zp][l];

                    corner =  _values [xm] [ym] [zm][l] +_values [xm] [ym] [zp][l] +_values [xm] [yp] [zm][l] +_values [xm] [yp] [zp][l]+_values [xp] [ym] [zm][l] +_values [xp] [ym] [zp][l] +_values [xp] [yp] [zm][l] +_values [xp] [yp] [zp][l];
                    
                    edge = _values [xm] [y0] [zm][l] +_values [xm] [ym] [z0][l] +_values [xm] [yp] [z0][l] +_values [xm] [y0] [zp][l] +_values [x0] [ym] [zm][l] +_values [x0] [yp] [zm][l] +_values [x0] [ym] [zp][l] +_values [x0] [yp] [zp][l] +_values [xp] [y0] [zm][l] +_values [xp] [ym] [z0][l] +_values [xp] [yp] [z0][l] +_values [xp] [y0] [zp][l];
    
                    g._values [i][j][k][l]=scale * (8.0 * _values [x0][y0][z0][l] +4.0 * face +2.0 * edge +corner);
                      }
            }
        }
    }
    return g;
}

void Grid::save(ofstream& nameFile)
{
    nameFile.precision(14);
    if(_domain.get_Nval()==1)
    {
        nameFile<<"Grid generated by CdftT "<<std::endl;
        nameFile<<"Density"<<std::endl;
    }
    else if(_domain.get_Nval()==4)
    {
        nameFile<<"Grid generated by CdftT "<<std::endl;
        nameFile<<"Gradient"<<std::endl;
    }
    else
    {
        nameFile<<"Grid generated by CdftT "<<std::endl;
        nameFile<<"Orbitals"<<std::endl;
    }
    nameFile<<scientific;
    nameFile<<_structure.number_of_atoms()<<" ";
    for(int i=0;i<3;i++)
    {
        nameFile<<_domain.get_origin()[i]<<" ";
    }
    nameFile<<_domain.get_Nval()<<" ";
    nameFile<<std::endl<<_domain.get_N1()<<" ";
    if(_domain.get_N1()>0)
    {
        for(int i=0;i<3;i++)
        {
            nameFile<<_domain.get_Tij(0,i)<<" ";
        }
        nameFile<<std::endl;
        nameFile<<_domain.get_N2()<<" ";
        for(int i=0;i<3;i++)
        {
            nameFile<<_domain.get_Tij(1,i)<<" ";
        }
        nameFile<<std::endl<<_domain.get_N3()<<" ";
        for(int i=0;i<3;i++)
        {
            nameFile<<_domain.get_Tij(2,i)<<" ";
        }
    }
    else
    {
        for(int i=0;i<3;i++)
        {
            nameFile << _domain.get_Tij(0, i) * Constants::BOHR_RADIUS_TO_ANGSTROM << " ";
        }
        nameFile<<std::endl;
        nameFile<<_domain.get_N2()<<" ";
        for(int i=0;i<3;i++)
        {
            nameFile << _domain.get_Tij(1, i) * Constants::BOHR_RADIUS_TO_ANGSTROM << " ";
        }
        nameFile<<std::endl<<_domain.get_N3()<<" ";
        for(int i=0;i<3;i++)
        {
            nameFile << _domain.get_Tij(2, i) * Constants::BOHR_RADIUS_TO_ANGSTROM << " ";
        }
    }
    nameFile<<std::endl;
    for(int i=0;i<_structure.number_of_atoms();i++)
    {
        nameFile<<_structure.atom(i).get_atomicNumber()<<" ";
        if(_structure.atom(i).get_charge()<1)
        {
            nameFile<<double(_structure.atom(i).get_atomicNumber())<<" ";
        }
        else
        {
            nameFile<<_structure.atom(i).get_charge()<<" ";
        }
        for(int j=0; j<3;j++)
        {
            nameFile<<_structure.atom(i).get_coordinates()[j]<<" ";
        }
        nameFile<<std::endl;
    }
    for(int i=0; i<_domain.get_N1();i++)
    {    
        
        for(int j=0; j<_domain.get_N2();j++)
        {    
            int R=0;
            for(int k=0; k<_domain.get_N3();k++)
            {
                for(int l=0; l<_domain.get_Nval();l++)
                {
                    nameFile<<_values[i][j][k][l]<<" ";
                    R++;
                    if(R%6==0)
                    {
                        nameFile<<std::endl;
                    }
                }
            }
            if(R%6!=0)
            {
                nameFile<<std::endl;
            }
        }
    }
}
double Grid::value(int i, int j, int k) const
{
    return _values[i][j][k][0];
}
double Grid::value(int i, int j, int k, int l) const
{
    return _values[i][j][k][l];
}

void Grid::next(int i, int j, int k, double& current, std::vector<std::vector<int>>& trajectory)
{
    std::vector<int> v=trajectory.back();
    int I[3];
    int i2 =i-1;
    int i1 =i+1;
    I[0] = i2;
    I[1] = i;
    I[2] = i1;
        
    int J[3];
    int j1 = j+1;
    int j2 = j-1;
    J[0] = j2;
    J[1] = j;
    J[2] = j1;
            
    int K[3];
    int k1 = k+1;
    int k2 = k-1;
    K[0] = k2;
    K[1] = k;
    K[2] = k1;
    if(i2<0) i2 = i;
    if(i1>_domain.get_N1()-1) i1 = i;
    if(j2<0) j2 = j;
    if(j1>_domain.get_N2()-1) j1 = j;
    if(k2<0) k2 = k;
    if(k1>_domain.get_N3()-1) k1 = k;
    for(int ic=0;ic<3;ic++)
    {
        for(int jc=0;jc<3;jc++)
        {    
                
            for(int kc=0;kc<3;kc++)
            {
                if(ic==1 and jc==1 and kc==1) continue;
                std::vector<double> ds(3);
                ds[0]=(I[ic]-I[1])*_domain.get_dx();
                ds[1]=(J[jc]-J[1])*_domain.get_dy();
                ds[2]=(K[kc]-K[1])*_domain.get_dz();
                double grad =0;
                double normds=sqrt(ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2]);
                for(int m=1;m<=3;m++)
                {
                    grad +=_values[I[ic]][J[jc]][K[kc]][m]*ds[m-1];
                }
                grad=grad/normds;
                if(grad>current)
                {
                    v[0]=I[ic];
                    v[1]=J[jc];
                    v[2]=K[kc];
                    current = grad;
                }
            }
        }
    }
    
    trajectory.push_back(v);
}

std::vector<double> Grid::atom_attract_diff(const std::vector<std::vector<int>>& attract)
{
    std::vector<double> v(_structure.number_of_atoms());
    double distance=0;
    double d1=100;
    std::vector<double> ds(3);
    ds[0]=_domain.get_dx();
    ds[1]=_domain.get_dy();
    ds[2]=_domain.get_dz();
    for(int j=0; j<_structure.number_of_atoms();j++)
    {
        for(int n=0;n<int(attract.size()); n++)
        {
            for(int i=0;i<3;i++)
            {
                distance += double((_structure.get_atoms()[j].get_coordinates()[i]-attract[n][i])*(_structure.get_atoms()[j].get_coordinates()[i]-attract[n][i]))*ds[i]*ds[i];
            }
            distance=sqrt(distance);
            if(distance<d1)
            {
                d1=distance;
            }    
        }
        v[j]=d1;
        
    }
    return v;
}


void Grid::addSurroundingEqualPoints(int i,int j,int k, std::vector<std::vector<int>>& equals, double& current)
{
    std::vector<int> v(3);
    int I[3];
    int i2 =i-1;
    int i1 =i+1;
    I[0] = i2;
    I[1] = i;
    I[2] = i1;
    
    int J[3];
    int j1 = j+1;
    int j2 = j-1;
    J[0] = j2;
    J[1] = j;
    J[2] = j1;
    
    int K[3];
    int k1 = k+1;
    int k2 = k-1;
    K[0] = k2;
    K[1] = k;
    K[2] = k1;
    for(int ic=0;ic<3;ic++)
    {    
        for(int jc=0;jc<3;jc++)
        {
            for(int kc=0;kc<3;kc++)
            {
                if(ic==1 and jc==1 and kc==1) continue;
                if(i2<0) i2 = i;
                if(i1>_domain.get_N1()-1) i1 = i;
                if(j2<0) j2 = j;
                if(j1>_domain.get_N2()-1) j1 = j;
                if(k2<0) k2 = k;
                if(k1>_domain.get_N3()-1) k1 = k;
                std::vector<double> ds(3);
                ds[0]=(I[ic]-I[1])*_domain.get_dx();
                ds[1]=(J[jc]-J[1])*_domain.get_dy();
                ds[2]=(K[kc]-K[1])*_domain.get_dz();
                double grad=0;
                double normds=sqrt((ds[0]*ds[0])+(ds[1]*ds[1])+(ds[2]*ds[2]));
                for(int m=1;m<3;m++)
                {
                    grad +=_values[I[1]][J[1]][K[1]][m]*ds[m-1];
                }
                grad=grad/normds;
                if(grad==current)    
                {
                    v[0]=I[ic];
                    v[1]=J[jc];
                    v[2]=K[kc];
                    equals.push_back(v);
                }    
            }
        }
    }
}
Grid Grid::aim_On_Grid(int nBound)
{
    try
    {
        if(_domain.get_Nval()==4)
        {
            Grid g(Domain(1,_domain.get_N1(), _domain.get_N2(), _domain.get_N3(), _domain.get_origin()));
            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    g._domain.set_Tij(_domain.get_Tij(i,j), i, j);
                }
            }
            g._structure=_structure;
            //make list instead
            std::vector<std::vector<int>> attractors(0, std::vector<int>(3));
            std::cout<<_values[36][29][32][2]<<std::endl;
            std::cout<<_values[36][30][32][2]<<std::endl;
            for(int i=nBound;i<_domain.get_N1()-nBound;i++)
            {
                for(int j=nBound;j<_domain.get_N2()-nBound;j++)
                {
                    for(int k=nBound;k<_domain.get_N3()-nBound;k++)
                    {
                        if(g._values[i][j][k][0]>PRECISION)
                        {
                            continue;
                        }
                        else
                        {    
                            double current=0;
                            bool KnownPt=false;
                            std::vector<std::vector<int>> trajectory={{i,j,k}};
                            std::vector<std::vector<int>> equals(0, std::vector<int>(3));
                            do
                            {
                                current=0;
                                int I=trajectory.back()[0];
                                int J=trajectory.back()[1];
                                int K=trajectory.back()[2];
                                next(I,J,K,current,trajectory);
                                
                                
                                std::cout<<trajectory.back()[0]<<std::endl;
                                std::cout<<trajectory.back()[1]<<std::endl;
                                std::cout<<trajectory.back()[2]<<std::endl;
                                if(g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0]>PRECISION)
                                {
                                    KnownPt=true;
                                    break;
                                }
                            } while(current>PRECISION);
                            for(int p=0;p<int(trajectory.size()-1);p++)
                            {
                                g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];
                            }
                            addSurroundingEqualPoints(trajectory.back()[0],trajectory.back()[1],trajectory.back()[2],equals,current);
                            for(int p=0;p<int(equals.size());p++)
                            {
                                g._values[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];    
                            }
                            //if point has been already indexed exit loop: new point
                            if(KnownPt)
                            {
                                std::cout<<"gone"<<std::endl;
                                continue;    
                            }
                            
                            //Do we know max?
                            bool WhoIsMax=true;
                            //compare max to list
                            std::cout<<i<<std::endl;
                            
                            for(int m=0;m<int(attractors.size());m++)
                            {
                                int I=1;
                                
                                if(attractors[m]==trajectory.back())
                                {
                                    WhoIsMax=false;
                                    for(int p=0;p<int(trajectory.size());p++)
                                    {
                                        g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(I);
                                    }
                                }
                                I++;
                            }
                            //dont know max
                            //add to known attractors
                            if(WhoIsMax)
                            {    
                                if(int(trajectory.size())<(_domain.get_N1()*_domain.get_N2()*_domain.get_N3())/1000)
                                {
                                    
                                    for(int p=0;p<int(equals.size());p++)
                                    {
                                        g._values[equals[p][0]][equals[p][1]][equals[p][2]][0]=double(attractors.size());        
                                    }
                                    for(int p=0;p<int(trajectory.size());p++)
                                    {
                                        g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
                                    }
                                }
                                attractors.push_back(trajectory.back());
                                //index all points from trajectory
                                std::cout<<"traj size "<<trajectory.size()<<std::endl;
                                std::cout<<-_values[36][29][32][2]<<std::endl;
                                std::cout<<"attract size "<<attractors.size()<<std::endl;
                                for(int p=0;p<int(equals.size());p++)
                                {
                                    g._values[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];        
                                }
                                for(int p=0;p<int(trajectory.size());p++)
                                {
                                    g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
                                }
                            }
                        
                        }                    
                    }
                }
            }
            /*
            std::vector<std::vector<int>> vpos(0,std::vector<int>(3));
            for(int p=0;p<int(attractors.size()); p++)
            {
                for(int m=0;m<int(attractors.size());m++)
                {
                    int I[2];
                    I[0]=-1;
                    I[1]=1;
                    for(int i=0;i<2;i++)
                    {
                        if(attractors[p][0]==attractors[m][0] +I[i] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2])
                        {
                            vpos.push_back(attractors[m]);
                            for(int l=nBound;l<_dom.N1()-nBound;l++)
                            {
                                for(int j=nBound;j<_dom.N2()-nBound;j++)
                                {
                                    for(int k=nBound;k<_dom.N3()-nBound;k++)
                                    {
                                        if(g._values[l][j][k][0]==m)
                                        {
                                            std::cout<<"okay"<<std::endl;
                                            g._values[l][j][k][0]=p;
                                            
                                        }
                                    }
                                }
                            }
                        }
                        if( attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] +I[i] and attractors[p][1]==attractors[m][1])
                        {
                            vpos.push_back(attractors[m]);
                            for(int l=nBound;l<_dom.N1()-nBound;l++)
                            {
                                for(int j=nBound;j<_dom.N2()-nBound;j++)
                                {
                                    for(int k=nBound;k<_dom.N3()-nBound;k++)
                                    {
                                        if(g._values[l][j][k][0]==m)
                                        {
                                            std::cout<<"okay"<<std::endl;
                                            g._values[l][j][k][0]=p;
                                            
                                        }
                                    }
                                }
                            }
                        }
                        if(attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2] +I[i])
                        {
                            vpos.push_back(attractors[m]);
                            for(int l=nBound;l<_dom.N1()-nBound;l++)
                            {
                                for(int j=nBound;j<_dom.N2()-nBound;j++)
                                {
                                    for(int k=nBound;k<_dom.N3()-nBound;k++)
                                    {
                                        if(g._values[l][j][k][0]==m)
                                        {
                                            std::cout<<"okay"<<std::endl;
                                            g._values[l][j][k][0]=p;
                                            
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for(int k=0;k<int(vpos.size());k++)
            {
                for(auto element:attractors)
                {    
                    auto it=find(attractors.begin(),attractors.end(), vpos[k]);
                    if(it!=attractors.end()) attractors.erase(it);
                }
            }*/
            for(int m=0; m<int(attractors.size());m++)
            {
                std::vector<double> ds(3);
                ds[0]=_domain.get_dx();
                ds[1]=_domain.get_dy();
                ds[2]=_domain.get_dz();
                std::cout<<"attractor "<<m<<std::endl;
                for(int a=0;a<3;a++)
                {
                    std::cout << (attractors[m][a] * ds[a] + _domain.get_origin()[a]) * Constants::BOHR_RADIUS_TO_ANGSTROM << std::endl;
                }
            }
            for(int m=0; m<_structure.number_of_atoms();m++)
            {
                std::cout<<atom_attract_diff(attractors)[m]<<std::endl;
                std::cout<<atom_attract_diff(attractors).size()<<std::endl;
            }
            
            
            return g;
        
        
        }
        else
        {
            throw std::string("Gradient grid must be used for AIM");
        }
    }
    
    catch(std::string error)
    {
        std::cout<<error<<std::endl;
        exit(1);
    }
}

void Grid::next_Density(int i, int j, int k, double& rhocenter, std::vector<std::vector<int>>& trajectory)
{
    //coordinates on grid of the latest point
    std::vector<int> v=trajectory.back();
    
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };
    
    double drhomax=0;
    std::vector<double> ds(3);
    
    // Iterate through neighboring cells
    for(int ic=0;ic<3;ic++)
    {
        for(int jc=0;jc<3;jc++)
        {    
                
            for(int kc=0;kc<3;kc++)
            {
                // Skip the center point
                if(ic==1 and jc==1 and kc==1) continue;
                
                // Calculate distances
                ds[0]=(I[ic]-I[1])*_domain.get_dx();
                ds[1]=(J[jc]-J[1])*_domain.get_dy();
                ds[2]=(K[kc]-K[1])*_domain.get_dz();
                
                // Calculate the norm of the distance std::vector
                double normds=sqrt(ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2]);
                
                // Get the density at the neighboring point
                double rho=_values[I[ic]][J[jc]][K[kc]][0];
                
                // Calculate the density gradient
                double drho = (rho-rhocenter)/normds;
                
                // Check if this is the maximum density gradient
                if(drho-drhomax>PRECISION)
                {
                    v[0]=I[ic];
                    v[1]=J[jc];
                    v[2]=K[kc];
                    drhomax = drho;
                }
            }
        }
    }
    
    //Update density for the new point
    rhocenter=_values[v[0]][v[1]][v[2]][0];
    
    //Add the new point to the ascent trajectory
    trajectory.push_back(v);
}

void Grid::addSurroundingDensity(int i,int j,int k, std::vector<std::vector<int>>& equals, double& rhocenter)
{
    //coordinates on grid of the latest point
    std::vector<int> v(3);
    
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };
    
    for(int ic=0;ic<3;ic++)
    {    
        for(int jc=0;jc<3;jc++)
        {
            for(int kc=0;kc<3;kc++)
            {
                // Skip the center point
                if(ic==1 and jc==1 and kc==1) 
                {
                    continue;
                }
                
                // Calculate distances
                std::vector<double> ds(3);
                ds[0]=(I[ic]-I[1])*_domain.get_dx();
                ds[1]=(J[jc]-J[1])*_domain.get_dy();
                ds[2]=(K[kc]-K[1])*_domain.get_dz();
                double normds=sqrt((ds[0]*ds[0])+(ds[1]*ds[1])+(ds[2]*ds[2]));
                
                // Get the density at the neighboring point
                double rho=_values[I[ic]][J[jc]][K[kc]][0];
                
                // Calculate the density gradient
                rho =(rho-rhocenter)/normds;
                
                // Check if the gradient is zero
                if(rho<PRECISION)
                {
                    v[0]=I[ic];
                    v[1]=J[jc];
                    v[2]=K[kc];
                    
                    // If it is, add point to list
                    equals.push_back(v);
                }    
            }
        }
    }
}

Grid Grid::aim_On_Grid_Density()
{

    try
    {
        if(_domain.get_Nval()==1)
        {
            // Initialise new grid for indices
            Grid g(Domain(1,_domain.get_N1(), _domain.get_N2(), _domain.get_N3(), _domain.get_origin()));
            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    g._domain.set_Tij(_domain.get_Tij(i,j), i, j);
                }
            }
            g._structure=_structure;
            //make list instead
            // Create a list of the attractors' coordinates
            std::vector<std::vector<int>> attractors(0, std::vector<int>(3));
            for(int i=0;i<_domain.get_N1();i++)
            {
                for(int j=0;j<_domain.get_N2();j++)
                {
                    for(int k=0;k<_domain.get_N3();k++)
                    {
                        // Check if the point has already been indexed
                        if(g._values[i][j][k][0]>PRECISION)
                        {
                            std::cout<<"already indexed"<<std::endl;
                            continue;
                        }
                        else
                        {    
                            // Initialise density at point i j k and the known point flag
                            double rhocenter=_values[i][j][k][0];
                            bool KnownPt=false;
                            
                            // Initialise ascent trajectory and equal density neighbours
                            std::vector<std::vector<int>> trajectory={{i,j,k}};
                            std::vector<std::vector<int>> equals(0, std::vector<int>(3));
                            
                            
                            // Start the loop to find the path of highest density gradient
                            do
                            {
                                //Initialize indices to the point
                                int I=trajectory.back()[0];
                                int J=trajectory.back()[1];
                                int K=trajectory.back()[2];

                                next_Density(I,J,K,rhocenter,trajectory);
                                
                                // Update indices to the new point
                                int newI = trajectory.back()[0];
                                int newJ = trajectory.back()[1];
                                int newK = trajectory.back()[2];
                                std::cout<<newI<<std::endl;
                                std::cout<<newJ<<std::endl;
                                std::cout<<newK<<std::endl;
                                std::cout<<rhocenter<<std::endl;
                                
                                // Check if a new point was added to the trajectory
                                if(I==newI and J==newJ and K==newK)
                                {
                                    std::cout<<"no new point found"<<std::endl;
                                    break;
                                }
                                
                                //check if the new point is already indexed
                                if(g._values[I][J][K][0]>PRECISION)
                                {
                                    KnownPt=true;
                                    break;
                                }
                            } while(true);
                            
                            // Iterate through the trajectory and set the index value of each point in the trajectory to the index of the associated maximum (last point)
                            for(int p=0;p<int(trajectory.size()-1);p++)
                            {
                                g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];
                            }
                            
                            
                            // Add points surrounding the maximum with equal density gradients to equals for indexing
                            addSurroundingDensity(trajectory.back()[0],trajectory.back()[1],trajectory.back()[2],equals,rhocenter);
                            
                            
                            // Iterate through the equal points and set the index value of each point in the the std::vector to the index of the associated maximum 
                            for(int p=0;p<int(equals.size());p++)
                            {
                                g._values[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];    
                            }
                            
                            
                            // If the point has been already indexed move on to a new point
                            if(KnownPt)
                            {
                                std::cout<<"gone"<<std::endl;
                                continue;    
                            }
                            
                            // Known attractor flag
                            bool WhoIsMax=true;

                            // Compare new attractor to known attractors. If the attractor is known index all points from trajectory and change flag value
                            for(int m=0;m<int(attractors.size());m++)
                            {
                                int I=1;
                                
                                if(attractors[m]==trajectory.back())
                                {
                                    WhoIsMax=false;
                                    for(int p=0;p<int(trajectory.size());p++)
                                    {
                                        g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(I);
                                    }
                                }
                                I++;
                            }
                            
                            // If the attractor is unknown
                            if(WhoIsMax)
                            {    
                            
                                //if the number of points leading to the attractor is insufficient, index trajectory and equals to previous attractor and disregard current
                                if(int(trajectory.size())<(_domain.get_N1()*_domain.get_N2()*_domain.get_N3())/1000)
                                {
                                    
                                    for(int p=0;p<int(equals.size());p++)
                                    {
                                        g._values[equals[p][0]][equals[p][1]][equals[p][2]][0]=double(attractors.size());        
                                    }
                                    for(int p=0;p<int(trajectory.size());p++)
                                    {
                                        g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
                                    }
                                    
                                }
                                
                                // Add the attractor to the list of attractors
                                attractors.push_back(trajectory.back());
                                
                                
                                std::cout<<"traj size "<<trajectory.size()<<std::endl;
                                std::cout<<"attract size "<<attractors.size()<<std::endl;
                                
                                // Index all points from trajectory and equals
                                for(int p=0;p<int(equals.size());p++)
                                {
                                    g._values[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._values[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];        
                                }
                                for(int p=0;p<int(trajectory.size());p++)
                                {
                                    g._values[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
                                }
                            }
                        
                        }                    
                    }
                }
            }
            /*
            std::vector<std::vector<int>> vpos(0,std::vector<int>(3));
            for(int p=0;p<int(attractors.size()); p++)
            {
                for(int m=0;m<int(attractors.size());m++)
                {
                    int I[2];
                    I[0]=-1;
                    I[1]=1;
                    for(int i=0;i<2;i++)
                    {
                        if(attractors[p][0]==attractors[m][0] +I[i] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2])
                        {
                            vpos.push_back(attractors[m]);
                            for(int l=nBound;l<_dom.N1()-nBound;l++)
                            {
                                for(int j=nBound;j<_dom.N2()-nBound;j++)
                                {
                                    for(int k=nBound;k<_dom.N3()-nBound;k++)
                                    {
                                        if(g._values[l][j][k][0]==m)
                                        {
                                            std::cout<<"okay"<<std::endl;
                                            g._values[l][j][k][0]=p;
                                            
                                        }
                                    }
                                }
                            }
                        }
                        if( attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] +I[i] and attractors[p][1]==attractors[m][1])
                        {
                            vpos.push_back(attractors[m]);
                            for(int l=nBound;l<_dom.N1()-nBound;l++)
                            {
                                for(int j=nBound;j<_dom.N2()-nBound;j++)
                                {
                                    for(int k=nBound;k<_dom.N3()-nBound;k++)
                                    {
                                        if(g._values[l][j][k][0]==m)
                                        {
                                            std::cout<<"okay"<<std::endl;
                                            g._values[l][j][k][0]=p;
                                            
                                        }
                                    }
                                }
                            }
                        }
                        if(attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2] +I[i])
                        {
                            vpos.push_back(attractors[m]);
                            for(int l=nBound;l<_dom.N1()-nBound;l++)
                            {
                                for(int j=nBound;j<_dom.N2()-nBound;j++)
                                {
                                    for(int k=nBound;k<_dom.N3()-nBound;k++)
                                    {
                                        if(g._values[l][j][k][0]==m)
                                        {
                                            std::cout<<"okay"<<std::endl;
                                            g._values[l][j][k][0]=p;
                                            
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for(int k=0;k<int(vpos.size());k++)
            {
                for(auto element:attractors)
                {    
                    auto it=find(attractors.begin(),attractors.end(), vpos[k]);
                    if(it!=attractors.end()) attractors.erase(it);
                }
            }*/
            for(int m=0; m<int(attractors.size());m++)
            {
                std::vector<double> ds(3);
                ds[0]=_domain.get_dx();
                ds[1]=_domain.get_dy();
                ds[2]=_domain.get_dz();
                std::cout<<"attractor "<<m<<std::endl;
                for(int a=0;a<3;a++)
                {
                    std::cout << (attractors[m][a] * ds[a] + _domain.get_origin()[a]) * Constants::BOHR_RADIUS_TO_ANGSTROM << std::endl;
                    std::cout<<attractors[m][a]<<std::endl;
                }
            }
            for(int m=0; m<_structure.number_of_atoms();m++)
            {
                std::cout<<atom_attract_diff(attractors)[m]<<std::endl;
                std::cout<<atom_attract_diff(attractors).size()<<std::endl;
            }
            //int nbrofindices=0;
            double value=0;
            std::cout<<"V0000"<<g._values[0][0][0][0]<<std::endl;
            for(int i=0;i<_domain.get_N1();i++)
            {
                for(int j=0;j<_domain.get_N2();j++)
                {
                    for(int k=0;k<_domain.get_N3();k++)
                    {
                        if(g._values[i][j][k][0]>value)
                        {
                            value =g._values[i][j][k][0];
                        }
                    }
                }
            }
            std::cout<<"number of of idices "<<int(value)<<std::endl;
            return g;
        
        
        }
        else
        {
            throw std::string("Gradient grid must be used for AIM");
        }
    }
    
    catch(std::string error)
    {
        std::cout<<error<<std::endl;
        exit(1);
    }
}

double Grid::value(double x, double y, double z) const 
{
    int i = _domain.i(x,y,z);
    int j = _domain.j(x,y,z);
    int k = _domain.k(x,y,z);
    
    if(i>=_domain.get_N1()-1 or j>=_domain.get_N2()-1 or k>=_domain.get_N3( )-1 or i<0 or j<0 or k<0)
    {
        return 0;
    }
    
    int I[2] = {i, i+1};
    int J[2] = {j, j+1};
    int K[2] = {k, k+1};
    
    double norm=0;
    double S=0;
    for(int ic=0;ic<2;ic++)
    {    
        for(int jc=0;jc<2;jc++)
        {
            for(int kc=0;kc<2;kc++)
            {
                double ex = exp(-sqrt((x-_domain.x(I[ic],J[jc],K[kc]))*(x-_domain.x(I[ic],J[jc],K[kc]))+(y-_domain.y(I[ic],J[jc],K[kc]))*(y-_domain.y(I[ic],J[jc],K[kc]))+(z-_domain.z(I[ic],J[jc],K[kc]))*(z-_domain.z(I[ic],J[jc],K[kc]))));
                S += _values[I[ic]][J[jc]][K[kc]][0]*ex;
                norm+=ex;
            }
        }
    }
    if(norm>0)
        S/= norm;
    return S;
}
