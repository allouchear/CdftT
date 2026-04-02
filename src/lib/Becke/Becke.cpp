#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Cube/Grid.h>
#include <Becke/Becke.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Becke::Becke()
{
	_molecule = Structure();
	_grid = GridPoints();
    _orbitals = Orbitals();
    _grid_points = std::vector<std::vector<std::vector<double>>> ();
    _grid_weights = std::vector<std::vector<double>> ();
    _grid_volumes = std::vector<std::vector<double>> ();
    _multigrid = false;
    
    _energy = 0.0;
}

Becke::Becke(const Structure& S)
{
	_molecule = S;
	_grid = GridPoints();
	_orbitals = Orbitals();
	_grid_points = std::vector<std::vector<std::vector<double>>> ();
	_grid_weights = std::vector<std::vector<double>> ();
	_grid_volumes = std::vector<std::vector<double>> ();
	_multigrid = false;
	_energy = 0.0;
}

Becke::Becke(const Grid& g)
{
	_molecule = g.get_structure();
	_grid = GridPoints();
	_orbitals = Orbitals();
	_grid_points = std::vector<std::vector<std::vector<double>>> ();
	_grid_weights = std::vector<std::vector<double>> ();
	_grid_volumes = std::vector<std::vector<double>> ();
	_multigrid = false;
	_energy = 0.0;
}

Becke::Becke(WFX& wfx, Binomial& bin, const PeriodicTable& table)
{
	_molecule = Structure(wfx, table);
	_grid = GridPoints();
    _orbitals = Orbitals(wfx, bin, table);
    _grid_points = std::vector<std::vector<std::vector<double>>> ();
    _grid_weights = std::vector<std::vector<double>> ();
    _grid_volumes = std::vector<std::vector<double>> ();
    _multigrid = false;
    _energy = wfx.Energy();
}

Becke::Becke(FCHK& fchk, Binomial& bin, const PeriodicTable& table)
{
    _molecule = Structure(fchk, table);
    _grid = GridPoints();
    _orbitals = Orbitals(fchk, bin, table);
    _grid_points = std::vector<std::vector<std::vector<double>>> ();
    _grid_weights = std::vector<std::vector<double>> ();
    _grid_volumes = std::vector<std::vector<double>> ();
    _multigrid = false;
    _energy = fchk.ScfEnergy();
}

Becke::Becke(MOLDENGAB& moldengab, Binomial& bin, const PeriodicTable& table)
{
    _molecule = Structure(moldengab, table);
    _grid = GridPoints();
    _orbitals = Orbitals(moldengab, bin, table);
    _grid_points = std::vector<std::vector<std::vector<double>>> ();
    _grid_weights = std::vector<std::vector<double>> ();
    _grid_volumes = std::vector<std::vector<double>> ();
    _multigrid = false;
    _energy = 0;
}

Becke::Becke(LOG& log, Binomial& bin, const PeriodicTable& table)
{
    _molecule = Structure(log, table);
    _grid = GridPoints();
    _orbitals = Orbitals(log, bin, table);
    _grid_points = std::vector<std::vector<std::vector<double>>> ();
    _grid_weights = std::vector<std::vector<double>> ();
    _grid_volumes = std::vector<std::vector<double>> ();
    _multigrid = false;
    _energy = log.Energy();
}


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

double Becke::get_energy() const
{
	return _energy;
}

const Structure& Becke::get_molecule() const
{
    return _molecule;
}

const Orbitals& Becke::get_orbitals() const
{
    return _orbitals;
}

const std::vector<double>& Becke::get_partial_charge() const
{
    return _partial_charge;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

double Becke::getHOMOEnergy()
{
    _orbitals.HOMO();
    return _orbitals.getHOMOEnergy();
}

double Becke::getLUMOEnergy()
{
    _orbitals.LUMO();
    return _orbitals.getLUMOEnergy();
}

int Becke::number_of_radial_points(int Z)
{
    /*
    select the number of radial grid points for the subintegral
    around an atom with atomic number Z
    */
    // Hydrogen atom receives an initial quota of 20 points
    int Nr = 20;
    // Each shell receives an additional 5 points
    if(Z >= 2)
        Nr += 5;
    if(Z >= 11)
        Nr += 5;
    if(Z >= 19)
        Nr += 5;
    if(Z >= 37)
        Nr += 5;
    if(Z >= 55)
        Nr += 5;
    if(Z >= 87)
        Nr += 5;

    return Nr;
}

GridPoints Becke::select_angular_grid(int lebedev_order)
{
	std::vector<int> L = _grid.Lebedev_Lmax();

    for(size_t i = 0; i < _grid.Lebedev_Lmax().size(); i++)
    {
        L[i] = std::abs(L[i] - lebedev_order);
    }
    
    int n = 0;
    int Ln = L[0];

	for(size_t i = 1; i < _grid.Lebedev_Lmax().size(); i++)
    {
		if(Ln > L[i])
        {
			Ln = L[i];
            n = i;
        }
    }

	int Lmax = _grid.Lebedev_Lmax()[n];
    if(lebedev_order != Lmax)
    {
        std::cout << "Error, no grid found." << std::endl;
        
        std::exit(1);
    }

	return GridPoints(Lmax);
}

double Becke::s(double mu, int k)
{
    // for nuclear weight functions
    double f = mu;
    for(int ik=0; ik<k; ik++)
      	f = 1.5 * f -0.5 * f*f*f;
    return 0.5*(1-f);
}

void Becke::multicenter_grids(int kmax, int lebedev_order, int radial_grid_factor)
{
    /*
    compute grid points and weights of the multicenter grids for visualization
   
    Parameters
    ----------
    atomlist           : list of tuples (Zat,(xI,yI,zI)) with atomic numbers and 
                         atom positions, which define the multicenter grid
    
    Optional
    --------
    kmax               : How fuzzy should the Voronoi polyhedrons be? Larger kmax
                         means borders are fuzzier.
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    grid_points        : list of tuples (x,y,z) with positions of points in each grid,
                         grid_points[I][0] contains the x-positions of the points
                         belonging to the grid around atom I
    grid_weights       : list of numpy arrays, grid_weights[I][k] contains the weight
                         of the k-th point in the grid around atom I due to the fuzzy
                         Voronoi decomposition.
    grid_volumes       : list of numpy arrays, grid_volumes[I][k] contains the volume
                         element around the k-th point in the grid at atom I.
    */

    // angular grid
    GridPoints Angular = select_angular_grid(lebedev_order);
    int Nang = Angular.Npts();
    std::vector<double> thang (Nang);
    std::vector<double> phiang (Nang);
    std::vector<double> wang (Nang);
    std::vector<double> sc (Nang);
    std::vector<double> ss (Nang);
    std::vector<double> c (Nang);

    for(int i=0; i<Nang; i++)
    {
        thang[i]=Angular.LebedevGridPoints()[i][0];
        phiang[i]=Angular.LebedevGridPoints()[i][1];
        wang[i] = Angular.LebedevGridPoints()[i][2];

        sc[i] = std::sin(thang[i]) * std::cos(phiang[i]);
        ss[i] = std::sin(thang[i]) * std::sin(phiang[i]);
        c[i]  = std::cos(thang[i]);
    }

    // declaration of all variables
    int Nat = _molecule.getNumberOfAtoms();
    double chi, uij, rm, mu, nu;
    std::vector<double> wr(Nat, 0.0);
    std::vector<std::vector<double>> R (Nat,std::vector<double> (Nat,0));     // distances between atoms i and j
    std::vector<std::vector<double>> a (Nat,std::vector<double> (Nat,0));     // scaling factor used in eqn. A2
    std::vector<std::vector<std::vector<double>>> list_coord_I (Nat);
    std::vector<std::vector<std::vector<double>>> grid_points;
    std::vector<std::vector<double>> grid_weights;
    std::vector<std::vector<double>> grid_volumes;
    int Nr, Npts;

    //cout<<endl;
    for(int i=0; i<Nat; i++)
        for(int j=i+1; j<Nat; j++)
        {
            R[i][j] = R[j][i] = _molecule.atom(i).computeDistance(_molecule.atom(j));
            // ratio of Slater radii
            chi = _molecule.atom(i).get_covalentRadius() / _molecule.atom(j).get_covalentRadius();
            uij = (chi-1)/(chi+1);
            a[i][j] = uij/(uij*uij - 1);
            a[j][i] = -a[i][j];
            //cout<<"a["<<i<<"]["<<j<<"] = "<<a[i][j]<<endl;
        }
    
    // atom-centered subintegral
    for(int I=0; I<Nat; I++)
    {
        // radial grid
        Nr = number_of_radial_points(_molecule.atom(I).get_atomicNumber());
        // increase number of grid points is requested
        Nr *= radial_grid_factor;
        rm = 0.5*_molecule.atom(I).get_covalentRadius();
        Npts = Nr*Nang;

        std::vector<double> k (Nr);
        std::vector<double> xr (Nr);
        std::vector<double> radial_weights (Nr);
        std::vector<double> g (Nr);
        std::vector<double> r (Nr);
        std::vector<double> wr(Npts, 0.0);
        std::vector<double> x(Npts), y(Npts) , z(Npts), weights(Npts);
        std::vector<double> Ptot(Npts, 0.0);
        std::vector<std::vector<double>> dist (Nat, std::vector<double> (Npts,0.0));
        std::vector<std::vector<double>> P (Nat, std::vector<double> (Npts,1.0));
        double div=0.0;

        for(int J=0; J<Nr; J++)
        {
            // grid points on interval [-1,1]
            xr[J] = std::cos((J + 1) / (Nr + 1.0) * M_PI);
            // weights
            radial_weights[J] = M_PI / (Nr + 1.0) * std::sin((J + 1) / (Nr + 1.0) * M_PI) * std::sin((J + 1) / (Nr + 1.0) * M_PI);
            // from variable transformation
            div=(1+xr[J])/((1-xr[J])*(1-xr[J])*(1-xr[J]));
            g[J] = 2 * rm * rm * rm * std::sqrt(div * div * div);
            radial_weights[J] *= g[J];
            // radial grid points on interval [0,infinity]
            r[J] = rm * (1+xr[J])/(1-xr[J]);
        }
        
        int n=0;
            // cartesian coordinates of grid
        for (int i=0; i<Nr; i++)
            for(int j=0; j<Nang; j++)
            {
                x[n] = r[i] * sc[j] + _molecule.atom(I).get_coordinates()[0];
                y[n] = r[i] * ss[j] + _molecule.atom(I).get_coordinates()[1];
                z[n] = r[i] * c[j] + _molecule.atom(I).get_coordinates()[2];
                weights[n] = radial_weights[i] * 4.0 * M_PI * wang[j];
                n++;
            }

        // distance between grid points and atom i
        for (int i=0; i<Nat; i++)
            for(int j=0; j<Npts; j++)
                dist[i][j] = std::sqrt((x[j] - _molecule.atom(i).get_coordinates()[0]) * (x[j] - _molecule.atom(i).get_coordinates()[0])
                                        + (y[j] - _molecule.atom(i).get_coordinates()[1]) * (y[j] - _molecule.atom(i).get_coordinates()[1])
                                        + (z[j] - _molecule.atom(i).get_coordinates()[2]) * (z[j] - _molecule.atom(i).get_coordinates()[2]));
                

        // P_i(r) as defined in eqn. (13)
        
        for(int i=0; i<Nat; i++)
            for(int j=0; j<Nat; j++)
            {
                // mu_ij as defined in eqn. (11)
                if(i!=j)
                    for(int k=0; k<Npts; k++)
                    {
                        mu = (dist[i][k]-dist[j][k])/R[i][j];
                        nu = mu + a[i][j]*(1-mu*mu);
                        P[i][k] *= s(nu, kmax);
                    }
            }

        for(int i=0; i<Npts; i++)
            for(int j=0; j<Nat; j++)
    	       Ptot[i] += P[j][i];
            
        // weight function due to partitioning of volume
        for(int i=0; i<Npts; i++)
            wr[i] = P[I][i]/Ptot[i];
        

        std::vector<std::vector<double>> coord_i (Npts);
        for(int i=0; i<Npts; i++)
            coord_i[i]={x[i], y[i], z[i]};

        list_coord_I[I]=coord_i;
        // The weights come from the fuzzy Voronoi partitioning 
        grid_weights.push_back( wr );
        // The naming is a little bit confusing, the `weights` are
        // actually the volume elements dV_i around each point.
        grid_volumes.push_back( weights );
    }

    grid_points=list_coord_I;

	_grid_points=grid_points;
    _grid_weights=grid_weights;
    _grid_volumes=grid_volumes;

    long int printnpts=0;
    for(int i=0; i<Nat; i++)
        printnpts+=_grid_points[i].size();

    std::cout << "Number of Becke grid points: " << printnpts << std::endl << std::endl;
}


std::vector<std::vector<double>> Becke::join_grids()
{   /*
    combine the multicenter grids into a single grid so that we get
    a quadrature rule for integration
         /
         | f(x,y,z) dV = sum  w  f(x ,y ,z ) 
         /                  i  i    i  i  i

    Parameters
    ----------
    points, weights, volumes:  return values of `multicenter_grids`

    Returns
    -------
    x,y,z     :  1d numpy arrays with cartesian coordinates of grid points
    w         :  1d numpy array with weights
    */
    // weights of quadrature rule
    std::vector<double> w;
    // sampling points of quadrature rule
    std::vector<double> x, y, z;
    // xI,yI,zI : grid points of spherical grid around atom I
    // wI : weights of spherical grid around atom I
    // vI : volume elements of spherical grid around atom I
    for(size_t i=0; i<_grid_points.size(); i++)
        for(size_t j=0; j<_grid_points[i].size(); j++)
        {
            x.push_back(_grid_points[i][j][0]);
            y.push_back(_grid_points[i][j][1]);
            z.push_back(_grid_points[i][j][2]);

            // The weights are the product of the weight function
            // (from fuzzy Voronoi decomposition of space) and the volume element.
            w.push_back(0);
            w[i]+=_grid_weights[i][j]*_grid_volumes[i][j];
        }

    std::vector<std::vector<double>> join_grid;
    join_grid.push_back(x);
    join_grid.push_back(y);
    join_grid.push_back(z);
    join_grid.push_back(w);

    return join_grid;
}

double Becke::multicenter_integration(std::function<double(const std::vector<GTF>& p, double x, double y, double z)> f, const std::vector<GTF>& p, int kmax, int lebedev_order, int radial_grid_factor)
{    /*
    compute the integral

             / 
         I = | f(x,y,z) dV
             / 

    numerically on a multicenter spherical grid using Becke's scheme 
   
    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    I       : value of the integral
    */

    if(_multigrid==false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid=true;
    }

    int Nat = _molecule.getNumberOfAtoms();

    double integral = 0.0;

    for(int I=0; I<Nat; I++)
    {
        double integ=0.0;
        #ifdef ENABLE_OMP
        #pragma omp parallel for reduction(+:integ)
        #endif
        for(size_t J=0; J<_grid_weights[I].size(); J++)
            integ += _grid_volumes[I][J] * _grid_weights[I][J] * f(p, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2]);       // evaluate function on the grid
        integral+=integ;
    }

    return integral;
}

double Becke::multicenter_integration(std::function<double(const Orbitals&, int i, int j, double x, double y, double z)> f, int i, int j, int kmax, int lebedev_order, int radial_grid_factor)
{    /*
    compute the integral

             / 
         I = | f(x,y,z) dV
             / 

    numerically on a multicenter spherical grid using Becke's scheme 
   
    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    I       : value of the integral
    */
    if(_multigrid==false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid=true;
    }

    int Nat = _molecule.getNumberOfAtoms();

    double integral = 0.0;

    for(int I=0; I<Nat; I++)
    {
        double integ=0.0;
        #ifdef ENABLE_OMP
        #pragma omp parallel for reduction(+:integ)
        #endif
        for(size_t J=0; J<_grid_weights[I].size(); J++)
            integ += _grid_volumes[I][J] * _grid_weights[I][J] * f(_orbitals, i, j, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2]);       // evaluate function on the grid
        integral+=integ;
    }

    return integral;
}

double Becke::multicenter_integration(std::function<double(const Orbitals&, int, int, double, double, double, SpinType)> f, int i, int j, int kmax, int lebedev_order, int radial_grid_factor, SpinType spinType)
{    /*
    compute the integral

             / 
         I = | f(x,y,z) dV
             / 

    numerically on a multicenter spherical grid using Becke's scheme 
   
    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    I       : value of the integral
    */

    if(_multigrid==false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid=true;
    }

    int Nat = _molecule.getNumberOfAtoms();

    double integral = 0.0;

    for(int I=0; I<Nat; I++)
    {
        double integ=0.0;
        #ifdef ENABLE_OMP
        #pragma omp parallel for reduction(+:integ)
        #endif
        for(size_t J=0; J<_grid_weights[I].size(); J++)
            integ += _grid_volumes[I][J] * _grid_weights[I][J] * f(_orbitals, i, j, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2], spinType);       // evaluate function on the grid

        integral+=integ;
    }

    return integral;
}

double Becke::multicenter_integration(std::function<double(Orbitals&, int, int, double, double, double, SpinType, const std::array<double, 3>&, double)> f, int i, int j, SpinType spinType, const std::array<double, 3>& chargePosition, double charge, int kmax, int lebedev_order, int radial_grid_factor)
{
    if (_multigrid == false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid = true;
    }

    int Nat = _molecule.getNumberOfAtoms();

    double integral = 0.0;

    for (int I = 0; I < Nat; I++)
    {
        double integ = 0.0;

        #ifdef ENABLE_OMP
        #pragma omp parallel for reduction(+ : integ)
        #endif
        for (size_t J = 0; J < _grid_weights[I].size(); ++J)
        {
            integ += _grid_volumes[I][J] * _grid_weights[I][J] * f(_orbitals, i, j, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2], spinType, chargePosition, charge); // evaluate function on the grid
        }
        
        integral += integ;
    }

    return integral;
}

std::vector<double> Becke::multicenter_sub_integration(std::function<double(const Orbitals&, double, double, double)> f, int kmax, int lebedev_order, int radial_grid_factor)
{   
    /*
    compute the integral

             / 
         I = | f(x,y,z) dV
             / 

    numerically on a multicenter spherical grid using Becke's scheme 
   
    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    I       : value of the integral
    */

    if(_multigrid==false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid=true;
    }

    int Nat = _molecule.getNumberOfAtoms();

    std::vector<double> sub_integral (Nat, 0.0);

    for(int I=0; I<Nat; I++)
    {
        double integ=0.0;
        #ifdef ENABLE_OMP
        #pragma omp parallel for reduction(+:integ)
        #endif
        for(size_t J=0; J<_grid_weights[I].size(); J++)
            integ += _grid_volumes[I][J] * _grid_weights[I][J] * f(_orbitals, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2]);       // evaluate function on the grid

        sub_integral[I] = integ;
    }

    return sub_integral;
}

double Becke::OverlapGTF(const GTF& A, const GTF& B, int kmax, int lebedev_order, int radial_grid_factor)
{
    /*
    overlap between two basis functions
    
        (a|b)

    Parameters
    ----------
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  callables, atomic basis functions
                      e.g. bfA(x,y,z) etc.

    Returns
    -------
    Sab            :  float, overlap integral    
    */
    // Now we compute the integrals numerically on a multicenter grid.
    // 1. define integrand s = a b
    // def s_integrand(x,y,z):
    //    return bfA(x,y,z).conjugate() * bfB(x,y,z)
    //                                    
    // 2. integrate density on a multicenter grid
    std::vector<GTF> p (2);
    p[0]=A;
    p[1]=B;

    double Sab = multicenter_integration(&prodGTF, p, kmax, lebedev_order, radial_grid_factor);

    return Sab;
}

double Becke::OverlapCGTF(int i, int j, int kmax, int lebedev_order, int radial_grid_factor)
{
    /*
    overlap between two basis functions
    
        (a|b)

    Parameters
    ----------
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  callables, atomic basis functions
                      e.g. bfA(x,y,z) etc.

    Returns
    -------
    Sab            :  float, overlap integral    
    */
    // Now we compute the integrals numerically on a multicenter grid.
    // 1. define integrand s = a b
    // def s_integrand(x,y,z):
    //    return bfA(x,y,z).conjugate() * bfB(x,y,z)
    //                                    
    // 2. integrate density on a multicenter grid

    double Sab = multicenter_integration(&CGTFstarCGTF, i, j, kmax, lebedev_order, radial_grid_factor);

    return Sab;
}

double Becke::overlap(int i, int j, int kmax, int lebedev_order, int radial_grid_factor, SpinType spinType)
{
    /*
    overlap between two basis functions
    
        (a|b)

    Parameters
    ----------
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  callables, atomic basis functions
                      e.g. bfA(x,y,z) etc.

    Returns
    -------
    Sab            :  float, overlap integral    
    */
    // Now we compute the integrals numerically on a multicenter grid.
    // 1. define integrand s = a b
    // def s_integrand(x,y,z):
    //    return bfA(x,y,z).conjugate() * bfB(x,y,z)
    //                                    
    // 2. integrate density on a multicenter grid

    double Sab = multicenter_integration(&phiStarPhi, i, j, kmax, lebedev_order, radial_grid_factor, spinType);

    return Sab;
}

void Becke::partial_charge(int kmax, int lebedev_order, int radial_grid_factor)
{
    int Nat=_molecule.getNumberOfAtoms();
    std::vector<double> qn (Nat);
    std::vector<double> In = multicenter_sub_integration(static_cast<double(*)(const Orbitals&, double, double, double)>(density), kmax, lebedev_order, radial_grid_factor); // must use the cast because density is an overloaded function

    for(int i=0; i<Nat ;i++)
        qn[i]=_molecule.atom(i).get_atomicNumber() - In[i];

    _partial_charge=qn;
}

void Becke::partial_charge(const Grid &g, int kmax, int lebedev_order, int radial_grid_factor)
{
    int Nat = _molecule.getNumberOfAtoms();
    std::vector<double> qn(Nat);
    std::vector<double> In = multicenter_sub_integration(g);
    for (int i = 0; i < Nat; i++)
    {
        qn[i] = _molecule.atom(i).get_atomicNumber() - In[i];
    }
    _partial_charge = qn;
}

std::vector<std::vector<std::vector<double>>> Becke::getIonicPotentialMatrix(const std::array<double, 3>& chargePosition, double charge, int kmax, int lebedev_order, int radial_grid_factor, bool debug, bool printAOMatrix, bool printMOMatrix)
{
    if (_multigrid == false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid = true;
    }

    int numberOfAo = _orbitals.get_numberOfAo();
    int numberOfMo = _orbitals.get_numberOfMo();
    std::vector<CGTF> vcgtf = _orbitals.get_vcgtf();
    const std::vector<std::vector<std::vector<double>>>& coefficients = _orbitals.get_coefficients();

    // debug
    if (__debug_AOMatrix.size() == 0)
    {
        __debug_AOMatrix = std::vector<std::vector<double>>(numberOfAo, std::vector<double>());
    }

    // Build ionic potential matrix in AO basis
    std::vector<std::vector<double>> ionicMatrixAO = std::vector<std::vector<double>>(numberOfAo, std::vector<double>());
    for (int i = 0; i < numberOfAo; ++i)
    {
        ionicMatrixAO[i].resize(i + 1, 0.0);

        // debug
        __debug_AOMatrix[i].resize(i + 1, 0.0);
    }

    // Compute ionic potential matrix elements in AO basis
    for (int i = 0; i < numberOfAo; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            for (int I = 0; I < _molecule.getNumberOfAtoms(); ++I)
            {
                for (size_t J = 0; J < _grid_weights[I].size(); J++)
                {
                    double x = _grid_points[I][J][0];
                    double y = _grid_points[I][J][1];
                    double z = _grid_points[I][J][2];

                    double distance = std::sqrt((chargePosition[0] - x) * (chargePosition[0] - x)
                                                + (chargePosition[1] - y) * (chargePosition[1] - y)
                                                + (chargePosition[2] - z) * (chargePosition[2] - z));
                    
                    if (distance > 1e-10)
                    {
                        ionicMatrixAO[i][j] += _grid_weights[I][J] * _grid_volumes[I][J] * vcgtf[i].func(x, y, z) * vcgtf[j].func(x, y, z) / distance;
                    }
                }
            }

            ionicMatrixAO[i][j] *= (- charge);
        }
    }

    // debug
    if (debug)
    {
        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        if (printAOMatrix) std::cout << "Ionic potential matrix in AO basis:" << std::endl;
        for (int ii = 0; ii < numberOfAo; ++ii)
        {
            for (int jj = 0; jj <= ii; ++jj)
            {
                __debug_AOMatrix[ii][jj] += ionicMatrixAO[ii][jj];
                __debug_totalSumAO += (ii == jj ? ionicMatrixAO[ii][jj] : 2.0 * ionicMatrixAO[ii][jj]);
                if (printAOMatrix) std::cout << std::right << std::setw(17) << __debug_AOMatrix[ii][jj] << '\t';
            }
            if (printAOMatrix) std::cout << std::endl;
        }
        if (printAOMatrix) std::cout << std::defaultfloat << "Total sum of AO matrix elements: " << std::setprecision(10) << __debug_totalSumAO << std::endl << std::endl;
    }

    // debug
    if (__debug_MOMatrix.size() == 0)
    {
        __debug_MOMatrix = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(numberOfMo, std::vector<double>()));
    }

    // Build ionic potential matrix in MO basis
    std::vector<std::vector<std::vector<double>>> ionicMatrixMO = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(numberOfMo, std::vector<double>()));
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int i = 0; i < numberOfMo; ++i)
        {
            ionicMatrixMO[spin][i].resize(i + 1, 0.0);

            // debug
            __debug_MOMatrix[spin][i].resize(i + 1, 0.0);
        }
    }

    // Compute ionic potential matrix elements in MO basis
    int spin, i, j;
    #ifdef ENABLE_OMP
    #pragma omp parallel for private(spin, i, j)
    #endif
    for (spin = 0; spin < 2; ++spin)
    {
        for (i = 0; i < numberOfMo; ++i)
        {
            for (j = 0; j <= i; ++j)
            {
                double sum = 0.0;

                for (size_t m = 0; m < coefficients[spin][i].size(); ++m)
                {
                    for (size_t n = 0; n < coefficients[spin][j].size(); ++n)
                    {
                        sum += coefficients[spin][i][m] * coefficients[spin][j][n] * (n <= m ? ionicMatrixAO[m][n] : ionicMatrixAO[n][m]);
                    }
                }

                ionicMatrixMO[spin][i][j] = sum;
            }
        }
    }

    // debug
    if (debug)
    {
        for (spin = 0; spin < 2; ++spin)
        {
            std::cout << std::scientific;
            std::cout << std::setprecision(10);

            if (printMOMatrix)
                std::cout << "Ionic potential matrix in MO basis for " << to_string(static_cast<SpinType>(spin)) << " spin:" << std::endl;
            for (int ii = 0; ii < numberOfMo; ++ii)
            {
                for (int jj = 0; jj <= ii; ++jj)
                {
                    __debug_MOMatrix[spin][ii][jj] += ionicMatrixMO[spin][ii][jj];
                    __debug_totalSumMO[spin] += (ii == jj ? ionicMatrixMO[spin][ii][jj] : 2.0 * ionicMatrixMO[spin][ii][jj]);

                    if (printMOMatrix)
                        std::cout << std::right << std::setw(17) << __debug_MOMatrix[spin][ii][jj] << '\t';
                }
                if (printMOMatrix)
                    std::cout << std::endl;
            }

            if (printMOMatrix)
                std::cout << std::defaultfloat << "Total sum of MO matrix elements for " << to_string(static_cast<SpinType>(spin)) << " spin: " << std::setprecision(10) << __debug_totalSumMO[spin] << std::endl;
        }
    }

    return ionicMatrixMO;
}

std::vector<std::vector<std::vector<std::vector<double>>>> Becke::getTripleOrbitalIntegralMatrix(int kmax, int lebedev_order, int radial_grid_factor)
{
    if (_multigrid == false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid = true;
    }

    int numberOfAo = _orbitals.get_numberOfAo();
    int numberOfMo = _orbitals.get_numberOfMo();
    std::vector<CGTF> vcgtf = _orbitals.get_vcgtf();
    const std::vector<std::vector<std::vector<double>>>& coefficients = _orbitals.get_coefficients();

    // Build triple-orbital-integral matrix in AO basis
    std::vector<std::vector<std::vector<double>>> toiMatrixAO = std::vector<std::vector<std::vector<double>>>(numberOfAo, std::vector<std::vector<double>>());
    for (int i = 0; i < numberOfAo; ++i)
    {
        toiMatrixAO[i] = std::vector<std::vector<double>>(i + 1, std::vector<double>());

        for (int j = 0; j <= i; ++j)
        {
            toiMatrixAO[i][j] = std::vector<double>(j + 1, 0.0);
        }
    }

    // Compute TOI matrix elements in AO basis
    for (int i = 0; i < numberOfAo; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            for (int k = 0; k <= j; ++k)
            {
                for (int I = 0; I < _molecule.getNumberOfAtoms(); ++I)
                {
                    for (size_t J = 0; J < _grid_weights[I].size(); J++)
                    {
                        double x = _grid_points[I][J][0];
                        double y = _grid_points[I][J][1];
                        double z = _grid_points[I][J][2];
                        
                        toiMatrixAO[i][j][k] += _grid_weights[I][J] * _grid_volumes[I][J] * vcgtf[i].func(x, y, z) * vcgtf[j].func(x, y, z) * vcgtf[k].func(x, y, z);
                    }
                }
            }
        }
    }

    /*
    std::cout << std::setprecision(10);
    for (int ii = 0; ii < numberOfAo; ++ii)
    {
        std::cout << "ii = " << ii << std::endl;

        for (int jj = 0; jj <= ii; ++jj)
        {
            for (int kk = 0; kk <= jj; ++kk)
            {
                std::cout << std::right << std::setw(17) << toiMatrixAO[ii][jj][kk] << ' ';
            }

            std::cout << std::endl;
        }
    }
    */
    
    // Build TOI matrix in MO basis
    std::vector<std::vector<std::vector<std::vector<double>>>> toiMatrixMO = std::vector<std::vector<std::vector<std::vector<double>>>>(2, std::vector<std::vector<std::vector<double>>>(numberOfMo, std::vector<std::vector<double>>()));
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int i = 0; i < numberOfMo; ++i)
        {
            toiMatrixMO[spin][i].resize(i + 1);

            for (int j = 0; j <= i; ++j)
            {
                toiMatrixMO[spin][i][j].resize(j + 1);
            }
        }
    }

    // Compute TOI matrix elements in MO basis
    int spin, i, j, k;
    #ifdef ENABLE_OMP
    #pragma omp parallel for private(spin, i, j, k)
    #endif
    for (spin = 0; spin < 2; ++spin)
    {
        for (i = 0; i < numberOfMo; ++i)
        {
            for (j = 0; j <= i; ++j)
            {
                for (k = 0; k <= j; ++k)
                {
                    double sum = 0.0;

                    for (size_t p = 0; p < coefficients[spin][i].size(); ++p)
                    {
                        for (size_t q = 0; q < coefficients[spin][j].size(); ++q)
                        {
                            for (size_t r = 0; r < coefficients[spin][k].size(); ++r)
                            {
                                // Get the correct AO matrix element (considering symmetry)
                                std::array<size_t, 3> indices = {p, q, r};
                                std::sort(indices.begin(), indices.end(), std::greater<size_t>());
                                double toiAOElement = toiMatrixAO[indices[0]][indices[1]][indices[2]];

                                sum += coefficients[spin][i][p] * coefficients[spin][j][q] * coefficients[spin][k][r] * toiAOElement;
                            }
                        }
                    }

                    toiMatrixMO[spin][i][j][k] = sum;
                }
            }
        }
    }

    // debug
    // if (debug)
    // {
    //     for (spin = 0; spin < 2; ++spin)
    //     {
    //         std::cout << std::scientific;
    //         std::cout << std::setprecision(10);

    //         if (printMOMatrix)
    //             std::cout << "Ionic potential matrix in MO basis for " << to_string(static_cast<SpinType>(spin)) << " spin:" << std::endl;
    //         for (int ii = 0; ii < numberOfMo; ++ii)
    //         {
    //             for (int jj = 0; jj <= ii; ++jj)
    //             {
    //                 __debug_MOMatrix[spin][ii][jj] += ionicMatrixMO[spin][ii][jj];
    //                 __debug_totalSumMO[spin] += (ii == jj ? ionicMatrixMO[spin][ii][jj] : 2.0 * ionicMatrixMO[spin][ii][jj]);

    //                 if (printMOMatrix)
    //                     std::cout << std::right << std::setw(17) << __debug_MOMatrix[spin][ii][jj] << '\t';
    //             }
    //             if (printMOMatrix)
    //                 std::cout << std::endl;
    //         }

    //         if (printMOMatrix)
    //             std::cout << std::defaultfloat << "Total sum of MO matrix elements for " << to_string(static_cast<SpinType>(spin)) << " spin: " << std::setprecision(10) << __debug_totalSumMO[spin] << std::endl;
    //     }
    // }

    return toiMatrixMO;
}

double Becke::ionic_potential(int i, int j, SpinType spinType, const std::array<double, 3>& chargePosition, double charge, int kmax, int lebedev_order, int radial_grid_factor)
{
    double sum = 0.0;

    if (spinType == SpinType::ALPHA)
    {
        sum = multicenter_integration(&phiStarVionicStarPhi, i, j, SpinType::ALPHA, chargePosition, charge, kmax, lebedev_order, radial_grid_factor);
    }
    else if (spinType == SpinType::BETA)
    {
        sum = multicenter_integration(&phiStarVionicStarPhi, i, j, SpinType::BETA, chargePosition, charge, kmax, lebedev_order, radial_grid_factor);
    }
    else
    {
        sum = (multicenter_integration(&phiStarVionicStarPhi, i, j, SpinType::ALPHA, chargePosition, charge, kmax, lebedev_order, radial_grid_factor)
                + multicenter_integration(&phiStarVionicStarPhi, i, j, SpinType::BETA, chargePosition, charge, kmax, lebedev_order, radial_grid_factor));
    }

    return sum;
}

void Becke::chiAtomic(std::vector<std::vector<double>>& chiAtomic, int kmax, int lebedev_order, int radial_grid_factor)
{
    // Build Becke grid if not already done
    if(_multigrid==false)
    {
        multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        _multigrid=true;
    }


    // Get necessary data from molecule and orbitals
    int numberOfAtoms = _molecule.getNumberOfAtoms();

    int numberOfAo = _orbitals.get_numberOfAo();

    std::vector<CGTF> vcgtf = _orbitals.get_vcgtf();
    const std::vector<std::vector<std::vector<double>>>& coefficients = _orbitals.get_coefficients();

    const std::vector<std::vector<double>>& orbitalEnergies = _orbitals.get_orbitalEnergy();

    std::vector<std::vector<int>> occupiedOrbitalNumbers;
    std::vector<std::vector<int>> virtualOrbitalNumbers;
    _orbitals.getOccupiedAndVirtualOrbitalNumbers(occupiedOrbitalNumbers, virtualOrbitalNumbers);


    // Compute the overlap matrix in AO basis (integrated on each atom)
    std::vector<std::vector<std::vector<double>>> overlapMatrixesAO(numberOfAtoms, std::vector<std::vector<double>>(numberOfAo));

    for(int I = 0; I < numberOfAtoms; ++I)
    {
        for (int i = 0; i < numberOfAo; ++i)
        {
            overlapMatrixesAO[I][i].resize(i + 1, 0.0);

            for (int j = 0; j <= i; ++j)
            {
                double sum = 0.0;

                #ifdef ENABLE_OMP
                #pragma omp parallel for reduction(+: sum)
                #endif
                for(size_t J = 0; J < _grid_weights[I].size(); ++J)
                {
                    double x = _grid_points[I][J][0];
                    double y = _grid_points[I][J][1];
                    double z = _grid_points[I][J][2];

                    sum += _grid_weights[I][J] * _grid_volumes[I][J] * vcgtf[i].func(x, y, z) * vcgtf[j].func(x, y, z);
                }

                overlapMatrixesAO[I][i][j] = sum;
            }
        }
    }

    // Build chiAtomic triangular matrix
    chiAtomic.resize(numberOfAtoms, std::vector<double>());
    for (int i = 0; i < numberOfAtoms; ++i)
    {
        chiAtomic[i].resize(i + 1, 0.0);
    }


    int I, J;
    #ifdef ENABLE_OMP
    #pragma omp parallel for private(I, J)
    #endif
    for (I = 0; I < numberOfAtoms; ++I)
    {
        for (J = 0; J <= I; ++J)
        {
            for (int spin = 0; spin < 2; ++spin)
            {
                size_t nbOccupiedOrbitals = occupiedOrbitalNumbers[spin].size();
                size_t nbVirtualOrbitals = virtualOrbitalNumbers[spin].size();

                for (size_t i = 0; i < nbOccupiedOrbitals; ++i)
                {
                    int occupiedOrbitalIndex = occupiedOrbitalNumbers[spin][i] - 1; // Convert to 0-based index

                    for (size_t j = 0; j < nbVirtualOrbitals; ++j)
                    {
                        int virtualOrbitalIndex = virtualOrbitalNumbers[spin][j] - 1; // Convert to 0-based index

                        double sum_I = 0.0;
                        double sum_J = 0.0;

                        for (size_t m = 0; m < coefficients[spin][occupiedOrbitalIndex].size(); ++m)
                        {
                            for (size_t n = 0; n < coefficients[spin][virtualOrbitalIndex].size(); ++n)
                            {
                                double coeffProduct = coefficients[spin][occupiedOrbitalIndex][m] * coefficients[spin][virtualOrbitalIndex][n];

                                sum_I += coeffProduct * (n <= m ? overlapMatrixesAO[I][m][n] : overlapMatrixesAO[I][n][m]);
                                sum_J += coeffProduct * (n <= m ? overlapMatrixesAO[J][m][n] : overlapMatrixesAO[J][n][m]);
                            }
                        }

                        chiAtomic[I][J] += sum_I * sum_J / (orbitalEnergies[spin][virtualOrbitalIndex] - orbitalEnergies[spin][occupiedOrbitalIndex]);
                    }
                }
            }

            chiAtomic[I][J] *= (-2.0);
        }
    }
}

double Becke::density(const Orbitals& orbitals, double x, double y, double z)
{
	return orbitals.density(x, y, z);
}

double Becke::density(const Orbitals& orbitals, const std::vector<int>& orbitalNumbers, const std::vector<SpinType>& orbitalSpins, double x, double y, double z)
{
    if (orbitalNumbers.size() != orbitalSpins.size())
    {
        print_error("Error in Becke::density(): orbitalNumbers and orbitalSpins vectors must have the same size.");
        std::exit(1);
    }

    return orbitals.density(orbitalNumbers, orbitalSpins, x, y, z);
}

double Becke::prodGTF(const std::vector<GTF>& p, double x, double y, double z)
{
    std::vector<double> c (3);
    c[0]=x;
    c[1]=y;
    c[2]=z;
    return p*c;
}

double Becke::CGTFstarCGTF(const Orbitals& orbitals, int i, int j, double x, double y ,double z)
{
    double c = 0.0;

    const std::vector<CGTF>& vcgtf = orbitals.get_vcgtf();

    if(i == j)
    {
        c = vcgtf[i].func(x, y, z);
        c *= c;
    }
    else
    {
        c = vcgtf[i].func(x, y, z) * vcgtf[j].func(x, y, z);
    }

    return c;
}

double Becke::phi(const Orbitals& orbitals, int i, double x, double y, double z, SpinType spinType)
{
    double phi = 0.0;
    
    for(size_t j = 0; j < orbitals.get_vcgtf().size(); ++j)
    {
        phi += orbitals.get_coefficients()[static_cast<int>(spinType)][i][j] * orbitals.get_vcgtf()[i].func(x, y, z);
    }

    return phi;
}

double Becke::phiStarPhi(const Orbitals& orbitals, int i, int j, double x, double y, double z, SpinType spinType)
{
    double phi_star_phi = 0.0;
    int spin = static_cast<int>(spinType);
    const std::vector<std::vector<std::vector<double>>>& coefficients = orbitals.get_coefficients();
    const std::vector<CGTF>& vcgtf = orbitals.get_vcgtf();

    if (i == j)
    {
        double phi_i = 0.0;

        for (size_t k = 0; k < orbitals.get_vcgtf().size(); k++)
        {
            phi_i += coefficients[spin][i][k] * vcgtf[k].func(x, y, z);
        }

        phi_star_phi = phi_i * phi_i;
    }
    else
    {
        double phi_i = 0.0;
        double phi_j = 0.0;

        for (size_t k = 0; k < orbitals.get_vcgtf().size(); k++)
        {
            double v_k = vcgtf[k].func(x, y, z);

            phi_i += coefficients[spin][i][k] * v_k;
            phi_j += coefficients[spin][j][k] * v_k;
        }

        phi_star_phi = phi_i * phi_j;
    }

    return phi_star_phi;
}

double Becke::phiStarVionicStarPhi(Orbitals& orbitals, int i, int j, double x, double y, double z, SpinType spinType, const std::array<double, 3>& chargePosition, double charge)
{
    double distance = std::sqrt((chargePosition[0] - x) * (chargePosition[0] - x)
                                + (chargePosition[1] - y) * (chargePosition[1] - y)
                                + (chargePosition[2] - z) * (chargePosition[2] - z));

    return (distance > 1e-10 ? - charge * phiStarPhi(orbitals, i, j, x, y, z, spinType) / distance : 0.0);
}

double Becke::multicenter_integration(const Grid& g, int kmax, int lebedev_order, int radial_grid_factor)
{  
	if(_multigrid==false)
	{
		multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        	_multigrid=true;
	}

	int Nat = _molecule.getNumberOfAtoms();

	double integral = 0.0;

	for(int I=0; I<Nat; I++)
	{
		double integ=0.0;
		#ifdef ENABLE_OMP
		#pragma omp parallel for reduction(+:integ)
		#endif
		for(size_t J=0; J<_grid_weights[I].size(); J++)
		{        
			integ += _grid_volumes[I][J] * _grid_weights[I][J] * g.value(_grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2]);       // evaluate function on the grid
		}
		integral+=integ;
	}
    return integral;
}

std::vector<double> Becke::multicenter_sub_integration(const Grid& g,int kmax , int lebedev_order, int radial_grid_factor)
{   
	if(_multigrid==false)
	{
		multicenter_grids(kmax, lebedev_order, radial_grid_factor);
		_multigrid=true;
	}
	int Nat = _molecule.getNumberOfAtoms();
	std::vector<double> sub_integral(Nat,0.0) ;
	for(int I=0; I<Nat; I++)
	{
		double integ=0.0;
		#ifdef ENABLE_OMP
		#pragma omp parallel for reduction(+:integ)
		#endif
		for(size_t J=0; J<_grid_weights[I].size(); J++)
		{
			integ  += _grid_volumes[I][J] * _grid_weights[I][J]*g.value(_grid_points[I][J][0] , _grid_points[I][J][1], _grid_points[I][J][2]);       // evaluate function on the grid
		}
		sub_integral[I] = integ;
	}
	return sub_integral;
}

std::vector<double> Becke::PartialChargeAndEnergy(const Grid& g, int kmax, int lebedev_order, int radial_grid_factor)
{
    partial_charge(g,kmax, lebedev_order, radial_grid_factor);
    std::vector<double> c = _partial_charge;
    c.insert(c.begin(), _energy);
    return c;
}
std::vector<std::vector<double>> Becke::PartialChargesAndEnergy(int kmax, int lebedev_order, int radial_grid_factor)
{
    partial_charge(kmax, lebedev_order, radial_grid_factor);
    std::vector<std::vector<double>> c(2);
    c[0].push_back({_energy});
    c[1] = _partial_charge;
    
    return c;
}



void Becke::printCharges()
{
    std::cout << "Number of atoms = " << _molecule.getNumberOfAtoms() << std::endl;
    for(int i = 0; i < _molecule.getNumberOfAtoms(); i++)
    {
        std::cout << " Atom = " << std::left << std::setw(10) << _molecule.atom(i).get_symbol() << ", " << std::setw(10) << " value = " << std::setw(15) << _partial_charge[i] << std::endl;
    }
}
