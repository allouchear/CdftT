#include<cmath>
#include<cstdlib>
#include<functional>
#include<iomanip>
#include<iostream>
#include<vector>

#include <Cube/Grid.h>
#include <Becke/Becke.h>



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
    _energy = fchk.TotalEnergy();
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
    int Nat = _molecule.number_of_atoms();
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

    std::cout << "Number of Becke  grid points = " << printnpts << std::endl << std::endl;
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

    int Nat = _molecule.number_of_atoms();

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

double Becke::multicenter_integration(std::function<double(Orbitals&, int i, int j, double x, double y, double z)> f, int i, int j, int kmax, int lebedev_order, int radial_grid_factor)
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

    int Nat = _molecule.number_of_atoms();

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

double Becke::multicenter_integration(std::function<double(Orbitals& Orb, int i, int j, double x, double y, double z, SpinType spinType)> f, int i, int j, int kmax, int lebedev_order, int radial_grid_factor, SpinType spinType)
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

    int Nat = _molecule.number_of_atoms();

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

std::vector<double> Becke::multicenter_sub_integration(std::function<double(Orbitals& Orb, double x, double y, double z)> f, int kmax, int lebedev_order, int radial_grid_factor)
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

    int Nat = _molecule.number_of_atoms();

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
    int Nat=_molecule.number_of_atoms();
    std::vector<double> qn (Nat);
    std::vector<double> In = multicenter_sub_integration(&density, kmax, lebedev_order, radial_grid_factor);

    for(int i=0; i<Nat ;i++)
        qn[i]=_molecule.atom(i).get_atomicNumber() - In[i];

    _partial_charge=qn;
}

void Becke::partial_charge(const Grid &g, int kmax, int lebedev_order, int radial_grid_factor)
{
    int Nat = _molecule.number_of_atoms();
    std::vector<double> qn(Nat);
    std::vector<double> In = multicenter_sub_integration(g);
    for (int i = 0; i < Nat; i++)
    {
        qn[i] = _molecule.atom(i).get_atomicNumber() - In[i];
    }
    _partial_charge = qn;
}

double Becke::density(Orbitals& Orb, double x, double y, double z)
{
	return Orb.density(x,y,z);
}

double Becke::prodGTF(const std::vector<GTF>& p, double x, double y, double z)
{
    std::vector<double> c (3);
    c[0]=x;
    c[1]=y;
    c[2]=z;
    return p*c;
}

double Becke::CGTFstarCGTF(Orbitals& Orb, int i, int j, double x, double y ,double z)
{
    double c=0.0;

    if(i==j)
    {
        c=Orb.get_vcgtf()[i].func(x,y,z);
        return c*c;
    }

    else
        c=Orb.get_vcgtf()[i].func(x,y,z)*Orb.get_vcgtf()[j].func(x,y,z);

    return c;
}

double Becke::phi(Orbitals& orbitals, int i, double x, double y, double z, SpinType spinType)
{
    double phi = 0.0;
    
    for(size_t j = 0; j < orbitals.get_vcgtf().size(); ++j)
    {
        phi += orbitals.get_coefficients()[static_cast<int>(spinType)][i][j] * orbitals.get_vcgtf()[i].func(x, y, z);
    }

    return phi;
}

double Becke::phiStarPhi(Orbitals& orbitals, int i, int j, double x, double y, double z, SpinType spinType)
{
    double phi_star_phi = 0.0;

    if (i == j)
    {
        double phi_i = 0.0;

        for (size_t k = 0; k < orbitals.get_vcgtf().size(); k++)
        {
            phi_i += orbitals.get_coefficients()[static_cast<int>(spinType)][i][k] * orbitals.get_vcgtf()[k].func(x, y, z);
        }

        phi_star_phi = phi_i * phi_i;
    }
    else
    {
        double phi_i = 0.0;
        double phi_j = 0.0;

        for (size_t k = 0; k < orbitals.get_vcgtf().size(); k++)
        {
            phi_i += orbitals.get_coefficients()[static_cast<int>(spinType)][i][k] * orbitals.get_vcgtf()[k].func(x, y, z);
            phi_j += orbitals.get_coefficients()[static_cast<int>(spinType)][j][k] * orbitals.get_vcgtf()[k].func(x, y, z);
        }

        phi_star_phi = phi_i * phi_j;
    }

    return phi_star_phi;
}

double Becke::phiStarVionicStarPhi(Orbitals& orbitals, int i, int j, double x, double y, double z, SpinType spinType, const std::array<double, 3>& chargePosition, const double charge)
{
    return phiStarPhi(orbitals, i, j, x, y, z, spinType) * charge / std::sqrt((x - chargePosition[0]) * (x - chargePosition[0])
                                                                              + (y - chargePosition[1]) * (y - chargePosition[1])
                                                                              + (z - chargePosition[2]) * (z - chargePosition[2]));
}

double Becke::multicenter_integration(const Grid& g, int kmax, int lebedev_order, int radial_grid_factor)
{  
	if(_multigrid==false)
	{
		multicenter_grids(kmax, lebedev_order, radial_grid_factor);
        	_multigrid=true;
	}

	int Nat = _molecule.number_of_atoms();

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
	int Nat = _molecule.number_of_atoms();
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
std::vector<double> Becke::get_Partial_Charge()
{
	return _partial_charge;
}
double Becke::get_Energy()
{
	return _energy;
}
void Becke::printCharges()
{
    std::cout << "Number of atoms = " << _molecule.number_of_atoms() << std::endl;
    for(int i = 0; i < _molecule.number_of_atoms(); i++)
    {
        std::cout << " Atom = " << std::left << std::setw(10) << _molecule.atom(i).get_symbol() << ", " << std::setw(10) << " value = " << std::setw(15) << _partial_charge[i] << std::endl;
    }
}
