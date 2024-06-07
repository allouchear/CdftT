#include<iostream>
#include<cmath>
#include<analytic/Becke/Becke.h>

using namespace std;

Becke::Becke()
{
	_molecule=Structure();
	_grid=GridPoints();
    _grid_points=vector<vector<vector<double>>> ();
    _grid_weights=vector<vector<double>> ();
    _grid_volumes=vector<vector<double>> ();
}


//A revoir
Becke::Becke(WFX& wfx, const PeriodicTable& Table)
{
	_molecule=Structure(wfx, Table);
	_grid=GridPoints();
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
	int n;
	vector<int> L=_grid.Lebedev_Lmax();
	for(size_t i=0; i<_grid.Lebedev_Lmax().size(); i++)
		L[i]=abs(L[i]-lebedev_order);
    
    int Ln=L[0];  
	n=0;

	for(size_t i=1; i<_grid.Lebedev_Lmax().size(); i++)
		if(Ln>L[i])
        {
			Ln=L[i];
            n=i;
        }

	int Lmax = _grid.Lebedev_Lmax()[n];
	if(lebedev_order != Lmax)
	{
		cout<<"Error, no grid found."<<endl;
		exit(1);
	}

	return GridPoints(Lmax);
}

double Becke::s(double mu, int k=3)
{
    // for nuclear weight functions
    double f = mu;
    for(int ik=0; ik<k; ik++)
      	f = 1.5 * f -0.5 * f*f*f;
    return 0.5*(1-f);
}

void Becke::multicenter_grids(int kmax=3, int lebedev_order=23, int radial_grid_factor=1)
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
    //cout<<endl<<"Check Point 1"<<endl;
    int Nang = Angular.Npts();
    vector<double> sc (Nang);
    vector<double> thang (Nang);
    vector<double> phiang (Nang);
    vector<double> ss (Nang);
    vector<double> c (Nang);
    vector<double> wang (Nang);

    for(int i=0; i<Nang; i++)
    {
        thang[i]=Angular.LebedevGridPoints()[i][0];
        phiang[i]=Angular.LebedevGridPoints()[i][1];
        wang[i] = Angular.LebedevGridPoints()[i][2];
    	sc[i] = sin(thang[i])*cos(phiang[i]);
        //cout<<endl<<"sc[i] = "<<sc[i]<<endl;
        //cout<<endl<<"Angular.LebedevGridPoints()[i][0] = "<<Angular.LebedevGridPoints()[i][0]<<endl;
        //cout<<"Angular.LebedevGridPoints()[i][1] = "<<Angular.LebedevGridPoints()[i][1]<<endl;
    	ss[i] = sin(thang[i])*sin(phiang[i]);
    	c[i]  = cos(thang[i]);
    }

    // declaration of all variables
    int Nat = _molecule.number_of_atoms();
    double chi, uij;
    vector<double> wr(Nat, 0.0);
    vector<vector<double>> R (Nat,vector<double> (Nat,0));     // distances between atoms i and j
    vector<vector<double>> a (Nat,vector<double> (Nat,0));     // scaling factor used in eqn. A2
    vector<vector<vector<double>>> list_coord_I (Nat);
    vector<vector<vector<double>>> grid_points;
    vector<vector<double>> grid_weights;
    vector<vector<double>> grid_volumes;
    int Nr, Npts;
    double rm, mu, nu;
    double Ptot=0;

    for(int i=0; i<Nat; i++)
        for(int j=i+1; j<Nat; j++)
        {
            R[i][j] = R[j][i] = _molecule.atoms(i).get_distance(_molecule.atoms(j));
            //for(int k=0; k<3; k++)
                //cout<<endl<<"Atom 1 | Atom 2 : "<<_molecule.atoms(i).coordinates()[k]<<" | "<<_molecule.atoms(j).coordinates()[k]<<endl;
            //cout<<endl<<"R[i][j] = "<<R[i][j]<<endl;
            // ratio of Slater radii
            chi = _molecule.atoms(i).covalent_radii() / _molecule.atoms(j).covalent_radii();
            uij = (chi-1)/(chi+1);
            a[i][j] = uij/(uij*uij - 1);
            a[j][i] = -a[i][j];
        }
    
    //cout<<"Check Point 3"<<endl;
    // atom-centered subintegral
    for(int I=0; I<Nat; I++)
    {
        //cout<<"Check Point Boucle I"<<endl;
        // radial grid
        Nr = number_of_radial_points(_molecule.atoms(I).atomic_number());
        // increase number of grid points is requested
        Nr *= radial_grid_factor;
        rm = 0.5*_molecule.atoms(I).covalent_radii();
        Npts = Nr*Nang;

        vector<double> k (Nr);
        vector<double> xr (Nr);
        vector<double> radial_weights (Nr);
        vector<double> g (Nr);
        vector<double> r (Nr);
        vector<double> x(Npts), y(Npts) , z(Npts), weights(Npts);
        vector<vector<double>> dist (Npts, vector<double> (Nat,0.0));
        vector<vector<double>> P (Npts, vector<double> (Nat,1));

        //cout<<"Check Point Boucle I 2"<<endl;
        for(int J=0; J<Nr; J++)
        {
            //cout<<"Check Point Boucle J"<<endl;
            k[J]=J+1;
            // grid points on interval [-1,1]
            xr[J] = cos(k[J]/(Nr+1.0) * M_PI);
            //cout<<endl<<"xr[J] = "<<xr[J]<<endl;
            // weights
            radial_weights[J] = (M_PI/(Nr+1.0) * sin(k[J]/(Nr+1.0) * M_PI))*(M_PI/(Nr+1.0) * sin(k[J]/(Nr+1.0) * M_PI));
            // from variable transformation
            g[J] = 2 * rm*rm*rm * sqrt(((1+xr[J])/((1-xr[J])*(1-xr[J])*(1-xr[J])))*((1+xr[J])/((1-xr[J])*(1-xr[J])*(1-xr[J])))*((1+xr[J])/((1-xr[J])*(1-xr[J])*(1-xr[J]))));
            radial_weights[J] *= g[J];
            // radial grid points on interval [0,infinity]
            r[J] = rm * (1+xr[J])/(1-xr[J]);
            //cout<<"Check Point Boucle J 2"<<endl;
        }
        
        int n=0;
            // cartesian coordinates of grid
        for (int i=0; i<Nr; i++)
            for(int j=0; j<Nang; j++)
            {
                //cout<<"Check Point Boucle ij"<<endl;
                x[n] = r[i] * sc[j] + _molecule.atoms(I).coordinates()[0];
                //cout<<endl<<"x[n] = "<<x[n]<<endl;
                //cout<<"Check Point Boucle ij 2"<<endl;
                y[n] = r[i] * ss[j] + _molecule.atoms(I).coordinates()[1];
                //cout<<"Check Point Boucle ij 3"<<endl;
                z[n] = r[i] * c[j] + _molecule.atoms(I).coordinates()[2];
                //cout<<"Check Point Boucle ij 4"<<endl;
                //cout<<"weights[i][j] = "<<weights[i][j]<<endl;
                //cout<<"radial_weights[I]"<<radial_weights[I]<<endl;
                //cout<<"Test access 5 "<<Angular.LebedevGridPoints()[j][2]<<endl;
                weights[n] = radial_weights[i] * 4.0 * M_PI * wang[j];
                //Angular.LebedevGridPoints()[j][2];
                //cout<<"Check Point Boucle ij 5"<<endl;
                n++;
            }

        //cout<<"Check Point 4"<<endl;

        // distance between grid points and atom i
        for (int i=0; i<Npts; i++)
            for(int j=0; j<Nat; j++)
            {
                dist[i][j] = sqrt((x[i] - _molecule.atoms(j).coordinates()[0])*(x[i] - _molecule.atoms(j).coordinates()[0]) + (y[i] - _molecule.atoms(j).coordinates()[1])*(y[i] - _molecule.atoms(j).coordinates()[1]) + (z[i] - _molecule.atoms(j).coordinates()[2])*(z[i] - _molecule.atoms(j).coordinates()[2]) );
                //cout<<endl<<"dist[i][j] = "<<dist[i][j]<<endl;
            }


        // P_i(r) as defined in eqn. (13)
        
        //cout<<"Check Point 5"<<endl;

        for(int i=0; i<Nat; i++)
            for(int j=0; j<Nat; j++)
            {
                // mu_ij as defined in eqn. (11)
                if(i!=j)
                    for(int k=0; k<Npts; k++)
                    {
                        mu = (dist[k][i]-dist[k][j])/R[i][j];
                        //cout<<endl<<"mu = "<<mu<<endl;
                        nu = mu + a[i][j]*(1-mu*mu);
                        P[k][i] *= s(nu);
                        //cout<<endl<<"P[k][I] = "<<P[k][I]<<endl;
                    }
            }

        //cout<<"Check Point 6"<<endl;

        for(int i=0; i<Npts; i++)
            for(int j=0; j<Nat; j++)
            {
    	       Ptot += P[i][j];
            //cout<<endl<<"Ptot[i] = "<<Ptot[i]<<endl; 
            }
    
        //cout<<"Check Point 7"<<endl;

        // weight function due to partitioning of volume
        for(int i=0; i<Npts; i++)
        {
            wr[I] += P[i][I]/Ptot;
                //cout<<endl<<"wr[i] = "<<wr[i]<<endl;
        }
        
        //cout<<"Check Point 8"<<endl;
        vector<vector<double>> coord_i (Npts);
        for(int i=0; i<Npts; i++)
            coord_i[i]={x[i], y[i], z[i]};

        list_coord_I[I]=coord_i;
        //cout<<"Test push_back 2"<<endl;
        // The weights come from the fuzzy Voronoi partitioning 
        grid_weights.push_back( wr );
        // The naming is a little bit confusing, the `weights` are
        // actually the volume elements dV_i around each point.
        //cout<<"Test push_back 3"<<endl;
        grid_volumes.push_back( weights );
        //cout<<"Test push_back 4"<<endl;
        //cout<<"Check Point 9"<<endl;
    }

    double test=0;
    for(int i=0; i<Nat; i++)
        test+=wr[i];
    cout<<endl<<"test sum wr_n = 1 : "<<test<<endl;

    grid_points=list_coord_I;

    //cout<<"Check Point 10"<<endl;

	_grid_points=grid_points;
    _grid_weights=grid_weights;
    _grid_volumes=grid_volumes;

    //cout<<"Check Point 11"<<endl;
}


vector<vector<double>> Becke::join_grids()
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
    vector<double> w;
    // sampling points of quadrature rule
    vector<double> x, y, z;
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

    vector<vector<double>> join_grid;
    join_grid.push_back(x);
    join_grid.push_back(y);
    join_grid.push_back(z);
    join_grid.push_back(w);

    return join_grid;
}

double Becke::multicenter_integration(function<double(vector<GTF> p, double x, double y, double z)> f, vector<GTF> p, double x=0, double y=0, double z=0, int lebedev_order=23, int radial_grid_factor=1)
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
    multicenter_grids();

    int Nat = _molecule.number_of_atoms();

    double integral = 0.0;
    double sub_integral =0.0;
    vector<double> fI (Nat, 0.0);


    // atom-centered subintegral
    for(int I=0; I<Nat; I++)
    {
        for(size_t J=0; J<_grid_weights[I].size(); J++)
        {
            // evaluate function on the grid
            fI[I] += _grid_weights[I][J] * f(p, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2]);
            //cout<<endl<<"f(p, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2]) = "<<f(p, _grid_points[I][J][0], _grid_points[I][J][1], _grid_points[I][J][2])<<endl;
            sub_integral += _grid_volumes[I][J] * fI[J];
        }
    integral += sub_integral;
    }

    return integral;
}

double Becke::overlap(GTF A, GTF B)
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
    vector<GTF> p (2);
    p[0]=A;
    p[1]=B;

    double Sab = multicenter_integration(&prod, p);

    return Sab;
}

double Becke::prod(vector<GTF> p, double x, double y, double z)
{
    vector<double> c (3);
    c[0]=x;
    c[1]=y;
    c[2]=z;
    return p*c;
}