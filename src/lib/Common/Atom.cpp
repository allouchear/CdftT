#include <cmath>
using namespace std;
#include <iostream>
#include <Common/Atom.h>
#include <Common/Constants.h>

/*constructeur par défaut: initialise les attributs cinématique de l'atome à 0*/
Atom::Atom()
{
    for(int i=0; i<3;i++)
    {
        _coordinates[i]=0;
        _gradient[i]=0;
        _velocity[i]=0;
    }
    _name = "none";
    _symbol = "none";
    _atomic_number = 0;
    _charge = 0;
    _charge_0 = 0;
    _hardness = 0;
    _width = 0;
    _covalent_radii=0;
    _e=Element();
}

//initialize par le nom/symbole de l'élément
Atom::Atom(const PeriodicTable& Table, const string& name)
{    
    for(int i=0; i<3;i++)
    {
        _coordinates[i]=0;
        _gradient[i]=0;
        _velocity[i]=0;
    }
    _e = Table.element(name);
    _name = _e.name();
    _symbol = _e.symbol();
    _atomic_number = _e.atomic_number();
    _covalent_radii=_e.covalent_radii();
}

//initialise par le numéro atomique
Atom::Atom(const PeriodicTable& Table, const int Z)
{    
    for(int i=0; i<3;i++)
    {
        _coordinates[i]=0;
        _gradient[i]=0;
        _velocity[i]=0;
    }
    _e = Table.element(Z);
    _name = _e.name();
    _symbol = _e.symbol();
    _atomic_number = _e.atomic_number();
    _covalent_radii=_e.covalent_radii();
}

Atom::~Atom(){}

const double* Atom::coordinates() const
{
    return _coordinates;
}

double* Atom::gradient() 
{
    return _gradient;
}

double* Atom::velocity() 
{
    return _velocity;
}

string Atom::name() const
{
    return _name;
}

string Atom::symbol() const
{
    return _symbol;
}

int Atom::atomic_number() const
{
    return _atomic_number;
}

double Atom::charge() const
{
    return _charge;
}

double Atom::charge_0() const
{
    return _charge_0;
}

double Atom::hardness() const
{
    return _hardness;
}

double Atom::width() const
{
    return _width;
}

double Atom::covalent_radii() const
{
    return _covalent_radii;
}

Element Atom::element()
{
    return _e;
}

void Atom::set_coordinates(const int i, const double d )
{
    _coordinates[i] = d;
}

void Atom::set_gradient(const int i, const double d )
{
    _gradient[i] = d;
}

void Atom::set_velocity(const int i, const double d )
{
    _velocity[i] = d;
}

void Atom::set_charge(const double c )
{
    _charge = c;
}



//renvoi la distance entre l'atome et un autre(a2)
double Atom::get_distance(const Atom& a2) const
{
    double C1[3];
    double C2[3];
    for(int i=0;i<3;i++)
    {
        C1[i] =_coordinates[i];
        C2[i] = a2._coordinates[i];
    }
    double x, y, z;
    
        x = C1[ 0 ] - C2[ 0 ];
           y = C1[ 1 ] - C2[ 1 ];
           z = C1[ 2 ] - C2[ 2 ];

    return sqrt( x * x + y * y + z * z );
}

//renvoi langle formée par l'atome et 2 autres(a2,a3)
double Atom::get_angle(const Atom& a2,const Atom& a3) const
{
    double C1[3];
    double C2[3];
    double C3[3];
    for(int i=0;i<3;i++)
    {
        C1[i] = _coordinates[i];
        C2[i] = a2._coordinates[i];
        C3[i] = a3._coordinates[i];
    }
    double x12, x32, y12, y32, z12, z32, l12, l32, dp;
    
      x12 = C1[ 0 ] - C2[ 0 ];
      y12 = C1[ 1 ] - C2[ 1 ];
      z12 = C1[ 2 ] - C2[ 2 ];
      x32 = C3[ 0 ] - C2[ 0 ];
      y32 = C3[ 1 ] - C2[ 1 ];
      z32 = C3[ 2 ] - C2[ 2 ];

      l12 = sqrt( x12 * x12 + y12 * y12 + z12 * z12 );
      l32 = sqrt( x32 * x32 + y32 * y32 + z32 * z32 );
      if( l12 == 0.0 )
    {
                   return 0.0;
      }
      if( l32 == 0.0 )
    {
                   return 0.0;
      }
    dp = ( x12 * x32 + y12 * y32 + z12 * z32 ) / (l12 * l32 );
    if ( dp < -1.0 )
        dp = -1.0;
    else if ( dp > 1.0 )
        dp = 1.0;
    return RADTODEG * acos(dp);
}


double Atom::get_torsion(const Atom& a2, const Atom& a3, const Atom& a4) const
{
    double C1[3];
    double C2[3];
    double C3[3];
    double C4[3];
    for(int i=0;i<3;i++)
    {
        C1[i] = _coordinates[i];
        C2[i] = a2._coordinates[i];
        C3[i] = a3._coordinates[i];
        C4[i] = a4._coordinates[i];
    }

    double   xij, yij, zij;
           double   xkj, ykj, zkj;
    double   xkl, ykl, zkl;
         double   dx, dy, dz;
    double   gx, gy, gz;
    double   bi, bk;
    double   ct, d, ap, app, bibk;

    xij = C1[ 0 ] - C2[ 0 ];
    yij = C1[ 1 ] - C2[ 1 ];
    zij = C1[ 2 ] - C2[ 2 ];
    xkj = C3[ 0 ] - C2[ 0 ];
    ykj = C3[ 1 ] - C2[ 1 ];
    zkj = C3[ 2 ] - C2[ 2 ];
    xkl = C3[ 0 ] - C4[ 0 ];
    ykl = C3[ 1 ] - C4[ 1 ];
    zkl = C3[ 2 ] - C4[ 2 ];

    dx = yij * zkj - zij * ykj;
    dy = zij * xkj - xij * zkj;
    dz = xij * ykj - yij * xkj;
    gx = zkj * ykl - ykj * zkl;
    gy = xkj * zkl - zkj * xkl;
    gz = ykj * xkl - xkj * ykl;

    bi = dx * dx + dy * dy + dz * dz;
    bk = gx * gx + gy * gy + gz * gz;
    ct = dx * gx + dy * gy + dz * gz;
    bibk = bi * bk;
    
    if ( bibk < 1.0e-6 )    
        return 0;
  
    ct = ct / sqrt( bibk );
  
    if( ct < -1.0 )
          ct = -1.0;
    else if( ct > 1.0 )
          ct = 1.0;

    ap = acos( ct );
  
    d  = xkj*(dz*gy-dy*gz) + ykj*(dx*gz-dz*gx) + zkj*(dy*gx-dx*gy);
  
    if( d < 0.0 )
          ap = -ap;
  
    ap = PI - ap;
    app = 180.0 * ap / PI;
    if( app > 180.0 )
             app = app - 360.0;
    return( app );
}

