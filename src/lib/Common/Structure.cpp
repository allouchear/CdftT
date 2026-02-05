#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../Common/Structure.h"


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Structure::Structure():
    _atoms()
{ }

Structure::Structure(const std::vector<Atom>& A):
    _atoms(A)
{ }


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

const std::vector<Atom>& Structure::get_atoms() const
{
    return _atoms;
}


void Structure::read_From_Cube(ifstream& nameFile, int Natoms,const PeriodicTable& Table )
{
    double input;
    _atoms.resize(Natoms);
    for(int i=0; i<Natoms; i++)
    {    
        nameFile>>input;
        Atom a(Table, int(input));
        _atoms[i]=a;
        nameFile>>input;
        _atoms[i].set_charge(input);
        for(int j=0; j<3; j++)
        {
            nameFile>>input;
            _atoms[i].setCoordinate(j,input);
        }
    }
}

Structure::Structure(ifstream& nameFile, const int Natoms, const PeriodicTable& Table)
{
    read_From_Cube(nameFile, Natoms, Table);    
}

Structure Structure::operator+(const Structure& S) const
{
    if(S._atoms.size()>_atoms.size())
    {
        return S;
    }
    else
    {
        return *this;
    }
}

Structure Structure::add(const Structure& S)
{
    try
    {
        for(int j=0; j<S.number_of_atoms(); j++)
        {    
            bool b=true;
            for(int k=0; k<number_of_atoms(); k++)
            {
                if( _atoms[k].computeDistance(S._atoms[j])==0 )
                {
                    if(_atoms[k].get_symbol()==S._atoms[j].get_symbol())
                    {
                        throw std::string("::add(const Structure& S): can't add two identical atoms at same coordinates");
                    }
                    else
                    {
                        b=false;
                        break;
                    }
                }
            }
            if(b)
            {
                    cout<<"added"<<endl;
                    _atoms.push_back(S._atoms[j]);
            }
        }
        return *this;
    }
    catch(std::string error)
    {
        cout<<error<<endl;
        exit(1);
    }
}

void Structure::read_from_wfx(WFX& wfx, const PeriodicTable& Table)
{
    int n=0;
    _atoms.resize(wfx.Number_of_Nuclei());
    for(int i=0; i<wfx.Number_of_Nuclei(); i++)
    {
        _atoms[i]=Atom(Table, wfx.Atomic_Number()[i]);
        _atoms[i].set_charge(wfx.Nuclear_Charges()[i]);
        n=3*i;
        for(int j=0; j<3; j++)
        {
            _atoms[i].setCoordinate(j, wfx.Nuclear_Cartesian_Coordinates()[n]);
            n++;
        }
    }
}

Structure::Structure(WFX& wfx, const PeriodicTable& Table)
{
    read_from_wfx(wfx, Table);
}

void Structure::read_from_fchk(FCHK& fchk, const PeriodicTable& Table)
{
    int n=0;
    _atoms.resize(fchk.NumberOfAtoms());
    for(int i=0; i<fchk.NumberOfAtoms(); i++)
    {
        _atoms[i]=Atom(Table, fchk.AtomicNumbers()[i]);
        _atoms[i].set_charge(fchk.NuclearCharges()[i]);
        n=3*i;
        for(int j=0; j<3; j++)
        {
            _atoms[i].setCoordinate(j, fchk.CurrentCartesianCoordinates()[n]);
            n++;
        }
    }
}

Structure::Structure(FCHK& fchk, const PeriodicTable& Table)
{
    read_from_fchk(fchk, Table);
}

void Structure::read_from_moldengab(MOLDENGAB& moldengab, const PeriodicTable& Table)
{
    _atoms.resize(moldengab.NumberOfAtoms());
    for(int i=0; i<moldengab.NumberOfAtoms(); i++)
    {
        _atoms[i]=Atom(Table, moldengab.AtomicNumbers()[i]);
        for(int j=0; j<3; j++)
            _atoms[i].setCoordinate(j, moldengab.Coordinates()[i][j]);
    }
}

Structure::Structure(MOLDENGAB& moldengab, const PeriodicTable& Table)
{
    read_from_moldengab(moldengab, Table);
}

void Structure::read_from_log(LOG& log, const PeriodicTable& Table)
{
    _atoms.resize(log.NumberOfAtoms());
    for(int i=0; i<log.NumberOfAtoms(); i++)
    {
        _atoms[i]=Atom(Table, log.AtomicNumbers()[i]);
        for(int j=0; j<3; j++)
            _atoms[i].setCoordinate(j, log.Coordinates()[i][j]);
    }
}

Structure::Structure(LOG& log, const PeriodicTable& Table)
{
    read_from_log(log, Table);
}