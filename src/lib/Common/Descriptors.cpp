#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Common/Descriptors.h>
#include <Common/PeriodicTable.h>
#include <Becke/Becke.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Descriptors::Descriptors()
{
    Structure str;
    _str=str;
    _okCharge=false;
    _mu=0.0;
    _mup=0.0;
    _mum=0.0;
    _xi=0.0;
    _hardness=0.0;
    _w=0.0;
    _wp=0.0;
    _wm=0.0;
    _S=0.0;
    _Qmax=0.0;
    _DEmin=0.0;
    reset();
}

Descriptors::Descriptors(const Structure& S, std::vector<double> Q0, std::vector<double> Qm, std::vector<double> Qp, double I, double A )
{
    _okCharge=true;
    _str=S;
    reset();
    std::vector<double> E(3,0);
    double i=0;
    double a=0;
    sortCharges(Q0, Qm, Qp, E, i, a);
    compute_All_From_Charge(I,A);
}

Descriptors::Descriptors(const Structure& S, std::vector<double> Q0, std::vector<double> Qm, std::vector<double> Qp, std::vector<double> E)
{
    compute_All_From_Charges(S, Q0, Qm, Qp, E);
}

Descriptors::Descriptors(const Grid &AIM1, const Grid &AIM2, const Grid &AIM3, double I, double A, PartitionMethod partitionMethod)
{
    _str=AIM1.get_structure();
    reset();
    compute_All_From_Grid(AIM1, AIM2, AIM3, I, A, partitionMethod);
}

Descriptors::Descriptors(std::ifstream &file0, std::ifstream &fileM, std::ifstream &fileP, double I, double A, PartitionMethod partitionMethod)
{
    compute_All_From_Cube(file0, fileM, fileP, I, A, partitionMethod);
}

Descriptors::Descriptors(std::ifstream &file0, std::ifstream &fileM, std::ifstream &fileP, std::vector<double> E, PartitionMethod partitionMethod)
{
    compute_All_From_Cube(file0, fileM, fileP, E, partitionMethod);

}

Descriptors::Descriptors(FCHK& fchk, const PeriodicTable& Table)
{
    _okCharge=false;
    _str=Structure(fchk, Table);
    _deltafk.resize(_str.getNumberOfAtoms());
    _wkm.resize(_str.getNumberOfAtoms());
    _wkp.resize(_str.getNumberOfAtoms());
    _Skm.resize(_str.getNumberOfAtoms());
    _Skp.resize(_str.getNumberOfAtoms());
    _Skfrac.resize(_str.getNumberOfAtoms());
    _hardnessk.resize(_str.getNumberOfAtoms());
    _hardnesskm.resize(_str.getNumberOfAtoms());
    _hardnesskp.resize(_str.getNumberOfAtoms());

    _fk0.resize(_str.getNumberOfAtoms());
}

Descriptors::Descriptors(LOG& log, const PeriodicTable& Table)
{
    _okCharge=false;
    _str=Structure(log, Table);
    _deltafk.resize(_str.getNumberOfAtoms());
    _wkm.resize(_str.getNumberOfAtoms());
    _wkp.resize(_str.getNumberOfAtoms());
    _Skm.resize(_str.getNumberOfAtoms());
    _Skp.resize(_str.getNumberOfAtoms());
    _Skfrac.resize(_str.getNumberOfAtoms());
    _hardnessk.resize(_str.getNumberOfAtoms());
    _hardnesskm.resize(_str.getNumberOfAtoms());
    _hardnesskp.resize(_str.getNumberOfAtoms());

    _fk0.resize(_str.getNumberOfAtoms());
}

Descriptors::Descriptors(MOLDENGAB& moldengab, const PeriodicTable& Table)
{
    _okCharge=false;
    _str=Structure(moldengab, Table);
    _deltafk.resize(_str.getNumberOfAtoms());
    _wkm.resize(_str.getNumberOfAtoms());
    _wkp.resize(_str.getNumberOfAtoms());
    _Skm.resize(_str.getNumberOfAtoms());
    _Skp.resize(_str.getNumberOfAtoms());
    _Skfrac.resize(_str.getNumberOfAtoms());
    _hardnessk.resize(_str.getNumberOfAtoms());
    _hardnesskm.resize(_str.getNumberOfAtoms());
    _hardnesskp.resize(_str.getNumberOfAtoms());

    _fk0.resize(_str.getNumberOfAtoms());
}

Descriptors::Descriptors(WFX& wfx, const PeriodicTable& Table)
{
    _okCharge=false;
    _str=Structure(wfx, Table);
    _deltafk.resize(_str.getNumberOfAtoms());
    _wkm.resize(_str.getNumberOfAtoms());
    _wkp.resize(_str.getNumberOfAtoms());
    _Skm.resize(_str.getNumberOfAtoms());
    _Skp.resize(_str.getNumberOfAtoms());
    _Skfrac.resize(_str.getNumberOfAtoms());
    _hardnessk.resize(_str.getNumberOfAtoms());
    _hardnesskm.resize(_str.getNumberOfAtoms());
    _hardnesskp.resize(_str.getNumberOfAtoms());

    _fk0.resize(_str.getNumberOfAtoms());
}


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

const std::vector<double>& Descriptors::get_deltafk() const
{
    return _deltafk;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void Descriptors::compute_All_From_Charge(double I, double A )
{
    compute_fk_From_Charge();
    set_all_mu(I, A);
    compute_all();
}

std::vector<double> Descriptors::compute_Charges_From_Becke(const Grid& grid)
{
    Becke B(grid);
    B.partial_charge(grid);
    std::vector<double> Bint=B.get_partial_charge();
    return B.get_partial_charge();
}

std::vector<double> Descriptors::compute_Charges_From_Grid(const Grid& AIM, PartitionMethod partitionMethod)
{
    GridCP gridcp;
    gridcp.buildBasins(AIM, partitionMethod);

    return gridcp.computeAIMCharges(AIM);
}

std::vector<double> Descriptors::compute_Charges_From_File(std::ifstream& file, PartitionMethod partitionMethod)
{
    PeriodicTable Table;
    Grid AIM(file, Table);
    _str=AIM.get_structure();
    reset();
    std::vector<double> Q;
    if(partitionMethod == PartitionMethod::BECKE)
    {
        Q=compute_Charges_From_Becke(AIM);
    }
    else
    {
        GridCP gridcp;
        gridcp.buildBasins(AIM, partitionMethod);
        Q=gridcp.computeAIMCharges(AIM);
    }
    return Q;
}

void Descriptors::compute_All_From_Grid(const Grid& AIM1,const Grid& AIM2,const Grid& AIM3 , double I, double A, PartitionMethod partitionMethod )
{
    _okCharge=true;
    _str=AIM1.get_structure();
    reset();
    std::vector<double> Q1 = compute_Charges_From_Grid(AIM1, partitionMethod );
    std::vector<double> Q2 = compute_Charges_From_Grid(AIM2, partitionMethod );
    std::vector<double> Q3 = compute_Charges_From_Grid(AIM3, partitionMethod );
    std::vector<double> E(3,0);
    double i=0;
    double a=0;
    sortCharges(Q1, Q2, Q3, E, i, a);
    compute_All_From_Charge(I,A);
}



void Descriptors::compute_All_From_Cube(std::ifstream& file1, std::ifstream& file2, std::ifstream& file3, double I, double A, PartitionMethod partitionMethod )
{
    _okCharge=true;
    std::vector<double> Q1 = compute_Charges_From_File(file1, partitionMethod);
    std::vector<double> Q2 = compute_Charges_From_File(file2, partitionMethod);
    std::vector<double> Q3 = compute_Charges_From_File(file3, partitionMethod);
    std::vector<double> E(3,0);
    double i=0;
    double a=0;
    sortCharges(Q1, Q2, Q3, E, i, a);
    compute_All_From_Charge(I,A);
}



void Descriptors::set_all_mu(const double& I,const double& A)
{
    _mum = -I;
    _mup = A;
    _mu = 0.5*(_mup+_mum);
}

void Descriptors::compute_all()
{
    _xi = -_mu;
    _hardness = _mup-_mum;
    _w = (_mu*_mu*0.5)/_hardness;
    _wp = (_mup*_mup*0.5)/_hardness;
    _wm = (_mum*_mum*0.5)/_hardness;
    _S = 1/_hardness;
    _Qmax = -_mu/_hardness;
    _DEmin = -0.5*(_mu*_mu)/_hardness;

    for(int i=0;i<_str.getNumberOfAtoms();i++)
    {
        _deltafk[i] = _fkp[i]-_fkm[i];
        _fk0[i]=0.5*(_fkm[i]+_fkp[i]);
        _wkm[i] = _fkm[i]*_w;
        _wkp[i] = _fkp[i]*_w;
        _Skm[i] = _fkm[i]*_S;
        _Skp[i] = _fkp[i]*_S;
        _Skfrac[i] = _Skm[i]/_Skp[i];
        _hardnessk[i] = _mup*_fkp[i]-_mum*_fkm[i];
        _hardnesskm[i] = _hardnessk[i]-((_mup-_mum)*(_fkp[i]-_fkm[i]));
        _hardnesskp[i] = _hardnessk[i]+((_mup-_mum)*(_fkp[i]-_fkm[i]));
    }
}

void Descriptors::compute_fk_From_Charge()
{    
    if(int(_Qm.size())!=_str.getNumberOfAtoms() or int(_Qp.size())!=_str.getNumberOfAtoms() or int(_Q0.size())!=_str.getNumberOfAtoms())
    {
        print_error("Error in Descriptors::compute_fk_From_Charge(): the number of atoms in _str is inconsistent with std::vector sizes... Please check the vectors.");
        
        std::exit(1);
    }
    else
    {
        for(int i=0; i<_str.getNumberOfAtoms();i++)
        {
            _fkp[i] = _Q0[i]-_Qm[i];
            _fkm[i] = _Qp[i]-_Q0[i];
        }
    }
}

void Descriptors::compute_fk()
{
    int atomNb = _str.getNumberOfAtoms();

    _fk0.resize(atomNb);
    _fkp.resize(atomNb);
    _fkm.resize(atomNb);

    for(int i = 0; i < atomNb ; ++i)
    {
        _fkm[i] = _Q0[i] - _Qp[i];
        _fkp[i] = _Qm[i] - _Q0[i];

        _fk0[i] = 0.5 * (_Qm[i] - _Qp[i]);
    }
}

void Descriptors::set_mu_fk_data(const std::vector<std::vector<double>>& f, double homoEnergy, double lumoEnergy)
{
    _mum = homoEnergy;
    _mup = lumoEnergy;

    _mu = 0.5 * (_mup + _mum);

    _fkm = f[0];
    _fkp = f[1];
}

void Descriptors::set_mu_fk_data(const std::vector<std::vector<double>>& data)
{
    int atomNb = _str.getNumberOfAtoms();

    _fkm.resize(atomNb);
    _fkp.resize(atomNb);
    _fk0.resize(atomNb);

    _mum = - (data[2][0] - data[0][0]);
    _mup = - (data[0][0] - data[1][0]);
    _mu = 0.5 * (_mup + _mum);

    for(int i = 0; i < atomNb; ++i)
    {
        _fkm[i] = data[0][i + 1] - data[1][i + 1];
        _fkp[i] = data[2][i + 1] - data[0][i + 1];
        _fk0[i] = 0.5 * (data[2][i + 1] - data[1][i + 1]);
    }

}



void Descriptors::reset()
{
    if (_str.getNumberOfAtoms()<1 )
    {
        _Q0 = std::vector<double>();
        _Qm = std::vector<double>();
        _Qp = std::vector<double>();
        _fk0 = std::vector<double>();
        _fkm = std::vector<double>();
        _fkp = std::vector<double>();
        _deltafk = std::vector<double>();
        _wkm = std::vector<double>();
        _wkp = std::vector<double>();
        _Skm = std::vector<double>();
        _Skp = std::vector<double>();
        _Skfrac = std::vector<double>();
        _hardnessk = std::vector<double>();
        _hardnesskm = std::vector<double>();
        _hardnesskp = std::vector<double>();
    }
    else 
    {
        _Q0 = std::vector<double>(_str.getNumberOfAtoms());
        _Qm = std::vector<double>(_str.getNumberOfAtoms());
        _Qp = std::vector<double>(_str.getNumberOfAtoms());
        _fk0 = std::vector<double>(_str.getNumberOfAtoms());
        _fkm = std::vector<double>(_str.getNumberOfAtoms());
        _fkp = std::vector<double>(_str.getNumberOfAtoms());
        _deltafk = std::vector<double>(_str.getNumberOfAtoms());
        _wkm = std::vector<double>(_str.getNumberOfAtoms());
        _wkp = std::vector<double>(_str.getNumberOfAtoms());
        _Skm = std::vector<double>(_str.getNumberOfAtoms());
        _Skp = std::vector<double>(_str.getNumberOfAtoms());
        _Skfrac = std::vector<double>(_str.getNumberOfAtoms());
        _hardnessk = std::vector<double>(_str.getNumberOfAtoms());
        _hardnesskm = std::vector<double>(_str.getNumberOfAtoms());
        _hardnesskp = std::vector<double>(_str.getNumberOfAtoms());
    }
}





/********************************************************************************************/

std::ostream& operator<<(std::ostream& flux, const Descriptors& desc)
{
    double HeV=27.21138469;
    
    flux<<std::scientific;
    flux<<std::setprecision(6);
    flux<<std::setw(15);
    if(desc._okCharge)
    {
        
        flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
        flux<<std::left<<std::setw(7)<<"Symbol"<<std::setw(4)<<"k"<<std::setw(15)<<std::right<<"Qk-"<<std::setw(15)<<std::right<<"Qk+"<<std::setw(15)<<std::right<<"Qk0"<<std::endl;
        for(int i=0; i<desc._str.getNumberOfAtoms(); i++)
            flux<<std::left<<std::setw(7)<<desc._str.atom(i).get_symbol()<<std::setw(4)<<i+1<<std::setw(15)<<std::right<<desc._Qm[i]<<std::setw(15)<<std::right<<desc._Qp[i]<<std::setw(15)<<std::right<<desc._Q0[i]<<std::endl;
        flux<<std::endl;
    }
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    flux<<std::left<<std::setw(7)<<"Symbol"<<std::setw(4)<<"k"<<std::setw(15)<<std::right<<"f-"<<std::setw(15)<< std::right<< "f+"<<std::setw(15)<<std::right<<"f0"<<std::setw(15)<<std::right<<"Delta f"<<std::endl;
    for(int i=0; i<desc._str.getNumberOfAtoms(); i++)
        flux<<std::left<<std::setw(7)<<desc._str.atom(i).get_symbol()<<std::setw(4)<<i+1<<std::setw(15)<<std::right<<desc._fkm[i]<<std::setw(15)<<std::right<<desc._fkp[i]<<std::setw(15)<<std::right<<desc._fk0[i]<<std::setw(15)<<std::right<<desc._deltafk[i]<<std::endl;
    flux<<std::endl;    
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    flux<<std::left<<std::setw(15)<<std::right<<"w-"<<std::setw(15)<<std::right<<"w+"<<std::setw(15)<<std::right<<"s-"<<std::setw(15)<<std::right<<"s+"<<std::setw(15)<<std::right<<"s-/s+"<<std::setw(15)<<std::right<<"hardness-"<<std::setw(15)<<std::right<<"hardness+"<<std::setw(15)<<std::right<<"hardness"<<std::endl;
    for(int i=0; i<desc._str.getNumberOfAtoms(); i++)
        flux<<std::setw(15)<<std::right<<desc._wkm[i]*HeV<<std::setw(15)<<std::right<<desc._wkp[i]*HeV<<std::setw(15)<<std::right<<desc._Skm[i]/HeV<<std::setw(15)<<std::right<<desc._Skp[i]/HeV<<std::setw(15)<<std::right<<desc._Skfrac[i]<<std::setw(15)<<std::right<<desc._hardnesskm[i]*HeV<<std::setw(15)<<std::right<<desc._hardnesskp[i]*HeV<<std::setw(15)<<std::right<<desc._hardnessk[i]*HeV<<std::endl;
    flux<<std::endl;
    
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    flux<<std::left<<std::setw(10)<<"mu+ "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._mup*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"mu- "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._mum*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"mu "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._mu*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"Xi "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._xi*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"hardness "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._hardness*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"w "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._w*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"S "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._S/HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"Qmax "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._Qmax<<std::endl;
    flux<<std::left<<std::setw(10)<<"DEmin "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._DEmin*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"w+ "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._wp*HeV<<std::endl;
    flux<<std::left<<std::setw(10)<<"w- "<<std::setw(2)<<"="<<std::setw(16)<<std::right<<desc._wm*HeV<<std::endl;
    flux<<std::endl;
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    flux<<"Energies (hardness, mu, w, Xi, DEmin, wk-, wk+, hardnessk-, hardnessk+, hardnessk) are given in eV"<<std::endl;
    flux<<"Softnesses (S, sk-, sk+) are given in eV^-1"<<std::endl;
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    flux<<std::left<<std::setw(12)<<"mu-"<<"= -I"<<std::endl;
    flux<<std::left<<std::setw(12)<<"mu+"<<"= A"<<std::endl;
    flux<<std::left<<std::setw(12)<<"mu"<<"= Chemical potential = (mu+ + mu-)/2"<<std::endl;
    flux<<std::left<<std::setw(12)<<"hardness"<<"= Chemical hardness = (mu+  -  mu-)"<<std::endl;
    flux<<std::left<<std::setw(12)<<"Xi"<<"= Electronegativity = -mu"<<std::endl;
    flux<<std::left<<std::setw(12)<<"w"<<"= Electrophilicity index = mu^2/(2 hardness)"<<std::endl;
    flux<<std::left<<std::setw(12)<<"w-"<<"= propensity to donate electron = mu-^2/(2 hardness)"<<std::endl;
    flux<<std::left<<std::setw(12)<<"w+"<<"= propensity to accept electron = mu+^2/(2 hardness)"<<std::endl; 
    flux<<std::left<<std::setw(12)<<"S"<<"= Global softness = 1/hardness"<<std::endl;
    flux<<std::left<<std::setw(12)<<"Qmax"<<"= Maximal electronic charge accepted by an electrophile = -mu/hardness"<<std::endl;
    flux<<std::left<<std::setw(12)<<"DEmin"<<"= Energy decrease if the electrophile take Qmax = -mu^2/(2 hardness)"<<std::endl; 
    flux<<std::left<<std::setw(12)<<"fk-"<<"= Local Fukui electrophilic attack"<<std::endl;
    flux<<std::left<<std::setw(12)<<"fk+"<<"= Local Fukui nucleophilic attack"<<std::endl;
    flux<<std::left<<std::setw(12)<<"sk-"<<"= Local softness electrophilic attack = S fk-"<<std::endl;
    flux<<std::left<<std::setw(12)<<"sk+"<<"= Local softness nucleophilic attack = S fk+"<<std::endl;
    flux<<std::left<<std::setw(12)<<"wk-"<<"= Local philicity index of electrophilic attack = w fk-"<<std::endl;
    flux<<std::left<<std::setw(12)<<"wk+"<<"= Local philicity index of nucleophilic attack = w fk+"<<std::endl;
    flux<<std::left<<std::setw(12)<<"hardnessk-"<<"= Local hardness = mu+ fk+ - mu- fk- - (mu+- mu-)*(fk+-fk-)"<<std::endl;
    flux<<std::left<<std::setw(12)<<"hardnessk+"<<"= Local hardness = mu+ fk+ - mu- fk- + (mu+- mu-)*(fk+-fk-)"<<std::endl;
    flux<<std::left<<std::setw(12)<<"hardnessk"<<"= Local hardness = mu+ fk+ - mu- fk-"<<std::endl;
    flux<<std::left<<std::setw(12)<<"Deltafk"<<"= Dual descripor = (fk+ - fk-) : "<<std::endl;
    flux<<std::left<<std::setw(9)<<" "<<">0 => site favored for a nucleophilic attack"<<std::endl;
    flux<<std::left<<std::setw(9)<<" "<<"<0 => site favored for an electrophilic attack"<<std::endl;
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    flux<<std::left<<std::setw(12)<<"References:"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"- Revisiting the definition of local hardness and hardness kernel"<<std::endl; 
    flux<<std::left<<std::setw(12)<<" "<<"C. A. Polanco-Ramrez et al"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Phys. Chem. Chem. Phys., 2017, 19, 12355-12364"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"DOI: 10.1039/c7cp00691h"<<std::endl;
    flux<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"- Applications of the Conceptual Density Functional Theory"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Indices to Organic Chemistry Reactivity"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Luis R. Domingo, Mar Ríos-Gutiérrez and Patricia Pérez"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Molecules 2016, 21, 748; doi:10.3390/molecules21060748"<<std::endl;
    flux<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"- Electrodonating and Electroaccepting Powers"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"José L. Gazquez, André Cedillo, and Alberto Vela"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"J. Phys. Chem. A 2007, 111, 1966-1970, DOI: 10.1021/jp065459f"<<std::endl;
    flux<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"- Introducing “UCA-FUKUI” software: reactivity-index calculations"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Jesús Sánchez-Márquez et al."<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"J Mol Model (2014) 20:2492, DOI 10.1007/s00894-014-2492-1"<<std::endl;
    flux<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"- Dual descriptor and molecular electrostatic potential:"<<std::endl; 
    flux<<std::left<<std::setw(12)<<" "<<"complementary tools for the study of the coordination"<<std::endl; 
    flux<<std::left<<std::setw(12)<<" "<<"chemistry of ambiphilic ligands"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"F.  Guégan et al."<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Phys.Chem.Chem.Phys., 2014, 16 , 15558-15569,"<<std::endl; 
    flux<<std::left<<std::setw(12)<<" "<<"DOI: 10.1039/c4cp01613k"<<std::endl;
    flux<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"- New Dual Descriptor for Chemical Reactivity"<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"Ch. Morell et al."<<std::endl;
    flux<<std::left<<std::setw(12)<<" "<<"J. Phys. Chem. A 2005, 109, 205-212, DOI: 10.1021/jp046577a"<<std::endl;
    flux<<"------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    return flux;
}

void Descriptors::compute_All_From_Cube(std::ifstream &file1, std::ifstream &file2, std::ifstream &file3, std::vector<double> E, PartitionMethod partitionMethod)
{
    _okCharge=true;
    std::vector<double> Q1 = compute_Charges_From_File(file1, partitionMethod);
    std::vector<double> Q2 = compute_Charges_From_File(file2, partitionMethod);
    std::vector<double> Q3 = compute_Charges_From_File(file3, partitionMethod);
    double I=0;
    double A=0;    
    sortCharges(Q1, Q2, Q3, E, I, A);
    compute_All_From_Charge(I,A);
}

void Descriptors::sortCharges(const std::vector<double>& Q1, const std::vector<double>& Q2, const std::vector<double>& Q3, std::vector<double>& E, double I, double A)
{
    std::vector<std::vector<double>> Q(3);
    Q[0] = Q1;
    Q[1] = Q2;
    Q[2] = Q3;

    std::vector<double> S(3, 0.0);
    for (size_t i = 0; i < Q1.size(); ++i)
    {
        for(size_t c = 0; c < 3; ++c)
        {
            S[c] += Q[c][i];
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        int k = i;
        
        for (int j = i + 1; j < 3; ++j)
        {
            if(S[j] < S[k])
            {
                k=j;
            }
        }

        if (k != i)
        {
            double s=S[k];
            S[k] = S[i];
            S[i] = s;
            double e=E[k];
            E[k] = E[i];
            E[i] = e;
            std::vector<double> q=Q[k];
            Q[k] = Q[i];
            Q[i] = q;
        }
    }

    _Qp=Q[2];
    _Qm=Q[0];
    _Q0=Q[1];
    I=E[2]-E[1];
    A=E[1]-E[0];
}

void Descriptors::compute_All_From_Charges(const Structure& structure, std::vector<double> Q0, std::vector<double> Qm, std::vector<double> Qp, std::vector<double> E)
{
    _okCharge = true;
    _str = structure;

    reset();

    double I=0;
    double A=0;

    sortCharges(Q0, Qm, Qp, E, I, A);
    compute_All_From_Charge(I, A);
}

