#include<iostream>
#include<iomanip>
#include<analytic/Orbitals/Orbitals.h>

using namespace std;

Density::Density()
{
	_gtf = vector<GTF> ();
	_occupation_number = vector<double> ();
	_orbital_coefficients = vector<vector<double>> ();
}

Density::Density(WFX& wfx, Binomial& Bin)
{
	int i,j;

	_gtf = vector<GTF> (wfx.Number_of_Primitives());
	_orbital_coefficients = vector<vector<double>> (wfx.Number_of_Occupied_Molecular_Orbital(), vector<double> (wfx.Number_of_Primitives()));
	vector<vector<double>> Coord(wfx.Number_of_Nuclei(), vector<double> (0));

	for(i=0; i<wfx.Number_of_Nuclei(); i++)
		for(j=i*3; j<3*(1+i); j++)
			Coord[i].push_back(wfx.Nuclear_Cartesian_Coordinates()[j]);

	for(i=0; i<wfx.Number_of_Primitives(); i++)		
		_gtf[i].push_back(wfx.Primitive_Exponents()[i], 1.0, Coord[wfx.Primitive_Centers()[i]-1], setLxyz(wfx.Primitive_Types()[i]), Bin);

	for(i=0; i<wfx.Number_of_Occupied_Molecular_Orbital(); i++)
		for(j=0; j<wfx.Number_of_Primitives(); j++)
			_orbital_coefficients[i][j]=wfx.Molecular_Orbital_Primitive_Coefficients()[i].Coefficients()[j];

	_occupation_number = wfx.Molecular_Orbital_Occupation_Numbers();
}

Orbitals::Orbitals()
{
	_all_f = vector<vector<double>> ();
	_lcao = vector<LCAO> ();
	_symbol = vector<string> ();
	_orbital_energy = vector<double> ();
	_primitive_centers = vector<int> ();
	_numOrb = vector<int> ();
	_numberOfFunctions=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_number_of_atoms=0;
}

Orbitals::Orbitals(WFX& wfx, Binomial& Bin)
{
	int i,j;

	GTF gtf;
	vector<vector<CGTF>> vcgtf (wfx.Number_of_Occupied_Molecular_Orbital(), vector<CGTF> (wfx.Number_of_Primitives()));
	_lcao = vector<LCAO> (wfx.Number_of_Occupied_Molecular_Orbital());
	vector<vector<double>> Coord(wfx.Number_of_Nuclei(), vector<double> (0));
	
	for(i=0; i<wfx.Number_of_Nuclei(); i++)
		for(j=i*3; j<3*(1+i); j++)
			Coord[i].push_back(wfx.Nuclear_Cartesian_Coordinates()[j]);

	for(i=0; i<wfx.Number_of_Occupied_Molecular_Orbital(); i++)
	{
		gtf=GTF();
		for(j=0; j<wfx.Number_of_Primitives(); j++)
		{			
			gtf.push_back(wfx.Primitive_Exponents()[j], 1.0, Coord[wfx.Primitive_Centers()[j]-1], setLxyz(wfx.Primitive_Types()[j]), Bin);
			vcgtf[i][j].push_back(gtf);
			_lcao[i].push_back(vcgtf[i][j], wfx.Molecular_Orbital_Primitive_Coefficients()[i].Coefficients()[j]);
		}
	}
	_density = Density(wfx, Bin);
	_primitive_centers=wfx.Primitive_Centers();
	_numberOfFunctions=wfx.Number_of_Occupied_Molecular_Orbital();
	_number_of_alpha_electrons=wfx.Number_of_Alpha_Electrons();
	_number_of_beta_electrons=wfx.Number_of_Beta_Electrons();
	_number_of_atoms=wfx.Number_of_Nuclei();
	_orbital_energy=wfx.Molecular_Orbital_Energies();
	_symbol=wfx.Nuclear_Names();
	_numOrb = vector<int> (2,0);
}

void Orbitals::PrintOverlap(int i, int j)
{
	cout<<"OverlapLCAO <"<<i<<"|"<<j<<"> = "<<_lcao[i].overlapLCAO(_lcao[j])<<endl;
}

void Orbitals::HOMO()
{
	if(_number_of_alpha_electrons>=_number_of_beta_electrons)
		_numOrb[0] = _number_of_alpha_electrons-1;
	else
		_numOrb[0] = _number_of_beta_electrons-1;
}

void Orbitals::LUMO()
{
	_numOrb[1] = _numOrb[0]+1;
}

vector<vector<double>> Orbitals::get_S() const
{
	int i,j,n;
	n=_lcao[0].numberOfFunctions();
	vector<vector<double>> S (n, vector<double> (n,0.0));
	vector<CGTF> cgtf =_lcao[0].cgtf();

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,j)
#endif
	for(i=0; i<n; i++)
		for(j=i; j<n; j++)
			S[i][j]=S[j][i]=cgtf[i].overlapCGTF(cgtf[j]);

	return S;
}

vector<double> Orbitals::get_f(int alpha) const
{
	int i;
	vector<vector<double>> S =get_S();
	int nu,xi;
	vector<double> f(_number_of_atoms,0.0);

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,nu,xi)
#endif
	for(i=0; i<_number_of_atoms; i++)
		for(nu=0; nu<_lcao[alpha].numberOfCoefficient(); nu++)
		{
			if(i+1 == _primitive_centers[nu])
				f[i]+=_lcao[alpha].coefficient()[nu]*_lcao[alpha].coefficient()[nu];

			for(xi=0; xi<_lcao[alpha].numberOfCoefficient(); xi++)
				if(xi!=nu && i+1 == _primitive_centers[nu])
					f[i]+=_lcao[alpha].coefficient()[xi]*_lcao[alpha].coefficient()[nu]*S[xi][nu];
		}

	return f;
}

void Orbitals::get_f()
{
	vector<vector<double>> S=get_S();
	int i,nu,xi;
	size_t j;
	vector<double> V(_number_of_atoms,0.0);
	vector<vector<double>> f(_numOrb.size(), V);

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,j,nu,xi)
#endif
	for(i=0; i<_number_of_atoms; i++)
		for(j=0; j<_numOrb.size(); j++)
			for(nu=0; nu<_lcao[_numOrb[j]].numberOfCoefficient(); nu++)
			{
				if(i+1 == _primitive_centers[nu])
					f[j][i]+=_lcao[_numOrb[j]].coefficient()[nu]*_lcao[_numOrb[j]].coefficient()[nu];

				for(xi=0; xi<_lcao[_numOrb[j]].numberOfCoefficient(); xi++)
					if(xi!=nu && i+1 == _primitive_centers[nu])
						f[j][i]+=_lcao[_numOrb[j]].coefficient()[xi]*_lcao[_numOrb[j]].coefficient()[nu]*S[xi][nu];
			}

	_all_f = f;
}

vector<double> Orbitals::fk0() const
{
	vector<double> f0 (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		f0[i]=(_all_f[0][i]+_all_f[1][i])/2;

	return f0;
}

vector<double> Orbitals::Deltafk() const
{
	vector<double> df (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		df[i]=_all_f[1][i]-_all_f[0][i];

	return df;
}

vector<double> Orbitals::skminus() const
{
	vector<double> skm (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		skm[i]=softness()*_all_f[0][i];

	return skm;
}

vector<double> Orbitals::skplus() const
{
	vector<double> skp (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		skp[i]=softness()*_all_f[1][i];

	return skp;
}

vector<double> Orbitals::skfrac() const
{
	vector<double> skf (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		skf[i]=skminus()[i]/skplus()[i];

	return skf;
}

vector<double> Orbitals::wkminus() const
{
	vector<double> wkm (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		wkm[i]=w()*_all_f[0][i];

	return wkm;
}

vector<double> Orbitals::wkplus() const
{
	vector<double> wkp (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		wkp[i]=w()*_all_f[1][i];

	return wkp;
}

vector<double> Orbitals::hardnesskminus() const
{
	vector<double> hkm (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		hkm[i]=muplus()*_all_f[1][i]-muminus()*_all_f[0][i]-(muplus()-muminus())*(_all_f[1][i]-_all_f[0][i]);

	return hkm;
}

vector<double> Orbitals::hardnesskplus() const
{
	vector<double> hkp (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		hkp[i]=muplus()*_all_f[1][i]-muminus()*_all_f[0][i]+(muplus()-muminus())*(_all_f[1][i]-_all_f[0][i]);

	return hkp;
}

vector<double> Orbitals::hardnessk() const
{
	vector<double> hk (_number_of_atoms);

	for(int i=0; i<_number_of_atoms; i++)
		hk[i]=muplus()*_all_f[1][i]-muminus()*_all_f[0][i];

	return hk;
}

void Orbitals::HOMO_LUMO()
{
	HOMO();
	LUMO();
}

void Orbitals::HOMO_LUMO(int i, int j)
{
	_numOrb[0]=i;
	_numOrb[1]=j;
}

ostream& operator<<(ostream& flux, const Orbitals& orb)
{
	double HeV=27.21138469;

	flux<<std::scientific;
	flux<<std::setprecision(6);
	flux<<std::setw(15);
	flux<<left<<setw(7)<<"Symbol"<<setw(4)<<"k"<<setw(15)<<"f-"<<setw(15)<<"f+"<<setw(15)<<"f0"<<setw(15)<<"Delta f"<<
	setw(15)<<"w-"<<setw(15)<<"w+"<<setw(15)<<"s-"<<setw(15)<<"s+"<<setw(15)<<"s-/s+"<<setw(15)
	<<"hardness-"<<setw(15)<<"hardness+"<<setw(15)<<"hardness"<<endl;
	for(int i=0; i<orb.NumberOfAtoms(); i++)
		flux<<left<<setw(7)<<orb.symbol()[i]<<setw(4)<<i+1<<setw(15)<<orb.fkminus()[i]<<setw(15)<<orb.fkplus()[i]<<
		setw(15)<<orb.fk0()[i]<<setw(15)<<orb.Deltafk()[i]<<setw(15)<<orb.wkminus()[i]*HeV<<setw(15)
		<<orb.wkplus()[i]*HeV<<setw(15)<<orb.skminus()[i]/HeV<<setw(15)<<orb.skplus()[i]/HeV<<setw(15)<<orb.skfrac()[i]
		<<setw(15)<<orb.hardnesskminus()[i]*HeV<<setw(15)<<orb.hardnesskplus()[i]*HeV<<setw(15)<<orb.hardnessk()[i]*HeV<<endl;

	flux<<endl;

	flux<<left<<setw(10)<<"mu+ "<<setw(2)<<"="<<setw(10)<<orb.muplus()*HeV<<endl;
	flux<<left<<setw(10)<<"mu- "<<setw(2)<<"="<<setw(10)<<orb.muminus()*HeV<<endl;
	flux<<left<<setw(10)<<"mu "<<setw(2)<<"="<<setw(10)<<orb.mu()*HeV<<endl;
	flux<<left<<setw(10)<<"Xi "<<setw(2)<<"="<<setw(10)<<orb.electronegativity()*HeV<<endl;
	flux<<left<<setw(10)<<"hardness "<<setw(2)<<"="<<setw(10)<<orb.hardness()*HeV<<endl;
	flux<<left<<setw(10)<<"w "<<setw(2)<<"="<<setw(10)<<orb.w()*HeV<<endl;
	flux<<left<<setw(10)<<"S "<<setw(2)<<"="<<setw(10)<<orb.softness()/HeV<<endl;
	flux<<left<<setw(10)<<"Qmax "<<setw(2)<<"="<<setw(10)<<orb.Qmax()<<endl;
	flux<<left<<setw(10)<<"DEmin "<<setw(2)<<"="<<setw(10)<<orb.DEmin()*HeV<<endl;
	flux<<left<<setw(10)<<"w+ "<<setw(2)<<"="<<setw(10)<<orb.wplus()*HeV<<endl;
	flux<<left<<setw(10)<<"w- "<<setw(2)<<"="<<setw(10)<<orb.wminus()*HeV<<endl;
	flux<<endl;
	flux<<"------------------------------------------------------------------------------------------------------------------------------"<<endl;
 	flux<<"Energies (hardness, mu, w, Xi, DEmin, wk-, wk+, hardnessk-, hardnessk+, hardnessk) are given in eV"<<endl;
 	flux<<"Softnesses (S, sk-, sk+) are given in eV^-1"<<endl;
	flux<<"------------------------------------------------------------------------------------------------------------------------------"<<endl;
   	flux<<left<<setw(12)<<"mu-"<<"= eHOMO"<<endl;
   	flux<<left<<setw(12)<<"mu+"<<"= eLUMO"<<endl;
   	flux<<left<<setw(12)<<"mu"<<"= Chemical potential = (mu+ + mu-)/2"<<endl;
   	flux<<left<<setw(12)<<"hardness"<<"= Chemical hardness = (mu+  -  mu-)"<<endl;
   	flux<<left<<setw(12)<<"Xi"<<"= Electronegativity = -mu"<<endl;
   	flux<<left<<setw(12)<<"w"<<"= Electrophilicity index = mu^2/(2 hardness)"<<endl;
  	flux<<left<<setw(12)<<"w-"<<"= propensity to donate electron = mu-^2/(2 hardness)"<<endl;
   	flux<<left<<setw(12)<<"w+"<<"= propensity to accept electron = mu+^2/(2 hardness)"<<endl; 
   	flux<<left<<setw(12)<<"S"<<"= Global softness = 1/hardness"<<endl;
   	flux<<left<<setw(12)<<"Qmax"<<"= Maximal electronic charge accepted by an electrophile = -mu/hardness"<<endl;
   	flux<<left<<setw(12)<<"DEmin"<<"= Energy decrease if the electrophile take Qmax = -mu^2/(2 hardness)"<<endl; 
   	flux<<left<<setw(12)<<"fk-"<<"= Local Fukui electrophilic attack"<<endl;
   	flux<<left<<setw(12)<<"fk+"<<"= Local Fukui nucleophilic attack"<<endl;
   	flux<<left<<setw(12)<<"sk-"<<"= Local softness electrophilic attack = S fk-"<<endl;
   	flux<<left<<setw(12)<<"sk+"<<"= Local softness nucleophilic attack = S fk+"<<endl;
   	flux<<left<<setw(12)<<"wk-"<<"= Local philicity index of electrophilic attack = w fk-"<<endl;
   	flux<<left<<setw(12)<<"wk+"<<"= Local philicity index of nucleophilic attack = w fk+"<<endl;
   	flux<<left<<setw(12)<<"hardnessk-"<<"= Local hardness = mu+ fk+ - mu- fk- - (mu+- mu-)*(fk+-fk-)"<<endl;
   	flux<<left<<setw(12)<<"hardnessk+"<<"= Local hardness = mu+ fk+ - mu- fk- + (mu+- mu-)*(fk+-fk-)"<<endl;
   	flux<<left<<setw(12)<<"hardnessk"<<"= Local hardness = mu+ fk+ - mu- fk-"<<endl;
   	flux<<left<<setw(12)<<"Deltafk"<<"= Dual descripor = (fk+ - fk-) : "<<endl;
    flux<<left<<setw(9)<<" "<<">0 => site favored for a nucleophilic attack"<<endl;
    flux<<left<<setw(9)<<" "<<"<0 => site favored for an electrophilic attack"<<endl;
	flux<<"------------------------------------------------------------------------------------------------------------------------------"<<endl;
 	flux<<left<<setw(12)<<"References:"<<endl;
  	flux<<left<<setw(12)<<" "<<"- Revisiting the definition of local hardness and hardness kernel"<<endl; 
    flux<<left<<setw(12)<<" "<<"C. A. Polanco-Ramrez et al"<<endl;
    flux<<left<<setw(12)<<" "<<"Phys. Chem. Chem. Phys., 2017, 19, 12355-12364"<<endl;
    flux<<left<<setw(12)<<" "<<"DOI: 10.1039/c7cp00691h"<<endl;
    flux<<endl;
  	flux<<left<<setw(12)<<" "<<"- Applications of the Conceptual Density Functional Theory"<<endl;
    flux<<left<<setw(12)<<" "<<"Indices to Organic Chemistry Reactivity"<<endl;
    flux<<left<<setw(12)<<" "<<"Luis R. Domingo, Mar Ríos-Gutiérrez and Patricia Pérez"<<endl;
    flux<<left<<setw(12)<<" "<<"Molecules 2016, 21, 748; doi:10.3390/molecules21060748"<<endl;
    flux<<endl;
  	flux<<left<<setw(12)<<" "<<"- Electrodonating and Electroaccepting Powers"<<endl;
    flux<<left<<setw(12)<<" "<<"José L. Gazquez, André Cedillo, and Alberto Vela"<<endl;
    flux<<left<<setw(12)<<" "<<"J. Phys. Chem. A 2007, 111, 1966-1970, DOI: 10.1021/jp065459f"<<endl;
    flux<<endl;
  	flux<<left<<setw(12)<<" "<<"- Introducing “UCA-FUKUI” software: reactivity-index calculations"<<endl;
    flux<<left<<setw(12)<<" "<<"Jesús Sánchez-Márquez et al."<<endl;
    flux<<left<<setw(12)<<" "<<"J Mol Model (2014) 20:2492, DOI 10.1007/s00894-014-2492-1"<<endl;
    flux<<endl;
  	flux<<left<<setw(12)<<" "<<"- Dual descriptor and molecular electrostatic potential:"<<endl; 
    flux<<left<<setw(12)<<" "<<"complementary tools for the study of the coordination"<<endl; 
    flux<<left<<setw(12)<<" "<<"chemistry of ambiphilic ligands"<<endl;
    flux<<left<<setw(12)<<" "<<"F.  Guégan et al."<<endl;
    flux<<left<<setw(12)<<" "<<"Phys.Chem.Chem.Phys., 2014, 16 , 15558-15569,"<<endl; 
    flux<<left<<setw(12)<<" "<<"DOI: 10.1039/c4cp01613k"<<endl;
    flux<<endl;
  	flux<<left<<setw(12)<<" "<<"- New Dual Descriptor for Chemical Reactivity"<<endl;
    flux<<left<<setw(12)<<" "<<"Ch. Morell et al."<<endl;
    flux<<left<<setw(12)<<" "<<"J. Phys. Chem. A 2005, 109, 205-212, DOI: 10.1021/jp046577a"<<endl;
	flux<<"------------------------------------------------------------------------------------------------------------------------------"<<endl;

	return flux;
}