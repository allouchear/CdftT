#include<iostream>
#include<iomanip>
#include<analytic/Orbitals.h>

using namespace std;

Orbitals::Orbitals()
{
	_all_f = vector<vector<double>> ();
	_lcao = vector<LCAO> ();
	_symbol = vector<string> ();
	_orbital_energy = vector<double> ();
	_primitive_centers = vector<int> ();
	_numberOfFunctions=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_number_of_atoms=0;
	_HeV=27.21138469;
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

	_primitive_centers=wfx.Primitive_Centers();
	_numberOfFunctions=wfx.Number_of_Occupied_Molecular_Orbital();
	_number_of_alpha_electrons=wfx.Number_of_Alpha_Electrons();
	_number_of_beta_electrons=wfx.Number_of_Beta_Electrons();
	_number_of_atoms=wfx.Number_of_Nuclei();
	_orbital_energy=wfx.Molecular_Orbital_Energies();
	_symbol=wfx.Nuclear_Names();
	vector<int> NumOrb (2);
	NumOrb[0]=HOMO();
	NumOrb[1]=LUMO();
	get_f(NumOrb);
}

void Orbitals::PrintOverlap(int i, int j)
{
	cout<<"OverlapLCAO <"<<i<<"|"<<j<<"> = "<<_lcao[i].overlapLCAO(_lcao[j])<<endl;
}

int Orbitals::HOMO() const
{
	if(_number_of_alpha_electrons>=_number_of_beta_electrons)
		return _number_of_alpha_electrons-1;
	else
		return _number_of_beta_electrons-1;
}

int Orbitals::LUMO() const
{
	return HOMO()+1;
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

void Orbitals::get_f(vector<int> NumOrb)
{
	vector<vector<double>> S=get_S();
	int i,nu,xi;
	size_t j;
	vector<double> V(_number_of_atoms,0.0);
	vector<vector<double>> f(NumOrb.size(), V);

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,j,nu,xi)
#endif
	for(i=0; i<_number_of_atoms; i++)
		for(j=0; j<NumOrb.size(); j++)
			for(nu=0; nu<_lcao[NumOrb[j]].numberOfCoefficient(); nu++)
			{
				if(i+1 == _primitive_centers[nu])
					f[j][i]+=_lcao[NumOrb[j]].coefficient()[nu]*_lcao[NumOrb[j]].coefficient()[nu];

				for(xi=0; xi<_lcao[NumOrb[j]].numberOfCoefficient(); xi++)
					if(xi!=nu && i+1 == _primitive_centers[nu])
						f[j][i]+=_lcao[NumOrb[j]].coefficient()[xi]*_lcao[NumOrb[j]].coefficient()[nu]*S[xi][nu];
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

	return flux;
}