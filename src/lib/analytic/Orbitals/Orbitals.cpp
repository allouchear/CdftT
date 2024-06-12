#include<iostream>
#include<iomanip>
#include<analytic/Orbitals/Orbitals.h>

using namespace std;

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
	_occupation_number=vector<double> ();
}

Orbitals::Orbitals(WFX& wfx, Binomial& Bin, const PeriodicTable& Table)
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
	_numOrb=vector<int> (2,0);
	_occupation_number=wfx.Molecular_Orbital_Occupation_Numbers();
	_descriptors=Descriptors(wfx, Table);
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

void Orbitals::PrintDescriptors()
{
	HOMO_LUMO();
	get_f();
	_descriptors.set_mu_fk_data(_all_f, eHOMO(), eLUMO());
	_descriptors.compute_all();
	cout<<_descriptors<<endl;
}

void Orbitals::PrintDescriptors(int i, int j)
{
	HOMO_LUMO(i,j);
	get_f();
	_descriptors.set_mu_fk_data(_all_f, eHOMO(), eLUMO());
	_descriptors.compute_all();
	cout<<_descriptors<<endl;
}