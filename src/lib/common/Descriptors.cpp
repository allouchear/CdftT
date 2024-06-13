#include <iomanip>
#include <common/Descriptors.h>

using namespace std;

Descriptors::Descriptors()
{
	_Qp=vector<double> ();
	_Q0=vector<double> ();
	_Qm=vector<double> ();
	_mu=0.0;
	_mup=0.0;
	_mum=0.0;
	_fk0=vector<double> ();
	_fkm=vector<double> ();
	_fkp=vector<double> ();
}

Descriptors::Descriptors(WFX& wfx, const PeriodicTable& Table)
{
	_str=Structure(wfx, Table);
	_Deltafk.resize(_str.number_of_atoms());
	_wkm.resize(_str.number_of_atoms());
	_wkp.resize(_str.number_of_atoms());
	_Skm.resize(_str.number_of_atoms());
	_Skp.resize(_str.number_of_atoms());
	_Skfrac.resize(_str.number_of_atoms());
	_hardnessk.resize(_str.number_of_atoms());
	_hardnesskm.resize(_str.number_of_atoms());
	_hardnesskp.resize(_str.number_of_atoms());

	_fk0.resize(_str.number_of_atoms());
}

void Descriptors::set_all_mu(const double& I,const double& A)
{
	_mum = -I;
	_mup = -A;
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
	_DEmin = -(_mu*_mu)/(2*_hardness);

	for(int i=0;i<_str.number_of_atoms();i++)
	{
		_Deltafk[i] = _fkm[i]-_fkp[i];
		_fk0[i] = 0.5 * (_fkp[i]+_fkm[i]);
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

void Descriptors::compute_fk()
{
	_fk0.resize(_str.number_of_atoms());
	_fkp.resize(_str.number_of_atoms());
	_fkm.resize(_str.number_of_atoms());
	for(int i=0; i<_str.number_of_atoms();i++)
	{
		_fk0[i] = 0.5*(_Qm[i]-_Qp[i]);
		_fkp[i] = _Qm[i]-_Q0[i];
		_fkm[i] = _Q0[i]-_Qp[i];
	}
}

void Descriptors::set_mu_fk_data(vector<vector<double>> f, double eH, double eL)
{
	_mum=eH;
	_mup=eL;
	_mu = 0.5*(_mup+_mum);
	_fkm=f[0];
	_fkp=f[1];
}

void Descriptors::set_mu_fk_data(vector<vector<double>> data)
{
	_fk0.resize(_str.number_of_atoms());
	_fkp.resize(_str.number_of_atoms());
	_fkm.resize(_str.number_of_atoms());

	_mum= -(data[2][0]-data[0][0]);
	_mup= -(data[0][0]-data[1][0]);
	_mu = 0.5*(_mup+_mum);

	for(int i=0; i<_str.number_of_atoms(); i++)
	{
		_fkm[i]=data[0][i+1]-data[1][i+1];
		_fkp[i]=data[2][i+1]-data[0][i+1];
		_fk0[i]=0.5*(data[2][i+1]-data[1][i+1]);
	}

}

/********************************************************************************************/

ostream& operator<<(ostream& flux, const Descriptors& desc)
{
	double HeV=27.21138469;
	
	flux<<scientific;
	flux<<setprecision(6);
	flux<<setw(15);
	flux<<left<<setw(7)<<"Symbol"<<setw(4)<<"k"<<setw(15)<<"f-"<<setw(15)<<"f+"<<setw(15)<<"f0"<<setw(15)<<"Delta f"
	<<setw(15)<<"w-"<<setw(15)<<"w+"<<setw(15)<<"s-"<<setw(15)<<"s+"<<setw(15)<<"s-/s+"<<setw(15)
	<<"hardness-"<<setw(15)<<"hardness+"<<setw(15)<<"hardness"<<endl;
	for(int i=0; i<desc._str.number_of_atoms(); i++)
		flux<<left<<setw(7)<<desc._str.atoms(i).symbol()<<setw(4)<<i+1<<setw(15)<<desc._fkm[i]<<setw(15)<<desc._fkp[i]<<
		setw(15)<<desc._fk0[i]<<setw(15)<<desc._Deltafk[i]<<setw(15)<<desc._wkm[i]*HeV<<setw(15)
		<<desc._wkp[i]*HeV<<setw(15)<<desc._Skm[i]/HeV<<setw(15)<<desc._Skp[i]/HeV<<setw(15)<<desc._Skfrac[i]
		<<setw(15)<<desc._hardnesskm[i]*HeV<<setw(15)<<desc._hardnesskp[i]*HeV<<setw(15)<<desc._hardnessk[i]*HeV<<endl;
	
	flux<<endl;
	
	flux<<left<<setw(10)<<"mu+ "<<setw(2)<<"="<<setw(10)<<desc._mup*HeV<<endl;
	flux<<left<<setw(10)<<"mu- "<<setw(2)<<"="<<setw(10)<<desc._mum*HeV<<endl;
	flux<<left<<setw(10)<<"mu "<<setw(2)<<"="<<setw(10)<<desc._mu*HeV<<endl;
	flux<<left<<setw(10)<<"Xi "<<setw(2)<<"="<<setw(10)<<desc._xi*HeV<<endl;
	flux<<left<<setw(10)<<"hardness "<<setw(2)<<"="<<setw(10)<<desc._hardness*HeV<<endl;
	flux<<left<<setw(10)<<"w "<<setw(2)<<"="<<setw(10)<<desc._w*HeV<<endl;
	flux<<left<<setw(10)<<"S "<<setw(2)<<"="<<setw(10)<<desc._S/HeV<<endl;
	flux<<left<<setw(10)<<"Qmax "<<setw(2)<<"="<<setw(10)<<desc._Qmax<<endl;
	flux<<left<<setw(10)<<"DEmin "<<setw(2)<<"="<<setw(10)<<desc._DEmin*HeV<<endl;
	flux<<left<<setw(10)<<"w+ "<<setw(2)<<"="<<setw(10)<<desc._wp*HeV<<endl;
	flux<<left<<setw(10)<<"w- "<<setw(2)<<"="<<setw(10)<<desc._wm*HeV<<endl;
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

