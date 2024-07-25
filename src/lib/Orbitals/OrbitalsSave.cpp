#include<iostream>
#include<iomanip>
#include<algorithm>
#include <Orbitals/Orbitals.h>
#include <Utils/LM.h>

using namespace std;


void Orbitals::Save(string& tag)
{
	if(tag.find(".wfx")!=string::npos)
		Save_wfx(tag);
	else if(tag.find(".molden")!=string::npos)
		Save_molden(tag);
	else if(tag.find(".gab")!=string::npos)
		Save_gab(tag);
	else
	{
		cout<<"Format not recognize, please choose a valide format"<<endl;
		exit(1);
	}
}

void Orbitals::Save_wfx(string& tag)
{
	int n, nMO;
	ofstream s;
	s.open(tag);

	s<<"<Title>"<<endl;
	s<<" Input file generated by cdftt..."<<endl;
	s<<"</Title>"<<endl;

	s<<"<Keywords>"<<endl;
	s<<" GTO"<<endl;
	s<<"</Keywords>"<<endl;

	s<<"<Number of Nuclei>"<<endl;
	s<<" "<<_number_of_atoms<<endl;
	s<<"</Number of Nuclei>"<<endl;

	s<<"<Number of Occupied Molecular Orbitals>"<<endl;
	if(!_alpha_and_beta)
		nMO=_numberOfMo*2;
	else
		nMO=_numberOfMo;
	s<<" "<<nMO<<endl;
	s<<"</Number of Occupied Molecular Orbitals>"<<endl;

	s<<"<Number of Perturbations>"<<endl;
	s<<" 0"<<endl;
	s<<"</Number of Perturbations>"<<endl;

	s<<"<Net Charge>"<<endl;
	s<<" 0"<<endl;
	s<<"</Net Charge>"<<endl;

	s<<"<Number of Electrons>"<<endl;
	s<<" "<<_number_of_alpha_electrons+_number_of_beta_electrons<<endl;
	s<<"</Number of Electrons>"<<endl;

	s<<"<Number of Alpha Electrons>"<<endl;
	s<<" "<<_number_of_alpha_electrons<<endl;
	s<<"</Number of Alpha Electrons>"<<endl;

	s<<"<Number of Beta Electrons>"<<endl;
	s<<" "<<_number_of_beta_electrons<<endl;
	s<<"</Number of Beta Electrons>"<<endl;

	s<<"<Nuclear Names>"<<endl;
	for(int i=0; i<_number_of_atoms; i++)
		s<<" "<<_symbol[i]<<i+1<<endl;
	s<<"</Nuclear Names>"<<endl;

	s<<"<Atomic Numbers>"<<endl;
	for(int i=0; i<_number_of_atoms; i++)
		s<<" "<<_atomic_numbers[i]<<endl;
	s<<"</Atomic Numbers>"<<endl;	

	s<<"<Nuclear Charges>"<<endl;
	for(int i=0; i<_number_of_atoms; i++)
		s<<right<<scientific<<setprecision(12)<<" "<<setw(20)<<double(_atomic_numbers[i])<<endl;
	s<<"</Nuclear Charges>"<<endl;

	s<<"<Nuclear Cartesian Coordinates>"<<endl;
	for(int i=0; i<_number_of_atoms; i++)
	{
		for(int j=i*3; j<(i+1)*3; j++)
			s<<right<<scientific<<setprecision(12)<<" "<<setw(20)<<_coordinates[j];
		s<<endl;
	}
	s<<"</Nuclear Cartesian Coordinates>"<<endl;

	s<<"<Number of Primitives>"<<endl;
	s<<" "<<_number_of_gtf<<endl;
	s<<"</Number of Primitives>"<<endl;

	s<<"<Primitive Centers>"<<endl;

	if(int(_vcgtf_non_normalise.size())==_number_of_gtf)
		for(int i=0; i<_number_of_gtf; i++)
		{
			s<<"\t"<<_primitive_centers[i];
			if((i+1)%5==0)
				s<<endl;
		}
	else
	{
		vector<int> pc;
		for(size_t i=0; i<_vcgtf.size(); i++)
			for(int j=0; j<_vcgtf[i].numberOfFunctions(); j++)
				pc.push_back(_vcgtf[i].NumCenter());

		for(int i=0; i<_number_of_gtf; i++)
		{
			s<<"\t"<<pc[i];
			if((i+1)%5==0)
				s<<endl;
		}
	}
	s<<endl;
	s<<"</Primitive Centers>"<<endl;

	n=0;
	s<<"<Primitive Types>"<<endl;
	for(size_t i=0; i<_vcgtf.size(); i++)
		for(int j=0; j<_vcgtf[i].numberOfFunctions(); j++)
		{
			s<<"\t"<<getwfxType(_vcgtf[i].gtf()[j].l());
			if((n+1)%5==0)
				s<<endl;
			n++;
		}
	s<<endl;
	s<<"</Primitive Types>"<<endl;

	n=0;
	s<<"<Primitive Exponents>"<<endl;
	for(size_t i=0; i<_vcgtf.size(); i++)
		for(int j=0; j<_vcgtf[i].numberOfFunctions(); j++)
		{
			s<<right<<scientific<<setprecision(12)<<"\t"<<setw(20)<<_vcgtf[i].gtf()[j].exposant();
			if((n+1)%5==0)
				s<<endl;
			n++;
		}
	s<<endl;
	s<<"</Primitive Exponents>"<<endl;

	int m;

	if(_alpha_and_beta)
		m=1;
	else
		m=2;

	s<<"<Molecular Orbital Occupation Numbers>"<<endl;
	for(int i=0; i<m; i++)
		for(int j=0; j<_numberOfMo; j++)
			s<<" "<<_occupation_number[i][j]<<endl;
	s<<"</Molecular Orbital Occupation Numbers>"<<endl;

	s<<"<Molecular Orbital Energies>"<<endl;
	for(int i=0; i<m; i++)
		for(size_t j=0; j<_orbital_energy[i].size(); j++)
			s<<right<<scientific<<setprecision(12)<<" "<<setw(20)<<_orbital_energy[i][j]<<endl;
	s<<"</Molecular Orbital Energies>"<<endl;

	s<<"<Molecular Orbital Spin Types>"<<endl;
	if(m==1)
		for(int i=0; i<_numberOfMo; i++)
			s<<" "<<"Alpha and Beta"<<endl;
	else
	{
		for(size_t i=0; i<_orbital_energy[0].size(); i++)
			s<<" "<<"Alpha"<<endl;
		for(size_t i=0; i<_orbital_energy[1].size(); i++)
			s<<" "<<"Beta"<<endl;
	}
	s<<"</Molecular Orbital Spin Types>"<<endl;

	s<<"<Molecular Orbital Primitive Coefficients>"<<endl;
	int a=0,b=0;
	for(int i=0; i<m; i++)
		for(int j=0; j<_numberOfMo; j++)
		{
			a++;
			s<<"<MO Number>"<<endl;
			s<<" "<<a<<endl;
			s<<"</MO Number>"<<endl;
		
			for(size_t k=0; k<_vcgtf.size(); k++)
				for(int l=0; l<_vcgtf[k].numberOfFunctions(); l++)
				{
					s<<right<<scientific<<setprecision(12)<<" "<<setw(20)<<_vcgtf[k].coefficients()[l]*_vcgtf[k].gtf()[l].coefficient()*_coefficients[i][j][k];

					if((b+1)%4==0)
						s<<endl;
					b++;
				}

			if(b%4!=0)
				s<<endl;
		}
	s<<"</Molecular Orbital Primitive Coefficients>"<<endl;

	s<<"<Energy = T + Vne + Vee + Vnn>"<<endl;
	s<<" "<<_energy<<endl;
	s<<"</Energy = T + Vne + Vee + Vnn>"<<endl;

	s<<"<Virial Ratio (-V/T)>"<<endl;
	s<<" "<<scientific<<setprecision(12)<<2<<endl;
	s<<"</Virial Ratio (-V/T)>"<<endl;

	s.close();
}
															// Moldengab faire attention au format sphe/cart !!!!
void Orbitals::Save_molden(string& tag)	
{
	if(int(_vcgtf_non_normalise.size())==_number_of_gtf)
	{
		cout<<"This option is nnot implemente."<<endl;
		return;
	}

	ofstream s;
	s.open(tag);

	s<<"[Molden Format]"<<endl;

	s<<"[Atoms] AU"<<endl;
	for(int i=0; i<_number_of_atoms; i++)
	{
		s<<_symbol[i]<<" "<<i+1<<" "<<_atomic_numbers[i];
		for(int j=i*3; j<(i+1)*3; j++)
			s<<std::fixed<<setprecision(6)<<right<<"  "<<setw(10)<<_coordinates[j];
		s<<endl;
	}

	s<<endl;
	
	bool d=false, f=false, g=false;

	for(size_t i=0; i<_vcgtf_non_normalise.size(); i++)
	{
		if(_vcgtf_non_normalise[i].Lformat()=="Sphe")
		{
			if((_vcgtf_non_normalise[i].Ltype()=="D" || _vcgtf_non_normalise[i].Ltype()=="d") && !d)
			{
				s<<"[5D] ";
				d=true;
			}
			else if((_vcgtf_non_normalise[i].Ltype()=="F" || _vcgtf_non_normalise[i].Ltype()=="f") &&!f)
			{
				s<<"[7F] ";
				f=true;
			}
			else if((_vcgtf_non_normalise[i].Ltype()=="G" || _vcgtf_non_normalise[i].Ltype()=="d") && !g)
			{
				s<<"[9G] ";
				g=true;
			}
		}
	}

	s<<endl<<"[GTO]"<<endl;

	int lt,q,m,k=0;
	double save_alpha;

	for(int i=0; i<_number_of_atoms; i++)
	{
		s<<"\t"<<_vcgtf_non_normalise[k].NumCenter()<<" "<<0<<endl;

		do{
			if(k==0)	//First shell
			{
				s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
				s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

				for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
					s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

				k++;
			}

			else if(k+1<int(_vcgtf_non_normalise.size()))	//Other shell
			{
				if(_vcgtf_non_normalise[k].Ltype()==_vcgtf_non_normalise[k+1].Ltype() && _vcgtf_non_normalise[k].gtf()[0].exposant()==_vcgtf_non_normalise[k+1].gtf()[0].exposant())	//Other format
				{
					lt=_vcgtf_non_normalise[k].gtf()[0].l()[0]+_vcgtf_non_normalise[k].gtf()[0].l()[1]+_vcgtf_non_normalise[k].gtf()[0].l()[2];

					if(_vcgtf_non_normalise[k].Lformat()=="Cart")
						m=(lt+1)*(lt+2)/2;
					else
						m=2*lt+1;

					s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
					s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

					for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
						s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

					k+=m;
				}

				else if(_vcgtf_non_normalise[k].Ltype()==_vcgtf_non_normalise[k+1].Ltype())   		//WFX
				{
					lt=_vcgtf_non_normalise[k].gtf()[0].l()[0]+_vcgtf_non_normalise[k].gtf()[0].l()[1]+_vcgtf_non_normalise[k].gtf()[0].l()[2];

					if(lt!=0)
					{
						m=(lt+1)*(lt+2)/2;
						q=k;
						save_alpha=_vcgtf_non_normalise[k].gtf()[0].exposant();

						do{
							s<<" "<<_vcgtf_non_normalise[q].Ltype()<<"\t"<<_vcgtf_non_normalise[q].numberOfFunctions()<<"  ";
							s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[q].FactorCoef()<<endl;

							for(int j=0; j<_vcgtf_non_normalise[q].numberOfFunctions(); j++)
								s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[q].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[q].coefficients()[j]<<endl;

							q++;
						}while(q+1<int(_vcgtf_non_normalise.size()) && _vcgtf_non_normalise[q].gtf()[0].exposant()!=save_alpha);
						//q=q-k+1;		don't work for h2otest.wfx
						q=q-k;
						//k+=q*(m-1)+1;      don't work for h2otest.wfx
						k+=q*m;
					}

					else
					{
						s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
						s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

						for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
							s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

						k++;
					}
				}

				else   			//Other case
				{
					s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
					s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

					for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
						s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

					k++;
				}
			}
			else   			//Last shell
			{
				s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
				s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

				for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
					s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

				k++;
			}
		}while(k<int(_vcgtf_non_normalise.size()) && i+1==_vcgtf_non_normalise[k].NumCenter());
		s<<endl;
	}

	s<<endl;

	s<<"[MO]"<<endl;

	int n;
	if(_alpha_and_beta)
		n=2;
	else
		n=1;

	for(int i=0; i<2; i++)
		for(size_t j=0; j<_coefficients[i].size(); j++)
		{
			s<<setprecision(6)<<" Ene= "<<_orbital_energy[i][j]<<endl;
			if(i==0)
				s<<" Spin= Alpha"<<endl;
			else
				s<<" Spin= Beta"<<endl;
			s<<setprecision(6)<<" Occup= "<<_occupation_number[i][j]/double(n)<<endl;
			s<<" Sym= unk"<<endl;

			for(size_t k=0; k<_coefficients[i][j].size(); k++)
			{
				s<<"\t"<<" "<<k+1<<"\t"<<"\t";
				s<<std::fixed<<setprecision(6)<<right<<setw(15)<<_coefficients[i][j][k]<<endl;
			}
		}

	s<<endl;
	s<<"[AO]"<<endl;
	s<<endl;

	s.close();
}

void Orbitals::Save_gab(string& tag)
{
	if(_mixte)
	{
		cout<<"Gabedit Format can't read mixte basis."<<endl;
		return;
	}

	if(int(_vcgtf_non_normalise.size())==_number_of_gtf)
	{
		cout<<"This option is nnot implemente."<<endl;
		return;
	}

	ofstream s;
	s.open(tag);

	s<<"[Gabedit Format] Cart"<<endl;

	s<<"[Atoms] AU"<<endl;
	for(int i=0; i<_number_of_atoms; i++)
	{
		s<<_symbol[i]<<" "<<i+1<<" "<<_atomic_numbers[i];
		for(int j=i*3; j<(i+1)*3; j++)
			s<<std::fixed<<setprecision(6)<<right<<"  "<<setw(10)<<_coordinates[j];
		s<<endl;
	}

	s<<"[Basis]"<<endl;

	int lt,q,m,k=0;
	double save_alpha;

	for(int i=0; i<_number_of_atoms; i++)
	{
		s<<"\t"<<_vcgtf_non_normalise[k].NumCenter()<<" "<<0<<endl;

		do{
			if(k==0)	//First shell
			{
				s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
				s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

				for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
					s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

				k++;
			}

			else if(k+1<int(_vcgtf_non_normalise.size()))	//Other shell
			{
				if(_vcgtf_non_normalise[k].Ltype()==_vcgtf_non_normalise[k+1].Ltype() && _vcgtf_non_normalise[k].gtf()[0].exposant()==_vcgtf_non_normalise[k+1].gtf()[0].exposant())	//Other format
				{
					lt=_vcgtf_non_normalise[k].gtf()[0].l()[0]+_vcgtf_non_normalise[k].gtf()[0].l()[1]+_vcgtf_non_normalise[k].gtf()[0].l()[2];

					if(_vcgtf_non_normalise[k].Lformat()=="Cart")
						m=(lt+1)*(lt+2)/2;
					else
						m=2*lt+1;

					s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
					s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

					for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
						s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

					k+=m;
				}

				else if(_vcgtf_non_normalise[k].Ltype()==_vcgtf_non_normalise[k+1].Ltype())   		//WFX
				{
					lt=_vcgtf_non_normalise[k].gtf()[0].l()[0]+_vcgtf_non_normalise[k].gtf()[0].l()[1]+_vcgtf_non_normalise[k].gtf()[0].l()[2];

					if(lt!=0)
					{
						m=(lt+1)*(lt+2)/2;
						q=k;
						save_alpha=_vcgtf_non_normalise[k].gtf()[0].exposant();

						do{
							s<<" "<<_vcgtf_non_normalise[q].Ltype()<<"\t"<<_vcgtf_non_normalise[q].numberOfFunctions()<<"  ";
							s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[q].FactorCoef()<<endl;

							for(int j=0; j<_vcgtf_non_normalise[q].numberOfFunctions(); j++)
								s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[q].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[q].coefficients()[j]<<endl;

							q++;
						}while(q+1<int(_vcgtf_non_normalise.size()) && _vcgtf_non_normalise[q].gtf()[0].exposant()!=save_alpha);
						//q=q-k+1;		don't work for h2otest.wfx
						q=q-k;
						//k+=q*(m-1)+1;      don't work for h2otest.wfx
						k+=q*m;
					}

					else
					{
						s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
						s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

						for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
							s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

						k++;
					}
				}

				else   			//Other case
				{
					s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
					s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

					for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
						s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

					k++;
				}
			}
			else   			//Last shell
			{
				s<<" "<<_vcgtf_non_normalise[k].Ltype()<<"\t"<<_vcgtf_non_normalise[k].numberOfFunctions()<<"  ";
				s<<std::fixed<<setprecision(1)<<_vcgtf_non_normalise[k].FactorCoef()<<endl;

				for(int j=0; j<_vcgtf_non_normalise[k].numberOfFunctions(); j++)
					s<<std::fixed<<setprecision(6)<<right<<setw(6)<<" "<<setw(15)<<_vcgtf_non_normalise[k].gtf()[j].exposant()<<"  "<<setw(15)<<_vcgtf_non_normalise[k].coefficients()[j]<<endl;

				k++;
			}
		}while(k<int(_vcgtf_non_normalise.size()) && i+1==_vcgtf_non_normalise[k].NumCenter());
		s<<endl;
	}

	s<<endl;

	s<<"[MO]"<<endl;

	int n;
	if(_alpha_and_beta)
		n=2;
	else
		n=1;

	for(int i=0; i<2; i++)
		for(size_t j=0; j<_coefficients[i].size(); j++)
		{
			s<<setprecision(6)<<" Ene= "<<_orbital_energy[i][j]<<endl;
			if(i==0)
				s<<" Spin= Alpha"<<endl;
			else
				s<<" Spin= Beta"<<endl;
			s<<setprecision(6)<<" Occup= "<<_occupation_number[i][j]/double(n)<<endl;
			s<<" Sym= unk"<<endl;

			for(size_t k=0; k<_coefficients[i][j].size(); k++)
			{
				s<<"\t"<<" "<<k+1<<"\t"<<"\t";
				s<<std::fixed<<setprecision(6)<<right<<setw(15)<<_coefficients[i][j][k]<<endl;
			}
		}

	s<<endl;
	s<<"[AO]"<<endl;
	s<<endl;

	s.close();
}

void Orbitals::Sorting()
{
	int k=0,q=0,pos,lt,m=0;

	if(_number_of_gtf==int(_vcgtf.size()))
	{
		for(int i=0; i<_number_of_atoms; i++)
		{
			do{
				if(k+1<int(_vcgtf.size()))
				{
					if((_vcgtf[k].Ltype()=="s" || _vcgtf[k].Ltype()=="S") || (_vcgtf[k+1].Ltype()=="s" || _vcgtf[k+1].Ltype()=="S"))
						q=0;

					else if(_vcgtf[k].Ltype()==_vcgtf[k+1].Ltype())
					{
						lt=_vcgtf[k].gtf()[0].l()[0]+_vcgtf[k].gtf()[0].l()[1]+_vcgtf[k].gtf()[0].l()[2];
						m=(lt+1)*(lt+2)/2;
						q=k;

						do{
							q++;
						}while(_vcgtf[q].Ltype()==_vcgtf[k].Ltype());

						q=q-k;
						k+=q+1;
					}

					if(q>m)
					{
						pos=k-q-1;
						q=pos;

						for(int ind=0; ind<m-1; ind++)
						{
							q++;
							swap(_vcgtf[q], _vcgtf[q+m-1]);
							q+=m;
						}

						q=pos+m-1;

						swap(_vcgtf[q], _vcgtf[q+m+1]);

						q=0;
						k++;
					}

					else
					{
						k++;
						continue;
					}
				}
				else
				{
					k++;
					continue;
				}
			}while(k<int(_vcgtf.size()) && i+1==_vcgtf[k].NumCenter());
		}
		DenormaliseAllBasis();
	}
}
