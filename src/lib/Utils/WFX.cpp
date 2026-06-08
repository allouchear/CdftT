#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <Utils/WFX.h>
#include <Utils/Utils.h>


//****************************************//
//************ MOPC class ****************//
//****************************************//

MOPC::MOPC()
{
    _MO_Number=0;
    _Coefficients.resize(0);
}

void MOPC::push_back(int a, std::vector<double> b)
{
    _MO_Number=a;
    _Coefficients=b;
}


//****************************************//
//************ NCEG class ****************//
//****************************************//

NCEG::NCEG()
{
    _symbol="None";
    _gradient.resize(0);
}

void NCEG::push_back(std::string a, std::vector<double> b)
{
    _symbol=a;
    _gradient=b;
}


//****************************************//
//*********** AEDF class *****************//
//****************************************//

AEDF::AEDF()
{
    _Number_of_EDF_Primitives=0;
    _EDF_Primitives_Centers.resize(0);
    _EDF_Primitives_Types.resize(0);
    _EDF_Primitives_Exponents.resize(0);
    _EDF_Primitives_Coefficients.resize(0);
}

void AEDF::push_back(int a, std::vector<int> b, std::vector<int> c, std::vector<double> d, std::vector<double> e)
{
    _Number_of_EDF_Primitives=a;
    _EDF_Primitives_Centers=b;
    _EDF_Primitives_Types=c;
    _EDF_Primitives_Exponents=d;
    _EDF_Primitives_Coefficients=e;
}

//****************************************//
//*********** WFX class ******************//
//****************************************//

WFX::WFX()
{
    _alpha_and_beta = false;
    _Title="None";                                                        
    _Keywords="None";                                                    
    _Number_of_Nuclei=0;                                                
    _Number_of_Primitives=0;    
    _Number_of_Occupied_Molecular_Orbital=0;    
    _Number_of_Perturbations=0;                                    
    _Nuclear_Names.resize(0);    
    _Atomic_Number.resize(0);
    _Nuclear_Charges.resize(0);    
    _Nuclear_Cartesian_Coordinates.resize(0);
    _Net_Charge=0.0;                                                    
    _Number_of_Electrons=0;                                            
    _Number_of_Alpha_Electrons=0;                                        
    _Number_of_Beta_Electrons=0;                                        
    _Electronic_Spin_Multiplicity=0;
    _Model="None";                            
    _Primitive_Centers.resize(0);                                        
    _Primitive_Types.resize(0);
    _Lxyz =std::vector<std::vector<int>>();                                        
    _Primitive_Exponents.resize(0);                                    
    _Molecular_Orbital_Occupation_Numbers.resize(0);                
    _Molecular_Orbital_Energies.resize(0);
    _Molecular_Orbital_Spin_Types.resize(0);                    
    _Molecular_Orbital_Primitive_Coefficients.resize(0);
    _Energy=0.0;
    _Virial_Ratio=0.0;
    _Nuclear_Cartesian_Energy_Gradients.resize(0);
    _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei=0.0;
    _Full_Virial_Ratio=0.0;
    _Number_of_Core_Electrons=0;
    _Additionnal_Electron_Density_Function=AEDF();
}

WFX::WFX(std::ifstream& file)
{
    _alpha_and_beta = false;
    read_file_wfx(file);
    std::vector<int> l(3);
    _Lxyz.resize(_Number_of_Primitives, l);
    int i;
    for(i=0; i<_Number_of_Primitives; i++)
        _Lxyz[i]=setLxyz(_Primitive_Types[i]);

}

std::vector<int> WFX::read_one_block_int(std::ifstream& f, std::string b, bool r, int nint)
{
    int data;
    std::vector<int> block_of_data(0);
    long int pos;
    int i,n;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return std::vector<int> (0);
    }
    f.seekg(pos);
    for(i=0; i<nint; i++)
    {
        f>>data;
        block_of_data.push_back(data);
    }
    return block_of_data;
}

std::vector<double> WFX::read_one_block_real(std::ifstream& f, std::string b, bool r)
{
    double data;
    std::vector<double> block_of_data(0);
    long int pos;
    int i,n;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return std::vector<double> (0);
    }
    f.seekg(pos);
    for(i=0; i<n; i++)
    {
        f>>data;
        block_of_data.push_back(data);
    }

    return block_of_data;
}

std::vector<double> WFX::read_one_block_real(std::ifstream& f, std::string b, bool r, int ndouble)
{
    double data;
    std::vector<double> block_of_data(0);
    long int pos;
    int i,n;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return std::vector<double> (0);
    }
    f.seekg(pos);
    for(i=0; i<ndouble; i++)
    {
        f>>data;
        block_of_data.push_back(data);
    }

    return block_of_data;
}

std::vector<std::string> WFX::read_one_block_string(std::ifstream& f, std::string b, bool r)
{
    std::string data;
    std::vector<std::string> block_of_data(0);
    long int pos;
    int i,n;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return std::vector<std::string> (0);
    }
    f.seekg(pos);
    for(i=0; i<n; i++)
    {
        getline(f,data);
        block_of_data.push_back(data);
    }

    return block_of_data;
}

int WFX::read_int(std::ifstream& f, std::string b, bool r)
{
    int data,i,n;
    long int pos;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return 0;
    }
    f.seekg(pos);
    for(i=0; i<n; i++)
        f>>data;

    return data;
}

double WFX::read_real(std::ifstream& f, std::string b, bool r)
{
    double data;
    int i,n;
    long int pos;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return 0.0;
    }
    f.seekg(pos);
    for(i=0; i<n; i++)
        f>>data;

    return data;
}

std::string WFX::read_string(std::ifstream& f, std::string b, bool r)
{
    std::string data;
    int i,n;
    long int pos;
    pos=LocaliseBlock(f, n, b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return "None";
    }
    f.seekg(pos);
    for(i=0; i<n; i++)
        getline(f,data);

    return data;
}

std::vector<MOPC> WFX::read_MOPC_block(std::ifstream& f, std::string b, bool r)
{
    int i,j,n,nMO,ndat;
    double c;
    std::vector<double> cdat(0);
    MOPC data;
    std::vector<MOPC> block_of_data(0);
    std::string MO="MO Number";
    std::string p;
    std::string test;
    long int pos;
    pos=LocaliseMO(f,n,nMO,b);

    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return std::vector<MOPC> (0);
    }

    f.clear();
    f.seekg(pos);
    n=_Number_of_Primitives;

    for(i=0; i<nMO/2; i++)
    {
        getline(f,p);
        f>>ndat;
        getline(f,p);
        getline(f,p);

        for(j=0; j<n; j++)
        {
            f>>c;
            cdat.push_back(c);
        }
        data.push_back(ndat,cdat);
        block_of_data.push_back(data);
        cdat.resize(0);
        getline(f,p);
    }

    return block_of_data;
}

std::vector<NCEG> WFX::read_NCEG_block(std::ifstream& f, std::string b, bool r)
{
    int i,j,n;
    long int pos;

    std::string s;
    double g;
    std::vector<double> gdat(0);
    NCEG data;
    std::vector<NCEG> block_of_data(0);

    pos=LocaliseBlock(f,n,b);
    if(pos==-1)
    {
        std::cout<<b+" : data not found"<<std::endl;
        if(r)
        {
            std::stringstream errorMessage;
            errorMessage << b << " : data not found. Data required, please check your file." << std::endl;
            print_error(errorMessage.str());

            std::exit(1);
        }
        else
            return std::vector<NCEG> (0);
    }

    f.seekg(pos);
    for(i=0; i<n/2; i++)
    {
        f>>s;
        for(j=0; j<3; j++)
        {
            f>>g;
            gdat.push_back(g);
        }
        data.push_back(s,gdat);
        block_of_data.push_back(data);
    }

    return block_of_data;
}

AEDF WFX::read_AEDF_block(std::ifstream& f, std::string b, bool r)
{
    std::string snp="Number of EDF Primitives";
    std::string spc="EDF Primitive Centers";
    std::string spt="EDF Primitive Types";
    std::string spe="EDF Primitive Exponents";
    std::string spcoef="EDF Primitive Coefficients";
    int n=read_int(f, snp, r);
    std::vector<int> pc=read_one_block_int(f, spc, r, n);
    std::vector<int> pt=read_one_block_int(f, spt, r, n);
    std::vector<double> pe=read_one_block_real(f, spe, r);
    std::vector<double> pcoef=read_one_block_real(f, spcoef, r);

    AEDF block_of_data;
    block_of_data.push_back(n, pc, pt, pe, pcoef);

    return block_of_data;
}

long int LocaliseBlock(std::ifstream& f, int& n, std::string b)
{
    n=0;
    long int position;
    f.clear();
    f.seekg(0,f.beg);
    std::string test;
    bool ok=false;
    while(!f.eof())
    {    
        getline(f, test);
        if(test.find(b)!=std::string::npos)
        {
            ok=true;
            position=f.tellg();
            break;
        }
    }

    if(!ok) 
        return -1;

    while(!f.eof())
    {    
        getline(f, test);
        if(test.find(b)!=std::string::npos)
        {
            ok=true;
            break;
        }
        else n++;
    }
    if(!ok) 
        return -1;    

    return position;
}

long int LocaliseMO(std::ifstream& f, int& n, int& n_MO, std::string b)
{
    n=0;
    n_MO=0;
    std::string m="MO Number";
    long int position;
    position=LocaliseBlock(f,n,b);
    f.clear();
    f.seekg(position);
    std::string test;

    while(!f.eof())
    {
        getline(f,test);
        if(test.find(m)!=std::string::npos)
            n_MO++;
        else if(test.find(b)!=std::string::npos)
            break;
    }
    if(n==0 || n_MO==0)
        position=-1;

    return position;
}

void WFX::read_file_wfx(std::ifstream& file)
{
    std::string a="Title";
    std::string b="Keywords";
    std::string c="Number of Nuclei";
    std::string d="Number of Primitives";
    std::string e="Number of Occupied Molecular Orbitals";
    std::string f="Number of Perturbations";
    std::string g="Nuclear Names";
    std::string h="Atomic Number";
    std::string i="Nuclear Charges";
    std::string j="Nuclear Cartesian Coordinates";
    std::string k="Net Charge";
    std::string l="Number of Electrons";
    std::string m="Number of Alpha Electrons";
    std::string n="Number of Beta Electrons";
    std::string o="Electronic Spin Multiplicity";
    std::string p="Model";
    std::string q="Primitive Centers";
    std::string r="Primitive Types";
    std::string s="Primitive Exponents";
    std::string t="Molecular Orbital Occupation Numbers";
    std::string u="Molecular Orbital Energies";
    std::string v="Molecular Orbital Spin Types";
    std::string w="Molecular Orbital Primitive Coefficients";
    std::string x="Energy = T + Vne + Vee + Vnn";
    std::string y="Virial Ratio (-V/T)";
    std::string z="Nuclear Cartesian Energy Gradients";
    std::string aa="Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W";
    std::string bb="Full Virial Ratio, -(V - W)/T";
    std::string cc="Number of Core Electrons";
    std::string dd="Additional Electron Density Function (EDF)";

    _Molecular_Orbital_Occupation_Numbers=std::vector<std::vector<double>> (2);
    _Molecular_Orbital_Energies=std::vector<std::vector<double>> (2);
    _Molecular_Orbital_Primitive_Coefficients=std::vector<std::vector<MOPC>> (2);

    _Title=read_string(file, a, true);
    _Keywords=read_string(file, b, true);
    _Number_of_Nuclei=read_int(file, c, true);
    _Number_of_Primitives=read_int(file, d, true);
    _Number_of_Occupied_Molecular_Orbital=read_int(file, e, true);
    _Number_of_Perturbations=read_int(file,f, true);
    _Nuclear_Names=read_one_block_string(file,g, true);
    _Atomic_Number=read_one_block_int(file,h, false,_Number_of_Nuclei);
    _Nuclear_Charges=read_one_block_real(file,i, true);
    _Nuclear_Cartesian_Coordinates=read_one_block_real(file,j, true, _Number_of_Nuclei*3);
    _Net_Charge=read_real(file,k, true); 
    _Number_of_Electrons=read_int(file,l, true);
    _Number_of_Alpha_Electrons=read_int(file,m, true);
    _Number_of_Beta_Electrons=read_int(file,n, true);
    _Electronic_Spin_Multiplicity=read_int(file,o, false);
    _Model=read_string(file,p, false);
    _Primitive_Centers=read_one_block_int(file,q, true, _Number_of_Primitives);
    _Primitive_Types=read_one_block_int(file,r, true, _Number_of_Primitives);
    _Primitive_Exponents=read_one_block_real(file,s, true, _Number_of_Primitives);
    _Molecular_Orbital_Occupation_Numbers[0]=read_one_block_real(file,t, true);
    _Molecular_Orbital_Energies[0]=read_one_block_real(file,u, true);
    _Molecular_Orbital_Spin_Types=read_one_block_string(file,v, true);
    _Molecular_Orbital_Primitive_Coefficients[0]=read_MOPC_block(file,w, true);
    _Energy=read_real(file,x, true);
    _Virial_Ratio=read_real(file,y, true);
    _Nuclear_Cartesian_Energy_Gradients=read_NCEG_block(file,z, false);
    _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei=read_real(file,aa, false);
    _Full_Virial_Ratio=read_real(file,bb, false);
    _Number_of_Core_Electrons=read_int(file,cc, false);
    _Additionnal_Electron_Density_Function=read_AEDF_block(file,dd, false);

    if(_Molecular_Orbital_Spin_Types[0]==" Alpha and Beta")
    {
        _alpha_and_beta = true;
        _Molecular_Orbital_Occupation_Numbers[1]=_Molecular_Orbital_Occupation_Numbers[0];
        _Molecular_Orbital_Energies[1]=_Molecular_Orbital_Energies[0];
        _Molecular_Orbital_Primitive_Coefficients[1]=_Molecular_Orbital_Primitive_Coefficients[0];
    }
    
    else
    {    
        _Molecular_Orbital_Occupation_Numbers[1]=std::vector<double> (0);
        _Molecular_Orbital_Energies[1]=std::vector<double> (0);
        _Molecular_Orbital_Primitive_Coefficients[1]=std::vector<MOPC> (0);

        for(size_t np=0; np<_Molecular_Orbital_Spin_Types.size(); np++)
        {
            if(_Molecular_Orbital_Spin_Types[np]==" Beta")
            {
                _Molecular_Orbital_Occupation_Numbers[1].push_back(_Molecular_Orbital_Occupation_Numbers[0][np]);
                _Molecular_Orbital_Energies[1].push_back(_Molecular_Orbital_Energies[0][np]);
                _Molecular_Orbital_Primitive_Coefficients[1].push_back(_Molecular_Orbital_Primitive_Coefficients[0][np]);
            }
        }

        int nb=_Molecular_Orbital_Primitive_Coefficients[1].size();
        int nc=0;

        while(nc<nb)
        {
            _Molecular_Orbital_Occupation_Numbers[0].pop_back();
            _Molecular_Orbital_Energies[0].pop_back();
            _Molecular_Orbital_Primitive_Coefficients[0].pop_back();
            nc++;
        }
    }
    for(int i=0; i<_Number_of_Nuclei; i++)
    {
        _Nuclear_Names[i].erase(0,1);
        std::stringstream s;
        s<<i+1;
        std::string p=s.str();
        if(_Nuclear_Names[i].find(p)!=std::string::npos)
            for(size_t j=0; j<p.size(); j++)
                _Nuclear_Names[i].pop_back();
    }
}

void WFX::write_one_block_int(std::ofstream& f, std::vector<int> v, std::string b, bool r, int nint)
{
    if(!r && v.size()==0)
        return;
    f<<b<<std::endl;
    int i;
    for(i=0; i<int(v.size()); i++)
        f<<" "<<v[i];
    b.insert(1,"\\");
    f<<std::endl<<b<<std::endl;;
}

void WFX::write_one_block_real(std::ofstream& f, std::vector<double> v, std::string b, bool r)
{
    if(!r && v.size()==0)
        return;
    f<<b<<std::endl;
    int i;
    for(i=0; i<int(v.size()); i++)
        f<<" "<<v[i]<<"\t";
    b.insert(1,"\\");
    f<<std::endl<<b<<std::endl;
}

void WFX::write_one_block_real(std::ofstream& f, std::vector<double> v, std::string b, bool r, int ndouble)
{
    if(!r && v.size()==0)
        return;
    f<<b<<std::endl;
    int i;
    for(i=0; i<ndouble; i++)
        f<<" "<<v[i]<<"\t";
    b.insert(1,"\\");
    f<<std::endl<<b<<std::endl;
}

void WFX::write_one_block_string(std::ofstream& f, std::vector<std::string> v, std::string b, bool r)
{
    if(!r && v.size()==0)
        return;
    f<<b<<std::endl;
    int i;
    for(i=0; i<int(v.size()); i++)
        f<<v[i]<<std::endl;
    b.insert(1,"\\");
    f<<b<<std::endl;
}

void WFX::write_one_matrix_real(std::ofstream& f, std::vector<std::vector<double>> v, std::string b, bool r)
{
    if(!r && v.size()==0)
        return;
    f<<b<<std::endl;
    int i,j;
    for(i=0; i<int(v.size()); i++)
        for(j=0; j<int(v[i].size()); j++)
            f<<" "<<v[i][j]<<"\t";
    b.insert(1,"\\");
    f<<std::endl<<b<<std::endl;
}

void WFX::write_int(std::ofstream& f, int i, std::string b, bool r)
{
    if(!r && i==0)
        return;
    f<<b<<std::endl;
    f<<" "<<i<<std::endl;
    b.insert(1,"\\");
    f<<b<<std::endl;
}

void WFX::write_real(std::ofstream& f, double d, std::string b, bool r)
{
    if(!r && d==0.0)
        return;
    f<<b<<std::endl;
    f<<" "<<d<<std::endl;
    b.insert(1,"\\");
    f<<b<<std::endl;
}

void WFX::write_string(std::ofstream& f, std::string s, std::string b, bool r)
{
    if(!r && s=="None")
        return;
    f<<b<<std::endl;
    f<<" "<<s<<std::endl;
    b.insert(1,"\\");
    f<<b<<std::endl;
}

void WFX::write_MOPC_block(std::ofstream& f, std::vector<std::vector<MOPC>> v, bool r)
{
    if(!r && v.size()==0)
        return;
    std::string b="<Molecular Orbital Primitive Coefficients>";
    std::string m1="<MO Number>";
    std::string m2=m1;
    m2.insert(1,"\\");
    f<<b<<std::endl;
    int i,j,k;
    int n=_Number_of_Primitives;
    std::vector<int> nMO (2);
    
    nMO[0]=v[0].size();
    nMO[1]=v[1].size();

    for(i=0; i<2; i++)
    {
        for(j=0; j<nMO[i]; j++)
        {
            f<<m1<<std::endl;
            f<<" "<<v[i][j].MO_Number()<<std::endl;
            f<<m2<<std::endl;
            for(k=0; k<n; k++)
            {
                f<<" "<<v[i][j].Coefficients()[k];
            }
            f<<std::endl;
        }
        if(_alpha_and_beta)
            break;
    }
    b.insert(1,"\\");
    f<<b<<std::endl;
}

void WFX::write_NCEG_block(std::ofstream& f, std::vector<NCEG> v, bool r)
{
    if(!r && v.size()==0)
        return;
    std::string b="<Nuclear Cartesian Energy Gradients>";
    f<<b<<std::endl;
    int i,j;
    for(i=0; i<int(v.size()); i++)
    {
        f<<" "<<v[i].symbol();
        for(j=0; j<3; j++)
            f<<" "<<v[i].gradient()[i];
    }
    b.insert(1,"\\");
    f<<std::endl<<b<<std::endl;
}

void WFX::write_AEDF_block(std::ofstream& f, AEDF a, bool r)
{
    if(!r && a.EDF_Primitives_Coefficients().size()==0)
        return;
    std::string b="<Additional Electron Density Function (EDF)>";
    std::string ub1="<Number of EDF Primitives>";
    std::string ub2="<EDF Primitive Centers>";
    std::string ub3="<EDF Primitive Types>";
    std::string ub4="<EDF Primitive Exponents>";
    std::string ub5="<EDF Primitive Coefficients>";
    f<<b<<std::endl;
    write_int(f, a.Number_of_EDF_Primitives(), ub1, r);
    write_one_block_int(f, a.EDF_Primitives_Centers(), ub2, r, a.Number_of_EDF_Primitives());
    write_one_block_int(f, a.EDF_Primitives_Types(), ub3, r, a.Number_of_EDF_Primitives());
    write_one_block_real(f, a.EDF_Primitives_Exponents(), ub4, r);
    write_one_block_real(f, a.EDF_Primitives_Coefficients(), ub5, r);
    b.insert(1,"\\");
    f<<b<<std::endl;
}

void WFX::write_file_wfx(std::ofstream& file)
{
    file<<std::scientific;
    file<<std::setprecision(15);
    //std::cout<<std::left<<std::setw(20);


    std::string a="<Title>";
    std::string b="<Keywords>";
    std::string c="<Number of Nuclei>";
    std::string d="<Number of Primitives>";
    std::string e="<Number of Occupied Molecular Orbitals>";
    std::string f="<Number of Perturbations>";
    std::string g="<Nuclear Names>";
    std::string h="<Atomic Number>";
    std::string i="<Nuclear Charges>";
    std::string j="<Nuclear Cartesian Coordinates>";
    std::string k="<Net Charge>";
    std::string l="<Number of Electrons>";
    std::string m="<Number of Alpha Electrons>";
    std::string n="<Number of Beta Electrons>";
    std::string o="<Electronic Spin Multiplicity>";
    std::string p="<Model>";
    std::string q="<Primitive Centers>";
    std::string r="<Primitive Types>";
    std::string s="<Primitive Exponents>";
    std::string t="<Molecular Orbital Occupation Numbers>";
    std::string u="<Molecular Orbital Energies>";
    std::string v="<Molecular Orbital Spin Types>";
    std::string w="<Energy = T + Vne + Vee + Vnn>";
    std::string x="<Virial Ratio (-V/T)>";
    std::string y="<Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W>";
    std::string z="<Full Virial Ratio, -(V - W)/T>";
    std::string aa="<Number of Core Electrons>";

    write_string(file, _Title, a, true);
    write_string(file, _Keywords, b, true);
    write_int(file, _Number_of_Nuclei, c, true);
    write_int(file, _Number_of_Primitives, d, true);
    write_int(file, _Number_of_Occupied_Molecular_Orbital, e, true);
    write_int(file, _Number_of_Perturbations, f, true);
    write_one_block_string(file, _Nuclear_Names, g, true);
    write_one_block_int(file, _Atomic_Number, h, false, _Number_of_Nuclei);
    write_one_block_real(file, _Nuclear_Charges, i, true);
    write_one_block_real(file, _Nuclear_Cartesian_Coordinates, j, true, _Number_of_Nuclei*3);
    write_real(file, _Net_Charge, k, true); 
    write_int(file, _Number_of_Electrons, l, true);
    write_int(file, _Number_of_Alpha_Electrons, m, true);
    write_int(file, _Number_of_Beta_Electrons, n, true);
    write_int(file, _Electronic_Spin_Multiplicity, o, false);
    write_string(file, _Model, p, false);
    write_one_block_int(file, _Primitive_Centers, q, true, _Number_of_Primitives);
    write_one_block_int(file, _Primitive_Types, r, true, _Number_of_Primitives);
    write_one_block_real(file, _Primitive_Exponents, s, true);
    write_one_matrix_real(file, _Molecular_Orbital_Occupation_Numbers, t, true);        //probleme
    write_one_matrix_real(file, _Molecular_Orbital_Energies, u, true);        // probleme
    write_one_block_string(file, _Molecular_Orbital_Spin_Types, v, true);
    write_MOPC_block(file, _Molecular_Orbital_Primitive_Coefficients, true); //probleme
    write_real(file, _Energy, w, true);
    write_real(file, _Virial_Ratio, x, true);
    write_NCEG_block(file, _Nuclear_Cartesian_Energy_Gradients, false);
    write_real(file, _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei, y, false);
    write_real(file, _Full_Virial_Ratio, z, false);
    write_int(file, _Number_of_Core_Electrons, aa, false);
    write_AEDF_block(file, _Additionnal_Electron_Density_Function, false);
}

std::vector<int> setLxyz(int iType)
{
    std::vector<int> l (3);
    l[0]=l[1]=l[2]=0;
    if(iType==1) 
        return l; // 1S
    else if(iType==2) 
        l[0]=1; // 2 PX
    else if(iType==3) 
        l[1]=1;// 2 PY
    else if(iType==4) 
        l[2]=1;// 4 PZ
    else if(iType==5)
        l[0]=2; // 5 DXX
    else if(iType==6)
        l[1]=2; // 6 DYY
    else if(iType==7) 
        l[2]=2; // 7 DZZ
    else if(iType==8)
    {
        l[0]=1;  
        l[1]=1;
    } // 8 DXY
    else if(iType==9)
    {
        l[0]=1; 
        l[2]=1;
    } // 9 DXZ
    else if(iType==10)
    {
        l[1]=1;
        l[2]=1;
    } // 10 DYZ
    else if(iType==11)
        l[0]=3; // 11 FXXX
    else if(iType==12)
        l[1]=3; // 12 FYYY
    else if(iType==13)
        l[2]=3; // 13 FZZZ
    else if(iType==14)
    {
        l[0]=2; 
        l[1]=1; // 14 FXXY
    }
    else if(iType==15)
    {
        l[0]=2;
        l[2]=1;
    } // 15 FXXZ
    else if(iType==16)
    {
        l[1]=2;
        l[2]=1;
    } // 16 FYYZ
    else if(iType==17)
    {
        l[0]=1;
        l[1]=2;
    } // 17 FXYY
    else if(iType==18)
    {
        l[0]=1;
        l[2]=2;
    } // 18 FXZZ
    else if(iType==19)
    {
        l[1]=1;
        l[2]=2;
    } // 19 FYZZ
    else if(iType==20)
    {
        l[0]=1;
        l[1]=1;
        l[2]=1;
    } // 20 FXYZ
    else if(iType==21)
        l[0]=4; // 21 GXXXX
    else if(iType==22)
        l[1]=4; // 22 GYYYY
    else if(iType==23)
        l[2]=4; // 23 GZZZZ
    else if(iType==24)
    {
        l[0]=3;
        l[1]=1;
        l[2]=0;
    } // 24 GXXXY
    else if(iType==25)
    {
        l[0]=3;
        l[1]=0;
        l[2]=1;
    } // 25 GXXXZ
    else if(iType==26)
    {
        l[0]=1;
        l[1]=3;
        l[2]=0;
    } // 26 GXYYY
    else if(iType==27)
    {
        l[0]=0;
        l[1]=3;
        l[2]=1;
    } // 27 GYYYZ
    else if(iType==28)
    {
        l[0]=1;
        l[1]=0;
        l[2]=3;
    } // 28 GXZZZ
    else if(iType==29)
    {
        l[0]=0;
        l[1]=1;
        l[2]=3;
    } // 29 GYZZZ
    else if(iType==30) 
    { 
        l[0]=2;
        l[1]=2;
        l[2]=0;
    } // 30 GXXYY
    else if(iType==31)
    {
        l[0]=2;
        l[1]=0;
        l[2]=2;
    } // 31 GXXZZ
    else if(iType==32)
    {
        l[0]=0;
        l[1]=2;
        l[2]=2;
    } // 32 GYYZZ
    else if(iType==33)
    {
        l[0]=2;
        l[1]=1;
        l[2]=1;
    } // 33 GXXYZ
    else if(iType==34)
    {
        l[0]=1;
        l[1]=2;
        l[2]=1;
    } // 34 GXYYZ
    else if(iType==35)
    {
        l[0]=1;
        l[1]=1;
        l[2]=2;
    } // 35 GXYZZ
    else
    {
        int it=35;
        int L, ix,iy;
        for(L=5;L<=30;L++)
            for(ix=0;ix<L;ix++)
                for(iy=0;iy<=L-ix;iy++)
                {
                    it++;
                    if(it==iType)
                    {
                        l[0] =ix; l[1]=iy; l[2]=L-ix-iy;
                        return l;;
                    }
                }
    }
    return l;
}
