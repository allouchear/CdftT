#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <Common/Constants.h>
#include <Utils/MOLDENGAB.h>


MOLDENGAB::MOLDENGAB() :
    _symbol(),
    _atomic_number(),
    _coord(),
    _shell_types(),
    _L_types(),
    _exposants(),
    _number_of_gtf(),
    _num_center(),
    _cgtf_coefs(),
    _factor_coefs(),
    _MO_energy(),
    _MO_coefs(),
    _occupation(),
    _spin_types(),
    _coord_type("None"),
    _basis_or_gto("None"),
    _number_of_atoms(0),
    _number_of_MO_coefs(0),
    _number_of_MO(0),
    _alpha_and_beta(true),
    _n_at_basis(),
    _alpha_occupation(),
    _beta_occupation(),
    _alpha_MO_coefs(),
    _beta_MO_coefs(),
    _alpha_energies(),
    _beta_energies(),
    _format("None"),
    _cart_sphe("none"),
    _mixte(false)
{ }

MOLDENGAB::MOLDENGAB(std::ifstream& file) :
    _symbol(),
    _atomic_number(),
    _coord(),
    _shell_types(),
    _L_types(),
    _exposants(),
    _number_of_gtf(),
    _num_center(),
    _cgtf_coefs(),
    _factor_coefs(),
    _MO_energy(),
    _MO_coefs(),
    _occupation(),
    _spin_types(),
    _coord_type("None"),
    _basis_or_gto("None"),
    _number_of_atoms(0),
    _number_of_MO_coefs(0),
    _number_of_MO(0),
    _alpha_and_beta(true),
    _n_at_basis(),
    _alpha_occupation(),
    _beta_occupation(),
    _alpha_MO_coefs(),
    _beta_MO_coefs(),
    _alpha_energies(),
    _beta_energies(),
    _format("None"),
    _cart_sphe("none"),
    _mixte(false)
{

    file.clear();
    file.seekg(0,file.beg);
    std::string p;
    getline(file,p);

    if(p.find("[Molden Format]") != std::string::npos)
    {
        _format = "molden";
        _basis_or_gto = "[GTO]";
    }
    else if(p.find("[Gabedit Format]") != std::string::npos)
    {
        _format = "gabedit";
        _basis_or_gto = "[Basis]";

        if(p.find("Sphe") != std::string::npos)
        {
            _cart_sphe = "sphe";
        }
        else if(p.find("Cart") != std::string::npos)
        {
            _cart_sphe = "cart";
        }
        else
        {
            std::cout << "Error, can't recognize data format (sphe/cart)." << std::endl;
            std::cout << "Please check your file." << std::endl;

            exit(1);
        }
    }
    else
    {
        std::cout << "Error, can't recognize file format." << std::endl;
        std::cout << "Please check your file." << std::endl;

        exit(1);
    }

    read_atom_data(file);

    _number_of_atoms = _atomic_number.size();
    _n_at_basis = std::vector<int>(_number_of_atoms);
    
    read_basis_data(file);
    read_MO_data(file);

    if(_coord_type == "Angs")
    {
        for(int i = 0; i < _number_of_atoms; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                _coord[i][j] *= Constants::ANGSTROM_TO_BOHR_RADIUS;
            }
        }
    }

    if(_alpha_and_beta)
    {
        _alpha_energies = _beta_energies = _MO_energy;
        _alpha_occupation = _beta_occupation = _occupation;
        _alpha_MO_coefs = _beta_MO_coefs = _MO_coefs;
    }
    else
    {
        for(int i = 0; i < _number_of_MO; ++i)
        {
            if(_spin_types[i] == "Alpha")
            {
                _alpha_energies.push_back(_MO_energy[i]);
                _alpha_occupation.push_back(_occupation[i]);
                _alpha_MO_coefs.push_back(_MO_coefs[i]);
            }
            
            else if(_spin_types[i] == "Beta")
            {
                _beta_energies.push_back(_MO_energy[i]);
                _beta_occupation.push_back(_occupation[i]);
                _beta_MO_coefs.push_back(_MO_coefs[i]);
            }
        }
    }
}

void MOLDENGAB::read_atom_data(std::ifstream& f)
{
    std::string p;
    int an;
    std::vector<double> c (3);
    long int pos=LocaliseDataMolGabBefore(f,"[Atoms]");

    if(pos==-1)
    {
        std::cout<<"Atoms data not found"<<std::endl;
        std::cout<<"Data required, please check your file"<<std::endl;
        exit(1);
    }

    f.seekg(pos);
    getline(f,p);
    std::stringstream t(p);
    t>>p;
    t>>_coord_type;
    
    getline(f,p);

    do{
        std::stringstream s(p);
        s>>p;    
        _symbol.push_back(p);
        s>>p;
        s>>an;
        _atomic_number.push_back(an);
        s>>c[0];
        s>>c[1];
        s>>c[2];
        _coord.push_back(c);
        getline(f,p);
    }while(p.find("[")==std::string::npos);
}

void MOLDENGAB::read_basis_data(std::ifstream& f)
{
    std::string p;

    long int pos=LocaliseDataMolGab(f, _basis_or_gto);

    if(pos==-1)
    {
        std::cout<<"Basis (GTO) data not found"<<std::endl;
        std::cout<<"Data required, please check your file"<<std::endl;
        exit(1);
    }

    f.seekg(pos);
    int n= _atomic_number.size();

    for(int i=0; i<n; i++)
        read_one_basis_data(f);

    int m=_shell_types.size();
    
    long int posd, posf, posg;

    if(_format=="molden")
    {
        posd = LocaliseDataMolGab(f, "[5D]");
        posf = LocaliseDataMolGab(f, "[7F]");
        posg = LocaliseDataMolGab(f, "[9G]");

        if((posd==-1 && (posf!=-1 || posg!=1)) || (posf==-1 && (posd!=-1 || posg!=1)) || (posg==-1 && (posf!=-1 || posd!=1)))
            _mixte=true;
    }

    else if(_format=="gabedit" && _cart_sphe=="sphe")
        posd = posf = posg = 0;

    else
        posd = posf = posg =-1;

    for(int i=0; i<m; i++)
    {
        if(_shell_types[i]=="s" || _shell_types[i]=="S")
        {
            _L_types.push_back(0);
            _number_of_MO_coefs+=1;
        }

        else if(_shell_types[i]=="p" || _shell_types[i]=="P")
        {
            _L_types.push_back(1);
            _number_of_MO_coefs+=3;
        }

        else if((_shell_types[i]=="d" || _shell_types[i]=="D") && posd==-1)
        {
            _L_types.push_back(2);
            _number_of_MO_coefs+=6;
        }
        
        else if((_shell_types[i]=="d" || _shell_types[i]=="D") && posd!=-1)
        {
            _L_types.push_back(-2);
            _number_of_MO_coefs+=5;
        }
        
        else if((_shell_types[i]=="f" || _shell_types[i]=="F") && posf==-1)
        {
            _L_types.push_back(3);
            _number_of_MO_coefs+=10;
        }

        else if((_shell_types[i]=="f" || _shell_types[i]=="F") && posf!=-1)
        {
            _L_types.push_back(-3);
            _number_of_MO_coefs+=7;
        }
        
        else if((_shell_types[i]=="g" || _shell_types[i]=="G") && posg==-1)
        {
            _L_types.push_back(4);
            _number_of_MO_coefs+=15;
        }

        else if((_shell_types[i]=="g" || _shell_types[i]=="G") && posg!=-1)
        {
            _L_types.push_back(-4);
            _number_of_MO_coefs+=9;
        }

        else if(_format=="molden")
        {
            int N=int(toupper(_shell_types[i][0]))-int('F')+3;
            _L_types.push_back(-N);
            _number_of_MO_coefs+= 2*abs(N)+1;
        }

        else if(_cart_sphe=="sphe")
        {
            int N=int(toupper(_shell_types[i][0]))-int('F')+3;
            _L_types.push_back(-N);
            _number_of_MO_coefs+= 2*abs(N)+1;
        }

        else if(_cart_sphe=="cart")
        {
            int N=int(toupper(_shell_types[i][0]))-int('F')+3;
            _L_types.push_back(N);
            _number_of_MO_coefs+= 2*N+1;
        }

        else
        {
            std::cout<<"Error, shell type no recognize."<<std::endl;
            std::cout<<"Please check your file"<<std::endl;
            exit(1);
        }
    }

    _number_of_MO=_number_of_MO_coefs;
}

void MOLDENGAB::read_one_basis_data(std::istream& f)
{
    std::string p, t;
    int n, m;
    int v = 0;
    double pc, c;

    getline(f, p);
    std::stringstream k(p);
    k >> m;
    getline(f, p);
    
    do{
        _num_center.push_back(m);
        std::stringstream s(p);
        s>>p;
        _shell_types.push_back(p);
        s>>n;
        _number_of_gtf.push_back(n);
        s>>c;

        for(int i=0; i<n; i++)
        {
            v++;
            _factor_coefs.push_back(c);
            getline(f,p);
            if(p.find("D")!=std::string::npos)
                p.replace(p.find("D"),1,"E");
            else if(p.find("d")!=std::string::npos)
                p.replace(p.find("d"),1,"e");

            if(p.find("D")!=std::string::npos)
                p.replace(p.find("D"),1,"E");
            else if(p.find("d")!=std::string::npos)
                p.replace(p.find("d"),1,"e");

            std::stringstream ss(p);
            ss>>pc;
            _exposants.push_back(pc);
            ss>>pc;
            _cgtf_coefs.push_back(pc);
        }

        getline(f,p);
        t=p;
        while(t.find(" ")!=std::string::npos)
            t.erase(t.find(" "),1);

    }while(!t.empty());
    m--;
    _n_at_basis[m]=v;
}

void MOLDENGAB::read_MO_data(std::ifstream& f)
{
    std::string p,t;
    double a;
    std::vector<double> aa;

    if(LocaliseDataMolGab(f,"Spin= Beta")!=-1)
    {
        _alpha_and_beta=false;
    }

    long int pos=LocaliseDataMolGab(f, "[MO]");

    if(pos==-1)
    {
        std::cout<<"Basis (GTO) data not found"<<std::endl;
        std::cout<<"Data required, please check your file"<<std::endl;
        exit(1);
    }

    f.seekg(pos);
    getline(f,p);

    do{
        for(int j=0; j<4; j++)
        {
            if(j!=0)
                getline(f,p);
            std::stringstream s(p);
            s>>p;

            if(p=="Ene=")
            {
                s>>a;
                _MO_energy.push_back(a);
            }
            else if(p=="Spin=")
            {
                s>>p;
                _spin_types.push_back(p);
            }
            else if(p=="Occup=")
            {
                s>>a;
                _occupation.push_back(a);
            }
        }

        for(int k=0; k<_number_of_MO_coefs; k++)
        {
            getline(f,p);
            std::stringstream ss(p);
            ss>>p;
            ss>>a;
            aa.push_back(a);
        }
        _MO_coefs.push_back(aa);
        aa=std::vector<double> ();
        getline(f,p);
        t=p;
        while(t.find(" ")!=std::string::npos)
            t.erase(t.find(" "),1);
    }while(!t.empty());

    _number_of_MO=_MO_energy.size();
}

void MOLDENGAB::PrintData()
{
    std::cout<<"Number of atoms = "<<_number_of_atoms<<std::endl;
    for(size_t i=0; i<_symbol.size(); i++)
        std::cout<<"Symbol "<<i<<" = "<<_symbol[i]<<std::endl;
    for(size_t i=0; i<_atomic_number.size(); i++)
        std::cout<<"Atomic number "<<i<<" = "<<_atomic_number[i]<<std::endl;
    for(size_t i=0; i<_coord.size(); i++)
        for(size_t j=0; j<_coord[i].size(); j++)
            std::cout<<"Coordinates "<<j<<" for atom "<<i<<" = "<<_coord[i][j]<<std::endl;
    for(size_t i=0; i<_num_center.size(); i++)
        std::cout<<"Num center "<<i<<" = "<<_num_center[i]<<std::endl;
    for(size_t i=0; i<_n_at_basis.size(); i++)
        std::cout<<"N at basis "<<i<<" = "<<_n_at_basis[i]<<std::endl;
    for(size_t i=0; i<_shell_types.size(); i++)
        std::cout<<"Shell type "<<i<<" = "<<_shell_types[i]<<std::endl;
    for(size_t i=0; i<_L_types.size(); i++)
        std::cout<<"L type "<<i<<" = "<<_L_types[i]<<std::endl;
    for(size_t i=0; i<_exposants.size(); i++)
        std::cout<<"Exposant "<<i<<" = "<<_exposants[i]<<std::endl;
    for(size_t i=0; i<_number_of_gtf.size(); i++)
        std::cout<<"Number of GTF "<<i<<" = "<<_number_of_gtf[i]<<std::endl;
    for(size_t i=0; i<_cgtf_coefs.size(); i++)
        std::cout<<"GTF coefficient "<<i<<" = "<<_cgtf_coefs[i]<<std::endl;
    for(size_t i=0; i<_factor_coefs.size(); i++)
        std::cout<<"CGTF coefficient "<<i<<" = "<<_factor_coefs[i]<<std::endl;
    for(size_t i=0; i<_alpha_energies.size(); i++)
        std::cout<<"Alpha MO energy "<<i<<" = "<<_MO_energy[i]<<std::endl;
    for(size_t i=0; i<_beta_energies.size(); i++)
        std::cout<<"Beta MO energy "<<i<<" = "<<_MO_energy[i]<<std::endl;
    for(size_t i=0; i<_alpha_MO_coefs.size(); i++)
        for(size_t j=0; j<_alpha_MO_coefs[i].size(); j++)
            std::cout<<"Alpha MO coefficient ["<<i<<"]["<<j<<"] = "<<_alpha_MO_coefs[i][j]<<std::endl;
    for(size_t i=0; i<_beta_MO_coefs.size(); i++)
        for(size_t j=0; j<_beta_MO_coefs[i].size(); j++)
            std::cout<<"Beta MO coefficient ["<<i<<"]["<<j<<"] = "<<_beta_MO_coefs[i][j]<<std::endl;
    for(size_t i=0; i<_alpha_occupation.size(); i++)
        std::cout<<"Alpha Occupation "<<i<<" = "<<_alpha_occupation[i]<<std::endl;
    for(size_t i=0; i<_beta_occupation.size(); i++)
        std::cout<<"Beta Occupation "<<i<<" = "<<_beta_occupation[i]<<std::endl;
    for(size_t i=0; i<_spin_types.size(); i++)
        std::cout<<"Spin types "<<i<<" = "<<_spin_types[i]<<std::endl;
    std::cout<<"Coordinates type = "<<_coord_type<<std::endl;
    std::cout<<"Number of MO = "<<_number_of_MO<<std::endl;
    std::cout<<"Number of MO coefficients = "<<_number_of_MO_coefs<<std::endl;

    if(_alpha_and_beta)
        std::cout<<"Alpha == Beta"<<std::endl;
    else
        std::cout<<"Alpha != Beta"<<std::endl;
}

long int LocaliseDataMolGab(std::ifstream& f, std::string b)
{
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

    return position;
}

long int LocaliseDataMolGabBefore(std::ifstream& f, std::string b)
{
    long int position;
    f.clear();
    f.seekg(0,f.beg);
    std::string test;
    bool ok=false;
    while(!f.eof())
    {    
        position=f.tellg();
        getline(f, test);
        if(test.find(b)!=std::string::npos)
        {
            ok=true;
            break;
        }
    }

    if(!ok) 
        return -1;    

    return position;
}
