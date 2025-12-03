#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Common/Constants.h>
#include <Utils/LOG.h>


LOG::LOG() :
    _number_of_atoms(0),
    _number_of_basis_functions(0),
    _number_of_cartesian_basis_functions(0),
    _number_of_primitive_gaussians(0),
    _number_of_alpha_electrons(0),
    _number_of_beta_electrons(0),
    _energy(0.0),
    _num_center(),
    _symbol(),
    _atomic_numbers(),
    _coordinates(),
    _mulliken_charges(),
    _shell_types(),
    _l_types(),
    _number_of_gtf(),
    _exposants(),
    _cgtf_coefficients(),
    _cgtf_sp_coefficients(),
    _factor_coefficients(),
    _alpha_occupation(),
    _beta_occupation(),
    _alpha_MO_coefs(),
    _beta_MO_coefs(),
    _alpha_energy(),
    _beta_energy(),
    _d_cart_sphe(),
    _f_cart_sphe(),
    _number_of_MO(0),
    _number_of_MO_coefs(0),
    _n_at_basis(),
    _alpha_and_beta(false),
    _mixte(false)
{ }

LOG::LOG(std::ifstream& file) :
    _number_of_atoms(0),
    _number_of_basis_functions(0),
    _number_of_cartesian_basis_functions(0),
    _number_of_primitive_gaussians(0),
    _number_of_alpha_electrons(0),
    _number_of_beta_electrons(0),
    _energy(0.0),
    _num_center(),
    _symbol(),
    _atomic_numbers(),
    _coordinates(),
    _mulliken_charges(),
    _shell_types(),
    _l_types(),
    _number_of_gtf(),
    _exposants(),
    _cgtf_coefficients(),
    _cgtf_sp_coefficients(),
    _factor_coefficients(),
    _alpha_occupation(),
    _beta_occupation(),
    _alpha_MO_coefs(),
    _beta_MO_coefs(),
    _alpha_energy(),
    _beta_energy(),
    _d_cart_sphe(),
    _f_cart_sphe(),
    _number_of_MO(0),
    _number_of_MO_coefs(0),
    _n_at_basis(),
    _alpha_and_beta(false),
    _mixte(false)
{
    read_atoms_data(file);
    read_basis_data(file);
    read_MO_data(file);
}

void LOG::read_atoms_data(std::ifstream& f)
{
    long int pos, pos2;
    int n;
    double d;
    std::string p;

    pos=LocaliseDataLogBefore(f, "NAtoms=");

    f.seekg(pos);
    getline(f,p);
    std::stringstream sp(p);
    sp>>p;
    sp>>_number_of_atoms;

    _coordinates=std::vector<std::vector<double>> (_number_of_atoms);

    pos=LocaliseDataLog(f, "Standard orientation:");

    f.clear();
    f.seekg(0,f.beg);        //Cas où pos=-1;

    if(pos==-1)
        do{
            pos2=pos;
            pos=LocaliseNextDataLog(f,"Input orientation:");
        }while(pos!=-1);

    else
        do{
            pos2=pos;
            pos=LocaliseNextDataLog(f,"Standard orientation:");
        }while(pos!=-1);
    
    f.clear();
    f.seekg(pos2);
    getline(f,p);
    getline(f,p);
    getline(f,p);
    getline(f,p);

    for(int i=0; i<_number_of_atoms; i++)
    {
        getline(f,p);
        std::stringstream s(p);
        s>>p;
        s>>n;
        _atomic_numbers.push_back(n);
        s>>p;

        for(int j=0; j<3; j++)
        {
            s>>d;
            d*=ANGTOBOHR;
            _coordinates[i].push_back(d);
        }
    }
}

void LOG::read_basis_data(std::ifstream& f)
{
    long int position, lastPosition;
    int n, l;
    double d, c;
    std::string line;

    position = LocaliseDataLog(f, "AO basis set in the form of general basis input (Overlap normalization):");
    
    do
    {
        lastPosition = position;
        position = LocaliseNextDataLog(f, "AO basis set in the form of general basis input (Overlap normalization):");
    } while (position != -1);

    f.clear();
    f.seekg(lastPosition);

    _n_at_basis = std::vector<int>(_number_of_atoms);
    int m = 0;
    int v;

    for(int i = 0; i < _number_of_atoms; ++i)
    {
        getline(f, line);

        std::stringstream nc(line);
        nc >> l;
        
        getline(f, line);

        v = 0;
        do
        {
            std::stringstream ss(line);
            ss >> line;
            _shell_types.push_back(line);
            
            ss >> n;
            _number_of_gtf.push_back(n);
            
            ss >> c;

            for(int j = 0; j < n; ++j)
            {
                v++;

                _num_center.push_back(l);
                _factor_coefficients.push_back(c);

                getline(f, line);


                // TODO (lgardre) : remplacer par des regex
                if (line.find("D") != std::string::npos)
                {
                    line.replace(line.find("D"), 1, "E");
                }
                else if (line.find("d") != std::string::npos)
                {
                    line.replace(line.find("d"), 1, "e");
                }

                if (line.find("D") != std::string::npos)
                {
                    line.replace(line.find("D"), 1, "E");
                }
                else if (line.find("d") != std::string::npos)
                {
                    line.replace(line.find("d"), 1, "e");
                }

                if(_shell_types[m] == "SP" || _shell_types[m] == "sp")
                {
                    if (line.find("D") != std::string::npos)
                    {
                        line.replace(line.find("D"), 1, "E");
                    }
                    else if (line.find("d") != std::string::npos)
                    {
                        line.replace(line.find("d"), 1, "e");
                    }
                }

                std::stringstream sss(line);
                sss >> d;
                _exposants.push_back(d);

                sss >> d;
                _cgtf_coefficients.push_back(d);

                if(_shell_types[m] == "SP" || _shell_types[m] == "sp")                    //A voir comment on organise les données !!!
                {
                    sss >> d;
                    _cgtf_sp_coefficients.push_back(d);
                }
                else
                {
                    _cgtf_sp_coefficients.push_back(0.0);
                }
            }

            getline(f, line);

            m++;
        } while (line.find("*") == std::string::npos);

        _n_at_basis[i] = v;
    }

    position = LocaliseNextDataLogBefore(f, "basis functions,");

    f.clear();
    f.seekg(position);

    getline(f, line);

    std::stringstream ss(line);
    ss >> _number_of_basis_functions;
    ss >> line;
    ss >> line;
    ss >> _number_of_primitive_gaussians;
    ss >> line;
    ss >> line;
    ss >> _number_of_cartesian_basis_functions;

    getline(f, line);

    std::stringstream sss(line);
    sss >> _number_of_alpha_electrons;
    sss >> line;
    sss >> line;
    sss >> _number_of_beta_electrons;

    do
    {
        lastPosition = position;
        position = LocaliseNextDataLogBefore(f, "E(");
    } while (position != -1);

    f.clear();
    f.seekg(lastPosition);

    getline(f,line);
    
    std::stringstream ssss(line);
    ssss >> line;
    ssss >> line;
    ssss >> line;
    ssss >> line;
    ssss >> _energy;

    position = LocaliseDataLogBefore(f, "Standard basis:", "General basis read from cards:");

    if(position == -1)
    {
        _d_cart_sphe = "sphe";
        _f_cart_sphe = "sphe";
    }
    else
    {
        do
        {
            lastPosition = position;
            position = LocaliseNextDataLogBefore(f, "Standard basis:", "General basis read from cards:");
        } while(position != -1);
    }

    f.clear();
    f.seekg(lastPosition);

    getline(f, line);
    
    if(line.find("5D") != std::string::npos || line.find("5d") != std::string::npos)
    {
        _d_cart_sphe = "sphe";
    }
    else if(line.find("6D") != std::string::npos || line.find("6d") != std::string::npos)
    {
        _d_cart_sphe = "cart";
    }
    else
    {
        _d_cart_sphe = "sphe";
    }

    if(line.find("7F") != std::string::npos || line.find("7f") != std::string::npos)
    {
        _f_cart_sphe = "sphe";
    }
    else if(line.find("10F") != std::string::npos || line.find("10f") != std::string::npos)
    {
        _f_cart_sphe = "cart";
    }
    else
    {
        _f_cart_sphe = "sphe";
    }

    for(size_t i = 0; i < _shell_types.size(); ++i)
    {
        if(_shell_types[i] == "s" || _shell_types[i] == "S")
        {
            _l_types.push_back(0);
            _number_of_MO += 1;
        }
        else if(_shell_types[i] == "p" || _shell_types[i] == "P")
        {
            _l_types.push_back(1);
            _number_of_MO += 3;
        }
        else if(_shell_types[i] == "sp" || _shell_types[i] == "SP")
        {
            _l_types.push_back(-1);
            _number_of_MO += 4;
        }
        else if((_shell_types[i] == "d" || _shell_types[i] == "D") && _d_cart_sphe == "cart")
        {
            _l_types.push_back(2);
            _number_of_MO += 6;
        }
        else if(_shell_types[i] == "d" || _shell_types[i] == "D")
        {
            _l_types.push_back(-2);
            _number_of_MO += 5;
        }
        else if((_shell_types[i] == "f" || _shell_types[i] == "F") && _f_cart_sphe == "cart")
        {
            _l_types.push_back(3);
            _number_of_MO += 10;
        }
        else if(_shell_types[i] == "f" || _shell_types[i] == "F")
        {
            _l_types.push_back(-3);
            _number_of_MO += 7;
        }
        else
        {
            int N = int(toupper(_shell_types[i][0])) - int('F') + 3;
            _l_types.push_back(-N);
            _number_of_MO += 2 * N + 1;

            if(_d_cart_sphe == "cart" || _f_cart_sphe == "cart")
            {
                _mixte = true;
            }
        }
    }
    
    _number_of_MO_coefs = _number_of_MO;
    _beta_MO_coefs = _alpha_MO_coefs = std::vector<std::vector<double>> (_number_of_MO_coefs);

    if(_d_cart_sphe != _f_cart_sphe)
    {
        _mixte = true;
    }
}

void LOG::read_MO_data(std::ifstream& f)
{
    long int pos, pos2;
    int n,m;
    double d;
    std::string p,p2,pp, name;

    pos=LocaliseDataLog(f, "Alpha Molecular Orbital Coefficients:");

    if(pos==-1)
    {
        _alpha_and_beta=true;
        pos=LocaliseDataLog(f, "Molecular Orbital Coefficients:");

        do{
            pos2=pos;
            pos=LocaliseNextDataLog(f, "Molecular Orbital Coefficients:");
        }while(pos!=-1);
    }

    else
    {
        _alpha_and_beta=false;
        do{
            pos2=pos;
            pos=LocaliseNextDataLog(f, "Alpha Molecular Orbital Coefficients:");
        }while(pos!=-1);
    }

    f.clear();
    f.seekg(pos2);

    getline(f,p);

    n=0;

    do{
        m=0;
        std::stringstream t(p);
        do{
            t>>p;
            m++;
        }while(!t.eof());

        getline(f,p);
        std::stringstream s(p);
        
        for(int i=0; i<m; i++)
        {
            s>>p;

            if(p.find("O")!=std::string::npos && _alpha_and_beta)
                _alpha_occupation.push_back(2.0);
            else if(p.find("O")!=std::string::npos)
                _alpha_occupation.push_back(1.0);
            else if(p.find("V")!=std::string::npos)
                _alpha_occupation.push_back(0.0);
        }

        getline(f,p);
        std::stringstream ss(p);
        ss>>p;
        ss>>p;

        for(int i=0; i<m; i++)
        {
            ss>>d;
            _alpha_energy.push_back(d);
        }

        for(int i=0; i<_number_of_basis_functions; i++)
        {
            getline(f,p);
            std::stringstream sss(p);

            sss>>p2;

            if(p.find("1S")!=std::string::npos || p.find("1s")!=std::string::npos)
            {
                sss>>p;
                sss>>p;
                if(n==0)
                    _symbol.push_back(p);
            }
            
            sss>>p;
            
            if(_d_cart_sphe=="sphe" and (p.find("d")!=std::string::npos or p.find("D")!=std::string::npos) and (p.size()==2 or p.size()==3))
                sss>>p;    
            
            if(_f_cart_sphe=="sphe" and (p.find("f")!=std::string::npos or p.find("F")!=std::string::npos) and (p.size()==2 or p.size()==3))
                sss>>p;
            
            for(int j=0; j<m; j++)
            {
                sss>>d;
                _alpha_MO_coefs[n+j].push_back(d);
            }
        }
        n+=m;
        getline(f,p);
        std::stringstream nn;
        nn<<n+1;
        pp=nn.str();
    }while(p.find(pp)!=std::string::npos);

    if(!_alpha_and_beta)
    {
        f.clear();
        f.seekg(pos2);
        pos=LocaliseNextDataLog(f, "Beta Molecular Orbital Coefficients:");

        f.clear();
        f.seekg(pos);

        getline(f,p);

        n=0;

        do{
            m=0;
            std::stringstream t(p);
            do{
                t>>p;
                m++;
            }while(!t.eof());

            getline(f,p);
            std::stringstream s(p);
        
            for(int i=0; i<m; i++)
            {
                s>>p;

                if(p.find("O")!=std::string::npos)
                    _beta_occupation.push_back(1.0);
                else if(p.find("V")!=std::string::npos)
                    _beta_occupation.push_back(0.0);
            }

            getline(f,p);
            std::stringstream ss(p);
            ss>>p;
            ss>>p;

            for(int i=0; i<m; i++)
            {
                ss>>d;
                _beta_energy.push_back(d);
            }

            for(int i=0; i<_number_of_MO; i++)
            {
                getline(f,p);
                std::stringstream sss(p);

                if(p.find("1S")!=std::string::npos || p.find("1s")!=std::string::npos)
                {
                    sss>>p;
                    sss>>p;
                }

                sss>>p;
                sss>>p;
        
                for(int j=0; j<m; j++)
                {
                    sss>>d;
                    _beta_MO_coefs[n+j].push_back(d);
                }
            }
            n+=m;
            getline(f,p);
            std::stringstream nn;
            nn<<n+1;
            pp=nn.str();
        }while(p.find(pp)!=std::string::npos);
    }

    else
    {
        _beta_occupation=_alpha_occupation;
        _beta_MO_coefs=_alpha_MO_coefs;
        _beta_energy=_alpha_energy;
    }

    pos=LocaliseNextDataLog(f, "Mulliken charges:");

    f.clear();
    f.seekg(pos);

    getline(f,p);

    for(int i=0; i<_number_of_atoms; i++)
    {
        getline(f,p);
        std::stringstream mc(p);
        mc>>p;
        mc>>p;
        mc>>d;
        _mulliken_charges.push_back(d);
    }
}

void LOG::PrintData()
{
    std::cout << "Number of atoms = " << _number_of_atoms << std::endl;
    std::cout << "Number of basis functions = " << _number_of_basis_functions << std::endl;
    std::cout << "Number of cartesian basis functions = " << _number_of_cartesian_basis_functions << std::endl;
    std::cout << "Number of primitive gaussians = " << _number_of_primitive_gaussians << std::endl;
    std::cout << "Number of alpha electrons = " << _number_of_alpha_electrons << std::endl;
    std::cout << "Number of beta electrons = " << _number_of_beta_electrons << std::endl;
    std::cout << "Energy = " << _energy << std::endl;

    for(size_t i = 0; i < _num_center.size(); ++i)
    {
        std::cout << "Number center [" << i << "] = " << _num_center[i] << std::endl;
    }

    for(size_t i = 0; i < _symbol.size(); ++i)
    {
        std::cout << "Symbol [" << i << "] = " << _symbol[i] << std::endl;
    }

    for(size_t i = 0; i < _atomic_numbers.size(); ++i)
    {
        std::cout << "Atomic number [" << i << "] = " << _atomic_numbers[i] << std::endl;
    }
    
    for(size_t i = 0; i < _coordinates.size(); ++i)
    {
        std::cout << "Coordinates atom [" << i << "] = ";
        for(size_t j = 0; j < _coordinates[i].size(); ++j)
        {
            std::cout << _coordinates[i][j] << "    ";
        }
        std::cout << std::endl;
    }

    for(size_t i = 0; i < _mulliken_charges.size(); ++i)
    {
        std::cout << "Mulliken charges " << i << " = " << _mulliken_charges[i] << std::endl;
    }

    for(size_t i = 0; i < _shell_types.size(); ++i)
    {
        std::cout << "Shell types " << i << " = " << _shell_types[i] << std::endl;
    }

    for(size_t i = 0; i < _l_types.size(); ++i)
    {
        std::cout << "Ltypes [" << i << "] = " << _l_types[i] << std::endl;
    }

    for(size_t i = 0; i < _number_of_gtf.size(); ++i)
    {
        std::cout << "Number of GTF in CGTF " << i << " = " << _number_of_gtf[i] << std::endl;
    }

    for(size_t i = 0; i < _exposants.size(); ++i)
    {
        std::cout << "Exposant " << i << " = " << _exposants[i] << std::endl;
    }

    for(size_t i = 0; i < _cgtf_coefficients.size(); ++i)
    {
        std::cout << "CGTF coefficient " << i << " = " << _cgtf_coefficients[i] << std::endl;
    }

    for(size_t i = 0; i < _cgtf_sp_coefficients.size(); ++i)
    {
        std::cout << "CGTF SP coefficient " << i << " = " << _cgtf_sp_coefficients[i] << std::endl;
    }

    for(size_t i = 0; i < _factor_coefficients.size(); ++i)
    {
        std::cout << "Factor coefficient " << i << " = " << _factor_coefficients[i] << std::endl;
    }

    for(size_t i = 0; i < _alpha_occupation.size(); ++i)
    {
        std::cout << "Alpha occupation " << i << " = " << _alpha_occupation[i] << std::endl;
    }

    for(size_t i = 0; i < _beta_occupation.size(); ++i)
    {
        std::cout << "Beta occupation " << i << " = " << _beta_occupation[i] << std::endl;
    }

    for(size_t i = 0; i < _alpha_MO_coefs.size(); ++i)
    {
        for(size_t j = 0; j < _alpha_MO_coefs[i].size(); ++j)
        {
            std::cout << "Alpha MO coefficients [" << i << "][" << j << "] = " << _alpha_MO_coefs[i][j] << std::endl;
        }
    }

    for(size_t i = 0; i < _beta_MO_coefs.size(); ++i)
    {
        for(size_t j = 0; j < _beta_MO_coefs[i].size(); ++j)
        {
            std::cout << "Beta MO coefficients [" << i << "][" << j << "] = " << _beta_MO_coefs[i][j] << std::endl;
        }
    }
    
    for(size_t i = 0; i < _alpha_energy.size(); ++i)
    {
        std::cout << "Alpha energy " << i << " = " << _alpha_energy[i] << std::endl;
    }

    for(size_t i = 0; i < _beta_energy.size(); ++i)
    {
        std::cout << "Beta energy " << i << " = " << _beta_energy[i] << std::endl;
    }

    std::cout << "D cart/sphe = " << _d_cart_sphe << std::endl;
    std::cout << "F cart/sphe = " << _f_cart_sphe << std::endl;

    std::cout << "Number of MO = " << _number_of_MO << std::endl;
    std::cout << "Number of MO coefficients = " << _number_of_MO_coefs << std::endl;

    for(size_t i = 0; i < _n_at_basis.size(); ++i)
    {
        std::cout << "n in basis " << i << " = " << _n_at_basis[i] << std::endl;
    }

    if(_alpha_and_beta)
    {
        std::cout << "Alpha == Beta" << std::endl;
    }
    else
    {
        std::cout << "Alpha != Beta" << std::endl;
    }
}

long int LocaliseDataLog(std::ifstream& f, std::string b)
{
    bool ok = false;
    
    f.clear();
    f.seekg(0, f.beg);

    std::string line;
    long int position;
    while(!f.eof())
    {    
        getline(f, line);
        if(line.find(b) != std::string::npos)
        {
            position = f.tellg();
            ok = true;
            break;
        }
    }

    return ok ? position : -1;
}

long int LocaliseDataLogBefore(std::ifstream& f, std::string b1, std::string b2)
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
        if(test.find(b1)!=std::string::npos || test.find(b2)!=std::string::npos)
        {    
            ok=true;
            break;
        }
    }

    if(!ok) 
        return -1;    

    return position;
}

long int LocaliseDataLogBefore(std::ifstream& f, std::string b)
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

long int LocaliseNextDataLog(std::ifstream& f, std::string b)
{
    bool ok = false;

    f.clear();

    std::string line;
    long int position;
    while(!f.eof())
    {    
        getline(f, line);
        if(line.find(b) != std::string::npos)
        {
            position = f.tellg();
            ok = true;
            break;
        }
    }

    return ok ? position : -1;
}

long int LocaliseNextDataLogBefore(std::ifstream& f, std::string b)
{
    long int position;
    f.clear();
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

long int LocaliseNextDataLogBefore(std::ifstream& f, std::string b1, std::string b2)
{
    long int position;
    f.clear();
    std::string test;
    bool ok=false;
    while(!f.eof())
    {
        position=f.tellg();    
        getline(f, test);
        if(test.find(b1)!=std::string::npos || test.find(b2)!=std::string::npos)
        {
            ok=true;
            break;
        }
    }

    if(!ok) 
        return -1;    

    return position;
}
