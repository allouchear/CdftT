using namespace std;

#include <cmath>
#include <iomanip>

#include <Common/Constants.h>
#include <Common/PeriodicTable.h>
#include <Common/Structure.h>
#include <Cube/GridCP.h>
#include <Utils/Enums.hpp>

#ifdef ENABLE_OMP
#include <omp.h>
#endif

/* Point with a volume number :
 * i>0 : point of volume # i
 * i<0 : the critical point of volume # -i
 * i= : point not yet assigned
 */
#define TOL 1e-12

/* See W. Tang et al J. Phys. Condens. Matter 21 (2009) 084204 */


/**************************************************************************/
void GridCP::resetKnown()
{
    for(size_t i=0;i< _known.size() ;i++)
    for(size_t j=0;j< _known[i].size() ;j++)
    for(size_t k=0;k< _known[i][j].size();k++)
        _known[i][j][k] = 0;
}
/**************************************************************************/
CriticalPoint  GridCP::newCriticalPoint(int i, int j, int k, int numV)
{
    CriticalPoint cp;
    cp.index[0] = i;
    cp.index[1] = j;
    cp.index[2] = k;
    cp.rank = 0;
    cp.signature = 0;
    cp.lambda[0] = 0;
    cp.lambda[1] = 0;
    cp.lambda[2] = 0;
    cp.integral = 0;
    cp.volume = 0;
    cp.nuclearCharge = 0;
    cp.numVolume = numV;
    cp.numCenter = 0;
    cp.value = 0;
    return cp;
}
/**************************************************************************/
void GridCP::reset()
{
    _volumeNumberOfPoints = vector< vector< vector<int>> > ();
    _known = vector< vector<vector<int>>  > ();
    _integral = 0;
    _nuclearCharge = 0;
    _criticalPoints=vector< CriticalPoint > ();
    _domain = Domain();
    _str = Structure();
}
/**************************************************************************/
GridCP::GridCP()
{
    reset();
}
/**************************************************************************/
vector< vector < vector<int>> > GridCP::get3DIntVector()
{
    vector<int> V(_domain.get_N3(),0);
    vector<vector<int>> M(_domain.get_N2(), V);
    vector<vector<vector<int>>> V3D(_domain.get_N1(), M);
    return V3D;
}
/**************************************************************************/
bool GridCP::okDomain()
{
    bool ok=false;
    if (_domain.get_N1()<2) return ok;
    if (_domain.get_N2()<2) return ok;
    if (_domain.get_N3()<2) return ok;
    if (_domain.get_Nval()<1) return ok;
    return true;
}
/**************************************************************************/
void GridCP::initGridCP(const Grid& grid, bool ongrid)
{
    reset();
    _domain = grid.get_domain();
    if (!okDomain()) return;
    _str = grid.get_structure();
    _V = grid.get_values();
    int nval = _V[0][0][0].size();
    if (nval<4 && !ongrid)
    {
        cout<<"begin grad"<<endl;
        Grid g= grid.gradient(1);
        cout<<"end grad"<<endl;
        _V = g.get_values();
        _domain = g.get_domain();
        _str = g.get_structure();
    }
    if (!ongrid)
    for(size_t i=0;i<_V.size();i++)
    for(size_t j=0;j<_V[i].size();j++)
    for(size_t k=0;k<_V[i][j].size();k++)
    {
        double c= abs(_V[i][j][k][1]);
        if (c < abs(_V[i][j][k][2])) c= abs(_V[i][j][k][2]);
        if (c < abs(_V[i][j][k][3])) c= abs(_V[i][j][k][3]);
        if (c>0)
        {
            c = 1.0/c;
            for(size_t l = 1;l<=3;l++)
                _V[i][j][k][l] *= c;
        }
    }
    if (!ongrid)
    for(size_t i=1;i<_V.size()-1;i++)
    for(size_t j=1;j<_V[i].size()-1;j++)
    for(size_t k=1;k<_V[i][j].size()-1;k++)
    {
        if ( _V[i-1][j][k][0]<_V[i][j][k][0] && _V[i+1][j][k][0]<_V[i][j][k][0]) _V[i][j][k][1]  = 0;
        if ( _V[i][j-1][k][0]<_V[i][j][k][0] && _V[i][j+1][k][0]<_V[i][j][k][0]) _V[i][j][k][2]  = 0;
        if ( _V[i][j][k-1][0]<_V[i][j][k][0] && _V[i][j][k+1][0]<_V[i][j][k][0]) _V[i][j][k][3]  = 0;
    }

    _volumeNumberOfPoints = get3DIntVector();
    _known = get3DIntVector();
}
/**************************************************************************/
bool GridCP::isVolumeEdge(int current[])
{
    if (!okDomain()) return 0;
    int i = current[0];
    int j = current[1];
    int k = current[2];
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };

    for(size_t ic=0;ic<3;ic++)
    for(size_t jc=0;jc<3;jc++)
    for(size_t kc=0;kc<3;kc++)
        if (2!=_known[I[ic]][J[jc]][K[kc]] && abs(_volumeNumberOfPoints[i][j][k]) != abs(_volumeNumberOfPoints[I[ic]][J[jc]][K[kc]]))
            return true;

    return false;
}
/**************************************************************************/
int GridCP::setArroundTo(int current[], int kn)
{
    if (!okDomain()) return 0;

    int i = current[0];
    int j = current[1];
    int k = current[2];
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };
    int n = 0;


    for(size_t ic=0;ic<3;ic++)
    for(size_t jc=0;jc<3;jc++)
    for(size_t kc=0;kc<3;kc++)
    {
        if (ic==1 && jc==1 && kc ==1) continue;
        if (_known[I[ic]][J[jc]][K[kc]] != 1 && _volumeNumberOfPoints[I[ic]][J[jc]][K[kc]]==_volumeNumberOfPoints[i][j][k]) 
        {
            _known[I[ic]][J[jc]][K[kc]] = kn;
            n++;
        }
    }
    return n;
}
/**************************************************************************/
bool GridCP::isMax(int current[])
{
    if (!okDomain()) return false;
    int i = current[0];
    int j = current[1];
    int k = current[2];
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };


    for(size_t ic=0;ic<3;ic++)
    for(size_t jc=0;jc<3;jc++)
    for(size_t kc=0;kc<3;kc++)
        if (_V[I[ic]][J[jc]][K[kc]][0]>_V[I[1]][J[1]][K[1]][0]) return false;

    return true;
}
/**************************************************************************/
bool GridCP::nextPointOnGrid(int current[], int next[])
{
    if (!okDomain()) return false;
    int i = current[0];
    int j = current[1];
    int k = current[2];
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };

    double dx;
    double dy;
    double dz;
    double rho;
    double drho;
    double dr;
    double x0,y0,z0;


    int im = 1;
    int jm = 1;
    int km = 1;
    double rhoCenter = _V[i][j][k][0];
    double drhoMax = 0;
    x0=_domain.x(i,j,k);
    y0=_domain.y(i,j,k);
    z0=_domain.z(i,j,k);
    for(size_t ic=0;ic<3;ic++)
    for(size_t jc=0;jc<3;jc++)
    for(size_t kc=0;kc<3;kc++)
    {
        if (ic==1 && jc==1 && kc==1) continue;
        if (_known[I[ic]][J[jc]][K[kc]] >1) continue;
        rho =_V[I[ic]][J[jc]][K[kc]][0];
        /*
        dx=(I[ic]-I[1])*_domain.dx();
        dy=(J[jc]-J[1])*_domain.dy();
        dz=(K[kc]-K[1])*_domain.dz();
        */
        dx=_domain.x(I[ic],J[jc],K[kc])-x0;
        dy=_domain.y(I[ic],J[jc],K[kc])-y0;
        dz=_domain.z(I[ic],J[jc],K[kc])-z0;
        dr = sqrt(dx*dx+dy*dy+dz*dz); 
        drho = (rho-rhoCenter)/dr;
        if (drho>drhoMax)
        {
            drhoMax = drho;
            im = ic;
            jm = jc;
            km = kc;
        }
    }
    next[0] = I[im];
    next[1] = J[jm];
    next[2] = K[km];
    return true;
}
/**************************************************************************/
bool GridCP::nextPoint(double deltaR[], int current[], int next[])
{
    if (!okDomain()) return false;
    int nval = _V[0][0][0].size();
    if (nval<4) return false;
    double gradrl[3];
    // Initialize indices for neighboring points
    int i = current[0];
    int j = current[1];
    int k = current[2];
    int N[3] ={ _domain.get_N1(), _domain.get_N2(), _domain.get_N3()};

    for(size_t c=0;c<3;c++) gradrl[c] = _V[i][j][k][c+1];
    if (int(rint(gradrl[0])) ==0 && int(rint(gradrl[1])) ==0 && int(rint(gradrl[2])) ==0)
    {
        if (isMax(current))
        {
            for(size_t c=0;c<3;c++) next[c] = current[c];
            for(size_t c=0;c<3;c++) deltaR[c] = 0.0;
            return true;
        }
        else
        {
            nextPointOnGrid(current, next);
            for(size_t c=0;c<3;c++) deltaR[c] = 0.0;
        }
    }
    else
    {
        for(size_t c=0;c<3;c++) next[c] = current[c] + int(rint(gradrl[c])); 
        for(size_t c=0;c<3;c++) deltaR[c] += gradrl[c]-int(rint(gradrl[c]));
        for(size_t c=0;c<3;c++) next[c] += int(rint(deltaR[c]));
        for(size_t c=0;c<3;c++) deltaR[c] -=  int(rint(deltaR[c]));
        for(size_t c=0;c<3;c++) if (next[c]<0 ) next[c] = 0;
        for(size_t c=0;c<3;c++) if (next[c]>N[c]-1) next[c] = N[c]-1;

        i = current[0];
        j = current[1];
        k = current[2];
        _known[i][j][k] = 1;
        i = next[0];
        j = next[1];
        k = next[2];
        if (_known[i][j][k]==1)
        {
            nextPointOnGrid(current, next);
            for(size_t c=0;c<3;c++) deltaR[c] = 0.0;
        }

    }
    return true;
}
/**************************************************************************/
bool GridCP::addSurroundingEqualPoints(int current[], vector<vector<int>>& listOfVisitedPoints)
{
    int i = current[0];
    int j = current[1];
    int k = current[2];
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };

    double rho0 = 0;
    double dRho = 0;

    if (!okDomain()) return false;

    rho0 = _V[I[1]][J[1]][K[1]][0];
    for(size_t ic=0;ic<3;ic++)
    for(size_t jc=0;jc<3;jc++)
    for(size_t kc=0;kc<3;kc++)
    {
        if (ic==1 && jc==1 && kc==1) continue;
        dRho =_V[I[ic]][J[jc]][K[kc]][0]-rho0;
        if (fabs(dRho)<TOL)
        {
            vector<int>  data = { I[ic], J[jc], K[kc]};
            listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);
            _known[I[ic]][J[jc]][K[kc]] = 1;
        }
    }
    return true;
}
/**************************************************************************/
vector<vector<int>> GridCP::assentTrajectory(int current[], bool ongrid, bool refine)
{
    vector<vector<int>> listOfVisitedPoints;
    if (!okDomain()) return listOfVisitedPoints;

    vector<int>  data = { current[0], current[1], current[2]};
    listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);

    size_t imax = _domain.get_N1()* _domain.get_N2()* _domain.get_N3();
    int next[3];
    double deltaR[3] = {0,0,0};
    for(size_t l=0;l<imax;l++)
    {

        if (ongrid)
            nextPointOnGrid(current, next);
        else
            nextPoint(deltaR, current, next);
        vector<int>  data = { next[0], next[1], next[2]};
        listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);

        if ((next[0] == current[0] && next[1] == current[1] && next[2] == current[2])
        )
        {
            addSurroundingEqualPoints(current, listOfVisitedPoints);
            break;
        }
        if (!refine && _volumeNumberOfPoints [next[0]][next[1]][next[2]] != 0)
        {
            for(size_t c=0;c<3;c++)current[c] = next[c];
            break;
        }
        for(size_t c=0;c<3;c++)current[c] = next[c];
    }
    return listOfVisitedPoints;
}
/**************************************************************************/
void GridCP::computeVolumes()
{
    for(int i=0;i<_domain.get_N1();i++)
    {
        for(int j=0;j<_domain.get_N2();j++)
        for(int k=0;k<_domain.get_N3();k++)
        {
            int n=abs(_volumeNumberOfPoints [i][j][k]);
            for(size_t ip=0;ip<_criticalPoints.size();ip++)
                if (_criticalPoints[ip].numVolume==n) 
                    _criticalPoints[ip].volume +=1 ;
        }
    }
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
            _criticalPoints[ip].volume *= _domain.get_dv() ;
}
/**************************************************************************/
void GridCP::computeNumCenters()
{
    double r[3];
    double dr[3];
    double r2;
    double r2old = 0;

    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        CriticalPoint  cp = _criticalPoints[ip];
        int i = cp.index[0];
        int j = cp.index[1];
        int k = cp.index[2];
        _criticalPoints[ip].numCenter=0;
        for(int c=0; c<3;c++)
            r[c] = _domain.get_origin()[ c ] + i*_domain.get_T()[c][0] + j*_domain.get_T()[c][1] +  k*_domain.get_T()[c][2]; 
        for(int ia=0;ia<_str.number_of_atoms();ia++)
        {
            for(int c=0; c<3;c++)
                dr[c] = r[c]-_str.atom(ia).get_coordinates()[c];
            r2 = 0;
            for(int c=0; c<3;c++)
                r2 += dr[c]*dr[c];
            
            if (ia==0 || r2<r2old )
            {
                _criticalPoints[ip].numCenter=ia;
                _criticalPoints[ip].nuclearCharge=_str.atom(ia).get_atomicNumber();
                r2old = r2;
            }
        }
    }
}
/**************************************************************************/
void GridCP::computeVolCenters()
{
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        _criticalPoints[ip].integral=0;
        int N=0;
        for(int i=0;i<_domain.get_N1();i++)
        for(int j=0;j<_domain.get_N2();j++)
        for(int k=0;k<_domain.get_N3();k++)
        for(int I=0;I<3;I++)
        {
            int n=abs(_volumeNumberOfPoints [i][j][k]);
            if (n==_criticalPoints[ip].numVolume)
            {
                _criticalPoints[ip].volCenter[0]+=i;
                _criticalPoints[ip].volCenter[1]+=j;
                _criticalPoints[ip].volCenter[2]+=k;
                N++;
            }
        }
        for(int i=0;i<3;i++)
            _criticalPoints[ip].volCenter[i]/=N;
    }
}
/**************************************************************************/
void GridCP::removeBasins0()
{
    vector< CriticalPoint > newCriticalPoints;
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        CriticalPoint  cp = _criticalPoints[ip];
        int i = cp.index[0];
        int j = cp.index[1];
        int k = cp.index[2];
        if (_V[i][j][k][0]<TOL)
            _volumeNumberOfPoints [i][j][k] = 0;
        else
            newCriticalPoints.push_back(cp);
    }
    _criticalPoints= newCriticalPoints;
}
/**************************************************************************/
void GridCP::removeNonSignificantBasins()
{
    int nAtoms=_str.number_of_atoms();
    int nCP=_criticalPoints.size();

    if (nCP==nAtoms) return;
    if (nCP<nAtoms)
    {
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"WARNING : Number of attractors < number of atoms    "<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        return;
    }
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    cout<<"WARNING : Number of attractors > number of atoms    "<<endl;
    cout<<"        : We will remove "<<nCP-nAtoms<<" attractor(s)"<<endl;
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

    vector< double > nPoints(_criticalPoints.size(),0);
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {

        CriticalPoint cp = _criticalPoints[ip];
        nPoints[ip] = cp.volume/_domain.get_dv();
    }
    // sort by nPoints
    for(size_t i=0;i<_criticalPoints.size();i++)
    {
        size_t k=0;
        for(size_t j=i+1;j<_criticalPoints.size();j++)
        {
            if (nPoints[j]>nPoints[k])
                k=j;
        }
        if (k!=i)
        {
            double t=nPoints[i];
            CriticalPoint cp = _criticalPoints[i];
            nPoints[i]=nPoints[k];
            nPoints[k]=t;
            _criticalPoints[i] = _criticalPoints[k];
            _criticalPoints[k] = cp;
        }
    }

    _criticalPoints.resize(nAtoms);
}
/**************************************************************************/
int GridCP::refineEdgeNearGrad()
{
    if (!okDomain()) return 0;
    size_t N[3] ={ size_t(_domain.get_N1()), size_t(_domain.get_N2()), size_t(_domain.get_N3())};

    string str ="Refining grid points adjacent to Bader surface... Please wait";
    int current[3];
    int ne = 0;

    resetKnown();

    for(size_t i=0;i<N[0];i++)
    for(size_t j=0;j<N[1];j++)
    for(size_t k=0;k<N[2];k++)
    {
        if (_V[i][j][k][0]<TOL) _known[i][j][k] = 2;
        if (_volumeNumberOfPoints [i][j][k]<0) _known[i][j][k] = 2;
    }
    for(size_t i=0;i<N[0];i++)
    for(size_t j=0;j<N[1];j++)
    for(size_t k=0;k<N[2];k++)
    {
        current[0] = i;
        current[1] = j;
        current[2] = k;
        if (_volumeNumberOfPoints[i][j][k]>0 && !isMax(current) && isVolumeEdge(current)) 
        {
            ne++;
            setArroundTo(current, 0);
            _known[i][j][k] = 1;
        }
        else _known[i][j][k] = 2; 
    }

    ne = 0;
    for(size_t i=0;i<N[0];i++)
    for(size_t j=0;j<N[1];j++)
    for(size_t k=0;k<N[2];k++)
    {
        if (_known[i][j][k]==1) 
        {
            ne++;
            _volumeNumberOfPoints [i][j][k] = 0;
            _known[i][j][k]=0;
        }
    }

    //double scal = 1.0/(N[0]-1);
    cout<<str<<endl;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
    for(size_t i=0;i<N[0];i++)
    {

        //cout<<setw(10)<<fixed<<setprecision(5)<<i*scal*100<<"%";
        //cout<<" ; Number of critical points = "<<_criticalPoints.size()<<endl;
        for(size_t j=0;j<N[1];j++)
        {
            for(size_t k=0;k<N[2];k++)
            {
                int current[3];
                vector<vector<int>> listOfVisitedPoints;
                current[0] = i;
                current[1] = j;
                current[2] = k;
                if (_volumeNumberOfPoints [i][j][k]!=0) continue;
                if (_known[i][j][k] !=0 ) continue;

                listOfVisitedPoints = assentTrajectory(current, false, true);
                _volumeNumberOfPoints [i][j][k] = abs(_volumeNumberOfPoints[current[0]][current[1]][current[2]]);
#ifdef ENABLE_OMP
#pragma omp critical
#endif
                for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _known[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = 0;
            }
        }
    }
    resetKnown();
    return ne;
}
/**************************************************************************/
void GridCP::assignPointsByGradient(bool ongrid)
{
    double scal;
    string str ="Assigning points to volumes... Please wait";

    if (!okDomain()) return;

    if (_criticalPoints.size()>0)
        _criticalPoints.clear();

    resetKnown();
    int numberOfCriticalPoints = 0;

    scal = 1.0/(_domain.get_N1()-1);
    cout<<str<<endl;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
    for(int i=0;i<_domain.get_N1();i++)
    {
        int current[3];
#ifdef ENABLE_OMP
#ifndef G_OS_WIN32
#pragma omp critical
#endif
#endif
        cout<<setw(10)<<fixed<<setprecision(5)<<i*scal*100<<"%";
        cout<<" ; Number of critical points = "<<_criticalPoints.size()<<endl;
        for(int j=0;j<_domain.get_N2();j++)
        {
            for(int k=0;k<_domain.get_N3();k++)
            {
                vector<vector<int>> listOfVisitedPoints;
                current[0] = i;
                current[1] = j;
                current[2] = k;
                if (_volumeNumberOfPoints[i][j][k] != 0) continue;
                if (_V[i][j][k][0]<TOL) continue;

                listOfVisitedPoints =assentTrajectory(current, ongrid,false);

#ifdef ENABLE_OMP
#pragma omp critical
#endif
                if (_volumeNumberOfPoints [current[0]][current[1]][current[2]] != 0)
                {
                    int icp = abs(_volumeNumberOfPoints[current[0]][current[1]][current[2]]);
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;

                }
                else
                {
                    numberOfCriticalPoints++;
                    vector<int> data;
                    int icp = numberOfCriticalPoints;
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;

                    _volumeNumberOfPoints[current[0]][current[1]][current[2]] = -icp;
                    CriticalPoint  cp = newCriticalPoint(current[0], current[1], current[2],numberOfCriticalPoints);
                    cp.value= _V[current[0]][current[1]][current[2]][0];
                    _criticalPoints.insert(_criticalPoints.begin(),cp);
                }

            }
        }
    }
}
/**************************************************************************/
void GridCP::buildBasins(const Grid& grid, PartitionMethod partitionMethod)
{
    bool ongrid = (partitionMethod == PartitionMethod::AIM_ON_GRID
                   || partitionMethod == PartitionMethod::VDD
                   || partitionMethod == PartitionMethod::BECKE);
    
    initGridCP(grid, ongrid);

    if (partitionMethod == PartitionMethod::AIM_ON_GRID
        || partitionMethod == PartitionMethod::AIM_NEAR_GRID
        || partitionMethod == PartitionMethod::AIM_NEAR_GRID_REFINEMENT)
    {
        assignPointsByGradient(ongrid);
    }

    if (partitionMethod == PartitionMethod::AIM_NEAR_GRID_REFINEMENT) 
    {
        int itermax=30;
        int N = _domain.get_N1() * _domain.get_N2() * _domain.get_N3() / 1000;
        int nold = refineEdgeNearGrad();
        int n = refineEdgeNearGrad();
        int iter = 0;
        //cout<<"n-nold="<<abs(n-nold)<<" => N max="<<N<<endl;
        while(abs(n - nold) > N)
        {
            iter++;
            if (iter > itermax)
            {
                break;
            }

            nold = n;
            n = refineEdgeNearGrad();
            //cout<<"n-nold="<<abs(n-nold)<<" => N max="<<N<<endl;
        }
    }

    if (partitionMethod == PartitionMethod::AIM_ON_GRID
        || partitionMethod == PartitionMethod::AIM_NEAR_GRID
        || partitionMethod == PartitionMethod::AIM_NEAR_GRID_REFINEMENT) 
    {
        resetKnown();
        //cout<<"Begin computeNumCenters"<<endl;
        computeNumCenters();
        //cout<<"End computeNumCenters"<<endl;
        computeVolumes();
        //cout<<"End computeVolumes"<<endl;
        removeBasins0();
        //cout<<"End removeBasins0"<<endl;
        removeNonSignificantBasins();
        computeVolCenters();
        //cout<<"End removeNonSignificantBasins"<<endl;
        
        int nt=0;
        for(size_t ip = 0; ip < _criticalPoints.size(); ++ip)
        {
            int i = int(_criticalPoints[ip].volume / _domain.get_dv());
            nt += i;
            cout << "Point n " << ip << " nPoints=" << i << endl;
        }

        cout << "ntot=" << nt << endl;
        cout << "N1*N2*n3=" << _domain.get_N1() * _domain.get_N2() * _domain.get_N3() << endl;
    }

    if (partitionMethod == PartitionMethod::VDD)
    {
        buildVDD();
    }

}
/**************************************************************************/
void GridCP::printCriticalPoints()
{
    cout<<"Number of critical points = "<<_criticalPoints.size()<<endl;
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        CriticalPoint cp = _criticalPoints[ip];
        int i = cp.index[0];
        int j = cp.index[1];
        int k = cp.index[2];
        int I = cp.volCenter[0];
        int J = cp.volCenter[1];
        int K = cp.volCenter[2];
        cout<<" Index = "<<left<<setw(10)<<i<<", "<<setw(10)<<j<<", "<<setw(10)<<k
        <<" value = "<<setw(15)<<cp.value<<" Z="<< setw(15)<<cp.nuclearCharge
        <<" Integral = "<<setw(15)<< cp.integral<<", "<<right<<"  Volume center = "<<left<<setw(10)<<I<<", "<<setw(10)<<J<<", "<<setw(10)<<K<<"; "<<" x = "<<setw(10)<<_domain.x(I,J,K)<<", "<<" y= "<<setw(10)<<_domain.y(I,J,K)<<", "<<" z = "<<setw(10)<<_domain.z(I,J,K)<<endl;
    }
}
/**************************************************************************/
void GridCP::computeIntegrals(const Grid& grid)
{
    double dv=_domain.get_dv() ;
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        _criticalPoints[ip].integral=0;
        for(int i=0;i<_domain.get_N1();i++)
        for(int j=0;j<_domain.get_N2();j++)
        for(int k=0;k<_domain.get_N3();k++)
        {
            int n=abs(_volumeNumberOfPoints [i][j][k]);
            if (n==_criticalPoints[ip].numVolume)
            {
                _criticalPoints[ip].integral+=dv*grid.value(i,j,k);
            }
        }
    }
}
/**************************************************************************/
vector<double> GridCP::computeAIMCharges(const Grid& grid)
{
    computeIntegrals(grid);
    cout<<"Partial charges (Pos in Angstrom) "<<endl;
    double r[3];
    double totalCharge = 0;
    vector<double> PartCharges(_criticalPoints.size());
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        CriticalPoint cp = _criticalPoints[ip];
        double PartCharge= -cp.integral+cp.nuclearCharge;
        int i = cp.index[0];
        int j = cp.index[1];
        int k = cp.index[2];
        double x = _domain.x(i, j, k) * Constants::BOHR_RADIUS_TO_ANGSTROM;
        double y = _domain.y(i, j, k) * Constants::BOHR_RADIUS_TO_ANGSTROM;
        double z = _domain.z(i, j, k) * Constants::BOHR_RADIUS_TO_ANGSTROM;
        int ia=_criticalPoints[ip].numCenter;
        for(int c=0; c<3;c++)
            r[c] = _str.atom(ia).get_coordinates()[c] * Constants::BOHR_RADIUS_TO_ANGSTROM;

        cout<<" Center = "<<fixed<<setprecision(6)<<right<<setw(10)<<x<<", "<<setw(10)<<y<<", "<<setw(10)<<z;
        cout<<left<<" Atom = "<<right<<setw(10)<<r[0]<<", "<<setw(10)<<r[1]<<", "<<setw(10)<<r[2]<<", "<< setw(10)<<cp.nuclearCharge;
        cout<<left<<" Integral = "<<setw(10)<<cp.integral;
        cout<<left<<" Charge = "<<fixed<<setprecision(10)<< setw(15)<<PartCharge <<endl;
        PartCharges[ia]=PartCharge;
        totalCharge +=  PartCharge;
    }
    cout<<"Total charge = "<<totalCharge<<endl;
    for(size_t ia=0;ia<PartCharges.size();ia++)
    {
        cout<<"ia "<<ia<<" Q "<<PartCharges[ia]<<endl;
    }
    return PartCharges;
    
}

Structure GridCP::str() const
{
    return _str;
}

// Ref : Fonseca Guerra et al. https://doi.org/10.1002/jcc.10351
/**************************************************************************/
void GridCP::buildVDD()
{
    PeriodicTable pTable;
    double scal;
    string str ="Assigning points to volumes... Please wait";
    int N[3]={_domain.get_N1(),_domain.get_N2(), _domain.get_N3()};

    if (!okDomain()) return;

    if (_criticalPoints.size()>0)
        _criticalPoints.clear();

    int nAtoms=_str.number_of_atoms();
    //cout<<"nAtoms="<<nAtoms<<endl;
    for(int ia=0;ia<nAtoms;ia++)
    {
            CriticalPoint  cp = newCriticalPoint(0,0,0,ia+1);
            cp.value= 0;
            _criticalPoints.push_back(cp);
    }
    vector< double > V(3,0);
    vector< vector<double> > coordinates(nAtoms,V);
    vector< double > rcov2(nAtoms,0);
    for(int ia=0;ia<nAtoms;ia++)
    {
            //double rcov =   _str.atom(ia).element().covalent_radii();
            double rcov =   _str.atom(ia).get_element().radii();
            rcov2[ia] = rcov*rcov;
            //rcov2[ia] = 0.5*0.5;
            for(int c=0;c<3;c++)
                coordinates[ia][c] = _str.atom(ia).get_coordinates()[c];
    }
    vector< double > r2minAtom(nAtoms,-1);

    scal = 1.0/(_domain.get_N1()-1);
    cout<<str<<endl;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
    for(int i=0;i<N[0];i++)
    {
        cout<<setw(10)<<fixed<<setprecision(5)<<i*scal*100<<"%"<<endl;
        for(int j=0;j<N[1];j++)
        {
            for(int k=0;k<N[2];k++)
            {
                int iamin=-1;
                double r2min=-1;
                _volumeNumberOfPoints[i] [j] [k] = 0;
                double r[] = {_domain.x(i,j,k), _domain.y(i,j,k), _domain.z(i,j,k)};
                for(int ia=0;ia<nAtoms;ia++)
                {
                    double r2 = 0;
                    for(int c=0;c<3;c++)
                        r2 += (r[c]-coordinates[ia][c])*(r[c]-coordinates[ia][c]);
                    r2 /= rcov2[ia];
                    if (iamin==-1 || r2<r2min)
                    {
                        r2min= r2;
                        iamin=ia;
                    }
                }
                if (iamin<0) 
                    continue;
                _volumeNumberOfPoints[i] [j] [k] = iamin+1;
                if (r2minAtom[iamin]<0 || r2min<r2minAtom[iamin])
                {
                    r2minAtom[iamin]=r2min;
                    _criticalPoints[iamin].value = _V[i][j][k][0];
                    _criticalPoints[iamin].index[0] = i;
                    _criticalPoints[iamin].index[1] = j;
                    _criticalPoints[iamin].index[2] = k;

                }
            }
        }
    }
    resetKnown();
    computeNumCenters();
    computeVolumes();
    computeVolCenters();
    int nt=0;
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        int i=int(_criticalPoints[ip].volume/_domain.get_dv());
        nt+=i;
        cout<<"Point n "<<ip<<" nPoints="<<i<<endl;
    }
    cout<<"ntot="<<nt<<endl;
    cout<<"N1*N2*n3="<<_domain.get_N1()*_domain.get_N2()*_domain.get_N3()<<endl;
}
/**************************************************************************/
bool GridCP::addSurroundingSignPoints(int current[], vector<vector<int>>& listOfVisitedPoints, double cutoff)
{
    int i = current[0];
    int j = current[1];
    int k = current[2];
    // Initialize indices for neighboring points
    int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
    int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
    int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };

    if (!okDomain()) return false;

    double v0 = _V[I[1]][J[1]][K[1]][0];
    _known[i][j][k]=1;
    for(size_t ic=0;ic<3;ic++)
    for(size_t jc=0;jc<3;jc++)
    for(size_t kc=0;kc<3;kc++)
    {
        if (ic==1 && jc==1 && kc==1) continue;
        if (_known[I[ic]][J[jc]][K[kc]]>0) continue;
        double v =_V[I[ic]][J[jc]][K[kc]][0];
        //if (abs(v)<TOL) continue;
        if (abs(v)>cutoff && v*v0>0)
        {
            vector<int>  data = { I[ic], J[jc], K[kc]};
            listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);
            _known[I[ic]][J[jc]][K[kc]] = 1;
        }
    }
    return true;
}
/**************************************************************************/
void GridCP::buildBasinsBySign(const Grid& grid, double cutoff)
{
    initGridCP(grid, 0);
    if (abs(cutoff)<TOL) cutoff=TOL;
    cutoff=abs(cutoff);
    //double scal;
    string str ="Assigning points to volumes... Please wait";

    if (!okDomain()) return;

    if (_criticalPoints.size()>0)
        _criticalPoints.clear();

    resetKnown();
    //scal = 1.0/(_domain.N1()-1);
    cout<<str<<endl;
    int numberOfCriticalPoints=0;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
    for(int i=0;i<_domain.get_N1();i++)
    {
        int current[3];
        //cout<<setw(10)<<fixed<<setprecision(5)<<i*scal*100<<"%";
        for(int j=0;j<_domain.get_N2();j++)
        {
            for(int k=0;k<_domain.get_N3();k++)
            {
                double v0 = _V[i][j][k][0];
                if (abs(v0)<cutoff) continue;
                vector<vector<int>> listOfVisitedPoints;
                current[0] = i;
                current[1] = j;
                current[2] = k;
                vector<int>  data = { current[0], current[1], current[2]};
                listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);
                addSurroundingSignPoints(current, listOfVisitedPoints,cutoff);
                //cout<<"Number of visited points ="<<listOfVisitedPoints.size()<<endl;
                if (_volumeNumberOfPoints [current[0]][current[1]][current[2]] != 0)
                {
                    int icp = _volumeNumberOfPoints[current[0]][current[1]][current[2]];
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;

                }
                else
                {
                    numberOfCriticalPoints++;
                    int icp = -numberOfCriticalPoints;
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;
                }

            }
        }
    }
    resetKnown();
    int n=0;
    for(int i=0;i<_domain.get_N1();i++)
            for(int k=0;k<_domain.get_N3();k++)
                for(int j=0;j<_domain.get_N2();j++)
            {
                double v0 = _V[i][j][k][0];
                if (abs(v0)<cutoff) continue;
                vector<vector<int>> listOfVisitedPoints;
                int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
                int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
                int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };
                for(size_t ic=0;ic<3;ic++)
                for(size_t jc=0;jc<3;jc++)
                for(size_t kc=0;kc<3;kc++)
                {
                    if (ic==1 && jc==1 && kc==1) continue;
                    double v =_V[I[ic]][J[jc]][K[kc]][0];
                    if (abs(v)>cutoff && v*v0>0 && _volumeNumberOfPoints[i][j][k] != _volumeNumberOfPoints[I[ic]][J[jc]][K[kc]])
                    {
                        vector<int>  data = { I[ic], J[jc], K[kc]};
                        listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);
                    }
                }
                if (listOfVisitedPoints.size()>0) 
                {
                    n+=listOfVisitedPoints.size();
                    int icp = _volumeNumberOfPoints[i][j][k];
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;
                }
            }
    cout<<"Number of changed points : ======>"<<n<<endl;
    n=0;
    for(int k=0;k<_domain.get_N3();k++)
        for(int i=0;i<_domain.get_N1();i++)
            for(int j=0;j<_domain.get_N2();j++)
            {
                double v0 = _V[i][j][k][0];
                if (abs(v0)<cutoff) continue;
                vector<vector<int>> listOfVisitedPoints;
                int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
                int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
                int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };
                for(size_t ic=0;ic<3;ic++)
                for(size_t jc=0;jc<3;jc++)
                for(size_t kc=0;kc<3;kc++)
                {
                    if (ic==1 && jc==1 && kc==1) continue;
                    double v =_V[I[ic]][J[jc]][K[kc]][0];
                    if (abs(v)>cutoff && v*v0>0 && _volumeNumberOfPoints[i][j][k] != _volumeNumberOfPoints[I[ic]][J[jc]][K[kc]])
                    {
                        vector<int>  data = { I[ic], J[jc], K[kc]};
                        listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);
                    }
                }
                if (listOfVisitedPoints.size()>0) 
                {
                    n+=listOfVisitedPoints.size();
                    //cout<<"Number of visited points======>"<<listOfVisitedPoints.size()<<endl;
                    int icp = _volumeNumberOfPoints[i][j][k];
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;
                }
            }
    cout<<"Number of changed points : ======>"<<n<<endl;
    n=0;
        for(int i=0;i<_domain.get_N1();i++)
            for(int j=0;j<_domain.get_N2();j++)
    for(int k=0;k<_domain.get_N3();k++)
            {
                double v0 = _V[i][j][k][0];
                if (abs(v0)<cutoff) continue;
                vector<vector<int>> listOfVisitedPoints;
                int I[3] = { max(0, i-1), i, min(i+1, _domain.get_N1()-1) };
                int J[3] = { max(0, j-1), j, min(j+1, _domain.get_N2()-1) };
                int K[3] = { max(0, k-1), k, min(k+1, _domain.get_N3()-1) };
                for(size_t ic=0;ic<3;ic++)
                for(size_t jc=0;jc<3;jc++)
                for(size_t kc=0;kc<3;kc++)
                {
                    if (ic==1 && jc==1 && kc==1) continue;
                    double v =_V[I[ic]][J[jc]][K[kc]][0];
                    if (abs(v)>cutoff && v*v0>0 && _volumeNumberOfPoints[i][j][k] != _volumeNumberOfPoints[I[ic]][J[jc]][K[kc]])
                    {
                        vector<int>  data = { I[ic], J[jc], K[kc]};
                        listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);
                    }
                }
                if (listOfVisitedPoints.size()>0) 
                {
                    n+=listOfVisitedPoints.size();
                    //cout<<"Number of visited points======>"<<listOfVisitedPoints.size()<<endl;
                    int icp = _volumeNumberOfPoints[i][j][k];
                    for(size_t ip=0;ip<listOfVisitedPoints.size();ip++)
                        _volumeNumberOfPoints[listOfVisitedPoints[ip][0]] [listOfVisitedPoints[ip][1]] [listOfVisitedPoints[ip][2]] = icp;
                }
            }
    cout<<"Number of changed points : ======>"<<n<<endl;



    cout<<"build CP "<<endl;
    for(int i=0;i<_domain.get_N1();i++)
        for(int j=0;j<_domain.get_N2();j++)
            for(int k=0;k<_domain.get_N3();k++)
            {
                int icp=_volumeNumberOfPoints[i][j][k];
                if (icp==0) continue;
                //cout<<"CP size = "<<_criticalPoints.size()<<endl;
                if (_criticalPoints.size()==0)
                {
                    CriticalPoint  cp = newCriticalPoint(i,j,k,icp);
                    cp.value= _V[i][j][k][0];
                    _criticalPoints.insert(_criticalPoints.begin(),cp);
                }
                else {
                    bool ok=true;
                    for(size_t ip=0;ip<_criticalPoints.size();ip++)
                    {
                        if (icp==_criticalPoints[ip].numVolume)
                        {
                            ok = false;
                            break;
                        }
                    }
                    if (ok)
                    {
                        CriticalPoint  cp = newCriticalPoint(i,j,k,icp);
                        cp.value= _V[i][j][k][0];
                        _criticalPoints.insert(_criticalPoints.begin(),cp);
                    }
                    
                }
            }
    // resort number
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
            int icp = _criticalPoints[ip].numVolume;
            int nv = ip+1;
            double maxv=0;
            for(int i=0;i<_domain.get_N1();i++)
                for(int j=0;j<_domain.get_N2();j++)
                    for(int k=0;k<_domain.get_N3();k++)
                        if (_volumeNumberOfPoints[i][j][k]==icp)
                        {
                            _volumeNumberOfPoints[i][j][k] = nv;
                            if (abs(_V[i][j][k][0])>abs(maxv))
                                maxv= _V[i][j][k][0];
                        }
             _criticalPoints[ip].numVolume=nv;
             _criticalPoints[ip].value=maxv;
    }
    resetKnown();
    computeNumCenters();
    computeVolumes();
    computeVolCenters();
    int nt=0;
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        int i=int(_criticalPoints[ip].volume/_domain.get_dv());
        nt+=i;
        cout<<"Point n "<<ip<<" nPoints="<<i<<endl;
    }
    cout<<"ntot="<<nt<<endl;
    cout<<"N1*N2*n3="<<_domain.get_N1()*_domain.get_N2()*_domain.get_N3()<<endl;
}
/**************************************************************************/
void GridCP::build2BasinSign(const Grid& grid)
{
    initGridCP(grid, 0);
    //double scal;
    string str ="Assigning points to volumes... Please wait";
    int N[3]={_domain.get_N1(),_domain.get_N2(), _domain.get_N3()};

    if (!okDomain()) return;

    if (_criticalPoints.size()>0)
        _criticalPoints.clear();

    //int nAtoms=_str.number_of_atoms();
    //cout<<"nAtoms="<<nAtoms<<endl;
    CriticalPoint  cp;
    cp = newCriticalPoint(0,0,0,1);
    cp.value= 0;
    _criticalPoints.push_back(cp);
    cp = newCriticalPoint(0,0,0,2);
    cp.value= 0;
    _criticalPoints.push_back(cp);

    //scal = 1.0/(_domain.N1()-1);
    cout<<str<<endl;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
    for(int i=0;i<N[0];i++)
    {
        //cout<<setw(10)<<fixed<<setprecision(5)<<i*scal*100<<"%"<<endl;
        for(int j=0;j<N[1];j++)
        {
            for(int k=0;k<N[2];k++)
            {
                if (_V[i][j][k][0]<0)
                {
                    _volumeNumberOfPoints[i][j][k] = 1;
                    if (_criticalPoints[0].value>_V[i][j][k][0]) 
                        _criticalPoints[0].value=_V[i][j][k][0]; 
                }
                else
                if (_V[i][j][k][0]>0)
                {
                    _volumeNumberOfPoints[i] [j] [k] = 2;
                    if (_criticalPoints[1].value<_V[i][j][k][0]) 
                        _criticalPoints[1].value=_V[i][j][k][0]; 
                }
            }
        }
    }
    resetKnown();
    computeNumCenters();
    computeVolCenters();
    computeVolumes();
    int nt=0;
    for(size_t ip=0;ip<_criticalPoints.size();ip++)
    {
        int i=int(_criticalPoints[ip].volume/_domain.get_dv());
        nt+=i;
        cout<<"Point n "<<ip<<" nPoints="<<i<<endl;
    }
    cout<<"ntot="<<nt<<endl;
    cout<<"N1*N2*n3="<<_domain.get_N1()*_domain.get_N2()*_domain.get_N3()<<endl;
}
