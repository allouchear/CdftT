using namespace std;
#include <numeric/GridCP.h>
#include <common/Constants.h>
#include <common/Structure.h>
#include <cmath>
#include <iomanip>
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
	cp.rho = 0;
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
	vector<int> V(_domain.N3(),0);
	vector<vector<int>> M(_domain.N2(), V);
	vector<vector<vector<int>>> V3D(_domain.N1(), M);
	return V3D;
}
/**************************************************************************/
bool GridCP::okDomain()
{
	bool ok=false;
	if(_domain.N1()<2) return ok;
	if(_domain.N2()<2) return ok;
	if(_domain.N3()<2) return ok;
	if(_domain.Nval()<1) return ok;
	return true;
}
/**************************************************************************/
void GridCP::initGridCP(const Grid& grid, bool ongrid)
{
	reset();
	_domain = grid.dom();
	if(!okDomain()) return;
	_str = grid.str();
	_V = grid.V();
	int nval = _V[0][0][0].size();
	if(nval<4 && !ongrid)
	{
		cout<<"begin grad"<<endl;
		Grid g= grid.gradient(1);
		cout<<"end grad"<<endl;
		_V = g.V();
		_domain = g.dom();
		_str = g.str();
	}
	if(!ongrid)
	for(size_t i=0;i<_V.size();i++)
	for(size_t j=0;j<_V[i].size();j++)
	for(size_t k=0;k<_V[i][j].size();k++)
	{
		double c= abs(_V[i][j][k][1]);
		if(c < abs(_V[i][j][k][2])) c= abs(_V[i][j][k][2]);
		if(c < abs(_V[i][j][k][3])) c= abs(_V[i][j][k][3]);
		if(c>0)
		{
			c = 1.0/c;
			for(size_t l = 1;l<=3;l++)
				_V[i][j][k][l] *= c;
		}
	}
	if(!ongrid)
	for(size_t i=1;i<_V.size()-1;i++)
	for(size_t j=1;j<_V[i].size()-1;j++)
	for(size_t k=1;k<_V[i][j].size()-1;k++)
	{
		if( _V[i-1][j][k][0]<_V[i][j][k][0] && _V[i+1][j][k][0]<_V[i][j][k][0]) _V[i][j][k][1]  = 0;
		if( _V[i][j-1][k][0]<_V[i][j][k][0] && _V[i][j+1][k][0]<_V[i][j][k][0]) _V[i][j][k][2]  = 0;
		if( _V[i][j][k-1][0]<_V[i][j][k][0] && _V[i][j][k+1][0]<_V[i][j][k][0]) _V[i][j][k][3]  = 0;
	}

	_volumeNumberOfPoints = get3DIntVector();
	_known = get3DIntVector();
}
/**************************************************************************/
bool GridCP::isVolumeEdge(int current[])
{
	if(!okDomain()) return 0;
	int i = current[0];
	int j = current[1];
	int k = current[2];
	// Initialize indices for neighboring points
	int I[3] = { max(0, i-1), i, min(i+1, _domain.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _domain.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _domain.N3()-1) };

	for(size_t ic=0;ic<3;ic++)
	for(size_t jc=0;jc<3;jc++)
	for(size_t kc=0;kc<3;kc++)
		if(2!=_known[I[ic]][J[jc]][K[kc]] && abs(_volumeNumberOfPoints[i][j][k]) != abs(_volumeNumberOfPoints[I[ic]][J[jc]][K[kc]]))
			return true;

	return false;
}
/**************************************************************************/
int GridCP::setArroundTo(int current[], int kn)
{
	if(!okDomain()) return 0;

	int i = current[0];
	int j = current[1];
	int k = current[2];
	// Initialize indices for neighboring points
	int I[3] = { max(0, i-1), i, min(i+1, _domain.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _domain.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _domain.N3()-1) };
	int n = 0;


	for(size_t ic=0;ic<3;ic++)
	for(size_t jc=0;jc<3;jc++)
	for(size_t kc=0;kc<3;kc++)
	{
		if(ic==1 && jc==1 && kc ==1) continue;
		if(_known[I[ic]][J[jc]][K[kc]] != 1 && _volumeNumberOfPoints[I[ic]][J[jc]][K[kc]]==_volumeNumberOfPoints[i][j][k]) 
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
	if(!okDomain()) return false;
	int i = current[0];
	int j = current[1];
	int k = current[2];
	// Initialize indices for neighboring points
	int I[3] = { max(0, i-1), i, min(i+1, _domain.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _domain.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _domain.N3()-1) };


	for(size_t ic=0;ic<3;ic++)
	for(size_t jc=0;jc<3;jc++)
	for(size_t kc=0;kc<3;kc++)
		if(_V[I[ic]][J[jc]][K[kc]][0]>_V[I[1]][J[1]][K[1]][0]) return false;

	return true;
}
/**************************************************************************/
bool GridCP::nextPointOnGrid(int current[], int next[])
{
	if(!okDomain()) return false;
	int i = current[0];
	int j = current[1];
	int k = current[2];
	// Initialize indices for neighboring points
	int I[3] = { max(0, i-1), i, min(i+1, _domain.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _domain.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _domain.N3()-1) };

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
		if(ic==1 && jc==1 && kc==1) continue;
		if(_known[I[ic]][J[jc]][K[kc]] >1) continue;
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
		if(drho>drhoMax)
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
	if(!okDomain()) return false;
	int nval = _V[0][0][0].size();
	if(nval<4) return false;
	double gradrl[3];
	// Initialize indices for neighboring points
	int i = current[0];
	int j = current[1];
	int k = current[2];
	int N[3] ={ _domain.N1(), _domain.N2(), _domain.N3()};

	for(size_t c=0;c<3;c++) gradrl[c] = _V[i][j][k][c+1];
	if(int(rint(gradrl[0])) ==0 && int(rint(gradrl[1])) ==0 && int(rint(gradrl[2])) ==0)
	{
		if(isMax(current))
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
		for(size_t c=0;c<3;c++) if(next[c]<0 ) next[c] = 0;
		for(size_t c=0;c<3;c++) if(next[c]>N[c]-1) next[c] = N[c]-1;

		i = current[0];
		j = current[1];
		k = current[2];
		_known[i][j][k] = 1;
		i = next[0];
		j = next[1];
		k = next[2];
		if(_known[i][j][k]==1)
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
	int I[3] = { max(0, i-1), i, min(i+1, _domain.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _domain.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _domain.N3()-1) };

	double rho0 = 0;
	double dRho = 0;

	if(!okDomain()) return false;

	rho0 = _V[I[1]][J[1]][K[1]][0];
	for(size_t ic=0;ic<3;ic++)
	for(size_t jc=0;jc<3;jc++)
	for(size_t kc=0;kc<3;kc++)
	{
		if(ic==1 && jc==1 && kc==1) continue;
		dRho =_V[I[ic]][J[jc]][K[kc]][0]-rho0;
		if(fabs(dRho)<TOL)
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
	if(!okDomain()) return listOfVisitedPoints;

	vector<int>  data = { current[0], current[1], current[2]};
	listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);

	size_t imax = _domain.N1()* _domain.N2()* _domain.N3();
	int next[3];
	double deltaR[3] = {0,0,0};
	for(size_t l=0;l<imax;l++)
	{

		if(ongrid)
			nextPointOnGrid(current, next);
		else
			nextPoint(deltaR, current, next);
		vector<int>  data = { next[0], next[1], next[2]};
		listOfVisitedPoints.insert (listOfVisitedPoints.begin(),data);

		if((next[0] == current[0] && next[1] == current[1] && next[2] == current[2])
		)
		{
			addSurroundingEqualPoints(current, listOfVisitedPoints);
			break;
		}
		if(!refine && _volumeNumberOfPoints [next[0]][next[1]][next[2]] != 0)
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
	for(int i=0;i<_domain.N1();i++)
	{
		for(int j=0;j<_domain.N2();j++)
		for(int k=0;k<_domain.N3();k++)
		{
			int n=abs(_volumeNumberOfPoints [i][j][k]);
			for(size_t ip=0;ip<_criticalPoints.size();ip++)
				if(_criticalPoints[ip].numVolume==n) 
					_criticalPoints[ip].volume +=1 ;
		}
	}
	for(size_t ip=0;ip<_criticalPoints.size();ip++)
			_criticalPoints[ip].volume *= _domain.dv() ;
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
			r[c] = _domain.O()[ c ] + i*_domain.T()[c][0] + j*_domain.T()[c][1] +  k*_domain.T()[c][2]; 
		for(int ia=0;ia<_str.number_of_atoms();ia++)
		{
			for(int c=0; c<3;c++)
				dr[c] = r[c]-_str.atom(ia).coordinates()[c];
			r2 = 0;
			for(int c=0; c<3;c++)
				r2 += dr[c]*dr[c];
			
			if(ia==0 || r2<r2old )
			{
				_criticalPoints[ip].numCenter=ia;
				_criticalPoints[ip].nuclearCharge=_str.atom(ia).atomic_number();
				r2old = r2;
			}
		}
	}
}
/**************************************************************************/
void GridCP::removeAttractor0()
{
	vector< CriticalPoint > newCriticalPoints;
	for(size_t ip=0;ip<_criticalPoints.size();ip++)
	{
		CriticalPoint  cp = _criticalPoints[ip];
		int i = cp.index[0];
		int j = cp.index[1];
		int k = cp.index[2];
		if(_V[i][j][k][0]<TOL)
			_volumeNumberOfPoints [i][j][k] = 0;
		else
			newCriticalPoints.push_back(cp);
	}
	_criticalPoints= newCriticalPoints;
}
/**************************************************************************/
void GridCP::removeNonSignificantAttractor()
{
	int n = _domain.N1()* _domain.N2()* _domain.N3()/100;
	vector< CriticalPoint > newCriticalPoints;
	for(size_t ip=0;ip<_criticalPoints.size();ip++)
	{

		CriticalPoint cp = _criticalPoints[ip];
		//cout<<" CP # "<<ip+1<<" npts="<< cp.volume/_domain.dv()<<" nmin="<<n<<" rho= "<<cp.rho<<endl;
		if(cp.volume/_domain.dv()>=n) 
			newCriticalPoints.push_back(cp);
	}
	_criticalPoints= newCriticalPoints;
}
/**************************************************************************/
int GridCP::refineEdge()
{
	if(!okDomain()) return 0;
	size_t N[3] ={ size_t(_domain.N1()), size_t(_domain.N2()), size_t(_domain.N3())};

	string str ="Refining grid points adjacent to Bader surface... Please wait";
	int current[3];
	int ne = 0;

	resetKnown();

	for(size_t i=0;i<N[0];i++)
	for(size_t j=0;j<N[1];j++)
	for(size_t k=0;k<N[2];k++)
	{
		if(_V[i][j][k][0]<TOL) _known[i][j][k] = 2;
		if(_volumeNumberOfPoints [i][j][k]<0) _known[i][j][k] = 2;
	}
	for(size_t i=0;i<N[0];i++)
	for(size_t j=0;j<N[1];j++)
	for(size_t k=0;k<N[2];k++)
	{
		current[0] = i;
		current[1] = j;
		current[2] = k;
		if(_volumeNumberOfPoints[i][j][k]>0 && !isMax(current) && isVolumeEdge(current)) 
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
		if(_known[i][j][k]==1) 
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
				if(_volumeNumberOfPoints [i][j][k]!=0) continue;
				if(_known[i][j][k] !=0 ) continue;

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
void GridCP::assignGridCP(bool ongrid)
{
	double scal;
	string str ="Assigning points to volumes... Please wait";

	if(!okDomain()) return;

	if(_criticalPoints.size()>0)
		_criticalPoints.clear();

	resetKnown();
	int numberOfCriticalPoints = 0;

	scal = 1.0/(_domain.N1()-1);
	cout<<str<<endl;
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
	for(int i=0;i<_domain.N1();i++)
	{
		int current[3];
#ifdef ENABLE_OMP
#ifndef G_OS_WIN32
#pragma omp critical
#endif
#endif
		cout<<setw(10)<<fixed<<setprecision(5)<<i*scal*100<<"%";
		cout<<" ; Number of critical points = "<<_criticalPoints.size()<<endl;
		for(int j=0;j<_domain.N2();j++)
		{
			for(int k=0;k<_domain.N3();k++)
			{
				vector<vector<int>> listOfVisitedPoints;
				current[0] = i;
				current[1] = j;
				current[2] = k;
				if(_volumeNumberOfPoints[i][j][k] != 0) continue;
				if(_V[i][j][k][0]<TOL) continue;

				listOfVisitedPoints =assentTrajectory(current, ongrid,false);

#ifdef ENABLE_OMP
#pragma omp critical
#endif
				if(_volumeNumberOfPoints [current[0]][current[1]][current[2]] != 0)
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
					cp.rho= _V[current[0]][current[1]][current[2]][0];
					_criticalPoints.insert(_criticalPoints.begin(),cp);
				}

			}
		}
	}
}
/**************************************************************************/
void GridCP::buildAttractors(const Grid& grid, int method)
{
	cout<<"done"<<endl;
	bool ongrid=(method==0);
	initGridCP(grid, ongrid);
	cout<<"init"<<endl;
	assignGridCP(ongrid);
	cout<<"assign"<<endl;
	if(method==2) 
	{
		cout<<"done"<<endl;
		int itermax=30;
		int N = _domain.N1()* _domain.N2()* _domain.N3()/1000;
		int nold = refineEdge();
		int n = refineEdge();
		int iter = 0;
		cout<<"done"<<endl;
		//cout<<"n-nold="<<abs(n-nold)<<" => N max="<<N<<endl;
		while(abs(n-nold)>N)
		{
			iter++;
			if(iter>itermax) break;
			nold = n;
			n = refineEdge();
			//cout<<"n-nold="<<abs(n-nold)<<" => N max="<<N<<endl;
		}
		cout<<"done"<<endl;
	}

	resetKnown();
	computeNumCenters();
	computeVolumes();
	removeAttractor0();
	removeNonSignificantAttractor();
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
		cout<<" Index = "<<i<<", "<<j<<", "<<k<<" rho = "<<cp.rho<<" Z="<< cp.nuclearCharge<<" Integral = "<< cp.integral<<endl;
	}
}
/**************************************************************************/
vector<double> GridCP::computeAIMCharges(const Grid& grid)
{
	double dv=_domain.dv() ;
	for(size_t ip=0;ip<_criticalPoints.size();ip++)
	{
		_criticalPoints[ip].integral=0;
		for(int i=0;i<_domain.N1();i++)
		for(int j=0;j<_domain.N2();j++)
		for(int k=0;k<_domain.N3();k++)
		{
			int n=abs(_volumeNumberOfPoints [i][j][k]);
			if(n==_criticalPoints[ip].numVolume)
			{
				_criticalPoints[ip].integral+=dv*grid.value(i,j,k);
			}
		}
	}
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
		double x = _domain.x(i,j,k)*BOHRTOANG;
		double y = _domain.y(i,j,k)*BOHRTOANG;
		double z = _domain.z(i,j,k)*BOHRTOANG;
		int ia=_criticalPoints[ip].numCenter;
		for(int c=0; c<3;c++)
			r[c]=_str.atom(ia).coordinates()[c]*BOHRTOANG;
		
		cout<<" Center = "<<fixed<<setprecision(6)<<x<<", "<<y<<", "<<z;
		cout<<" Atom = "<<r[0]<<", "<<r[1]<<", "<<r[2]<<", "<< cp.nuclearCharge;
		cout<<" Charge = "<<fixed<<setprecision(10)<< PartCharge <<endl;
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

