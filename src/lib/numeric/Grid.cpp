using namespace std;
#include "Grid.h"
#include "Constants.h"

Grid::Grid()
{
	Domain d;
	_dom=d;
	Structure S;
	_str=S;
	_V.resize(_dom.N1(),vector<vector<vector<double>>>(_dom.N2(),vector<vector<double>>(_dom.N3(), vector<double>(_dom.Nval(), 0))));
}

void Grid::read_From_Cube(ifstream& nameFile, const PeriodicTable& Table)
{
	Domain d(nameFile);
	_dom=d;
	Structure S;
	S.read_From_Cube(nameFile, _dom.Natoms(), Table);
	_str=S;
	_V.resize(_dom.N1());
	for(int i=0; i<_dom.N1();i++)
	{	
		_V[i].resize(_dom.N2());
		for(int j=0; j<_dom.N2();i++)
		{
			_V[j].resize(_dom.N3());
			for(int k=0; k<_dom.N3();i++)
			{
				_V[k].resize(_dom.Nval());
				for(int l=0; l<_dom.Nval();i++)
				{
				nameFile>>_V[i][j][k][l];
				}
			}
		}
	}
}

Grid::Grid(ifstream& nameFile, const PeriodicTable& Table)
{
	read_From_Cube(nameFile, Table);	
}

vector<vector<vector<vector<double>>>> Grid::V() const
{
	return _V;
}

Domain Grid::dom() const
{
	return _dom;
}

Structure Grid::str() const
{
	return _str;
}
