#include<iostream>
#include<cmath>
#include"GTF.h"

using namespace std;

GTF::GTF(double exposant, double coefficient, vector<double> coord, vector<double> l) : _exposant(exposant), _coefficient(coefficient), _coord(coord), _l(l) {}

GTF::operator+ (vector<double> left, vector<double> right)
{
	vector<double> result(left.size());
	for(int i=0; i<result.size(); i++)
		result[i]=left[i]+right[i];
	return result;
}

GTF::operator* (GTF left, GTF right)
{
	return GTF(left.exposant()+right.exposant(), left.coefficient()*right.coefficient(), left.coord(), left.l()+right.l());
}

double GTF::normeGTF(Factorial FTable)
{
	return 2*_exposant/M_PI * sqrt(2*_exposant/M_PI * pow(4*_exposant, _l[0]+_l[1]+_l[2])/(FTable.double_factorial(2*_l[0]+1) * FTable.double_factorial(2*_l[1]+1) * FTable;double_factorial(2*_l[2]+1)));
}

void GTF::normalise_radialGTF(Factorial FTable)
{
	vector<double> l_bis(_l.size());
	l_bis=_l[0]+_l[1]+_l[2];
	_l[0]=l_bis;
	_l[1]=0;
	_l[2]=0;
	_coefficient*=normeGTF(FTable);
}

void GTF::normaliseGTF(Factorial FTable)
{
	_coefficient*=normeGTF(Ftable);
}

GTF GTF::overlapGTF(GTF left, GTF right)
{
	return //Voir avec le prof et vérifier la classe avec lui
}