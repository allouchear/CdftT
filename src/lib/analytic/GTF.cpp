#include<iostream>
#include"GTF.h"

using namespace std;

GTF::GTF(double exposant, double coefficient, vector<double> coord, vector<double> l, Factorial& F) : _exposant(exposant), 
		_coefficient(coefficient), _coord(coord), _l(l), _fact(F) {}

GTF::operator* (GTF left, GTF right)
{
	int i,j;
	vector<double> sum(3);
	double t;
	vector<double> PA(3);
	vector<double> PB(3);
	double gama=left.exposant()+right.exposant();
	double R2=0.0;
	double c=0;

	for(j=0; j<3; j++)
	{
		t=(left.exposant()*left.coord()[j] + right.exposant()*right.coord()[j])/gama;
		PA[j]=left.coord()[j]-t;
		PB[j]=right.coord()[j]-t;
		R2+=(left.coord()[j]-right.coord()[j])*(left.coord()[j]-right.coord()[j]);
	}

	c = (M_PI/gama)*sqrt(PI/gama)*exp(-left.exposant()*right.exposant/gama*R2);

	for(j=0; j<3; j++)
	{
		sum[j]=0.0;
		for(i=0; i<=(left.l()[j]+right.l()[j]/2); i++)
			sum[j]+=f(2*i, left.l()[j], right.l()[j], PA[j], PB[j])*_fact.double_factorial(2*i+1)/(power(gama,i));
	}
	return c*sum[0]*sum[1]*sum[2];
}

double GTF::normeGTF()
{
	return sqrt(2*_exposant/M_PI*sqrt(2*_exposant/M_PI)*power(4*_exposant, _l[0]+_l[1]+_l[2])
			/(_fact.double_factorial(_l[0])*_fact.double_factorial(_l[1])*_fact.double_factorial(_l[2])));
}

double GTF::normeGTF(GTF& q)
{
	return sqrt(2*q.exposant()/M_PI*sqrt(2*q.exposant()/M_PI)*power(4*q.exposant(), q.l()[0]+q.l()[1]+q.l()[2])
			/(_fact.double_factorial(q.l()[0])*_fact.double_factorial(q.l()[1])*_fact.double_factorial(q.l()[2])));

void GTF::normalise_radialGTF()
{
	GTF q(_exposant, _coefficient, _coord, _l)
	int l_bis=_l[0]+_l[1]+_l[2];
	q.l()[0]=l_bis;
	q.l()[1]=0;
	q.l()[2]=0;
	_coefficient*=normeGTF(GTF& q, _fact);
}

void GTF::normaliseGTF()
{
	_coefficient*=normeGTF(_fact);
}

GTF GTF::overlapGTF(GTF left, GTF right)
{
	return left.coefficient()*right.coefficient()*left*right;
}

GTF GTF::overlap3GTF(GTF left, GTF mid, GTF right)
{
	return left.coefficient()*mid.coefficient()*right.coefficient()*left*mid*right;
}

GTF GTF::overlap4GTF(GTF A, GTF B, GTF C, GTF D)
{
	return A.coefficient()*B.coefficient()*C.coefficient()*D.coefficient()*A*B*C*D;
}

double GTFxyzGTF(GTF p, GTF q, int ix, int iy, int iz)
{
	vector<double> C(3,0);
	vector<double> l {ix, iy, iz};
	GTF m(0.0, 1.0, l, C);
	return overlap3GTF(p, m, q);
}

double GTF::kineticGTF(GTF left, GTF right)
{
	int j;
	GTF a,b;
	vector<double> Ti(7);
	vector<double> sum(3);
	double T=0.0;

	for(j=0; j<7; j++)
		Ti[j]=0.0;
	a=left;
	b=right;
	Ti[0]=a*b;

	b.l()[0]=right.l()[0]+2;
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2];
	Ti[1]=a*b;

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1]+2;
	b.l()[2]=right.l()[2];
	Ti[2]=a*b;

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2]+2;
	Ti[3]=a*b;

	b.l()[0]=right.l()[0]-2;
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2];
	Ti[4]=a*b;

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1]-2;
	b.l()[2]=right.l()[2];
	Ti[5]=a*b;

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2]-2;
	Ti[6]=a*b;

	sum[0]=right.exposant()*(2*(right.l()[0]+right.l()[1]+right.l()[2])+3)*Ti[0];
	sum[1]=0.0;
	for(j=1; j<=3; j++)
		sum[1]+=Ti[j];
	sum[1]=-2*right.exposant()*right.exposant()*sum[1];

	sum[2]=0.0;
	for(j=4; j<=6; j++)
		sum[2]+=right.l()[j-4]*(right.l()[j-4]-1)*Ti[j];
	sum[2]*=-0.5;

	for(j=0; j<3; j++)
		T+=sum[j];
	return T*left.coefficient()*right.coefficient();

}

double ionicPotentialGTF(GTF left, GTF right, vector<double> C, double Z)
{
	int i,r,u;
	int j,s,n;
	int k,t,w;
	double Sx,Sy,Sz;
	double temp;
	vector<double> PA(3);
	vector<double> PB(3);
	vector<double> PC(3);
	double gama=left.exposant()+right.exposant();
	double R2=0.0;
	double PC2=0.0;
	double sum;
	//double FTable=NULL;										???????

	for(j=0; j<3; j++)
	{
		temp=(left.exopsant()*left.coord()[j]+right.exposant()*right.coord()[j])/gama;
		PA[j]=left.coord()[j]-temp;
		PB[j]=right.coord()[j]-temp;
		PC[j]=-C[j]+temp;
		R2+=(left.coefficient()[j]-right.coefficient()[j])*(left.coefficient()[j]-right.coefficient()[j]);
		PC2+=PC[j]*PC[j];
	}
	//FTable=getFTable(left.l()[0]+right.l()[0]+left.l()[1]+right.l()[1]+left.l()[2]+right.l()[2], gama*PC2);		????

	sum=0.0;
	for(i=0; i<=left.l()[0]+right.l()[0]; i++)
		for(r=0; r<=i/2; r++)
			for(u=0; u<=(i-2*r)/2; u++)
			{
				//Sx=A(i,r,u,left.l()[0],right.l()[0],PA[0],PB[0],PC[0],gama);
				for(j=0; j<=left.l()[1]+right.l()[1]; j++)
					for(s=0; s<=j/2; s++)
						for(n=0; n<=(j-2*s)/2;n++)
						{
							//Sy=A(j,s,n,left.l()[1],right.l()[1],PA[1],PB[1],PC[1],gama);
							for(k=0; k<=left.l()[2]+right.l()[2]; k++)
								for(t=0; t<=k/2; t++)
									for(w=0; w<=(k-2*t)/2; w++)
									{
										//Sz=A(k,t,w,left.l()[2],right.l()[2],PA[2],PB[2],PC[2],gama);
										//sum+=Sx*Sy*Sz*FTable[i+j+k-2*(r+s+t)-u-n-w];
									}
						}
			}
// Il faut coder la fonction A() et comprendre ce qu'est FTable (je n'ai pas l'impression que ce soit mon factorial table)
	sum *=2*PI/gama*exp(-left.exposant()*right.exposant()/gama*R2)*left.coefficient()*right.coefficient();
	return -Z*sum;
}

double GTF::ERIGTF(GTF p, GTF q, GTF r, GTF s)
{
	int I,Ip,R,Rp,U;
	int J,Jp,S,Sp,N;
	int K,Kp,T,Tp,W;
	double Sx,Sy,Sz;
	vector<vector<double>> Te;
	Te.resize(2, vector<double>(3));
	double temp1,temp2;
	vector<double> PA(3);
	vector<double> PB(3);

	vector<double> QC(3);
	vector<double> QD(3);
	vector<double> PQ(3);

	
	double g1=p.exposant()+q.exposant();
	double g2=r.exposant()+s.exposant();
	double d=(1.0/g1+1.0/g2)/4;
	double RAB2=0.0;
	double RCD2=0.0;
	double RPQ2=0.0;
	int j;
	
	double sum;

	for(j=0; j<3; j++)
	{
		temp1=(p.exposant()*p.coord()[j]+q.exposant()*q.coord()[j])/g1;
		PA[j]=p.coord()[j]-temp1;
		PB[j]=q.coord()[j]-temp1;
		
		temp2=(r.exposant()*r.coord()[j]+s.exposant()*s.coord()[j])/g2;
		QC[j]=r.coord()[j]-temp2;
		QD[j]=s.coord()[j]-temp2;

		PQ[j]=temp2-temp1;

		RAB2+=(p.coord()[j]-q.coord()[j])*(p.coord()[j]-q.coord()[j]);
		RCD2+=(r.coord()[j]-s.coord()[j])*(r.coord()[j]-s.coord()[j]);
		RPQ2+=PQ[j]*PQ[j];
	}

	sum=0.0;
	for(I=0; I<=p.l()[0]+q.l()[0]; I++)
		for(R=0; R<=I/2; R++)
		{
			Te[0][0]=Theta(I,R,p.l()[0],q.l()[0],PA[0],PB[0],g1);
			for(Ip=0; Ip<=r.l()[0]+s.l()[0]; Ip++)
				for(Rp=0; Rp<=Ip/2; Rp++)
				{
					Te[1][0]=Theta(Ip,Rp,r.l()[0],s.l()[0],QC[0],QD[0],g2);
					for(U=0; U<=(I+2*Ip)/2-R-Rp; U++)
					{
						Sx=B(I,Ip,R,Rp,U,PQ[0],d,Te[0][0],Te[1][0]);
						for(J=0; J<=p.l()[1]+q.l()[1]; J++)
							for(S=0; S<=J/2; S++)
							{
								Te[0][1]=Theta(J,S,p.l()[1],q.l()[1],PA[1],PB[1],g1);
								for(Jp=0; Jp<=r.l()[1]+s.l()[1]; Jp++)
									for(Sp=0; Sp<=Jp/2; Sp++)
									{
										Te[1][1]=Theta(Jp,Sp,r.l()[1],s.l()[1],QC[1],QD[1],g2);
										for(N=0; N<=(J+Jp)/2-S-Sp; N++)
										{
											Sy=B(J,Jp,S,Sp,N,PQ[1],d,Te[0][1],Te[1][1]);
											for(K=0; K<=p.l()[2]+q.l()[2]; K++)
												for(T=0; T<=K/2; T++)
												{
													Te[0][2]=Theta(K,T,p.l()[2],q.l()[2],PA[2],PB[2],g1);
													for(Kp=0; Kp<=r.l()[2]+s.l()[2]; Kp++)
														for(Tp=0; Tp<=Kp/2; Tp++)
														{
															Te[1][2]=Theta(Kp,Tp,r.l()[2],s.l()[2],QC[2],QD[2],g2);
															for(W=0; W<=(K+Kp)/2-T-Tp; W++)
															{
																Sz=B(K,Kp,T,Tp,W,PQ[2],d,Te[0][2],Te[1][2]);
																sum +=Sx*Sy*Sz*F(I+Ip+J+Jp+K+Kp-2*(R+Rp+S+Sp+T+Tp)-U-N-W,RPQ2/4/d);
															}
														}

												}
										}
									}
							}
					}
				}
		}

	sum*=2*PI*PI*sqrt(M_PI)/g1/g2/sqrt(g1+g2)*exp(-p.exposant()*q.exposant()*RAB2/g1)*exp(-r.exposant()*s.exposant()*RCD2/g2)*
			p.coefficient()*q.coefficient()*r.coefficient()*s.coefficient();
	return sum;
}















	//Recopier et remplacer -> par des . (changer aussi les noms des valeurs pointés par les miennes)
	// Revoir GTF*GTF (faire pareil que pour overlap)
	// Rajouter des fonctions (boucle) pour plusieurs star (si besoin)
	// Puis faire ERIGTF, ionicPotentialGTF, kineticGTF et GTFxyzGTF
	// OK mais il faut coder /Theta(), A(), B(), F() (et FTable ??) dans MathFunction.h
