#include<iostream>
#include<analytic/GTF.h>

using namespace std;

GTF::GTF(const double& exposant, const double& coefficient, const vector<double>& coord, const vector<int>& l, 
	Binomial& B) : _exposant(exposant), _coefficient(coefficient), _coord(coord), _l(l), _bino(B){}

double GTF::GTFstarGTF (GTF& left, GTF& right)
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

	c = (M_PI/gama)*sqrt(M_PI/gama)*exp(-left.exposant()*right.exposant()/gama*R2);

	for(j=0; j<3; j++)
	{
		sum[j]=0.0;
		for(i=0; i<=(left.l()[j]+right.l()[j]/2); i++)
			sum[j]+=f(2*i, left.l()[j], right.l()[j], PA[j], PB[j], _bino)*_bino.fact().double_factorial(2*i+1)/(power(gama,i));
	}
	return c*sum[0]*sum[1]*sum[2];
}

double GTF::GTFstarGTFstarGTF (GTF& left, GTF& mid, GTF& right)
{
	vector<double> sum(3);
	double t;
	vector<double> PA(3);
	vector<double> PB(3);
	vector<double> P(3);
	vector<double> QP(3);
	vector<double> QC(3);
	double gama1=left.exposant()+right.exposant();
	double gama=gama1+mid.exposant();
	double R2AB=0.0;
	double R2PC=0.0;
	double c = 0;
	int iAB;
	int i,j;

	for(j=0;j<3;j++)
	{
		t=(left.exposant()*left.coord()[j]+right.exposant()*right.coord()[j])/gama1;
		P[j]=t;
		PA[j]=left.coord()[j]-t;
		PB[j]=right.coord()[j]-t;
		R2AB += (left.coord()[j]-right.coord()[j])*(left.coord()[j]-right.coord()[j]);
	}
	for(j=0;j<3;j++)
	{
		t=(gama1*P[j]+mid.exposant()*mid.coord()[j])/gama;
		QP[j]=P[j]-t;
		QC[j]=mid.coord()[j]-t;
		R2PC += (P[j]-mid.coord()[j])*(P[j]-mid.coord()[j]);
	}
	c = (M_PI/gama)*sqrt(M_PI/gama)*exp(-left.exposant()*right.exposant()/gama1*R2AB)*exp(-gama1*mid.exposant()/gama*R2PC);

	for(j=0;j<3;j++)
	{
		sum[j]=0.0;
		for(iAB=0;iAB<=(left.l()[j]+right.l()[j]);iAB++)
		{
			double fiAB = f(iAB,left.l()[j],right.l()[j],PA[j],PB[j], _bino);
			for(i=0;i<=(iAB+mid.l()[j])/2;i++)
			{
				sum[j] +=
				fiAB*
				f(2*i,iAB,mid.l()[j],QP[j],QC[j], _bino)*
				_bino.fact().double_factorial(2*i-1)/(power(2.0,i)*power(gama,i));
	 		}
		}
	}
	return  c*sum[0]*sum[1]*sum[2];
}

double GTF::GTFstarGTFstarGTFstarGTF(GTF& A,GTF& B, GTF& C, GTF& D)
{
	vector<double> sum(3);
	double t;
	vector<double> PA(3);
	vector<double> PB(3);
	vector<double> QC(3);
	vector<double> QD(3);
	vector<double> P(3);
	vector<double> Q(3);
	vector<double> GP(3);
	vector<double> GQ(3);
	double gama1=A.exposant()+B.exposant();
	double gama2=C.exposant()+D.exposant();
	double gama=gama1+gama2;
	double R2AB=0.0;
	double R2CD=0.0;
	double R2PQ=0.0;
	double c = 0;
	int iAB;
	int iCD;
	int i,j;

	for(j=0;j<3;j++)
	{
		t=(A.exposant()*A.coord()[j]+B.exposant()*B.coord()[j])/gama1;
		P[j]=t;
		PA[j]=A.coord()[j]-t;
		PB[j]=B.coord()[j]-t;
		R2AB += (A.coord()[j]-B.coord()[j])*(A.coord()[j]-B.coord()[j]);
	}
	for(j=0;j<3;j++)
	{
		t=(C.exposant()*C.coord()[j]+D.exposant()*D.coord()[j])/gama2;
		Q[j]=t;
		QC[j]=C.coord()[j]-t;
		QD[j]=D.coord()[j]-t;
		R2CD += (C.coord()[j]-D.coord()[j])*(C.coord()[j]-D.coord()[j]);
	}
	for(j=0;j<3;j++)
	{
		t=(gama1*P[j]+gama2*Q[j])/gama;
		GP[j]=P[j]-t;
		GQ[j]=Q[j]-t;
		R2PQ += (P[j]-Q[j])*(P[j]-Q[j]);
	}
	c = (M_PI/gama)*sqrt(M_PI/gama)
		*exp(-A.exposant()*B.exposant()/gama1*R2AB)
		*exp(-C.exposant()*D.exposant()/gama2*R2CD)
		*exp(-gama1*gama2/gama*R2PQ);


	for(j=0;j<3;j++)
	{
		sum[j]=0.0;
		for(iAB=0;iAB<=(A.l()[j]+B.l()[j]);iAB++)
		{
			double fiAB = f(iAB,A.l()[j],B.l()[j],PA[j],PB[j],_bino);
			for(iCD=0;iCD<=(C.l()[j]+D.l()[j]);iCD++)
			{
				double fiCD = f(iCD,C.l()[j],D.l()[j],QC[j],QD[j],_bino);
				for(i=0;i<=(iAB+iCD)/2;i++)
				{
					sum[j] +=
					fiAB*
					fiCD*
					f(2*i,iAB,iCD,GP[j],GQ[j],_bino)*
					_bino.fact().double_factorial(2*i-1)/(power(2.0,i)*power(gama,i));
	 			}
			}
		}
	}
	return  c*sum[0]*sum[1]*sum[2];
}

double GTF::normeGTF()
{
	return sqrt(2*_exposant/M_PI*sqrt(2*_exposant/M_PI)*power(4*_exposant, _l[0]+_l[1]+_l[2])
			/(_bino.fact().double_factorial(_l[0])*_bino.fact().double_factorial(_l[1])*_bino.fact().double_factorial(_l[2])));
}

double GTF::normeGTF(GTF& q)
{
	return sqrt(2*q.exposant()/M_PI*sqrt(2*q.exposant()/M_PI)*power(4*q.exposant(), q.l()[0]+q.l()[1]+q.l()[2])
			/(_bino.fact().double_factorial(q.l()[0])*_bino.fact().double_factorial(q.l()[1])*_bino.fact().double_factorial(q.l()[2])));
}

void GTF::normalise_radialGTF()
{
	GTF q(_exposant, _coefficient, _coord, _l, _bino);
	int l_bis=_l[0]+_l[1]+_l[2];
	q.l()[0]=l_bis;
	q.l()[1]=0;
	q.l()[2]=0;
	_coefficient*=normeGTF(q);
}

void GTF::normaliseGTF()
{
	_coefficient*=normeGTF();
}

double GTF::overlapGTF(GTF& left, GTF& right)
{
	return left.coefficient()*right.coefficient()*GTFstarGTF(left, right);
}

double GTF::overlap3GTF(GTF& left, GTF& mid, GTF& right)
{
	return left.coefficient()*mid.coefficient()*right.coefficient()*GTFstarGTFstarGTF(left, mid, right);
}

double GTF::overlap4GTF(GTF& A, GTF& B, GTF& C, GTF& D)
{
	return A.coefficient()*B.coefficient()*C.coefficient()*D.coefficient()*GTFstarGTFstarGTFstarGTF(A,B,C,D);
}

double GTF::GTFxyzGTF(GTF& p, GTF& q, int ix, int iy, int iz)
{
	vector<double> C(3,0);
	vector<int> l {ix, iy, iz};
	GTF m(0.0, 1.0, C, l, _bino);
	return overlap3GTF(p, m, q);
}

double GTF::kineticGTF(GTF& left, GTF& right)
{
	int j;
	GTF a=left;
	GTF b=right;
	vector<double> Ti(7);
	vector<double> sum(3);
	double T=0.0;

	for(j=0; j<7; j++)
		Ti[j]=0.0;
	Ti[0]=GTFstarGTF(a, b);

	b.l()[0]=right.l()[0]+2;
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2];
	Ti[1]=GTFstarGTF(a, b);

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1]+2;
	b.l()[2]=right.l()[2];
	Ti[2]=GTFstarGTF(a, b);

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2]+2;
	Ti[3]=GTFstarGTF(a, b);

	b.l()[0]=right.l()[0]-2;
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2];
	Ti[4]=GTFstarGTF(a, b);

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1]-2;
	b.l()[2]=right.l()[2];
	Ti[5]=GTFstarGTF(a, b);

	b.l()[0]=right.l()[0];
	b.l()[1]=right.l()[1];
	b.l()[2]=right.l()[2]-2;
	Ti[6]=GTFstarGTF(a, b);

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

double ionicPotentialGTF(GTF& left, GTF& right, vector<double> C, double Z)
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
	vector<double> FTable;

	for(j=0; j<3; j++)
	{
		temp=(left.exposant()*left.coord()[j]+right.exposant()*right.coord()[j])/gama;
		PA[j]=left.coord()[j]-temp;
		PB[j]=right.coord()[j]-temp;
		PC[j]=-C[j]+temp;
		R2+=(left.coord()[j]-right.coord()[j])*(left.coord()[j]-right.coord()[j]);
		PC2+=PC[j]*PC[j];
	}
	FTable=getFTable(left.l()[0]+right.l()[0]+left.l()[1]+right.l()[1]+left.l()[2]+right.l()[2], gama*PC2);

	sum=0.0;
	for(i=0; i<=left.l()[0]+right.l()[0]; i++)
		for(r=0; r<=i/2; r++)
			for(u=0; u<=(i-2*r)/2; u++)
			{
				Sx=A(i,r,u,left.l()[0],right.l()[0],PA[0],PB[0],PC[0],gama,left.bino());
				for(j=0; j<=left.l()[1]+right.l()[1]; j++)
					for(s=0; s<=j/2; s++)
						for(n=0; n<=(j-2*s)/2;n++)
						{
							Sy=A(j,s,n,left.l()[1],right.l()[1],PA[1],PB[1],PC[1],gama,left.bino());
							for(k=0; k<=left.l()[2]+right.l()[2]; k++)
								for(t=0; t<=k/2; t++)
									for(w=0; w<=(k-2*t)/2; w++)
									{
										Sz=A(k,t,w,left.l()[2],right.l()[2],PA[2],PB[2],PC[2],gama,left.bino());
										sum+=Sx*Sy*Sz*FTable[i+j+k-2*(r+s+t)-u-n-w];
									}
						}
			}
	sum *=2*M_PI/gama*exp(-left.exposant()*right.exposant()/gama*R2)*left.coefficient()*right.coefficient();
	return -Z*sum;
}

double GTF::ERIGTF(GTF& p, GTF& q, GTF& r, GTF& s)
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
			Te[0][0]=Theta(I,R,p.l()[0],q.l()[0],PA[0],PB[0],g1,p.bino());
			for(Ip=0; Ip<=r.l()[0]+s.l()[0]; Ip++)
				for(Rp=0; Rp<=Ip/2; Rp++)
				{
					Te[1][0]=Theta(Ip,Rp,r.l()[0],s.l()[0],QC[0],QD[0],g2,p.bino());
					for(U=0; U<=(I+2*Ip)/2-R-Rp; U++)
					{
						Sx=B(I,Ip,R,Rp,U,PQ[0],d,Te[0][0],Te[1][0],p.bino());
						for(J=0; J<=p.l()[1]+q.l()[1]; J++)
							for(S=0; S<=J/2; S++)
							{
								Te[0][1]=Theta(J,S,p.l()[1],q.l()[1],PA[1],PB[1],g1,p.bino());
								for(Jp=0; Jp<=r.l()[1]+s.l()[1]; Jp++)
									for(Sp=0; Sp<=Jp/2; Sp++)
									{
										Te[1][1]=Theta(Jp,Sp,r.l()[1],s.l()[1],QC[1],QD[1],g2,p.bino());
										for(N=0; N<=(J+Jp)/2-S-Sp; N++)
										{
											Sy=B(J,Jp,S,Sp,N,PQ[1],d,Te[0][1],Te[1][1],p.bino());
											for(K=0; K<=p.l()[2]+q.l()[2]; K++)
												for(T=0; T<=K/2; T++)
												{
													Te[0][2]=Theta(K,T,p.l()[2],q.l()[2],PA[2],PB[2],g1,p.bino());
													for(Kp=0; Kp<=r.l()[2]+s.l()[2]; Kp++)
														for(Tp=0; Tp<=Kp/2; Tp++)
														{
															Te[1][2]=Theta(Kp,Tp,r.l()[2],s.l()[2],QC[2],QD[2],g2,p.bino());
															for(W=0; W<=(K+Kp)/2-T-Tp; W++)
															{
																Sz=B(K,Kp,T,Tp,W,PQ[2],d,Te[0][2],Te[1][2],p.bino());
																sum +=Sx*Sy*Sz*F(I+Ip+J+Jp+K+Kp-2*(R+Rp+S+Sp+T+Tp)-U-N-W,RPQ2/4/d,p.bino().fact());
															}
														}

												}
										}
									}
							}
					}
				}
		}

	sum*=2*M_PI*M_PI*sqrt(M_PI)/g1/g2/sqrt(g1+g2)*exp(-p.exposant()*q.exposant()*RAB2/g1)*exp(-r.exposant()*s.exposant()*RCD2/g2)*
			p.coefficient()*q.coefficient()*r.coefficient()*s.coefficient();
	return sum;
}
