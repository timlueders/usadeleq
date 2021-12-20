#include "nr3.h"			//As per usual
#include "Solvde.h"			//Solvde.h contains the include of Difeq.h
# define M_PI           3.14159265358979323846  //Defining Pi, as we will need it

struct Solving {		//using structures simplifies solving the equation for many sets of parameters and saves time  
    complex<double> J;				
    complex<double> delta;				
    Doub ds,df,dw;
    Int M,MM,MMM;		
    Int itmax;			
    Doub conv;			
    Doub slowc;	
	MatDoub &N;
	Doub h;
	VecDoub x;
	Doub spinup;		
	Doub spindown;
	Doub phi;
	Doub epsilon;
	Doub sigmaS,sigmaF;
	Doub EL;
	Doub EM;
    //Solvde &solvde;
    Solving(complex<double> JJ, complex<double> deltaa, Doub dss,  Doub dff ,  Doub dww,
    Int Mm, Int MMm, Int MMMm, Int itmaxx, Doub convv, Doub slowcc, MatDoub_IO &NN, Doub hh, VecDoub xx,
	Doub spinupp, Doub spindownn, Doub phii, Doub epsilonn, Doub sigmaSS, Doub sigmaFF, Doub ELL, Doub EMM);
};

Solving::Solving(complex<double> JJ, complex<double> deltaa, Doub dss,  Doub dff ,  Doub dww,
    Int Mm, Int MMm, Int MMMm, Int itmaxx, Doub convv, Doub slowcc, MatDoub_IO &NN, Doub hh, VecDoub xx,
	Doub spinupp, Doub spindownn, Doub phii, Doub epsilonn, Doub sigmaSS, Doub sigmaFF, Doub ELL, Doub EMM) : 
	J(JJ), delta(deltaa), ds(dss), df(dff), dw(dww), M(Mm),MM(MMm),MMM(MMMm), itmax(itmaxx), 
	conv(convv),slowc(slowcc), N(NN), h(hh), x(xx), spinup(spinupp), spindown(spindownn), phi(phii), epsilon(epsilonn),
	sigmaS(sigmaSS), sigmaF(sigmaFF), EL(ELL), EM(EMM)
{	
    double N_0 = 1, Ntot = 0;				//Ntot is the density of states, N_0 is just a scaling factor for N.
	complex<double> Gupup;					//These are defined in the thesis, eq.(1),(35)
	complex<double> Gdowndown;				//These are defined in the thesis, eq.(1),(35)
	const Int NE=16,NB=8,NYJ=NE,NYK=M+1;	//NE is the number of equations, NB is the number of edgeconditions at the first edge
	Int mpt= M+1;							//Is often needed, which is why it is given an name aswell
	VecInt indexv(NE);						//This vector is there for ordering the columns of the matrix s in Difeq.
	VecDoub scalv(NE);						//scalv measures the typical size of each variable
	MatComplex y(NYJ,NYK);					//y is the matrix we search, the first entry is the equation, the second the grid point
	Int n1,n2;								//this will be the inner edge points
	Doub scalvestimation;					//for an approximate size of the variables or scalv, respectively.
	indexv[0]=8;							//the entries of indexv come from the outer edges in Difeq.h
	indexv[1]=9;						
	indexv[2]=10;						
	indexv[3]=11;						
	indexv[4]=0;							
	indexv[5]=1;						 	
	indexv[6]=2;						
	indexv[7]=3;						
	indexv[8]=12;						
	indexv[9]=13;					
	indexv[10]=14;							
	indexv[11]=15;							
	indexv[12]=4;							
	indexv[13]=5;						 	
	indexv[14]=6;							
	indexv[15]=7;							
	for(int k=0;k<=M;k++){
		for(int q8 = 0; q8<=15;q8++){
			y[q8][k] = 0;					//initial guess, setting everything to zero works usually fine
		}
	}
	for(int q7=0; q7<=15; q7++){
		scalv[q7] = 1.0;		//initial size of the variables, setting everything to one works usually fine
	}
	for (Int k=0;k<M;k++) {		//to determine the inner edges
		if( x[k]<=ds && x[k+1] > ds){
			n1 = k+1;		//n1 is the first inner edge
		}
		if( x[k]<= ds+df && x[k+1] > ds +df){
			n2 = k+1;		//n2 is the second inner edge
		}
	}
	for(int a=-MMM; a<=MMM ; a++){
		complex<double> E = (EL/MMM) * a + EM + epsilon * i;		//Here we start the loop over the energy
		for(int q=0;q<=15;q++){
			scalvestimation = 0;
			for(int k=0;k<=M;k++){
				scalvestimation += abs(y[q][k]);
			}
			scalvestimation /= M;
			scalv[q] = scalvestimation + 0.001;			//this ensures senseful values for scalv, as it measures the size of the y[][] in a certain scope.
		}
		for (Int j=0;j<MM;j++) {		//solving the differential equation with Solvde.h and Difeq.h
			Difeq difeq(mpt,n1,n2,h,ds,df,E,delta,x,dw,J,phi,sigmaS,sigmaF); 		
			Solvde solvde(itmax,conv,slowc,scalv,indexv,NB,y,difeq);
		}
		for(int k=0;k<=M;k++){			//Calculating the Density of states.
			complex <double> zeta1 = y[0][k]*y[8][k] +  y[1][k]*y[10][k] ;
			complex <double> zeta2 = y[0][k]*y[9][k] +  y[1][k]*y[11][k] ;
			complex <double> zeta3 = y[2][k]*y[8][k] +  y[3][k]*y[10][k] ;
			complex <double> zeta4 = y[2][k]*y[9][k] +  y[3][k]*y[11][k] ;
			complex <double> zeta5 = -M_PI * i * pow((one-zeta1)*(one-zeta4)-zeta2*zeta3 ,-1);
			Gupup = spinup * ((one-zeta4)*(one+zeta1)+(zeta2*zeta3)) * zeta5 ;
			Gdowndown = spindown * ((one-zeta1)*(one+zeta4)+(zeta2*zeta3)) * zeta5 ;
			Ntot = -N_0 * pow(2*M_PI,-1) * imag(Gupup + Gdowndown);			//the total density of states
			N[a+MMM][k] = Ntot;												//defines a matrix that contains the density of states on the energy and spatial domain
		}																	//With this matrix we can compute the width of the gap etc. 

		//cout << real(E) << endl;		//If there are convergence problems etc. for certain energies, this is the way to find them.
	}
}