//The order ist nr3.h -> function.h -> Difeq.h -> Solvde.h -> Solving.h -> main.cpp
//The files are included with respect to this order
#include "nr3.h"		
#include "Solving.h"   

ofstream out1("zdata1.dat");     //these are for outputs
ofstream out2("zdata2.dat");     //these are for outputs
ofstream out3("zdata3.dat");     //these are for outputs
ofstream out4("zdata4.dat");     //these are for outputs

int main(){

    complex<double> J(0,0);         //this is for the magnetic domain wall, the absolut value of J. complex only for convenience.
    complex<double> delta(1,0);     //this is the pair potential for the superconductor
    Doub phi = 0.25 * M_PI;         //phase difference between the superconducting elements
    Doub spinup = 1;                //these are for the spin resolution, specifically the up part.
    Doub spindown = 1;              //these are for the spin resolution, specifically the down part.
    Doub ds=5.0;                    //length of the superconducting part of the System
    Doub df=2.0;                    //length of the ferromagnetic(normal) part of the System
    Doub dw=0.1;                    //width of the magnetic domain wall
    Int M=600,MM=1,MMM=100;         //z-grid points, iterations after convergence, E-grid points
    Int itmax = 4000;               //maximal number of iterations
    Doub conv =  10e-8;             //error for convergence
    Doub slowc = 0.5;               //stepsize for convergence 
    Doub epsilon = 10e-5;           //small imaganary part for numerical stability
    MatDoub N(2*MMM+1,M+1);         //Density of states Ntot(E,z)
    Doub sigmaS = 1;                //conducting parameter in superconducting element
    Doub sigmaF = 1;                //conducting parameter in ferromagnetic element
    Doub EL = 0.01;                 //defines which energy interval were looking at: (EM - EL , EM + EL)
    Doub EM = -0.12;                //defines which energy interval were looking at: (EM - EL , EM + EL)
    


    /*    //This yields an image of the 3D-Density of states
    VecDoub x(M+1);                 //The z-grid
    Doub h=(2*ds+df)/M;             //The distance of two neighbored points on the z-grid
    for(int k=0; k<=M; k++){
		x[k] = h *k ;   
    }
    Solving solving(J,delta,ds,df,dw,M,MM,MMM,itmax,conv,slowc,N,h,x,spinup,spindown,phi,epsilon,sigmaS,sigmaF,EL,EM);
    for(int a=-MMM;a<=MMM;a++){
        complex<double> E = (EL/MMM) * a + EM + epsilon * i;
        for(int k=0;k<=M;k++){
            out1 << x[k] << " " << real(E) << " " << N[a+MMM][k] << endl;
        }
        out2 << real(E) << " " << N[a+MMM][M/2] << " " << 3 <<  endl;   
    } 
    */
    

    
    /* //Energy gap with respect to some parameter, width, height and plot
    Int n3=40;                     //number of different paramters we compare in one execution of the program
    Doub deltaN = 0.05;            //acceptance for minigap. If the density of states is at any point smaller then this, it has a minigap.
    VecDoub Nmid(2*MMM+1);         
    VecDoub DF(n3+1);
    VecDoub DW(n3+1);
    VecDoub DS(n3+1);   
    VecDoub PHI(n3+1);    
    for(int k=0;k<=n3;k++){
        DF[k] = 2.0;
    }
    for(int k=0;k<=n3;k++){
        //DW[k] = 0.5 * k;
        DW[k] = 0.05 * k;
    }
    for(int k=0;k<=n3;k++){
        //DS[k] = 0.5 * k;
        DS[k] = 5.0;
    }
    for(int k=0;k<=n3;k++){
        PHI[k] = 0;
        //PHI[k] = 0;
    }
    VecDoub x(M+1);
    Doub h;
    Doub Nmidmax ;
    Doub Emidmax ;
    Doub Nmidmin ;
    Doub Emidmin ;
    Doub Nmax;
    Doub widthEmin,widthEmax;
    Doub widthE;
    for(int j=0;j<=n3;j++){
        df = DF[j];
        ds = DS[j];
        dw = DW[j];
        phi = PHI[j];
        h=(2*ds+df)/M;
        for(int k=0; k<=M; k++){
		    x[k] = h * k ;
        }
        Solving solving(J,delta,ds,df,dw,M,MM,MMM,itmax,conv,slowc,N,h,x,spinup,spindown,phi,epsilon);
        for(int a=-MMM;a<=MMM;a++){
            complex<double> E = (1.5/MMM) * a + epsilon * i;
            out1 << real(E) << " " << N[a+MMM][M/2] << " " << j << endl;        //plotting the density of states at z=M/2 (In the middle of the domain)
        }
        Nmidmax = 0;
        Emidmax = 0;
        Nmidmin = 0;
        Emidmin = 0;
        for(int a=-MMM;a<MMM;a++){
            complex<double> E = (1.5/MMM) * a + epsilon * i;
            Nmid[a+MMM] = N[a+1+MMM][M/2] - N[a+MMM][M/2];
            if(Nmid[a+MMM]>Nmidmax){
                Nmidmax = Nmid[a+MMM];
                Emidmax = real(E);
            }
            if(Nmid[a+MMM]<Nmidmin){
                Nmidmin = Nmid[a+MMM];
                Emidmin = real(E);
            }
        }                                                                       //computing the maximal/minimal slope yiels the edges of the minigap (!) ->
        out2 << abs(Emidmax-Emidmin) << " " << PHI[j] << endl;                  //width of the minigaps (1.method) (doesn't decide if there is a minigap)
        for(int a=-MMM;a<MMM;a++){
            complex<double> E = (1.5/MMM) * a + epsilon * i;
            if(N[a+MMM][M/2] > deltaN && N[a+1+MMM][M/2] <=deltaN){
                widthEmin = real(E);
            }
            if(N[a+MMM][M/2] <= deltaN && N[a+1+MMM][M/2] > deltaN){
                widthEmax = real(E);
            }    
            widthE = widthEmax - widthEmin;
        }                                                                       //deciding whether a point is under the threshold for a minigap or not
        out3 << abs(widthE) << " " << PHI[j] << endl;                           //width of the minigaps (2.method) (decides if there is a minigap)
        Nmax = 0;
        for(int a=-MMM;a<MMM;a++){
            if(N[a+MMM][M/2] > Nmax){
                Nmax = N[a+MMM][M/2];
            }
        }
        //out4 << N[MMM][M/2] << " " << DW[j] << endl;                          //maximal height of the density of states at ds+df/2
    }
    */
    


    /* Illustration of the magentic structure (on the whole domain). Works independent from the parameters above!
     //plotting Jx and Jz
    complex<double> J(0.5,0);   //for the magentic domain wall, J is the effective exchange field.
	Doub ds=5,df=2;			    //as described earlier, these are the region lengths		
    const Int M=1200;		    //number of grid points
	VecDoub x(M+1);				//grid vector
	Doub h=(2*ds+df)/M;			//step size 
	int n1,n2;                  //inner edges
    
    for(int k=0;k<=M;k++){
		x[k] = h *k ;						
	}
    for(int k=0;k<=M;k++){
        if(x[k] <= ds && x[k+1] >ds){
            n1 = k+1;
        }
        if(x[k] <= ds+df && x[k+1] >ds+df){
            n2 = k+1;
        }
    }
    VecDoub DW(7);
    DW[0]=0, DW[1] =0.01, DW[2] =0.1, DW[3] =0.25, DW[4] =0.5, DW[5] =1.0, DW[6] =2.0;
    for(int j=1;j<=6;j++){
        for(int k=n1;k<n2-1;k++){
            Doub dw = DW[j];
            out1 << x[k] << " " << real(Jx(x[k],ds,df,dw,J)/J) << " " << j <<  endl;        //the x-component of J (not to mistake with the grid x)
            out2 << x[k] << " " << real(Jz(x[k],ds,df,dw,J)/J) << " " << j <<  endl;        //the z-component of J (not to mistake with the spatial coordinate z)
        }
    }
    */    

}
