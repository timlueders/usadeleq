#include "nr3.h"		//As mentioned, this file is needed and provided by numerical recipes

//ofstream out("test.dat");  //this is for outputs, especially to check definitions here or in Difeq

complex<double>i(0,1);     //imaginary unit
complex<double>one(1,0);   //complex one, appears to be useful sometimes

//These are some generic operations to avoid getting problems by adding integers to complex numbers for example.
complex<double> dtc(double x, complex<double> y){      //multiplication of a double with a complex number, order is important!
    complex<double> z(x,0);
    return z*y ;
}complex<double> dpc(double x, complex<double> y){     //addition of a double with a complex number, order is important!
    complex<double> z(x,0);
    return z+y ;
}
complex<double> itc(int x, complex<double> y){         //multiplication of a integer with a complex number, order is important!
    complex<double> z(x,0);
    return z*y ;
}complex<double> ipc(int x, complex<double> y){        //addition of a integer with a complex number, order is important!
    complex<double> z(x,0);
    return z+y ;
}

//The pair potential with phase difference (as described in the thesis eq.(34))
//z is the spatial coordinate, ds is the length of each superconducting elements, df is the lengt of the ferromagnetic(normal) metall,
//delta is the gap parameter, phi the phase difference
complex<double> Delta(double z, double ds, double df, complex<double> delta, Doub phi){ 
    if(z<=ds){
        return delta * exp(i*phi);       //implement a phase difference via exp(i*phi) at one superconductor
    } 
    else if(z >= ds+df){
        return delta;                    //the other superconductor has only the delta part.
    }
    else{                                //for BCS-Systems we have Delta=delta on the whole domain.
        return 0;                        //for the SNS/SFS system we have a zero here.
    }
}

//Definition of the magnetic structure (as described in the thesis eq.(36))
//dw is the length of the magnetic domain wall
complex<double> Theta(double z,double ds,double df,double dw){
    double z_0 = (2*ds+df) * 0.5;
    return -atan((z-z_0) * pow(dw,-1));
}   
complex<double> Jx(double z, double ds, double df, double dw, complex<double> J){       //Jx for the SFS-Systems
    if(z > ds && z < ds + df){
        return J * cos(Theta(z,ds,df,dw));             //magnetic domain wall
        //return zero;                                 //fully polarized or in the SNS-System
    }
    else{
        return 0;                                      //J neq 0 only in the F part of SFS, as per definition
    }
    
}
complex<double> Jy(double z, double ds, double df, double dw, complex<double> J){       //Jy for the SFS-Systems                       
    return 0;                                          //always zero; SNS, fully polarized and magnetic domain wall
}
complex<double> Jz(double z, double ds, double df, double dw, complex<double> J){       //Jz for the SFS-Systems               
    if(z > ds && z < ds + df){
        return J * sin(Theta(z,ds,df,dw));             //magnetic domain wall
        //return J;                                    //fully polarized or in the SNS-System(J=0)
    }
    else{
        return 0;                                      //J neq 0 only in the F part of SFS, as per definition
    }      
}   

//matrix operations, only for the present case n=2
MatComplex Matadd(MatComplex A, MatComplex B, MatComplex &C){      //Addition
    for(int j=0; j<=1 ;j++){
        for(int k=0; k<=1; k++){
            C[j][k] = A[j][k] + B[j][k];
        }
    }
    return C;
}
MatComplex Skalmult(Complex Z, MatComplex B, MatComplex &C){      //2x2-Scalarmultiplication
    for(int j=0; j<=1 ;j++){
        for(int k=0; k<=1; k++){
            C[j][k] = Z * B[j][k];
        }
    }
    return C;
}
MatComplex Matmult(MatComplex A, MatComplex B, MatComplex &C){      //2x2-Matrixmultiplication
    for(int m=0; m<=1; m++){
       for(int n=0; n<=1; n++){
            C[m][n] = 0;
        }
    }
    for(int j=0; j<=1 ;j++){
        for(int k=0; k<=1; k++){
            for(int n=0; n<=1; n++){
                C[j][k] += A[j][n]*B[n][k];
            }
        }
    }
    return C;
}
MatComplex Matinv(MatComplex A, MatComplex &C){                     //Inversion
    if(abs(A[0][0]*A[1][1] - A[0][1]*A[1][0]) != 0){
        C[0][0] = A[1][1];
        C[0][1] = -A[0][1];
        C[1][0] = -A[1][0];
        C[1][1] = A[0][0];
        for(int m=0; m<=1; m++){
        	 for(int n=0; n<=1; n++){
                C[m][n] /= A[0][0]*A[1][1] - A[0][1]*A[1][0];
            }
        }
        return C;
    }
    else{
        cout << "the matrix is not invertible" << endl;
        return A;
    }
}
