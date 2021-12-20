#include "nr3.h"
#include "function.h"


struct Difeq {				
	const Int &mpt, &n1, &n2;
	const Doub &h, &ds, &df;	
	const complex<double> &E, &delta;
	const VecDoub &x;		
	const Doub &dw;
	const complex<double> &J;	
	const Doub &phi;					
	Difeq(const Int &mptt, const Int &n11, const Int &n22, const Doub &hh, const Doub &dss, const Doub &dff, 
		const complex<double> &EE, const complex<double> &deltaa, VecDoub_I &xx, const Doub &dww, 
		const complex<double> &JJ, const Doub &phii) :
		mpt(mptt), n1(n11), n2(n22), h(hh), ds(dss), df(dff), E(EE), delta(deltaa), x(xx),dw(dww),J(JJ),
		phi(phii){}
	void smatrix(const Int k, const Int k1, const Int k2, const Int jsf,
		const Int is1, const Int isf, VecInt_I &indexv, MatComplex_O &s,
		MatComplex_I &y)
	{	
		for(int q2=0; q2 <= 15; q2 ++){		
			for(int q3=0; q3 <= 15; q3 ++){			
				s[q2][indexv[q3]] = 0;				
				s[q2][16+indexv[q3]] = 0;			
			}
		}

		if (k==k1) {						//left edge				
			s[8][16+indexv[4]] = 1.0;
			s[9][16+indexv[5]] = 1.0;
			s[10][16+indexv[6]] = 1.0;
			s[11][16+indexv[7]] = 1.0;
			s[12][16+indexv[12]] = 1.0;
			s[13][16+indexv[13]] = 1.0;
			s[14][16+indexv[14]] = 1.0;
			s[15][16+indexv[15]] = 1.0;
			s[8][jsf] = y[4][0];
			s[9][jsf] = y[5][0];
			s[10][jsf] = y[6][0];
			s[11][jsf] = y[7][0];
			s[12][jsf] = y[12][0];
			s[13][jsf] = y[13][0];
			s[14][jsf] = y[14][0];
			s[15][jsf] = y[15][0];
		}

		else if (k > k2-1) {				//right edge	
			s[0][16+indexv[4]] = 1.0;
			s[1][16+indexv[5]] = 1.0;
			s[2][16+indexv[6]] = 1.0;
			s[3][16+indexv[7]] = 1.0;
			s[4][16+indexv[12]] = 1.0;
			s[5][16+indexv[13]] = 1.0;
			s[6][16+indexv[14]] = 1.0;
			s[7][16+indexv[15]] = 1.0;
			s[0][jsf] = y[4][mpt-1];
			s[1][jsf] = y[5][mpt-1];
			s[2][jsf] = y[6][mpt-1];
			s[3][jsf] = y[7][mpt-1];
			s[4][jsf] = y[12][mpt-1];
			s[5][jsf] = y[13][mpt-1];
			s[6][jsf] = y[14][mpt-1];
			s[7][jsf] = y[15][mpt-1];
		}

		else if ( k == n1 || k == n2){		//inner edges	
			for(int q1=0; q1<= 15 ; q1++){
				s[q1][indexv[q1]]= -1 ;
				s[q1][16+indexv[q1]] = 1;
				s[q1][jsf] = y[q1][k] - y[q1][k-1] ;
			}
		}
		
		else{								//interior
			//first the simple derivatives/equations for gamma and gamma tilde.	
			for(int q4=0; q4 <= 3; q4 ++){	
				s[q4][indexv[q4]] = -1;
				s[q4][indexv[4+q4]] = -0.5*h;
				s[q4][16+indexv[q4]] = 1;
				s[q4][16+indexv[4+q4]] = -0.5*h;
				s[8+q4][indexv[8+q4]] = -1;
				s[8+q4][indexv[8+4+q4]] = -0.5*h;
				s[8+q4][16+indexv[8+q4]] = 1;
				s[8+q4][16+indexv[8+4+q4]] = -0.5*h;
			}
			s[0][jsf] = y[0][k] - y[0][k-1] - 0.5 * h * (y[4][k]+y[4][k-1]);
			s[1][jsf] = y[1][k] - y[1][k-1] - 0.5 * h * (y[5][k]+y[5][k-1]);
			s[2][jsf] = y[2][k] - y[2][k-1] - 0.5 * h * (y[6][k]+y[6][k-1]);
			s[3][jsf] = y[3][k] - y[3][k-1] - 0.5 * h * (y[7][k]+y[7][k-1]);
			s[8][jsf] = y[8][k] - y[8][k-1] - 0.5 * h * (y[12][k]+y[12][k-1]);
			s[9][jsf] = y[9][k] - y[9][k-1] - 0.5 * h * (y[13][k]+y[13][k-1]);
			s[10][jsf] = y[10][k] - y[10][k-1] - 0.5 * h * (y[14][k]+y[14][k-1]);
			s[11][jsf] = y[11][k] - y[11][k-1] - 0.5 * h * (y[15][k]+y[15][k-1]);

			//Now the complicated derivatives/equations for gamma and gamma tilde. For this we use the matrix operations from function.h
			Complex two(2,0);								
			Doub Xk = 0.5 * (x[k] + x[k-1]);
			MatComplex gamma(2,2);							//Setting up all the matrices we will need
			MatComplex gamma_tilde(2,2);
			MatComplex gamma_prime(2,2);
			MatComplex gamma_tilde_prime(2,2);
			MatComplex Delta_mat(2,2);
			MatComplex Delta_conj_mat(2,2);
			MatComplex J_sigma(2,2);
			MatComplex J_sigma_conj(2,2);
			gamma[0][0]= 0.5 * (y[0][k]+y[0][k-1]);			//translating the notation of the book to our gamma matrices
			gamma[0][1]= 0.5 * (y[1][k]+y[1][k-1]);
			gamma[1][0]= 0.5 * (y[2][k]+y[2][k-1]);
			gamma[1][1]= 0.5 * (y[3][k]+y[3][k-1]);
			gamma_prime[0][0]= 0.5 * (y[4][k]+y[4][k-1]);
			gamma_prime[0][1]= 0.5 * (y[5][k]+y[5][k-1]);
			gamma_prime[1][0]= 0.5 * (y[6][k]+y[6][k-1]);
			gamma_prime[1][1]= 0.5 * (y[7][k]+y[7][k-1]);
			gamma_tilde[0][0]= 0.5 * (y[8][k]+y[8][k-1]);
			gamma_tilde[0][1]= 0.5 * (y[9][k]+y[9][k-1]);
			gamma_tilde[1][0]= 0.5 * (y[10][k]+y[10][k-1]);
			gamma_tilde[1][1]= 0.5 * (y[11][k]+y[11][k-1]);
			gamma_tilde_prime[0][0]= 0.5 * (y[12][k]+y[12][k-1]);
			gamma_tilde_prime[0][1]= 0.5 * (y[13][k]+y[13][k-1]);
			gamma_tilde_prime[1][0]= 0.5 * (y[14][k]+y[14][k-1]);
			gamma_tilde_prime[1][1]= 0.5 * (y[15][k]+y[15][k-1]);
			Delta_mat[0][0] = 0;							//defining Delta and J*sigma in the natural way
			Delta_mat[0][1] = Delta(Xk,ds,df,delta,phi);
			Delta_mat[1][0] = -Delta(Xk,ds,df,delta,phi);
			Delta_mat[1][1] = 0;
			Delta_conj_mat[0][0] = 0;
			Delta_conj_mat[0][1] = conj(Delta(Xk,ds,df,delta,phi));
			Delta_conj_mat[1][0] = -conj(Delta(Xk,ds,df,delta,phi));
			Delta_conj_mat[1][1] = 0;
			J_sigma[0][0] = Jz(Xk,ds,df,dw,J);
			J_sigma[0][1] = Jx(Xk,ds,df,dw,J) - i * Jy(Xk,ds,df,dw,J);
			J_sigma[1][0] = Jx(Xk,ds,df,dw,J) + i * Jy(Xk,ds,df,dw,J);
			J_sigma[1][1] = -Jz(Xk,ds,df,dw,J);
			J_sigma_conj[0][0] = Jz(Xk,ds,df,dw,J);
			J_sigma_conj[0][1] = Jx(Xk,ds,df,dw,J) + i * Jy(Xk,ds,df,dw,J);
			J_sigma_conj[1][0] = Jx(Xk,ds,df,dw,J) - i * Jy(Xk,ds,df,dw,J);
			J_sigma_conj[1][1] = -Jz(Xk,ds,df,dw,J);

			MatComplex J_sigma_gamma(2,2);										//calculating J_sigma_gamma 
			Matmult(J_sigma,gamma,J_sigma_gamma);
			
			MatComplex minus_gamma_J_sigma_conj(2,2);							//calculating minus_gamma_J_sigma_conj
			MatComplex gamma_J_sigma_conj(2,2);
			Matmult(gamma,J_sigma_conj,gamma_J_sigma_conj);
			Skalmult(-one,gamma_J_sigma_conj,minus_gamma_J_sigma_conj);

			MatComplex J_sigma_conj_gamma_tilde(2,2);							//calculating J_sigma_conj_gamma_tilde
			Matmult(J_sigma_conj,gamma_tilde,J_sigma_conj_gamma_tilde);

			MatComplex minus_gamma_tilde_J_sigma(2,2);							//calculating minus_gamma_tilde_J_sigma
			MatComplex gamma_tilde_J_sigma(2,2);
			Matmult(gamma_tilde,J_sigma,gamma_tilde_J_sigma);
			Skalmult(-one,gamma_tilde_J_sigma,minus_gamma_tilde_J_sigma);

			MatComplex gamma_Delta_conj_gamma(2,2);			//calculating gamma*conj(Delta)*gamma
			MatComplex gamma_Delta_conj(2,2);
			MatComplex Delta_conj_gamma(2,2);
			Matmult(gamma,Delta_conj_mat,gamma_Delta_conj);
			Matmult(Delta_conj_mat,gamma,Delta_conj_gamma);
			Matmult(gamma_Delta_conj,gamma,gamma_Delta_conj_gamma);
			
			MatComplex gamma_tilde_Delta_gamma_tilde(2,2);	//calculating gammatilde*Delta*gammatilde
			MatComplex gamma_tilde_Delta(2,2);
			MatComplex Delta_gamma_tilde(2,2);
			Matmult(gamma_tilde,Delta_mat,gamma_tilde_Delta);
			Matmult(Delta_mat,gamma_tilde,Delta_gamma_tilde);
			Matmult(gamma_tilde_Delta,gamma_tilde,gamma_tilde_Delta_gamma_tilde);

			MatComplex minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime(2,2); //calculating  minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime
			MatComplex gamma_tilde_gamma(2,2);
			MatComplex one_minus_gamma_tilde_gamma(2,2);
			MatComplex inv_one_minus_gamma_tilde_gamma(2,2);
			MatComplex gamma_prime_inv_one_minus_gamma_tilde_gamma(2,2);
			MatComplex gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde(2,2);
			MatComplex gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime(2,2);
			Matmult(gamma_tilde,gamma,gamma_tilde_gamma);
			one_minus_gamma_tilde_gamma[0][0] = one - gamma_tilde_gamma[0][0]; 
			one_minus_gamma_tilde_gamma[1][1] = one - gamma_tilde_gamma[1][1]; 
			one_minus_gamma_tilde_gamma[0][1] = gamma_tilde_gamma[0][1]; 
			one_minus_gamma_tilde_gamma[1][0] = gamma_tilde_gamma[1][0]; 
			Matinv(one_minus_gamma_tilde_gamma,inv_one_minus_gamma_tilde_gamma);
			Matmult(gamma_prime,inv_one_minus_gamma_tilde_gamma,gamma_prime_inv_one_minus_gamma_tilde_gamma);
			Matmult(gamma_prime_inv_one_minus_gamma_tilde_gamma,gamma_tilde, gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde);
			Matmult(gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde, gamma_prime, gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime);
			Skalmult(-two,gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime,minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime);

			MatComplex minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime(2,2); //calculating  minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime
			MatComplex gamma_gamma_tilde(2,2);
			MatComplex one_minus_gamma_gamma_tilde(2,2);
			MatComplex inv_one_minus_gamma_gamma_tilde(2,2);
			MatComplex gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde(2,2);
			MatComplex gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma(2,2);
			MatComplex gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime(2,2);
			Matmult(gamma,gamma_tilde,gamma_gamma_tilde);
			one_minus_gamma_gamma_tilde[0][0] = one - gamma_gamma_tilde[0][0]; 
			one_minus_gamma_gamma_tilde[1][1] = one - gamma_gamma_tilde[1][1]; 
			one_minus_gamma_gamma_tilde[0][1] = gamma_gamma_tilde[0][1]; 
			one_minus_gamma_gamma_tilde[1][0] = gamma_gamma_tilde[1][0]; 
			Matinv(one_minus_gamma_gamma_tilde,inv_one_minus_gamma_gamma_tilde);
			Matmult(gamma_tilde_prime,inv_one_minus_gamma_gamma_tilde,gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde);
			Matmult(gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde,gamma,gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma);
			Matmult(gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma,gamma_tilde_prime,gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime);
			Skalmult(-two,gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime,minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime);

			compMat3DDoub gamma_deriv(2,2,16);						//calculating all derivatives of gamma
			compMat3DDoub gamma_tilde_deriv(2,2,16);				//calculating all derivatives of gammatilde
			compMat3DDoub gamma_prime_deriv(2,2,16);				//calculating all derivatives of gamma_prime
			compMat3DDoub gamma_tilde_prime_deriv(2,2,16);			//calculating all derivatives of gammatilde_prime
				for(int q1=0;q1<=1;q1++){
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=15;q3++){
							gamma_deriv[q1][q2][q3] = 0;
							gamma_tilde_deriv[q1][q2][q3] = 0;
							gamma_prime_deriv[q1][q2][q3] = 0;
							gamma_tilde_prime_deriv[q1][q2][q3]=0;
						}
					}
				}
				gamma_deriv[0][0][0] = 0.5;
				gamma_deriv[0][1][1] = 0.5;
				gamma_deriv[1][0][2] = 0.5;
				gamma_deriv[1][1][3] = 0.5;
				gamma_tilde_deriv[0][0][8] = 0.5;
				gamma_tilde_deriv[0][1][9] = 0.5;
				gamma_tilde_deriv[1][0][10] = 0.5;
				gamma_tilde_deriv[1][1][11] = 0.5;
				gamma_prime_deriv[0][0][4] = 0.5;
				gamma_prime_deriv[0][1][5] = 0.5;
				gamma_prime_deriv[1][0][6] = 0.5;
				gamma_prime_deriv[1][1][7] = 0.5;
				gamma_tilde_prime_deriv[0][0][12] = 0.5;
				gamma_tilde_prime_deriv[0][1][13] = 0.5;
				gamma_tilde_prime_deriv[1][0][14] = 0.5;
				gamma_tilde_prime_deriv[1][1][15] = 0.5;
			
			compMat3DDoub J_sigma_gamma_deriv(2,2,16);						//calculating all derivatives of J_sigma_gamma
			compMat3DDoub J_sigma_conj_gamma_tilde_deriv(2,2,16);			//calculating all derivatives of J_sigma_conj_gamma_tilde
			compMat3DDoub minus_gamma_J_sigma_conj_deriv(2,2,16);			//calculating all derivatives of minus_gamma_J_sigma_conj
			compMat3DDoub minus_gamma_tilde_J_sigma_deriv(2,2,16);			//calculating all derivatives of minus_gamma_tilde_J_sigma
				MatComplex X1(2,2),X2(2,2),X3(2,2),X4(2,2),X5(2,2),X6(2,2);
				for(int q3=0;q3<=15;q3++){	
					for(int q4=0;q4<=1;q4++){
						for(int q5=0;q5<=1;q5++){
							X1[q4][q5] = gamma_deriv[q4][q5][q3];
							X2[q4][q5] = gamma_tilde_deriv[q4][q5][q3];
						}
					}	
					Matmult(J_sigma,X1,X3);
					Matmult(J_sigma_conj,X2,X4);
					Matmult(X1,J_sigma_conj,X5);
					Matmult(X2,J_sigma,X6);
					for(int q4=0;q4<=1;q4++){
						for(int q5=0;q5<=1;q5++){
							J_sigma_gamma_deriv[q4][q5][q3] = X3[q4][q5];
							J_sigma_conj_gamma_tilde_deriv[q4][q5][q3] = X4[q4][q5];
							minus_gamma_J_sigma_conj_deriv[q4][q5][q3] = - X5[q4][q5];
							minus_gamma_tilde_J_sigma_deriv[q4][q5][q3] = - X6[q4][q5];
						}
					}	
				}

			compMat3DDoub gamma_Delta_conj_gamma_deriv(2,2,16);			//calculating all derivaties of gamma_Delta_conj_gamma
			compMat3DDoub gamma_tilde_Delta_gamma_tilde_deriv(2,2,16); 	//calculating all derivaties of gamma_tilde_Delta_gamma_tilde
				MatComplex A1(2,2),B1(2,2),C1(2,2),A2(2,2),B2(2,2),C2(2,2);
				for(int q3=0;q3<=15;q3++){	
					for(int q4=0;q4<=1;q4++){
						for(int q5=0;q5<=1;q5++){
							A1[q4][q5] = gamma_deriv[q4][q5][q3];
							A2[q4][q5] = gamma_tilde_deriv[q4][q5][q3];
						}
					}
					Matmult(A1,Delta_conj_gamma,B1);
					Matmult(gamma_Delta_conj,A1,C1);
					Matmult(A2,Delta_gamma_tilde,B2);
					Matmult(gamma_tilde_Delta,A2,C2);	
					for(int q4=0;q4<=1;q4++){
						for(int q5=0;q5<=1;q5++){
							gamma_Delta_conj_gamma_deriv[q4][q5][q3] = B1[q4][q5] + C1[q4][q5];
							gamma_tilde_Delta_gamma_tilde_deriv[q4][q5][q3] = B2[q4][q5] + C2[q4][q5];
						}
					}
				}
			
			compMat3DDoub minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime_deriv(2,2,16);	//calculating all derivaties of minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime
			compMat3DDoub minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime_deriv(2,2,16);  //calculating all derivaties of  minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime
				MatComplex inv_one_minus_gamma_tilde_gamma_gamma_tilde(2,2);
				Matmult(inv_one_minus_gamma_tilde_gamma,gamma_tilde,inv_one_minus_gamma_tilde_gamma_gamma_tilde);
				MatComplex X7(2,2),X8(2,2),X9(2,2),X10(2,2),X11(2,2);
				for(int q1=4;q1<=7;q1++){			//derivatives after gamma_prime
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							X7[q2][q3] = gamma_prime_deriv[q2][q3][q1];
						}
					}
					Matmult(X7, inv_one_minus_gamma_tilde_gamma_gamma_tilde,X8);
					Matmult(X8, gamma_prime,X9);
					Matmult(inv_one_minus_gamma_tilde_gamma_gamma_tilde,X7,X10);
					Matmult(gamma_prime,X10,X11);
					for(int q4=0;q4<=1;q4++){
						for(int q5=0;q5<=1;q5++){
							minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime_deriv[q4][q5][q1] = 0;
							minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime_deriv[q4][q5][q1] = -two * (X9[q4][q5] + X11[q4][q5]);
						}
					}	
				}
				MatComplex inv_one_minus_gamma_gamma_tilde_gamma(2,2);
				Matmult(inv_one_minus_gamma_gamma_tilde,gamma,inv_one_minus_gamma_gamma_tilde_gamma);
				MatComplex X12(2,2),X13(2,2),X14(2,2),X15(2,2),X16(2,2);
				for(int q1=12;q1<=15;q1++){			//derivatives after gamma_tilde_prime
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							X12[q2][q3] = gamma_tilde_prime_deriv[q2][q3][q1];
						}
					}
					Matmult(X12,inv_one_minus_gamma_gamma_tilde_gamma ,X13);
					Matmult(X13, gamma_tilde_prime,X14);
					Matmult(inv_one_minus_gamma_gamma_tilde_gamma,X12,X15);
					Matmult(gamma_tilde_prime,X15,X16);
					for(int q4=0;q4<=1;q4++){
						for(int q5=0;q5<=1;q5++){
							minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime_deriv[q4][q5][q1] = 0;
							minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime_deriv[q4][q5][q1] = -two * (X14[q4][q5] + X16[q4][q5]);
						}
					}	
				}
				MatComplex X17(2,2),X18(2,2),X19(2,2),X20(2,2),X21(2,2),X22(2,2),X23(2,2),X24(2,2),
				X25(2,2),X26(2,2),X27(2,2),X28(2,2),X29(2,2),X30(2,2),X31(2,2);
				for(int q1=0;q1<=3;q1++){		//derivatives after gamma	
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							X17[q2][q3] = gamma_deriv[q2][q3][q1];
						}
					}
					Matmult(inv_one_minus_gamma_tilde_gamma,gamma_tilde,X18);
					Matmult(X18,X17,X19);
					Matmult(X19,inv_one_minus_gamma_tilde_gamma,X20);
					Matmult(X20,gamma_tilde,X21);
					Matmult(X21,gamma_prime,X22);
					Matmult(gamma_prime,X22,X23);
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime_deriv[q2][q3][q1] = -two * X23[q2][q3];
						}
					}
					Matmult(inv_one_minus_gamma_gamma_tilde,X17,X24);		
					Matmult(inv_one_minus_gamma_gamma_tilde,X17,X25);		
					Matmult(X25,gamma_tilde,X26);
					Matmult(X26,inv_one_minus_gamma_gamma_tilde,X27);
					Matmult(X27,gamma,X28);
					Matadd(X24,X28,X29);
					Matmult(X29,gamma_tilde_prime,X30);							
					Matmult(gamma_tilde_prime,X30,X31);							
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime_deriv[q2][q3][q1] = - two * X31[q2][q3];
						}
					}
				}
				MatComplex X32(2,2),X33(2,2),X34(2,2),X35(2,2),X36(2,2),X37(2,2),X38(2,2),X39(2,2),
				X40(2,2),X41(2,2),X42(2,2),X43(2,2),X44(2,2),X45(2,2),X46(2,2);
				for(int q1=8;q1<=11;q1++){		//derivatives after gamma tilde
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							X32[q2][q3] = gamma_tilde_deriv[q2][q3][q1];
						}
					}
					Matmult(inv_one_minus_gamma_gamma_tilde,gamma,X33);
					Matmult(X33,X32,X34);
					Matmult(X34,inv_one_minus_gamma_gamma_tilde,X35);		
					Matmult(X35,gamma,X36);
					Matmult(X36,gamma_tilde_prime,X37);
					Matmult(gamma_tilde_prime,X37,X38);
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime_deriv[q2][q3][q1] = -two * X38[q2][q3];
						}
					}
					Matmult(inv_one_minus_gamma_tilde_gamma,X32,X39);		
					Matmult(inv_one_minus_gamma_tilde_gamma,X32,X40);		
					Matmult(X40,gamma,X41);
					Matmult(X41,inv_one_minus_gamma_tilde_gamma,X42);
					Matmult(X42,gamma_tilde,X43);
					Matadd(X39,X43,X44);
					Matmult(X44,gamma_prime,X45);			
					Matmult(gamma_prime,X45,X46);			
					for(int q2=0;q2<=1;q2++){
						for(int q3=0;q3<=1;q3++){
							minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime_deriv[q2][q3][q1] = -two * X46[q2][q3];
						}
					}
				}


			//Implementation of the calculated derivatives in the s-matrix
			int qder,pder;
			MatComplex ONE(16,16);
			for(int j1=0;j1<=15;j1++){
				for(int j2=0;j2<=15;j2++){
					ONE[j1][j2] = 0;
				}
			}
			for(int zder=0;zder<=15;zder++){
				ONE[zder][indexv[zder]] = one;
			}
			for(int jder=4;jder<=7;jder++){
				if(jder==4){qder=0; pder=0;}
				if(jder==5){qder=0; pder=1;}
				if(jder==6){qder=1; pder=0;}
				if(jder==7){qder=1; pder=1;}

				for(int kder=0;kder<=15;kder++){			
					s[jder][indexv[kder]] = - ONE[jder][indexv[kder]]
											- h * minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime_deriv[qder][pder][kder]
					 						- h * i * gamma_Delta_conj_gamma_deriv[qder][pder][kder]
					 						- h * i * J_sigma_gamma_deriv[qder][pder][kder] 
					 						- h * i * minus_gamma_J_sigma_conj_deriv[qder][pder][kder]
					 						- h * i * (-two)*E*gamma_deriv[qder][pder][kder];												 
				}
				for(int kder=0;kder<=15;kder++){
					if(jder != kder){
						s[jder][16+indexv[kder]] = s[jder][indexv[kder]];
					}
					else{
						s[jder][16+indexv[kder]] = two + s[jder][indexv[kder]];
					}		
				}
			}
			for(int jder=12;jder<=15;jder++){
				if(jder==12){qder=0; pder=0;}
				if(jder==13){qder=0; pder=1;}
				if(jder==14){qder=1; pder=0;}
				if(jder==15){qder=1; pder=1;}
				for(int kder=0;kder<=15;kder++){
					s[jder][indexv[kder]] = - ONE[jder][indexv[kder]]
											- h * minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime_deriv[qder][pder][kder]
										    - h * (-i) * gamma_tilde_Delta_gamma_tilde_deriv[qder][pder][kder]
											- h * (-i) * J_sigma_conj_gamma_tilde_deriv[qder][pder][kder]
											- h * (-i) * minus_gamma_tilde_J_sigma_deriv[qder][pder][kder]
											- h * (-i) * two *E*gamma_tilde_deriv[qder][pder][kder];								 
				}
				for(int kder=0;kder<=15;kder++){
					if(jder != kder){
						s[jder][16+indexv[kder]] = s[jder][indexv[kder]];
					}
					else{
						s[jder][16+indexv[kder]] = two + s[jder][indexv[kder]];
					}		
				}
			}




			//now the implementation in the original equations
			int qor,por;
			for(int jor=4;jor<=7;jor++){
				if(jor==4){qor=0; por=0;}
				if(jor==5){qor=0; por=1;}
				if(jor==6){qor=1; por=0;}
				if(jor==7){qor=1; por=1;}
				s[jor][jsf] = y[jor][k] - y[jor][k-1] 
					- h * minus_two_gamma_prime_inv_one_minus_gamma_tilde_gamma_gamma_tilde_gamma_prime[qor][por]
					- h * i * (gamma_Delta_conj_gamma[qor][por] - Delta_mat[qor][por] 
					- two * E * gamma[qor][por] + J_sigma_gamma[qor][por] + minus_gamma_J_sigma_conj[qor][por]);
			}
			for(int jor=12;jor<=15;jor++){
				if(jor==12){qor=0; por=0;}
				if(jor==13){qor=0; por=1;}
				if(jor==14){qor=1; por=0;}
				if(jor==15){qor=1; por=1;}
				s[jor][jsf] = y[jor][k]- y[jor][k-1] 
					- h * minus_two_gamma_tilde_prime_inv_one_minus_gamma_gamma_tilde_gamma_gamma_tilde_prime[qor][por]
					- h * (-i) * (gamma_tilde_Delta_gamma_tilde[qor][por] - Delta_conj_mat[qor][por]
					+ two * E * gamma_tilde[qor][por] + J_sigma_conj_gamma_tilde[qor][por] + minus_gamma_tilde_J_sigma[qor][por]);
			}
		}	
	}
};
