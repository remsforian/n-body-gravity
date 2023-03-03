/*
Teddy Masters Copyright 2023
A numerical simulation of the n body problem using gravitational forces
This project was inspired by the novel "the three body problem" by Cixin Liu

Things I still want to do:
parallelization
*/
#include<stdio.h>
#include<cmath>
#include<iostream>
#include<cstdlib>
#define dt 0.000001
#define dtao 0.05
//N: number of bodies, I: number of diffeqs D: number of dimentions
#define N 3
#define I 2
#define D 3
//mass + g (made up units)
#define M 1
#define G 1
//1 for recentering 0 for not
#define CENTER 1
//for loops
#define FORI for(int i=0;i<I;i++)
#define FORD for(int d=0;d<D;d++)
#define FORN for(int n=0;n<N;n++)
#define FORP for(int p=0;p<N;p++)
//random stuff for initial conditions

void derivs(long double[N][I][D], long double[N][I][D]);
void dy(long double[N][I][D], long double[N][I][D],long double[N][I][D]);
void recenter(long double [N][I][D]);
long double randfunc(long double,long double,int);
long double POW(long double);
int main(){
	//create arrays
	long double y[N][I][D];
	long double dydt[N][I][D];
	long double deltay[N][I][D];
	//random stuff
	int RANDSEED = 19258;
	//create time variables
	double t;
	double tau;
	double tf = 100; //final time
	FORN{
		FORD{
			int seed = (RANDSEED)*(n+1)*(d+1);
			y[n][0][d] = randfunc(-10.,10.,seed);
			y[n][1][d] = randfunc(-1.,1.,seed);
		}
	}
	//file setup
	FILE* fptr;
	fptr = fopen("body_cords.dat","w+");
	//beginning output
	std::cout << "Program begining. Final t: " << tf << "\n"; //prints output for user 
	for(t = 0;t<tf;t+=dt){
		//update the values
		dy(y,dydt,deltay);
		FORN{
			FORI{
				FORD{
					y[n][i][d] += deltay[n][i][d];
				}
			}
		}
		if(t>=tau){//prints to file
			tau += dtao;
			fprintf(fptr,"%lf ",t); //prints time
			if(CENTER == 1){recenter(y);} //recenters if desired setting is there
			FORN{
				FORD{
					fprintf(fptr,"%Lf ",y[n][0][d]); //prints positions of all bodies
				}
			}
			fprintf(fptr,"\n");
			std::cout << "Current Program Time:" << t << "\n"; //prints current time to user
		}
	}
}
long double POW(long double x){
	return(x*x);
}
long double randfunc(long double lbound, long double ubound,int seed){
	srandom(seed);//random number seed
	const long max_rand = 1000000;
	return lbound + (ubound-lbound) * (random()%max_rand)/max_rand;
}
void recenter(long double y[N][I][D]){
	long double sum[D];
	long double avg[D];
	FORD{
		sum[d] = 0;
		avg[d] = 0;
	}
	FORN{
		FORD{
			sum[d] += y[n][0][d]; //finds the sum in each direction
		}
	}
	FORD{
		avg[d] = sum[d]/N;//calculates the center of mass for the system
	}
	FORD{
		FORN{
			y[n][0][d] -= avg[d]; //moves the bodies so they are in the center of the cord system
		}
	}
}
void dy(long double y[N][I][D], long double dydt[N][I][D], long double deltay[N][I][D]){ //fhte fuction where we numerically integrate
	//temporary variables for RK4 use
	long double f1[N][I][D]; 
	long double f2[N][I][D];
	long double f3[N][I][D];
	long double f4[N][I][D];
	long double df1[N][I][D];
	long double df2[N][I][D];
	long double df3[N][I][D];
	long double df4[N][I][D];
	derivs(y,dydt); //first step calc
	FORN{
		FORI{
			FORD{
				f1[n][i][d] = y[n][i][d]; //first RK
			}
		}
	}
	derivs(f1,df1); //second step calc
	FORN{
		FORI{
			FORD{
				f2[n][i][d] = y[n][i][d] + df1[n][i][d] * dt/2; //second RK
			}
		}
	}
	derivs(f2,df2); //third step calc
	FORN{
		FORI{
			FORD{
				f3[n][i][d] = y[n][i][d] + df2[n][i][d] * dt/2;//third RK
			}
		}
	}
	derivs(f3,df3); //fourth step calc
	FORN{
		FORI{
			FORD{
				f4[n][i][d] = y[n][i][d] + df3[n][i][d] * dt;//final RK
			}
		}
	}
	derivs(f4,df4); //final calc
	FORN{
		FORI{
			FORD{
				deltay[n][i][d] = (1./6.) * (df1[n][i][d] + 2 * df2[n][i][d] + 2 * df3[n][i][d] + df4[n][i][d]) * dt; //weighted avrage
			}
		}
	}
}
void derivs(long double y[N][I][D],long double dydt[N][I][D]){ //the function where the physics lives 
	long double r[N][D][D];//vector between bodies
	long double sqrmag[N][N];//square magnitude of r
	long double rmag[N][N]; //magnitude of r
	long double rhat[N][N][D];//unit vectors of r
	FORN{
		FORI{
			FORD{
				dydt[n][i][d] = 0; // clears the dydt array
			}
		}
	}
	FORN{
		FORP{
			sqrmag[n][p] = 0; //clears sqrmag variable
		}
	}
	FORN{
		FORD{
			dydt[n][0][d] = y[n][1][d]; //change in position from velocity
		}
	}
	FORN{
		FORP{
			FORD{
				r[n][p][d] = y[p][0][d] - y[n][0][d]; //finds the position vectors between pairs of bodies
			}
		}
	}
	FORN{
		FORP{
			FORD{
				sqrmag[n][p] += POW(r[n][p][d]); //calculates the square magnitude of the r vector
			}
		}
	}
	FORN{
		FORP{
			rmag[n][p] = sqrt(sqrmag[n][p]); //calculates the magnitude of the r vector
		}
	}
	FORN{
		FORP{
			FORD{
				if(n==p){continue;} //dosen't let the program devide by zero (which wouldn't be ideal)
				rhat[n][p][d] = r[n][p][d] / rmag[n][p]; //creates the appropreate unit vector
			}
		}
	}
	FORN{
		FORD{
			FORP{
				if(p==n){continue;} //this is a safety guard for dividing by zero
				dydt[n][1][d] += ((-G*M*M)/(POW(rmag[n][p]))*(rhat[n][p][d])); //adds in new accelerations
			}
		}
	}
}