#include <iostream>
#include <fstream>
#include <cstdlib>
#include "lib.h"
#include <cmath>
#include <omp.h>
#include <QElapsedTimer>
using namespace std;

double int_function_spherical(double theta1, double theta2, double phi1, double phi2, double rho1, double rho2);
double int_function_cartesian(double x1, double y1, double z1, double x2, double y2, double z2);
void gauss_laguerre(double *x, double *w, int n, double alf);
int main()
{
    double const pi = 3.14159265359;
    double const exact = 5*pi*pi/(16*16);
    //Task a)
    int n = 20;
    double a,b;
    a = -2;
    b = 2;
    double *x = new double [n];
    double *w = new double [n];


    gauleg(a, b, x ,w, n);
    double int_gauss = 0.;
    int i1,j1,k1,i2,j2,k2;
//    QElapsedTimer elapsedTimer;
//    elapsedTimer.start();

//#pragma omp for private(k2)
    for (i1 = 0; i1 < n; i1++){
        for (j1 = 0; j1 < n; j1++){
            for (k1 = 0; k1 < n; k1++){
                for (i2 = 0; i2 < n; i2++){
                    for (j2 = 0; j2 < n; j2++){
                        for (k2 = 0; k2 < n; k2++){
                            //cout << int_function_cartesian(x[i1], x[j1], x[k1], x[i2], x[j2], x[k2]) << endl;
                            int_gauss += w[i1]*w[j1]*w[k1]*w[i2]*w[j2]*w[k2]
                                    *int_function_cartesian(x[i1], x[j1], x[k1], x[i2], x[j2], x[k2]);
                        }
                    }
                }
            }
        }
        //cout << i1 << endl;
    }
    //cout << "Time: " << elapsedTimer.elapsed()/1000. << "s" << endl;
    cout << "Exact integral: " << exact << endl;
    cout << "Integral Gauss-Legendre: " << int_gauss << endl;
    cout << "Relative error Gauss-Legendre: " << (exact-int_gauss)/exact << endl;
    delete [] w;
    //Task b)
    double *rho = new double[n+1];
    double *wr = new double[n+1];
    double *theta = new double[n];
    double *wt = new double[n];
    double *phi = new double[n];
    double *wp = new double[n];

    gauleg(0.,pi,theta,wt,n);
    gauleg(0.,2*pi,phi,wp,n);
    gauss_laguerre(rho, wr, n, 2.0);
    for (i1 = 0; i1 < n; i1++){
            for (j1 = 0; j1 < n; j1++){
                for (k1 = 0; k1 < n; k1++){
                    for (i2 = 0; i2 < n; i2++){
                        for (j2 = 0; j2 < n; j2++){
                            for (k2 = 0; k2 < n; k2++){
                                int_gauss += wt[i1]*wt[j1]*wp[k1]*wp[i2]*wr[j2+1]*wr[k2+1]*
                                        int_function_spherical(theta[i1],theta[j1],phi[k1],phi[i2],rho[j2+1],rho[k2+1]);
                            }
                        }
                    }
                }
            }
            //cout << i1 << endl;
        }
    cout << "Integral Gauss-Leguerre: " << int_gauss/1024. << endl;
    cout << "Relative error Gauss-Leguerre: " << (exact-int_gauss/1024.)/exact << endl;
    delete []rho;
    delete []wr;
    delete []theta;
    delete []wt;
    delete []phi;
    delete []wp;
//    //Task c)
    int N = 1e6;
    long idum;
    double jacobiDeterminant = pow(4,6);
    idum = -1;
    double brute_integral = 0;
    double sigma_sum = 0;
    double f;
    //Brute force
    for(int i = 0; i < N; i++){
        for(int j = 0; j < 6; j++){
            x[j]= -2+4*ran0(&idum);
        }
        f = int_function_cartesian(x[0],x[1],x[2],x[3],x[4],x[5]);
        sigma_sum += f*f;
        brute_integral += f;
    }

    //Task d)
    //Importance sampling
    double sigma_sum_importance = 0;
    double importance_integral = 0;
    double theta1, theta2, phi1, phi2, rho1, rho2;
    double JacobiDeterminantSpherical = 4*pow(pi,4)/1024;
    for(int i = 0; i < N; i++){
        theta1 = pi*ran0(&idum);
        theta2 = pi*ran0(&idum);
        phi1 = 2*pi*ran0(&idum);
        phi2 = 2*pi*ran0(&idum);
        rho1 = -log(ran0(&idum));
        rho2 = log(ran0(&idum));
        f = rho1*rho1*rho2*rho2*int_function_spherical(theta1, theta2, phi1, phi2, rho1, rho2);
        importance_integral += f;
        sigma_sum_importance += f*f;
    }

    double var_brute;
    double var_importance;
    sigma_sum = sigma_sum/N;
    sigma_sum_importance = sigma_sum_importance/N;
    brute_integral = brute_integral/N;
    importance_integral = importance_integral/N;

    var_brute = (sigma_sum-brute_integral*brute_integral)/N;
    var_importance = (sigma_sum_importance-importance_integral*importance_integral)/N;
    cout << "Brute force MC: " << brute_integral*jacobiDeterminant << endl;
    cout << "Standard deviation brute force MC: " << jacobiDeterminant*sqrt(var_brute) << endl;
    cout << "Importance sampling MC: " << importance_integral*JacobiDeterminantSpherical << endl;
    cout << "Standard deviation importance MC: " << JacobiDeterminantSpherical*sqrt(var_importance) << endl;
    delete [] x;
    return 0;
}

double int_function_cartesian(double x1, double y1, double z1, double x2, double y2, double z2){
    double epsilon = 1e-12;

    double distance = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    double r1 = sqrt(x1*x1+y1*y1+z1*z1);
    double r2 = sqrt(x2*x2+y2*y2+z2*z2);

    if(distance > epsilon){
        return exp(-4*(r1+r2))/sqrt(distance);
    }else{return 0;}
}

double int_function_spherical(double theta1, double theta2, double phi1, double phi2, double rho1, double rho2){
    double epsilon = 1e-12;

    double cosB = (cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2));
    double distance2 = rho1*rho1 + rho2*rho2- 2*rho1*rho2*cosB;

    if(distance2 > epsilon){
        return sin(theta1)*sin(theta2)/sqrt(distance2);
    }else{return 0;}
}

