#include<R.h>
#include<Rmath.h>
#include<Rinternals.h>

double gammln(double xx){
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  
  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for(j = 0; j <=5; j++) ser += cof[j]/++y;
  return(-tmp+log(2.5066282746310005*ser/x));
}

double factrl(int n){
  //void nrerror(char error_text[]);
  static int ntop = 4;
  static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
  int j;
  
  if (n < 0) printf("Negative factorial in routine factrl");
  if (n > 32) return (exp(gammln(n + 1.0)));
  while (ntop < n){
    j = ntop++;
    a[ntop] = a[j]*ntop;
  }
  return(a[n]);
}

//funcao log-fatorial
double factln(int n){
  static double a[101];
  
  if (n < 0) printf("Negative factorial in routine factln");
  if (n <= 1) return 0.0;
  if (n <= 100) return (a[n] ? a[n]:(a[n] = gammln(n + 1.0)));
  else return(gammln(n+1.0));
}

//funcao log-choose
double Lchoose(int n, int k){
  return(factln(n) - factln(k) - factln((n) - (k)));
}


double fun1(int n, int J, double m, double x){ //alonsoA5
  double gama, nu, lambda, lny;
  gama = (m)*(J - 1)/(1 - (m));
  nu = J + gama*(1 - (x));
  lambda = gama*(x);
  //printf("lchoose e %f\n", Lchoose(J, n));
  //printf("gama e %f\n", gama);
  //printf("nu e %f\n", nu);
  //printf("lambda e %f\n", lambda);
  lny = Lchoose(J, n) + gammln(n + lambda) - gammln(lambda) + gammln(nu - (n)) - gammln(nu - (J)) + gammln(lambda + nu - (J)) - gammln(lambda + nu);
  //printf("lny e %f\n", lny);
  //y = y*pow(1-x, theta-1)/x
  return exp(lny);
}

double SomaAreas(int n, int J, double m, double theta, double dx){
  double x, s;
  
  for(x = (dx/2), s = 0.0; x < 1; x += dx)
    s += fun1(n, J, m, x)*pow(1-x, theta-1)*dx/x;
    
  return s;
}

void Intgrl1(int *n, int *J, double *m, double *theta, double *precisao, double *result){
  double dx, s0, s1;
  dx = 1/100.;
  s0 = 0;
  do {
    //printf("integrando: %g passos, dx = %g\n", 1/dx, dx);
    s1 = s0;
    s0 = SomaAreas(*n, *J, *m, *theta, dx);
    dx /= 2;
  } while(abs(s1 - s0) > *precisao);
  *result = s0;
}

double Intgrl2(int n, int J, double m, double theta, double precisao){
  double dx, s0, s1;
  dx = 1/100.;
  s0 = 0;
  do {
    //printf("integrando: %g passos, dx = %g\n", 1/dx, dx);
    s1 = s0;
    s0 = SomaAreas(n, J, m, theta, dx);
    dx /= 2;
  } while(abs(s1 - s0) > precisao);
  return s0;
}

void sn(int *n, int *lengn, int *J, double *m, double *theta, double *precisao, double *result){
  int i;
  int *med;
  for (i=0; i < *lengn; i++){
    med = &n[i];
    //printf("%d\n", *med);
    result[i] = *theta*Intgrl2(*med, *J, *m, *theta, *precisao); //precisao 0.00001
    //printf("%g\n", result[i]);
  }
}

void msn(int *n, int *lengn, int *J, double *theta, double *result){
  int i;
  int *med;
  for (i=0; i < *lengn; i++){
    med = &n[i];
    result[i] = (*theta/(*med))*pow((1-(*med/(*J))), *theta-1) + *theta*(*theta-1)*(*theta-2)*pow((1-(*med/(*J))), *theta-3)/(*J*(*J)*2);
  }
}