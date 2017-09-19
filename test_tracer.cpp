#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <vector>
#include <cfloat>
#include <iostream>
#include <fstream>

#define HEIGHT 168.7//cm
#define RADIUS  72.8
#define TrFOIL 0.//66802
double pmtRad[4] = {0.627,0.720,0.619,0.709}; //outer then inner
double gridsA[5] = { 0.20, 0.20, 0.20, 0.20, 0.20 };
double gridsZ[5] = { 1.00, 14.8, 160.4,161.7,167.7};
double gridsD[5] = { 75.0, 100., 75.0, 100., 75.0 }; //um
double gridsS[5] = { 0.50, 0.50, 0.50, 0.25, 0.50 };
double gridsT[5], mm2[2][60];
void histogram ( double milliMSq, int which );
#define BORDER 160.9
//#define PMTliveFrac 0.548 //Al ring counted as dead
#define TOL 1e-9 //precision 1 nm

using namespace std;

typedef enum {
  GXe = 0,
  LXe = 1,
  Tef = 2,
  Qtz = 3,
} MATERIALS;

typedef enum {
  up = 1,
  dn =-1,
} DIRECTIONS;

double QE; void Initialize ( ) {
  
  QE=.25;
  
}

double rand_gauss();
double QuantumEfficiency ( int array, double angle, unsigned index, double fudge ) {
  
  //if ( array == 0 ) return 0.25;//(0.27+rand_gauss()*0.05)*(1.-exp(-1.7/cos(angle)))/(1-exp(-1.7));
  //             else return 0.25;//(0.27+rand_gauss()*0.05)*(1.-exp(-1.7/cos(angle)))/(1-exp(-1.7));

  return QE;
  
}

double wavelength_nm ( double energy ) { //E, eV; lambda, nm
  
  return 6.626e-34 * 299792458./
    ( energy * 1e-9 * 1.602e-19 );

}

double direction[6]; void RaylScat ( bool nonRndm, double check );

double rand_uniform ( ) {

  return ( double ) rand ( ) / ( double ) RAND_MAX;

}

double rand_exp ( double mfp ) {
  
  return - mfp * log ( rand_uniform ( ) );
  
}

double rand_gauss ( ) {
  
  double u = rand_uniform(), v = rand_uniform();
  return sqrt(-2.*log(u))*cos(2.*M_PI*v);
  
}

double AngleToNormal(double origin[6]); //unsigned int FindWithinPMT(double PMTsize);
double dotP ( const double vec1[6], const double vec2[6] ) {
  
  return (vec1[3]-vec1[0])*(vec2[3]-vec2[0])+
         (vec1[4]-vec1[1])*(vec2[4]-vec2[1])+
         (vec1[5]-vec1[2])*(vec2[5]-vec2[2]);
  
}

double distanceSq ( double vect[6] ) {
  
  return (vect[3]-vect[0])*(vect[3]-vect[0])+
         (vect[4]-vect[1])*(vect[4]-vect[1])+
         (vect[5]-vect[2])*(vect[5]-vect[2]);
  
}

void truncation ( double level ) {
  
  double length2 = distanceSq ( direction );
  double scatter = ((level-direction[2])*sqrt(length2))/
    (direction[5]-direction[2]);
  direction[3] = direction[0] + scatter * (direction[3]-direction[0]) / sqrt(length2);
  direction[4] = direction[1] + scatter * (direction[4]-direction[1]) / sqrt(length2);
  direction[5] = direction[2] + scatter * (direction[5]-direction[2]) / sqrt(length2);
  
  return;
  
}

void reflection ( const double normal[6], bool diffuse ) {
  
  const double origin[6] =
    
    { direction[0],
      direction[1],
      direction[2],
      direction[3],
      direction[4],
      direction[5] }; int arrow;

  if ( (direction[5]-direction[2]) > 0. ) arrow = up;
                                     else arrow = dn;
  
  direction[0] = origin[3];
  direction[1] = origin[4];
  direction[2] = origin[5];

  if ( diffuse ) {
    
    //vector<double> pos(3);
  RE_TRY: RaylScat(false,0);
    //pos = randomDirection ( RADIUS );
    
    //direction[3] = pos[0];
    //direction[4] = pos[1];
    //direction[5] = pos[2];
    
    if ( fabs(sqrt(pow(direction[3],2.)+pow(direction[4],2.)) - RADIUS) > TOL ||
	 (fabs(direction[5]-direction[2]) < TOL &&
	  fabs(direction[4]-direction[1]) < TOL &&
	  fabs(direction[3]-direction[0]) < TOL) ) goto RE_TRY;
    
    if ( (direction[2] < (0.0000+TrFOIL+TOL) && direction[5] < direction[2]) ||
	 (direction[2] > (HEIGHT-TrFOIL-TOL) && direction[5] > direction[2]) ) goto RE_TRY;
    
    //if ( fabs(sqrt(pow(direction[0],2.)+pow(direction[1],2.)) - RADIUS) < TOL &&
    // ( (direction[2] < direction[5] && arrow == dn) || (direction[2] > direction[5] && arrow == up) )
    // && direction[2] > (0.0000+TrFOIL+TOL) && direction[2] < (HEIGHT-TOL-TrFOIL) ) { if(rand_uniform()<0.75) goto RE_TRY; }
    
    //if ( fabs(sqrt(pow(direction[0],2.)+pow(direction[1],2.)) - RADIUS) < TOL &&
    //if ( ((direction[2] > BORDER && direction[5] < BORDER) || (direction[2] < BORDER && direction[5] > BORDER)) && TIR ) goto RE_TRY;
    //&& direction[2] > (0.0000+TrFOIL+TOL) && direction[2] < (HEIGHT-TOL-TrFOIL) ) goto RE_TRY;
    
    return;
    
  }
  
  const double cosI = -dotP(normal, origin);
  
  const double transmit[3] = { (origin[3]-origin[0])+2.*cosI*(normal[3]-normal[0]),
			       (origin[4]-origin[1])+2.*cosI*(normal[4]-normal[1]),
			       (origin[5]-origin[2])+2.*cosI*(normal[5]-normal[2]) };
  
  double X = (sqrt((transmit[0]*transmit[0]+transmit[1]*transmit[1])*RADIUS*RADIUS-origin[3]*origin[3]*transmit[1]*transmit[1]+2.*origin[3]*origin[4]*transmit[0]*
		   transmit[1]-origin[4]*origin[4]*transmit[0]*transmit[0])-origin[3]*transmit[0]-origin[4]*transmit[1])/(transmit[0]*transmit[0]+transmit[1]*transmit[1]);
  
  direction[3] = origin[3] + X * transmit[0];
  direction[4] = origin[4] + X * transmit[1];
  direction[5] = origin[5] + X * transmit[2];
  
  return; // new dir
  
}

double refractiveIndex ( int material, double wavelength ) {
  
  switch ( material ) {
    
  case LXe:
    return 1.3805-9.5730/wavelength+12850./pow(wavelength,2.)-3.7518e6/pow(wavelength,3.)+6.4012e8/pow(wavelength,4.);
  case Qtz:
    return 1.5038+44.569/wavelength-17601./pow(wavelength,2.)+3.8406e6/pow(wavelength,3.)-1.9606e8/pow(wavelength,4.);
  case Tef:
    return 1.5;
  case GXe:
    return 1.000702;
    
  }
  
  return 1.0000; //if calling anything else
  
}

void refraction ( const double normal[6], int materials[2], double wavelength ) {
  
  double n1 = refractiveIndex ( materials[0], wavelength );
  double n2 = refractiveIndex ( materials[1], wavelength );
  
  const double magnitude = sqrt(distanceSq(direction));
  
  double origin[6] = { direction[0]/magnitude, direction[1]/magnitude, direction[2]/magnitude,
		       direction[3]/magnitude, direction[4]/magnitude, direction[5]/magnitude };
  
  const double n = n1 / n2;
  const double cosI = -dotP(normal, origin);
  const double sinT2 = n * n * ( 1. - cosI * cosI );
  if ( sinT2 > 1.0 ) {
    reflection ( normal, false ); return;
  }
  const double cosT = sqrt ( 1. - sinT2 );
  
  direction[0] = origin[3] * magnitude;
  direction[1] = origin[4] * magnitude;
  direction[2] = origin[5] * magnitude;
  
  const double transmit[3] = { n*(origin[3]-origin[0])+(n*cosI-cosT)*(normal[3]-normal[0]),
			       n*(origin[4]-origin[1])+(n*cosI-cosT)*(normal[4]-normal[1]),
			       n*(origin[5]-origin[2])+(n*cosI-cosT)*(normal[5]-normal[2]) };
  
  origin[3] *= magnitude;
  origin[4] *= magnitude;
  origin[5] *= magnitude;
  
  double X = (sqrt((transmit[0]*transmit[0]+transmit[1]*transmit[1])*RADIUS*RADIUS-origin[3]*origin[3]*transmit[1]*transmit[1]+2.*origin[3]*origin[4]*transmit[0]*
		   transmit[1]-origin[4]*origin[4]*transmit[0]*transmit[0])-origin[3]*transmit[0]-origin[4]*transmit[1])/(transmit[0]*transmit[0]+transmit[1]*transmit[1]);
  
  direction[3] = origin[3] + X * transmit[0];
  direction[4] = origin[4] + X * transmit[1];
  direction[5] = origin[5] + X * transmit[2];
  
  return; // the direction global variable has been updated
  
}

double reflectance ( int material );
double reflectivity ( int material, double wavelength,
		      int material2,const double normal[6] ) {
  
  double n1 = refractiveIndex ( material, wavelength );
  double n2 = refractiveIndex ( material2,wavelength );
  
  const double magnitude = sqrt(distanceSq(direction));
  
  double origin[6] = { direction[0]/magnitude, direction[1]/magnitude, direction[2]/magnitude,
                       direction[3]/magnitude, direction[4]/magnitude, direction[5]/magnitude };
  
  const double n = n1 / n2;
  const double cosI = -dotP(normal, origin);
  const double sinT2 = n * n * ( 1. - cosI * cosI );
  if ( sinT2 > 1.0 ) return 1.0; // total internal reflection
  const double cosT = sqrt ( 1. - sinT2 );
  
  const double r_s = ( n1 * cosI - n2 * cosT ) / ( n1 * cosI + n2 * cosT );
  const double r_p = ( n2 * cosI - n1 * cosT ) / ( n2 * cosI + n1 * cosT );
  
  return (r_s*r_s+r_p*r_p)/2.0; //averaging: orth and par
  
}

double reflectivity_Schlick ( int material, double wavelength,
			      int material2,const double normal[6] ) {
  
  double n1 = refractiveIndex ( material, wavelength );
  double n2 = refractiveIndex ( material2,wavelength );
  
  const double magnitude = sqrt(distanceSq(direction));
  
  double origin[6] = { direction[0]/magnitude, direction[1]/magnitude, direction[2]/magnitude,
                       direction[3]/magnitude, direction[4]/magnitude, direction[5]/magnitude };
  
  double r0 = ( n1 - n2 ) / ( n1 + n2 );
  r0 *= r0;
  double cosX = -dotP(normal, origin);
  
  if ( n1 > n2 ) {
    const double n = n1 / n2;
    const double sinT2 = n * n * ( 1. - cosX * cosX );
    if ( sinT2 > 1.0 ) return 1.0; // TIR
    cosX = sqrt ( 1. - sinT2 );
  }
  const double x = 1. - cosX;
  
  return r0 + ( 1. - r0 ) * x * x * x * x * x; //faster than pow func
  
}

double absorptionLength ( int material, double wavelength ) {
  
  switch ( material ) {
    
  case LXe:
    return 3e3; //in centimeters
  case Qtz:
    return 30.03082 + ( -0.301817 - 30.03082 ) / pow ( 1. + pow(wavelength/177.0802,10.5777), 0.0864376 ); //synthetic
  case Tef:
    return 0.;
  case GXe:
    return 500e2;
    
  }
  
  return 1e4; //default of "infinity" (100m)
  
}

double RaylScatLength ( int material, double wavelength ) {
  
  double temp = 732.79-13.494*wavelength+0.092174*pow(wavelength,2.)-0.00033424*pow(wavelength,3.)+6.6757e-07*pow(wavelength,4.);
  
  switch ( material ) {
    
  case LXe: return 30.;
    if ( temp < 20. ) return 20.0;
                 else return temp;
  default:
    return 500e2;
    
  } // the units are cm
  
}
double PMTmap[2][1456][1456];
double FindWithinPMT();
void Cartographer( double PMTmap[2][1456][1456]);


int main ( int argc, char** argv ) {
  
  int k, arrow, volume, surf[2]; vector<double> pos(3); double alpha, kludge; double dt[2];
  long numPhotons = (long)1e6, numSurvivors[2] = { 0, 0 }, numAbsorbed[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double E, nm, length2, perp[6], absorb, r, phi, scatter, X[8], minX, check; bool firstTime; double dz[2];
  unsigned index; double fudge;
   
  ofstream output1;
  output1.open("output1.txt");

  ofstream output2;
  output2.open("output2.txt"); 

  ofstream output3;
  output3.open("output3.txt"); 
 
  gridsT[0] = 1. - (2*gridsD[0]*1e-4) / gridsS[0];
  gridsT[1] = 1. - (2*gridsD[1]*1e-4) / gridsS[1];
  gridsT[2] = 1. - (2*gridsD[2]*1e-4) / gridsS[2];
  gridsT[3] = 1. - (2*gridsD[3]*1e-4) / gridsS[3];
  gridsT[4] = 1.;// - (1*gridsD[4]*1e-4) / gridsS[4];
  
  Initialize();
  Cartographer(PMTmap);
  
  for ( long i = 0; i < numPhotons; i++ ) { firstTime = true; //cout << i << endl;
    E = 6.97 + 0.23 * rand_gauss(); nm = wavelength_nm(E); //totalD[0]=0.;totalD[1]=0.;//absorb[0]=rand_exp(absorptionLength(GXe,nm));absorb[1]=rand_exp(absorptionLength(LXe,nm));
    phi = 2.*M_PI*rand_uniform(); r = 68.8 * sqrt ( rand_uniform() ); //double InitrValue = r;
    direction[0] = r * cos(phi);
    direction[1] = r * sin(phi);
    //dt[0] = 290.; dt[1] = 310.;
    //dz[0] = BORDER - dt[1] * 0.151; dz[1] = BORDER - dt[0] * 0.151;
    direction[2] = 16.3 + (146.9-16.3)*rand_uniform(); //atof(argv[1]);//BORDER + rand_uniform() * ( gridsZ[3] - BORDER );
    if ( direction[2] > BORDER ) { kludge = 2.; fudge = 1.000; }
                            else { kludge = 1.18; fudge = 1.000; }
    //printf("init %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    if ( direction[2] < BORDER ) volume = LXe;
                            else volume = GXe;
  RAYLEIGH: check = rand_uniform();
    if ( firstTime ) { RaylScat(false,0.000); firstTime = false; }
                  else RaylScat(true ,check);
    //printf("ray2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
  BOUNCE:
    if ( (direction[5]-direction[2]) > 0. ) arrow = up;
                                       else arrow = dn;
    X[0] = direction[2] - 0.0000;
    X[1] = direction[2] - gridsZ[0];
    X[2] = direction[2] - gridsZ[1];
    X[3] = direction[2] - gridsZ[2];
    X[4] = direction[2] - BORDER;
    X[5] = direction[2] - gridsZ[3];
    X[6] = direction[2] - gridsZ[4];
    X[7] = direction[2] - HEIGHT;
    minX = 1e100; k = -1;
    for ( int j = 0; j < 8; j++ ) {
      if(fabs(X[j]) < TOL) X[j]=0.;
      if ( fabs(X[j]) < minX && fabs(X[j]) > 0. ) {
	if ( arrow == up && X[j] < 0. ) { minX = X[j]; k = j; }
	if ( arrow == dn && X[j] > 0. ) { minX = X[j]; k = j; }
      }
      else { ; }
    } alpha = AngleToNormal ( direction ); //cout << alpha << endl;
  DO_OVER:
    switch ( k ) {
    case 0:
      if ( direction[2] > (0.0000+TOL) && direction[5] < (0.0000+TOL) ) goto PMTb;
      else goto REFL;
    case 1:
      if ( ((direction[2] > gridsZ[0] && direction[5] < gridsZ[0]) || (direction[2] < gridsZ[0] && direction[5] > gridsZ[0])) && rand_uniform() > gridsT[0]*1.01 ) {
	goto BOT;
      }
      if ( arrow == up ) { k = 2; goto DO_OVER; } else { k = 0; goto DO_OVER; }
    case 2:
      if ( ((direction[2] > gridsZ[1] && direction[5] < gridsZ[1]) || (direction[2] < gridsZ[1] && direction[5] > gridsZ[1])) && rand_uniform() > gridsT[1]*1.00 ) {
	goto CTH;
      }
      if ( arrow == up ) { k = 3; goto DO_OVER; } else { k = 1; goto DO_OVER; }
    case 3:
      if ( ((direction[2] > gridsZ[2] && direction[5] < gridsZ[2]) || (direction[2] < gridsZ[2] && direction[5] > gridsZ[2])) && rand_uniform() > gridsT[2]*1.00 ) {
	goto GAT;
      }
      if ( arrow == up ) { k = 4; goto DO_OVER; } else { k = 2; goto DO_OVER; }
    case 4:
      if ( (direction[2] < BORDER && direction[5] > BORDER) || (direction[2] > BORDER && direction[5] < BORDER) ) { //cout << (180./M_PI)*alpha << endl;
	goto REFR;
      }
      if ( arrow == up ) { k = 5; goto DO_OVER; } else { k = 3; goto DO_OVER; }
    case 5:
      if ( ((direction[2] > gridsZ[3] && direction[5] < gridsZ[3]) || (direction[2] < gridsZ[3] && direction[5] > gridsZ[3])) && rand_uniform() > gridsT[3]*pow(cos(alpha/2./kludge),2) ) {
	goto ANE;
      }
      if ( arrow == up ) { k = 6; goto DO_OVER; } else { k = 4; goto DO_OVER; }
    case 6:
      if ( ((direction[2] > gridsZ[4] && direction[5] < gridsZ[4]) || (direction[2] < gridsZ[4] && direction[5] > gridsZ[4])) && rand_uniform() > gridsT[4]*1.00 ) {
	goto TOP;
      }
      if ( arrow == up ) { k = 7; goto DO_OVER; } else { k = 5; goto DO_OVER; }
    case 7:
      if ( direction[2] < (HEIGHT-TOL) && direction[5] > (HEIGHT-TOL) ) goto PMTt;
      else goto REFL;
    default:
      goto REFL;
    }
  PMTt:
    truncation ( HEIGHT ); length2 = distanceSq ( direction ); //totalD[0] += length2;
    //printf("top1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("top2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(GXe,nm));
    scatter = rand_exp(RaylScatLength(GXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    } PMTtt:
    perp[3] = direction[3];
    perp[4] = direction[4];
    perp[5] = 0.;
    perp[0] = direction[3];
    perp[1] = direction[4];
    perp[2] = 1.;
    if ( FindWithinPMT() > 0.1 ) {
      if ( rand_uniform() < reflectivity(GXe,nm,Qtz,perp) )
	{ reflection ( perp, false ); goto BOUNCE; }
      if ( FindWithinPMT() > 0.7 ) { //if(alpha>0.35)alpha=0.35;
	alpha = asin ( (refractiveIndex(GXe,nm)/refractiveIndex(Qtz,nm)) * sin ( alpha ) );
	if ( rand_uniform() < QuantumEfficiency(0,alpha,index,fudge) ){
          numSurvivors[0]++;           
          output1 << (10*direction[3]) + 728 << "    ";
          output1 << (10*direction[4]) + 728 << endl;
          output3 << direction[5] << endl;
        }     
        else numAbsorbed[1]++;
      }
      else numAbsorbed[2]++;
    }
    else {
      if( rand_uniform() < 0.80 ) {
	direction[5] -= TrFOIL;
	reflection ( perp, true );
	//printf("top3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
	//printf("top4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
	goto BOUNCE;
      }
      else numAbsorbed[3]++;
    }
    continue;
  PMTb:
    truncation ( 0.0000 ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("bot1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("bot2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    } PMTbb:
    perp[3] = direction[3];
    perp[4] = direction[4];
    perp[5] = 1.;
    perp[0] = direction[3];
    perp[1] = direction[4];
    perp[2] = 0.;
    if ( FindWithinPMT() > 0.1 ) {
      if ( rand_uniform() < reflectivity(LXe,nm,Qtz,perp) )
	{ reflection ( perp, false ); goto BOUNCE; }
      if ( FindWithinPMT() > 0.7 ) { //if(alpha>0.35)alpha=0.35;
        alpha = asin ( (refractiveIndex(LXe,nm)/refractiveIndex(Qtz,nm)) * sin ( alpha ) );
	if ( rand_uniform() < QuantumEfficiency(1,alpha,index,fudge) ){
          numSurvivors[1]++;
          output2 << (10*direction[3]) + 728 << "    ";
          output2 << (10*direction[4]) + 728 << endl;
          output3 << direction[5] << endl;
      }
       else numAbsorbed[5]++;
      }
      else numAbsorbed[6]++;
    }
    else {
      if( rand_uniform() < reflectance(LXe) ) {
	direction[5] += TrFOIL;
	reflection ( perp, true );
	//printf("bot3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
        //printf("bot4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
	goto BOUNCE;
      }
      else numAbsorbed[7]++;
    }
    continue;
  BOT:
    truncation ( gridsZ[0] ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("bot1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("bot2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[0] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[8]++; continue; }
  CTH:
    truncation ( gridsZ[1] ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("cth1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("cth2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[1] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[9]++; continue; }
  GAT:
    truncation ( gridsZ[2] ); length2 = distanceSq ( direction ); //totalD[1] += length2;
    //printf("gat1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("gat2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(LXe,nm));
    scatter = rand_exp(RaylScatLength(LXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[2] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[10]++; continue; }
  ANE:
    truncation ( gridsZ[3] ); length2 = distanceSq ( direction ); //totalD[0] += length2;
    //printf("ano1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("ano2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(GXe,nm));
    scatter = rand_exp(RaylScatLength(GXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[3] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[11]++; continue; }
  TOP:
    truncation ( gridsZ[4] ); length2 = distanceSq ( direction ); //totalD[0] += length2;
    //printf("gtp1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("gtp2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(GXe,nm));
    scatter = rand_exp(RaylScatLength(GXe,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < gridsA[4] ) {
      perp[0] = direction[0]; perp[1] = direction[1]; perp[2] = direction[2];
      perp[3] = direction[3]; perp[4] = direction[4]; perp[5] = direction[5];
      while ( (direction[2] < direction[5] && arrow == up) || (direction[5] < direction[2] && arrow == dn) ) {
	direction[0] = perp[0]; direction[1] = perp[1]; direction[2] = perp[2]; direction[3] = perp[3]; direction[4] = perp[4]; direction[5] = perp[5];
	reflection ( perp, true );
      } goto BOUNCE;
    }
    else { numAbsorbed[12]++; continue; }
  REFR:
    truncation ( BORDER ); length2 = distanceSq ( direction ); //totalD[volume] += length2;
    //printf("rer1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("rer2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(volume,nm));
    scatter = rand_exp(RaylScatLength(volume,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { if ( volume == GXe ) numAbsorbed[0]++; else numAbsorbed[4]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ray1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( direction[2] < BORDER ) {
      perp[3] = direction[3];
      perp[4] = direction[4];
      perp[5] = 0.;
      perp[0] = direction[3];
      perp[1] = direction[4];
      perp[2] = 1.;
      surf[0] = LXe; surf[1] = GXe; refraction(perp, surf, nm);
      //printf("rer3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      //printf("rer4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
      if ( direction[5] > BORDER ) volume = GXe;
    }
    else {
      perp[3] = direction[3];
      perp[4] = direction[4];
      perp[5] = 1.;
      perp[0] = direction[3];
      perp[1] = direction[4];
      perp[2] = 0.;
      surf[0] = GXe; surf[1] = LXe; refraction(perp, surf, nm);
      //printf("rer3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      //printf("rer4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
      if ( direction[5] < BORDER ) volume = LXe;
    }
    goto BOUNCE;
  REFL:
    length2 = distanceSq ( direction ); //totalD[volume] += length2;
    //printf("ref1 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
    //printf("ref2 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
    absorb = rand_exp(absorptionLength(volume,nm));
    scatter = rand_exp(RaylScatLength(volume,nm));
    if ( absorb*absorb < length2 && absorb < scatter ) { if ( volume == LXe ) numAbsorbed[4]++; else numAbsorbed[0]++; continue; }
    if ( scatter < absorb && scatter*scatter < length2 ) {
      direction[0] += scatter * (direction[3]-direction[0]) / sqrt(length2);
      direction[1] += scatter * (direction[4]-direction[1]) / sqrt(length2);
      direction[2] += scatter * (direction[5]-direction[2]) / sqrt(length2);
      //printf("ry1' %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      goto RAYLEIGH;
    }
    if ( rand_uniform() < reflectance(volume) ) {
      r = sqrt(direction[3]*direction[3]+direction[4]*direction[4]);
      perp[3] = 0.;
      perp[4] = 0.;
      perp[5] = direction[5] / r;
      perp[0] = direction[3] / r;
      perp[1] = direction[4] / r;
      perp[2] = direction[5] / r;
      if(rand_uniform()<0.5)reflection ( perp, true ); else reflection ( perp, false );
      //printf("ref3 %.6f\t%.6f\t%.6f\n",direction[0],direction[1],direction[2]);
      //printf("ref4 %.6f\t%.6f\t%.6f\n",direction[3],direction[4],direction[5]);
      goto BOUNCE;
    }
    else {
      if ( volume == GXe ) numAbsorbed[13]++; else numAbsorbed[14]++; continue;
    }
  } //end of photons loop
 
  output1.close();
  output2.close();
  output3.close(); 

  cout << "gas\t\t\t" << numAbsorbed[0]  << endl;
  cout << "photo-cathode (T)\t" << numAbsorbed[1]  << endl;
  cout << "window (T)\t\t" << numAbsorbed[2]  << endl;
  cout << "trifoil (T)\t\t" << numAbsorbed[3]  << endl;
  cout << "liquid\t\t\t" << numAbsorbed[4]  << endl;
  cout << "photo-cathode (B)\t" << numAbsorbed[5]  << endl;
  cout << "window (B)\t\t" << numAbsorbed[6]  << endl;
  cout << "trifoil (B)\t\t" << numAbsorbed[7]  << endl;
  cout << "bot grid\t\t" << numAbsorbed[8]  << endl;
  cout << "cathode\t\t\t" << numAbsorbed[9]  << endl;
  cout << "gate\t\t\t" << numAbsorbed[10] << endl;
  cout << "anode\t\t\t" << numAbsorbed[11] << endl;
  cout << "top grid\t\t" << numAbsorbed[12] << endl;
  cout << "PTFE wall in gas\t" << numAbsorbed[13] << endl;
  cout << "PTFE wall in liquid\t" << numAbsorbed[14] << endl;
  
  cout << "#sur" << "\t" << "#abs" << "\t" << "#tot" << "\t" << "g1" << "\t\t" << "[%]" << "\t" << "TBA" << endl;
  double total = 1. * double(numSurvivors[0]) + 1. * double(numSurvivors[1]);
  long noAbs = numAbsorbed[0] + numAbsorbed[1] + numAbsorbed[2] + numAbsorbed[3] + numAbsorbed[4] + numAbsorbed[5] + numAbsorbed[6] +
    numAbsorbed[7] + numAbsorbed[8] + numAbsorbed[9] + numAbsorbed[10]+ numAbsorbed[11]+ numAbsorbed[12]+ numAbsorbed[13] +numAbsorbed[14];
  cout << total << "\t" << noAbs << "\t" << numSurvivors[0]+numSurvivors[1]+noAbs << "\t" << total
    /double(numPhotons) << "\t" << 100.*total
    /double(numPhotons) << "\t" <<
    ( 1. * double(numSurvivors[0]) - 1. * double(numSurvivors[1]) ) / total << endl;
  
  //for ( int pq = 0; pq < 60; pq++ ) cout << mm2[0][pq] << "\t" << mm2[1][pq] << endl;
  
  return 1;
  
}

double reflectance ( int material ) {
  
  switch ( material ) {
    
  case LXe:
    return 0.95;
  case GXe:
    return 0.20;
    
  }
  
  return 0.; //if calling anything else make it non-reflective
  
}

double AngleToNormal ( double origin[6] ) {
  
  double normal[6];
  
  normal[0] = origin[3]; normal[1] = origin[4]; normal[2] = 0.;
  normal[3] = origin[3]; normal[4] = origin[4]; normal[5] = 1.;
  
  double numer = fabs(dotP(normal,origin));
  double denom = sqrt(distanceSq(normal)*distanceSq(origin));
  
  return acos ( numer / denom ); //in radians!

}

void RaylScat ( bool nonRndm, double check ) {
  
  double CosTheta=rand_uniform(); double portion=.5;
  if ( nonRndm && check < 0.5 ) {
    double dist = rand_uniform() + 1.;
    CosTheta = sqrt ( dist - 1. ); portion = 0.499;
    //CosTheta = 1.;
    //CosTheta = 2. * rand_uniform() - 1.;
    //fcostheta = ( 1. + CosTheta*CosTheta)/2.;
  }
  double SinTheta = sqrt(1.-CosTheta*CosTheta);
  // consider for the angle 90-180 degrees
  if (rand_uniform() < portion ) CosTheta = -CosTheta; //if(nonRndm) cout << CosTheta << endl;

  // simulate the phi angle
  double rand = 2.*M_PI*rand_uniform();
  double SinPhi = sin(rand);
  double CosPhi = cos(rand);
  
  // start constructing the new momentum direction
  double unit_x = SinTheta * CosPhi;
  double unit_y = SinTheta * SinPhi;
  double unit_z = CosTheta;
  
  double alpha = (sqrt(-(direction[1]*direction[1]-RADIUS*RADIUS)*unit_x*unit_x+2.*direction[0]*direction[1]*unit_x*unit_y-(direction[0]*direction[0]-RADIUS*RADIUS)*unit_y*
		       unit_y)-direction[0]*unit_x-direction[1]*unit_y)/(unit_x*unit_x+unit_y*unit_y); //extend the unit vector until ray of light hits the wall at RADIUS
  
  direction[3] = alpha * unit_x + direction[0];
  direction[4] = alpha * unit_y + direction[1];
  direction[5] = alpha * unit_z + direction[2];
  
  return; //void because dir is global var
  
}

/*unsigned int FindWithinPMT ( double PMTsize ) {
  
  double PMTy[5] = { 0., 5.19615, 10.3923, 15.5885, 20.7846 };
  double PMTx[9] = { 0., 3., 6., 9., 12., 15., 18., 21., 24 };
  
  int i, minX, minY; double minV = DBL_MAX; double sign[2];
  for ( i = 0; i < 5; i++ ) {
    if ( fabs ( direction[4] - PMTy[i] ) < minV ) { minV = fabs ( direction[4] - PMTy[i] ); minY = i; sign[1] = 1.; }
  }
  for ( i = 0; i < 5; i++ ) {
    if ( fabs ( direction[4] + PMTy[i] ) < minV ) { minV = fabs ( direction[4] + PMTy[i] ); minY = i; sign[1] =-1.; }
  }
  if ( minV > PMTsize ) return 0;
  
  if ( minY == 0 ) { PMTx[1] = DBL_MAX; PMTx[3] = DBL_MAX; PMTx[5] = DBL_MAX; PMTx[7] = DBL_MAX; }
  else if ( minY == 1 ) { PMTx[0] = DBL_MAX; PMTx[2] = DBL_MAX; PMTx[4] = DBL_MAX; PMTx[6] = DBL_MAX; PMTx[8] = DBL_MAX; }
  else if ( minY == 2 ) { PMTx[1] = DBL_MAX; PMTx[3] = DBL_MAX; PMTx[5] = DBL_MAX; PMTx[7] = DBL_MAX; PMTx[8] = DBL_MAX; }
  else if ( minY == 3 ) { PMTx[0] = DBL_MAX; PMTx[2] = DBL_MAX; PMTx[4] = DBL_MAX; PMTx[6] = DBL_MAX; PMTx[7] = DBL_MAX; PMTx[8] = DBL_MAX; }
  else { PMTx[1] = DBL_MAX; PMTx[3] = DBL_MAX; PMTx[5] = DBL_MAX; PMTx[6] = DBL_MAX; PMTx[7] = DBL_MAX; PMTx[8] = DBL_MAX; }
  
  minV = DBL_MAX;
  for ( i = 0; i < 9; i++ ) {
    if ( fabs ( direction[3] - PMTx[i] ) < minV ) { minV = fabs ( direction[3] - PMTx[i] ); minX = i; sign[0] = 1.; }
  }
  for ( i = 0; i < 9; i++ ) {
    if ( fabs ( direction[3] + PMTx[i] ) < minV ) { minV = fabs ( direction[3] + PMTx[i] ); minX = i; sign[0] =-1.; }
  }
  if ( minV > PMTsize ) return 0;
  
  if ( ( (direction[3]-sign[0]*PMTx[minX])*(direction[3]-sign[0]*PMTx[minX]) +
	 (direction[4]-sign[1]*PMTy[minY])*(direction[4]-sign[1]*PMTy[minY]) ) < PMTsize*PMTsize )
    return int(floor(sign[0]*PMTx[minX]+1.5*sign[1]*PMTy[minY])+45.);
  else
    return 0;
  */
//}

void histogram ( double milliMSq, int which ) {
  
  if ( milliMSq > 0.0000 && milliMSq <= 1000.0 ) mm2[which][0]++ ;
  if ( milliMSq > 1000.0 && milliMSq <= 2000.0 ) mm2[which][1]++ ;
  if ( milliMSq > 2000.0 && milliMSq <= 3000.0 ) mm2[which][2]++ ;
  if ( milliMSq > 3000.0 && milliMSq <= 4000.0 ) mm2[which][3]++ ;
  if ( milliMSq > 4000.0 && milliMSq <= 5000.0 ) mm2[which][4]++ ;
  if ( milliMSq > 5000.0 && milliMSq <= 6000.0 ) mm2[which][5]++ ;
  if ( milliMSq > 6000.0 && milliMSq <= 7000.0 ) mm2[which][6]++ ;
  if ( milliMSq > 7000.0 && milliMSq <= 8000.0 ) mm2[which][7]++ ;
  if ( milliMSq > 8000.0 && milliMSq <= 9000.0 ) mm2[which][8]++ ;
  if ( milliMSq > 9000.0 && milliMSq <= 10000. ) mm2[which][9]++ ;
  if ( milliMSq > 10000. && milliMSq <= 11000. ) mm2[which][10]++;
  if ( milliMSq > 11000. && milliMSq <= 12000. ) mm2[which][11]++;
  if ( milliMSq > 12000. && milliMSq <= 13000. ) mm2[which][12]++;
  if ( milliMSq > 13000. && milliMSq <= 14000. ) mm2[which][13]++;
  if ( milliMSq > 14000. && milliMSq <= 15000. ) mm2[which][14]++;
  if ( milliMSq > 15000. && milliMSq <= 16000. ) mm2[which][15]++;
  if ( milliMSq > 16000. && milliMSq <= 17000. ) mm2[which][16]++;
  if ( milliMSq > 17000. && milliMSq <= 18000. ) mm2[which][17]++;
  if ( milliMSq > 18000. && milliMSq <= 19000. ) mm2[which][18]++;
  if ( milliMSq > 19000. && milliMSq <= 20000. ) mm2[which][19]++;
  if ( milliMSq > 20000. && milliMSq <= 21000. ) mm2[which][20]++;
  if ( milliMSq > 21000. && milliMSq <= 22000. ) mm2[which][21]++;
  if ( milliMSq > 22000. && milliMSq <= 23000. ) mm2[which][22]++;
  if ( milliMSq > 23000. && milliMSq <= 24000. ) mm2[which][23]++;
  if ( milliMSq > 24000. && milliMSq <= 25000. ) mm2[which][24]++;
  if ( milliMSq > 25000. && milliMSq <= 26000. ) mm2[which][25]++;
  if ( milliMSq > 26000. && milliMSq <= 27000. ) mm2[which][26]++;
  if ( milliMSq > 27000. && milliMSq <= 28000. ) mm2[which][27]++;
  if ( milliMSq > 28000. && milliMSq <= 29000. ) mm2[which][28]++;
  if ( milliMSq > 29000. && milliMSq <= 30000. ) mm2[which][29]++;
  if ( milliMSq > 30000. && milliMSq <= 31000. ) mm2[which][30]++;
  if ( milliMSq > 31000. && milliMSq <= 32000. ) mm2[which][31]++;
  if ( milliMSq > 32000. && milliMSq <= 33000. ) mm2[which][32]++;
  if ( milliMSq > 33000. && milliMSq <= 34000. ) mm2[which][33]++;
  if ( milliMSq > 34000. && milliMSq <= 35000. ) mm2[which][34]++;
  if ( milliMSq > 35000. && milliMSq <= 36000. ) mm2[which][35]++;
  if ( milliMSq > 36000. && milliMSq <= 37000. ) mm2[which][36]++;
  if ( milliMSq > 37000. && milliMSq <= 38000. ) mm2[which][37]++;
  if ( milliMSq > 38000. && milliMSq <= 39000. ) mm2[which][38]++;
  if ( milliMSq > 39000. && milliMSq <= 40000. ) mm2[which][39]++;
  if ( milliMSq > 40000. && milliMSq <= 41000. ) mm2[which][40]++;
  if ( milliMSq > 41000. && milliMSq <= 42000. ) mm2[which][41]++;
  if ( milliMSq > 42000. && milliMSq <= 43000. ) mm2[which][42]++;
  if ( milliMSq > 43000. && milliMSq <= 44000. ) mm2[which][43]++;
  if ( milliMSq > 44000. && milliMSq <= 45000. ) mm2[which][44]++;
  if ( milliMSq > 45000. && milliMSq <= 46000. ) mm2[which][45]++;
  if ( milliMSq > 46000. && milliMSq <= 47000. ) mm2[which][46]++;
  if ( milliMSq > 47000. && milliMSq <= 48000. ) mm2[which][47]++;
  if ( milliMSq > 48000. && milliMSq <= 49000. ) mm2[which][48]++;
  if ( milliMSq > 49000. && milliMSq <= 50000. ) mm2[which][49]++;
  if ( milliMSq > 50000. && milliMSq <= 51000. ) mm2[which][50]++;
  if ( milliMSq > 51000. && milliMSq <= 52000. ) mm2[which][51]++;
  if ( milliMSq > 52000. && milliMSq <= 53000. ) mm2[which][52]++;
  if ( milliMSq > 53000. && milliMSq <= 54000. ) mm2[which][53]++;
  if ( milliMSq > 54000. && milliMSq <= 55000. ) mm2[which][54]++;
  if ( milliMSq > 55000. && milliMSq <= 56000. ) mm2[which][55]++;
  if ( milliMSq > 56000. && milliMSq <= 57000. ) mm2[which][56]++;
  if ( milliMSq > 57000. && milliMSq <= 58000. ) mm2[which][57]++;
  if ( milliMSq > 58000. && milliMSq <= 59000. ) mm2[which][58]++;
  if ( milliMSq > 59000. && milliMSq <= 60000. ) mm2[which][59]++;
  
}



double FindWithinPMT () {
  if ( direction[5] > 10 ) {
    int xo = (int) (10*direction[3]) + 728 + 0.5; int yo = (int) (10*direction[4]) + 728 + 0.5;
    return PMTmap[0][xo][yo];
  }
  if  ( direction[5] <= 10 ) {
    int xo = (int) (10*direction[3]) + 728 + 0.5; int yo = (int) (10*direction[4]) + 728 + 0.5;
    return PMTmap[1][xo][yo];
  }
}      //Will return 0 if xy is not a PMT location; 0.5 if inside only outer radius; 1.0 if inside inner radius


void Cartographer ( double PMTmap[2][1456][1456] ) {
  //double PMTmap[2][1456][1456];
  //int top = 0; int bot = 1;
  int PMTrad[2] = {38,32}; // Outer and Inner radii of the PMTs

  for (int i = 0; i < 1456; i++ ) {
    for (int  j = 0; j < 1456; j++ ) {
      PMTmap[0][i][j] = 0.0;
      PMTmap[1][i][j] = 0.0;           // fills the entire array with zeros
    }
  }


  int data_top[][2] = {{728,728},{820,728},{774,808},{682,808},{636,728},{682,648},{774,648},{866,808},{728,887},{590,808},{590,648},{728,569},{866,648},{912,728},{820,887},{636,887},{544,728},{636,569},{820,569},{958,808},{912,887},{774,967},{682,967},{544,887},{498,808},{498,648},{544,569},{682,489},{774,489},{912,569},{958,648},{1004,728},{866,967},{590,967},{452,728},{590,489},{866,489},{1004,887},{728,1047},{452,887},{452,569},{728,409},{1004,569},{1050,808},{958,967},{820,1047},{636,1047},{498,967},{406,808},{406,648},{498,489},{636,409},{820,409},{958,489},{1050,648},{1096,728},{912,1047},{544,1047},{360,728},{544,409},{912,409},{1095,887},{1049,966},{774,1125},{682,1125},{407,966},{361,887},{361,569},{407,490},{682,331},{774,331},{1049,490},{1095,569},{1136,805},{999,1043},{866,1120},{590,1120},{457,1043},{320,805},{320,651},{457,413},{590,336},{866,336},{999,413},{1136,651},{1202,728},{965,1138},{491,1138},{254,728},{491,318},{965,318},{1141,966},{728,1204},{315,966},{315,490},{728,252},{1141,490},{1185,891},{1097,1042},{816,1205},{640,1205},{359,1042},{271,891},{271,565},{359,414},{640,251},{816,251},{1097,414},{1185,565},{1223,813},{1049,1114},{902,1199},{554,1199},{407,1114},{233,813},{233,643},{407,342},{554,257},{902,257},{1049,342},{1223,643},{1232,969},{1189,1044},{772,1285},{684,1285},{267,1044},{224,969},{224,487},{267,412},{684,171},{772,171},{1189,412},{1232,487},{1280,877},{1133,1132},{875,1281},{581,1281},{323,1132},{176,877},{176,579},{323,324},{581,175},{875,175},{1133,324},{1280,579},{1299,771},{1051,1201},{976,1244},{480,1244},{405,1201},{157,771},{157,685},{405,255},{480,212},{976,212},{1051,255},{1299,685},{1319,973},{1236,1117},{811,1362},{645,1362},{220,1117},{137,973},{137,483},{220,339},{645,94},{811,94},{1236,339},{1319,483},{1374,814},{1125,1245},{977,1331},{479,1331},{331,1245},{82,814},{82,642},{331,211},{479,125},{977,125},{1125,211},{1374,642},{1364,899},{1194,1193},{898,1364},{558,1364},{262,1193},{92,899},{92,557},{262,263},{558,92},{898,92},{1194,263},{1364,557},{1387,728},{1299,1058},{1058,1299},{728,1387},{398,1299},{157,1058},{69,728},{157,398},{398,157},{728,69},{1058,157},{1299,398},{1459,776},{1447,871},{1422,964},{1385,1052},{1337,1135},{1279,1211},{1211,1279},{1135,1337},{1052,1385},{964,1422},{871,1447},{776,1459},{680,1459},{585,1447},{492,1422},{404,1385},{321,1337},{245,1279},{177,1211},{119,1135},{71,1052},{34,964},{9,871},{-3,776},{-3,680},{9,585},{34,492},{71,404},{119,321},{177,245},{245,177},{321,119},{404,71},{492,34},{585,9},{680,-3},{776,-3},{871,9},{964,34},{1052,71},{1135,119},{1211,177},{1279,245},{1337,321},{1385,404},{1422,492},{1447,585},{1459,680}}; 
  //data[i][0] are x-coordinates
  //data[i][1] are y-coordinates
 
 int data_bot[][2] = {{728,728},{799,769},{728,811},{657,769},{657,687},{728,645},{799,687},{871,728},{799,852},{657,852},{585,728},{657,604},{799,604},{871,811},{728,893},{585,811},{585,645},{728,563},{871,645},{942,769},{871,893},{799,934},{657,934},{585,893},{514,769},{514,687},{585,563},{657,522},{799,522},{871,563},{942,687},{942,852},{728,976},{514,852},{514,604},{728,480},{942,604},{1014,728},{871,976},{585,976},{442,728},{585,480},{871,480},{1014,811},{942,934},{799,1017},{657,1017},{514,934},{442,811},{442,645},{514,522},{657,439},{799,439},{942,522},{1014,645},{1014,893},{728,1058},{442,893},{442,563},{728,398},{1014,563},{1085,769},{942,1017},{871,1058},{585,1058},{514,1017},{371,769},{371,687},{514,439},{585,398},{871,398},{942,439},{1085,687},{1085,852},{1014,976},{799,1099},{657,1099},{442,976},{371,852},{371,604},{442,480},{657,357},{799,357},{1014,480},{1085,604},{1085,934},{728,1141},{371,934},{371,522},{728,315},{1085,522},{1157,728},{942,1099},{514,1099},{299,728},{514,357},{942,357},{1157,811},{1014,1058},{871,1141},{585,1141},{442,1058},{299,811},{299,645},{442,398},{585,315},{871,315},{1014,398},{1157,645},{1157,893},{1085,1017},{799,1182},{657,1182},{371,1017},{299,893},{299,563},{371,439},{657,274},{799,274},{1085,439},{1157,563},{1157,976},{728,1223},{299,976},{299,480},{728,233},{1157,480},{1228,769},{1014,1141},{942,1182},{514,1182},{442,1141},{228,769},{228,687},{442,315},{514,274},{942,274},{1014,315},{1228,687},{1228,852},{1085,1099},{871,1223},{585,1223},{371,1099},{228,852},{228,604},{371,357},{585,233},{871,233},{1085,357},{1228,604},{1228,934},{1157,1058},{799,1265},{657,1265},{299,1058},{228,934},{228,522},{299,398},{657,191},{799,191},{1157,398},{1228,522},{1300,728},{1014,1223},{442,1223},{156,728},{442,233},{1014,233},{1300,811},{1228,1017},{1085,1182},{942,1265},{728,1306},{514,1265},{371,1182},{228,1017},{156,811},{156,645},{228,439},{371,274},{514,191},{728,150},{942,191},{1085,274},{1228,439},{1300,645},{1300,893},{1157,1141},{871,1306},{585,1306},{299,1141},{156,893},{156,563},{299,315},{585,150},{871,150},{1157,315},{1300,563},{1300,976},{1228,1099},{799,1347},{657,1347},{228,1099},{156,976},{156,480},{228,357},{657,109},{799,109},{1228,357},{1300,480},{1371,769},{1085,1265},{1014,1306},{442,1306},{371,1265},{85,769},{85,687},{371,191},{442,150},{1014,150},{1085,191},{1371,687},{1371,852},{1157,1223},{942,1347},{514,1347},{299,1223},{85,852},{85,604},{299,233},{514,109},{942,109},{1157,233},{1371,604},{1300,1058},{728,1388},{156,1058},{156,398},{728,68},{1300,398}}; 
 
 int dump =0;
   
  for (int i = 0; i < 253; i++ ) {
    for (int j = -PMTrad[0]; j <= PMTrad[0]; j++ ) {
      int xo1 = data_top[i][0] + j;
      double ymax1 = sqrt(pow(PMTrad[0],2) - pow(j,2)) + 0.5;  //plus 1/2 to round properly during truncation in next line
      int y_max1 = (int) ymax1;
      for (int k = -y_max1; k <= y_max1; k++ ) {
        int yo1 = data_top[i][1] + k;
        if (sqrt( pow(xo1-728,2) + pow(yo1-728,2) ) > 728) {
          dump++;
        }        
        else { if ( sqrt( pow(j,2) + pow(k,2)) <= PMTrad[1] ) {
          PMTmap[0][xo1][yo1] = 1.;
          }
          else { PMTmap[0][xo1][yo1] = 0.5; }
        }
      }
    }
  }

  for (int i = 0; i < 241; i++) {
    for (int j = -PMTrad[0]; j <= PMTrad[0]; j++ ) {
      int xo2 = data_bot[i][0] + j;
      double ymax2 = sqrt(pow(PMTrad[0],2) - pow(j,2)) + 0.5;
      int y_max2 = (int) ymax2;
      for (int k = -y_max2; k <= y_max2; k++ ) { 
        int yo2 = data_bot[i][1] + k;
        if ( sqrt( pow(xo2-728,2) + pow(yo2-728,2)) > 728 ) {
          dump++;
        }
        else { if ( sqrt( pow(j,2) + pow(k,2)) <= PMTrad[1] )  {
          PMTmap[1][xo2][yo2] = 1.; 
          }
          else { PMTmap[1][xo2][yo2] = 0.5; }
        }
      }
    } 
  } 

}

