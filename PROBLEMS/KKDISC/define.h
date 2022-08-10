//#define NONRELMHD
//#define MHD
/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/ 
#define RESTART
#define RESTARTGENERALINDICES
//#define RESCALEDENSITY 2.
#define RESTARTNUM -1

/************************************/
//magnetic choices
/************************************/
/**/
#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
//from NSDISK
/*
#define MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN
#define DYNAMORADIUS 15.
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28
#define MAXBETA .01 //target pmag/pgas int the midplane
*/
/**/
/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM 0.5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//YSO
/************************************/
#define MASS 10.//MSUNCM
#define BHSPIN 0.

/************************************/
//blackhole
/************************************/
//#define MASS 1.e1//6
//#define BHSPIN 0.10

/************************************/
//coordinates
/************************************/
//#define myMKS3COORDS  // use for 3D
#define myMKS2COORDS // good for 2D
#define MKSR0 0.
#define MKSH0 0.6
#define MKSMY1 0.0025
#define MKSMY2 0.025
#define MKSMP0 1.5
#define METRICAXISYMMETRIC
#define RH 0.6//1.348
#define RMIN 2.0 //1.85
#define RMAX 50
//1.e5
/**/ //setup in MKS2 coords
#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS//6 for spherical
#define MINX (log(RMIN-MKSR0))//(RMIN)
#define MAXX (log(RMAX-MKSR0))//(RMAX)
#define MINY (0.1)
#define MAXY (1.-MINY)//M_PI-MINY)//(1.-MINY)//M_PI-MINY)
#endif
#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS2COORDS
#define MINX (RMIN)
#define MAXX (RMAX)
//#define Y_OFFSET 0.009
#define MINY (0.)
#define MAXY (M_PI-MINY)
#endif

#define PHIWEDGE (M_PI/2.)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)
/**/ //end of setup in MKS2 coords

/* //setup in spherical coords
//#define MKS2COORDS
//#define mySPHCOORDS
#define MYCOORDS 6
//#define MYCOORDS MKS2COORDS 
#define METRICAXISYMMETRIC

//#define MINX 0.
//#define MAXX 1.
#define MINX (RMIN)
#define MAXX (RMAX)

//#define Y_OFFSET 0.009
#define MINY (0.1)
#define MAXY (M_PI-MINY)

//#define MKSMY1 0.0025
//#define MKSMY2 0.025
//#define MKSMP0 1.2

//#endif //End of whateverCOORDS
#define PHIWEDGE (0.5*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)
*/ //end of setup in spherical coords
/************************************/
//resolution 
/************************************/
//total resolution
//#define TNX 384 //16*17
//#define TNY 384 //16*12
#define TNX 124 //16*17
#define TNY 124 //16*12
#define TNZ 1//32//128 //16*8

//number of tiles
#define NTX 1//32//16//32//16//17
#define NTY 1//16//16//32//8//12
//#define NTX 8//16//32//16//17
//#define NTY 8//16//32//8//12
#define NTZ 1//2//8

/************************************/
//boundary conditions 
/************************************/
#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define DTOUT1 1. //res - files
//#define DTOUT2 1000.  //avg - files
//#define DTOUT3 1. //box,var - files
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
//#define OUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000//5000
#define SILOOUTPUT 1
#define OUTOUTPUT 0

#if(TNZ==1)//if 2D
#define SILO2D_XZPLANE
#else
#define FULLPHI
#endif

// Inner edge of the disc
#define RINNER 6.
// Max density in the disc centre
#define RHO_DISC_MAX 1.0
// max atm. density = RHO_EPS * RHO_DISC_MAX (at horizon)
#define RHO_EPS 1.0e-2

/************************************/
//physics 
/************************************/
#define GAMMA (5./3.)

/************************************/
//initial disk 
/************************************/
/**/ //This for later choices
#define TDISK 2

#if(TDISK==1) //Kluzniak&Kita(2000) thin disk in HD
#endif

#if(TDISK==2) //Kluzniak&Kita(2000) thin disk with hourglass Bfield
//in TDISK
#define EPSS 0.1//thin disk height ratio

#ifdef RADIATION
//#define RGAMMA 4./3.
//#else
//#define RGAMMA 5./3.
#endif

#endif//ending TDISK

#if(TDISK==3) //Kluzniak&Kita(2000) thin disk in RMHD
#endif
/**/
//#define BETANORMFULL

#define RHO0 1

/************************************/
//rmhd floors
/************************************/
/**/
#define RHOFLOOR 1.e-50
//#define UURHORATIOMIN 1.e-8//1K: uu/rho = 7.259162e+12
//#define UURHORATIOMAX 1.e20
//NSDISK
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.
