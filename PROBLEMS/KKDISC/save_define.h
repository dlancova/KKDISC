
/************************************/
//general
/************************************/
//#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/ 
#define RESTART
//#define RESTARTGENERALINDICES
//#define RESCALEDENSITY 2.
#define RESTARTNUM -1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
//#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM 0.9
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.0

/************************************/
//coordinates 
/************************************/
#define MKS3COORDS

#define METRICAXISYMMETRIC

#define RH 2.//1.348
#define RMIN 1.3
#define RMAX 50

#define MYCOORDS BLCOORDS

#define ROUT 1000.
#define MKSR0 0.
#define MKSH0 0.6
#define MKSMY1 0.0025
#define MKSMY2 0.025
#define MKSMP0 1.2

#define MINX (log(RMIN-MKSR0))
#define MAXX (log(ROUT-MKSR0))
#define MINY 0.
#define MAXY 1.

#define PHIWEDGE (0.5*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

/************************************/
//resolution 
/************************************/
//total resolution
#define TNX 128//28*9
#define TNY 64//128//26*9
#define TNZ 1//32//32 //2*8
//number of tiles
#define NTX 1
#define NTY 1
#define NTZ 1//4

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
#define DTOUT1 100. //res - files
//#define DTOUT2 1000.  //avg - files
//#define DTOUT3 1. //box,var - files
#define OUTCOORDS BLCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 10.e10
#define NOUTSTOP 5000//5000
#define SILOOUTPUT 1
#define SIMOUTPUT 2
#define COORDOUTPUT 0
/*
#if(TNZ==1)
#define SILO2D_XZPLANE
#else
#define FULLPHI
#endif
*/
/************************************/
//physics 
/************************************/
#define GAMMA (5./3.)//(4./3.)
/************************************/
//polar axis
/************************************/
#ifndef myCYLCOORDS
#define CORRECT_POLARAXIS
#endif
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2
