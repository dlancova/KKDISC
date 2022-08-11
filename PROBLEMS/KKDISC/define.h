
/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/ 
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN

//#define MULTIPLELOOPS
#define HLOOPS 0.5
#define MAXBETA 1e-3 //target pmag/pgas int the midplane

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
//blackhole
/************************************/
#define MASS 10.//MSUNCM
#define BHSPIN 0.

/************************************/
//physics 
/************************************/
#define GAMMA (5./3.)


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
#define RMIN 1.85
#define RMAX 50

/************************************/
//setup for MKS2 coords
/************************************/
#ifdef myMKS2COORDS 
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.1)  // 0.001 for bigger domin closer to polar axis 
#define MAXY (1.-MINY)
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

/************************************/
//resolution 
/************************************/

#define TNX 124 //16*17
#define TNY 124 //16*12
#define TNZ 1//32//128 //16*8

//number of tiles for MPI
#define NTX 1//32//16//32//16//17
#define NTY 1//16//16//32//8//12
#define NTZ 1//2//8

/************************************/
//boundary conditions 
/************************************/
#define SPECIFIC_BC //in bc.c
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define DTOUT1 1. //res - files
#define DTOUT2 1000.  //avg - files

#define OUTCOORDS BLCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define SILOOUTPUT 1
#define SIMOUTPUT 1
#define RADOUTPUT 0 
#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define THOUTPUT 0
#define THPROFRAD1US 30

#if(TNZ==1)//if 2D
#define SILO2D_XZPLANE
#else
#define FULLPHI
#endif

/************************************/
//initial disk 
/************************************/
/**/ //This for later choices
#define TDISK 1

/************************************/
//disc choices
/************************************/

#if(TDISK==1) //Kluzniak&Kita(2000) thin disk in HD
    // Inner edge of the disc
    #define RINNER r_ISCO_BL(BHSPIN)
    // Max density in the disc centre
    #define RHO_DISC_MAX 1.0
    // max atm. density = RHO_EPS * RHO_DISC_MAX (at horizon)
    #define RHO_EPS 1.0e-4
    #define EPSS 0.1//thin disk height ratio
    #define ALPHA_DISC 0.5 //viscous alpha
    #define HR_INIT 0.1 //initial disc thickness 
#endif //TDISK 1


#if(TDISK==2) //Kluzniak&Kita(2000) thin disk with hourglass Bfield

    #define EPSS 0.1//thin disk height ratio

#endif//ending TDISK

#if(TDISK==3) //Kluzniak&Kita(2000) thin disk in RMHD
    // TBD
#endif
//#define BETANORMFULL

/************************************/
//rmhd floors
/************************************/
#define RHOFLOOR 1.e-50
//#define UURHORATIOMIN 1.e-8//1K: uu/rho = 7.259162e+12
//#define UURHORATIOMAX 1.e20
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
