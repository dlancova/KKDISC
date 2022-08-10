int init_KKdisk(ldouble r, ldouble th, ldouble *rho,ldouble *uint);


//MC240322--setup of KK00 disk, first HD version, then slowly moving towards
//a complete MHD version. Idea is to set it for NS and then shift to SMBH case.
ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,uT,uphi,uPhi;
ldouble rcyl, pres, eps2, coeff, lambda1, rhoc, pc, rd, alphav;
ldouble mub, rminb, mmb;
int loops, openn;//to be given by hand later, for choice of loops model
//geometries

//geometry in MYCOORDS
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

//geometry in BLCOORDS
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

//r and theta in BL
ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

init_KKdisk(r, th, &rho, &uint);

if (rho<0){
  // for now go with the 0 atm. 
  set_hdatmosphere(pp,geom.xxvec, geom.gg, geom.GG,5);
}
else{
  set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,5);
  eps2=HR_INIT*HR_INIT;
  rcyl=r*sin(th);
  rd=RINNER;
  alphav=ALPHA_DISC;
  rhoc=RHO_EPS * RHO_DISC_MAX;
  //notice that rhoc and coeff are multiplied with rho0, in KK00 rho0=1.
  coeff=RHO_DISC_MAX*2./5./eps2*(1./r-(1.-5./2.*eps2)/rcyl);
  lambda1=11./5./(1.+64./25.*alphav*alphav);
  pp[RHO] = pow(coeff,3./2.);
  pres= 2./5.*rhoc*pow(r,-5./2.);
  pp[UU] = pres*pp[RHO]/(GAMMA-1);
  pc=pres;
  ldouble ucon[4];

  ucon[1] = 0;
  ucon[2] = 0;
  ucon[3] = 0;

  ucon[0] = sqrt(-1.0/geomBL.gg[0][0]);	


  pres=eps2*pow(coeff,5./2.);
       
  ucon[1] = -alphav/sin(th)*eps2*(10.-32./3.*lambda1*alphav*alphav       -lambda1*(5.-1./(eps2*tan(th)*tan(th))))/sqrt(rcyl*pow(sin(th),2.0));
//	ucon[1] = 0;       
	ucon[3] = 5.0e-3 * (sqrt(1.-5./2.*eps2)+2./3.*eps2*alphav*alphav *lambda1*(1.-6./(5.*eps2*tan(th)*tan(th))))/sqrt(rcyl)/r; 
	//ucon[3] = sqrt(1.0/pow(r,3.0));      
fill_utinucon(ucon,geomBL.gg, geomBL.GG);
	//ucon[0] = sqrt(-1.0/(geomBL.gg[0][0] + geomBL.gg[3][3]*ucon[3]*ucon[3] + geomBL.gg[1][1]*ucon[1]*ucon[1]));
  	//ucon[0] = sqrt((-1.0-geomBL.gg[3][3]*ucon[3]*ucon[3])/geomBL.gg[0][0]);   


 pp[UU]=GAMMA/GAMMAM1*pres/pp[RHO];

ucon[1]*= ucon[0];
ucon[2]*= ucon[0];
ucon[3]*= ucon[0];

fill_utinucon(ucon,geomBL.gg, geomBL.GG);
 //ucon[0]=1./sqrt(-geomBL.GG[0][0])i
 conv_vels(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG); 

//printf("%lf %lf %lf %lf\n",ucon[0],ucon[1],ucon[2],ucon[3]);

 pp[VX] = ucon[1];
 pp[VY] = ucon[2];
 pp[VZ] = ucon[3];

//----finish of the HD setup by KK00-------- 

#ifdef MAGNFIELD
//MC090522-choose here by hand the shape of initial field:
loops=0;//0 for hourglass shape, 1 for various poloidal loops, 2 for quadrupole
openn=1;//0 for open lines, 1 for quadrupole
if (loops == 0) {//beginning of "loops" loop
if(openn == 0){
 //large scale hourglass field Zhu&Stone (2018), Mishra et al.(2019)
 //mmb BMishra -5/4, Zhu&Stone -9/4
mub=1.;//free parameter for strength of mag.field
rminb=1.*rd;//free parameter for reach of inner B
mmb=-5./4.;//parameter for shape of B.
     if(rcyl < rminb){
      pp[B1]=mub*cos(th)*pow(rminb,mmb)*(1.+sin(th));
      pp[B2]=-mub*sin(th)*pow(rminb,mmb);
      pp[B3]=0.;  
    }else{
      pp[B1]=mub*pow(r*sin(th),mmb)*cos(th)*(1.+sin(th));
      pp[B2]=-mub*pow(r*sin(th),mmb)*sin(th);
      pp[B3]=0.;}
    }else{
//quadrupole
    if(rcyl < rminb){
      pp[B1]=0;
      pp[B2]=0;
      pp[B3]=0.;  
    }else{
      pp[B1]=3.0/2.0*mub*(3.0*cos(th)*cos(th)-1.0)/(r*r*r*r);
      pp[B2]=3.0*mub*cos(th)*sin(th)/(r*r*r*r);
      pp[B3]=0.;
    }
}
/* quadrupole 
  B0[0] = 3.0/2.0*g_inputParam[MU]*(3.0*cos(x2)*cos(x2)-1.0)/(x1*x1*x1*x1);
  B0[1] = 3.0*g_inputParam[MU]*cos(x2)*sin(x2)/(x1*x1*x1*x1);
  B0[2] = 0.0;
*/
    
//    
//--with poloidal field loops:    
    }else{//vector potential to calculate Bp loops
    if(loops == 1) {
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;

#ifdef NONELOOP//no loop, some extended loop in fact
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*sqrt(geomBL.xx)/1.e-5,2.)-0.0001,0.)*
      pow(sin(fabs(geomBL.yy)),4.);
#endif

#ifdef SINGLELOOP//standard single poloidal loop
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/(4.*rho0),2.)-0.02,0.)*sqrt(1.e-23);//*pow(sin(fabs(geomBL.yy)),4.);
#endif

#ifdef QUADLOOPS//mirrored quadrupole loops, one above, one below disk equator
//    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/(4.*rho0),2.)-0.002,0.)/sqrt(rho0*1.e3);    
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/(3000.),GAMMA)-0.02,0.)/sqrt(rho0)/rho0/1000.;    
    Acov[3]*=sin(geomBL.xx/5.); 
//   Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*sqrt(geomBL.xx)/1.e-5,2.)-0.0001,0.)
#endif

#ifdef MULTIPLELOOPS//sane-like loops
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/(3000.),GAMMA)-0.02,0.)/pow(rho0,1.5)/1000.;
    //Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/(4.*rho0),2.)-0.02,0.)*sqrt(rho0*1.e1);    
    //Acov[3]*=sin((M_PI/2.-geomBL.yy)/0.1);//sin(geomBL.xx/3.);
#endif

if(openn == 1){
}else{
    pp[B1]=0.;
    pp[B2]=0.;
    pp[B3]=Acov[3];
}

}//end of "loops=1" sub-loop

/*
else{//other types of loops 
pp[B1]=pp[B2]=pp[B3]=0.; 
   ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;

#if(NTORUS==3)
    //LIMOFIELD from a=0 SANE harm init.c
    ldouble lambda = 2.5;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 350.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vector potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.5*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }

    Acov[3]=vpot;

    //if(iy==TNY/2) {printf("%d %d %e %e %e %e %e %e\n",ix,iy,vpot,r,q,fr,uchop,uchopmid);getch();}

#elif (NTORUS==4)
    //LIMOFIELD from a=0 SANE harm init.c + denser loops
    ldouble lambda = 1.;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 550.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.5*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }

    Acov[3]=vpot;

#elif (NTORUS==5)

   //LIMOFIELD from a=0 SANE harm init.c for mimic_dynamo
    ldouble lambda = 2.5;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 550.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.5*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
        
    Acov[2]=vpot*sin((M_PI/2.-geomBL.yy));;

    

#elif (NTORUS==7) 

   //quadrupolar
    ldouble lambda = 1.5;//2.5;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 350.; //outer boundary of field loops
    ldouble u_av = uintorg;
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 1.25*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    Acov[3]=vpot*sin((M_PI/2.-geomBL.yy));;

#elif (NTORUS==17) 

    //single loop
    ldouble lambda = 15.;//2.5;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 350.; //outer boundary of field loops
    ldouble u_av = uintorg;
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 1.25*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    Acov[3]=vpot;

#elif (NTORUS==6)

   //LIMOFIELD from a=0 SANE harm init.c with flipping polarity in theta
    ldouble lambda = 2.5;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 350.; //outer boundary of field loops
    ldouble u_av = uintorg;
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 1.25*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    Acov[3]=vpot*sin((M_PI/2.-geomBL.yy));;

#elif (NTORUS==77 || NTORUS==78)


  
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*sqrt(geomBL.xx)/1.e-5,2.)-0.0001,0.)*
      pow(sin(fabs(geomBL.yy)),4.);

#elif (NTORUS==79) //a'la adaf paper - center too close
  
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*sqrt(geomBL.xx)/1.e-5,2.)-0.01,0.)*
      pow(sin(fabs(geomBL.yy)),4.);

#elif (NTORUS==80) //a'la adaf paper but ~ RHO
  
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/1.e-5,2.)-0.1,0.)*
      pow(sin(fabs(geomBL.yy)),4.);


#elif (NTORUS==81) //a'la adaf paper but ~ UU

  //LIMOFIELD from a=0 MAD harm init.c
    ldouble lambda = 25.;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 800.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.5*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 8); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }

    Acov[3]=vpot;

#else //standard single poloidal loop

    
  
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*sqrt(geomBL.xx)/1.e-5,2.)-0.0001,0.)*
      pow(sin(fabs(geomBL.yy)),4.);
				     //*step_function(-(geomBL.xx-350.),10.);
#endif

    pp[B1]=Acov[1];
    pp[B2]=Acov[2];
    pp[B3]=Acov[3];
#endif
}else{//quadrupole loops
p[B1]=pp[B2]=pp[B3]=0.; 
   ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;
*/
}//end of "loops" if-then loop

#ifdef PERTMAGN //perturb to break axisymmetry
//pp[UU]*=1.+((double)rand()/(double)RAND_MAX-0.5)*2.*PERTMAGN;
pp[UU]*=1.+PERTMAGN*sin(10.*2.*M_PI*(MAXZ-geomBL.zz)/(MAXZ-MINZ));
#endif

#endif //end of MAGNFIELD

  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
}

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1],ix,iy,iz);
//to conserved
p2u(pp,uu,&geom);

/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy_cell(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);


