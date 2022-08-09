//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,
//	ldouble *uu,ldouble *pp,int ifinit,int BCtype)

/**********************/
//geometries
ldouble rcyl, pres, eps, eps2, coeff, lambda, rhoc, pc, rd, alphav;
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv,iiya;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

/**********************/
ldouble r=geomBL.xx;
ldouble th=geomBL.yy;


eps=0.1;
eps2=eps*eps;
rcyl=r*sin(th);
rd=MKSR0;
alphav=1.0;
rhoc=0.01;

coeff=2./5./eps2*(1./r-(1.-5./2.*eps2)/rcyl);
lambda=11./5./(1.+64./25.*alphav*alphav);

//radius

//if(BCtype==XBCHI) //outflow in magn, atm in rad., atm. in HD
//ccm130422--keep the initial disk values at the rmax 
 if(ix>=NX) //analytical solution at rout only
//ccm220422--set the initial disk values at outer disk radius
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;

    //copying everything
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }
    
    //!! begin rescale
    //first transform to BL
    trans_pmhd_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);
    
    struct geometry geombdBL;
    fill_geometry_arb(iix,iiy,iiz,&geombdBL,BLCOORDS);
    ldouble rghost = geomBL.xx;
    ldouble rbound = geombdBL.xx;
    ldouble deltar = rbound - rghost;
    ldouble scale1 =  rbound*rbound/rghost/rghost;
    ldouble scale2 = rbound/rghost;

    pp[RHO]*=scale1;
    pp[UU] *=scale1;

//ccm220422--here starts the part Miki added to Debora's setup:    

// initial non-rotating adiabatic corona in hydrostatic equilibrium
 pp[RHO] = rhoc*pow(r,-3./2.);
 pres= 2./5.*rhoc*pow(r,-5./2.);
 pp[UU] = pres*pp[RHO]/(GAMMA-1);
 pc=pres;
 
//  pp[VX] = 0.0;
//  pp[VY] = 0.0;  
//  pp[VZ] = 0.0;
  
// Keplerian adiabatic disk in vertical pressure equilibrium with the
//   adiabatic corona, as given by Kluzniak & Kita (2000)

pres=eps2*pow(coeff,5./2.);
   
    if (pres >= pc && rcyl > rd)
      {pp[RHO] = pow(coeff,3./2.);      
       pp[VX] = -alphav/sin(th)*eps2*(10.-32./3.*lambda*alphav*alphav
       -lambda*(5.-1./(eps2*tan(th)*tan(th))))/sqrt(rcyl);
       pp[VZ] = (sqrt(1.-5./2.*eps2)+2./3.*eps2*alphav*alphav
       *lambda*(1.-6./(5.*eps2*tan(th)*tan(th))))/sqrt(rcyl); 
      }
    else
      {
       pres=pc;
      }  

 pp[UU]=GAMMA/GAMMAM1*pres/pp[RHO];

// Save conserved and primitives over domain + ghost (no corners)
// ccm--140422--all at t=1., could be t=0. also, why 1. here works and 0. not?
//  copyi_u(1.,u,upreexplicit); //conserved quantities before explicit update
//  copyi_u(1.,p,ppreexplicit); //primitive quantities before explicit update

/**/
//ccm--now we copy rho, vx and vz from T=1 into 2 ghostzones
// It should be done only in the disk part, but for now we take all.
// A loop in iiy should be done to limit only in the disk
//for(iiya=39; iiya<=65; iiya++){//41-63 if 10 above and below the middle cell
//   for(iv=0;iv<NV;iv++)
//      {
//      for(iiya = 0; iiya <= NY-1; iiy++){
//  copyi_u(1.,u,upreexplicit); //conserved quantities before explicit update
//  copyi_u(1.,p,ppreexplicit); //primitive quantities before explicit update

//    set_u(p,RHO,iix+1,iiy,iiz,get_u(ppreexplicit, 1, iix+1,iiya,iiz));
//    set_u(u,RHO,iix+1,iiya,iiz,get_u(upreexplicit, 1, iix+1,iiya,iiz));
    set_u(u,iv,iix,iiya,iiz,pp[RHO]);
    set_u(p,iv,iix,iiya,iiz,pp[RHO]);
//    set_u(p,VX,iix+1,iiya,iiz,get_u(ppreexplicit, 1, iix+1,iiya,iiz));
//    set_u(u,VX,iix+1,iiya,iiz,get_u(upreexplicit, 1, iix+1,iiya,iiz));
    set_u(u,iv,iix,iiya,iiz,pp[VX]);
    set_u(p,iv,iix,iiya,iiz,pp[VX]);
//    set_u(p,VZ,iix+1,iiya,iiz,get_u(ppreexplicit, 1, iix+1,iiya,iiz));
//    set_u(u,VZ,iix+1,iiya,iiz,get_u(upreexplicit, 1, iix+1,iiya,iiz));
    set_u(u,iv,iix,iiya,iiz,pp[VZ]);
    set_u(p,iv,iix,iiya,iiz,pp[VZ]);
/*
//    set_u(p,RHO,iix+2,iiya,iiz,get_u(ppreexplicit, 1, iix+2,iiya,iiz));
//    set_u(u,RHO,iix+2,iiya,iiz,get_u(upreexplicit, 1, iix+2,iiya,iiz));
     set_u(u,iv,iix,iiya,iiz,pp[RHO]);
    set_u(p,iv,iix,iiya,iiz,pp[RHO]);
//    set_u(p,VX,iix+2,iiya,iiz,get_u(ppreexplicit, 1, iix+2,iiya,iiz));
//    set_u(u,VX,iix+2,iiya,iiz,get_u(upreexplicit, 1, iix+2,iiya,iiz));
    set_u(u,iv,iix,iiya,iiz,pp[VX]);
    set_u(p,iv,iix,iiya,iiz,pp[VX]);
//    set_u(p,VZ,iix+2,iiya,iiz,get_u(ppreexplicit, 1, iix+2,iiya,iiz));
//    set_u(u,VZ,iix+2,iiya,iiz,get_u(upreexplicit, 1, iix+2,iiya,iiz));
    set_u(u,iv,iix,iiya,iiz,pp[VZ]);
    set_u(p,iv,iix,iiya,iiz,pp[VZ]);
*/
//}   
//}
/**/
    //transform back after rescaling
    trans_pmhd_coco(pp, pp,BLCOORDS, MYCOORDS, geomBL.xxvec,&geomBL, &geom);
    //!! end rescale


    //checking for the gas inflow
    ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
    conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
    if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
    if(ucon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in rho,uint and velocities and zero magn. field
	//set_hdatmosphere(pp,xxvec,gg,GG,4);
	//set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);

	//.!! begin reset to floor        
	//	pp[UU] = pp[RHO]*UURHORATIOMIN*3.;
	
	//!! end reset to floor

	ucon[1]=0.;
	/**
	#ifdef MAGNFIELD
	pp[B2]=pp[B3]=0.;
	#endif
	**/

	if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	pp[VX]=ucon[1];
	pp[VY]=ucon[2];
	pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
    
#ifdef RADIATION
    ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
    conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
    if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
    if(urfcon[1]<0.) //inflow, resseting to atmosphere
      {
	//!! begin reset to floor
	//pp[EE0] = ERADATMMIN;
	//!! end reset to floor

	//atmosphere in radiation
	//set_radatmosphere(pp,xxvec,gg,GG,0);
	urfcon[1]=0.;
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	pp[FX0]=urfcon[1];
	pp[FY0]=urfcon[2];
	pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
#endif

    
    p2u(pp,uu,&geom);
    return 0;   
  }
 else if(BCtype==XBCLO) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

      for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }

      if(RMIN>rhorizonBL)
	{
	  //checking for the gas inflow
	  ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
	  conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	  trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
	  if(ucon[1]<0.) //inflow, resseting to atmosphere
	    {
	      //atmosphere in rho,uint and velocities and zero magn. field
	      //set_hdatmosphere(pp,xxvec,gg,GG,4);
	      ucon[1]=0.;
	      trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	      conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	      pp[VX]=ucon[1];
	      pp[VY]=ucon[2];
	      pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
	    }

#ifdef RADIATION
	  ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
	  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
	  trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
	  if(urfcon[1]<0.) //inflow, resseting to atmosphere
	    {
	      //atmosphere in radiation
	      //set_radatmosphere(pp,xxvec,gg,GG,0);
	      urfcon[1]=0.;
	      trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	      conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	      pp[FX0]=urfcon[1];
	      pp[FY0]=urfcon[2];
	      pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
	    }
#endif
	}

     p2u(pp,uu,&geom);
     return 0;
   }

//reflections/outflow in theta 
//in 3D polar cells overwritten with #define CORRECT_POLARAXIS_3D
if(BCtype==YBCLO) //upper spin axis 
  {      
    
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//v_theta
	if(iv==VY || iv==B2 || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
     

    p2u(pp,uu,&geom);
    return 0;
  }

if(BCtype==YBCHI) //lower spin axis
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
  	  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==VY || iv==B2 || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	
      }
 
    p2u(pp,uu,&geom); 
    return 0; 
  }
   
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

