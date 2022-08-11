int init_KKdisk(ldouble r, ldouble th, ldouble *rho,ldouble *uint);


//MC240322--setup of KK00 disk, first HD version, then slowly moving towards
//a complete MHD version. Idea is to set it for NS and then shift to SMBH case.
ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,uT,uphi,uPhi;
ldouble rcyl, pres, eps2, coeff, lambda1, rhoc, rho0, pc, rd, alphav;
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
  rho0= RHO_DISC_MAX;
  rhoc = RHO_EPS * rho0;
  //notice that rhoc and coeff are multiplied with rho0, in KK00 rho0=1.
  
  coeff= (GAMMAM1/GAMMA/eps2)*(RINNER/r-(1.-GAMMA*eps2/GAMMAM1)*RINNER/rcyl);
  lambda1=(11./5.)/(1.+(64./25.)*alphav*alphav);
  pp[RHO] = rho0*pow(coeff,1./GAMMAM1);

  ldouble ucon[4];

  ucon[1] = 0;
  ucon[2] = 0;
  ucon[3] = 0;

  ucon[0] = sqrt(-1.0/geomBL.gg[0][0]);	

  pres=(1./RINNER)*eps2*pow(coeff,GAMMA/GAMMAM1);
       
        ucon[1] =  -(alphav/sin(th))*eps2*(10.-(32./3.)*lambda1*alphav*alphav-lambda1*(5.-1./(eps2*tan(th)*tan(th))))/sqrt(rcyl);//
	//ucon[1] = 0;       
	ucon[3] = (sqrt(1.-5.*eps2/2.)+(2./3.)*eps2*alphav*alphav*lambda1*(1.-6./(5.*eps2*tan(th)*tan(th))))/sqrt(rcyl)/r; //
	//ucon[3] = sqrt(1.0/pow(r,3.0));      
	fill_utinucon(ucon,geomBL.gg, geomBL.GG);

  	//ucon[0] = sqrt((-1.0-geomBL.gg[3][3]*ucon[3]*ucon[3])/geomBL.gg[0][0]);   


   pp[UU] = pres/GAMMAM1;

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

  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);

  #ifdef MAGNFIELD
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;

    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx/4.e-22,2.)-0.02,0.)*sqrt(1.e-23);
	  #ifdef QUADLOOPS
		  Acov[3]*=sin((M_PI/2.-geomBL.yy)/0.1);
	  #endif

	  #ifdef MULTIPLELOOPS
		  Acov[3]*=sin(geomBL.xx/3.);
	  #endif
    
     
    //Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/4.e-20,2.)-0.02,0.)*sqrt(1.e-23)*pow(sin(fabs(geomBL.yy)),4.);
    
    //printf("%lf\n",Acov[3]);
     

    pp[B1]=Acov[1];
    pp[B2]=Acov[2];
    pp[B3]=Acov[3];

    #endif
 
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


