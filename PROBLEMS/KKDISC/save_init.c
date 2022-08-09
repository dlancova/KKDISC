/*MC271021--thin disc problem, setup with Kluzniak&Kita(2000) HD disc */
int init_dsandvels_KKHDthindisc(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

ldouble pre, R, L, uint;
ldouble eps, eps2, coeff, lambda, rhoc, pca, rd, alphav;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

init_dsandvels_KKHDthindisc(r, th, BHSPIN, &rho, &uint, &ell); 
uintorg=uint;

if(rho<0.) //outside donut
{
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif

#ifdef MAGNFIELD
#ifdef BATMZ //put field in the atmosphere and torus (Better for MAD with small torus?)
    if(r*sin(th) < BATMZ_MAXR && r*sin(th) > BATMZ_MINR)
    {
      pp[B1]=0.;//BATMZ_B0*cos(th)/r/r; //will this get auto rescaled?
      pp[B2]=0.;//BATMZ_B0*sin(th)/r/r;  
      pp[B3]=0.;
    }
#endif
#endif
 }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

    pgas = GAMMAM1 * uint;
    ell*=-1.;


    ldouble ult,ulph,ucov[4],ucon[4];
    ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
         //  (sqrt(1.-5./2.*eps2)+2./3.*eps2*alphav*alphav
      // *lambda*(1.-6./(5.*eps2*tan(th)*tan(th))))/sqrt(R);

    ult = ulph / ell;

    ucov[0]=ult;
    ucov[1]=0.;
    ucov[2]=0.;
    ucov[3]=ulph;

    ucov[0]=ult;
    ucov[1]=0.;//-alphav/sin(th)*eps2*(10.-32./3.*lambda*alphav*alphav
       //-lambda*(5.-1./(eps2*tan(th)*tan(th))))/sqrt(R);
    ucov[2]=0.;
    ucov[3]=ulph;
    
    indices_12(ucov,ucon,geomBL.GG);

    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
   
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
//#if (NTORUS==2)
    //standard single poloidal loop
    ldouble rin=KT_R0;
    ldouble rout = KT_ROUT;//5.e3; //cover full range of torus?
    ldouble lambda = 2.*(rout-rin);

      Acov[3] = my_max( pow(r, RADIUS_POWER) * (pp[RHO] - RHO_CUT_FACTOR * KT_RHO0) * pow(sin(th), 3) * pow(sin(r/lambda), SIN_POWER) , 0.);
//#endif
   
   
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
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
