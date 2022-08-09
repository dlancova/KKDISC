
int init_dsandvels_KKdisk(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell)
{
  ldouble rho,pre,phi,phi0,uint,phit0,phit,L,L0,cs0,uint0,rho0,n,R;
  ldouble eps, eps2, coeff, lambda, rhoc, pca, rd, alphav;
  
//  if(r<KT_RMIN) {*rhoout=-1.;return 0;}

  R=r*sin(th); //cylindrical radius 
  eps=0.1;
  eps2=eps*eps;
  rd=1.7;
  alphav=1.0;
  rhoc=0.01;
  
  coeff=2./5./eps2*(1./r-(1.-5./2.*eps2)/R);
  lambda=11./5./(1.+64./25.*alphav*alphav); 

//ccm--271021
/* initial non-rotating adiabatic corona in hydrostatic equilibrium  */
 rho = rhoc*pow(r,-3./2.);
 pre = 2./5.*rhoc*pow(r,-5./2.);
 
 pca=pre;

/* Keplerian adiabatic disk in vertical pressure equilibrium with the
   adiabatic corona, as given by Kluzniak & Kita (2000) */

  pre=eps2*pow(coeff,5./2.);

    if (pre >= pca && R > rd)
      {rho = pow(coeff,3./2.);
      }
    else
      {
       pre=2./5.*rhoc*pow(r,-5./2.);
      }  
//ccm
  uint = pre / GAMMAM1;

/*
  #ifdef KT_RHO_FLOOR //This keeps the same Bernoulli number but reduces mass in disk
  if (rho < KT_RHO_FLOOR) {*rhoout=-1.;return 0;}
  #endif
*/
  *rhoout = rho;
  *uuout = uint;
  *ell = L;

  return 0;

}

