
int init_dsandvels_KKHDthindisc(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell)
{
//  ldouble rho,pre,phi,phi0,uint,phit0,phit,L,L0,cs0,uint0,rho0,n,R;
  ldouble rho, pre, R, L, uint;
  ldouble eps, eps2, coeff, lambda, rhoc, pca, rd, alphav;

  R=r*sin(th); //cylindrical radius 
  eps=0.1;
  eps2=eps*eps;
  rd=1.7;
  alphav=1.0;
  rhoc=0.01;
  
   coeff=RHO_DISC_MAX*2./5./eps2*(RINNER/r-(1.-5./2.*eps2)*RINNER/rcyl);
   lambda=11./5./(1.+64./25.*alphav*alphav); 
/* initial non-rotating adiabatic corona in hydrostatic equilibrium  */
 rho = rhoc*pow(r,-3./2.);
 pre = 2./5.*rhoc*pow(r,-5./2.);
 
 pca=pre;

/* 
  v[VX1] = 0.0;
  v[VX2] = 0.0;  
  v[VX3] = 0.0;
*/  
/* Keplerian adiabatic disk in vertical pressure equilibrium with the
   adiabatic corona, as given by Kluzniak & Kita (2000) */
   
  pre=eps2*pow(coeff,5./2.);
  
    if (pre >= pca && R > rd)
      {rho = pow(coeff,3./2.);
//       ucov[1] = -alphav/sin(th)*eps2*(10.-32./3.*lambda*alphav*alphav
 //      -lambda*(5.-1./(eps2*tan(th)*tan(th))))/sqrt(R);
 //      ucov[3] = (sqrt(1.-5./2.*eps2)+2./3.*eps2*alphav*alphav
//       *lambda*(1.-6./(5.*eps2*tan(th)*tan(th))))/sqrt(R); 
      }
    else
      {
       pre=2./5.*rhoc*pow(r,-5./2.);
      }  

  *rhoout = rho;
  *uuout = uint;
  *ell = L;

  return 0;

}

