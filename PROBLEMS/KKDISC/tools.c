
int init_KKdisk(ldouble r, ldouble th, ldouble *rhoout,ldouble *uintout){
  ldouble eps2=HR_INIT*HR_INIT;
  ldouble rcyl=r*sin(th);
  ldouble rd=RINNER;
  ldouble alphav=ALPHA_DISC;
  ldouble rho0=RHO_DISC_MAX;//free param. multipl. factor to KK00=1 max disk density
  ldouble rhoc=RHO_EPS * RHO_DISC_MAX;

  ldouble coeff=RHO_DISC_MAX*2./5./eps2*(RINNER/r-(1.-5./2.*eps2)*RINNER/rcyl);
  ldouble lambda1=11./5./(1.+64./25.*alphav*alphav);

  ldouble rho = rhoc*pow(r,-3./2.);
  ldouble pres= 2./5.*rhoc*pow(r,-5./2.);
  ldouble uint = pres/(GAMMA-1);
  ldouble pc=pres;

  pres=(1./RINNER)*eps2*pow(coeff,5./2.);

  if (pres >= pc && rcyl > rd){
      rho = pow(coeff,3./2.);     
  }
  else{
    rho = -1.0;
  }  

  *rhoout = rho;
  *uintout = pres/GAMMAM1;

}

