
int init_KKdisk(ldouble r, ldouble th, ldouble *rhoout,ldouble *uintout){
  ldouble eps2=EPSS*EPSS;
  ldouble rcyl=r*sin(th);
  ldouble rd=RINNER;
  ldouble alphav=1.0;
  ldouble rho0=RHO0;//free param. multipl. factor to KK00=1 max disk density
  ldouble rhoc=0.01*rho0;

  ldouble coeff=rho0*2./5./eps2*(1./r-(1.-5./2.*eps2)/rcyl);
  ldouble lambda1=11./5./(1.+64./25.*alphav*alphav);

  ldouble rho = rhoc*pow(r,-3./2.);
  ldouble pres= 2./5.*rhoc*pow(r,-5./2.);
  ldouble uint = pres*rho/(GAMMA-1);
  ldouble pc=pres;

  pres=eps2*pow(coeff,5./2.);

  if (pres >= pc && rcyl > rd){
      rho = pow(coeff,3./2.);     
  }
  else{
    rho = -1.0;
  }  

  *rhoout = rho;
  *uintout = GAMMA/GAMMAM1*pres/rho;

}

