
int init_KKdisk(ldouble r, ldouble th, ldouble *rhoout,ldouble *uintout){
  ldouble eps2=HR_INIT*HR_INIT;
  ldouble rcyl=r*sin(th);
  ldouble rd=RINNER;
  ldouble alphav=ALPHA_DISC;
  ldouble rho0=RHO_DISC_MAX;//free param. multipl. factor to KK00=1 max disk density
  ldouble rhoc=RHO_EPS * RHO_DISC_MAX;

  ldouble coeff=(GAMMAM1/GAMMA/eps2)*(RINNER/r-(1.-GAMMA*eps2/GAMMAM1)*RINNER/rcyl);
  ldouble lambda1=11./5./(1.+64./25.*alphav*alphav);

  ldouble rho = rhoc*pow(r/RINNER,-1./GAMMAM1);
  ldouble pres= (GAMMAM1/GAMMA)*rhoc*pow(r/RINNER,-GAMMA/GAMMAM1)/RINNER;
  ldouble uint = pres/(GAMMA-1);
  ldouble pc=pres;

  pres=(1./RINNER)*eps2*rho0*pow(coeff,GAMMA/GAMMAM1);

  if (pres >= pc && rcyl > rd){
      rho = RHO_DISC_MAX*pow(coeff,1./GAMMAM1);     
  }
  else{
    rho = -1.0;
  }  

  *rhoout = rho;
  *uintout = pres/GAMMAM1;

}

