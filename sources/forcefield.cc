#define MINCUT 2.232 // minimum 2D LJ cutoff

function_t u,ufull,f,ffull;

char *FFinfo(int verbose) /****************************************** FFinfo */
{
  static char info[64];

  if (verbose) switch (ss.ff) {
      case LJ: sprintf(info,"2D Lennard-Jones cutoff c=%g",ss.C2); break;
      case WCA: sprintf(info,"2D WCALJ"); break;
      case PD: sprintf(info,"Penetrable disks a=%g c=%g",ss.a,ss.C2); break;
      case DW: sprintf(info,"Double Well a=%g b=%g c=%g",ss.a,ss.b,ss.C2); break;
      default: strcpy(info,"ERROR"); }
  else switch (ss.ff) {
      case LJ: sprintf(info,"LJ c=%g",ss.C2); break;
      case WCA: sprintf(info,"WCALJ"); break;
      case PD: sprintf(info,"PD a=%g c=%g",ss.a,ss.C2); break;
      case DW: sprintf(info,"DW%g,%g,%g",ss.a,ss.b,ss.C2); break;
      default: strcpy(info,"ERROR"); }

  return info;
}

/*
  the potentials and forces as functions of rr=r^2
  u=u(rr) is the pair energy at distance r
  f(rr) = -2 d u(rr)/d rr = pair_force(r)/r
*/

/* 2D Lennard-Jones */
inline
double uLJfull(double rr) /***************************************** uLJfull */
{
  rr=1/Sqr(rr);
  return 4*(Sqr(rr)-rr);
}

inline
double fLJfull(double rr) /***************************************** fLJfull */
{
  double r4=Sqr(rr);
  return (32/r4-16)/(rr*r4);
}

inline
double uLJ(double rr) /***************************************************** u */
{
  if (rr>ss.C2q)
    return 0;
  else if (rr<ss.C1q) {
    rr=1/Sqr(rr);
    return 4*(Sqr(rr)-rr); }
  else {
    return ss.A*Sqr(rr-ss.C2q); }
}

inline
double fLJ(double rr) /***************************************************** f */
{
  if (rr>ss.C2q)
    return 0;
  else if (rr<ss.C1q) {
    double r4=Sqr(rr);
    return (32/r4-16)/(rr*r4); }
  else {
    return ss.A4*(rr-ss.C2q); }
}


/* 2D WCALJ */
inline
double uWCALJ(double rr) /******************************************* uWCALJ */
{
  if (rr>1.414213562373095)
    return 0;
  else {
    rr=1/Sqr(rr);
    return 4*(Sqr(rr)-rr)+1; }
}

inline
double fWCALJ(double rr) /******************************************* fWCALJ */
{
  if (rr>1.414213562373095)
    return 0;
  else {
    double r4=Sqr(rr);
    return (32/r4-16)/(rr*r4); }
}


/* penetrable disks */
inline
double uPD(double rr) /***************************************************** u */
{
  if (rr>ss.C2q)
    return 0;
  else
    return ss.A*Sqr(rr-ss.C2q);
}

inline
double fPD(double rr) /***************************************************** f */
{
  if (rr>ss.C2q)
    return 0;
  else
    return ss.A4*(ss.C2q-rr);
}


/* double well, for Penrose */
inline
double uDW(double rr) /************************************************* uDW */
{
  if (rr>ss.C2q)
    return 0;
  else
    return ss.A*(1-rr)*(ss.aq-rr)*(ss.bq-rr)*Sqr(ss.C2q/rr-1);
}

double fDW(double rr) /************************************************* fDW */
{
  if (rr>ss.C2q)
    return 0;
  else
    return (((ss.P[0]/rr+ss.P[1])/rr+ss.P[2])/rr + ss.P[3] + ss.P[4]*rr)*(ss.C2q-rr);
}


/* smooth cutoff setup */
void setss(ff_e ff,double cutoff) /************************************** setss */
{
  int n;
  double oldC1;

  ss.ff=ff;
  sliders.c->value(cutoff);
  ss.C2=cutoff; // also denoted as c
  ss.C2q=Sqr(ss.C2);
  ss.aq=Sqr(ss.a);
  ss.bq=Sqr(ss.b);
  //  fprintf(stderr,"setss: ff=%d cutoff=%g\n",(int)ff,cutoff);

  switch (ff) {
    case LJ: 

      
    again:
      if (cutoff<MINCUT) {
        fl_alert("\
Cannot determine smooth cutoff,\n\
the minimum c=%g will be used.",(double)MINCUT);
        ss.C2=MINCUT; 
        ss.C2q=Sqr(ss.C2); }

      /* 2DLJ smoothing */
      ss.C2q=ss.C1q=1e99; // in case of problems - will be changed
      u=uLJ; f=fLJ;
      ufull=uLJfull; ffull=fLJfull;
      
      ss.C1=0.7*ss.C2;
      
      n=0;
      do {
        double den=f(Sqr(ss.C1))*ss.C1;
        
        oldC1=ss.C1;
        ss.C1=1+4*ufull(Sqr(ss.C1))/ss.C1/den;
        if (fabs(den)<1e-50 || n++>2000 || ss.C1<0) {
          // should never be activated ..
          cutoff=1;
          goto again; }
        
        ss.C1=ss.C2/sqrt(ss.C1);
      } while (fabs(1-ss.C1/oldC1) > 1e-15);
      
      if (ss.C1>0.999*ss.C2 || ss.C1<0.5*ss.C2) {
        // should never be activated ..
        cutoff=1; // will cause error
        goto again; }

      sliders.c->value(ss.C2);
      ss.C1q=Sqr(ss.C1);
      ss.C2q=Sqr(ss.C2);
      ss.A=ufull(ss.C1q)/Sqr(ss.C2q-ss.C1q);
      // MACSIMUS: ss->A4= -SS_F(C1)/(C2*C2-C1*C1)/C1;
      ss.A4=-ffull(ss.C1q)/(ss.C2q-ss.C1q);
      break;
    
      case DW:  // Penrose
        // see Maple, two minima
        // NB: C2 is negative
        // integral ∫_1^C2 u(rr)^2 drr, see Penrose.mw
        /* from Maple:
           ss.A=log(c^2)*c^2*(-4*a^4*b^4
            +c^2*(-12*((a^2 + 1)*b^2 + a^2)*a^2*b^2
            +c^2*((-4*b^4 - 16*b^2 - 4)*a^4 + (-16*b^4 - 16*b^2)*a^2 - 4*b^4
            +c^2*((-2*b^2 - 2)*a^4 + (-2*b^4 - 8*b^2 - 2)*a^2 - 2*b^4 - 2*b^2))))
            -1/105 + ((-10*b^4 + 5*b^2 - 1)*a^4)/30 + ((5*b^4 - 4*b^2 + 1)*a^2)/30 - b^4/30 + b^2/30
            +c^2*(1/15 + ((-28*b^4 - 8*b^2 + 1)*a^4)/3 + 4*(-10*b^4 + 5*b^2 - 1)*a^2/15 + b^4/3 - (4*b^2)/15
            +c^2*(-1/5 + 2*(-9*b^2 - 1)*a^4 + (-18*b^4 - 8*b^2 + 1)*a^2 - 2*b^4 + b^2
            +c^2*(1/3 + 4*(7*b^4 + 10*b^2 - 2)*a^4/3 + 8*(5*b^4 - 4*b^2 - 1)*a^2/3 - (8*b^4)/3 - (8*b^2)/3
            +c^2*(((2*b^4 + 43*b^2 + 25)*a^4)/6 + ((43*b^4 + 100*b^2 + 7)*a^2)/6 + (25*b^4)/6 + (7*b^2)/6 - 1/3
            +c^2*(a^4/5 + ((4*b^2 + 4)*a^2)/5 + b^4/5 + (4*b^2)/5 + 1/5
            +c^2*(-a^2/15 - b^2/15 - 1/15
            +c^2/105))))));  */
        ss.A=log(ss.C2q)*ss.C2q*(-4*Sqr(ss.aq)*Sqr(ss.bq)
          +ss.C2q*(-12*((ss.aq + 1)*ss.bq + ss.aq)*ss.aq*ss.bq
          +ss.C2q*((-4*Sqr(ss.bq) - 16*ss.bq - 4)*Sqr(ss.aq) + (-16*Sqr(ss.bq) - 16*ss.bq)*ss.aq - 4*Sqr(ss.bq)
          +ss.C2q*((-2*ss.bq - 2)*Sqr(ss.aq) + (-2*Sqr(ss.bq) - 8*ss.bq - 2)*ss.aq - 2*Sqr(ss.bq) - 2*ss.bq))))
          -1./105 + ((-10*Sqr(ss.bq) + 5*ss.bq - 1)*Sqr(ss.aq))/30 + ((5*Sqr(ss.bq) - 4*ss.bq + 1)*ss.aq)/30 - Sqr(ss.bq)/30 + ss.bq/30
          +ss.C2q*(1./15 + ((-28*Sqr(ss.bq) - 8*ss.bq + 1)*Sqr(ss.aq))/3 + 4*(-10*Sqr(ss.bq) + 5*ss.bq - 1)*ss.aq/15 + Sqr(ss.bq)/3 - (4*ss.bq)/15
          +ss.C2q*(-1./5 + 2*(-9*ss.bq - 1)*Sqr(ss.aq) + (-18*Sqr(ss.bq) - 8*ss.bq + 1)*ss.aq - 2*Sqr(ss.bq) + ss.bq
          +ss.C2q*(1./3 + 4*(7*Sqr(ss.bq) + 10*ss.bq - 2)*Sqr(ss.aq)/3 + 8*(5*Sqr(ss.bq) - 4*ss.bq - 1)*ss.aq/3 - (8*Sqr(ss.bq))/3 - (8*ss.bq)/3
          +ss.C2q*(((2*Sqr(ss.bq) + 43*ss.bq + 25)*Sqr(ss.aq))/6 + ((43*Sqr(ss.bq) + 100*ss.bq + 7)*ss.aq)/6 + (25*Sqr(ss.bq))/6 + (7*ss.bq)/6 - 1./3
          +ss.C2q*(Sqr(ss.aq)/5 + ((4*ss.bq + 4)*ss.aq)/5 + Sqr(ss.bq)/5 + (4*ss.bq)/5 + 1./5
          +ss.C2q*(-ss.aq/15 - ss.bq/15 - 1./15
          +ss.C2q/105))))));
    
        ss.A=512./1155/ss.A; // normalizing factor to the same integral (512/1155) as for 2DLJ
        ss.P[0]=(4*ss.aq*ss.bq*ss.C2q)*ss.A;
        ss.P[1]=(-2*((ss.bq+1)*ss.aq + ss.bq)*ss.C2q)*ss.A;
        ss.P[2]=((-2*ss.bq-2)*ss.aq - 2*ss.bq)*ss.A;
        ss.P[3]=(4*ss.aq + 4*ss.bq + 2*ss.C2q + 4)*ss.A;
        ss.P[4]=-6*ss.A;

        u=ufull=uDW; f=ffull=fDW;
        break;
   
    case PD:  // penetrable disks, normally c=ss.C2=2
      ss.A=ss.a/Sqr(ss.C2);
      ss.A4=ss.A*4;
      u=ufull=uPD; f=ffull=fPD;
      break;
      
    case WCA: // WCALJ, cutoff=2^{1/4}
      ss.C2q=1.414213562373095;
      ss.C2=sqrt(ss.C2q);
      u=ufull=uWCALJ; f=ffull=fWCALJ;
      break;

    default:
      
      fl_alert("Wrong force field.");
      usleep(4000000);
      exit(0); }
      
  calculateB2(); // change vs rg
}

/*** second virial coefficient ***/
inline
double Mayer(double r) /********************************************** Mayer */
{
  return (exp(-u(r*r)/T)-1)*r;
}

/*
  Gauss 2 point formula for numerical integration, O(dx^4)
*/
double Gauss4(double (*f)(double),double a,double b,int n) /********* Gauss4 */
{
  double h=(b-a)/(n*2),sum=0,q=sqrt((double)1/3);
  int i;

  if (n<=0) return 0;

  for (i=n*2-1; i>0; i-=2) sum += (*f)(a+(i+q)*h) + (*f)(a+(i-q)*h);

  return h*sum;
}

void calculateB2(void) /**************************************** calculateB2 */
{
  int n=(int)(fabs(ss.C2)*128)+32;

  if (ss.C2>MINCUT) // 2DLJ, split at C1
    ss.B2=-PI*(Gauss4(Mayer,0,ss.C1,n*2/3)+Gauss4(Mayer,ss.C1,ss.C2,n/3));
  else if (ss.C2<-1) // DW (double well)
    ss.B2=-PI*Gauss4(Mayer,0,-ss.C2,n);
  else if (ss.C2>=1 && ss.C2<=2) // WCALJ and Penetrable disks
    ss.B2=-PI*Gauss4(Mayer,0,ss.C1,n);
  else
    ss.B2=0;
}
