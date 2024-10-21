// Simulation code
// This code was originally written in Pascal (with Turbo graphics),
// then rewritten to C (with Turbo graphics),
// and then Turbo-graphics was emulated by X11.
// This version includes graphics rewritten to FLTK.

#define MAXN 2000     // max # of atoms: static arrays used
#define MINRHO 0.03   // NUCLEATION: 0.04
#define MAXRHO 1.1
#define PROBDMAX 0.03 // of long MC move, default
#define NPOW 3        // power function for the N-slider

#define MINT 0.1 // minimum temperature
#define MAXT 5 // maximum temperature

#define MINP 0 // minimum pressure
#define MAXP 2 // maximum pressure

#define MING -0.2 // minimum gravity
#define MAXG 0.2 // maximum gravity

#define CUTOFF 4 // default cutoff

#define DSCALE 7 // MC displacement d (in units of Lh) d=L/2*exp(min MC displacement d, in the unit of L/2, for log-based slider

// List of variables exported to the 1st line of sim-files;
// if changed, also the formats must be updated.
#define SIMVARS N,bc,walls,T,L,walldens,dt,d,thermostat,P,dV
#define STR(X) #X

typedef struct {
  double x,y;
} vector;

double gravity;  /* acceleration of gravity */
double walldens=0.75*PI; /* rho_wall = walldens/PI = density of smoothed 
                            atoms on walls; 0.75 ~ triple point liquid */
int N=300; /* # of atoms */
double L=20,Lh; /* box size, Lh=L/2 */
double T=3; /* temperature as NVT parameter */
double P=1; /* pressure as NPT parameter */
double qtau=5; /* tauP/tau, for Berendsen and MTK barostat */
double Tk; /* Tk=kinetic T (MD) or demon T (Creutz MC) */
double bag; /* Creutz daemon bag: static */
double d; /* displacement, in Lh=L/2 */
double dV; /* volume displacement, relative (=in ln(L)) */
double t=0; /* MD time or MC steps (for CP only) */
vector r[MAXN],v[MAXN],a[MAXN]; /* configuration, velocities, accelerations */
vector vlast[MAXN]; /* velocity(t-3*dt/2) for Nose-Hoover + TRVP */
double xi,xiv,xivlast,xivpred; /* Nose-Hoover variable, xivlast=xiv(t-3*dt/2) */
double lambdav,lambdavlast,lambdavpred; /* MTK ln(L) */
double Econserved,Ekin; /* total energy incl. extended degrees of freedom */
int justreset=1; /* reset scaling of energy convergence profile */
int debug=0; /* verbose debug mode */
int circlemethod=2; /* 0=fl_pie  1=fl_circle  2=custom of fl_line */
int iblock,block; // stride moved to speed.stride
int delayed_y_color=0;

// it's called thermostat, but it includes NVE and NPT
// the order is the same as in the menu; AUTO is not used (directly)
enum thermostat_t     {AUTO,       CREUTZ,         METROPOLIS,          MCNPT,             NVE,     BERENDSEN,         NOSE_HOOVER,  ANDERSEN,         MAXWELL,         LANGEVIN,         BUSSI,     MDNPT,             MTK, NTH}
  thermostat=BUSSI,mauto;
const char *thinfo[NTH]={"NVT MC→MD","MC/NVE/Creutz","MC/NVT/Metropolis","MC/NPT/Metropolis","MD/NVE","MD/NVT/Berendsen","MD/NVT/Nosé–Hoover","MD/NVT/Andersen","MD/NVT/Maxwell–Bol.","MD/NVT/Langevin","MD/Bussi CSVR","MD/NPT/Berendsen","MD/NPT/Martyna et al." };
int mcstart; // # of MC steps (Metropolis/NPT) to equilibrate before the selected MD starts
int mymessage;

#define isMC(thermostat) (thermostat>=CREUTZ && thermostat<=MCNPT)
#define isMD(thermostat) (thermostat>=NVE)
#define isNPT(thermostat) (thermostat==MCNPT || thermostat==MDNPT || thermostat==MTK)

enum bctype {BOX,SLIT,PERIODIC} bc=BOX; /* boundary conditions */

// NB: order important - >= used
enum cfg_t {GAS,DIFFUSION,GRAVITY,
            VLE,SLAB,NUCLEATION,
            LIQUID,CAPILLARY,CAVITY,
            CRYSTAL,DEFECT,VACANCY,INTERSTITIAL,
            SENTINEL} cfg=GAS; /* initial cfg. */

int accepted=0,Vaccepted=0; /* # of accepted MC moves */

double acc=0.3; /* acceptance ratio for auto set */
double tau=1; /* MD thermostat time constant */
double range=Sqr(RNBR); /* range^2 for counting neighbors (hotkey 'N') */

double dt=0.02,dtadj=0.02,dtfixed; /* MD timestep, adaptive version, fixed */
struct P_s {
  double xx,yy; /* diagonal components of the pressure tensor */
  double fwx,fwxL,fwy,fwyL; /* direct pressure on walls */
} Pcfg;
double Upot; /* potential energy */

/* the last 2 enum types must be RPROFILE,CPROFILE */
enum measure_t { NONE,QUANTITIES,ENERGY,TEMPERATURE,PRESSURE,VOLUME,INTMOTION,RDF,YPROFILE,RPROFILE,CPROFILE,NMEASURE } measure;
const char *minfo[NMEASURE]={"None","Quantities","Energy","Temperature","Pressure","Volume","Integral of motion","RDF","Yprofile","Rprofile","Cprofile" };

struct sum_s {
  double Pxx,Pyy,U,Ekin,Tk,V,H,Econserved;
  double fwx,fwxL,fwy,fwyL;
  double ar,Var;
} sum;

int boxsize=600;
double scale;

char ch;
char st[32];
int hot=1;
int ndelay=1;

FILE *plb;

#include <sys/time.h>
/* real time in better resolution (10^-6 s if provided by the system) */
double mytime(void) /************************************************ mytime */
{
  struct timeval tv;
  struct timezone tz;

  gettimeofday(&tv,&tz);

  return (unsigned)tv.tv_usec/1e6+tv.tv_sec;
}

inline
double sqr(double x) /************************************************** sqr */
{
  return x*x;
}

double cutoff=CUTOFF;

struct ss_s {
  double C1q; // C1^2
  double C2q; // C2^2 (C2=cutoff)
  double A,A4; // constants in the smoothed part
} ss;

/* truncated potential as a function of rr=r^2, rr<cutoff^2 assumed */
inline
double u(double rr) /***************************************************** u */
{
  if (rr>ss.C2q)
    return 0;
  else if (rr<ss.C1q) {
    rr=1/Sqr(rr);
    return 4*(Sqr(rr)-rr); }
  else {
    return ss.A*Sqr(rr-ss.C2q); }
}

/* truncated forces/r as a function of rr=r^2, rr<cutoff^2 assumed */
inline
double f(double rr) /***************************************************** f */
{
  if (rr>ss.C2q)
    return 0;
  else if (rr<ss.C1q) {
    double r4=Sqr(rr);
    return (32/r4-16)/(rr*r4); }
  else {
    return ss.A4*(rr-ss.C2q); }
}

/* smooth cutoff setup */
void setss(void) /**************************************************** setss */
{
  double oldC1,C1,C2;
  int n;

  ss.C2q=1e99;
  ss.C1q=1e99;

 again:
  C2=cutoff;
  C1=0.7*C2;

  n=0;
  do {
    double den=f(C1*C1)*C1;

    oldC1=C1;
    C1=1+4*u(C1*C1)/C1/den;
    if (fabs(den)<1e-50 || n++>2000 || C1<0 || C2<1.2) {
      cutoff=CUTOFF;
      fl_alert("Cannot determine smooth cutoff.\n\
The default cutoff=%d is used.",CUTOFF);
      goto again; }
    C1=C2/sqrt(C1);
  } while (fabs(1-C1/oldC1) > 1e-15);

  if (C1>0.999*C2 || C1<0.5*C2) {
    cutoff=CUTOFF;
    fl_alert("Cannot determine smooth cutoff.\n\
The default cutoff=%d is used.",CUTOFF);
    goto again; }

  ss.A=u(C1*C1)/Sqr(C2*C2-C1*C1);
  ss.A4=-f(C1*C1)/(C2*C2-C1*C1);
  ss.C2q=C2*C2;
  ss.C1q=C1*C1;
}

int degrees_of_freedom() /******************************* degrees_of_freedom */
{
  if (thermostat==ANDERSEN || thermostat==LANGEVIN || thermostat==MAXWELL) return 2*N;
  else return 2*N-(int)bc;
}

/* very approximate EOS for barostat */
double EOS(double T,double rho) /*************************************** EOS */
{
  return T*rho + (exp(2*pow(rho,1.2))-1);
}

/* attractive walls: site-wall potential */
double uwalla(double d) /******************************************** uwalla */
{
  if (d<=0)
    return -1e30*d;
  else {
    double dd=d*d;
    return walldens*(5./24/Sqr(dd)-1)/dd; }
}

/* attractive walls: site-wall force */
double fwalla(double d) /******************************************** fwalla */
{
  double dd=d*d;

  return walldens*(5./4/Sqr(dd)-2)/(dd*d);
}

/* repulsive walls: site-wall potential */
double uwallr(double d) /******************************************** uwallr */
{
  if (d<=0)
    return -1e30*d;
  else {
    double dd=d*d;
    return walldens*5./24/(Sqr(dd)*dd); }
}

/* repulsive walls: site-wall force */
double fwallr(double d) /******************************************** fwallr */
{
  double dd=Sqr(d);

  return walldens*5./4/(Sqr(dd)*dd*d);
}

double L2rho(double L) /********************************************** L2rho */
{
  switch (bc) {
    case BOX: return N/Sqr(L-1.2);
    case SLIT: return N/(L*(L-1.2));
      // case PERIODIC:
    default: return N/Sqr(L); }
}

double rho2L(double rho) /******************************************** rho2L */
{
  switch (bc) {
    case BOX: return 1.2+sqrt(N/rho);
    case SLIT: return sqrt(0.36+N/rho)+0.6;
      //    case PERIODIC:
    default: return sqrt(N/rho); }
}

typedef double (*walltype)(double r);
walltype uwallx,fwallx,uwally,fwally,uwallxL,fwallxL,uwallyL,fwallyL;
int walls; // flags: 1=x,2=xL,4=y,8=yL; set=attractive, unset=repulsive

void setwalls() /************************************************** setwalls */
{
  if (walls&1) { uwallx=uwalla; fwallx=fwalla; } else { uwallx=uwallr; fwallx=fwallr; }
  if (walls&2) { uwallxL=uwalla; fwallxL=fwalla; } else { uwallxL=uwallr; fwallxL=fwallr; }
  if (walls&4) { uwally=uwalla; fwally=fwalla; } else { uwally=uwallr; fwally=fwallr; }
  if (walls&8) { uwallyL=uwalla; fwallyL=fwalla; } else { uwallyL=uwallr; fwallyL=fwallr; }
}

/* squared vector in nearest image convention */
double Sqrni(double r) /********************************************** Sqrni */
{
  return Sqr(Lh-fabs(Lh-fabs(r)));
}

/* nearest image convention (x or y) */
double ni(double dx) /*************************************************** ni */
{
  if (dx>Lh) dx-=L;
  else if (dx<-Lh) dx+=L;
  return dx;
}

// WARNING: some constants are hardwired and not derived from HISTGRID,HISTMAX
// e.g., axis labeling 0 L/4 L/2 L etc.
#define HISTGRID 50 // 50 histogram bins per unity for RDF
//#define HISTMAX 276 // 5.5*HISTGRID+1, should be even for profiles
#define HISTMAX ((PANELW/2)*2-48) // PANELW=329 => HISTMAX=280
unsigned hist[HISTMAX+1];
double etot[HISTMAX];
int itot=HISTMAX+2; /* wait 2, then initialize; then itot=0..HISTMAX-1 */
double rhosum[HISTMAX];

void randomv(int i,double vv) /************************************* randomv */
{
  v[i].x=vv*rndgauss();
  v[i].y=vv*rndgauss();
}

void put_data() /************************************************** put_data */
{
  // copy information back to sliders
  if (bc==PERIODIC) gravity=0;
  sliders.g->value(gravity);
  sliders.d->value(log(d)/DSCALE);
  sliders.rho->value(log(L2rho(L)));
  sliders.P->value(P);
  sliders.tau->value(log(tau));
  sliders.T->value(log(T));
  sliders.N->value(pow(N,1./NPOW));
  buttons.wallx->value(walls&1);
  buttons.wallxL->value(!!(walls&2));
  buttons.wally->value(!!(walls&4));
  buttons.wallyL->value(!!(walls&8));
}

int molcol[MAXN];

void molblack() /************************************************** molblack */
{
  int i;
  loop (i,0,N) molcol[i]=FL_BLACK;
}

static unsigned int nbrcol[8]={
  OI_ORANGE,
  OI_CYAN,
  OI_GREEN,
  OI_LILAC,
  OI_BLUE,
  OI_RED,
  OI_BLACK, /* 6 */
  OI_YELLOW };

void y_split() /**************************************************** y_split */
// two colors by y-coordinate
{
  double ymedian=Lh,dy=Lh;
  int i,nred;

  for (;;) {
    nred=0;
    loop (i,0,N) nred+=(r[i].y>ymedian);
    if (nred>=N/2 && nred<=N/2+1) break;
    dy*=0.5;
    if (nred>N/2) ymedian+=dy; else ymedian-=dy; }

  if (debug) printf("%d red  %d blue  N=%d  ymedian/L=%g\n",nred,N-nred,N,ymedian/L);

  loop (i,0,N) molcol[i]=r[i].y>ymedian ? OI_RED : OI_BLUE;
}

void randomcolors() /****************************************** randomcolors */
{
  int i;

  loop (i,0,N) molcol[i]=256*irnd(0x1000000);
}

void dtadjset() /************************************************** dtadjset */
/* automatic determination of timestep */
{
  if (dtfixed>0)
    dtadj=dtfixed;
  else {
    dtadj=0.018/sqrt(T);
    if (dtadj>0.03) dtadj=0.03;
    if (thermostat==NOSE_HOOVER || thermostat==MTK || thermostat==BUSSI)
      dtadj=(1-(thermostat==MTK)*0.2)/sqrt(1/Sqr(dtadj)+Sqr(8/tau));
    else if (thermostat==NVE) {
      if (dtadj>0.014) dtadj=0.014; } }
}

void initcfg(enum cfg_t cfg) /************************************** initcfg */
{
  int i,j,k,ii;
  double ff;
  vector rt;

  mauto=BUSSI;
  mcstart=40;

  measure=RDF;
  tau=1;
  colormode->value(0); // black
  drawmode->value(0); // movie
  molsize->value(0); // full size

  thermostat=METROPOLIS;
  d=0.1;
  gravity=0;
  walls=0;
  bc=BOX;
  L=rho2L(0.2); Lh=L/2;
  T=3;

  switch (cfg) {

    case DIFFUSION:
      delayed_y_color++;
      molblack();
      colormode->value(0); // black for now
    case GAS:
      /* initial random configuration in a box */
      loop (i,0,N) {
        r[i].x=0.5+rnd()*(L-1); r[i].y=0.5+rnd()*(L-1); }
      break;

    case GRAVITY:
      /* initial random configuration in a box + gravity */
      bc=SLIT;
      measure=YPROFILE;
      L=rho2L(0.15); Lh=L/2;
      d=1; /* longest */
      loop (i,0,N) {
        /* slight gradient - not exact */
        r[i].x=rnd()*L; r[i].y=0.5+pow(rnd(),0.7)*(L-1); }
      gravity=-0.15;
      break;

    case VLE:
      /* layer of liquid at bottom */
      bc=SLIT;
      L=rho2L(0.2); Lh=L/2;
      walls=8;
      measure=YPROFILE;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=0.8+(rnd()*0.35+0.65)*(L-1); }
      T=0.65;
      break;

    case SLAB:
      /* layer of liquid as a slab */
      bc=PERIODIC;
      measure=YPROFILE;
      L=rho2L(0.3); Lh=L/2;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=(rnd()*0.46+0.27)*L; }
      mcstart=100;
      d=0.3/Lh;
      setd->value(1);
      T=0.6;
      break;

    case NUCLEATION:
      /* initial random configuration, low density */
      bc=PERIODIC;
      measure=ENERGY;
      L=rho2L(0.04); Lh=L/2;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=rnd()*L; }
      T=0.5;
      tau=5; /* higher: coalescence, lower: Ostwald ripening */
      mauto=LANGEVIN;
      break;

    case CAPILLARY:
      /* box, layer at bottom */
      L=rho2L(0.35); Lh=L/2;
      loop (i,0,N) {
        r[i].x=0.5+rnd()*(L-1); r[i].y=0.5+(rnd()*0.5+0.5)*(L-1); }
      T=0.6;
      gravity=-0.01;
      walls=15-4;
      mcstart=60;
      measure=YPROFILE;
      break;

    case CAVITY:
      /* cavity in attractive box */
      measure=CPROFILE;
      L=rho2L(0.5); Lh=L/2;
      walls=15;
      mcstart=60;
      loop (i,0,N) {
        do {
          r[i].x=rnd()*L; r[i].y=rnd()*L; }
        while (Sqr(r[i].x-Lh)+Sqr(r[i].y-Lh)<N/6); }
      T=0.6;
      break;

    case LIQUID:
      /* big droplet */
      bc=PERIODIC;
      L=rho2L(0.15); Lh=L/2;
      measure=RPROFILE;
      loop (i,0,N) {
        do {
          r[i].x=rnd()*L; r[i].y=rnd()*L; }
        while (Sqr(r[i].x-Lh)+Sqr(r[i].y-Lh)>N/2); }
      T=0.6;
      mcstart=50;
      break;

    default: // crystals
      L=rho2L(0.25); Lh=L/2;
      k=sqrt(N/3)+1;
      do {
        k--;
        N=3*k*(k+1)+1; // ideal hexagon of edge k+1
        switch (cfg) {
          case DEFECT: T=0.12; N-=k+1; break;
          case VACANCY: T=0.24; N--; break;
          case INTERSTITIAL: T=0.12; N++; break;
          default: T=0.2; }
      } while (N>=MAXN);

      ii=0;
      rt.x=1.155; /* optimized for cutoff=4 */
      rt.y=0.8660254*rt.x;

      for (i=k; i>=0; i--) {
        for (j=-k; j<=k-i; j++) {
          if (cfg==VACANCY && i==0 && j==1) ii--;
          r[ii].x=Lh+(i*0.5+j)*rt.x;
          r[ii].y=Lh+i*rt.y;
          if (i!=0) {
            ii++;
            r[ii].x=r[ii-1].x;
            r[ii].y=L-r[ii-1].y; }
          ii++; } }

      if (cfg==INTERSTITIAL) {
        r[ii].x=Lh+0.5;
        r[ii].y=Lh+0.5; }

      ff=sqrt(T)*0.3;

      loop (i,0,N)  {
        r[i].x+=rnd()*ff; r[i].x+=rnd()*ff; }

      if (cfg==DEFECT) loop (i,0,N) if (r[i].x>Lh) {
          ff=sqrt((r[i].x-Lh)/k)*0.43;
          if (r[i].y>Lh) { r[i].x+=ff; r[i].y+=ff; }
          else { r[i].x-=ff; r[i].y+=ff; } }
      break;
  }

  if (debug) {
    printf("# initial cfg: N=%d cfg=%d\n",N,cfg);
    loop (i,0,N) printf("%g %g\n",r[i].x,r[i].y); }

  dtadjset();
  if (cfg>=CRYSTAL) colormode->value(2);
  put_data();
}

void insdel(int newN) /********************************************** insdel */
// change the number of particles
{
  int i,j;

  if (N==newN) return;

  if (newN<1) newN=1;
  if (newN>MAXN) newN=MAXN;

  if (newN>N) // insert
    loop (i,N,newN) {
      if (bc==BOX) r[i].x=0.8+rnd()*(L-1.6);
      else r[i].x=rnd()*L;
      if (bc==BOX || bc==SLIT) r[i].y=0.8+rnd()*(L-1.6);
      else r[i].y=rnd()*L; }
  else if (newN<N) // delete
    for (i=N; i>newN; i--) {
      j=(int)(i*rnd());
      if (j<i-1) r[j]=r[i-1]; }
  N=newN;
}

void addPonly(double dx,double dy) /******************************* addPonly */
{
  double rr=Sqr(dx)+Sqr(dy);
  double ff=f(rr);

  //  P+=rr*f(rr);
  Pcfg.xx+=dx*dx*ff;
  Pcfg.yy+=dy*dy*ff;
  Upot+=u(rr);
}

void addPRDF(double dx,double dy) /********************************* addPRDF */
{
  unsigned ir;
  double rr=Sqr(dx)+Sqr(dy);
  double ff=f(rr);

  //  P+=rr*f(rr);
  Pcfg.xx+=dx*dx*ff;
  Pcfg.yy+=dy*dy*ff;
  Upot+=u(rr);
  ir=(unsigned)(sqrt(rr)*(double)HISTGRID);
  if (ir<=HISTMAX) hist[ir]++;
}

void (*addP)(double dx,double dy);

vector shift;

void cfgcenter() /************************************************ cfgcenter */
/* returns:
   zero if not periodic
   center of void if periodic
   (not good for measure=CPROFILE)
*/
{
  double c,s,q;
  int i;

  shift.x=shift.y=L;

  if (bc==BOX) return;

  // x-periodic mass center
  c=s=0;
  q=2*PI/L;
  loop (i,0,N) s+=sin(r[i].x*q),c+=cos(r[i].x*q);
  shift.x=atan2(s,c)/q+L; // mass center
  shift.x+=Lh; // void center

  if (bc==SLIT) return;

  // y-periodic mass center
  c=s=0;
  loop (i,0,N) s+=sin(r[i].y*q),c+=cos(r[i].y*q);
  shift.y=atan2(s,c)/q+L; // mass center
  shift.y+=Lh; // void center
}

void MDreset(thermostat_t th) /************************************* MDreset */
{
  double vv=sqrt(T);
  int i;

  thermostat=th;
  loop (i,0,N) randomv(i,vv);

  if (thermostat==NOSE_HOOVER || thermostat==MTK) {
    xivlast=xiv=xi=0;
    loop (i,0,N) vlast[i]=v[i];
    if (thermostat==MTK)
      lambdavlast=lambdav=0; }

  dtadjset();
}

void writesim(const char *fn) /************************************ writesim */
{
  FILE *out=fl_fopen(fn,"wt");
  int i;

  if (out) {
    fprintf(out,"%d %d %d %g %g %g %g %g %d %g %g # " STR(SIMVARS) "\n", SIMVARS);
    fprintf(out,"#  x         y\n");
    loop (i,0,N) fprintf(out,"%9.5f %9.5f\n",r[i].x,r[i].y);
    fclose(out); }
}

void MDerror(void) /************************************************ MDerror */
{
  if (debug) writesim("debug.sim");
  if (!mauto) mymessage=1;
  mcstart=256;
  if (N>500) mcstart=128000/N;
  d=0.1;
  setd->value(1);
  mauto=thermostat;
  thermostat = isNPT(mauto) ? MCNPT : METROPOLIS;

  if (std::isnan(L) || L>sqrt(N)*100) {
    // NPT failed (e.g., P<0) => total reset
    int i;

    L=sqrt(N)*10;
    sliders.rho->value(log(L2rho(L)));
    loop (i,0,N) {
      r[i].x=rnd()*(L-1)+0.5;
      r[i].y=rnd()*(L-1)+0.5;
      v[i].x=0;
      v[i].y=0; } }
}

void debugtimer(const char *info) /****************************** debugtimer */
{
  static int ns=0;
  static unsigned last;
  static double mylast;

  if (ns==0) {
    fprintf(stderr,"\n\
timing (all times are in s):\n\
* speed: slider \'simulation speed\' value\n\
* timer (frame delay): use slider \'simulation speed\'\n\
* stride=MDsteps/display: use slider \'simulation speed\' (to right)\n\
* block: use slider \'measurement block\'\n");
    fprintf(stderr,"SIM nsteps clockstep tvstep  dispstep  timer stride block speed   abstime\n");
    mylast=mytime();
    last=clock(); }
  else if (ns%debug==0) {
    unsigned t=clock();
    double myt=mytime();

    fprintf(stderr,"%s%6d  %9.6f %9.6f %7.5f   %7.5f  %2d %3d %5.2f %.4f\n",
            info,
            ns,
            (double)(t-last)/((double)CLOCKS_PER_SEC*debug),
            (myt-mylast)/debug,
            (myt-mylast)*speed.stride/debug,
            speed.timerdelay,
            speed.stride,
            block,
            sliders.speed->value(),
            myt);
    last=t;
    mylast=myt; }
  ns++;
}

void MDstep(void) /************************************************** MDstep */
{
  vector momentum,ri,rt;
  double ff,vv,xx,maxvv;
  double vf;
  int i,j,k,l;
  int ndegf=degrees_of_freedom();
  int nimg=cutoff/L+1; // # of replicas if cutoff<Lh

  dtadjset();
  dt=dtadj;

  /* calculate forces (accelerations) caused by walls and gravity */
  loop (i,0,N)
    if (bc==PERIODIC) {
      a[i].x=a[i].y=0; }
    else {
      ri=r[i];
      a[i].y=fwally(ri.y)-fwallyL(L-ri.y)-gravity;
      if (bc==BOX) a[i].x=fwallx(ri.x)-fwallxL(L-ri.x);
      else a[i].x=0; }

  /* calculate atom-atom pair forces */
  switch (bc) {
    case BOX:
      loop (i,0,N) {
        ri=r[i];
        loop (j,i+1,N) {
          rt.x=ri.x-r[j].x; rt.y=ri.y-r[j].y;
          ff=f(Sqr(rt.x)+Sqr(rt.y));
          a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
          a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
      break;
    case SLIT:
      if (cutoff<Lh)
        loop (i,0,N) {
          ri=r[i];
          loop (j,i+1,N) {
            rt.x=ri.x-r[j].x;
            if (rt.x>Lh) rt.x=rt.x-L;
            if (rt.x<-Lh) rt.x=rt.x+L;
            rt.y=ri.y-r[j].y;
            ff=f(Sqr(rt.x)+Sqr(rt.y));
            a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
            a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
        else
          loop (i,0,N) {
            ri=r[i];
            loop (j,i+1,N)
              loopto (k,-nimg,nimg) {
                rt.x=ri.x-r[j].x+L*k;
                rt.y=ri.y-r[j].y;
                ff=f(Sqr(rt.x)+Sqr(rt.y));
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
      break;
    case PERIODIC:
      if (cutoff<Lh)
        loop (i,0,N) {
          ri=r[i];
          loop (j,i+1,N) {
            rt.x=ri.x-r[j].x;
            if (rt.x>Lh) rt.x=rt.x-L;
            if (rt.x<-Lh) rt.x=rt.x+L;
            rt.y=ri.y-r[j].y;
            if (rt.y>Lh) rt.y=rt.y-L;
            if (rt.y<-Lh) rt.y=rt.y+L;
            ff=f(Sqr(rt.x)+Sqr(rt.y));
            a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
            a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
      else
        loop (i,0,N) {
          ri=r[i];
          loop (j,i+1,N) {
            loopto (k,-nimg,nimg) {
              rt.x=ri.x-r[j].x+L*k;
              loopto (l,-nimg,nimg) {
                rt.y=ri.y-r[j].y+L*l;
                ff=f(Sqr(rt.x)+Sqr(rt.y));
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } } }}
  }

  if (thermostat==LANGEVIN) {
    /* random force by the fluctuation-dissipation theorem */
    vv=sqrt(T/(tau*dt));
    loop (i,0,N) { a[i].x+=vv*rndgauss(); a[i].y+=vv*rndgauss(); } }

  if (thermostat==NOSE_HOOVER || thermostat==MTK) {
    /* velocity predictor TRVP(k=1), and acceleration modification
       see J. Chem. Theory Comput. 2011, 7, 3596–3607 
       @-codes below refer to the 3D MTK algorithm in the MACSIMUS manual */

    vf=xivpred=1.5*xiv-0.5*xivlast; // @4
    xivlast=xiv;

    if (thermostat==MTK) {
      lambdavpred=1.5*lambdav-0.5*lambdavlast; // @4
      vf+=(1+2./ndegf)*lambdavpred; // @41
      lambdavlast=lambdav; }

    loop (i,0,N) {
      a[i].x-=(1.5*v[i].x-0.5*vlast[i].x)*vf; // @6
      a[i].y-=(1.5*v[i].y-0.5*vlast[i].y)*vf; // @6
      vlast[i]=v[i]; } }

  /* leap-frog integrator and kinetic energy */
  vv=0;
  maxvv=0;
  momentum.x=momentum.y=0;

  loop (i,0,N) {
    /* crash integrator test: this form detects NaN, too */
    if (!(Sqr(a[i].x)+Sqr(a[i].y)<1e7)) {
      if (debug)
        fprintf(stderr,"MDerror: i=%d r=[%g %g] v=[%g %g] a=[%g %g]\n",i,r[i].x,r[i].y,v[i].x,v[i].y,a[i].x,a[i].y);
      MDerror();
      return; }

    xx=Sqr(v[i].x)+Sqr(v[i].y);
    vv+=xx;  /* leap-frog Ekin (at t-dt/2) */

    /* leap-frog, momentum */
    momentum.x+=v[i].x+=dt*a[i].x;
    momentum.y+=v[i].y+=dt*a[i].y;

    xx=Sqr(v[i].x)+Sqr(v[i].y);
    vv+=xx;  /* leap-frog Ekin (at t+dt/2) */
    if (xx>maxvv) maxvv=xx;

    r[i].x+=dt*v[i].x;
    r[i].y+=dt*v[i].y;

    //    if (r[i].x<0 || r[i].x>L || r[i].y<0 || r[i].y>L ) fprintf(stderr,"%8.4f %8.4f\n",r[i].x/L, r[i].y/L);
  } /*N*/

  // remove total momentum - does not have to be done every step...
  if ((thermostat>=NVE && thermostat<=NOSE_HOOVER) || thermostat>=BUSSI) {
    momentum.x/=N; momentum.y/=N;
    if (bc==PERIODIC) loop (i,0,N) { v[i].x-=momentum.x; v[i].y-=momentum.y; }
    else if (bc==SLIT) loop (i,0,N) { v[i].x-=momentum.x; } }

  if (ndegf<=0) {
    fl_alert("\
No degree of freedom!\n\
(Periodic boundary conditions, N=1,\n\
and momentum conservation.)\n\
Switching to BOX boundary conditions...");
    ndegf=2;
    bc=BOX; }

  /* vv=4*Ekin(t), leap-frog (VERLET=3 in MACSIMUS) style */
  Ekin=Econserved=vv/4;
  Tk=2*Ekin/ndegf; /* kinetic temperature */

  double xikin,xipot,lakin; // DEBUG - erase

  /* MTK and NOSE_HOOVER extended masses: */
  double M_T=ndegf*T*Sqr(tau);
  double M_P=2*(ndegf+2)*T*Sqr(tau*qtau);

  switch (thermostat) {

    case MDNPT: {
      /* BERENDSEN barostat */
      double Ptr=(Pcfg.xx+Pcfg.yy)/2;
      double rho=L2rho(L);
      /* rough estimate of bulk modulus: */
      double B=(EOS(T,rho+1e-2)-EOS(T,rho-1e-2))/2e-2;

      ff=exp(dt*(Ptr-P)/(B*tau*qtau)); /* tau.P = tau*qtau, cf. MTK */
      if (L*ff>1e10) ff=1; // added in 12/2021
      L*=ff; Lh=L/2;
      sliders.rho->value(log(L2rho(L)));
      loop (i,0,N) { r[i].x*=ff; r[i].y*=ff; } }
      /* continue by the BERENDSEN thermostat */

    case BERENDSEN:
      /* friction thermostat: rescale velocities, correlation time ~ tau */
      ff=exp((T-Tk)*dt/tau);
      loop (i,0,N) { v[i].x*=ff; v[i].y*=ff; }
      break;

    case BUSSI: {
      /* canonical sampling through velocity rescaling
         DOI: 10.1016/j.cpc.2008.01.006 */
      double c=exp(-dt/tau); // c
      double Rt=rndgauss(); // R(t)
      double Si=0; // S_{Nf-1}

      loop (i,1,ndegf) Si+=sqr(rndgauss()); /* S_{Nf-1} (can be optimized) */

      ff=T/(ndegf*Tk); // av(K)/Nf/K(t)
      ff=sqrt(c + (1-c)*(Rt*Rt+Si)*ff + 2*Rt*sqrt(ff*c*(1-c))); // alpha(t)
      loop (i,0,N) { v[i].x*=ff; v[i].y*=ff; } }
      // NB: change of sign not considered because improbablu
      break;

    case MAXWELL:
      /* Maxwell-Boltzmann thermostat: random hits */
      vv=sqrt(T);
      static double counter=0;
      counter-=dt;
      if (counter<0) {
        loop (i,0,N) randomv(i,vv);
        counter+=tau; }
      break;

    case ANDERSEN:
      /* Andersen thermostat: random hits */
      vv=sqrt(T);
      loop (i,0,N) if (rnd()<dt/tau) randomv(i,vv);
      break;

    case LANGEVIN:
      /* friction part by the fluctuation-dissipation theorem
         - the simplest version giving best Epot */
      ff=1-0.5*dt/tau;
      loop (i,0,N) { v[i].x*=ff; v[i].y*=ff; }
      break;

    case MTK:
      xikin=M_T*Sqr(xiv/2); // 1/2 of leap-frog kin. energy
      xiv+=dt/M_T*(2*(Ekin+M_P*Sqr(lambdavpred))-(ndegf+1)*T); // @A1
      xikin+=M_T*Sqr(xiv/2); // +1/2 of leap-frog kin. energy
      xipot=(ndegf+1)*T*xi; // xi(t)
      xi+=dt*xiv; // xi(t+dt) @B1
      Econserved+=xikin+xipot; // Ekin included, +U later...
      lakin=M_P*Sqr(lambdav)/2; // NB: 2 degrees of freedom have energy M*v^2
      lambdav+=dt*((Sqr(L)*((Pcfg.xx+Pcfg.yy)/2-P)+Tk)/M_P-xivpred*lambdavpred);
      lakin+=M_P*Sqr(lambdav)/2;
      Econserved+=lakin+P*Sqr(L);
      // lambda=log(L);
      ff=exp(dt*lambdav);
      if (L*ff>1e10) ff=1; // added in 12/2021
      L*=ff; Lh=L/2;
      sliders.rho->value(log(L2rho(L)));
      loop (i,0,N) { r[i].x*=ff; r[i].y*=ff; }
      // printf("%g %g %g %g %g %g %g\n",Econserved+Upot,Ekin,Upot,xikin,xipot,lakin,P*Sqr(L));
      break;

    case NOSE_HOOVER:
      xikin=M_T*Sqr(xiv/2);
      xiv+=dt*(Tk/T-1)/Sqr(tau);
      xikin+=M_T*Sqr(xiv/2);
      xipot=ndegf*T*xi;
      Econserved+=xikin+xipot; // Ekin included, +U later...
      xi+=dt*xiv; // xi(t+dt)
      //      printf("%g %g %g %g %g\n",Econserved+Upot,Ekin,Upot,xikin,xipot);
      break;

    default:; /* to make the compiler happy */
  }

  t+=dt;

  if (debug) debugtimer("MD");
} /* MD */

void MCsweep(void) /************************************************ MCsweep */
{
  vector rt,ri,incirc,drt,dri; /* trial position, r[i], random vector */
  int i,j,k,l,nimg;
  double deltaU; /* MC: Utrial-U */
  int isdauto=setd->value();

  /* max displacement = Lh=L/2 */
  Lh=L/2;
  nimg=cutoff/L+1; // # of replicas if cutoff>Lh

  Tk=0; /* Creutz temperature */

  loop (i,0,N) {
    if (d>1) d=1;

    /* trial move of atom i */

    ri=r[i]; /* old position */

    /* rnd() displacement in a 1-circle */
    do {
      incirc.x=rndcos();
      incirc.y=rndcos(); }
    while (Sqr(incirc.x)+Sqr(incirc.y)>1);

    /* new (trial) position (not too much efficient but simple) */
    /* long moves removed */
    rt.x=ri.x+incirc.x*(d*Lh);
    /* normalize to the box, fool-proof
       even for box+slit b.c. which prevents out-of-box problems */
    while (rt.x<0) rt.x+=L;
    while (rt.x>L) rt.x-=L;

    /* ... and the same for y coord */
    rt.y=ri.y+incirc.y*(d*Lh);
    while (rt.y<0) rt.y=rt.y+L;
    while (rt.y>L) rt.y=rt.y-L;

    /* energy difference caused by the trial move: walls and gravity */
    deltaU=0;

    if (bc<=SLIT) {
      deltaU=uwally(rt.y)-uwally(ri.y)
        +uwallyL(L-rt.y)-uwallyL(L-ri.y)
        -gravity*(ri.y-rt.y);
      if (bc==BOX)
        deltaU+=
          uwallx(rt.x)-uwallx(ri.x)
          +uwallxL(L-rt.x)-uwallxL(L-ri.x);
    }

    /* : all other atoms ... not efficient, old energies calculated
       again and again though they might have been stored */
    switch (bc) {
      case BOX:
        loop (j,0,N) if (i!=j)
          deltaU+=u(Sqr(rt.x-r[j].x)+Sqr(rt.y-r[j].y))-u(Sqr(ri.x-r[j].x)+Sqr(ri.y-r[j].y));
        break;
      case SLIT:
        if (cutoff<Lh) {
          loop (j,0,N) if (i!=j)
            deltaU+=u(Sqrni(rt.x-r[j].x)+Sqr(rt.y-r[j].y))-u(Sqrni(ri.x-r[j].x)+Sqr(ri.y-r[j].y)); }
        else
          loop (j,0,N) if (i!=j)
            loopto (k,-nimg,nimg) {
              dri.x=r[j].x+k*L-ri.x;
              drt.x=r[j].x+k*L-rt.x;
              deltaU+=u(Sqr(drt.x)+Sqr(rt.y-r[j].y))-u(Sqr(dri.x)+Sqr(ri.y-r[j].y)); }
        break;
      case PERIODIC:
        if (cutoff<Lh) {
          loop (j,0,N) if (i!=j)
            deltaU+=u(Sqrni(rt.x-r[j].x)+Sqrni(rt.y-r[j].y))-u(Sqrni(ri.x-r[j].x)+Sqrni(ri.y-r[j].y)); }
        else
          loop (j,0,N) if (i!=j)
            loopto (k,-nimg,nimg) {
              dri.x=ri.x-r[j].x+k*L;
              drt.x=rt.x-r[j].x+k*L;
              loopto (l,-nimg,nimg) {
                dri.y=ri.y-r[j].y+l*L;
                drt.y=rt.y-r[j].y+l*L;
                deltaU+=u(Sqr(drt.x)+Sqr(drt.y))-u(Sqr(dri.x)+Sqr(dri.y)); } }
    }

    if (thermostat!=CREUTZ) bag=-T*log(rnd());

    /* Metropolis/Creutz test */
    if (deltaU<bag) {
      /* move accepted */
      bag-=deltaU; // not needed for METROPOLIS
      accepted++;
      r[i]=rt;

      /* automatic adjustment of displacement size to reach
         acceptance ratio = acc; tiny violation of
         microreversibility is usually acceptable */
      if (isdauto) d*=1+(1-acc)*(0.003/sqrt(N)); }
    else {
      /* move rejected */
      /* automatic adjustment of displacement size */
      if (isdauto) d*=1-acc*(0.003/sqrt(N)); }

    if (d>1) d=1;

    Tk+=bag;
  } /*N*/
  Tk/=N;

  if (isdauto) sliders.d->value(log(d)/DSCALE);

  if (thermostat==MCNPT) {
    double lf=rndcos()*dV,f=exp(lf),ff=f*f,rr;
    int isni=cutoff<Lh && cutoff<Lh*f;

    nimg=cutoff/fmin(L,L*f)+1; // # of replicas if !isni

    deltaU=0;

    loop (i,0,N) {

      /* energy difference caused by the trial volume change: walls and gravity */
      if (bc<=SLIT) {
        deltaU+=
          uwally(r[i].y*f)-uwally(r[i].y)
          +uwallyL((L-r[i].y)*f)-uwallyL(L-r[i].y)
          -gravity*(r[i].y*(f-1));
        if (bc==BOX)
          deltaU+=
            uwallx(r[i].x*f)-uwallx(r[i].x)
            +uwallxL((L-r[i].x)*f)-uwallxL(L-r[i].x); }

      /* : inefficient, see above */
      switch (bc) {
        case BOX:
          loop (j,0,i) {
            rr=Sqr(r[i].x-r[j].x)+Sqr(r[i].y-r[j].y);
            deltaU+=u(rr*ff)-u(rr); }
          break;
        case SLIT:
          if (isni)
            loop (j,0,i) {
              rr=Sqrni(r[i].x-r[j].x)+Sqr(r[i].y-r[j].y);
              deltaU+=u(rr*ff)-u(rr); }
          else
            loop (j,0,i)
              loopto (k,-nimg,nimg) {
                drt.x=r[i].x-r[j].x+k*L;
                drt.y=r[i].y-r[j].y;
                rr=Sqr(drt.x)+Sqr(drt.y);
                deltaU+=u(rr*ff)-u(rr); }
          break;
        case PERIODIC:
          if (isni)
            loop (j,0,i) {
              rr=Sqrni(r[i].x-r[j].x)+Sqrni(r[i].y-r[j].y);
              deltaU+=u(rr*ff)-u(rr); }
          else
            loop (j,0,i)
              loopto (k,-nimg,nimg) {
                drt.x=r[i].x-r[j].x+k*L;
                loopto (l,-nimg,nimg) {
                  drt.y=r[i].y-r[j].y+l*L;
                  rr=Sqr(drt.x)+Sqr(drt.y);
                  deltaU+=u(rr*ff)-u(rr); } }
      } /* bc */
    }

    if (rnd()<exp(lf*(2*N+2)-(P*L*L*(ff-1)+deltaU)/T)) {
      Vaccepted++;
      if (isdauto) {
        if (dV<0.1) dV*=1+(1-acc)*0.01; }
      rr=L2rho(L);
      if (rr<1e-5) {
        mymessage=2;
        f=0.8; /* shrink box */
        rr/=f*f;
        if (debug) writesim("debug.sim");
        thermostat=METROPOLIS;
        d=0.3/Lh; }
      if (L*f>1e10) ff=1; // added in 12/2021
      L*=f; Lh=L/2;
      sliders.rho->value(log(rr));
      loop (i,0,N) { r[i].x*=f; r[i].y*=f; } }
    else {
      /* move rejected */
      if (isdauto) {
        if (dV>0.001/sqrt(N)) dV*=1-acc*0.01; } }
  } /* NPT */

  if (debug) debugtimer("MC");

  t+=1;
} /* MC */

int shownbrs;

void neighbors() /************************************************ neighbors */
{
  int i,j;
  int nbr;
  vector rt;

  loop (i,0,N) {
    rt=r[i];
    nbr=0;
    switch (bc) {
      case BOX:
        loop (j,0,N) if (i!=j)
          if (Sqr(rt.x-r[j].x)+Sqr(rt.y-r[j].y)<range) nbr++;
        break;
      case SLIT:
        loop (j,0,N) if (i!=j)
          if (Sqrni(rt.x-r[j].x)+Sqr(rt.y-r[j].y)<range) nbr++;
        break;
      case PERIODIC:
        loop (j,0,N) if (i!=j)
          if (Sqrni(rt.x-r[j].x)+Sqrni(rt.y-r[j].y)<range) nbr++;
        break; }

    if (nbr>7) nbr=7;
    molcol[i]=nbrcol[nbr]; }
}

// return the last single extension .ext (last .)
// FILE. => zero-length string
// BUG: \ is treated as directory seperator in linux
const char *getext(const char *fn) /********************************* getext */
{
  const char *c,*d;

  for (c=fn,d=NULL; *c; c++) {
    if (*c=='.') d=c;
    // I do not know how the directory separator works over platforms...
    if (*c=='/') d=NULL;
    if (*c=='\\') d=NULL; }

  return d;
}

void loadsim(const char *fn) /************************************** loadsim */
{
  FILE *in=fl_fopen(fn,"rt");
  char line[256];
  int i;

  if (in) {
    if (fgets(line,256,in)) {
      if (!strstr(line,STR(SIMVARS))) {
        fl_alert("read \"%s\": is not a SIMOLANT file",fn);
        fclose(in);
        return; }
      sscanf(line,"%d %d %d %lf %lf %lf %lf %lf %d %lf %lf",
             // see macso SIMVARS for the list of variables
             &N,(int*)&bc,&walls,&T,&L,&walldens,&dt,&d,(int*)&thermostat,&P,&dV); }
    if (!fgets(line,256,in)) {
      fl_alert("read \"%s\": bad header or bad format",fn);
      fclose(in);
      return; }

    put_data();

    loop (i,0,N)
      if (fgets(line,256,in))
        sscanf(line,"%lf %lf",&r[i].x,&r[i].y);
      else {
        fl_alert("read \"%s\": too short file or bad format",fn);
        fclose(in);
        return ; }

    fclose(in); }
  else {
    fl_alert("file \"%s\" unreadable or not found",fn);
    fclose(in);
    return; }
}
