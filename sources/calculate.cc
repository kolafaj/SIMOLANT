/*
   Simulation code.  It originates in old DOS Pascal+Turbo graphics,
   then it was rewritten to C (with Turbo C graphics),
   then Turbo-graphics was emulated by X11.
   This version includes graphics rewritten to C++/FLTK, but the code is still
   rather C than C++ in design.
*/

#define MAXN 3000     // max # of atoms: static arrays used
#define MINRHO 0.01   // for slider (NUCLEATION 0.04, isotherms 0.01)
#define MAXRHO 1.13   // for slider
#define MINRHONPT 1e-5// for MCNPT and MDNPT(Berendsen)
#define PROBDMAX 0.03 // of long MC move, default
#define MAXDV 0.25    // max dV: box scaling is max exp(MAXDV)
// the N-slider is cubic (not logarithmic nor linear)
#define NtoSlider(N) cbrt((double)N)
#define SlidertoN(V) Cub(V) // rounded to N in insdel()

#define MINT 0.1 // minimum temperature
#define MAXT 5 // maximum temperature

#define MINP 0.001 // minimum pressure
#define MAXP 2 // maximum pressure

#define MING -0.2 // minimum gravity
#define MAXG 0.2 // maximum gravity

#define MINC 1.0 // minimum cutoff
#define MAXC 100.0 // maximum cutoff

#define DSCALE 7 // MC displacement d (in units of Lh) d=L/2*exp(min MC displacement d), in the unit of L/2, for log-based slider

#define TRVP 2 // time-reversible velocity predictor length
#if TRVP<1 || TRVP>2
#  error "unsupported TRVP (1 or 2 allowed, 2 recommended)"
#endif

// List of variables exported to the 1st line of sim-files;
// if changed, also the formats must be updated.
#define SIMVARS   N,bc,walls,T,L,walldens,dt,d,method,P,dV,tau,qtau,gravity,ss.a,ss.b,ss.C2
#define SIMVARSQ "N bc walls T L wall dt d method P dV tau qtau g a b c"

vector zero={0,0};

double gravity;  /* acceleration of gravity */
double walldens=0.75*PI; /* rho_wall = walldens/PI = density of smoothed
                            atoms on walls; 0.75 ~ triple point liquid */
int N=300; /* # of atoms */
double L=20,Lh; /* box size, Lh=L/2 */
double T=3; /* temperature as NVT parameter */
double P=1; /* pressure as NPT parameter */
double qtau=5; /* tauP/tau, for Berendsen and MTK barostat */
double Tk; /* Tk = averaged kinetic T (MD) or demon T (Creutz MC) */
double bag; /* Creutz daemon bag */
double d=0.5; /* MC displacement, in Lh=L/2 */
double dV=0.1; /* MC volume displacement, relative (=in ln(L)) */
double t=0; /* MD time or MC steps (for CP only) */
vector r[MAXN],v[MAXN],a[MAXN]; /* configuration, velocities, accelerations */
vector r0[MAXN]; /* configuration, at recording start, divided by L */
vector rpbc[MAXN]; /* configuration, followed periodic b.c., divided by L */
int MSDskip; /* skip MSD calculation (set after error cannot follow periodic b.c.) */
vector vlast[MAXN]; /* velocity(t-3*dt/2) for Nose-Hoover + TRVP */
vector vllast[MAXN]; /* velocity(t-5*dt/2) for Nose-Hoover + TRVP */
double wnbrs[MAXN]; /* for SWARM: weight of neighbors */
double xi,xiv,xivlast,xivllast,xivpred; /* Nose-Hoover variable, xivlast=xiv(t-3*dt/2), xivllast=xiv(t-5*dt/2) */
double lambdav,lambdavlast,lambdavllast,lambdavpred; /* MTK ln(L) */
double Econserved; /* total energy incl. extended degrees of freedom */
int justreset=1; /* reset scaling of energy convergence profile */
int debug=0; /* verbose debug mode */
struct circle_s {
  int method=2; /* 0=fl_pie  1=fl_circle  2=custom of fl_line */
  int opt=2; /* locally may use faster method for speed */
} circle;
int iblock,block; // also slider (stride moved to speed.stride)
int nimg; // for cutoff>Lh: number of images
double trace=64; // length for draw mode Traces
struct drop_s { /* for function TWODROPS */
  double v=0.5; /* initial velocity */
  double T=0.4; /* initial temperature */
} drop;

// method: thermostat/barostat/MC/MD
// the order is the same as in the menu; AUTO is not used (directly)
enum method_e                {AUTO,       CREUTZ,         METROPOLIS,         MCNPT,              NVE,     BERENDSEN,         NOSE_HOOVER,         ANDERSEN,         MAXWELL,              LANGEVIN,         BUSSI,          MDNPT,             MTK,                   VICSEK,NTH}
  method=BUSSI,mauto;
//UNICODE problem const char *methodinfo[NTH]={"NVT MC→MD","MC/NVE/Creutz","MC/NVT/Metropolis","MC/NPT/Metropolis","MD/NVE","MD/NVT/Berendsen","MD/NVT/Nosé–Hoover","MD/NVT/Andersen","MD/NVT/Maxwell–Bol.","MD/NVT/Langevin","MD/NVT/Bussi CSVR","MD/NPT/Berendsen","MD/NPT/Martyna et al.","MD/NVT/VICSEK" };
const char *methodinfo[NTH]={"NVT MC->MD","MC/NVE/Creutz","MC/NVT/Metropolis","MC/NPT/Metropolis","MD/NVE","MD/NVT/Berendsen","MD/NVT/Nosé–Hoover","MD/NVT/Andersen","MD/NVT/Maxwell–Bol.","MD/NVT/Langevin","MD/NVT/Bussi CSVR","MD/NPT/Berendsen","MD/NPT/Martyna et al.","MD/NVT/VICSEK" };
int mcstart; // # of MC steps (Metropolis/NPT) to equilibrate before the selected MD starts
#if 1
#define NOERROR 0
#define MDFAILED 1
#define NPTFAILED 2
#define MSDFAILED 3
int errmessage=NOERROR;
#else
// this line causes ERROR "expected identifier before {" on Windows:
enum errmessage_e {NOERROR,MDFAILED,NPTFAILED,MSDFAILED} errmessage=NOERROR;
#endif
int lasterrmessage=NOERROR;

typedef double (*function_t)(double r); // for potential() and forces()
void MDlcPneeded(linklist_t *l);
void MDlc(linklist_t *l);
void MDlcVicsek(linklist_t *l);
void MCNPTlc(linklist_t *l);

int wasMDerror; // previously MDFAILED - to break "Switching temporarily to Monte Carlo…" loop

enum colormode_e {CM_BLACK,CM_ONERED,CM_YSPLIT,CM_NEIGHBORS,CM_RANDOM,CM_ART,CM_KEEP};

#define isMC(method) (method>=CREUTZ && method<=MCNPT)
#define isMD(method) (method>=NVE)
#define isNPT(method) (method==MCNPT || method==MDNPT || method==MTK)
#define isNVE(method) (method==CREUTZ || method==NVE)

enum bctype {BOX,SLIT,PERIODIC} bc=BOX; // boundary conditions

// NB: order important - >= used
enum cfg_t {GAS,DIFFUSION,GRAVITY,
            VLE,SLAB,NUCLEATION,
            LIQUIDDROP,TWODROPS,CAVITY,CAPILLARY,PERIODICLIQUID,
            CRYSTAL,DEFECT,VACANCY,INTERSTITIAL,PERIODICCRYSTAL,
            INITVICSEK,
            SENTINEL} cfg=GAS; /* initial cfg. */

int accepted=0,Vaccepted=0; /* # of accepted MC moves */

double acc=0.3; /* acceptance ratio for auto set */
double tau=1; /* MD thermostat time constant */

double dt=0.02,dtadj=0.02,dtfixed; /* MD timestep, adaptive version, fixed */

struct En_s {
  vector P;  // diagonal components of the pressure tensor, P_cfg
             // (first, the virial only, then the kin part is added)
  vector Ek; // MD: ∑ vx², ∑ vy², then corrected by qx,qy and /2
             // MC: sum of bag energy over N MC steps (x,y the same)
  double Ekin; // MD: total kinetic energy = (Ek.x+Ek.y)/2
               // MC: sum of bag energy over N MC steps
  vector q; // kinetic pressure corrections:
            //   En.Ek*=En.q; P=(virial+Ek)/L² (diagonal terms)
  double fwx,fwxL,fwy,fwyL; // direct pressure on walls
  double Upot; // potential energy
  double H; // enthalpy
  vector momentum; // averaged total momentum (velocity), good for Vicsek
  vector MSD; // mean square displacement
  int Nf; // see degrees_of_freedom();
} En; // similar as in MACSIMUS

/*
  What to show in the right top panel, see menu Show.
  Option NONE was removed in V03/2025, former name measure was changed to Show;
  note that lowercase show collides with a member function show().
  The last 3 enum types must be RPROFILE,CPROFILE,NSHOW,
*/
enum Show_e { QUANTITIES,ENERGY,TEMPERATURE,PRESSURE,VOLUME,MOMENTUM,INTMOTION,RDF,YPROFILE,RPROFILE,CPROFILE,NSHOW }
  Show=QUANTITIES, lastShow=NSHOW;
const char *Showinfo[NSHOW]={"Quantities","Energy","Temperature","Pressure","Volume","Momentum","Integral of motion","RDF","Yprofile","Rprofile","Cprofile" };

/* for mean values over blocks - sums accumulated using addM() */
struct sum_s {
  vector P; // diagonal of pressure tensor
  double U,Ekin,Tk,V,H,Econserved;  // Tk → Tbag for MC
  double fwx,fwxL,fwy,fwyL; // forces on walls
  double accr,Vaccr; // acceptance ratios in MC
  double MSD; // MSD from r0, averaged over N
  double momentum; // from |En.momentum| averaged total momentum (velocity)
} sum; // sums accumulated

/* for variances of selected quantities */
struct sumVar_s {
  double Upot,Ekin,V,H;
} sum0, // 1st value
  sums, // sum (from 1st value)
  sumq; // sums of squares (from 1st value)
int measureVar; // number of measurements in sumq

#include <sys/time.h>

/* real time in s with better resolution (10^-6 s if provided by the system) */
double mytime(void) /************************************************ mytime */
{
  struct timeval tv;
  struct timezone tz;
  static double t0;
  double t;

  gettimeofday(&tv,&tz);

  t=(unsigned)tv.tv_usec/1e6+tv.tv_sec;
  if (!t0) t0=t;

  return t-t0;
}

inline
double sqr(double x) /************************************************** sqr */
{
  return x*x;
}

#include "forcefield.cc"

/* the number of mechanical degrees of freedom, conserved momenta not included */
void degrees_of_freedom() /****************************** degrees_of_freedom */
{
  En.q.x=En.q.y=1;

  En.Nf=2*N;
  if (isMC(method) || method==ANDERSEN || method==LANGEVIN || method==VICSEK || method==MAXWELL) return;

  switch (bc) {
    case BOX: break;
    case SLIT:
      En.q.x=(double)N/(N-1);
      En.Nf=2*N-1; break;
    case PERIODIC:
      En.q.x=En.q.y=(double)N/(N-1);
      En.Nf=2*N-2; break; }

  if (method==MTK) En.q.x=En.q.y=1;
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

/* box size L → number density ρ conversion */
double L2rho(double L) /********************************************** L2rho */
{
  switch (bc) {
    case BOX: return N/Sqr(L-1.2);
    case SLIT: return N/(L*(L-1.2));
      // case PERIODIC:
    default: return N/Sqr(L); }
}

/* number density ρ → box size L conversion */
double rho2L(double rho) /******************************************** rho2L */
{
  switch (bc) {
    case BOX: return 1.2+sqrt(N/rho);
    case SLIT: return sqrt(0.36+N/rho)+0.6;
      //    case PERIODIC:
    default: return sqrt(N/rho); }
}

function_t uwallx,fwallx,uwally,fwally,uwallxL,fwallxL,uwallyL,fwallyL;
int walls; // flags: 1=x,2=xL,4=y,8=yL; set=attractive, unset=repulsive

void setwalls() /************************************************** setwalls */
{
  if (walls&1) { uwallx=uwalla; fwallx=fwalla; } else { uwallx=uwallr; fwallx=fwallr; }
  if (walls&2) { uwallxL=uwalla; fwallxL=fwalla; } else { uwallxL=uwallr; fwallxL=fwallr; }
  if (walls&4) { uwally=uwalla; fwally=fwalla; } else { uwally=uwallr; fwally=fwallr; }
  if (walls&8) { uwallyL=uwalla; fwallyL=fwalla; } else { uwallyL=uwallr; fwallyL=fwallr; }
}

/* squared vector in nearest image convention */
inline
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
#define HISTMAX ((PANELW/2)*2-51) // should be even for profiles
unsigned hist[HISTMAX+1];
double etot[HISTMAX];
int itot=HISTMAX+2; /* wait 2, then initialize; then itot=0..HISTMAX-1 */
double rhosum[HISTMAX];

/* Maxwell-Boltzmann velocity assignment */
void randomv(int i,double vv) /************************************* randomv */
{
  v[i].x=vv*rndgauss();
  v[i].y=vv*rndgauss();
}

/* copy information back to sliders */
void put_data() /************************************************** put_data */
{
  if (bc==PERIODIC) gravity=0;
  sliders.g->value(gravity);
  sliders.d->value(log(d)/DSCALE);
  sliders.rho->value(log(L2rho(L)));
  sliders.P->value(P);
  sliders.tau->value(log(tau));
  sliders.T->value(log(T));
  sliders.N->value(NtoSlider(N));
  buttons.wallx->value(walls&1);
  buttons.wallxL->value(!!(walls&2));
  buttons.wally->value(!!(walls&4));
  buttons.wallyL->value(!!(walls&8));
}

unsigned molcol[MAXN];

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

/* two colors by y-coordinate */
void y_split() /**************************************************** y_split */
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

  loop (i,0,N) molcol[i]=irnd(0x1000000)<<8;
}

/* automatic determination of timestep */
void dtadjset() /************************************************** dtadjset */
{
  if (dtfixed>0)
    dtadj=dtfixed;
  else {
    dtadj=0.018/sqrt(T);
    /* low temp limit; not that for NVE, T may be irrelevant */
    if (dtadj>0.025) dtadj=0.025;
    if (u==uDW)
      if (dtadj>0.02) dtadj=0.02; // DW needs a shorter timestep
    if (method==NOSE_HOOVER || method==MTK || method==BUSSI)
      dtadj=(1-(method==MTK)*0.2)/sqrt(1/Sqr(dtadj)+Sqr(8/tau));
    else if (method==NVE) {
      if (dtadj>0.014) dtadj=0.014; } }
}

void MDreset(method_e th) /***************************************** MDreset */
{
  double vv=sqrt(T);
  int i;

  method=th;

  loop (i,0,N) randomv(i,vv);

  if (method==NOSE_HOOVER || method==MTK) {
    xivlast=xivllast=xiv=xi=0;
    loop (i,0,N) vlast[i]=vllast[i]=v[i];
    lambdavlast=lambdavllast=lambdav=0; }

  degrees_of_freedom();
  dtadjset();
}

void MCreset(method_e th) /***************************************** MDreset */
{
  method=th;
  degrees_of_freedom();
  bag=T;
  if (wasMDerror) mauto=AUTO; // kill the "switching temporarily to Monte Carlo…" loop
}

enum cfg_t lastcfg; // for DIFFUSION (color half AFTER MC) and TWODROPS

void preopt(int from,int to,double xy0,double RR) /****************** preopt */
/* remove worst overlaps for particles in [from,to)
   in disk centered (xy0,xy0) and radius sqrt(rr) */
{
  int i,j;
  vector rt;
  double d=0.3,rr,ff;
  // printf("%d -3\n",N);

  while (d>0.05) {

    loop (i,from,to)
      a[i].x=a[i].y=0;

    loop (i,from,to) {
      rt.x=r[i].x-xy0;
      rt.y=r[i].y-xy0;
      rr=Sqr(rt.x)+Sqr(rt.y);
      ff=rr-RR;
      if (ff>0) { a[i].x-=rt.x*ff; a[i].y-=rt.y*ff; }
      
      loop (j,i+1,to) {
        rt.x=r[i].x-r[j].x; rt.y=r[i].y-r[j].y;
        ff=f(Sqr(rt.x)+Sqr(rt.y));
        a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
        a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
    // printf("%g %g 1\n",L,L);

    loop (i,from,to) {
      ff=Sqr(a[i].x)+Sqr(a[i].y);
      if (ff>1) {
        ff=sqrt(ff);
        a[i].x/=ff;
        a[i].y/=ff; }
      r[i].x+=a[i].x*d;
      r[i].y+=a[i].y*d;
      // printf("%g %g 0\n",r[i].x,r[i].y);
    }
    
    d*=0.9; }
}

void initcfg(enum cfg_t cfg) /************************************** initcfg */
{
  int i,j,k,ii;
  double ff;
  vector rt;
  double RR;
  
  lastcfg=cfg;

  mauto=BUSSI;
  mcstart=40;

  lastShow=NSHOW; // to reset the graph

  Show=RDF;
  tau=1;
  colormode->value(CM_BLACK); // black
  drawmode->value(0); // movie
  molsize->value(0); // full size

  MCreset(METROPOLIS);
  d=1./sqrt(N);
  gravity=0;
  walls=0;
  bc=BOX;
  L=rho2L(0.2); Lh=L/2;
  T=3;

  switch (cfg) {

    case DIFFUSION:
      molblack();
      colormode->value(CM_BLACK); // black for now
    case GAS:
      /* initial random configuration in a box */
      loop (i,0,N) {
        r[i].x=0.5+rnd()*(L-1); r[i].y=0.5+rnd()*(L-1); }
      break;

    case GRAVITY:
      /* initial random configuration in a box + gravity */
      bc=SLIT;
      Show=YPROFILE;
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
      Show=YPROFILE;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=0.8+(rnd()*0.35+0.65)*(L-1); }
      T=0.65;
      break;

    case SLAB:
      /* layer of liquid as a slab */
      bc=PERIODIC;
      Show=YPROFILE;
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
      Show=ENERGY;
      L=rho2L(0.04); Lh=L/2;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=rnd()*L; }
      T=0.5;
      tau=5; /* higher: coalescence, lower: Ostwald ripening */
      mauto=LANGEVIN;
      break;

    case CAVITY:
      /* cavity in attractive box */
      Show=CPROFILE;
      L=rho2L(0.5); Lh=L/2;
      walls=15;
      mcstart=60;
      loop (i,0,N) {
        do {
          r[i].x=rnd()*L; r[i].y=rnd()*L; }
        while (Sqr(r[i].x-Lh)+Sqr(r[i].y-Lh)<N/6); }
      T=0.6;
      break;

    case LIQUIDDROP:
      /* big droplet */
      bc=PERIODIC;
      L=rho2L(0.15); Lh=L/2;
      Show=RPROFILE;
      RR=0.48*N; // RR=0.5*N ⇒ ρ=2/π=0.6366 (better slightly < ρ(l)
      loop (i,0,N) {
        do { r[i].x=rnd()*L; r[i].y=rnd()*L; }
        while (Sqr(r[i].x-Lh)+Sqr(r[i].y-Lh)>RR); }
      preopt(0,N,Lh,RR);
      T=0.6;
      mcstart=20;
      break;

    case TWODROPS:
      /* two droplets */
      bc=PERIODIC;
      L=rho2L(0.15); Lh=L/2;
      Show=TEMPERATURE;
      RR=0.23*N;
      mauto=NVE;
      loop (i,0,N/2) {
        do { r[i].x=rnd()*L; r[i].y=rnd()*L; }
        while (Sqr(r[i].x-Lh/2)+Sqr(r[i].y-Lh/2)>RR);
        molcol[i]=OI_RED; }
      preopt(0,N/2,Lh/2,RR);
      loop (i,N/2,N) {
        do { r[i].x=rnd()*L; r[i].y=rnd()*L; }
        while (Sqr(r[i].x-Lh*1.5)+Sqr(r[i].y-Lh*1.5)>RR);
        molcol[i]=OI_BLUE; }
      preopt(N/2,N,Lh*1.5,RR); 

      d=0.3/sqrt(N);
      mcstart=50; // causes also set velocity of the droplets
      colormode->value(CM_KEEP);
      T=drop.T;
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
      Show=YPROFILE;
      break;

    case PERIODICLIQUID:
      /* initial random configuration, low density */
      bc=PERIODIC;
      Show=RDF;
      L=rho2L(0.73); Lh=L/2;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=rnd()*L; }
      T=0.7;
      P=0.5;
      colormode->value(CM_NEIGHBORS);
      mauto=MDNPT;
      break;

    case PERIODICCRYSTAL: {
      /* 0.87 = experimental density for T=0.3, P=0.5, c=4, N=209
         (accidentally RHOCRYST ~ \0.75) */
#define RHOCRYST 0.87
      /* magic periodic N and crystal: see hexinsquare.c (\ = sqrt)
         Resulting magical N_n = [(2+\3)^n-(2-\3)^n]/2\3 ~ (2+\3)^n/2\3 */
      int n=log(N*3.464101615137754)/1.316957896924817+0.5; // 2\3, ln(2+\3)
      do {
        N=exp(n*1.316957896924817)/3.464101615137754+0.5;
        if (N>MAXN) n--;
        if (N<15) n++; }
      while (N<15 || N>MAXN);
      bc=PERIODIC;
      Show=RDF;
      L=rho2L(RHOCRYST);
      Lh=L/2;

      // n=4, 56, 780: (i+j/2,j*\(3/4)) or i*(1,0) + j*(1/2,\0.75)
      // n=3,5  15, 209: i*(1/\2,1/\2) + j*(a,b); a=(1/2-q)/\2, b=(1/2+q)/\2
      const vector v1t[2]={{1,0},{0.707106781186547,0.707106781186547}};
      const vector v2t[2]={{0.5,0.866025403784439},{-0.258819045102521,0.965925826289068}};
      double q=1/sqrt(RHOCRYST*0.866025403784439); // atom-atom distance = 1 for density rho=1/\0.75
      vector v1={v1t[n&1].x*q,v1t[n&1].y*q};
      vector v2={v2t[n&1].x*q,v2t[n&1].y*q};
      int i,j,nn,maxn=1.5*sqrt(N);

      n=0;
      loop (i,-maxn,maxn) {
        loop (j,-maxn,maxn) {
          r[n].x=v1.x*i+v2.x*j;
          if (r[n].x<0 || r[n].x>=L) continue;
          r[n].y=v1.y*i+v2.y*j;
          if (r[n].y<0 || r[n].y>=L) continue;
          loop (nn,0,n)
            if (Sqrni(r[n].x-r[nn].x) + Sqrni(r[n].y-r[nn].y) < 0.25) goto again;
          n++;
          if (n==N) break;
        again:; }
        if (n==N) break; }

      T=0.3;
      P=0.5;
      colormode->value(CM_NEIGHBORS);
      mauto=MDNPT; }
      break;

    case INITVICSEK:
      // Vicsek implemented if (ss.vicsek) only in periodic+Langevin:
      bc=PERIODIC;
      //  method=VICSEK; // shows [MD/NVT/VICSEK] instead of MD/NVT/VICSEK - why?
      mauto=VICSEK;
      Show=MOMENTUM;
      mcstart=10;
      L=rho2L(0.03); Lh=L/2;
      loop (i,0,N) {
        r[i].x=rnd()*L; r[i].y=rnd()*L;
        // random velocities, T=1
        v[i].x=rndgauss(); v[i].y=rndgauss(); }
      // slightly repulsive potential for Vicsek:
      ss.a=0.2; setss(PD,CUTOFF);
      T=1;
      drawmode->value(1); // traces
      break;

    default: // nonperiodic crystals: CRYSTAL,DEFECT,VACANCY,INTERSTITIAL
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
      colormode->value(CM_NEIGHBORS);
      break;
  }

  if (cfg!=INITVICSEK) setss(LJ,CUTOFF);

  if (debug) {
    printf("# initial cfg: N=%d cfg=%d\n",N,cfg);
    loop (i,0,N) printf("%g %g XY\n",r[i].x,r[i].y); }

  dtadjset();
  put_data();
}

void insdel(double dN) /********************************************* insdel */
// change the number of particles
{
  int i,j,newN=dN+0.5;
  double minrr;

  if (N==newN) return;

  if (newN<1) newN=1;
  if (newN>MAXN) newN=MAXN;

  if (newN>N) // insert
    loop (i,N,newN) {
      minrr=0.68;
    again:
      // "Infinite" loop to insert a particle to a cavity, but the distance
      // limit minrr decreases so that the insertion should succeed eventually
      minrr*=0.999;
      if (bc==BOX) r[i].x=0.8+rnd()*(L-1.6);
      else r[i].x=rnd()*L;
      if (bc==BOX || bc==SLIT) r[i].y=0.8+rnd()*(L-1.6);
      else r[i].y=rnd()*L;
      loop (j,0,i)
        if (sqr(Lh-fabs(Lh-fabs(r[i].x-r[j].x))) +
            sqr(Lh-fabs(Lh-fabs(r[i].y-r[j].y))) < minrr) goto again;
      randomv(i,sqrt(T));
      a[i].x=a[i].y=0; }
  else if (newN<N) // delete
    for (i=N; i>newN; i--) {
      j=(int)(i*rnd());
      if (j<i-1) r[j]=r[i-1]; }
  N=newN;
  degrees_of_freedom();
} // insdel()

void addPonly(double dx,double dy) /******************************* addPonly */
{
  double rr=Sqr(dx)+Sqr(dy);
  double ff=f(rr);

  //  P+=rr*f(rr);
  En.P.x+=dx*dx*ff;
  En.P.y+=dy*dy*ff;
  En.Upot+=u(rr);
} // addPonly

void addPRDF(double dx,double dy) /********************************* addPRDF */
{
  unsigned ir;
  double rr=Sqr(dx)+Sqr(dy);
  double ff=f(rr);

  //  P+=rr*f(rr);
  En.P.x+=dx*dx*ff;
  En.P.y+=dy*dy*ff;
  En.Upot+=u(rr);
  ir=(unsigned)(sqrt(rr)*(double)HISTGRID);
  if (ir<=HISTMAX) hist[ir]++;
} // addPRDF

void (*addP)(double dx,double dy);

vector shift;

void cfgcenter() /************************************************ cfgcenter */
/* returns:
   zero if not periodic
   center of void if periodic
   (not good for Show=CPROFILE)
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
} // cfgcenter()

void savesim(const char *fn) /************************************** savesim */
{
  FILE *out=fl_fopen(fn,"wt");
  int i;

  if (out) {//  N bc walls T L wall dt d method
    //                                      P dV tau qτ g a  b  c"
    fprintf(out,"%d %d %d %g %g %g %g %g %d %g %g %g %g %g %g %g %g # " SIMVARSQ "\n", SIMVARS);
    fprintf(out,"#  x         y        vx        vy\n");
    loop (i,0,N) fprintf(out,"%9.5f %9.5f %9.5f %9.5f\n",r[i].x,r[i].y,v[i].x,v[i].y);
    fclose(out); }
} // savesim()

void MDerror(void) /************************************************ MDerror */
{
  if (debug) savesim("debug.sim");
  if (!mauto) errmessage=MDFAILED;
  mcstart=256;
  if (N>500) mcstart=128000/N;
  d=0.1;
  setd->value(1);
  mauto=method;
  method = isNPT(mauto) ? MCNPT : METROPOLIS;

  if (std::isnan(L) || L>sqrt(N)*100 || L<sqrt(N)*0.1) {
    // NPT failed (e.g., P<0) => total reset
    int i;

    L=sqrt(N)*10;
    sliders.rho->value(log(L2rho(L)));
    loop (i,0,N) {
      r[i].x=rnd()*(L-1)+0.5;
      r[i].y=rnd()*(L-1)+0.5;
      v[i].x=0;
      v[i].y=0; } }
} // MDerror()

void debugtimer(const char *info) /****************************** debugtimer */
{
  static int step=0;
  static clock_t last;
  static double mylast;

  if (step==0) {
    fprintf(stderr,"\n\
timing (all times are in s):\n\
* MM = MC|MD\n\
* st = stride (slider simulation speed)\n\
* bk = block\n\
* speed = speed slider value (is recalculated to delay/stride)\n\
* step = counter\n\
* clock-prev = clock() - previous value of clock()\n\
* wall-prev = real (wall) time - previous value\n\
  with repeat_timeout: exp-like running average\n\
* to = timeout\n\
* avclk = averaged clock-prev over 1 block\n\
* abstime = tome (from 0) of repeat_timeout\n\
* wall/timer = wall_time / timer_value\n");
    fprintf(stderr,"MMst.bk speed step  clock-prev wall-prev  to|avclk  abstime wall/timer\n");
    mylast=mytime();
    last=clock(); }
  else if (step%debug==0) {
    clock_t t=clock();
    double myt=mytime(),cl=(double)(t-last)/((double)CLOCKS_PER_SEC*debug);
    static double avclock=0;

    avclock=avclock*0.999+cl*0.001;

    fprintf(stderr,"%s%2d.%02d %5.2f%6d  %9.6f %9.6f %9.6f %9.6f\n",
            info, speed.stride, block,
            sliders.speed->value(),
            step,
            cl,
            (myt-mylast)/debug,
            avclock,
            myt);
    last=t;
    mylast=myt; }
  step++;
}

void MDstep(int Pneeded) /******************************************* MDstep */
{
  vector ri,rt;
  double ff,maxvv,rr;
  double Vf;
  int i,j,k,l;

  if (method!=MDNPT) Pneeded=0;
  if (method==MTK) Pneeded=1;
  /* Pneeded: pressure is calculated here because it will be needed for a barostat:
     - MTK needs Ptensor calculated every step
     - MDNPT needs Ptensor in the last step of a sweep unless it is calculated
       in this step by meas() due to measure condition */

  nimg=ss.C2/L+1; // # of replicas if cutoff<Lh

  dtadjset();
  dt=dtadj;

  if (method==VICSEK) bc=PERIODIC; // should issue a warning...

  if (Pneeded) En.P=zero; // pressure tensor every step for MDNPT,MTK

  /* calculate forces (accelerations) caused by walls and gravity */
  loop (i,0,N)
    if (bc==PERIODIC) {
      a[i].x=a[i].y=0; }
    else {
      double fw,fwL;
      ri=r[i];
      fw=fwally(ri.y);
      fwL=fwallyL(L-ri.y);
      a[i].y=fw-fwL-gravity;
      // gravity not part of the virial
      if (Pneeded) En.P.y+=ri.y*fw+(L-ri.y)*fwL;
      if (bc==BOX) {
        fw=fwallx(ri.x);
        fwL=fwallxL(L-ri.x);
        a[i].x=fw-fwL;
       if (Pneeded) En.P.x+=ri.x*fw+(L-ri.x)*fwL; }
      else a[i].x=0; }

  if (method==VICSEK) {
    /* vllast used for the new velocity */
    memset(vllast,0,sizeof(vllast[0])*N);
    memset(wnbrs,0,sizeof(wnbrs[0])*N); }

  if (ss.C2>Lh) ss.ncell=1;

  /* calculate atom-atom pair forces; for MTK also the virial of force */
  if (ss.ncell>1) { // linked-cell list method
    DO_f *DO = Pneeded ? MDlcPneeded : MDlc;
    ss.Clc=ss.C2; // might be a tiny bit < C2 to neglect corners
    if (method==VICSEK) DO=MDlcVicsek;
    ss.lcfunc=f; // ugly patch, to be polished
    lc(DO); }

  else switch (bc) { // ss.ncell==1 means direct pair sum
    case BOX:
      loop (i,0,N) {
        ri=r[i];
        loop (j,i+1,N) {
          rt.x=ri.x-r[j].x; rt.y=ri.y-r[j].y;
          ff=f(Sqr(rt.x)+Sqr(rt.y));
          if (Pneeded) {
            En.P.x+=Sqr(rt.x)*ff;
            En.P.y+=Sqr(rt.y)*ff; }
          a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
          a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
      break;

    case SLIT:
      if (ss.C2<Lh)
        loop (i,0,N) {
          ri=r[i];
          loop (j,i+1,N) {
            rt.x=ri.x-r[j].x;
            if (rt.x>Lh) rt.x=rt.x-L;
            if (rt.x<-Lh) rt.x=rt.x+L;
            rt.y=ri.y-r[j].y;
            ff=f(Sqr(rt.x)+Sqr(rt.y));
            if (Pneeded) {
              En.P.x+=Sqr(rt.x)*ff;
              En.P.y+=Sqr(rt.y)*ff; }
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
                if (Pneeded) {
                  En.P.x+=Sqr(rt.x)*ff;
                  En.P.y+=Sqr(rt.y)*ff; }
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } }
      break;

    case PERIODIC:
      if (ss.C2<Lh)
        loop (i,0,N) {
          ri=r[i];
          loop (j,i+1,N) {
            rt.x=ri.x-r[j].x;
            if (rt.x>Lh) rt.x=rt.x-L;
            if (rt.x<-Lh) rt.x=rt.x+L;
            rt.y=ri.y-r[j].y;
            if (rt.y>Lh) rt.y=rt.y-L;
            if (rt.y<-Lh) rt.y=rt.y+L;
            rr=Sqr(rt.x)+Sqr(rt.y);
            if ( (ff=f(rr)) ) {
              if (Pneeded) {
                En.P.x+=Sqr(rt.x)*ff;
                En.P.y+=Sqr(rt.y)*ff; }
              if (method==VICSEK) {
                rr=ss.C2q-rr;
                vllast[i].x+=rr*v[j].x;
                vllast[i].y+=rr*v[j].y;
                vllast[j].x+=rr*v[i].x;
                vllast[j].y+=rr*v[i].y;
                wnbrs[i]+=rr;
                wnbrs[j]+=rr; }
              a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
              a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } } }
      else
        loop (i,0,N) {
          ri=r[i];
          loop (j,i+1,N) {
            loopto (k,-nimg,nimg) {
              rt.x=ri.x-r[j].x+L*k;
              loopto (l,-nimg,nimg) {
                rt.y=ri.y-r[j].y+L*l;
                ff=f(Sqr(rt.x)+Sqr(rt.y));
                if (Pneeded) {
                  En.P.x+=Sqr(rt.x)*ff;
                  En.P.y+=Sqr(rt.y)*ff; }
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; } } } } }

  //  fprintf(stderr,"wnbrs=%g %g %g ..\n",wnbrs[0],wnbrs[1],wnbrs[2]); exit(0);

  if (method==LANGEVIN || method==VICSEK) {
    /* random force by the fluctuation-dissipation theorem */
    double vv=sqrt(T/(tau*dt));

#define MIXPARM ss.C2q
    if (method==VICSEK) loop (i,0,N) {
      v[i].x=(vllast[i].x+v[i].x*MIXPARM)/(wnbrs[i]+MIXPARM);
      v[i].y=(vllast[i].y+v[i].y*MIXPARM)/(wnbrs[i]+MIXPARM);
      double rr=sqrt(T/(Sqr(v[i].x)+Sqr(v[i].y)));
      v[i].x*=rr;
      v[i].y*=rr; }

    loop (i,0,N) { a[i].x+=vv*rndgauss(); a[i].y+=vv*rndgauss(); } }

  if (method==NOSE_HOOVER || method==MTK) {
    /* velocity predictor TRVP(k=1,2)
       see J. Chem. Theory Comput. 2011, 7, 3596 [doi.org/10.1021/ct200108g]
       Codes @* below refer to the MACSIMUS code and manual */

#if TRVP==1
    Vf=xivpred=1.5*xiv-0.5*xivlast; // @4
#else
    Vf=xivpred=5./3*xiv-5./6*xivlast+xivllast/6.; // @4
#endif
    xivllast=xivlast; xivlast=xiv;

    if (method==MTK) {
#if TRVP==1
      lambdavpred=1.5*lambdav-0.5*lambdavlast; // @4
#else
      lambdavpred=5./3*lambdav-5./6*lambdavlast+lambdavllast/6.; // @4
#endif
      Vf+=(1+2./En.Nf)*lambdavpred; // @41
      lambdavllast=lambdavlast; lambdavlast=lambdav; }

    loop (i,0,N) {
#if TRVP==1
      a[i].x-=(1.5*v[i].x-0.5*vlast[i].x)*Vf; // @6
      a[i].y-=(1.5*v[i].y-0.5*vlast[i].y)*Vf; // @6
#else
      a[i].x-=(5./3*v[i].x-5./6*vlast[i].x+vllast[i].x/6.)*Vf; // @6
      a[i].y-=(5./3*v[i].y-5./6*vlast[i].y+vllast[i].y/6.)*Vf; // @6
#endif
      vllast[i]=vlast[i]; vlast[i]=v[i]; } }

  /* leap-frog integrator and kinetic energy */
  En.Ek=zero; // x,y separately (for diag P tensor)
              // En.Ekin = total kin. energy or Creutz's bag
  En.momentum=zero;
  maxvv=0;

  loop (i,0,N) {
    /* crash integrator test: this form detects NaN, too */
    //    fprintf(stderr,"%g %g\n",a[i].x,a[i].y);

    if (Sqr(a[i].x)+Sqr(a[i].y)>1e7) {
      if (debug)
        fprintf(stderr,"MDerror: i=%d r=[%g %g] v=[%g %g] a=[%g %g]\n",i,r[i].x,r[i].y,v[i].x,v[i].y,a[i].x,a[i].y);
      MDerror();
      return; }

    /* leap-frog Ekin components at t-dt/2 */
    En.Ek.x+=Sqr(v[i].x);
    En.Ek.y+=Sqr(v[i].y);

    /* leap-frog, momentum */
    En.momentum.x+=v[i].x+=dt*a[i].x;
    En.momentum.y+=v[i].y+=dt*a[i].y;

    /* leap-frog Ekin components at t+dt/2 */
    En.Ek.x+=Sqr(v[i].x);
    En.Ek.y+=Sqr(v[i].y);

    if (Sqr(v[i].x)+Sqr(v[i].y)>maxvv) maxvv=Sqr(v[i].x)+Sqr(v[i].y);

    r[i].x+=dt*v[i].x;
    r[i].y+=dt*v[i].y;

    //    if (r[i].x<0 || r[i].x>L || r[i].y<0 || r[i].y>L ) fprintf(stderr,"%8.4f %8.4f\n",r[i].x/L, r[i].y/L);
  } /*N*/

  En.momentum.x/=N; En.momentum.y/=N; // = averaged velocity, also for Vicsek
  // remove total momentum - does not have to be done every step...
  if (method>=NVE && method<=NOSE_HOOVER || method>=BUSSI && method<=MTK) {
    if (bc==PERIODIC) loop (i,0,N) {
        v[i].x-=En.momentum.x;
        v[i].y-=En.momentum.y; }
    else if (bc==SLIT) loop (i,0,N) {
        v[i].x-=En.momentum.x; } }

  if (En.Nf<=0) {
    fl_alert("\
No degree of freedom!\n\
(Periodic boundary conditions, N=1,\n\
and momentum conservation.)\n\
Switching to Metropolis Monte Carlo...");
    En.Nf=2;
    r[0].x=r[1].x=Lh;
    method=METROPOLIS; }

  /*
     En.Ek.x+En.Ek.y = 4*(total leap-frog kinetic energy), kin-corrected
     En.Ekin=total leap-frog kinetic energy
     see VERLET=3 in MACSIMUS
  */
  En.Ekin=Econserved=(En.Ek.x+En.Ek.y)/4; // total kinetic energy
  En.Ek.x*=En.q.x/2; // corrected for kinetic pressure tensor
  En.Ek.y*=En.q.y/2; // corrected for kinetic pressure tensor
  Tk=2*En.Ekin/En.Nf; // kinetic temperature, Tk=2*Ekin/Nf

  if (Pneeded) {
    /* En.Ek.x,En.Ek.y already scaled above, for MTK En.q.x=En.q.x=1 anyway */
    En.P.x=(En.Ek.x+En.P.x)/Sqr(L);
    En.P.y=(En.Ek.y+En.P.y)/Sqr(L); }

  double xikin,xipot,lakin; // mostly for debugging (but does not hurt)

  /* MTK and NOSE_HOOVER extended masses: */
  double M_T=En.Nf*T*Sqr(tau);
  double M_P=2*(En.Nf+2)*T*Sqr(tau*qtau);

  switch (method) {

    case MDNPT: {
      /* BERENDSEN barostat */
      double rho=L2rho(L);
      /* rough estimate of bulk modulus:
       * from T=3 the EOS is p=rho*T+rho^2+19*rho^5 */
      double B=rho*(T+rho*(1+50*rho*Sqr(rho)));

      ff=dt*((En.P.x+En.P.y)/2-P)/(B*tau*qtau); /* tau.P = tau*qtau, cf. MTK */
      //      fprintf(stderr,"MDNPT: L=%g rho=%g ncell=%d B=%g ff=%g\n",L,rho,ss.ncell,B,ff);

      BRACKETVAL(ff,-0.01,0.01);
      if (L2rho(L)<MINRHONPT) {
        errmessage=NPTFAILED;
        if (debug) savesim("debug.sim");
        method=METROPOLIS;
        d=0.3/Lh;
        ff=-0.1; /* shrink */ }

      ff=exp(ff);

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

      loop (i,1,En.Nf) Si+=sqr(rndgauss()); /* S_{Nf-1} (can be optimized) */

      ff=T/(En.Nf*Tk); // av(K)/Nf/K(t)
      ff=sqrt(c + (1-c)*(Rt*Rt+Si)*ff + 2*Rt*sqrt(ff*c*(1-c))); // alpha(t)
      loop (i,0,N) { v[i].x*=ff; v[i].y*=ff; } }
      // NB: change of sign not considered because improbablu
      break;

    case MAXWELL:
      /* Maxwell-Boltzmann thermostat: random hits */
      ff=sqrt(T);
      static double counter=0;
      counter-=dt;
      if (counter<0) {
        loop (i,0,N) randomv(i,ff);
        counter+=tau; }
      break;

    case ANDERSEN:
      /* Andersen thermostat: random hits */
      ff=sqrt(T);
      loop (i,0,N) if (rnd()<dt/tau) randomv(i,ff);
      break;

    case LANGEVIN:
    case VICSEK:
      /* friction part by the fluctuation-dissipation theorem
         - the simplest version giving the best Epot */
      ff=1-0.5*dt/tau;
      loop (i,0,N) { v[i].x*=ff; v[i].y*=ff; }
      break;

    case MTK:
      /* Pcfg calculated every step - finished here */

      xikin=M_T*Sqr(xiv/2); // 1/2 of leap-frog kin. energy of xi
      // 2*M_P was wrong (or M_P rescaled), see macsimus/history.txt
      xiv+=dt/M_T*(2*(En.Ekin+M_P*Sqr(lambdavpred))-(En.Nf+1)*T); // @A1
      xikin+=M_T*Sqr(xiv/2); // +1/2 of leap-frog kin. energy of xi
      xipot=(En.Nf+1)*T*xi; // xi=xi(t)
      xi+=dt*xiv; // xi=xi(t+dt) @B1
      Econserved+=xikin+xipot; // Ekin included, +U later...
      lakin=M_P*Sqr(lambdav)/2; // 1/2 of leap-frog kin. energy of lambda (*2 deg.f.)
      lambdav+=dt*((Sqr(L)*((En.P.x+En.P.y)/2-P)+Tk)/M_P-xivpred*lambdavpred); // @C1
      lakin+=M_P*Sqr(lambdav)/2; // 1/2 of leap-frog kin. energy of lambda (*2 deg.f.)
      Econserved+=lakin+P*Sqr(L);
      ff=exp(dt*lambdav); // equivalent to lambda+=dt*lambdav with L=exp(lambda) @D1
      if (L*ff>1e10) ff=1; // added in 12/2021
      L*=ff; Lh=L/2;
      sliders.rho->value(log(L2rho(L)));
      loop (i,0,N) { r[i].x*=ff; r[i].y*=ff; }
      // printf("Nf=%d mom=%g %g\n",En.Nf,momentum.x,momentum.y);
      // printf("%g %g %g %g %g %g %g\n",Econserved+Upot,Ekin,Upot,xikin,xipot,lakin,P*Sqr(L));
      break;

    case NOSE_HOOVER:
      xikin=M_T*Sqr(xiv/2);
      xiv+=dt*(Tk/T-1)/Sqr(tau);
      xikin+=M_T*Sqr(xiv/2);
      xipot=En.Nf*T*xi;
      Econserved+=xikin+xipot; // Ekin included, +U later...
      xi+=dt*xiv; // xi(t+dt)
      //      printf("%g %g %g %g %g\n",Econserved+Upot,Ekin,Upot,xikin,xipot);
      break;

    default:; /* to make the compiler happy */
  }

  t+=dt;

  if (debug<0) lc(debuglc);

  if (debug) debugtimer("MD");
} // MDstep()

void MCsweep(void) /************************************************ MCsweep */
{
  vector rt,ri,incirc,drt,dri; /* trial position, r[i], random vector */
  int i,j,k,l;
  double deltaU; /* MC: Utrial-U */
  int isdauto=setd->value();

  /* max displacement = Lh=L/2 */
  Lh=L/2;
  nimg=ss.C2/L+1; // # of replicas if cutoff>Lh

  En.Ekin=0; // = sum of Creutz's bag over N attempted steps

  loop (i,0,N) { // 1 sweep = N attempted MC moves
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
        if (ss.C2<Lh) {
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
        if (ss.C2<Lh) {
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

    if (method!=CREUTZ) bag=-T*log(rnd());

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

    En.Ekin+=bag;
  } /*N*/

  Tk=En.Ekin/N; // consistent with Tk=2*Ekin/Nf in MD (Nf=2N)

  if (isdauto) sliders.d->value(log(d)/DSCALE);

  if (method==MCNPT) {
    double Lf=rndcos()*dV,f=exp(Lf),rr;
    int isni=ss.C2<Lh && ss.C2<Lh*f;

    //    nimg=ss.C2/fmin(L,L*f)+1; // # of replicas if !isni

    ss.rrscale=f*f;
    ss.Clc=ss.C2*fmax(1,1/f); // max cutoff incl. rescaling
    nimg=ss.Clc/L+1; // # of replicas if !isni

    if (ss.Clc>Lh) ss.ncell=1;

    ss.DeltaU=0;

    if (bc!=PERIODIC)
      /* energy difference caused by the trial volume change: walls and gravity */
      loop (i,0,N) {
        ss.DeltaU+=
          uwally(r[i].y*f)-uwally(r[i].y)
          +uwallyL((L-r[i].y)*f)-uwallyL(L-r[i].y)
          -gravity*(r[i].y*(f-1));
        if (bc==BOX)
          ss.DeltaU+=
            uwallx(r[i].x*f)-uwallx(r[i].x)
            +uwallxL((L-r[i].x)*f)-uwallxL(L-r[i].x); }

    if (ss.ncell>1) {
      DO_f *DO = MCNPTlc;
      ss.lcfunc=u; // function called inside MCNPTlc
      lc(DO); }

    else loop (i,0,N) {

      switch (bc) {
        case BOX:
          loop (j,0,i) {
            rr=Sqr(r[i].x-r[j].x)+Sqr(r[i].y-r[j].y);
            ss.DeltaU+=u(rr*ss.rrscale)-u(rr); }
          break;
        case SLIT:
          if (isni)
            loop (j,0,i) {
              rr=Sqrni(r[i].x-r[j].x)+Sqr(r[i].y-r[j].y);
              ss.DeltaU+=u(rr*ss.rrscale)-u(rr); }
          else
            loop (j,0,i)
              loopto (k,-nimg,nimg) {
                drt.x=r[i].x-r[j].x+k*L;
                drt.y=r[i].y-r[j].y;
                rr=Sqr(drt.x)+Sqr(drt.y);
                ss.DeltaU+=u(rr*ss.rrscale)-u(rr); }
          break;
        case PERIODIC:
          if (isni)
            loop (j,0,i) {
              rr=Sqrni(r[i].x-r[j].x)+Sqrni(r[i].y-r[j].y);
              ss.DeltaU+=u(rr*ss.rrscale)-u(rr); }
          else
            loop (j,0,i)
              loopto (k,-nimg,nimg) {
                drt.x=r[i].x-r[j].x+k*L;
                loopto (l,-nimg,nimg) {
                  drt.y=r[i].y-r[j].y+l*L;
                  rr=Sqr(drt.x)+Sqr(drt.y);
                  ss.DeltaU+=u(rr*ss.rrscale)-u(rr); } }
      } /* bc */
    }

    if (rnd()<exp(Lf*(2*N)-(P*Sqr(L)*(ss.rrscale-1)+ss.DeltaU)/T)) {
      /*
         This is for NPT ensemble with volume measure mu=1/V,
         which is compatible with MTK (with conserved momentum).
         See https://doi.org/10.1063/5.0193281 for details.
         (For mu=const (e.g., mu=P/T), use Lf*(2*N+2) instead of Lf*(2*N)
         in the formula above.)
      */
      Vaccepted++;
      if (isdauto) {
        if (dV<MAXDV) dV*=1+(1-acc)*0.01; }
      rr=L2rho(L);
      if (rr<MINRHONPT) {
        errmessage=NPTFAILED;
        f=0.8; /* shrink box */
        rr/=f*f;
        if (debug) savesim("debug.sim");
        method=METROPOLIS;
        d=0.3/Lh; }
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
} // MCsweep()

void art(void) /******************************************************** art */
{
  double phase=t/trace;
  int i;

  if (isMD(method)) phase/=dt;

  loop (i,0,N) {
    phase=fmod(phase+2*PI*i/N,2*PI);

    molcol[i]=((unsigned)(cos(phase)*127.9999+128)<<24)
            | ((unsigned)(cos(phase+2*PI/3)*127.9999+128)<<16)
            | ((unsigned)(cos(phase+4*PI/3)*127.9999+128)<<8); }
}

int shownbrs;

void neighbors(void) /******************************************** neighbors */
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
          if (Sqr(rt.x-r[j].x)+Sqr(rt.y-r[j].y)<ss.rnbrq) nbr++;
        break;
      case SLIT:
        loop (j,0,N) if (i!=j)
          if (Sqrni(rt.x-r[j].x)+Sqr(rt.y-r[j].y)<ss.rnbrq) nbr++;
        break;
      case PERIODIC:
        loop (j,0,N) if (i!=j)
          if (Sqrni(rt.x-r[j].x)+Sqrni(rt.y-r[j].y)<ss.rnbrq) nbr++;
        break; }

    if (nbr>7) nbr=7;
    molcol[i]=nbrcol[nbr]; }
} // neighbors()

// return the last single extension .ext (last .)
// FILE. => zero-length string
// BUG: \ is treated as directory separator in linux

const char *getext(const char *fn) /********************************* getext */
{
  const char *c,*d;

  for (c=fn,d=NULL; *c; c++) {
    if (*c=='.') d=c;
    // I do not know how the directory separator works over platforms...
    if (*c=='/') d=NULL;
    if (*c=='\\') d=NULL; }

  return d;
} // getext()

void loadsim(const char *fn) /************************************** loadsim */
{
  FILE *in=fl_fopen(fn,"rt");
  char line[256];
  int i;

  if (in) {
    if (fgets(line,256,in)) {
      if (!strstr(line,SIMVARSQ)) {
        fl_alert("read \"%s\": is not a SIMOLANT file",fn);
        fclose(in);
        return; }
      sscanf(line,"%d %d %d %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf",
             // see macro SIMVARS for the list of variables
             &N,(int*)&bc,&walls,&T,&L,&walldens,&dt,&d,(int*)&method,&P,&dV,&tau,&qtau,&gravity,&ss.a,&ss.b,&ss.C2); }
    if (!fgets(line,256,in)) {
      fl_alert("read \"%s\": bad header or bad format",fn);
      fclose(in);
      return; }

    put_data();

    loop (i,0,N)
      if (fgets(line,256,in))
        sscanf(line,"%lf %lf %lf %lf",&r[i].x,&r[i].y,&v[i].x,&v[i].y);
      else {
        fl_alert("read \"%s\": too short file or bad format",fn);
        fclose(in);
        return ; }

    fclose(in); }
  else
    fl_alert("file \"%s\" unreadable or not found",fn);
} // loadsim()

/* mean squared displacement, following the periodic b.c. if applicable */
void MSD(vector *MSDp) /************************************************ MSD */
{
  double D,maxD=0;
  int i;

  /* enum bctype {BOX,SLIT,PERIODIC}
     NB: because of compatibility with NPT,
     the b.c. calculations are done with the box rescaled to [1,1] */

  if (bc>=SLIT)
    loop (i,0,N) {
      D=r[i].x/L-rpbc[i].x;
      //    fprintf(stderr,"%d dx=%g ",i,D);
      rpbc[i].x=r[i].x/L;
      while (D>0.5) {
        rpbc[i].x-=1; D-=1; }
      while (D<-0.5) {
        rpbc[i].x+=1; D+=1; }
      if (fabs(D)>maxD) maxD=fabs(D);
      MSDp->x+=Sqr(rpbc[i].x-r0[i].x); }
  else
    loop (i,0,N) {
      rpbc[i].x=r[i].x/L; // in fact, need not be stored
      MSDp->x+=Sqr(rpbc[i].x-r0[i].x); }
  MSDp->x*=Sqr(L)/N;

  if (bc==PERIODIC)
    loop (i,0,N) {
      D=r[i].y/L-rpbc[i].y;
      //    fprintf(stderr,"%d dy=%g\n",i,D);
      rpbc[i].y=r[i].y/L;
      while (D>0.5) {
        rpbc[i].y-=1; D-=1; }
      while (D<-0.5) {
        rpbc[i].y+=1; D+=1; }
      if (fabs(D)>maxD) maxD=fabs(D);
      MSDp->y+=Sqr(rpbc[i].y-r0[i].y); }
  else
    loop (i,0,N) {
      rpbc[i].y=r[i].y/L; // in fact, need not be stored
      MSDp->y+=Sqr(rpbc[i].y-r0[i].y); }
  MSDp->y*=Sqr(L)/N;

  if (!MSDskip && maxD>0.45) {
    errmessage=MSDFAILED;
    MSDskip=1; }
} // MSD()
