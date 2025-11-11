// fltk-config --use-images --compile simolant.cc
// g++ -O2 -o simolant simolant.cc -lfltk -lfltk_images
// g++ -g -o simolant simolant.cc -lfltk -lfltk_images

#define VERSION "11/2025"

// cmath must be first (otherwise windows problems)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <unistd.h>
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Box.H>
#include <FL/fl_ask.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Help_Dialog.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Fill_Slider.H>
#include <FL/Fl_Native_File_Chooser.H>

#include "include.h"

#define PI M_PI // needed in rndetc.c (#included from rndgeni.c)
#include "rndgeni.c"

#define PANELW 480       // right panel width (needs tinkering if changed)
#define MENUH 24         // top menu height (reasonable range 20..30)
#define INFOPANELH 92    // height of the top panel with N=, ensemble, parameters
                         // given by font size - do not change
#define RESULTPANELH 250 // height of the results panel and graphs
                         // given by font size and graph - do not change
#define BUTTONPANELH 370 // height of the bottom panel with buttons,sliders
                         // partly flexible, but cannot shrink too much
#define BORDER 7         // box border (can change)
// square box size (initial, wall-wall area, can change):
#define BOXSIZE (INFOPANELH+RESULTPANELH+BUTTONPANELH-2*BORDER)
// default geometry
#define INITGEOMETRY_Y (INFOPANELH+RESULTPANELH+BUTTONPANELH+MENUH)
#define INITGEOMETRY_X (INFOPANELH+RESULTPANELH+BUTTONPANELH+PANELW)
#define MAXNCELL 34

// Okabe and Ito color-blind-safe palette
// https://www.nceas.ucsb.edu/sites/default/files/2022-06/Colorblind%20Safe%20Color%20Schemes.pdf

#define OI_BLACK  0x00000000
#define OI_GREEN  0x009E7300
#define OI_BLUE   0x0072B200
#define OI_CYAN   0x56B4E900
#define OI_YELLOW 0xF0E44200
#define OI_ORANGE 0xE69F0000
#define OI_RED    0xD55E0000
#define OI_LILAC  0xCC79A700
#define LIGHTYELLOW 0xFFF08000 // for graphs vs. OI_GREEN

// force field stuff moved here because of ss.rnbr used in panel
#define RNBR 1.55 // range for neighbors, good for 2D LJ
#define CUTOFF 4 // default cutoff (applies for 2D LJ)

// parameter value normalized to range [FROM,TO]
#define BRACKETVAL(VAL,FROM,TO) do { \
  if (VAL<FROM) VAL=FROM; \
  if (VAL>TO) VAL=TO; } while(0)

enum ff_e {LJ,WCA,PD,DW,NFF};
/*
   LJ  = 2D Lennard-Jones, 4/r^8-4/r^4, smoothly truncated
   WCA = WCA version of 2D Lennard-Jones
   PD  = Penetrable Disks a*(r^2-c^)^2
   DW  = Double Well (Penrose)
   NFF = undef, sentinel
*/

typedef double (*function_t)(double r); // for potential() and forces()

struct ss_s {
  double C1;        // smooth cutoff from C1 to C2=c
  double C2=CUTOFF; // cutoff=c, from cmd: and -P
  double C1q;       // C1^2
  double C2q;       // C2^2
  double A,A4;      // constants in the smoothed part
  double a=1.25;    // a from cmd: and -P
  double aq;        // a^2
  double b=1.4;     // b from cmd: and -P
  double bq;        // b^2
  double rnbr=RNBR; // thereshold to count neighbors
  double rnbrq=Sqr(RNBR);
  double P[5];      // DW forces expansion, see Penrose.mw - two minima
  double B2;        // 2nd virial coefficient for given thermostat parameter T
  double Clc;       // linked-cell list: see MCNPT (might be a tiny bit < C2)
  double rrscale;   // for MCNPT
  double DeltaU;    // for MCNPT
  function_t lcfunc;// linked-cell list: function passed to lc() 
  enum ff_e ff=NFF; // force field
  int ncell=0;      // for linked-cell list
                    // ncell=0 automatic setup
                    // ncell=1 direct pair sum
                    // ncell>1 linke-cell list, ss.Clc<=L/2 required
  int autoncell=1;  // 1 if ncell should be automatically determined
} ss;

static
int autoncell(int N) /******************************************** sutoncell */
{
  if (N<57) return 1;
  else return sqrt(N)/3+0.5;
}

// some global variables collected here
struct files_s {
  int record=0;         // state of recording
  int record0=0;        // last state of recording
  int nmeas=0;          // counter of blocked measurements
  char sep=',';         // CSV field separator (global)
  FILE *csv=NULL;       // output CSV file with blocked neasurements
  char protocol[256]="simolant"; // .txt or .csv will be added, cf. cb_protocol()
  char *ext=protocol+8; // . of the extension
  char info[1024]="";   // blocked measurement to record
} files;

// KM: menu-redraw must be explicitly requested, so Simulate must know it

struct Sliders {
  Fl_Fill_Slider *T, *tau,*d, *g,*c,  *rho,*P, *N,
    *speed,*block;
} sliders;

struct Parms_s {
  Fl_Button *plus,*minus,*gzero;
} Parms;

struct Buttons {
  Fl_Light_Button *wallx,*wally,*wallxL,*wallyL;
  Fl_Button *shift_up,*shift_down,*shift_left,*shift_right;
  Fl_Button *invertwalls;
  Fl_Light_Button *record,*csv,*comma;
  Fl_Input *cmd;
} buttons;

Fl_Help_Dialog *help;
Fl_Choice *molsize,*drawmode,*colormode,*incl;
Fl_Button *resetgraph,*resetview,*setd;
Fl_Light_Button *run; // run/hold

void calculateB2(void); // needed in setss()

// replacements for deprecated fl_ask
static int fl_ask1(const char *fmt,const char *arg) /*************** fl_ask1 */
{
  return !fl_choice(fmt,"Yes","No",0,arg);
} // fl_ask1()

/* MESS IN DECLARATIONS -- to be cleaned */
typedef struct {
  double x,y;
} vector;

typedef struct linklist_s { /* 64 B = cache line */
  vector r;      /* site position physically here (periodic b.c. optimizing) */
  vector *f;     /* forces (indirected - see MACSIMUS for other options) */
  vector *v;     /* original velocities for Vicsek */
  vector *vnew;  /* new velocities for Vicsek = vllast+i */
  double *wnbrs; /* for Vicsek */
  int i;         /* atom number (debug only) */
  int padding;   /* sizeof(linklist_t)=64 (padded anyway on 64b architecture) */
  struct linklist_s *next;
} linklist_t;

typedef void DO_f(linklist_t *l);
DO_f debuglc;
void lc(DO_f DO);

#include "speed.cc"
#include "calculate.cc"
#include "linklist.cc"
#include "callbacks.cc"
#include "record.cc"

Fl_Menu_Item molsizeitems[]={
  { "&Real" },
  { "&Small" },
  { "&Dot" },
  { "&Pixel" },
  { 0 } };
Fl_Menu_Item drawmodeitems[]={
  { "&Movie" },
  { "&Traces" },
  { "&Lines" },
  { "&Nothing" },
  { 0 } };
Fl_Menu_Item colormodeitems[]={
  { "&Black" },
  { "&One" },
  { "&y-split" },
  { "&Neighbors" },
  { "&Random" },
  { "&Art" },
  { "&Keep" },
  { 0 } };

Fl_Menu_Item inclitems[]={
  { "&Nothing" },
  { "&Convergence profile" },
  { "&Density profile" },
  { "&Both" },
  { 0 } };


#include "panel.cc"
Fl_Menu_Bar *menu;
#include "simulate.cc"

Fl_Menu_Item menuitems[] = {
  // NAME KEYBOARD_SHORTCUT CALLBACK *USERDATA
  { "&File ", FL_ALT|'f', 0, 0, FL_SUBMENU },
  { "&Protocol name...", 0, cb_protocol, 0, FL_MENU_DIVIDER },
  { "&Load configuration...", 0, cb_load },
  { "&Save configuration as...", 0, cb_save, 0, FL_MENU_DIVIDER },
  { "&Export force field", 0, cb_ff, 0, FL_MENU_DIVIDER },
  { "&Quit", 0, cb_exit },
  { 0 },
  { "&Prepare system  ", FL_ALT|'p', 0, 0, FL_SUBMENU },
  { "&Gas in box", 0, cb_gas },
  { "+ &Diffusion", 0, cb_diffusion },
  { "+ G&ravity", 0, cb_gravity, 0, FL_MENU_DIVIDER },
  { "Vapor-liquid e&quilibrium", 0, cb_vle },
  { "Horizontal &Slab", 0, cb_slab },
  { "&Nucleation", 0, cb_nucleation },
  { "&Liquid droplet", 0, cb_liquid },
  { "&Two droplets", 0, cb_twodrops },
  { "&Bubble (cavity)", 0, cb_cavity },
  { "&Capillary", 0, cb_capillary },
  { "&Periodic liquid", 0, cb_periodicliquid, 0, FL_MENU_DIVIDER },
  { "&Hexagonal crystal", 0, cb_crystal },
  { "+ &Edge defect", 0, cb_defect },
  { "+ &Vacancy", 0, cb_vacancy },
  { "+ &Intersticial", 0, cb_intersticial },
  { "Pe&riodic crystal", 0, cb_periodiccrystal, 0, FL_MENU_DIVIDER },
  { "Vicse&k active matter", 0, cb_vicsek },
  { 0 },
  { "F&orce field  ", FL_ALT|'o', 0, 0, FL_SUBMENU },
  { "&Lennard-Jones c=4 (recommended)", 0, cb_LJ, 0, FL_MENU_DIVIDER},
  { "&Repulsive (WCALJ, c=1.1892)", 0, cb_WCA},
  { "&Penetrable disks (a=1,c=2)", 0, cb_PD},
  { "&Ideal gas (a=0,c=2)", 0, cb_IG},
  { "&Attractive penetrable disks (a=-1,c=2)", 0, cb_APD},
  { "&Double well (a=1.25,b=1.4,c=-2)", 0, cb_DW },
  { 0 },
  { "&Method  ", FL_ALT|'m', 0, 0, FL_SUBMENU },
  //UNICODE problem { "&Automatic (MC→MD/Bussi CSVR)", 0, cb_AUTO, 0, FL_MENU_DIVIDER },
  { "&Automatic (MC->MD/Bussi CSVR)", 0, cb_AUTO, 0, FL_MENU_DIVIDER },
  { "Monte Carlo NVE (&Creutz)", 0, cb_MC_C },
  { "Monte Carlo  NVT (&Metropolis)", 0, cb_MC },
  { "Monte Carlo   N&PT (Metropolis)", 0, cb_MC_NPT, 0, FL_MENU_DIVIDER },
  { "Molecular Dynamics NV&E (microcanonical)", 0, cb_MD },
  { "Molecular Dynamics  NVT (&Berendsen)", 0, cb_MD_B },
  { "Molecular Dynamics  NVT (&Nosé-Hoover)", 0, cb_MD_NH },
  { "Molecular Dynamics  NVT (An&dersen)", 0, cb_MD_A },
  { "Molecular Dynamics  NVT (Ma&xwell-Boltzmann)", 0, cb_MD_M },
  { "Molecular Dynamics  NVT (&Langevin)", 0, cb_MD_L },
  { "Molecular Dynamics  NVT (B&ussi CSVR)", 0, cb_MD_BUSSI },
  { "Molecular Dynamics   NPT (Be&rendsen)", 0, cb_MD_MDNPT },
  { "Molecular Dynamics   NPT (Mar&tyna et al.)", 0, cb_MD_MTK,0, FL_MENU_DIVIDER },
  { "Vicse&k with NVT (Langevin)", 0, cb_VICSEK },
  { 0 },
  { "&Boundary conditions  ", FL_ALT|'b', 0, 0, FL_SUBMENU },
  { "&Box (soft walls)", 0, cb_box },
  { "&Slit (x-periodic, y-walls)", 0, cb_slit },
  { "&Periodic", 0, cb_periodic },
  { 0 },
  { "&Show  ", FL_ALT|'s', 0, 0, FL_SUBMENU },
  { "&Quantities", 0, cb_quantities, 0, FL_MENU_DIVIDER },
  { "&Energy/enthalpy convergence profile", 0, cb_energy },
  { "&Temperature convergence profile", 0, cb_temperature },
  { "&Pressure convergence profile", 0, cb_pressure },
  { "&Volume convergence profile", 0, cb_volume },
  { "&Momentum (e.g., for Vicsek)", 0, cb_momentum },
  { "&Integral of motion convergence profile", 0, cb_intmotion, 0, FL_MENU_DIVIDER },
  { "&Radial distribution function", 0, cb_rdf },
  { "Vertical densit&y profile", 0, cb_zprofile },
  { "&Droplet radial density profile", 0, cb_dprofile },
  { "&Cavity radial density profile", 0, cb_cprofile },
  { 0 },
  { "&Tinker  ", FL_ALT|'t', 0, 0, FL_SUBMENU },
  { "FPS=&15", 0, cb_fps15 },
  { "FPS=&30", 0, cb_fps30 },
  { "FPS=&60", 0, cb_fps60 },
  { "FPS=12&0", 0, cb_fps120 , 0, FL_MENU_DIVIDER },
  { "Timeout &99% of simulation (use if shaky)", 0, cb_timeout99 },
  { "Timeout &50% of simulation", 0, cb_timeout50 },
  { "Timeout &20% of simulation (faster)", 0, cb_timeout20 },
  { "Timeout &7% of simulation (fastest, may stuck)", 0, cb_timeout7 },
  { 0 },
  { "&Help", FL_ALT|'h', 0, 0, FL_SUBMENU  },
  { "&Manual", 0, cb_help, 0, FL_MENU_DIVIDER },
  { "&About", 0, cb_about },
  { 0 },
  { 0 } };

int main(int narg,char **arg) /**************************************** main */
{
  Simulate *simulate;
  typedef char *charptr;
  char **newarg=new charptr[narg];
  const char *fn="",*ext;
  int iarg,newiarg=0,fltkoptions=0,SEED=0;
  enum cfg_t init=GAS;
  char *optionP=NULL;

  //  fprintf(stderr,"sizeof(linklist_t)=%ld\n",sizeof(linklist_t));
  
  if (narg==2 && arg[1][0]!='-') fn=arg[1]; // single argument
  else for (iarg=0; iarg<narg; iarg++) {
    if (arg[iarg][0]=='-' && isupper(arg[iarg][1]))
      switch (arg[iarg][1]) {
        case 'D': debug=atoi(arg[iarg]+2); if (!debug) debug=1; break;
        case 'F': speed.FPS=atof(arg[iarg]+2); break;
        case 'I': fn=arg[iarg]+2; init=(enum cfg_t)atoi(arg[iarg]+2); break;
        case 'N': N=atoi(arg[iarg]+2); break; // legacy
        case 'O':
          if (strlen(files.protocol)>250) exit(4);
          strcpy(files.protocol,arg[iarg]+2);
          files.ext=files.protocol+strlen(files.protocol);
          break;
        case 'P': optionP=arg[iarg]+2; break;
        case 'T': speed.MINTIMEOUT=atof(arg[iarg]+2); break;
        case 'U': speed.extradelay=(int)(atof(arg[iarg]+2)*1e6+0.5); break;
        case 'Z': SEED=atoi(arg[iarg]+2); break;
        default: fprintf(stderr,"\
SIMOLANT " VERSION " options (NO SPACE BEFORE NUMBER):\n\
  <FILE>     (NAME.sim) open sim-file on start (single argument only)\n\
  -D<CYCLES> debug mode (timing and technical info printed every CYCLES)\n\
             disabled in Windows (recompile with -mconsole)\n\
  -D         the same as -D1\n\
  -D-1       debug linked-cell list method\n\
  -h         get help on lowercase FLTK options\n\
  -I<INIT>   start, see menu Prepare system for numbers [default=0=Gas]\n\
  -N<N>      number of particles (the same as -PN=<N>)\n\
  -P<LIST>   comma-separated parameter list, the same as for input field \"cmd:\"\n\
             performed after -I, -N\n\
  -U<DELAY>  delays in s added, use if the widgets are lazy [default=%g]\n\
  -Z<SEED>   random number seed; default=0=time\n\
In case of multiple instances of the same option, the last one applies.\n\
FLTK options are in lowercase and with a space before a number, get them by:\n\
  simolant -x\n\
Do not change too much the default geometry -g %dx%d\n\
EXAMPLE:\n\
  simolant -g 1150x700 -I7 -N600 -Pg=-0.02,tau=0.5\n\
",speed.extradelay,INITGEOMETRY_X,INITGEOMETRY_Y);
          if (!fltkoptions) exit (0);
      }
    else {
      newarg[newiarg++]=arg[iarg]; } }

  if (N<1) N=1;
  if (N>=MAXN) N=MAXN;

  rndinit(0,SEED);
  //  Fl::set_color(FL_BACKGROUND_COLOR,0x05558500);

  Fl_Double_Window win(BOXSIZE+2*BORDER+PANELW,
                       BOXSIZE+2*BORDER+MENUH,"SIMOLANT " VERSION);
  win.color(FL_WHITE);

  if (debug) fprintf(stderr,"geometry:\n  -g %dx%d\n",BOXSIZE+2*BORDER+PANELW,BOXSIZE+2*BORDER+MENUH);

  // the upper menu
  menu=new Fl_Menu_Bar(0,0,BOXSIZE+2*BORDER+PANELW,MENUH);
  menu->menu(menuitems);

  // right panel with info, graphs, and many widgets
  panel = new Panel(BOXSIZE+2*BORDER, MENUH, PANELW, BOXSIZE+2*BORDER);

  /*
    Simulate window is the simulation box to which the configuration is drawn.
    The draw() method (actually accessed as redraw() which "schedules
    redrawing") contains also drawing all panels and graphs and MC/MD
    calculations.  At the same time, it initiates the timer_callback() loop.
    with redrawing everything and refreshing the timer.
    The simulation box is resizable, the right panel is resized accordingly,
    cf. the note in panel.cc: Panel()
  */
  simulate = new Simulate(0, MENUH, BOXSIZE+2*BORDER, BOXSIZE+2*BORDER);
  win.resizable(simulate);

  help = new Fl_Help_Dialog();
  help->textsize(16);
  help->resize(0,0,800,600);

  if (init<GAS || init>=SENTINEL) init=GAS;

  Fl::visual(FL_DOUBLE|FL_INDEX);

  if (newiarg>1) fprintf(stderr,"FLTK ");
  win.show(newiarg,newarg); // only FLTK options passed
  loop (iarg,1,newiarg) if (!memcmp(newarg[iarg],"-h",2)) {
    fprintf(stderr,"Use option -H to get simolant options\n");
    exit(0); }

  if (ss.ff==NFF) ss.ff=LJ; // this is because INITVICSEK sets ss.ff=PD
  setss(ss.ff,CUTOFF); // if -P contains T, calculateB2 must have parameters set

  ext=getext(fn);
  if (ext && !strcmp(ext,".sim")) loadsim(fn); // name.sim => load it
  else initcfg(init); // no .ext=.sim => is a number

  if (optionP) {
    // reading option -P: decimal separator=period, command separator=comma
    char *c,*e;

    optionP=strdup(optionP);
    for (c=optionP; ; ) {
      e=strchr(c,',');
      parse_cmd(c);
      if (e) c=e+1;
      else break; }

    free(optionP); }

  setss(ss.ff,ss.C2); // in case parameters have changed

  return Fl::run();
}
