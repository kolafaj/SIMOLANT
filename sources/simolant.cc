// fltk-config --use-images --compile simolant.cc
// g++ -O2 -o simolant simolant.cc -lfltk -lfltk_images

#define VERSION "02/2024"

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

#if FL_MAJOR_VERSION==1 && FL_MINOR_VERSION<3
// old version:
//   no <FL/Fl_Native_File_Chooser.H>
//   UTF8 characters not treated in V<1.3 (ugly anyway)
#  define fl_fopen fopen
#  include <FL/Fl_File_Chooser.H>
#else
#  define NATIVE
#  include <FL/Fl_Native_File_Chooser.H>
#endif

#include "include.h"

#define PI M_PI // needed in rndetc.c
#include "rndgeni.c"

#define PANELW 329 // right panel width (needs tinkering if changed)
#define MENUH 25  // top menu height (reasonable range 20..30)
#define INFOPANELH 92 // height of the top panel with N=, ensemble, parameters
                      // given by font size - do not change
#define RESULTPANELH 250 // height of the results panel and graphs
                         // given by font size and graph - do not change
#define BUTTONPANELH 370 // height of the bottom panel with buttons,sliders
                         // partly flexible, but cannot shrink too much
#define BORDER 7 // box border (can change)
// square box size (initial, wall-wall area, can change):
#define BOXSIZE (INFOPANELH+RESULTPANELH+BUTTONPANELH-2*BORDER)

static double TIMERDELAY=0.05; // timer default value in s, option -T
static double timerdelay=TIMERDELAY; // to be changed by sliders.speed
static double initspeed=1.1; // nit=3
static int extradelay=10000; // in us - delay added to avoid Xorg jam, see -U
//static double initspeed=0.7; // nit=2
#define RNBR 1.55 // range for neighbors

//KM menu-redraw must be explicitly requested, so Simulate must know it

struct Sliders {
  Fl_Fill_Slider *T,*tau,*d,*g,*rho,*P,*N, *speed,*block;
} sliders;

struct Buttons {
  Fl_Light_Button *wallx,*wally,*wallxL,*wallyL, *meas,*comma;
  Fl_Input *parms;
} buttons;

Fl_Help_Dialog *help;
Fl_Choice *molsize,*drawmode,*colormode,*incl;
Fl_Button *resetchoices,*setd;
Fl_Light_Button *run; // run/hold
char *protocol=(char*)"simolant.txt";
#ifdef INCL
int writeCP=0; // full convergence profile will be written, changed to incl->value()&1
int writeDF=0; // density profiles, RDF will be written, changed to incl->value()&2
#endif

#include "calculate.cc"
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
  { "&y-split" },
  { "&Neighbors" },
  { "&Random" },
  { "&Keep" },
  { 0 } };

Fl_Menu_Item inclitems[]={
  { "&Nothing" },
  { "&Conv.prof." },
  { "&Dens.prof." },
  { "&Both" },
  { 0 } };

#include "panel.cc"
Fl_Menu_Bar *menu;
#include "simulate.cc"

Fl_Menu_Item menuitems[] = {
  { "&File", FL_ALT|'f', 0, 0, FL_SUBMENU },
  { "&Open configuration...", 0, cb_open },
  { "&Save configuration as...", 0, cb_save, 0, FL_MENU_DIVIDER },
  { "&Protocol name...", 0, cb_record, 0, FL_MENU_DIVIDER },
#ifdef INCL  
  { "Include &Convergence profile...", 0, cb_cp },
  { "Include &Density functions...", 0, cb_df, 0, FL_MENU_DIVIDER },
#endif
  { "&Quit", 0, cb_exit },
  { 0 },
  { "&Prepare system  ", FL_ALT|'p', 0, 0, FL_SUBMENU },
  { "&Gas in box", 0, cb_gas },
  { "+ &Diffusion", 0, cb_diffusion },
  { "+ G&ravity", 0, cb_gravity, 0, FL_MENU_DIVIDER },
  { "Vapor-liquid e&quilibrium", 0, cb_vle },
  { "Horizontal &Slab", 0, cb_slab },
  { "&Nucleation", 0, cb_nucleation, 0, FL_MENU_DIVIDER },
  { "&Liquid droplet", 0, cb_liquid },
  { "&Capillary", 0, cb_capillary },
  { "&Bubble (cavity)", 0, cb_cavity, 0, FL_MENU_DIVIDER },
  { "&Hexagonal crystal", 0, cb_crystal },
  { "+ &Edge defect", 0, cb_defect },
  { "+ &Vacancy", 0, cb_vacancy },
  { "+ &Intersticial", 0, cb_intersticial },
  { 0 },
  { "&Method  ", FL_ALT|'m', 0, 0, FL_SUBMENU },
  { "&Automatic (MC→MD/Bussi CSVR)", 0, cb_AUTO, 0, FL_MENU_DIVIDER },
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
  { "Molecular Dynamics   NPT (B&erendsen)", 0, cb_MD_MDNPT },
  { "Molecular Dynamics   NPT (Ma&rtyna et al.)", 0, cb_MD_MTK },
  { 0 },
  { "&Boundary conditions  ", FL_ALT|'b', 0, 0, FL_SUBMENU },
  { "&Box (soft walls)", 0, cb_box },
  { "&Slit (x-periodic, y-walls)", 0, cb_slit },
  { "&Periodic", 0, cb_periodic },
  { 0 },
  { "&Show  ", FL_ALT|'s', 0, 0, FL_SUBMENU },
  { "&Minimum", 0, cb_nomeas },
  { "&Quantities", 0, cb_quantities, 0, FL_MENU_DIVIDER },
  { "&Energy/enthalpy convergence profile", 0, cb_energy },
  { "&Temperature convergence profile", 0, cb_temperature },
  { "&Pressure convergence profile", 0, cb_pressure },
  { "&Volume convergence profile", 0, cb_volume },
  { "&Integral of motion convergence profile", 0, cb_intmotion, 0, FL_MENU_DIVIDER },
  { "&Radial distribution function", 0, cb_rdf },
  { "Vertical densit&y profile", 0, cb_zprofile },
  { "&Droplet radial density profile", 0, cb_dprofile },
  { "&Cavity radial density profile", 0, cb_cprofile },
  { 0 },
  { "&Help", FL_ALT|'h', 0, 0, FL_SUBMENU  },
  { "&Manual", 0, cb_help, 0, FL_MENU_DIVIDER },
  { "&About", 0, cb_about },
  { 0 },
  { 0 } };

int main(int narg,char **arg) {
  Simulate *simulate;
  typedef char *charptr;
  char **newarg=new charptr[narg];
  const char *fn="";
  int iarg,newiarg=0;
  enum cfg_t init=GAS;
  char *optionP=NULL;

  if (narg==2 && arg[1][0]!='-') fn=arg[1]; // single argument
  else for (iarg=0; iarg<narg; iarg++) {
    if (arg[iarg][0]=='-' && isupper(arg[iarg][1]))
      switch (arg[iarg][1]) {
        case 'C': circlemethod=atoi(arg[iarg]+2); break;
        case 'D': debug=atoi(arg[iarg]+2); if (!debug) debug=100; break;
        case 'I': fn=arg[iarg]+2; init=(enum cfg_t)atoi(arg[iarg]+2); break;
        case 'N': N=atoi(arg[iarg]+2); break;
        case 'P': optionP=arg[iarg]+2; break;
        case 'S': initspeed=atof(arg[iarg]+2); break;
        case 'T': TIMERDELAY=atof(arg[iarg]+2); break;
        case 'U': extradelay=atoi(arg[iarg]+2)*1000; break;
        default: fprintf(stderr,"\
SIMOLANT " VERSION " options (NO SPACE BEFORE NUMBER):\n\
  -C<CIRCLE>   circle drawing method:\n\
               0=fl_pie  1=fl_circle  2=custom of fl_line [default=2]\n\
  -D<MDSTEPS>  debug mode (timing and technical info printed every MDSTEPS)\n\
               disabled in Windows (recompile with -mconsole)\n\
  -I<INIT>     (number) start with system, see the Prepare system menu for\n\
               numbers [default=0=Gas]\n\
  -I<FILE>     (name.sim) open sim-file on start\n\
  -N<N>        initial number of particles\n\
  -P<LIST>     comma-separated parameter list, performed after -I, -N,\n\
               the same as for input field \"cmd:\"\n\
  -S<SPEED>    initial value of the simulation speed slider [default=0.7]\n\
               affects the timer and (for SPEED>0) the skipped frames (stride)\n\
  -T<TIME>     timer delay in s for speed=0 [default=0.05]\n\
  -U<DELAY>    delay in ms added to avoid display jam [default=10]\n\
  <FILE>       (NAME.sim) open sim-file on start (single argument only)\n\
In case of multiple instances of the same option, the last one applies\n\
FLTK options are in lowercase and with a space before a number, get them by:\n\
  simolant -x\n\
EXAMPLE:\n\
  simolant -g 1235x960 -I7 -N600 -Pg=-0.02,tau=0.5\n\
");
      }
    else {
      newarg[newiarg++]=arg[iarg]; } }

  if (N<1) N=1;
  if (N>=MAXN) N=MAXN;

  rndinit(0,0); /* not C++ style :( */
  //  Fl::set_color(FL_BACKGROUND_COLOR,0x05558500);

  Fl_Double_Window win(BOXSIZE+2*BORDER+PANELW,
                       BOXSIZE+2*BORDER+MENUH,"SIMOLANT " VERSION);
  win.color(FL_WHITE);

  if (debug) fprintf(stderr,"geometry:\n  -g %dx%d\n",BOXSIZE+2*BORDER+PANELW,BOXSIZE+2*BORDER+MENUH);

  menu=new Fl_Menu_Bar(0,0,BOXSIZE+2*BORDER+PANELW,MENUH);
  menu->menu(menuitems);
  panel = new Panel(BOXSIZE+2*BORDER, MENUH, PANELW, BOXSIZE+2*BORDER);
  simulate = new Simulate(0, MENUH, BOXSIZE+2*BORDER, BOXSIZE+2*BORDER);
  win.resizable(simulate);

  help = new Fl_Help_Dialog();
  help->textsize(16);
  help->resize(0,0,640,480);

  if (init<GAS || init>=SENTINEL) init=GAS;

  Fl::visual(FL_DOUBLE|FL_INDEX);

  win.show(newiarg,newarg);

  if (!strcmp(getext(fn),".sim")) loadcfg(fn); // name.sim => load it
  else initcfg(init); // no .ext=.sim => is a number
  
  if (optionP) {
    char *c,*e;

    optionP=strdup(optionP);
    for (c=optionP; ; ) {
      e=strchr(c,',');
      read_parm(c);
      if (e) c=e+1;
      else break; }

    free(optionP); }

  setss();
  if (debug) {
    int i;
    printf("# r   u(r)   f(r)   -du(r)/dr_numeric\n");
    loop (i,150,650) {
      double r=i/200.;
      printf("%5.3f %.9g %.9g %.9g\n",
             r,u(r*r),f(r*r)*r,(u(Sqr(r-5e-6))-u(Sqr(r+5e-6)))/1e-5); } }
  
  return Fl::run();
}
