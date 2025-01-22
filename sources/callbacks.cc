// Callbacks of most widgets are here

// replacements for deprecated fl_ask
static int fl_ask1(const char *fmt,const char *arg) /*************** fl_ask1 */
{
  return !fl_choice(fmt,"Yes","No",0,arg);
} // fl_ask1()

// File menu
// choose and load a sim-file (a configuration)
static void cb_load(Fl_Widget* w,void *data) /********************** cb_load */
{
  const char *fn=NULL;
  Fl_Native_File_Chooser fnfc;

  fnfc.title("Open");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  fnfc.filter("SIMOLANT\t*.sim");
  fnfc.directory(".");
  switch (fnfc.show()) {
    case -1: fl_alert("%s\n", fnfc.errmsg()); break;
    case  1: break;  // CANCEL
    default: fn=fnfc.filename(); }

  if (debug) fprintf(stderr,"open file: \"%s\"\n",fn);

  if (fn) loadsim(fn);
} // cb_load()

// choose/type filename and save a sim-file (a configuration)
static void cb_save(Fl_Widget* w,void *data) /********************** cb_save */
{
  const char *fn=NULL;
  char *fnext=NULL;
  Fl_Native_File_Chooser fnfc;

  fnfc.title("Save as");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
  fnfc.filter("SIMOLANT\t*.sim");
  fnfc.directory(".");
  switch (fnfc.show()) {
    case -1: fl_alert("%s\n", fnfc.errmsg()); break;
    case  1: break;  // CANCEL
    default: fn=fnfc.filename(); }

  if (debug) fprintf(stderr,"save file: \"%s\"\n",fn);

  if (fn) {
    const char *ext=getext(fn);
    if (!ext || strcmp(ext,".sim")) {
      // add extension .sim (ALWAYS case-sensitive)
      fnext=(char*)malloc(strlen(fn)+5);
      sprintf(fnext,"%s.sim",fn);
      fn=fnext; }

    FILE *out=fl_fopen(fn,"rt");

    if (out) {
      if (!fl_ask1("File \"%s\" exists. Overwrite?",fn)) {
        fclose(out);
        goto ret_cb_save; }
      fclose(out); }

    savesim(fn);
  ret_cb_save:
    if (fnext) free(fnext);
  }
} // cb_save()

// type protocol name or choose in *.txt files (extension .txt will be removed)
static void cb_protocol(Fl_Widget* w,void *data) /************** cb_protocol */
{
  const char *fn=NULL;
  Fl_Native_File_Chooser fnfc;

  fnfc.title("Name of the measurement protocol (start/stop [record])");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
  fnfc.filter("TEXT\t*.txt");
  fnfc.directory(".");
  switch (fnfc.show()) {
    case -1: fl_alert("%s\n", fnfc.errmsg()); break;
    case  1: break;  // CANCEL
    default: fn=fnfc.filename(); }

  if (debug) fprintf(stderr,"record file: \"%s\"\n",fn);

  if (fn) {
    int sl=strlen(fn);

    if (sl>243) {
      //     456789012345
      // NAME-RDF####.csv accepted
      fl_alert("filename too long");
      return; }

    strcpy(files.protocol,fn);
    // const removed because files.protocol is not const
    files.ext=(char*)getext(files.protocol);
    if (files.ext && !strcmp(files.ext,".txt"))
      *files.ext=0;
    else
      files.ext=files.protocol+sl;
    files.nmeas=0;
    fl_alert("\
The protocol name has changed!\n\
Measurements will be numbered from 1.\n\
If the following files:\n\
 %s.txt\n\
 %s.csv\n\
exist, they will be overwitten on recording.",files.protocol,files.protocol);
  }
} // cb_protocol()

static void cb_exit(Fl_Widget* w,void *data) /********************** cb_exit */
{
  //KM Close all windows - better than exit exit()
  while (Fl::first_window())
    Fl::first_window()->hide();
} // cb_exit()

FILE *ffout;
void comma(char *s); // defined in record.cc
char ffsep; // files.sep for csv, tabulator if txt output

void ffprtcsv(double d) /****************************************** ffprtcsv */
// similar to prtcsv in record()
{
  char x[24];

  sprintf(x,"%.10g",d);
  comma(x);
  fprintf(ffout,"%c%s",ffsep,x);
} // ffprtcsv()

static void cb_ff(Fl_Widget* w,void *data) /************************** cb_ff */
{
  // print force field (potential+forces) file
  int i;

  if (buttons.csv->value()) {
    strcpy(files.ext,"-ff.csv");
    ffsep=files.sep; }
  else {
    strcpy(files.ext,"-ff.txt");
    ffsep='\t'; }

  ffout=fl_fopen(files.protocol,"wt");
  fprintf(ffout,"#r%cufull(r)%cu(r)%cffull(r)%cf(r)%c-du(r)/dr_numeric\n",
          ffsep,ffsep,ffsep,ffsep,ffsep);
  loop (i,0,(int)(fabs(ss.C2)*100+10.5)) {
    double r=i/100.;
    if (u(r*r)<1e5) {
      char x[16];
      sprintf(x,"%4.2f",r);
      comma(x);
      fprintf(ffout,"%s",x);
      ffprtcsv(ufull(r*r));
      ffprtcsv(u(r*r));
      ffprtcsv(ffull(r*r)*r);
      ffprtcsv(f(r*r)*r);
      ffprtcsv((u(Sqr(r-5e-6))-u(Sqr(r+5e-6)))/1e-5);
      fprintf(ffout,"\n"); } }
  fprintf(ffout,"#params:%ca%cb%cc\n#",ffsep,ffsep,ffsep);
  ffprtcsv(sqrt(ss.aq));
  ffprtcsv(sqrt(ss.bq));
  ffprtcsv(ss.C2);
  fprintf(ffout,"\n");
  fclose(ffout);
} // cb_ff()

// Force field menu

static void cb_LJ4(Fl_Widget*, void *data) { setss(CUTOFF); }
static void cb_WCALJ(Fl_Widget*, void *data) { setss(1); }
static void cb_PD(Fl_Widget*, void *data) { ss.a=0.5; setss(2); }
static void cb_IG(Fl_Widget*, void *data) { ss.a=0; setss(2); }
static void cb_APD(Fl_Widget*, void *data) { ss.a=-0.5; setss(2); }
static void cb_DW(Fl_Widget*, void *data) { ss.a=1.25; ss.b=1.4; setss(-2); }

// Method menu
static void cb_AUTO(Fl_Widget*, void *data) { MCreset(METROPOLIS); mauto=BUSSI; mcstart=60; }

static void cb_MC_C(Fl_Widget*, void *data) { MCreset(CREUTZ); }
static void cb_MC(Fl_Widget*, void *data) { MCreset(METROPOLIS); }
static void cb_MC_NPT(Fl_Widget*, void *data) { MCreset(MCNPT); dV=0.02; }

static void cb_MD(Fl_Widget*, void *data) { MDreset(NVE); }
static void cb_MD_B(Fl_Widget*, void *data) { MDreset(BERENDSEN); }
static void cb_MD_NH(Fl_Widget*, void *data) { MDreset(NOSE_HOOVER); }
static void cb_MD_A(Fl_Widget*, void *data) { MDreset(ANDERSEN); }
static void cb_MD_M(Fl_Widget*, void *data) { MDreset(MAXWELL); }
static void cb_MD_L(Fl_Widget*, void *data) { MDreset(LANGEVIN); }
static void cb_MD_BUSSI(Fl_Widget*, void *data) { MDreset(BUSSI); }
static void cb_MD_MDNPT(Fl_Widget*, void *data) { MDreset(MDNPT); }
static void cb_MD_MTK(Fl_Widget*, void *data) { MDreset(MTK); }

// Prepare menu
static void cb_gas(Fl_Widget*, void *data) { initcfg(GAS); }
static void cb_diffusion(Fl_Widget*, void *data) { initcfg(DIFFUSION); }
static void cb_gravity(Fl_Widget*, void *data) { initcfg(GRAVITY); }

static void cb_vle(Fl_Widget*, void *data) { initcfg(VLE); }
static void cb_slab(Fl_Widget*, void *data) { initcfg(SLAB); }
static void cb_nucleation(Fl_Widget*, void *data) { initcfg(NUCLEATION); }

static void cb_liquid(Fl_Widget*, void *data) { initcfg(LIQUID); }
static void cb_twodrops(Fl_Widget*, void *data) { initcfg(TWODROPS); }
static void cb_cavity(Fl_Widget*, void *data) { initcfg(CAVITY); }
static void cb_capillary(Fl_Widget*, void *data) { initcfg(CAPILLARY); }

static void cb_crystal(Fl_Widget*, void *data) { initcfg(CRYSTAL); }
static void cb_defect(Fl_Widget*, void *data) { initcfg(DEFECT); }
static void cb_vacancy(Fl_Widget*, void *data) { initcfg(VACANCY); }
static void cb_intersticial(Fl_Widget*, void *data) { initcfg(INTERSTITIAL); }

// Boundary conditions menu
static void cb_box(Fl_Widget*, void *data) { bc=BOX; }
static void cb_periodic(Fl_Widget*, void *data) { bc=PERIODIC; }
static void cb_slit(Fl_Widget*, void *data) { bc=SLIT; }

// Measure menu
static void cb_nomeas(Fl_Widget*, void *data) { measure=NONE; }
static void cb_quantities(Fl_Widget*, void *data) { measure=QUANTITIES; }
static void cb_energy(Fl_Widget*, void *data) { measure=ENERGY; }
static void cb_temperature(Fl_Widget*, void *data) { measure=TEMPERATURE; }
static void cb_intmotion(Fl_Widget*, void *data) { measure=INTMOTION; }
static void cb_pressure(Fl_Widget*, void *data) { measure=PRESSURE; }
static void cb_volume(Fl_Widget*, void *data) { measure=VOLUME; }
static void cb_rdf(Fl_Widget*, void *data) { measure=RDF; }
static void cb_zprofile(Fl_Widget*, void *data) { measure=YPROFILE; }
static void cb_dprofile(Fl_Widget*, void *data) { measure=RPROFILE; }
static void cb_cprofile(Fl_Widget*, void *data) { measure=CPROFILE; }

// Help menu
static void cb_help(Fl_Widget* w,void *data) /********************** cb_help */
{
  FILE *f=fl_fopen("simolant.html","rt");

  if (f) {
    fclose(f);
    help->load("simolant.html");
    help->show(); }
  else
    fl_alert("\
Missing manual \"simolant.html\".\n\
Run SIMOLANT from the directory where\n\
it has been installed, or download it from\n\
https://github.com/kolafaj/SIMOLANT.");
} // cb_help()

static void cb_about(Fl_Widget* w,void *data) /******************** cb_about */
{
  fl_message("SIMOLANT " VERSION ", © Jiří Kolafa\n\
License:\n\
    GNU General Public License 3\n\
Where is:\n\
    https://github.com/kolafaj/SIMOLANT\n\
    http://old.vscht.cz/fch/software/simolant/index-en.html\n\
Supported by:\n\
    1985–2001 Institute of Chemical Process Fundamentals\n\
    2001–2020 University of Chemistry and Technology, Prague\n\
    2012–2013 TUL Liberec, EU: esf (evropský sociální fond\n\
        v ČR), MŠMT: OP Vzdělání pro konkurenceschopnost\n\
Acknowledgements:\n\
    Ivo Nezbeda (pioneering simulation workshop in ~1985)\n\
    Karel Matas (FLTK advise)\n\
    students");
} // cb_about()

// Tinker menu
static void cb_fps15(Fl_Widget*, void *data) { speed.FPS=15; }
static void cb_fps30(Fl_Widget*, void *data) { speed.FPS=30; }
static void cb_fps60(Fl_Widget*, void *data) { speed.FPS=60; }
static void cb_fps120(Fl_Widget*, void *data) { speed.FPS=120; }
static void cb_timeout99(Fl_Widget*, void *data) { speed.MINTIMEOUT=0.99; }
static void cb_timeout50(Fl_Widget*, void *data) { speed.MINTIMEOUT=0.5; }
static void cb_timeout20(Fl_Widget*, void *data) { speed.MINTIMEOUT=0.2; }
static void cb_timeout7(Fl_Widget*, void *data) { speed.MINTIMEOUT=0.07; }

// parameter value normalized to range [FROM,TO]
#define BRACKETVAL(VAL,FROM,TO) do { \
  if (VAL<FROM) VAL=FROM; \
  if (VAL>TO) VAL=TO; } while(0)

// read in a parameter from a string
static void read_parm(char *cmd) /******************************** read_parm */
{
  char *nr;

  /* case insensitive */
  for (nr=cmd; *nr; nr++) *nr=tolower(*nr);

  for (nr=cmd; *nr!='='; nr++)
    if (nr>cmd+6) {
      /* max variable length */
      buttons.parms->value("? variable");
      return; }
  *nr++=0;

  if (strchr("0123456789-+.",nr[0])) {
    double val=atof(nr);
    int ival=atoi(nr);
    buttons.parms->value("");

    if (!strcmp(cmd,"t")) {
      BRACKETVAL(val,1e-9,99);
      sliders.T->value(log(val));
      calculateB2(); }
    else if (!strcmp(cmd,"tau") || !strcmp(cmd,"τ")) {
      BRACKETVAL(val,0.01,1e9);
      sliders.tau->value(log(val)); }
    else if (!strcmp(cmd,"d")) {
      /* in the units of L */
      BRACKETVAL(val,1e-4,1);
      sliders.d->value(log(val)/DSCALE); }
    else if (!strcmp(cmd,"dv")) {
      /* MC volume change */
      BRACKETVAL(val,1e-4,MAXDV);
      dV=val; }
    else if (!strcmp(cmd,"g")) {
      BRACKETVAL(val,-9,9);
      sliders.g->value(val);}
    else if (!strcmp(cmd,"rho") || !strcmp(cmd,"ρ")) {
      BRACKETVAL(val,0.0001,2);
      sliders.rho->value(log(val)); }
    else if (!strcmp(cmd,"p")) {
      BRACKETVAL(val,-9,99);
      sliders.P->value(val); }
    else if (!strcmp(cmd,"n")) {
      BRACKETVAL(ival,1,MAXN);
      sliders.N->value(pow(ival,1./NPOW)); }
    else if (!strcmp(cmd,"dt") || !strcmp(cmd,"h")) {
      dtfixed=val;
      if (val) BRACKETVAL(val,1e-9,1);
      else dtadjset(); }
    else if (!strcmp(cmd,"stride")) {
      speed.stride=val;
      BRACKETVAL(val,-10,100);
      slider2speed(val);
      sliders.speed->value(val); }
    else if (!strcmp(cmd,"block")) {
      BRACKETVAL(ival,1,1000);
      block=ival;
      sliders.block->value(log10(block)); }
    /* no sliders from here: */
    else if (!strcmp(cmd,"qtau") || !strcmp(cmd,"qτ")) {
      BRACKETVAL(val,0.5,1e9);
      qtau=val; }
    else if (!strcmp(cmd,"wall")) {
      BRACKETVAL(val,0.1,9);
      walldens=PI*val; }
    else if (!strcmp(cmd,"c")) {
      //      BRACKETVAL(val,2.4,1000);
      // CHANGE BRACKETING?
      setss(val); }
    else if (!strcmp(cmd,"a")) {
      //      BRACKETVAL(val,?,?)
      ss.a=val; 
      setss(ss.C2); }
    else if (!strcmp(cmd,"b")) {
      //      BRACKETVAL(val,?,?)
      ss.b=val; 
      setss(ss.C2); }
    else if (!strcmp(cmd,"nbr")) {
      if (val<=0) val=RNBR; // the default
      BRACKETVAL(val,0,4);
      ss.rnbr=val;
      ss.rnbrq=Sqr(val); }
    else if (!strcmp(cmd,"l")) {
      double L0=L;
      L=val; val=L2rho(L);
      BRACKETVAL(val,0.0001,2); // this is rho
      L=L0;
      sliders.rho->value(log(val)); }
    else if (!strcmp(cmd,"circle"))
      circlemethod=ival;
    else if (!strcmp(cmd,"method")) {
      method=(method_e)ival;
      mauto=AUTO;
      if (method<=AUTO || method>=NTH) cb_AUTO(NULL,NULL); }
    else if (!strcmp(cmd,"bc")) {
      bc=(bctype)ival;
      if (bc<BOX || bc>PERIODIC) bc=BOX; }
    else if (!strcmp(cmd,"measure") || !strcmp(cmd,"show")) {
      measure=(measure_t)ival;
      if (measure<NONE || measure>=NMEASURE) measure=NONE; }
    else
      buttons.parms->value("ERROR variable"); }
  else
    buttons.parms->value("ERROR number");
} // read_parm()

/* read in parameters from input area [cmd:] */
static void cb_parms(Fl_Widget* w,void *data) /******************** cb_parms */
{
  char *cmd=(char*)buttons.parms->value(),*c,*cc;

  // change , into . so that decimal comma is accepted and remove whites
  for (c=cc=cmd; *c; c++)
    if (*c>' ') {
      if (*c==',') *cc='.';
      else *cc=*c;
      cc++; }
  *cc=0;

  fprintf(stderr,"%s\n",cmd);

  read_parm(cmd);

  buttons.parms->show();
} // cb_parms()

static void cb_resetview(Fl_Widget* w,void *data) /************ cb_resetview */
{
  drawmode->value(0);
  colormode->value(CM_BLACK);
  molsize->value(0);
  justreset=1;
  itot=HISTMAX+1;
} // cb_resetview()

static void cb_invertwalls(Fl_Widget* w,void *data) /******** cb_invertwalls */
{
  buttons.wallx->value(!buttons.wallx->value());
  buttons.wallxL->value(!buttons.wallxL->value());
  buttons.wally->value(!buttons.wally->value());
  buttons.wallyL->value(!buttons.wallyL->value());

  walls=(buttons.wallx->value())
    + 2*(buttons.wallxL->value())
    + 4*(buttons.wally->value())
    + 8*(buttons.wallyL->value());

  setwalls();
} // cb_invertwalls
