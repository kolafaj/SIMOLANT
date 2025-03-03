// Callbacks of most widgets are here

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

static void cb_LJ(Fl_Widget*, void *data) { setss(LJ,CUTOFF); }
static void cb_WCA(Fl_Widget*, void *data) { setss(WCA,1); }
static void cb_PD(Fl_Widget*, void *data) { ss.a=1; setss(PD,2); }
static void cb_IG(Fl_Widget*, void *data) { ss.a=0; setss(PD,2); }
static void cb_APD(Fl_Widget*, void *data) { ss.a=-1; setss(PD,2); }
static void cb_DW(Fl_Widget*, void *data) { ss.a=1.25; ss.b=1.4; setss(DW,2); }

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

static void cb_VICSEK(Fl_Widget*, void *data) { MDreset(VICSEK); }

// Prepare menu
static void cb_gas(Fl_Widget*, void *data) { initcfg(GAS); }
static void cb_diffusion(Fl_Widget*, void *data) { initcfg(DIFFUSION); }
static void cb_gravity(Fl_Widget*, void *data) { initcfg(GRAVITY); }

static void cb_vle(Fl_Widget*, void *data) { initcfg(VLE); }
static void cb_slab(Fl_Widget*, void *data) { initcfg(SLAB); }
static void cb_nucleation(Fl_Widget*, void *data) { initcfg(NUCLEATION); }

static void cb_liquid(Fl_Widget*, void *data) { initcfg(LIQUIDDROP); }
static void cb_twodrops(Fl_Widget*, void *data) { initcfg(TWODROPS); }
static void cb_cavity(Fl_Widget*, void *data) { initcfg(CAVITY); }
static void cb_capillary(Fl_Widget*, void *data) { initcfg(CAPILLARY); }
static void cb_periodicliquid(Fl_Widget*, void *data) { initcfg(PERIODICLIQUID); }

static void cb_crystal(Fl_Widget*, void *data) { initcfg(CRYSTAL); }
static void cb_defect(Fl_Widget*, void *data) { initcfg(DEFECT); }
static void cb_vacancy(Fl_Widget*, void *data) { initcfg(VACANCY); }
static void cb_intersticial(Fl_Widget*, void *data) { initcfg(INTERSTITIAL); }
static void cb_periodiccrystal(Fl_Widget*, void *data) { initcfg(PERIODICCRYSTAL); }

static void cb_vicsek(Fl_Widget*, void *data) { initcfg(INITVICSEK); }

// Boundary conditions menu
static void cb_box(Fl_Widget*, void *data) { bc=BOX; }
static void cb_periodic(Fl_Widget*, void *data) { bc=PERIODIC; }
static void cb_slit(Fl_Widget*, void *data) { bc=SLIT; }

// Measure menu
static void cb_quantities(Fl_Widget*, void *data) { Show=QUANTITIES; }
static void cb_energy(Fl_Widget*, void *data) { Show=ENERGY; }
static void cb_temperature(Fl_Widget*, void *data) { Show=TEMPERATURE; }
static void cb_intmotion(Fl_Widget*, void *data) { Show=INTMOTION; }
static void cb_pressure(Fl_Widget*, void *data) { Show=PRESSURE; }
static void cb_volume(Fl_Widget*, void *data) { Show=VOLUME; }
static void cb_momentum(Fl_Widget*, void *data) { Show=MOMENTUM; }
static void cb_rdf(Fl_Widget*, void *data) { Show=RDF; }
static void cb_zprofile(Fl_Widget*, void *data) { Show=YPROFILE; }
static void cb_dprofile(Fl_Widget*, void *data) { Show=RPROFILE; }
static void cb_cprofile(Fl_Widget*, void *data) { Show=CPROFILE; }

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

// the same as strcpm() but case insensitive
int nocasecmp(char const *a, char const *b) /********************* nocasecmp */
{
  for (;; a++, b++) {
    int cmp = toupper((unsigned char)*a) - toupper((unsigned char)*b);
    if (cmp || !*a) return cmp; }
}

char printcmd[64],*printcmd0=printcmd;

/* only double version -- not overloaded for int val */
void append(const char *var, double val) /*************************** append */
{
  if (val<=-9e99) printcmd0+=sprintf(printcmd0,"%s ",var);  
  else printcmd0+=sprintf(printcmd0,"%s=%g ",var,val);
  while (fl_width(printcmd)>PANELW/2-21) {
    char *x=strchr(printcmd,' ');
    if (!x) x=printcmd;
    memmove(printcmd,x+1,strlen(x)+1);
    printcmd0=printcmd+strlen(printcmd); }
}

// read in a parameter from a string
static void parse_cmd(char *cmd) /******************************** parse_cmd */
{
  char *c,*cc;
  const char *err=NULL;
  int assign=0;
  double val=0;
  int ival=0;

  /* remove whites (the string shrinks) */
  for (c=cc=cmd; *c; c++) {
    *cc=*c;
    // unsigned char for all strings would be a better solution to treat UTF8
    if ((unsigned)*c>' ') cc++; }
  *cc=0;

  if ( (c=strchr(cmd,'=')) ) {
    /* assignment; accepts garbage after number */
    *c++=0;
    if (strchr("+-.,0123456789",*c)) {
      assign++;
      val=atof(c);
      ival=atoi(c); }
    else
      err="Error:syntax"; }

  if (!err) {
  
    /* variables with sliders: */
    if (!nocasecmp(cmd,"T")) { /* temperature */
      if (assign) {
        BRACKETVAL(val,1e-9,100);
        sliders.T->value(log(val));
        calculateB2(); }
      append("T",exp(sliders.T->value())); }
  
    else if (!nocasecmp(cmd,"tau") || !strcmp(cmd,"τ")) { /* thermostat time */
      if (assign) {
        BRACKETVAL(val,0.01,1e9);
        sliders.tau->value(log(val)); }
      append("τ",exp(sliders.tau->value())); }
  
    else if (!nocasecmp(cmd,"d")) { /* MC step in the units of L */
      if (assign) {
        BRACKETVAL(val,1e-4,1);
        sliders.d->value(log(val)/DSCALE); }
      append("d",exp(DSCALE*sliders.d->value())); }
  
    else if (!nocasecmp(cmd,"g")) { /* gravity */
      if (assign) {
        BRACKETVAL(val,-10,10);
        sliders.g->value(val); }
      append("g",sliders.g->value()); }
  
    else if (!nocasecmp(cmd,"rho") || !strcmp(cmd,"ρ")) { /* density */
      if (assign) {
        BRACKETVAL(val,0.0001,2);
        sliders.rho->value(log(val)); }
      append("ρ",exp(sliders.rho->value())); }
  
    else if (!nocasecmp(cmd,"p")) { /* pressure */
      if (assign) {
        BRACKETVAL(val,-10,100);
        sliders.P->value(val); }
      append("P",sliders.P->value()); }
  
    else if (!nocasecmp(cmd,"l")) { /* box size */
      if (assign) {
        L=val; val=L2rho(L);
        BRACKETVAL(val,0.0001,2); // val=rho
        sliders.rho->value(log(val));
        L=rho2L(val); }
      append("L",L); }
  
    else if (!nocasecmp(cmd,"n")) { /* number of particles N, see insdel() */
      if (assign) {
        BRACKETVAL(ival,1,MAXN);
        sliders.N->value(NtoSlider(ival)); }
      append("N",SlidertoN(sliders.N->value())); }
  
    else if (!nocasecmp(cmd,"stride")) { /* speed, duplicated by speed.stride */
      if (assign) {
        BRACKETVAL(val,-10,100);
        speed.stride=val; // IS THIS NECESSARY?
        slider2speed(val);
        sliders.speed->value(val); }
      append("stride",(double)speed.stride); }
  
    else if (!nocasecmp(cmd,"block")) {
      if (assign) {
        BRACKETVAL(ival,1,1000);
        block=ival; // IS THIS NECESSARY?
        sliders.block->value(log10(block)); }
      append("block",(double)block); }
  
    /* variables without sliders: */
    else if (!nocasecmp(cmd,"dt") || !nocasecmp(cmd,"h")) { /* MD time step */
      if (assign) {
        dtfixed=val;
        if (dtfixed) BRACKETVAL(dtfixed,1e-5,0.1);
        else dtadjset(); }
      append("dt",dtfixed); }
  
    else if (!nocasecmp(cmd,"dV")) { /* MC volume change */
      if (assign) {
        BRACKETVAL(val,1e-4,MAXDV);
        dV=val; }
      append("dV",dV); }
  
    else if (!nocasecmp(cmd,"qtau") || !nocasecmp(cmd,"qτ")) { /* MD qτ=τP/τ */
      if (assign) {
        BRACKETVAL(val,0.5,1e9);
        qtau=val; }
      append("qτ",dV); }
  
    else if (!nocasecmp(cmd,"wall")) { /* wall density */
      if (assign) {
        BRACKETVAL(val,0.1,9);
        walldens=PI*val; }
      append("wall",walldens/PI); }
  
    else if (!nocasecmp(cmd,"ff")) { /* force field */
      if (assign) {
        BRACKETVAL(val,(int)LJ,(int)NFF-1);
        setss((ff_e)(int)val,ss.C2); }
      append("ff",val); }
  
    else if (!nocasecmp(cmd,"c")) { /* cutoff c */
      if (assign) {
        BRACKETVAL(val,MINC,MAXC);
        setss(ss.ff,val); }
      append("c",ss.C2); }
  
    else if (!nocasecmp(cmd,"a")) { /* ff parameter a */
      //      BRACKETVAL(val,?,?)
      if (assign) {
        ss.a=val;
        setss(ss.ff,ss.C2); }
      append("a",ss.a); }
  
    else if (!nocasecmp(cmd,"trace")) { /* lenght of Traces */
      BRACKETVAL(val,4,2048);
      if (assign) trace=val;
      append("trace",trace); }
  
    else if (!nocasecmp(cmd,"b")) { /* ff parameter b */
      //      BRACKETVAL(val,?,?)
      if (assign) {
        ss.b=val;
        setss(ss.ff,ss.C2); }
      append("b",ss.b); }
  
    else if (!nocasecmp(cmd,"nbr")) { /* neighbor distance limit for coloring */
      if (assign) {
        if (val<=0) val=RNBR; // the default
        BRACKETVAL(val,0,4);
        ss.rnbr=val;
        ss.rnbrq=Sqr(val); }
      append("nbr",ss.rnbr); }
  
    else if (!nocasecmp(cmd,"v")) { /* velocity for Two drops */
      if (assign) {
        BRACKETVAL(val,0,10);
        dropvel=val; }
      append("v",dropvel); }
  
    else if (!nocasecmp(cmd,"circle")) { /* circle drawing method */
      if (assign) circlemethod=ival%3;
      append("circle",(double)circlemethod); }
  
    else if (!nocasecmp(cmd,"method")) { /* simulation method, thermostat etc. */
      if (assign) {
        method=(method_e)ival;
        mauto=AUTO;
        if (method<=AUTO || method>=NTH) cb_AUTO(NULL,NULL); }
      append("method",(double)method); }
  
    else if (!nocasecmp(cmd,"bc")) { /* boundary conditions */
      if (assign) {
        bc=(bctype)ival;
        if (bc<BOX || bc>PERIODIC) bc=BOX; }
      append("bc",(double)bc); }
  
    else if (!nocasecmp(cmd,"show")) { /* what to show */
      if (assign) {
        Show=(Show_e)ival;
        if (Show<QUANTITIES || Show>=NSHOW) Show=QUANTITIES; }
      append("Show",(double)Show); }
  
    else
      err="Error:variable"; }

  if (err) append(err,-9e99);
  
  buttons.cmd->value("");

  if (debug) fprintf(stderr,"printcmd[%s]\n",printcmd);

} // parse_cmd()

/* read in parameters from input area [cmd:] */
static void cb_cmd(Fl_Widget* w,void *data) /************************ cb_cmd */
{
  char *cmd=(char*)buttons.cmd->value(),*c;

  // change ,->. so that decimal comma is accepted (not in -P)
  for (c=cmd; *c; c++) if (*c==',') *c='.';

  fprintf(stderr,"%s\n",cmd);

  parse_cmd(cmd);

  buttons.cmd->show();
} // cb_cmd()

static void cb_plus(Fl_Widget* w,void *data) /********************** cb_plus */
{
  int N=SlidertoN(sliders.N->value())+1.5;
  sliders.N->value(NtoSlider(N));
} // cb_plus()

static void cb_gzero(Fl_Widget* w,void *data) /******************** cb_gzero */
{
  sliders.g->value(0);
} // cb_gzero()

static void cb_minus(Fl_Widget* w,void *data) /******************** cb_minus */
{
  int N=SlidertoN(sliders.N->value())-0.5;
  sliders.N->value(NtoSlider(N));
} // cb_minus()

static void cb_resetview(Fl_Widget* w,void *data) /************ cb_resetview */
{
  drawmode->value(0);
  colormode->value(CM_BLACK);
  molsize->value(0);
} // cb_resetview()

static void cb_resetgraph(Fl_Widget* w,void *data) /********** cb_resetgraph */
{
  lastShow=NSHOW;
} // cb_resetgraph

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

static void cb_shift_down(Fl_Widget* w,void *data) /********** cb_shift_down */
{
  int i;
  loop (i,0,N) {
    r[i].y+=L/8;
    if (r[i].y>=L) r[i].y-=L; }
}
  
static void cb_shift_up(Fl_Widget* w,void *data) /************** cb_shift_up */
{
  int i;
  loop (i,0,N) {
    r[i].y-=L/13;
    if (r[i].y<0) r[i].y+=L; }
}

static void cb_shift_left(Fl_Widget* w,void *data) /********** cb_shift_left */
{
  int i;
  loop (i,0,N) {
    r[i].x-=L/13;
    if (r[i].x<0) r[i].x+=L; }
}
  
static void cb_shift_right(Fl_Widget* w,void *data) /******** cb_shift_right */
{
  int i;
  loop (i,0,N) {
    r[i].x+=L/8;
    if (r[i].x>=L) r[i].x-=L; }
}
