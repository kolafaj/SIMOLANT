// File menu (also save measurements)

// replacements for deprecated fl_ask
static int fl_ask1(const char *fmt,const char *arg)
{
  return !fl_choice(fmt,"Yes","No",0,arg);
}

#ifdef NATIVE
static void cb_open(Fl_Widget* w,void *data)
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

  if (fn) loadcfg(fn);
}

static void cb_save(Fl_Widget* w,void *data)
{
  const char *fn=NULL;

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
      char *fnx=(char*)malloc(strlen(fn)+5);
      strcpy(fnx,fn); strcat(fnx,".sim");
      fn=fnx; /* (small) memory leak... */ }

    FILE *out=fl_fopen(fn,"rt");

    if (out) {
      if (!fl_ask1("File \"%s\" exists. Overwrite?",fn)) {
        fclose(out);
        return; }
      fclose(out); }

    writefile(fn); }
}

static void cb_record(Fl_Widget* w,void *data)
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
    const char *ext=getext(fn);
    if (!ext || strcmp(ext,".txt")) {
      // add extension .txt (ALWAYS case-sensitive)
      char *fnx=(char*)malloc(strlen(fn)+5);
      strcpy(fnx,fn); strcat(fnx,".txt");
      fn=fnx; /* (small) memory leak... */ }
    else
      fn=strdup(fn);

    protocol=(char*)fn;
  }
}

#else

static void cb_open(Fl_Widget* w,void *data)
{
  Fl_File_Chooser *fb=new Fl_File_Chooser(".","*.sim",Fl_File_Chooser::SINGLE,"Open");

  fb->show();
  while (fb->shown()) { Fl::wait(); }

  const char *fn=fb->value();

  if (debug) fprintf(stderr,"open file: \"%s\"\n",fn);

  if (fn) loadcfg(fn);
}

static void cb_save(Fl_Widget* w,void *data)
{
  Fl_File_Chooser *fb=new Fl_File_Chooser(".","*.sim",Fl_File_Chooser::CREATE,"Save");

  //  fb->type(FL_NORMAL_BROWSER);
  fb->preview(0);
  fb->show();
  while (fb->shown()) { Fl::wait(); }

  const char *fn=fb->value();
  const char *ext=getext(fn);

  if (debug) fprintf(stderr,"save file: \"%s\"\n",fn);

  if (fn) {
    if (!ext || strcmp(ext,".sim")) {
      // add extension .sim (ALWAYS case-sensitive)
      char *fnx=(char*)malloc(strlen(fn)+5);
      strcpy(fnx,fn); strcat(fnx,".sim");
      fn=fnx; /* (small) memory leak... */ }

    FILE *out=fl_fopen(fn,"rt");

    if (out) {
      if (!fl_ask1("File \"%s\" exists. Overwrite?",fn)) {
        fclose(out);
        return; }
      fclose(out); }

    out=fl_fopen(fn,"wt");
    int i;

    /* see also writefile (NATIVE) */
    if (out) {
      fprintf(out,"%d %d %d %g %g %g %g %g %d %g 0 # N,bc,walls,T,L,walldens,dt,d,thermostat, reserved,reserved\n",
              N,bc,walls,T,L,walldens,dt,d,thermostat,P);
      fprintf(out,"#  x         y\n");
      loop (i,0,N) fprintf(out,"%9.5f %9.5f\n",r[i].x,r[i].y);
      fclose(out); } }
}

static void cb_record(Fl_Widget* w,void *data)
{
  Fl_File_Chooser *fb=new Fl_File_Chooser(".","*.txt",Fl_File_Chooser::CREATE,"Name of the measurement protocol (start/stop [record])");

  //  fb->type(FL_NORMAL_BROWSER);
  fb->preview(0);
  fb->show();
  while (fb->shown()) { Fl::wait(); }

  const char *fn=fb->value();
  const char *ext=getext(fn);

  if (fn) {
    const char *ext=getext(fn);
    if (!ext || strcmp(ext,".txt")) {
      // add extension .txt (ALWAYS case-sensitive)
      char *fnx=(char*)malloc(strlen(fn)+5);
      strcpy(fnx,fn); strcat(fnx,".txt");
      fn=fnx; /* (small) memory leak... */ }
    else
      fn=strdup(fn);

    protocol=(char*)fn;
  }
}
#endif

#ifdef INCL
static void cb_cp(Fl_Widget* w,void *data)
{
  writeCP=fl_ask1("Include full convergence profile in the recorded protocol?",protocol);
}

static void cb_df(Fl_Widget* w,void *data)
{
  writeDF=fl_ask1("Include RDF, y- or r-density functions in the recorded protocol?",protocol);
}
#endif

static void cb_exit(Fl_Widget* w,void *data)
{
  //KM Close all windows - better than exit exit()
  while (Fl::first_window())
    Fl::first_window()->hide();
}

static void cb_resetchoices(Fl_Widget* w,void *data)
{
  drawmode->value(0);
  colormode->value(0);
  molsize->value(0);
  justreset=1;
  itot=HISTMAX+1;
}

// Method menu
static void cb_AUTO(Fl_Widget*, void *data) { thermostat=METROPOLIS; mauto=BUSSI; mcstart=60; }

static void cb_MC_C(Fl_Widget*, void *data) { thermostat=CREUTZ; }
static void cb_MC(Fl_Widget*, void *data) { thermostat=METROPOLIS; }
static void cb_MC_NPT(Fl_Widget*, void *data) { thermostat=MCNPT; dV=0.02; }

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
static void cb_capillary(Fl_Widget*, void *data) { initcfg(CAPILLARY); }
static void cb_cavity(Fl_Widget*, void *data) { initcfg(CAVITY); }

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
static void cb_help(Fl_Widget* w,void *data)
{
  FILE *f=fl_fopen("simolant.html","rt");

  if (f) {
    fclose(f);
    help->load("simolant.html");
    help->show(); }
  else
    fl_alert("\
Missing manual \"simolant.html\".\n\
Downolad it from\n\
http://www.vscht.cz/fch/software/simolant/index-en.html\n\
and copy to the directory of the simolant executable.");
}

static void cb_about(Fl_Widget* w,void *data)
{
  fl_message("SIMOLANT " VERSION ", © Jiří Kolafa\n\
License:\n\
    GNU General Public License 3\n\
Main page:\n\
    http://www.vscht.cz/fch/software/simolant/index-en.html\n\
Supported by:\n\
    1985–2001 Institute of Chemical Process Fundamentals\n\
    2001–2020 University of Chemistry and Technology, Prague\n\
    2012–2013 TUL Liberec, EU: esf (evropský sociální fond\n\
        v ČR), MŠMT: OP Vzdělání pro konkurenceschopnost\n\
Acknowledgements:\n\
    Ivo Nezbeda (pioneering simulation workshop in ~1985)\n\
    Karel Matas (FLTK advise)");
}

#define bracketval(VAL,FROM,TO) do { \
  if (VAL<FROM) VAL=FROM; \
  if (VAL>TO) VAL=TO; } while(0)

/* read in a parameter from a string */
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
      bracketval(val,1e-9,99);
      sliders.T->value(log(val)); }
    else if (!strcmp(cmd,"tau") || !strcmp(cmd,"τ")) {
      bracketval(val,0.01,1e9);
      sliders.tau->value(log(val)); }
    else if (!strcmp(cmd,"d")) {
      /* in the units of L */
      bracketval(val,1e-3,1);
      sliders.d->value(log(val)/DSCALE); }
    else if (!strcmp(cmd,"dv")) {
      /* MC volume change */
      bracketval(val,1e-3,0.1);
      dV=val; }
    else if (!strcmp(cmd,"g")) {
      bracketval(val,-9,9);
      sliders.g->value(val);}
    else if (!strcmp(cmd,"rho") || !strcmp(cmd,"ρ")) {
      bracketval(val,0.0001,2);
      sliders.rho->value(log(val)); }
    else if (!strcmp(cmd,"p")) {
      bracketval(val,-9,99);
      sliders.P->value(val); }
    else if (!strcmp(cmd,"n")) {
      bracketval(ival,1,MAXN);
      sliders.N->value(pow(ival,1./NPOW)); }
    else if (!strcmp(cmd,"dt") || !strcmp(cmd,"h")) {
      dtfixed=val;
      if (val) bracketval(val,1e-9,1);
      else dtadjset(); }
    else if (!strcmp(cmd,"stride")) {
      nit=ival;
      bracketval(ival,1,99);
      timerdelay=TIMERDELAY*exp(-1.609438*nit);
      sliders.speed->value(log(nit)/1.1+1e-5); }
    else if (!strcmp(cmd,"block")) {
      bracketval(ival,1,1000);
      block=ival;
      sliders.block->value(log10(block)); }
    /* no sliders from here: */
    else if (!strcmp(cmd,"qtau") || !strcmp(cmd,"qτ")) {
      bracketval(val,0.5,1e9);
      qtau=val; }
    else if (!strcmp(cmd,"wall")) {
      bracketval(val,0.1,9);
      walldens=PI*val; }
    else if (!strcmp(cmd,"cutoff") || !strcmp(cmd,"cut")) {
      bracketval(val,2.4,1000);
      cutoff=val;
      setss(); }
    else if (!strcmp(cmd,"l")) {
      double L0=L;
      L=val; val=L2rho(L);
      bracketval(val,0.0001,2); // this is rho
      L=L0;
      sliders.rho->value(log(val)); }
    else if (!strcmp(cmd,"th")) {
      thermostat=(thermostat_t)ival;
      mauto=AUTO;
      if (thermostat<=AUTO || thermostat>=NTH) cb_AUTO(NULL,NULL); }
    else if (!strcmp(cmd,"bc")) {
      bc=(bctype)ival;
      if (bc<BOX || bc>PERIODIC) bc=BOX; }
    else
      buttons.parms->value("? variable"); }
  else
    buttons.parms->value("? number");
}

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
}
