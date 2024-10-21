// Data recording - activated by radio buttom [▯ record]

struct res_s {
  struct res_s *next;
  const char *name;
  int n;
  double sum,sumq;
} *head;

// to store the convergence profiles CP)
struct CP_s {
  struct CP_s *next;
  char line[1]; // variable length
} *CPhead=NULL;

// change spaces to underscores, this helps simple plots to use CSV
static char *space2_(const char *s) /******************************* space2_ */
{
  static char conv[16];
  char *c;

  if (strlen(s)>15) return (char *)s;
  strcpy(conv,s);
  for (c=conv; *c; c++) if (*c==' ') *c='_';
  
  return conv;
}

/* add measurement (for averages; if ((incl->value()&1)), then also append to CP */
void addM(const char *name,double x) /********************************* addM */
{
  struct res_s *a;
  static char line[1024]; // WARNING, fixed size!
  char item[32];
  int litem;

  if ((incl->value()&1)) {
    if (!name) {
      struct CP_s *CP;
      static struct CP_s *CPlast;

      litem=sprintf(item,"%g\t",t);
      memmove(line+litem,line,strlen(line)+1);
      memcpy(line,item,litem);

      CP=(struct CP_s *)malloc(sizeof(CP_s)+strlen(line));
      CP->next=NULL;
      strcpy(CP->line,line);

      if (CPhead) {
        CPlast->next=CP;
        CPlast=CP; }
      else
        CPhead=CPlast=CP;

      line[0]=0;

      return; }

    litem=sprintf(item,"%g\t",x);
    memmove(line+litem,line,strlen(line)+1);
    memcpy(line,item,litem); }

  if (!name) return;

  looplist (a,head) if (!strcmp(name,a->name)) {
    a->sum+=x;
    a->sumq+=x*x;
    a->n++;
    return; }

  a=(struct res_s*)malloc(sizeof(struct res_s));
  // a=new struct res_s;
  a->next=head;
  a->name=name;
  head=a;
  head->n=1;
  a->sum=x;
  a->sumq=x*x;
}

// change decimal points to commas if requested
void comma(char *s) /************************************************* comma */
{
  if (buttons.comma->value())
    while (*s) {
      if (*s=='.') *s=',';
      s++; }
}

// change tabulators to the csv separator (, if not comma, ; if comma)
void csvline(char *s) /********************************************* csvline */
{
  for (char *c=s; *c; c++) if (*c=='\t') *c=files.sep;
}

// print to csv file, overloaded, separator first
void prtcsv(int i) /************************************************* prtcsv */
{
  fprintf(files.csv,"%c%d",files.sep,i);
}

void prtcsv(double d) /********************************************** prtcsv */
{
  char x[16];
  
  sprintf(x,"%g",d);
  comma(x);
  fprintf(files.csv,"%c%s",files.sep,x);
}

void prtcsv(const char *s) /***************************************** prtcsv */
{
  fprintf(files.csv,"%c%s",files.sep,s); // on purpose without ""
}

void info2csv(const char *key) /*********************************** info2csv */
{
  char *c=strstr(files.info,key);

  if (c) {
    c=strchr(c,'=');
    if (!c) goto bad;
    c+=2;
    fputc(files.sep,files.csv);
    while ((unsigned)*c>' ') { fputc(*c,files.csv); c++; }
    c++;
    while ((unsigned)*c>' ') c++;
    c++;
    fputc(files.sep,files.csv);
    while ((unsigned)*c>' ') { fputc(*c,files.csv); c++; }
    return; }
 bad:
  prtcsv("n.a.");
  prtcsv("n.a.");
}

const char *wallsg(int key) /**************************************** wallsg */
{
  return key?"attractive":"repulsive";
}

// print one item to csv

// show the measurements and record them
// returns the record status (1=continue, 0=erased)
int showclearM() /*********************************************** showclearM */
{
  struct res_s *a,*next;
  char *end;
  char butt0[80];
  static int overwrite=1,i;
  int choice;
  double shift=0.5;
  FILE *out;
  struct CP_s *CP;

  if (buttons.csv->value()) {

  }
  else
    files.csv=NULL;

  files.sep=",;"[(int)buttons.comma->value()]; // global

  if (!head) return 0; // recording already off

  // recording just has turned off

  // generate info string to pop up and write to a file
  // also, csv generated from this (cumbersome...)
  end=files.info;

  looplist (a,head) {
    double av=a->sum/a->n;
    double stder=(a->sumq/a->n-Sqr(av))/(a->n-1);

    if (a->n<2 || stder<=0) stder=0;
    else stder=sqrt(stder);

    if (a==head) end+=sprintf(end,"\
%2d blocks (block length=%d, stride=%d)\n\
----------------------------------------\n\
quantity = average ± standard error (relative error)\n\
",a->n,block,speed.stride);

    end+=sprintf(end,"%s = %g ± %.3g (%.3g %%)\n",
                 a->name,
                 av,stder,fabs(stder/av*100)); }
  comma(files.info);

  if (head->n<2)
    choice=!fl_choice("Not enough blocks have been recorded.","Stop recording","Continue",NULL);
  else {
    if (overwrite) {
      if (strlen(files.protocol)<26) sprintf(butt0,"save (overwrite \"%s.{txt,csv}\") and clear",files.protocol);
      else strcpy(butt0,"save (overwrite protocol) and clear"); }
    else {
      if (strlen(files.protocol)<26) sprintf(butt0,"append to \"%s.{txt,csv}\" and clear",files.protocol);
      else strcpy(butt0,"append to protocol and clear"); }

    choice=fl_choice("%s","continue",butt0,"clear (data lost)",files.info);
    // this is warning:
    //   choice=fl_choice(info,"continue",butt0,"clear (data lost)");
  }

  switch (choice) {
    case 0: // [continue]
      return 1;

    case 1: // [save]
      strcpy(files.ext,".txt");
      out=fl_fopen(files.protocol,overwrite?"wt":"at");
      if (buttons.csv->value()) {
        strcpy(files.ext,".csv");
        files.csv=fl_fopen(files.protocol,overwrite?"wt":"at");
        if (!files.csv) fl_alert("cannot write the summary files.csv file"); }
      else
        files.csv=NULL;

      if (!out)
        fl_alert("cannot write the protocol");
      else {
        char lout[512]; // 2* safety size

        files.nmeas++;
        overwrite=0;

        if (files.nmeas==1)
          fprintf(out,"SIMOLANT VERSION %s\n",VERSION);

        sprintf(lout,"\
\n\
=========== MEASUREMENT =========== # %d ===========\n\
Parameters:\n\
  N=%d number of particles\n\
  bc=%s boundary conditions\n\
  method=\"%s\" MC/MD, thermostat/barostat\n\
  T=%g temperature\n\
  P=%g pressure\n\
  g=%g gravity\n",
                files.nmeas,
                N,
                bc==BOX?"Box":bc==SLIT?"Slit":"Periodic",
                thinfo[thermostat],T,gravity,P);
        comma(lout);
        fputs(lout,out);

        if (files.csv) {
          if (files.nmeas==1)
            fprintf(files.csv,"#,N,bc,method,T,g,P,walls,rho_wall,L,rho,cutoff,tau,qtau,dt|d,dV,speed.stride,block,Etot,err,Tkin,err,Epot,err,V,err,Z,err,Pvir,err,Pxx,err,Pyy,err,γ,err,P(right_wall),err,P(left_wall),err,P(top_wall),err,P(bottom_wall),err,enthalpy,err,Econserved,err\n");
          fprintf(files.csv,"%d",files.nmeas);
          prtcsv(N);
          prtcsv(bc);
          fprintf(files.csv,"%c\"%s\"",files.sep,space2_(thinfo[thermostat])); // do not use prtcsv(char*)
          prtcsv(T);
          prtcsv(gravity);
          prtcsv(P);
          prtcsv(walls);
          prtcsv(walldens/PI);
          prtcsv(L);
          prtcsv(L2rho(L));
          prtcsv(cutoff);
          prtcsv(tau);
          prtcsv(qtau);
          if (isMC(thermostat)) prtcsv(dt);
          else prtcsv(d);
          prtcsv(dV);
          prtcsv(speed.stride);
          prtcsv(block); }

        switch (bc) {
          case BOX:
            sprintf(lout,"\
  walls:\n\
    left=%s right=%s bottom=%s top=%s\n\
    rho_wall=%g (number density of smoothed wall atoms)\n",
                    wallsg(walls&1),wallsg(walls&2),wallsg(walls&4),wallsg(walls&8),
                    walldens/PI);
            break;
          case SLIT:
            sprintf(lout,"\
  walls:\n\
    bottom=%s top=%s\n\
    rho_wall=%g (number density of smoothed wall atoms)\n",
                    wallsg(walls&4),wallsg(walls&8),
                    walldens/PI);            break;
          default:
            sprintf(lout,"  walls: not present (rho_wall=%g)\n",
                    walldens/PI); }
        comma(lout);
        fputs(lout,out);

        sprintf(lout,"\
  L=%g box side\n\
  rho=%g  number density N/L²\n\
  cutoff=%g potential cutoff (is also smoothed)\n",L,L2rho(L),cutoff);
        comma(lout); fputs(lout,out);

        sprintf(lout,"\
  tau=%g thermostat characteristic time (MD)\n\
  qtau=%g =Ptau/tau, where Ptau = barostat characteristic time (NPT MD)\n\
  dt=%g leap-frog integration step (MD)\n",tau,qtau,dt);
        comma(lout);
        fputs(lout,out);

        sprintf(lout,"\
  d=%g(%s) trial displacement in L (MC)\n\
  dV=%g trial volume change (NPT MC)\n",d,setd->value()?"auto":"fixed",dV);
        comma(lout);
        fputs(lout,out);

        sprintf(lout,"\
  stride=%d quantities are measured every %d%s MD step or MC sweep\n\
  block=%d quantities are averaged in blocks by %d measured points\n",
                speed.stride,speed.stride,
                speed.stride==1?"st":speed.stride==2?"nd":speed.stride==3?"rd":"th",
                block,block);
        //        comma(lout);
        fputs(lout,out);

        fprintf(out,"\n----------------------------------------\n");
        fputs(files.info,out);
        //        fprintf(stderr,"strlen(info)=%d\n",(int)(strlen(info)));
        // <600 bytes detected, char info[2048] = safe size

        if (files.csv) {
          // writing csv from re-parsed info string
          // the order in info may differ
          // should match the CSV header
          info2csv("Etot");
          info2csv("Tkin");
          info2csv("Epot");
          info2csv("V");
          info2csv("Z");
          info2csv("Pvir");
          info2csv("Pxx");
          info2csv("Pyy");
          info2csv("γ");
          info2csv("P(right_wall)");
          info2csv("P(left_wall)");
          info2csv("P(top_wall)");
          info2csv("P(bottom_wall)");
          info2csv("H");
          info2csv("Econserved");
          fprintf(files.csv,"\n"); }

        if (CPhead && (incl->value()&1)) { // convergence profiles

          FILE *cpcsv=NULL;

          if (buttons.csv->value()) {
            sprintf(files.ext,"-CP%d.csv",files.nmeas);
            cpcsv=fl_fopen(files.protocol,"wt");
            if (!cpcsv) fl_alert("cannot write to CSV file (will use TXT)"); }

          if (cpcsv) {
            fprintf(cpcsv,"#t");
            looplist (a,head) {
              fprintf(cpcsv,"%c",files.sep);
              fprintf(cpcsv,"%s",a->name); }
            fprintf(cpcsv,"\n"); }
          else {
            fprintf(out,"\n\
----------------------------------------------------\n\
Convergence profile (time development of quantities)\n\
----------------------------------------------------\n");

            fprintf(out,"time");
            i='A';
            looplist (a,head) {
              fprintf(out,"\t%c=%s",i,a->name);
              i++; }
            fprintf(out,"\tKEY\n"); }

          for (CP=CPhead; CP; ) {
            struct CP_s *next=CP->next;

            comma(CP->line);

            if (cpcsv) {
              csvline(CP->line);
              fprintf(cpcsv,"%s\n",CP->line); }
            else
              fprintf(out,"%sCP%d\n",CP->line,files.nmeas);

            free(CP);
            CP=next; }

          CPhead=NULL;
          if (cpcsv) fclose(cpcsv);
          else fputs("\n",out); }

        if (measure>=RDF && (incl->value()&2)) { // density profiles
          int i;
          double q,cumul=0,Q;
          const char *dpkey=NULL;

          FILE *dpcsv=NULL;

          if (!files.csv) fputs("------------------------------\n",out);

          switch (measure) {
            case RDF:
              dpkey="RDF";
              q=1./50;
              if (isNPT(thermostat)) Q=L2rho(sum.V/iblock);
              else Q=L2rho(L);
              Q*=PI/5000; /* 50^2*2 */
              if (!files.csv) fprintf(out,"\
Radial distribution function\n\
------------------------------\n\
r\tg(r)\tN(r)\tKEY\n");
              break;
            case YPROFILE:
              dpkey="VDP";
              if (isNPT(thermostat)) {
                Q=(sum.V/iblock)/(HISTMAX*2);
                q=sqrt(sum.V/iblock)/HISTMAX; }
              else {
                Q=L*L/(HISTMAX*2);
                q=L/HISTMAX; }
              shift=(1-HISTMAX)*0.5; // y-profile is centered
              if (!files.csv) fprintf(out,"\
Vertical density profile\n\
------------------------------\n\
y\tρ(y)\tN(y)\tKEY\n");
              break;
            case CPROFILE:
              dpkey="CRDP";
              goto qqq;
            default:
              dpkey="DRDP";
            qqq:
              if (isNPT(thermostat)) {
                q=sqrt(sum.V/iblock)/(HISTMAX*2);
                Q=sum.V/iblock*(PI/Sqr(HISTMAX)/8); }
              else {
                q=L/(HISTMAX*2);
                Q=PI*Sqr(L/HISTMAX)/8; }
              if (!files.csv) fprintf(out,"\
Radial density profile\n\
------------------------------\n\
r\tρ(r)\tN(r)\tKEY\n"); }

          if (files.csv) {
            sprintf(files.ext,"-%s%d.csv",dpkey,files.nmeas);
            dpcsv=fl_fopen(files.protocol,"wt");
            if (!dpcsv) {
              fl_alert("cannot write to CSV file (will use TXT)");
              return 0; }
            fprintf(dpcsv,"#%c,%s,cumul\n",dpkey[0]=='V'?'y':'r',dpkey); }

          loop (i,0,HISTMAX) {
            if (measure==YPROFILE) cumul+=(rhosum[i]/head->n);
            else cumul+=(rhosum[i]/head->n)*(2*i+1);
            if (dpcsv) sprintf(lout,"%g\t%g\t%g\n",(i+shift)*q,rhosum[i]/head->n,cumul*Q);
            else sprintf(lout,"%g\t%g\t%g\t%s%d\n",(i+shift)*q,rhosum[i]/head->n,cumul*Q,dpkey,files.nmeas);
            if (measure==YPROFILE) cumul+=(rhosum[i]/head->n);
            else cumul+=(rhosum[i]/head->n)*(2*i+1);
            comma(lout);

            if (dpcsv) {
              csvline(lout);
              fprintf(dpcsv,"%s",lout); }
            else
              fputs(lout,out); }
          fclose(dpcsv); } // density profiles

        if (files.csv) fclose(files.csv);
        fclose(out); }

    case 2: // [clear]
      // erasing recorded data
      for (a=head; a; a=next) {
        next=a->next;
        // delete a;
        free(a); }

      int i;
      loop (i,0,HISTMAX) rhosum[i]=0;

      head=NULL;
      return 0; } /* switch (choice) */

  fprintf(stderr,"fl_choice: bad return value\n");

  return 0;
}
