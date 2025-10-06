// Data recording - activated by radio buttom [▯ record]

// hold the variable names and accumulated sums and sums of squared over blocks
struct res_s {
  struct res_s *next;
  const char *name;
  int n;
  double sum,sumq;
} *head;

// to store the convergence profiles (CP)
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

/* Add 1 measurement to a registered variable for to calculate a block average
   'name' should not contain , ; Tab space
   if ((incl->value()&1)), then also append a line to CP:
     the names are printed to CP in reverse order
   MDS is an exception, it is calculated at block end
*/
int addM(const char *name,double x) /********************************** addM */
{
  struct res_s *a;
  static char line[1024]; // WARNING, fixed size!
  char item[32];
  int litem;

  if ((incl->value()&1)) {
    // convergence profiles
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

      return 0; }

    litem=sprintf(item,"%g\t",x);
    memmove(line+litem,line,strlen(line)+1);
    memcpy(line,item,litem); }

  if (!name) return 0;

  looplist (a,head) if (!strcmp(name,a->name)) {
    if (memcmp(name,"MSD",3)) {
      a->sum+=x;
      a->sumq+=x*x;
      a->n++; }
    else {
      a->sum=x; // exceptional for MSDx,MSDy
      a->sumq=t; 
      a->n=1; }

    goto ret; }

  a=(struct res_s*)malloc(sizeof(struct res_s));
  // a=new struct res_s;
  a->next=head;
  a->name=name;
  head=a;
  a->n=1;
  a->sum=x;
  a->sumq=x*x;

 ret:
  //  fprintf(stderr,"%d %s %g\n",a->n,a->name,x);
  return a->n;
} // addM()

// change decimal points to commas in place if requested
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
// + set switch measureVar
int showclearM() /*********************************************** showclearM */
{
  struct res_s *a,*next;
  char *end;
  char butt0[80];
  int choice,i;
  double shift=0.5;
  FILE *out;
  struct CP_s *CP;

  if (!buttons.csv->value()) files.csv=NULL;

  files.sep=",;"[(int)buttons.comma->value()]; // global

  if (!head) return 0; // recording already off

  // recording just has turned off

  // generate info string to pop up and write to a file
  // also, csv generated from this (cumbersome...)
  end=files.info;

  looplist (a,head) {
    double av=a->sum/a->n;
    double err=(a->sumq/a->n-Sqr(av))/(a->n-1);
    char serr[32];

    if (!memcmp("MSD",a->name,3))
      sprintf(serr," at %g (last time)",a->sumq); // do not change - parsed
    else if (a->n>1 && err>0) {
      err=sqrt(err);
      if (fabs(err)>=fabs(av))
        sprintf(serr," ± %.3g",err);
      else
        sprintf(serr," ± %.3g (%.3g %%)",err,fabs(err/av*100)); }
    else
      serr[0]=0;

    if (a==head) end+=sprintf(end,"\
%2d blocks (block length=%d, stride=%d)\n\
----------------------------------------\n\
quantity = average ± standard error (relative standard error)\n\
",a->n,block,speed.stride);

    end+=sprintf(end,"%s = %g%s\n", a->name, av,serr); }

  comma(files.info);

  if (head->n<2)
    choice=!fl_choice("Not enough blocks have been recorded.","Stop recording","Continue",NULL);
  else {
    if (files.nmeas==0) {
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
      files.record0=1;
      return 1;

    case 1: // [save]
      strcpy(files.ext,".txt");
      out=fl_fopen(files.protocol,files.nmeas==0?"wt":"at");
      if (buttons.csv->value()) {
        strcpy(files.ext,".csv");
        files.csv=fl_fopen(files.protocol,files.nmeas==0?"wt":"at");
        if (!files.csv) fl_alert("cannot write the summary files.csv file"); }
      else
        files.csv=NULL;

      if (!out)
        fl_alert("cannot write the protocol");
      else {
        char lout[512]; // 2* safety size

        files.nmeas++;

        if (files.nmeas==1)
          fprintf(out,"SIMOLANT VERSION %s\n",VERSION);

        sprintf(lout,"\
\n\
=========== MEASUREMENT =========== # %d ===========\n\
%s\n\
Parameters:\n\
  N=%d number of particles\n\
  bc=%s boundary conditions\n\
  method=\"%s\"\n\
  Nf=%d degrees of freedom\n\
  qx,qy=%g,%g kinetic pressure correction factors\n\
  T=%g temperature -> B=%g (2nd virial coefficient)\n\
  P=%g pressure\n\
  g=%g gravity\n\
  drop.v=%g drop.T=%g v and T for two droplets\n",
                files.nmeas,
                FFinfo(1),
                N,
                bc==BOX?"Box":bc==SLIT?"Slit":"Periodic",
                methodinfo[method],
                En.Nf,
                En.q.x,En.q.y,
                T,ss.B2,
                P,
                gravity,
                drop.v,drop.T);
        comma(lout);
        fputs(lout,out);

        if (files.csv) {
          if (files.nmeas==1)
            // cf. info2csv below
            fprintf(files.csv,"#,N,bc,method,T,P,g,walls,rho_wall,L,rho,c,a,b,tau,qtau,dt|d,dV,speed.stride,block,Etot,err,Tkin,err,Epot,err,V,err,Z,err,Pvir,err,Pxx,err,Pyy,err,γ,err,enthalpy,err,P(right_wall),err,P(left_wall),err,P(top_wall),err,P(bottom_wall),err,Econserved,err,MSDx,at_t,MSDy,at_t,momentum\n");
          fprintf(files.csv,"%d",files.nmeas);
          prtcsv(N);
          prtcsv(bc);
          fprintf(files.csv,"%c\"%s\"",files.sep,space2_(methodinfo[method])); // do not use prtcsv(char*)
          prtcsv(T);
          prtcsv(P);
          prtcsv(gravity);
          prtcsv(walls);
          prtcsv(walldens/PI);
          prtcsv(L);
          prtcsv(L2rho(L));
          prtcsv(ss.C2);
          prtcsv(ss.a);
          prtcsv(ss.b);
          prtcsv(tau);
          prtcsv(qtau);
          if (isMC(method)) prtcsv(dt);
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
  c=%g a=%g b=%g potential cutoff and parameters\n",L,L2rho(L),ss.C2,ss.a,ss.b);
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
        // fprintf(stderr,"strlen(files.info)=%d\n",(int)(strlen(files.info)));
        // ~650 bytes detected ⇒ files.info: char info[1024] = safe size

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
          info2csv("H");
          info2csv("P(right_wall)");
          info2csv("P(left_wall)");
          info2csv("P(top_wall)");
          info2csv("P(bottom_wall)");
          info2csv("Econserved");
          info2csv("MSDx");
          info2csv("MSDy");
          info2csv("momentum");
          fprintf(files.csv,"\n"); }

        if (measureVar>1) {
          fprintf(out,"\n----------------------------------------\nVariances (%d points, stride=%d)\n----------------------------------------\n",measureVar,speed.stride);
          fprintf(out,"Var(V) = %g \n",(sumq.V-sums.V/measureVar)/(measureVar-1));
          fprintf(out,"Var(H) = %g \n",(sumq.H-sums.H/measureVar)/(measureVar-1));
          fprintf(out,"Var(Upot) = %g \n",(sumq.Upot-sums.Upot/measureVar)/(measureVar-1));
          fprintf(out,"Var(%s) = %g \n",isMC(method)?"Tbag":"Ekin",(sumq.Ekin-sums.Ekin/measureVar)/(measureVar-1)); }

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

        if (Show>=RDF && (incl->value()&2)) { // density profiles
          int i;
          double q,cumul=0,Q;
          const char *dpkey=NULL;

          FILE *dpcsv=NULL;

          if (!files.csv) fputs("------------------------------\n",out);

          switch (Show) {
            case RDF:
              dpkey="RDF";
              q=1./50; // grid per 1
              if (isNPT(method)) Q=L2rho(sqrt(sum.V/iblock));
              else Q=L2rho(L);
              Q*=PI/5000; /* 50^2*2 */
              if (!files.csv) fprintf(out,"\
Radial distribution function\n\
------------------------------\n\
r\tg(r)\tN(r)\tKEY\n");
              break;
            case YPROFILE:
              dpkey="VDP";
              if (isNPT(method)) {
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
              if (isNPT(method)) {
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

          /* print RDF and density profiles to string lout */
          loop (i,0,HISTMAX) {
            if (Show==YPROFILE) cumul+=(rhosum[i]/head->n);
            else cumul+=(rhosum[i]/head->n)*(2*i+1);
            if (dpcsv) sprintf(lout,"%g\t%g\t%g\n",(i+shift)*q,rhosum[i]/head->n,cumul*Q);
            else sprintf(lout,"%g\t%g\t%g\t%s%d\n",(i+shift)*q,rhosum[i]/head->n,cumul*Q,dpkey,files.nmeas);
            if (Show==YPROFILE) cumul+=(rhosum[i]/head->n);
            else cumul+=(rhosum[i]/head->n)*(2*i+1);

            comma(lout); // change . to ,

            if (dpcsv) {
              csvline(lout); // change tab CSV separator
              fprintf(dpcsv,"%s",lout); }
            else
              fputs(lout,out); }
          fclose(dpcsv); } // density profiles

        if (files.csv) fclose(files.csv);
        fclose(out); }

    case 2: // [clear]
      measureVar=0;
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
