// fltk-config --use-images --compile simolant.cc
// g++ -O2 -o simolant simolant.cc -lfltk -lfltk_images

int record;
const char *key;

struct res_s {
  struct res_s *next;
  const char *name;
  int n;
  double sum,sumq;
} *head;

#define looplist(PTR,HEAD) for ((PTR)=(HEAD); (PTR); (PTR)=(PTR)->next)

// to store all the CP
struct CP_s {
  struct CP_s *next;
  char line[1]; // variable length
} *CPhead=NULL;

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

void comma(char *s) /************************************************* comma */
{
  if (buttons.comma->value())
    while (*s) {
      if (*s=='.') *s=',';
      s++; }
}

int showclearM() /*********************************************** showclearM */
/* pops the results, asks, returns the record status (1=continue, 0=erased) */
{
  struct res_s *a,*next;
  char info[2048],*end;
  char butt0[64];
  static int overwrite=1,i,firstout=1;
  int choice;
  double shift=0.5;
  FILE *out;
  struct CP_s *CP;

  if (!head) return 0; // recording already off

  // recording just turned off

  // generate info string to pop up and write to a file
  end=info;

  looplist (a,head) {
    double av=a->sum/a->n;
    double stder=(a->sumq/a->n-Sqr(av))/(a->n-1);

    if (a->n<2 || stder<=0) stder=0;
    else stder=sqrt(stder);

    if (a==head) end+=sprintf(end,"\
%2d blocks (block length=%d, stride=%d)\n\
----------------------------------------\n\
quantity = average ± standard error (relative std.err.)\n\
",a->n,block,nit);

    end+=sprintf(end,"%s = %g ± %.3g (%.3g %%)\n",
                 a->name,
                 av,stder,fabs(stder/av*100)); }
  comma(info);

  if (head->n<2)
    choice=!fl_choice("Not enough blocks have been recorded.","Stop recording","Continue",NULL);
  else {
    if (overwrite) {
      if (strlen(protocol)<32) sprintf(butt0,"save (overwrite \"%s\") and clear",protocol);
      else strcpy(butt0,"save (overwrite protocol) and clear"); }
    else {
      if (strlen(protocol)<32) sprintf(butt0,"append to \"%s\" and clear",protocol);
      else strcpy(butt0,"append to protocol and clear"); }

    // It is not possible to use info as the 1st parameter
    // because the paranoic compiler shouts that
    // "format not a string literal and no format arguments"
    choice=fl_choice("%s","continue",butt0,"clear (data lost)",info);
  }

  switch (choice) {
    case 0:
      return 1;

    case 1:
      if (overwrite) {
        char *bak=strdup(protocol),*ext=(char*)getext(bak);
        if (!strcmp(ext,".txt")) strcpy(ext,".bak");
        rename(protocol,bak);
        free(bak);
        out=fl_fopen(protocol,"wt"); }
      else
        out=fl_fopen(protocol,"at");
      if (out) {
        static int BLOCK=0;
        static char lout[1024];

        BLOCK++;
        overwrite=0;

        if (firstout) {
          fprintf(out,"SIMOLANT VERSION %s\n",VERSION);
          firstout=0; }

        sprintf(lout,"\
\n\
=========== MEASUREMENT =========== # %d ===========\n\
\n\
N=%d bc=%s method=%s T=%g g=%g P=%g\n",
                BLOCK,
                N,bc==BOX?"Box":bc==SLIT?"Slit":"Periodic",
                thinfo[thermostat],T,gravity,P);
        comma(lout); fputs(lout,out);

        sprintf(lout,"walls: left=%s right=%s bottom=%s top=%s rho_wall=%g\n",
                walls&1?"attr":"rep",
                walls&2?"attr":"rep",
                walls&4?"attr":"rep",
                walls&8?"attr":"rep",walldens/PI);
        comma(lout); fputs(lout,out);

        sprintf(lout,"L=%g  rho=%g  cutoff=%g\n",L,L2rho(L),cutoff);
        comma(lout); fputs(lout,out);

        sprintf(lout,"tau=%g  qtau=%g  dt=%g\n",tau,qtau,dt);
        comma(lout); fputs(lout,out);

        sprintf(lout,"d=%g(%s)  dV=%g  stride*block=%d*%d\n",d,setd->value()?"auto":"fixed",dV,nit,block);
        comma(lout); fputs(lout,out);


        fprintf(out,"\n----------------------------------------\n");
        fputs(info,out);
        //        fprintf(stderr,"strlen(info)=%d\n",(int)(strlen(info)));
        // <600 bytes detected, char info[2048] = safe size

        if (CPhead && (incl->value()&1)) {

          fprintf(out,"\n\
----------------------------------------------------\n\
Convergence profile (time development of quantities)\n\
----------------------------------------------------\n");

          fprintf(out,"A=time");
          i='A';
          looplist (a,head) {
            i++;
            fprintf(out,"\t%c=%s",i,a->name); }
          fprintf(out,"\tKEY\n");

          for (CP=CPhead; CP; ) {
            struct CP_s *next=CP->next;
            comma(CP->line);
            fprintf(out,"%sCP%d\n",CP->line,BLOCK);
            free(CP);
            CP=next; }

          CPhead=NULL; }

        fputs("\n",out);

        if (measure>=RDF && (incl->value()&2)) {
          int i;
          double q,cumul=0,Q;

          fputs("------------------------------\n",out);
          switch (measure) {
            case RDF:
              key="RDF";
              q=1./50;
              if (isNPT(thermostat)) Q=L2rho(sum.V/iblock);
              else Q=L2rho(L);
              Q*=PI/5000; /* 50^2*2 */
              fprintf(out,"\
Radial distribution function\n\
------------------------------\n\
r\tg(r)\tN(r)\tKEY\n");
              break;
            case YPROFILE:
              key="VDP";
              if (isNPT(thermostat)) {
                Q=(sum.V/iblock)/(HISTMAX*2);
                q=sqrt(sum.V/iblock)/HISTMAX; }
              else {
                Q=L*L/(HISTMAX*2);
                q=L/HISTMAX; }
              shift=(1-HISTMAX)*0.5; // y-profile is centered
              fprintf(out,"\
Vertical density profile\n\
------------------------------\n\
y\tρ(y)\tN(y)\tKEY\n");
              break;
            case CPROFILE:
              key="CRDP";
              goto qqq;
            default:
              key="DRDP";
            qqq:
              if (isNPT(thermostat)) {
                q=sqrt(sum.V/iblock)/(HISTMAX*2);
                Q=sum.V/iblock*(PI/Sqr(HISTMAX)/8); }
              else {
                q=L/(HISTMAX*2);
                Q=PI*Sqr(L/HISTMAX)/8; }
              fprintf(out,"\
Radial density profile\n\
------------------------------\n\
r\tρ(r)\tN(r)\tKEY\n"); }

          loop (i,0,HISTMAX) {
            if (measure==YPROFILE) cumul+=(rhosum[i]/head->n);
            else cumul+=(rhosum[i]/head->n)*(2*i+1);
            sprintf(lout,"%g\t%g\t%g\t%s%d\n",(i+shift)*q,rhosum[i]/head->n,cumul*Q,key,BLOCK);
            if (measure==YPROFILE) cumul+=(rhosum[i]/head->n);
            else cumul+=(rhosum[i]/head->n)*(2*i+1);
            comma(lout);
            fputs(lout,out); } }

        fclose(out); }

    case 2:
      for (a=head; a; a=next) {
        next=a->next;
        // delete a;
        free(a); }

      int i;
      loop (i,0,HISTMAX) rhosum[i]=0;

      head=NULL;
      return 0; }

  fprintf(stderr,"fl_choice: bad return value\n");
  return 0;
}
