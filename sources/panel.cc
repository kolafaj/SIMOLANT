// fltk-config --use-images --compile simolant.cc
// g++ -O2 -o simolant simolant.cc -lfltk -lfltk_images

/* shorten the info string */
void shorten(char *s) /********************************************* shorten */
{
  char *c;
  int i;

  loop (i,0,8) {
    if (fl_width(s)<PANELW-14) break;
    if ( (c=strstr(s,"N/⟨V⟩=")) ) memcpy(c,c+strlen("N/⟨V⟩="),strlen(c));
    if ( (c=strstr(s,"N/<V>=")) ) memcpy(c,c+strlen("N/<V>="),strlen(c));
    if ( (c=strstr(s,"acc.r.")) ) memcpy(c+1,c+3,strlen(c));
    if ( (c=strstr(s,"gamma")) ) { memcpy(c+1,c+4,strlen(c)); c[1]='m'; }
    if (fl_width(s)<PANELW-14) break;
    if ( (c=strstr(s,"Tbag")) ) { memcpy(c+2,c+3,strlen(c)); }
    if ( (c=strstr(s,"  ")) ) memcpy(c+1,c+2,strlen(c));
    if ( (c=strstr(s,"e-0")) ) memcpy(c+2,c+3,strlen(c));
    if ( (c=strstr(s,"e+0")) ) memcpy(c+2,c+3,strlen(c));
    if ( (c=strstr(s,"=-0.")) ) memcpy(c+2,c+3,strlen(c));
    if ( (c=strstr(s,"=0.")) ) memcpy(c+1,c+2,strlen(c)); }
}

/* top: parameters and run info */
class InfoPanel : public Fl_Box
{
private:
  void draw () {
    char s[256]; /* not to crash for huge number and fixed format */
    static char info[32];
    int atx=x()+6;
#define TSIZE 20 // font size for the top line, otherwise 16
    int aty=y()+TSIZE;
    static struct last_s {
      double T,P,tau,d;
      double L,rho;
      double gravity;
      int N,nit,block;
    } last;

    static struct lastc_s {
      int N,Td,d,Tt,L,g;
    } lastc;

    Fl_Box::draw();

    fl_font(FL_HELVETICA,TSIZE);
    if (thermostat==CREUTZ||thermostat==NVE) fl_color(0,127,0);
    else if (isNPT(thermostat)) fl_color(150,0,160);
    else fl_color(0,100,150);

    /* basic info */
    if (mauto) {
      if (thermostat==METROPOLIS)
        strcpy(info,thinfo[AUTO]);
      else
        sprintf(info,"[%s]",thinfo[thermostat]); }
    else
      strcpy(info,thinfo[thermostat]);

    // something wrong here - FL_ALIGN_BOTTOM out of order, always y-centered
    //    fl_draw(info,atx,aty,PANELW-9+3*(info[0]=='['),TSIZE,FL_ALIGN_RIGHT|FL_ALIGN_BOTTOM);
    // using fl_width(info) is much easier
    fl_draw(info,atx+PANELW-11+3*(info[0]=='[')-fl_width(info),aty);

    fl_color(FL_BLACK);

    if (N!=last.N) lastc.N=10;
    if (lastc.N) fl_color(FL_RED),lastc.N--;
    last.N=N;
    sprintf(s,"N=%d",N); fl_draw(s,atx,aty);
    fl_color(FL_BLACK);
    fl_font(FL_HELVETICA,16);

    aty+=TSIZE+4;

    switch (thermostat) {
      case CREUTZ:
        if (d!=last.d) lastc.d=10;
        if (lastc.d) fl_color(FL_RED),lastc.d--;
        sprintf(s,"E=%.3f  d*Lh=%5.3f",Econserved,d*Lh);
        break;
      case METROPOLIS:
        if (T!=last.T || d!=last.d) lastc.Td=10;
        if (lastc.Td) fl_color(FL_RED),lastc.Td--;
        sprintf(s,"T=%6.4f d*Lh=%5.3f",T,d*L);
        break;
      case MCNPT:
        if (T!=last.T || d!=last.d || P!=last.P) lastc.Td=10;
        if (lastc.Td) fl_color(FL_RED),lastc.Td--;
        sprintf(s,"T=%6.4f d=%.3f dV=%.3f P=%g",T,d,dV,P);
        break;
      case NVE:
        sprintf(s,"Epot+Ekin=%.3f",(sum.Ekin+sum.U)/iblock);
        //        sprintf(s,"Epot+Ekin=%.3f",Econserved); last step
        break;
      case MDNPT:
      case MTK:
        if (T!=last.T || tau!=last.tau || P!=last.P) lastc.Td=10;
        sprintf(s,"T=%6.4f P=%5.3f τ=%5.3f qτ=%g",T,P,tau,qtau);
        break;
      default:
        if (T!=last.T || tau!=last.tau) lastc.Tt=10;
        if (lastc.Tt) fl_color(FL_RED),lastc.Tt--;
        sprintf(s,"T=%6.4f  τ=%5.3f",T,tau); }

    last.T=T;
    last.P=P;
    last.tau=tau;
    last.d=d;
    fl_draw(s,atx,aty); aty+=20;
    fl_color(FL_BLACK);

    if (thermostat==MCNPT || thermostat==MDNPT || thermostat==MTK) {
      fl_color(FL_RED);
      if (bc==PERIODIC) sprintf(s,"L=%6.3f",L);
      else sprintf(s,"L=%6.3f wall=%.3g",L,walldens/PI); }
    else {
      if (L!=last.L) lastc.L=10;
      if (lastc.L) fl_color(FL_RED),lastc.L--;
      last.L=L;
      if (bc==PERIODIC) sprintf(s,"L=%-.4g ρ=%.4g",L,L2rho(L));
      else sprintf(s,"L=%-.4g ρ=%.4g wall=%.3g",L,L2rho(L),walldens/PI); }

    fl_draw(s,atx,aty); aty+=20;
    fl_color(FL_BLACK);

    if (h()>84) {
      if (gravity!=last.gravity || nit!=last.nit || block!=last.block) lastc.g=10;
      if (lastc.g) fl_color(FL_RED),lastc.g--;
      last.gravity=gravity;
      last.nit=nit;
      last.block=block;
      sprintf(s,"g=%5.3f    stride*block=%s%d*%d",gravity,sliders.speed->value()<=0.48?"(slow)":"",nit,block); fl_draw(s,atx,aty); }
    fl_color(FL_BLACK);
  }

public:
  InfoPanel (int X, int Y, int W, int H) : Fl_Box(X,Y,W,H) {
    color(FL_WHITE);
    box(FL_UP_BOX);
  }
  ~InfoPanel() { }
};

/* results: numbers and graph */
class ResultPanel : public Fl_Box
{
private:
  void draw () {
    /* (blocked) results print/draw */
    if (iblock>=block) {
      char s[64],ss[64];
      const char *yaxis;
      int atx=x()+6,aty=y()+18;
      int i;
      double CPval=0;

      Fl_Box::draw();

      //      fprintf(stderr,"ResultPanel:draw %g record=%d\n",t,record);

      if (record) {
        if (bc<PERIODIC) {
          addM("P(top_wall)",sum.fwy/iblock);
          addM("P(bottom_wall)",sum.fwyL/iblock); }
        if (bc==BOX) {
          addM("P(left_wall)",sum.fwx/iblock);
          addM("P(right_wall)",sum.fwxL/iblock); }
        addM("γ",(sum.Pyy-sum.Pxx)/iblock*L);
        addM("Pyy",sum.Pyy/iblock);
        addM("Pxx",sum.Pxx/iblock);
        addM("Pvir",(sum.Pxx+sum.Pyy)/2/iblock);
        addM("Z",(sum.Pxx+sum.Pyy)/(2*iblock*T*L2rho(L)));
        addM("V",sum.V/iblock);
        addM("Epot",sum.U/iblock); }

      fl_font(FL_HELVETICA,16);

      fl_color(FL_BLUE);

      const char *eq=dtfixed?"=":"⇒";

      if (thermostat==MDNPT || thermostat==MTK) {
        if (record) sprintf(s,"Tkin=%5.3f  dt%s%5.4f  ρ=N/⟨V⟩=%.4f  n=%d",sum.Tk/iblock,eq,dt,N/(sum.V/iblock),head->n);
        else sprintf(s,"Tkin=%5.3f  dt%s%5.4f  ρ=N/⟨V⟩=%.4f",sum.Tk/iblock,eq,dt,N/(sum.V/iblock));
        if (record) addM("Tkin",sum.Tk/iblock); }
      else if (isMD(thermostat)) {
        if (record) sprintf(s,"Tkin=%5.3f  dt%s%5.4f  n=%d",sum.Tk/iblock,eq,dt,head->n);
        else sprintf(s,"Tkin=%5.3f  dt%s%5.4f",sum.Tk/iblock,eq,dt);
        if (record) addM("Tkin",sum.Tk/iblock); }
      else if (thermostat==MCNPT) {
        if (record) sprintf(s,"acc.r.=%.4f  ρ=N/⟨V⟩=%.4f  n=%d",sum.Var/iblock,N/(sum.V/iblock),head->n);
        else sprintf(s,"acc.r.=%5.3f  ρ=N/⟨V⟩=%.4f",sum.ar/iblock,N/(sum.V/iblock)); }
      else { /* Metropolis, CREUTZ */
        if (record) sprintf(s,"Tbag=%5.3f  acc.r.=%5.3f  n=%d",sum.Ekin/iblock,sum.ar/iblock,head->n);
        else sprintf(s,"Tbag=%5.3f  acc.r.=%5.3f",sum.Ekin/iblock,sum.ar/iblock);
        if (record) addM("Tbag",sum.Ekin/iblock); }

      /* conserved energy  */
      if (thermostat==NOSE_HOOVER || thermostat==MTK)
        if (record) addM("Econserved",(sum.Econserved+sum.U)/iblock);

      if (isNPT(thermostat))
        if (record) addM("H",sum.H/iblock);

      shorten(s);
      fl_draw(s,atx,aty);
      aty+=20;

      // MC: equivalent kinetic energy to print internal energy
      double eqEkin=T*degrees_of_freedom()*iblock/2.;

      /* convergence profile: interpretation of menu items */
      yaxis=NULL;
      switch (measure) {
        
        case ENERGY:
          if (isNPT(thermostat)) {
            yaxis="enthalpy";
            CPval=(sum.U+eqEkin)/iblock+sum.V/iblock*P; }
          else if (thermostat==NVE || thermostat==CREUTZ) {
            yaxis="Epot";
            CPval=sum.U/iblock; }
          else if (isMC(thermostat)) {
            yaxis="Epot+f/2*kT";
            CPval=(sum.U+eqEkin)/iblock; }
          else {
            yaxis="Epot+Ekin";
            CPval=(sum.U+sum.Ekin)/iblock; }
          break;

        case TEMPERATURE:
          if (isMC(thermostat)) {
            yaxis="Tbag";
            CPval=sum.Ekin/iblock; /* this is the bag */ }
          else {
            yaxis="Tkin";
            CPval=sum.Ekin/degrees_of_freedom()*2/iblock; }
          break;

        case PRESSURE:
          yaxis="Pvir";
          CPval=(sum.Pxx+sum.Pyy)/2/iblock;
          break;

        case VOLUME:
          if (isNPT(thermostat)) {
            yaxis="volume";
            CPval=sum.V/iblock; }
          break;

        case INTMOTION:
          switch (thermostat) {
            case MTK:
              yaxis="ext.enthalpy";
              CPval=sum.Econserved/iblock;
              break;
            case NOSE_HOOVER:
              yaxis="ext.energy";
              CPval=sum.Econserved/iblock;
              break;
            case NVE:
              yaxis="Ekin+Epot";
              CPval=(sum.Ekin+sum.U)/iblock;
              break;
            default:; /* to make the compiler happy */
          }
          
        default:; /* to make the compiler happy */
      }

      if (itot<HISTMAX) {
        etot[itot++]=CPval;
        itot%=HISTMAX; }
      else if (itot>HISTMAX) itot--;
      else if (itot==HISTMAX) {
         etot[--itot]=CPval;
         while (itot--) etot[itot]=-9e99; }

      if (record) {
        addM("Etot",(sum.U+sum.Ekin)/iblock);
        addM(NULL,0); // end of line - must be the last of all addM's
      }

      double Pav=(sum.Pxx+sum.Pyy)/2/iblock,Pid=T*L2rho(L),Z=Pav/Pid;

      if (measure==NONE)
        sprintf(ss,"Pid=%.4g",Pid);
      else {
        sprintf(ss,"Pvir=%.4g  Pid=%.4g  Z=%.4g",Pav,Pid,Z);
        shorten(ss); }

      if (measure==YPROFILE) {
        sprintf(s,"Pxx=%.3g  Pyy=%.3g  γ=%.3f",
                sum.Pxx/iblock,sum.Pyy/iblock,
                (sum.Pyy-sum.Pxx)/iblock/2*L);
        shorten(s);
        fl_draw(s,atx,aty); }
      else if (measure==QUANTITIES) {
        fl_draw(ss,atx,aty);
        aty+=20;
        if (bc==BOX) {
          sprintf(s,"P(left wall)=%.4g",sum.fwx/iblock);
          fl_draw(s,atx,aty);
          aty+=20;
          sprintf(s,"P(right wall)=%g",sum.fwxL/iblock);
          fl_draw(s,atx,aty);
          aty+=20; }
        if (bc<=SLIT) {
          sprintf(s,"P(top wall)=%.4g",sum.fwy/iblock);
          fl_draw(s,atx,aty);
          aty+=20;
          sprintf(s,"P(bottom wall)=%g",sum.fwyL/iblock);
          fl_draw(s,atx,aty);
          aty+=20; }

        /* WARNING: do not use fix format - overflow if P huge */
        sprintf(s,"Pxx=%6.3g  Pyy=%.3g  γ=%.3g",
                sum.Pxx/iblock,sum.Pyy/iblock,
                (sum.Pyy-sum.Pxx)/iblock/2*L);
        shorten(s);
        fl_draw(s,atx,aty); }
      else
        fl_draw(ss,atx,aty);

      aty+=20;

      if (measure!=NONE) {
        sprintf(s,"Etot=%g  Epot=%g",(sum.U+sum.Ekin)/iblock,sum.U/iblock);
        fl_draw(s,atx,aty);
        aty+=20;
      }

      if (shownbrs) {
        fl_color(FL_BLACK);
        aty+=4;
        fl_draw("neighbors:",atx,aty);
        atx+=86;
        int dx=26;

        loop (i,0,8) {
          fl_color(nbrcol[i]);
          fl_pie(atx+i*dx,aty-18,24,24,0.0,360.0); }
        fl_color(FL_WHITE);
        loop (i,0,7) {
          char s[8];
          sprintf(s,"%d",i);
          fl_draw(s,atx+i*dx+7,aty); }
        fl_color(FL_BLACK);
        fl_draw("7+",atx+i*dx+2,aty); }
      else if (isNPT(thermostat)) {
        if (thermostat==MCNPT) sprintf(s,"H=%g acc.r.(V)=%.2f Tbag=%.3f",sum.H/iblock,sum.Var/iblock,sum.Ekin/iblock);
        else sprintf(s,"H=%g",sum.H/iblock);
        shorten(s);
        fl_draw(s,atx,aty); }

      aty=y()+106;
      atx=x()+41;

      fl_color(FL_BLUE);

      if (h()>180) switch (measure) {

          case PRESSURE:
          case INTMOTION:
          case VOLUME:
          case ENERGY:
          case TEMPERATURE:
            if (!yaxis) {
              fl_color(FL_RED);
              fl_draw(minfo[measure],atx-20,aty+20);
              fl_draw("is not available",atx-20,aty+40);
              break; }

            if (itot<HISTMAX) {
              int i,netot=0;
              double emin=1e99,emax=-1e99;
              double sume=0,sq=0;
              double etot0=0;
              static double lmin,lmax;


              loop (i,0,HISTMAX) if (etot[i]>-8e99) { etot0=etot[0]; break; }

              loop (i,0,HISTMAX) if (etot[i]>-8e99) {
                netot++;
                sume+=etot[i]-etot0;
                sq+=Sqr(etot[i]-etot0);
                if (etot[i]<emin) emin=etot[i];
                if (etot[i]>emax) emax=etot[i]; }

              if (justreset) lmin=emin,lmax=emax,justreset=1;

              if (lmin<emin) lmin=0.9*lmin+0.1*emin; else lmin=emin;
              if (lmax>emax) lmax=0.9*lmax+0.1*emax; else lmax=emax;
              emin=lmin;
              emax=lmax;

              fl_color(FL_BLACK);
              sprintf(s,"%.6g",emax);
              fl_draw(s,atx-30,aty-5);
              sprintf(s,"%.6g",emin);
              fl_draw(s,atx-30,aty+120);
              if (emax-emin<1e-9) sprintf(s,"max-min=%.2g",emax-emin);
              else if (emax-emin<1e-2) sprintf(s,"max-min=%.3g",emax-emin);
              else sprintf(s,"max-min=%.4g",emax-emin);
              fl_draw(s,atx+100,aty+120);

              // quantity to show: should correspond to CPval above
              fl_draw(90,yaxis,atx-18,aty+55+4*strlen(yaxis));

              if (emax<=emin) {
                fl_color(FL_RED);
                fl_draw(minfo[measure],atx,aty+40);
                fl_draw("nothing to show – zero range",atx,aty+60);
                break; }

              sq=sqrt((sq-Sqr(sume)/netot)/(netot-1));
              sprintf(s,"stdev=%.3g",sq);
              fl_draw(s,atx+100,aty-5);

              //          fl_color(FL_DARK_GREEN); =(0,145,0)
              fl_color(fl_rgb_color(0,120,0));
              fl_line_style(FL_SOLID,0,NULL);
              fl_rectf(atx,aty,HISTMAX,102);

              /* draw grid */
              double demax=102/(emax-emin);
              double ex=log(emax-emin)/log(2.);
              double dl=floor(ex);
              double fr=ex-floor(ex);
              dl=pow(2.,dl-1.5); // -> number of shaded lines
              fl_color(fl_rgb_color(0,190,0));

              ex=floor(emin/dl)*dl;
              dl/=2;
              while (ex<=emax) {
                int y=102.5-(ex-emin)*demax;
                if (ex>=emin) {
                  fl_color(fl_rgb_color(0,220,0));
                  fl_line(atx,aty+y,atx+HISTMAX,aty+y); }
                ex+=dl;
                if (ex>emax) break;
                if (ex>=emin) {
                  y=102.5-(ex-emin)*demax;
                  fl_color(fl_rgb_color(0,(int)(220-100*fr),0));
                  fl_line(atx,aty+y,atx+HISTMAX,aty+y); }
                ex+=dl; }

              /* draw convergence profile */
              fl_color(FL_YELLOW);

              loop (i,0,HISTMAX) {
                int ii=(i+itot)%HISTMAX;
                if (etot[ii]>-8e99) {
                  int y=102.5-(etot[ii]-emin)*demax;
                  fl_line(atx+i,aty+y,atx+i,aty+102); } }

              fl_line_style(FL_SOLID,0,NULL); }

            break;

        case RDF: {
          double vv,rr;
          unsigned int ll,l,smooth;

          fl_color(FL_BLACK);

          fl_draw(90,"g(r)",x()+18,aty+70);

          loopto (i,0,5) {
            sprintf(s,"%d",i);
            fl_draw(s,atx-5+50*i,aty+115); }
          if (shownbrs) fl_draw("neighbor limit",atx+72,aty+130);
          else fl_draw("rmin",atx+54,aty+130);
          fl_draw("r",atx+220,aty+125);
          fl_draw("0",atx-13,aty+107);
          fl_draw("max",atx-37,aty+9);

          vv=Sqr(50*L/N)*2/PI/iblock;
          if (record) loop (i,0,HISTMAX) rhosum[i]+=hist[i]*vv/(2*i+1);

          loop (smooth,0,1+(block==1)) {
            l=hist[0];
            loop (i,0,HISTMAX) {
              ll=hist[i];
              hist[i]=l+ll+hist[i+1];
              l=ll; } }

          // determine maximum
          vv=0;
          loop (i,0,HISTMAX) {
            rr=(double)hist[i]/(2*i+1);
            if (rr>vv) vv=rr; }
          if (vv<=0) vv=1;
          vv=100.0/vv;

          fl_color(FL_DARK_GREEN);
          fl_line_style(FL_SOLID,0,NULL);
          fl_rectf(atx,aty,HISTMAX,100);

          fl_color(FL_WHITE);
          fl_draw("RDF",atx+3,aty+20);

          fl_color(FL_YELLOW);
          loop (i,1,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1);
            //            if (l)  // causes to draw (leave green) maximum, because l==0 at maximum
            fl_line(atx+i,aty+l,atx+i,aty+100);
            /* omit r=0,1,2,3,4,5: */
            if (i%50==48) i+=2; }

          // NB: it is very slow to switch colors - loops unrolled

          // vertical lines finished
          fl_color(FL_GREEN);
          loop (i,0,HISTMAX) {
            fl_line(atx+i,aty,atx+i,aty+100);
            if (i%50==0) i+=48; }
          if (shownbrs) i=(int)(RNBR*50); else i=59; /* rmin */
          fl_line(atx+i,aty,atx+i,aty+100);

          fl_color(FL_DARK_YELLOW);
          loop (i,0,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1);
            if (l<100) fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%50==0) i+=48; }

          if (shownbrs) {
            i=(int)(RNBR*50);
            l=100.5-vv*hist[i]/(2*i+1); }
          else {
            i=59;
            l=100.5-vv*hist[i]/(2*i+1); }

          fl_line(atx+i,aty+l,atx+i,aty+117);
          // fl_line(atx+i+2,aty+115,atx+i-1,aty+115);

          fl_line_style(FL_SOLID,0,NULL); }
          break;

        case YPROFILE: {
          double vv;
          unsigned int smooth,div=1;
          int l;

          fl_color(FL_BLACK);
          //          fl_draw("ρ(y)",x()+4,aty+35);
          fl_draw(90,"ρ(y)",x()+18,aty+75);
          fl_draw("0",atx-13,aty+107);
          //          fl_draw("1",atx-13,aty+8);
          fl_draw("0",atx-3,aty+115);
          fl_draw("L/2",atx+HISTMAX/2-10,aty+115);
          fl_draw("L",atx+HISTMAX-6,aty+115);
          fl_draw("y",atx+HISTMAX*14/19,aty+125);

          vv=HISTMAX/Sqr(L)/iblock;

          if (record) loop (i,0,HISTMAX) rhosum[i]+=hist[i]*vv;

          /* periodic smoothing */
          loop (smooth,0,1+(block==1)) {
            int h[HISTMAX];
            loopto (i,1,HISTMAX) h[i%HISTMAX]=hist[i-1]+hist[i%HISTMAX]+hist[(i+1)%HISTMAX];
            loop (i,0,HISTMAX) hist[i]=h[i];
            div*=3; }

          vv*=100.0/div;

          double maxh=0;
          loop (i,0,HISTMAX) if (hist[i]>maxh) maxh=hist[i];

          static double avzmax=1;
          static int izmax=9;
          double z=maxh*vv/100.;
          avzmax=0.8*avzmax+z*0.2;
          double zser[10]={0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,9e9};

          if (avzmax>zser[izmax]-.01) izmax++;
          if (z>zser[izmax]) izmax++;
          if (avzmax<zser[izmax-1]-0.01) izmax--;
          if (izmax>8) izmax=8;
          if (izmax<1) izmax=1;

          vv/=zser[izmax];

          sprintf(s,"%3.1f",zser[izmax]);
          fl_draw(s,atx-28,aty+8);

          fl_color(FL_DARK_GREEN);
          fl_line_style(FL_SOLID,0,NULL);
          fl_rectf(atx,aty,HISTMAX,100);

          fl_color(FL_WHITE);
          fl_draw("Vertical density profile, smoothed",atx+3,aty+20);

          fl_color(FL_YELLOW);
          loop (i,1,HISTMAX) {
            l=100.5-vv*hist[i];
            if (l<0) { fl_color(255,200,0); fl_line(atx+i,aty,atx+i,aty+100); fl_color(FL_YELLOW); }
            else fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%(HISTMAX/2)==HISTMAX/2-2) i+=2; }

          // NB: it is very slow to switch colors - loops unrolled
          fl_color(FL_GREEN);
          for (i=0; i<HISTMAX; i++) {
            fl_line(atx+i,aty,atx+i,aty+100);
            if (i%(HISTMAX/2)==0) i+=HISTMAX/2-2; }

          fl_color(FL_DARK_YELLOW);
          for (i=0; i<HISTMAX; i++) {
            l=100.5-vv*hist[i];
            if (l<0) { fl_color(200,180,0); fl_line(atx+i,aty,atx+i,aty+100); fl_color(FL_DARK_YELLOW); }
            else fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%(HISTMAX/2)==0) i+=HISTMAX/2-2; }

          fl_line_style(FL_SOLID,0,NULL);
          fl_color(FL_WHITE); }
          break;

        case CPROFILE:
        case RPROFILE: {
          double vv;
          int l;

          fl_color(FL_BLACK);
          //          fl_draw("ρ(r)",x()+4,aty+30);
          fl_draw(90,"ρ(r)",x()+18,aty+75);
          fl_draw("0",atx-13,aty+108);
          fl_draw("1",atx-13,aty+8);
          fl_draw("0",atx-3,aty+115);
          fl_draw("L/4",atx+HISTMAX/2-10,aty+115);
          fl_draw("L/2",atx+HISTMAX-15,aty+115);
          fl_draw("r",atx+HISTMAX*14/19,aty+125);

          vv=Sqr(HISTMAX/Lh)/PI/iblock;
          if (record) loop (i,0,HISTMAX) rhosum[i]+=hist[i]*vv/(2*i+1);

          fl_color(FL_DARK_GREEN);
          fl_line_style(FL_SOLID,0,NULL);
          fl_rectf(atx,aty,HISTMAX,100);

          fl_color(FL_WHITE);
          fl_draw(measure==CPROFILE?"radial density profile (cavity)":"radial density profile (droplet)",atx+3,aty+20);

          vv*=100.0;
          fl_color(FL_YELLOW);
          loop (i,1,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1);
            if (l<0) { fl_color(255,200,0); fl_line(atx+i,aty,atx+i,aty+100); fl_color(FL_YELLOW); }
            else fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%(HISTMAX/2)==HISTMAX/2-2) i+=2; }

          // NB: it is very slow to switch colors - loops unrolled
          fl_color(FL_GREEN);
          loop (i,0,HISTMAX) {
            fl_line(atx+i,aty,atx+i,aty+100);
            if (i%(HISTMAX/2)==0) i+=HISTMAX/2-2; }

          fl_color(FL_DARK_YELLOW);
          loop (i,0,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1)+0.5;
            if (l<0) { fl_color(200,180,0); fl_line(atx+i,aty,atx+i,aty+100); fl_color(FL_DARK_YELLOW); }
            else fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%(HISTMAX/2)==0) i+=HISTMAX/2-2; }

          fl_line_style(FL_SOLID,0,NULL); }
          break;
          
          default:; /* to make the compiler happy */ 
      }

      // for next cumulative sums:
      memset(hist,0,(HISTMAX+1)*sizeof(hist[0]));
      memset(&sum,0,sizeof(sum));
      iblock=0;
    } // iblock>=block
  } // draw()

public:
  ResultPanel (int X, int Y, int W, int H) : Fl_Box(X,Y,W,H) {
    color(FL_WHITE);
    box(FL_UP_BOX);
  }
  ~ResultPanel() { }
};

/* left top of ButtonPanel */
class SliderPanel : public Fl_Box
{
public:
  SliderPanel (int X, int Y, int W, int H,char const *LL=0)
    : Fl_Box(X,Y,W,H,LL) {

    int atx=X+8;
    int aty=Y+6;
    int dx=(W-6)/5;
    int sh=H-26;

    color(FL_GRAY);
    align(FL_ALIGN_TOP);
    box(FL_UP_BOX);

#define SLIDERCOLOR fl_gray_ramp(8)
    // WARNING: vertical sliders have maximum=bottom, minimum=top
    sliders.T=new Fl_Fill_Slider(atx, aty, 20, sh,"T");
    atx+=dx;
    sliders.T->type(FL_VERT_FILL_SLIDER);
    sliders.T->minimum(log(MAXT));
    sliders.T->maximum(log(MINT));
    sliders.T->value(0);
    sliders.T->color(SLIDERCOLOR);
    // see MIN    
    sliders.T->tooltip("\
THERMOSTAT TEMPERATURE T\n\
\n\
T ∈ [0.1, 5] (logarithmic scale);\n\
outside this range, use cmd: T=<value>");

    sliders.tau=new Fl_Fill_Slider(atx, aty, 20, sh,"τ");
    sliders.tau->type(FL_VERT_FILL_SLIDER);
    sliders.tau->minimum(log(10));
    sliders.tau->maximum(log(0.1));
    sliders.tau->color(SLIDERCOLOR);
    sliders.tau->value(0);
    sliders.tau->tooltip("\
THERMOSTAT COUPLING TIME τ\n\
\n\
τ ∈ [0.1, 10].\n\
Long τ means that the system is well insulated\n\
from the thermostat and the system temperature\n\
approaches the predefined value T slowly,\n\
and vice versa.\n\
Outside the slider range, use cmd: tau=<value>");

    // the same position as tau
    sliders.d=new Fl_Fill_Slider(atx, aty, 20, sh,"d");
    atx+=dx;
    sliders.d->type(FL_VERT_FILL_SLIDER);
    sliders.d->minimum(0);
    sliders.d->maximum(-1);
    sliders.d->value(-0.5);
    sliders.d->color(SLIDERCOLOR);
    sliders.d->tooltip("\
MONTE CARLO TRIAL DISPLACEMENT d\n\
\n\
The range depends on N, logarithmic scale.\n\
See also button [set d] for automatic adjustment.\n\
Applies also to volume-change in NPT.\n\
Maximum displacement = L/2.");

    sliders.g=new Fl_Fill_Slider(atx, aty, 20, sh,"g");
    atx+=dx;
    sliders.g->type(FL_VERT_FILL_SLIDER);
    sliders.g->minimum(MAXG);
    sliders.g->maximum(MING);
    sliders.g->value(0);
    sliders.g->color(SLIDERCOLOR);
    sliders.g->tooltip("GRAVITY g");

    sliders.rho=new Fl_Fill_Slider(atx, aty, 20, sh,"ρ");
    sliders.rho->type(FL_VERT_FILL_SLIDER);
    sliders.rho->minimum(log(MAXRHO));
    sliders.rho->maximum(log(MINRHO));
    sliders.rho->value(log(L2rho(L)));
    sliders.rho->color(SLIDERCOLOR);
    sliders.rho->tooltip("\
NUMBER DENSITY ρ=N/Area\n\
\n\
N is kept constant, therefore increasing ρ shrinks\n\
the box size L, and vice versa.\n\
The simulation box is displayed in the same size\n\
so that the particle diameters are rescaled instead.");

    sliders.P=new Fl_Fill_Slider(atx, aty, 20, sh,"P");
    atx+=dx;
    sliders.P->type(FL_VERT_FILL_SLIDER);
    sliders.P->minimum(MAXP);
    sliders.P->maximum(MINP);
    sliders.P->value(P);
    sliders.P->color(SLIDERCOLOR);
    sliders.P->tooltip("\
PRESSURE P\n\
\n\
for Method = MD/NPT (Berendsen) only.\n\
The barostat time constant is derived from\n\
τ for thermostat and estimated compressibility.\n\
N is kept constant, therefore increasing density\n\
shrinks the box size L, and vice versa.\n\
The simulation box is displayed in the same size\n\
so that the particle diameters are rescaled instead.");

    sliders.N=new Fl_Fill_Slider(atx, aty, 20, sh,"N");
    atx+=dx;
    sliders.N->type(FL_VERT_FILL_SLIDER);
    //    sliders.N->minimum(log(MAXN));
    //    sliders.N->maximum(0);
    sliders.N->minimum(pow(MAXN,1./NPOW));
    sliders.N->maximum(1);
    sliders.N->value(0);
    sliders.N->color(SLIDERCOLOR);
    sliders.N->tooltip("\
NUMBER OF PARTICLES N\n\n\
• decrease N: randomly selected particles are removed\n\
• increase N: particles are inserted at random places\n\
\n\
The box is rescaled to keep constant density.\n\
MD will likely fail and will be temporarily replaced by MC.");
  }
};

/* right top of ButtonPanel */
class WallPanel : public Fl_Box
{
public:
  WallPanel (int X, int Y, int W, int H,char const *LL=0)
  : Fl_Box(X,Y,W,H,LL) {
    color(FL_GRAY);
    align(FL_ALIGN_TOP);
    box(FL_UP_BOX);

    int atx=X+6;
    int aty=Y+6;
    const char *winfo="toggle attractive (green)\nand repulsive (red) walls";

    buttons.wally =new Fl_Light_Button(atx+W/2-39,aty,70,20,"   top");
    buttons.wally->selection_color(FL_GREEN);
    //    buttons.wally->align(FL_ALIGN_CENTER); does not center :(
    buttons.wally->tooltip(winfo);

    buttons.wallx =new Fl_Light_Button(atx,aty+27,53,20," left");
    buttons.wallx->selection_color(FL_GREEN);
    // no function to select the unselected state color ?!
    buttons.wallx->tooltip(winfo);

    buttons.wallxL=new Fl_Light_Button(atx+W-66,aty+27,53,20,"right");
    buttons.wallxL->selection_color(FL_GREEN);
    buttons.wallxL->tooltip(winfo);

    buttons.wallyL=new Fl_Light_Button(atx+W/2-39,aty+54,70,20,"bottom");
    buttons.wallyL->selection_color(FL_GREEN);
    buttons.wallyL->tooltip(winfo);
  }
};

/* below WallPanel */
class ExpertPanel : public Fl_Box
{
public:
  ExpertPanel (int X, int Y, int W, int H,char const *LL=0)
  : Fl_Box(X,Y,W,H,LL) {
    color(FL_GRAY);
    align(FL_ALIGN_TOP);
    box(FL_UP_BOX);

    int atx=X+6;
    int aty=Y+6;

    buttons.meas=new Fl_Light_Button(atx,aty,70,20,"record");
    buttons.meas->selection_color(FL_GREEN);
    // no function to select the unselected state color ?!
    buttons.meas->tooltip("\
RECORDING\n\
\n\
Start/stop recording of measured data.\n\
When recording stops, averages with estimated\n\
standard errors and optional graphs/convergence\n\
profiles are printed to a file.\n\
Use menu [Show] to select quantites/graphs\n\
Use menu [File] to change file name");

    buttons.comma=new Fl_Light_Button(atx+W-43,aty,30,20," ,");
    buttons.comma->selection_color(FL_GREEN);
    buttons.comma->tooltip("\
DECIMAL POINT=COMMA\n\
\n\
If selected, comma (,) will be used\n\
in decimal numbers instead of period (.)\n\
in the recorded file.");

    incl=new Fl_Choice(atx+59,aty+28,W-70,25,"include:");
    incl->menu(inclitems);
    incl->tooltip("\
INCLUDE TO THE RECORDED FILE\n\
\n\
Nothing = only statistics of basic quantities recorded\n\
Conv.prof. = + convergence profiles of basic quantities\n\
Dens.prof. = + the active density/distribution function\n\
     as selected in menu Show and shown in graph:\n\
     RDF, y-profile, droplet and cavity radial profiles\n\
Both = + both Conv.prof. and Dens.prof.");

    buttons.parms=new Fl_Input(atx+38,aty+60,W-49,20,"cmd:");
    buttons.parms->callback(cb_parms);
    buttons.parms->tooltip("\
PARAMETERS\n\
\n\
Enter parameter as command:\n\
  VARIABLE=number\n\
Both dec. comma (,) and period (.) accepted\n\
Available VARIABLEs are:\n\
  N = number of atoms\n\
  rho = number density\n\
  L = box size\n\
  T = temperature (not NVE)\n\
  tau = thermostat time constant (MD)\n\
  d = trial displacement size (MC), in L/2\n\
  dV = relative trial volume change (NPT MC)\n\
  dt = h = MD timestep; 0 = based on T,τ\n\
  g = gravity\n\
  P = pressure (NPT only)\n\
  qtau = for MD/NPT: tauP=qtau*tau\n\
  stride = every stride-th config. shown\n\
  block = measurement block\n\
  wall = wall number density\n\
Example:\n\
  rho=0.01\n\
The values may be out of slider range.");

  }
};

// lower half of the right panel
class ButtonPanel : public Fl_Box
{
public:
  ButtonPanel (int X, int Y, int W, int H)
    : Fl_Box(X,Y,W,H) {

    /*
      SliderPanel *sliderpanel=new SliderPanel(X+3,Y+20,W/2,H*3/5-25,"Parameters");
      WallPanel *wallpanel=new WallPanel(X+W/2+9,Y+20,W/2-12,91,"Walls");
      ExpertPanel *expertpanel=new ExpertPanel(X+W/2+9,Y+H*3/5-71,W/2-12,66,"Expert");
    */

    new SliderPanel(X+3,Y+20,W/2,H*3/5-20,"Parameters");
    new WallPanel(X+W/2+9,Y+20,W/2-12,H/4-7,"Walls");
    //    new ExpertPanel(X+W/2+9,Y+H*3/5-71,W/2-12,66,"Expert");
    new ExpertPanel(X+W/2+9,Y+H/4+37,W/2-12,H/4,"Expert");
    
#define CHOICEW 102

    // bottom four lines
    
    int atx=X+CHOICEW+7,aty=Y+H*2/3-10;

    molsize=new Fl_Choice(atx,aty,CHOICEW,25,"molecule size:");
    molsize->menu(molsizeitems);
    molsize->tooltip("\
MOLECULE SIZE\n\
\n\
Real = real molecule size (collision diameter) is shown\n\
Small = half the real size\n\
Dot = 5 pixels\n\
Pixel = 1 pixel");

    resetchoices=new Fl_Button(atx+CHOICEW+14,aty,88,25,"reset view");
    resetchoices->callback(cb_resetchoices);
    resetchoices->tooltip("\
RESET VIEW\n\
\n\
• Reset molecule size, draw mode,\n\
  and color to the defaults\n\
• Clear red HINTs\n\
• Reset convergence profile\n");

    setd=new Fl_Light_Button(atx+CHOICEW+4,aty+30,107,25,"set MC move");
    setd->selection_color(FL_GREEN);
    setd->tooltip("\
AUTOMATIC DETERMINATION OF MC MOVE SIZE\n\
\n\
checked (green) = the MC move size d is automatically\n\
   determined to reach acceptance ratio about 0.3\n\
   • microreversibility is slightly violated\n\
unchecked = the last value of d applies\n\
   • exactly microreversible\n\
In NPT MC, the same rules apply to volume change dV.\n\
Move length d can be also set by slider \"d\" or by\n\
\"d=number\" (in units of L/2) from input field cmd:.\n\
The value of dV (max relative change of L in a step)\n\
can be set from cmd:, but not from a slider.");
    setd->value(1);

    run=new Fl_Light_Button(atx+CHOICEW+32,aty+60,48,25,"run");
    run->selection_color(FL_GREEN);
    run->tooltip("RUN CONTROL\n\n\
checked (green) = simulation is running\n\
unchecked = simulation is temporarily stopped\n\
Useful to chage parameters if the current setup\n\
keeps launching errors.");
    run->value(1);

    aty+=30;
    drawmode=new Fl_Choice(atx,aty,CHOICEW,25,"draw mode:");
    drawmode->menu(drawmodeitems);
    drawmode->tooltip("\
DRAW MODE\n\
\n\
Movie = show moving molecules\n\
Traces = draw trajectories gradually fading out\n\
Lines = draw trajectories (best with Dot or Pixel)\n\
Nothing = do not draw (use to speed up)");

    aty+=30;
    colormode=new Fl_Choice(atx,aty,CHOICEW,25,"color mode:");
    colormode->menu(colormodeitems);
    colormode->tooltip("\
COLOR MODE\n\
\n\
Black = all molecules are black\n\
y-split = top half molecules are colored blue,\n\
    bottom half red, and these colors are kept\n\
Neighbors = color atoms by the number of neighbors\n\
    (color code shown top right)\n\
Random = colorize atoms randomly (and keep)\n\
Keep = keep the colors (e.g., after Neighbors)");

    atx=X+4;
    aty=Y+H-40;

    sliders.speed=new Fl_Fill_Slider(atx, aty, PANELW/2-8, 20,"simulation speed");
    sliders.speed->type(FL_HOR_FILL_SLIDER);
    sliders.speed->minimum(-0.5); // changed from -1 (too slow)
    sliders.speed->maximum(2.5); // nit=15
    sliders.speed->value(initspeed);
    sliders.speed->color(SLIDERCOLOR);
    sliders.speed->tooltip("\
SIMULATION SPEED\n\
\n\
• left: slower - delays are added\n\
• right: faster - not all configurations are\n\
   displayed (default=every 3th, max 15th)\n");

    sliders.block=new Fl_Fill_Slider(atx+PANELW/2-1, aty,PANELW/2-8, 20,"measurement block");
    sliders.block->type(FL_HOR_FILL_SLIDER);
    sliders.block->minimum(0); // log(1)
    sliders.block->maximum(2); // log(100)
    sliders.block->value(1); // log (10)
    sliders.block->color(SLIDERCOLOR);
    sliders.block->tooltip("MEASUREMENT BLOCKING\n\
\n\
block ∈ [1, 100] (log-scale)\n\
Measurements (of kinetic temperature Tkin,\n\
pressure P, RDF, etc.) are averaged in blocks\n\
and shown after a block is completed.");

  }
  ~ButtonPanel() { }
};


//JK Panel is Fl_Box, because:
// 1) system redraws the box
// 2) win.resizable does not work
class Panel : public Fl_Box {
public:
  InfoPanel *infopanel;
  ResultPanel *resultpanel;
  ButtonPanel *buttonpanel;

  Panel(int X,int Y,int W,int H) : Fl_Box(X,Y,W,H) {
    color(FL_GRAY);
    box(FL_UP_BOX);
    infopanel = new InfoPanel(X,Y,W,INFOPANELH); // constants, parameters
    resultpanel = new ResultPanel(X,Y+INFOPANELH,W,RESULTPANELH);
    buttonpanel = new ButtonPanel(X,Y+(INFOPANELH+RESULTPANELH),W,BUTTONPANELH);
  };
  ~Panel() { };
};

class Panel *panel;
