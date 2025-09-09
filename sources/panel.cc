/* the right panel is drawn here */

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
    if ( (c=strstr(s,"=0.")) ) memcpy(c+1,c+2,strlen(c));
  }
} // shorten()

/* top: parameters and run info */
class InfoPanel : public Fl_Box ///////////////////////////////////// InfoPanel
{
  //private:
protected:
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
      int N,stride,block;
    } last;

    static struct lastc_s {
      int N,Td,d,Tt,L,g;
    } lastc;

    Fl_Box::draw();

    fl_font(FL_HELVETICA,TSIZE);
    if (method==CREUTZ||method==NVE) fl_color(0,127,0);
    else if (isNPT(method)) fl_color(150,0,160);
    else fl_color(0,100,150);

    /* basic info */
    if (mauto) {
      if (method==METROPOLIS)
        strcpy(info,methodinfo[AUTO]);
      else
        sprintf(info,"[%s]",methodinfo[method]); }
    else
      strcpy(info,methodinfo[method]);

    // local short error message with timer
    static const char *shortmsg;

    if (lasterrmessage || shortmsg) {
      static double t0=0,t1;

      if (!shortmsg) t0=mytime();

      switch (lasterrmessage) { // max 31 chars!
        case MDFAILED: shortmsg="MD failed: switching to MC"; break;
        case NPTFAILED: shortmsg="NPT failed: switching to MC"; break;
        case MSDFAILED: shortmsg="MSD undefined: can't follow"; break;
        default:; }

      if (shortmsg) {
        strcpy(info,shortmsg);
        fl_color(OI_RED); }

      t1=mytime(); // show short message for 4 s
      if ((t1-t0)>4) {
        lasterrmessage=0;
        shortmsg=NULL; } }

    fl_draw(info,atx+PANELW-11+3*(info[0]=='[')-fl_width(info),aty);
    fl_color(FL_BLACK);

    // the variable that just have changed is RED
    if (N!=last.N) lastc.N=10;
    if (lastc.N) fl_color(OI_RED),lastc.N--;
    last.N=N;
    sprintf(s,"N=%d %s",N,FFinfo(0));
    fl_draw(s,atx,aty);
    fl_color(FL_BLACK);
    fl_font(FL_HELVETICA,16);

    aty+=TSIZE+4;

    switch (method) {
      case CREUTZ:
        if (d!=last.d) lastc.d=10;
        if (lastc.d) fl_color(OI_RED),lastc.d--;
        sprintf(s,"E=%.3f  d*Lh=%5.3f",Econserved,d*Lh);
        break;
      case METROPOLIS:
        if (T!=last.T || d!=last.d) lastc.Td=10;
        if (lastc.Td) fl_color(OI_RED),lastc.Td--;
        sprintf(s,"T=%6.4f d*Lh=%5.3f",T,d*L);
        break;
      case MCNPT:
        if (T!=last.T || d!=last.d || P!=last.P) lastc.Td=10;
        if (lastc.Td) fl_color(OI_RED),lastc.Td--;
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
        if (lastc.Tt) fl_color(OI_RED),lastc.Tt--;
        sprintf(s,"T=%6.4f  τ=%5.3f",T,tau); }

    last.T=T;
    last.P=P;
    last.tau=tau;
    last.d=d;
    fl_draw(s,atx,aty); aty+=20;
    fl_color(FL_BLACK);

    if (method==MCNPT || method==MDNPT || method==MTK) {
      fl_color(OI_RED);
      if (bc==PERIODIC) sprintf(s,"L=%6.3f",L);
      else sprintf(s,"L=%6.3f wall=%.3g",L,walldens/PI); }
    else {
      if (L!=last.L) lastc.L=10;
      if (lastc.L) fl_color(OI_RED),lastc.L--;
      last.L=L;
      if (bc==PERIODIC) sprintf(s,"L=%-.4g ρ=%.4g",L,L2rho(L));
      else sprintf(s,"L=%-.4g ρ=%.4g wall=%.3g",L,L2rho(L),walldens/PI); }

    fl_draw(s,atx,aty); aty+=20;
    fl_color(FL_BLACK);

    { // this is unnecessarily shown too often (every step)
      static double t0,tav;
      static int n=0,ncell; // 1st value will be wrong
      double t=mytime();

      if (gravity!=last.gravity || speed.stride!=last.stride || block!=last.block) lastc.g=10;
      if (lastc.g) fl_color(OI_RED),lastc.g--;
      last.gravity=gravity;
      last.stride=speed.stride;
      last.block=block;

      /* averages cycle time over at least 2 s */
      sprintf(s,tav>=10?
              "g=%5.3f  stride*block=%d*%d  %.3f%+.4g ms  ncell=%d":
              "g=%5.3f  stride*block=%d*%d  %.3f%+.3f ms  ncell=%d",
              gravity,
              speed.stride,block,
              speed.timerdelay*1000,tav,ncell);

      if (t-t0>2) { // update in 2 s or a bit more
        // this will be shown later
        tav=((t-t0)/n-speed.timerdelay)*1000; ncell=ss.ncell;
        n=0;
        t0=t; }
      n++;

      if (h()>84)
        fl_draw(s,atx,aty+(h()-92)/2);
    } // measure wall time and ncell autoset

    fl_color(FL_BLACK);
  } // draw()

public:
  InfoPanel (int X, int Y, int W, int H) : Fl_Box(X,Y,W,H) {
    color(FL_WHITE);
    box(FL_UP_BOX);
  }
  ~InfoPanel() { }
};

/* results: numbers and graph */
class ResultPanel : public Fl_Box ///////////////////////////////// ResultPanel
{
private:
  void draw () {
    /* (blocked) results print/draw */
    if (iblock>=block) {
      char s[64],s2[64];
      const char *yaxis;
      int atx=x()+6,aty=y()+18;
      int i;
      double CPval=0;
      // de-optimized to make the code a bit cleaner
      double Vav=sum.V/iblock;
      double momentumav=sum.momentum/iblock;
      double rhoav=N/Vav;
      double Tav=sum.Tk/iblock; // will be replaced by T for Pvir
      vector Pav={sum.P.x/iblock,sum.P.y/iblock};
      double Pvir=(Pav.x+Pav.y)/2;
      vector MSDav=zero; // calc. here = last block
      double Hav=sum.H/iblock,Z;
      int n=-1;

      if (isMC(method)) {
        Tav=T; // to be used in calculations of Pvir etc.
        if (method==CREUTZ) Tav=sum.Ekin/iblock; }

      Fl_Box::draw();

      if (files.record) {
        // the output order in CP is reversed
        // applies to SIMANME-CP1.csv or SIMANME.txt CP1

        // conserved energy
        if (method==NOSE_HOOVER || method==MTK)
          addM("Econserved",(sum.Econserved+sum.U)/iblock);

        // enthalpy
        if (isNPT(method)) addM("H",Hav);

        if (bc<PERIODIC) {
          addM("P(top_wall)",sum.fwy/iblock);
          addM("P(bottom_wall)",sum.fwyL/iblock); }
        if (bc==BOX) {
          addM("P(left_wall)",sum.fwx/iblock);
          addM("P(right_wall)",sum.fwxL/iblock); }
        // ? errmessage=MSDFAILED;
        if (!MSDskip) MSD(&MSDav);
        if (method==VICSEK) addM("momentum",momentumav);
        addM("MSDy",MSDav.y);
        addM("MSDx",MSDav.x);
        addM("γ",(sum.P.y-sum.P.x)/iblock*L);
        addM("Pyy",Pav.y);
        addM("Pxx",Pav.x);
        addM("Pvir",Pvir);
        addM("Z",(sum.P.x+sum.P.y)/(2*iblock*T*L2rho(L)));
        addM("V",Vav);
        n=addM("Epot",sum.U/iblock); }

      fl_font(FL_HELVETICA,16);

      fl_color(OI_BLUE);

      const char *eq=dtfixed?"=":"⇒"; // ⇒ not shown?

      if isMD(method) {
        if (method==MDNPT || method==MTK)
          sprintf(s,"Tkin=%5.3f  dt%s%5.4f  ρ=N/⟨V⟩=%.4f",Tav,eq,dt,rhoav);
        else
          sprintf(s,"Tkin=%5.3f  dt%s%5.4f",Tav,eq,dt);
        if (files.record) addM("Tkin",Tav); }
      else {
        if (method==MCNPT)
          sprintf(s,"Tbag=%5.3f  acc.r.=%.4f (V %.4f)  ρ=N/⟨V⟩=%.4f",
                  sum.Tk/iblock,sum.accr/iblock,sum.Vaccr/iblock,rhoav);
        else /* Metropolis, CREUTZ */
          sprintf(s,"Tbag=%5.3f  acc.r.=%5.3f",
                  sum.Tk/iblock,sum.accr/iblock);
        if (files.record) addM("Tbag",sum.Ekin/iblock); }

      shorten(s);
      fl_draw(s,atx,aty);

      if (files.record) {
        sprintf(s,"n=%d",n);
        fl_draw(s,atx+PANELW-11-fl_width(s),aty); }

      aty+=20;

      // MC: equivalent kinetic energy to print internal energy
      double eqEkin=T*En.Nf*iblock/2.;

      /* convergence profile: interpretation of menu items */
      yaxis=NULL;
      switch (Show) {

        case ENERGY:
          if (isNPT(method)) {
            yaxis="enthalpy";
            CPval=(sum.U+eqEkin)/iblock+Vav*P; }
          else if (method==NVE || method==CREUTZ) {
            yaxis="Epot";
            CPval=sum.U/iblock; }
          else if (isMC(method)) {
            yaxis="Epot+f/2*kT";
            CPval=(sum.U+eqEkin)/iblock; }
          else {
            yaxis="Epot+Ekin";
            CPval=(sum.U+sum.Ekin)/iblock; }
          break;

        case TEMPERATURE:
          if (isMC(method)) {
            yaxis="Tbag";
            CPval=sum.Ekin/iblock; /* this is the bag */ }
          else {
            yaxis="Tkin";
            CPval=sum.Ekin/En.Nf*2/iblock; }
          break;

        case PRESSURE:
          yaxis="Pvir";
          CPval=Pvir;
          break;

        case VOLUME:
          if (isNPT(method)) {
            yaxis="volume";
            CPval=Vav; }
          break;

        case MOMENTUM:
          if (isMD(method)) {
            yaxis="momentum";
            CPval=momentumav; }
          break;

        case INTMOTION:
          switch (method) {
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

      if (files.record) {
        addM("Etot",(sum.U+sum.Ekin)/iblock);
        addM(NULL,0); // end of line - must be the last of all addM's
      }

      // second line of cyan-blue info
      if (isNPT(method)) Z=P/(rhoav*T);
      else if (isNVE(method)) Z=Pvir/(rhoav*Tav);
      else Z=Pvir/(rhoav*T);
      sprintf(s2,"V=%.4g  Pvir=%.4g  Z=%.4g",Vav,Pvir,Z);
      shorten(s2);

      if (Show==YPROFILE) {
        sprintf(s,"Pxx=%.3g  Pyy=%.3g  γ=%.3f",
                Pav.x,Pav.y,(Pav.y-Pav.x)/2*L); // repeated below - why?
        shorten(s);
        fl_draw(s,atx,aty); }
      else if (Show==QUANTITIES) {
        fl_draw(s2,atx,aty);
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
                Pav.x,Pav.y,(Pav.y-Pav.x)/2*L); // this is doubled - why?
        shorten(s);
        fl_draw(s,atx,aty); }
      else
        fl_draw(s2,atx,aty);

      aty+=20;

      if (fabs(momentumav)>1e-9) sprintf(s,"Etot=%g  Epot=%g  p=%.4f",(sum.U+sum.Ekin)/iblock,sum.U/iblock,momentumav);
      else sprintf(s,"Etot=%g  Epot=%g",(sum.U+sum.Ekin)/iblock,sum.U/iblock);
      fl_draw(s,atx,aty);
      aty+=20;

      if (shownbrs) {
        fl_color(FL_BLACK);
        aty+=4;
        fl_draw("neighbors:",atx,aty);
        atx+=86;
        int dx=28;

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
      else if (isNPT(method)) {
        if (method==MCNPT) sprintf(s,"H=%g acc.r.(V)=%.2f",Hav,sum.Vaccr/iblock);
        else sprintf(s,"H=%g",Hav);
        shorten(s);
        fl_draw(s,atx,aty); }

      aty=y()+106;
      atx=x()+41;

      fl_color(OI_BLUE);

      if (Show!=lastShow) {
        // reset graph
        justreset=1;
        itot=HISTMAX+1;
        lastShow=Show; }

      if (h()>180) switch (Show) {

        case ENERGY:
        case TEMPERATURE:
        case PRESSURE:
        case VOLUME:
        case MOMENTUM:
        case INTMOTION:
          if (!yaxis) {
            fl_color(OI_RED);
            fl_draw(Showinfo[Show],atx-20,aty+20);
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

            if (justreset) lmin=emin,lmax=emax,justreset=0;

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
            fl_draw(90,yaxis,atx-18,aty+53+fl_width(yaxis)/2);

            if (emax<=emin) {
              fl_color(OI_RED);
              fl_draw(Showinfo[Show],atx,aty+40);
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
            // fl_color(OI_YELLOW); // too dark
            fl_color(LIGHTYELLOW);

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

          // print r=0,1,...
          loopto (i,0,8) {
            sprintf(s,"%d",i);
            fl_draw(s,atx-5+50*i,aty+115); }
          if (shownbrs) fl_draw("neighbor limit",atx+72,aty+130);
          else fl_draw("rmin",atx+54,aty+130);
          fl_draw("r",atx+220,aty+125);
          fl_draw("0",atx-13,aty+107);
          fl_draw("max",atx-37,aty+9);

          vv=Sqr(50*L/N)*2/PI/iblock;
          if (files.record) loop (i,0,HISTMAX) rhosum[i]+=hist[i]*vv/(2*i+1);

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
          fl_draw("RDF",atx+5,aty+20);

          fl_color(LIGHTYELLOW);
          loop (i,1,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1);
            //            if (l)  // causes to draw (leave green) maximum, because l==0 at maximum
            fl_line(atx+i,aty+l,atx+i,aty+100);
            /* omit r=0,1,2,3,4,5: */
            if (i%50==48) i+=2; }

          // NB: it is very slow to switch colors - loops unrolled

          // vertical lines finished
          fl_color(OI_GREEN);
          loop (i,0,HISTMAX) {
            fl_line(atx+i,aty,atx+i,aty+100);
            if (i%50==0) i+=48; }
          if (shownbrs) i=(int)(ss.rnbr*50+0.5); else i=59; /* rmin */
          fl_line(atx+i,aty,atx+i,aty+100);

          fl_color(FL_DARK_YELLOW);
          loop (i,0,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1);
            if (l<100) fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%50==0) i+=48; }

          if (shownbrs) {
            i=(int)(ss.rnbr*50+0.5);
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

          vv=HISTMAX/Sqr(L)/iblock; // NVT (NVE) assumed

          if (files.record) loop (i,0,HISTMAX) rhosum[i]+=hist[i]*vv;

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
          fl_draw("Vertical density profile",atx+5,aty+20);

          fl_color(LIGHTYELLOW);
          loop (i,1,HISTMAX) {
            l=100.5-vv*hist[i];
            if (l<0) { fl_color(255,200,0); fl_line(atx+i,aty,atx+i,aty+100); fl_color(LIGHTYELLOW); }
            else fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%(HISTMAX/2)==HISTMAX/2-2) i+=2; }

          // NB: it is very slow to switch colors - loops unrolled
          fl_color(OI_GREEN);
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

          vv=Sqr(HISTMAX/Lh)/PI/iblock; // NVT (NVE) assumed
          if (files.record) loop (i,0,HISTMAX) rhosum[i]+=hist[i]*vv/(2*i+1);

          fl_color(FL_DARK_GREEN);
          fl_line_style(FL_SOLID,0,NULL);
          fl_rectf(atx,aty,HISTMAX,100);

          fl_color(FL_WHITE);
          if (Show==CPROFILE)
            fl_draw("Cavity radial density profile",atx+5,aty+20);
          else
            fl_draw("Droplet radial density profile",atx+195,aty+20);

          vv*=100.0;
          fl_color(LIGHTYELLOW);
          loop (i,1,HISTMAX) {
            l=100.5-vv*hist[i]/(2*i+1);
            if (l<0) { fl_color(255,200,0); fl_line(atx+i,aty,atx+i,aty+100); fl_color(LIGHTYELLOW); }
            else fl_line(atx+i,aty+l,atx+i,aty+100);
            if (i%(HISTMAX/2)==HISTMAX/2-2) i+=2; }

          // NB: it is very slow to switch colors - loops unrolled
          fl_color(OI_GREEN);
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
        } // show


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
class ParameterPanel : public Fl_Box /////////////////////////// ParameterPanel
{
public:
  ParameterPanel (int X,int Y,int W,int H,char const *LL=0) : Fl_Box(X,Y,W,H,LL) {

    int atx=X+6;
    int aty=Y+6;
    int dx=(W-35)/5; // 5 vertical sliders and N+, N-, g=0
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
T ∈ [0.1, 5]\n\
Outside this range, use cmd: T=<value>.");

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
from the thermostat. The system temperature\n\
approaches the predefined value T slowly.\n\
Short τ means poor insulation and fast change of T.\n\
Outside the slider range, use cmd: tau=<value>.");

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
The range depends on N.\n\
See also button [set d] for automatic adjustment.\n\
Applies also to volume-change in NPT.\n\
Maximum displacement = L/2.");

    /* Vicsek only (implies periodic b.c.) */
    sliders.c=new Fl_Fill_Slider(atx, aty, 20, sh,"c");
    sliders.c->type(FL_VERT_FILL_SLIDER);
    sliders.c->minimum(6);
    sliders.c->maximum(2);
    sliders.c->value(0);
    sliders.c->color(SLIDERCOLOR);
    sliders.c->tooltip("\
CUTOFF c\n\
\n\
Defines the range for determining neighbors\n\
in the Vicsek model of active matter as well\n\
as the underlying penetrable disk radius.\n\
Equivalent to variable c in [cmd:]; however,\n\
with a narrower allowed range.\n\
This slider is not active for other models\n\
to prevent accidental change.");

    /* not periodic b.c. */
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
    sliders.N->type(FL_VERT_FILL_SLIDER);
    //    sliders.N->minimum(log(MAXN));
    //    sliders.N->maximum(0);
    sliders.N->minimum(NtoSlider(1000)); // slider max N=1000, command max = MAXN
    sliders.N->maximum(NtoSlider(2)); // slider min N=2, command min N=1 from cmd:
    sliders.N->value(0);
    sliders.N->color(SLIDERCOLOR);
    sliders.N->tooltip("\
NUMBER OF PARTICLES N\n\n\
• decrease N: randomly selected particles are removed\n\
• increase N: particles are inserted at random places\n\
\n\
The box is rescaled to keep constant density.\n\
MD may fail and will be temporarily replaced by MC.");

    atx=X+W-42; // from the "next slider"

    Parms.plus=new Fl_Button(atx,aty,35,20,"N+");
    Parms.plus->callback(cb_plus);
    Parms.plus->tooltip("\
INSERT PARTICLE\n\
\n\
One particle is inserted: N:=N+1.\n\
Useful if the N-slider does not have\n\
sufficient resolution.");

    Parms.gzero=new Fl_Button(atx,aty+sh/2-8,35,20,"g=0");
    Parms.gzero->callback(cb_gzero);
    Parms.gzero->tooltip("\
ZERO GRAVITY\n\
\n\
Set the gravity to zero, g:=0.");

    Parms.minus=new Fl_Button(atx,aty+sh-21,35,20,"N-"); // "N−" not shown
    Parms.minus->callback(cb_minus);
    Parms.minus->tooltip("\
REMOVE PARTICLE\n\
\n\
One particle is removed: N:=N-1.\n\
Useful if the N-slider does not have\n\
sufficient resolution.");
    // N:=N−1 not shown
  }
};

/* right top of ButtonPanel */
class WallPanel : public Fl_Box ///////////////////////////////////// WallPanel
{
public:
  WallPanel (int X,int Y, int W,int H,char const *LL=0) : Fl_Box(X,Y,W,H,LL) {
    color(FL_GRAY);
    align(FL_ALIGN_TOP);
    box(FL_UP_BOX);

    int atx=X+6;
    int aty=Y+6;
    const char *winfo="toggle attractive (green)\nand repulsive (red) walls";
    const char *shiftinfo="shift the configuration\nperiodically";

    buttons.wally =new Fl_Light_Button(atx+W/2-39,aty,70,20,"   top");
    buttons.wally->selection_color(OI_GREEN);
    //    buttons.wally->align(FL_ALIGN_CENTER); does not center :(
    buttons.wally->tooltip(winfo);
    buttons.shift_up=new Fl_Button(atx+W/2-43,aty,75,20,"shift up");
    buttons.shift_up->callback(cb_shift_up);
    buttons.shift_up->tooltip(shiftinfo);

    buttons.wallx =new Fl_Light_Button(atx,aty+26,53,20," left");
    buttons.wallx->selection_color(OI_GREEN);
    // no function to select the unselected state color ?!
    buttons.wallx->tooltip(winfo);
    buttons.shift_left=new Fl_Button(atx,aty+26,75,20,"shift left");
    buttons.shift_left->callback(cb_shift_left);
    buttons.shift_left->tooltip(shiftinfo);

    buttons.wallxL=new Fl_Light_Button(atx+W-66,aty+26,54,20,"right");
    buttons.wallxL->selection_color(OI_GREEN);
    buttons.wallxL->tooltip(winfo);
    buttons.shift_right=new Fl_Button(atx+W-86,aty+26,75,20,"shift right");
    buttons.shift_right->callback(cb_shift_right);
    buttons.shift_right->tooltip(shiftinfo);

    buttons.wallyL=new Fl_Light_Button(atx+W/2-39,aty+52,70,20,"bottom");
    buttons.wallyL->selection_color(OI_GREEN);
    buttons.wallyL->tooltip(winfo);
    buttons.shift_down=new Fl_Button(atx+W/2-45,aty+52,77,20,"shift down");
    buttons.shift_down->callback(cb_shift_down);
    buttons.shift_down->tooltip(shiftinfo);

    buttons.invertwalls=new Fl_Button(atx+W/2-34,aty+25,61,20,"invert");
    buttons.invertwalls->callback(cb_invertwalls);
    buttons.invertwalls->tooltip("\
INVERT WALLS\n\
\n\
Change attractive (green) walls\n\
to repulsive (red) and vice versa.");
  }
};

#define PRINTHEIGHT 23

class ExpertPanel : public Fl_Box ///////////////////////////////// ExpertPanel
{
protected:
  void draw () {
    /* print previous command of [cmd:] aggregated string */
    int atx=x()+6,aty=y()+6;

    //    Fl_Box::draw();
    // I'm missing function drawing a frame using FL_UP_BOX
    fl_frame("XXAA",x(),y(),w(),h());

    fl_font(FL_HELVETICA,16);
    fl_draw_box(FL_FLAT_BOX,atx,aty,PANELW/2-24,PRINTHEIGHT,FL_WHITE);
    fl_color(FL_BLACK);
    fl_draw(printcmd,atx+1,aty+17);
  }

public:
  ExpertPanel (int X,int Y,int W,int H,char const *LL=0) : Fl_Box(X,Y,W,H,LL) {
    color(FL_GRAY);
    align(FL_ALIGN_TOP);
    box(FL_UP_BOX);

    int atx=X+6;
    int aty=Y+6+PRINTHEIGHT+8;

    fl_font(FL_HELVETICA,16);

    buttons.cmd=new Fl_Input(atx+38,aty,W-49,24,"cmd:");
    buttons.cmd->callback(cb_cmd);
    buttons.cmd->tooltip("\
PARAMETERS\n\
\n\
Enter/print parameter:\n\
  VARIABLE                   (print the value)\n\
  VARIABLE=NUMBER (assign)\n\
Both . and , allowed as dec. separator in NUMBER.\n\
Selected available VARIABLEs are:\n\
  block = measurement block\n\
  d = trial displacement size (MC), in L/2\n\
  dt = h = MD timestep; 0 = based on T,τ\n\
  dV = relative trial volume change (NPT MC)\n\
  g = gravity\n\
  L = box size (ρ recalculated)\n\
  N = number of atoms (L recalculated)\n\
  P = pressure (NPT only)\n\
  qtau, qτ = for MD/NPT: tauP=qtau*tau\n\
  rho, ρ = number density (L recalculated)\n\
  T = temperature (not NVE)\n\
  tau, τ = thermostat time constant (MD)\n\
  wall = wall number density\n\
More VARIABLEs (see the manual):\n\
  a b bc c circle ff fps method nbr\n\
  ncell show speed trace v\n\
Example:\n\
  ρ=0.01\n\
Values out of slider range may be accepted.");

    aty+=30;

    incl=new Fl_Choice(atx+59,aty,W-70,25,"include:");
    incl->menu(inclitems);
    incl->tooltip("\
WHAT TO MEASURE AND EXPORT\n\
\n\
Nothing = only statistics of basic quantities\n\
Convergence profile = also convergence profiles\n\
  of basic quantities\n\
Density profile = also the active distribution function\n\
  as selected in menu Show and shown in the graph:\n\
  RDF, y-profile, droplet or cavity, radial profiles\n\
Both = also both above profiles\n\
• If CSV is selected, there is a separate CSV file for each output\n\
• In no CSV, the output goes to the TXT-file");

    aty+=32;

    buttons.record=new Fl_Light_Button(atx,aty,68,20,"record");
    buttons.record->selection_color(OI_GREEN);
    // no function to select the unselected state color ?!
    buttons.record->tooltip("\
RECORDING\n\
\n\
Start/stop recording of measured data.\n\
When recording stops, averages with estimated\n\
standard errors and optional graphs/convergence\n\
profiles are printed to a file.\n\
Use menu [Show] to select quantites/graphs\n\
Use menu [File] to change file name\n\
NB: [Show]->Minimum is changed to Quantities.");

    buttons.csv=new Fl_Light_Button(atx+W/2-34,aty,50,20,"CSV");
    buttons.csv->selection_color(OI_GREEN);
    buttons.csv->value(1);
    buttons.csv->tooltip("\
OUTPUT FORMAT\n\
\n\
Affects convergence profiles and distribution functions.\n\
• If selected, the outputs are exported\n\
   in the CSV format (.csv).\n\
   If also comma [,] is selected, semicolon [;]\n\
   is used as the separator.\n\
• If deselected, the outputs are in the protocol (.txt),\n\
   separated by tabulators.");

    buttons.comma=new Fl_Light_Button(atx+W-87,aty,73,20,"comma");
    //    buttons.comma=new Fl_Light_Button(atx+W-53,aty,40,20,"  , ");
    buttons.comma->selection_color(OI_GREEN);
    buttons.comma->tooltip("\
DECIMAL POINT=COMMA\n\
\n\
If selected, comma (,) will replace period (.)\n\
in the recorded files. Also, the separator\n\
in CSV files is changed into semicolon (;).\n\
In the cmd: field, both . and , are accepted.");
  }
};


// lower half of the right panel
class ButtonPanel : public Fl_Box ///////////////////////////////// ButtonPanel
{
public:
  ExpertPanel *expertpanel;

  ButtonPanel (int X, int Y, int W, int H) : Fl_Box(X,Y,W,H) {

    // ParameterPanel, WallPanel not redrawn ⇒ pointers not needed
    new ParameterPanel(X+4,Y+20,W/2-16,H*3/5-20,"Parameters");
    new WallPanel(X+W/2+4,Y+20,W/2-12,H/4-9,"Walls + shifts");
    // expertpanel contains previous cmd: to be redrawn
    expertpanel=new ExpertPanel(X+W/2+4,Y+H/4+36,W/2-12,H/4+33,"Expert");

    // bottom four lines (expertpanel covers 1st line right)
    // color(FL_YELLOW);
    // box(FL_UP_BOX);

    int atx=X+4, aty=Y+H*2/3-10;

    // 1st line
    molsize=new Fl_Choice(atx+W/2-130,aty,W/4-5,25,"molecule size:");
    molsize->menu(molsizeitems);
    molsize->tooltip("\
MOLECULE SIZE\n\
\n\
Real = real molecule size (collision diameter) is shown\n\
Small = half the real size\n\
Dot = 5 pixels\n\
Pixel = 1 pixel");

    // 2nd line
    aty+=30;
    drawmode=new Fl_Choice(atx+W/2-130,aty,W/4-5,25,"draw mode:");
    drawmode->menu(drawmodeitems);
    drawmode->tooltip("\
DRAW MODE\n\
\n\
Movie = show moving molecules\n\
Traces = draw trajectories gradually fading out\n\
         change length by variable trace (from cmd:)\n\
Lines = draw trajectories (best with Dot or Pixel)\n\
Nothing = do not draw (use to speed up)");

    resetview=new Fl_Button(atx+W/2+4,aty,88,25,"reset view");
    resetview->callback(cb_resetview);
    resetview->tooltip("\
RESET VIEW\n\
\n\
Reset molecule size, draw mode,\n\
and color to the defaults.\n");

    resetgraph=new Fl_Button(atx+W-114,aty,97,25,"reset graph");
    resetgraph->callback(cb_resetgraph);
    resetgraph->tooltip("\
RESET GRAPH\n\
\n\
Erase the convergence profile data\n\
in the green/yellow graph above\n\
and set the min-max range anew.\n");

    // 3rd line
    aty+=30;
    colormode=new Fl_Choice(atx+W/2-130,aty,W/4-5,25,"color mode:");
    colormode->menu(colormodeitems);
    colormode->tooltip("\
COLOR MODE\n\
\n\
Black = All atoms are black\n\
One = One atom black, the rest orange\n\
y-split = Top half of the configuration is blue,\n\
    bottom half red, and these colors are kept\n\
Neighbors = Colorize atoms by the number of\n\
    neighbors, the color code is shown top right\n\
Random = Colorize atoms randomly (and keep)\n\
Art = Cycle through rainbow colors, best with\n\
    draw mode=Line or Trace.\n\
    Change length by variable trace (from cmd:)\n\
Keep = Keep the colors (e.g., after Neighbors)");

    setd=new Fl_Light_Button(atx+W/2+4,aty,118,25,"set MC moves");
    setd->selection_color(OI_GREEN);
    setd->tooltip("\
AUTOMATIC DETERMINATION OF MC MOVE SIZE\n\
\n\
Should be off for productive runs!\n\
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

    run=new Fl_Light_Button(atx+W-65,aty,48,25,"run");
    run->selection_color(OI_GREEN);
    run->tooltip("RUN CONTROL\n\n\
checked (green) = simulation is running\n\
unchecked = simulation is temporarily stopped\n\
Useful to change parameters if the current setup\n\
keeps launching errors.");
    run->value(1);

    // 3rd line
    atx=X+4;
    aty=Y+H-40;

    //    sliders.speed=new Fl_Fill_Slider(atx, aty, PANELW/2-38, 20,"simulation speed");
    sliders.speed=new Fl_Fill_Slider(atx, aty, W/2-14, 20,"simulation speed");
    sliders.speed->type(FL_HOR_FILL_SLIDER);
    sliders.speed->minimum(MINSPEED);
    sliders.speed->maximum(MAXSPEED+1.5); // note that >MAXSPEED+1⇒ max speed
    sliders.speed->value(SPEEDINIT);
    slider2speed(SPEEDINIT);
    sliders.speed->color(SLIDERCOLOR);
    sliders.speed->tooltip("\
SIMULATION SPEED\n\
\n\
• leftmost: very slow (FPS decreased)\n\
• left: every frame shown with the selected FPS\n\
• center: every 3rd frame shown with the selected FPS\n\
• right: every 10th frame shown with the selected FPS\n\
• rightmost: max speed: as above, no delays, faster disk\n\
   drawing method for L>10; however, disks may overdraw text\n\
→watch stride*block and delay in the 4th line\n");

    //    sliders.block=new Fl_Fill_Slider(atx+W/2-26, aty,W/2-38, 20,"measurement block");
    sliders.block=new Fl_Fill_Slider(atx+PANELW/2-1, aty,PANELW/2-8, 20,"measurement block");
    sliders.block->type(FL_HOR_FILL_SLIDER);
    sliders.block->minimum(0); // log(1)
    sliders.block->maximum(2); // log(100)
    sliders.block->value(1);   // log(10)
    sliders.block->color(SLIDERCOLOR);
    sliders.block->tooltip("MEASUREMENT BLOCKING\n\
\n\
block ∈ [1, 100]\n\
Measurements (of kinetic temperature Tkin,\n\
pressure P, RDF, etc.) are averaged in blocks\n\
and shown/exported after a block has completed."); }

  ~ButtonPanel() { }
};


// Panel is inherited from Fl_Box, because:
// 1) system redraws the box
// 2) undocumented method win.resizable() is out of order or deprecated

class Panel : public Fl_Box ///////////////////////////////////////////// Panel
{
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
