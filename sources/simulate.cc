// The main simulation loop:
// Simulate.draw() calculates MC/MD erc ans shows moleciles
// Also, the timer control is here
class Simulate : public Fl_Box {
  void get_data () {
    int i;

    // get information from sliders; cf. cb_parms()
    slider2speed(sliders.speed->value());

    //    printf("%d %d %g\n",time(NULL),stride,timerdelay);
    block=pow(10,sliders.block->value())+0.5;

    T=exp(sliders.T->value());
    d=exp(DSCALE*sliders.d->value());
    tau=exp(sliders.tau->value());

    if (bc==PERIODIC) gravity=0;
    else gravity=sliders.g->value();

    if (mauto && thermostat<=MCNPT) {
      if (mcstart<=0) {
        // initial MC finished, switching to MD

        if (delayed_y_color) { // hack for DIFFUSION (color half AFTER MC)
          y_split();
          delayed_y_color=0; }

        MDreset(mauto);
        mauto=AUTO; /* = 0 */ }
      else
        mcstart-=speed.stride; }

    //    insdel((int)(exp(sliders.N->value())+0.5));
    insdel((int)(pow(sliders.N->value(),NPOW)+0.5));

    if (isNPT(thermostat)) {
      P=sliders.P->value();
      if (P<-1) P=-1;
      if (P>100) P=100; }
    else {
      double oldL=L;
      L=rho2L(exp(sliders.rho->value()));

      if (L2rho(L)>MAXRHO) {
        rho2L(MAXRHO);
        sliders.rho->value(log(L2rho(L))); }
      if (L2rho(L)<MINRHO) {
        rho2L(MINRHO);
        sliders.rho->value(log(L2rho(L))); }
      Lh=L/2;

      if (fabs(L/oldL-1)>1e-9) {
        if (debug) fprintf(stderr,"L: %g→%g (bc=%d)\n",oldL,L,bc);
        loop (i,0,N) {
          r[i].x*=L/oldL;
          r[i].y*=L/oldL; } } }

    walls=(buttons.wallx->value())
      + 2*(buttons.wallxL->value())
      + 4*(buttons.wally->value())
      + 8*(buttons.wallyL->value());
    setwalls();
  }

  void norm() {
    int i;

    /* out-of-box normalization, error for MD+wall */
    loop (i,0,N) {
      if ((r[i].x<=0) || (r[i].x>=L)) {
        if (isMD(thermostat) && (bc==BOX)) {
          if (debug) fprintf(stderr,"MDerror: r[%d] = [%g %g]\n",i,r[i].x,r[i].y);
          MDerror();
          r[i].x=rnd()*(L-1)+0.5; }
        while (r[i].x<0) r[i].x+=L;
        while (r[i].x>L) r[i].x-=L; }
      if ((r[i].y<=0) || (r[i].y>=L)) {
        if (isMD(thermostat) && (bc!=PERIODIC)) {
          if (debug) fprintf(stderr,"MDerror: r[%d] = [%g %g]\n",i,r[i].x,r[i].y);
          MDerror();
          r[i].y=rnd()*(L-1)+0.5; }
        while (r[i].y<0) r[i].y+=L;
        while (r[i].y>L) r[i].y-=L; } }
  }

  void meas() {
    int i,j;
    unsigned ih;
    vector ri;

    /* calculate virial pressure and RDF
       NB: in MD, Ekin will be known AFTER a MD step */
    Upot=0;
    Pcfg.xx=Pcfg.yy=0; /* P=sum r*f = -virial */
    Pcfg.fwx=Pcfg.fwxL=Pcfg.fwy=Pcfg.fwyL=0;
    addP=measure==RDF?addPRDF:addPonly;
    loop (i,0,N) {
      ri=r[i];
      Upot+=gravity*r[i].y; // potential of gravity
      if (bc<PERIODIC) {
        double fwy=fwally(ri.y),fwyL=fwallyL(L-ri.y);
        Pcfg.yy+=ri.y*fwy+(L-ri.y)*fwyL;
        Pcfg.fwy+=fwy;
        Pcfg.fwyL+=fwyL;
        Upot+=uwally(ri.y)+uwallyL(L-ri.y);
        if (bc==BOX) {
          double fwx=fwallx(ri.x),fwxL=fwallxL(L-ri.x);
          Pcfg.xx+=ri.x*fwx+(L-ri.x)*fwxL;
          Pcfg.fwx+=fwx;
          Pcfg.fwxL+=fwxL;
          Upot+=uwallx(ri.x)+uwallxL(L-ri.x); } }
      switch (bc) {
        case BOX:
          loop (j,i+1,N)
            addP(ri.x-r[j].x,ri.y-r[j].y);
          break;
        case SLIT:
          loop (j,i+1,N)
            addP(ni(ri.x-r[j].x),ri.y-r[j].y);
          break;
        case PERIODIC:
          loop (j,i+1,N)
            addP(ni(ri.x-r[j].x),ni(ri.y-r[j].y));
      } }

    if (measure==YPROFILE) {
      cfgcenter();
      loop (i,0,N) {
        unsigned ir=(shift.y-r[i].y)*(HISTMAX/L);
        hist[ir%HISTMAX]++; } }

    if (measure>=RPROFILE) {
      cfgcenter();
      if (measure==RPROFILE) { shift.x+=Lh; shift.y+=Lh; }
      else if (bc==SLIT) shift.y=Lh;
      else if (bc==BOX) shift.y=shift.x=Lh;
      while (shift.x>Lh) shift.x-=L;
      while (shift.y>Lh) shift.y-=L;
      loop (i,0,N) {
        double rr=0;

        switch (bc) {
          case BOX: rr=Sqr(r[i].x-shift.x)+Sqr(r[i].y-shift.y); break;
          case SLIT: rr=Sqrni(r[i].x-shift.x)+Sqr(r[i].y-shift.y); break;
          case PERIODIC: rr=Sqrni(r[i].x-shift.x)+Sqrni(r[i].y-shift.y); }

        ih=sqrt(rr)*(HISTMAX/Lh);
        if (ih<=HISTMAX) hist[ih]++; } }

    //    if (thermostat==METROPOLIS) Tk=T;
    //    Pcfg.xx=(N*Tk+Pcfg.xx)/Sqr(L);
    //    Pcfg.yy=(N*Tk+Pcfg.yy)/Sqr(L);
    // because of MTK, Ekin used instead of N*Tk (no "kinetic pressure correction")

    if (thermostat==CREUTZ) Ekin=0; // this is not fully correct...
    if (thermostat==METROPOLIS || thermostat==MCNPT) Ekin=degrees_of_freedom()*T/2;
    Pcfg.xx=(Ekin+Pcfg.xx)/Sqr(L);
    Pcfg.yy=(Ekin+Pcfg.yy)/Sqr(L);

    Pcfg.fwx/=L;
    Pcfg.fwxL/=L;
    Pcfg.fwy/=L;
    Pcfg.fwyL/=L;
  }

  void draw() { // in addition to drawing, all calculations are called from here
    int i;
    double V=0; // to suppress warning uninitilized

    get_data();
    if (run->value()) {
      loop (i,0,speed.stride) {

        if (i==speed.stride-1) {
          /* measurements BEFORE the last cycle so that Ekin is in sync 
             with Upot AFTER this loop has finished */
          if (measure) {
            meas();
            sum.Pxx+=Pcfg.xx;
            sum.Pyy+=Pcfg.yy;
            sum.fwx+=Pcfg.fwx;
            sum.fwxL+=Pcfg.fwxL;
            sum.fwy+=Pcfg.fwy;
            sum.fwyL+=Pcfg.fwyL;
            V=N/L2rho(L);
            sum.V+=V; /* some corrections for not PERIODIC */
            sum.U+=Upot;
            if (thermostat<=MCNPT) {
              /* the current bag + Upot is constant for CREUTZ */
              Econserved=Upot+bag;
              sum.Econserved+=Econserved;
              sum.Ekin+=bag; } } }

        if (isMD(thermostat)) {
          MDstep();
          sum.Tk+=Tk/speed.stride; }
        else {
          accepted=Vaccepted=0;
          MCsweep();
          sum.ar+=(double)accepted/(N*speed.stride);
          sum.Var+=(double)Vaccepted/speed.stride; }
      } // speed.stride

      /* last cycle finished, both Upot and Ekin are in sync */
      if (thermostat==NVE) {
        Econserved=Upot+Ekin;
        sum.Econserved+=Econserved; }

      if (isNPT(thermostat))
        sum.H+=Upot+Ekin+P*V;
      if (isMD(thermostat))
        sum.Ekin+=Ekin;
      if (thermostat==NOSE_HOOVER || thermostat==MTK) 
        sum.Econserved+=Econserved+Upot;
      iblock++; }

    else
      // simulation not running due to button [run] in off state
      usleep(speed.extradelay+5000);

    norm();

    if (drawmode->value()<3) {
      // draw background
      int cmode=colormode->value();
      if (cmode==0) { molblack(); cmode=4; colormode->value(cmode); }

      if (drawmode->value()==0) Fl_Box::draw();

      if (cmode==2) {
        shownbrs=1;
        neighbors(); }
      else {
        shownbrs=0;
        if (cmode==1) {
          y_split(); cmode=4; colormode->value(cmode); }
        else if (cmode==3) {
          randomcolors(); cmode=4; colormode->value(cmode); }  }

      fl_color(FL_WHITE);
      if (drawmode->value()==1)
        loop (i,0,Sqr(boxsize)/64) fl_point(rnd()*boxsize+BORDER,rnd()*boxsize+MENUH+BORDER);

      boxsize=w()<h() ? w()-2*BORDER : h()-2*BORDER;
      scale=boxsize/L;

      int rad=2;
      double drad=1; // to suppress warning uninitilized
      switch (molsize->value()) {
        case 0: drad=scale/2.; rad=int(drad+.99); break;
        case 1: drad=scale/4.; rad=int(drad+.99); break;
        case 2: rad=1; break;
        case 3: rad=0; break; }

      // DRAW WALLS
      // WARNING: fl_line_style must be after fl_color
      // (because of Windows implementation caveat)

      // left wall
      if (bc==BOX) {
        if (uwallx==uwalla) fl_color(OI_GREEN);
        else fl_color(OI_RED);
        fl_line_style(FL_SOLID,BORDER,NULL); }
      else {
        fl_color(FL_GRAY);
        fl_line_style(FL_DASH,BORDER,NULL); }
      fl_line(BORDER/2,MENUH+BORDER,BORDER/2,boxsize+MENUH+BORDER);

      // right wall
      if (bc==BOX) {
        if (uwallxL==uwalla) fl_color(OI_GREEN);
        else fl_color(OI_RED); }
      // else GRAY+DASH set above
      fl_line(boxsize+3*BORDER/2,MENUH+BORDER,boxsize+3*BORDER/2,boxsize+MENUH+BORDER);

      // top wall
      if (bc!=PERIODIC) {
        if (uwally==uwalla) fl_color(OI_GREEN);
        else fl_color(OI_RED);
        fl_line_style(FL_SOLID,BORDER,NULL); }
      else {
        fl_color(FL_GRAY);
        fl_line_style(FL_DASH,BORDER,NULL); }
      fl_line(BORDER,MENUH+BORDER/2,boxsize+BORDER,MENUH+BORDER/2);

      // bottom wall
      if (bc!=PERIODIC) {
        if (uwallyL==uwalla) fl_color(OI_GREEN);
        else fl_color(OI_RED);
        fl_line_style(FL_SOLID,BORDER,NULL); }
      // else GRAY+DASH set above
      fl_line(BORDER,boxsize+MENUH+3*BORDER/2,boxsize+BORDER,boxsize+MENUH+3*BORDER/2);

      fl_line_style(FL_SOLID,0,NULL);

      fl_color(FL_BLACK);
      // color switching is sllooowww, so optimize if all black
      cmode=0; // the same color
      loop (i,1,N) if (molcol[i]!=molcol[0]) { cmode=1; break; }
      if (cmode==0) fl_color(molcol[0]);

      // draw molecules
      if (rad==0)
        // pixel
        loop (i,0,N) {
          int xx=r[i].x*scale+BORDER+0.5;
          int yy=r[i].y*scale+BORDER+MENUH+0.5;
          if (cmode) fl_color(molcol[i]);
          fl_point(xx,yy); }
      else if (rad==1)
        // 5 pixels
        loop (i,0,N) {
          int xx=r[i].x*scale+BORDER+0.5;
          int yy=r[i].y*scale+BORDER+MENUH+0.5;
          if (cmode) fl_color(molcol[i]);
          fl_point(xx+1,yy);
          fl_point(xx-1,yy);
          fl_point(xx,yy-1);
          fl_point(xx,yy+1);
          fl_point(xx,yy); }
      else {
        if (circlemethod==0) {
          // use library disk functions, integer positions
          loop (i,0,N) {
            int xx=r[i].x*scale+BORDER+0.5;
            int yy=r[i].y*scale+BORDER+MENUH+0.5;
            if (cmode) fl_color(molcol[i]);
            fl_pie(xx-rad,yy-rad,rad*2,rad*2,0.0,360.0); } }
        else if (circlemethod==1) {
          // library circle, double float positions (more precise, slower)
          loop (i,0,N) {
            double xx=r[i].x*scale+BORDER;
            double yy=r[i].y*scale+BORDER+MENUH;
            if (cmode) fl_color(molcol[i]);
            fl_circle(xx,yy,drad); } }
        else {
          // custom made, periodicity obeyed
          if (bc==BOX)
            loop (i,0,N) {
              double xx=r[i].x*scale+BORDER;
              double yy=r[i].y*scale+BORDER+MENUH;
              double ss;
              int j;
              if (cmode) fl_color(molcol[i]);
              loopto (j,xx-drad+1,xx+drad) {
                ss=sqrt(Sqr(drad)-Sqr(j-xx));
                fl_line(j,yy-ss+1, j,yy+ss); } }
          else if (bc==SLIT)
            loop (i,0,N) {
              double xx=r[i].x*scale+BORDER;
              double yy=r[i].y*scale+BORDER+MENUH;
              double ss;
              int j,jj;

              if (cmode) fl_color(molcol[i]);

              loopto (j,xx-drad+1,xx+drad) {
                jj=j;
                if (jj<BORDER) jj+=boxsize;
                if (jj>BORDER+boxsize) jj-=boxsize;
                ss=sqrt(Sqr(drad)-Sqr(j-xx));
                fl_line(jj,yy-ss+1, jj,yy+ss); } }
          else
            loop (i,0,N) {
              double xx=r[i].x*scale+BORDER;
              double yy=r[i].y*scale+BORDER+MENUH;
              double ss;
              int j,jj;
              if (cmode) fl_color(molcol[i]);

              loopto (j,xx-drad+1,xx+drad) {
                jj=j;
                if (jj<BORDER) jj+=boxsize;
                if (jj>=BORDER+boxsize) jj-=boxsize;
                ss=sqrt(Sqr(drad)-Sqr(j-xx));

                if (yy-ss+1<BORDER+MENUH) {
                  fl_line(jj,yy-ss+1+boxsize, jj,BORDER+MENUH+boxsize-1);
                  fl_line(jj,BORDER+MENUH, jj,yy+ss); }
                else if (yy+ss>=BORDER+MENUH+boxsize) {
                  fl_line(jj,yy-ss+1, jj,BORDER+MENUH+boxsize-1);
                  fl_line(jj,BORDER+MENUH, jj,yy+ss-boxsize); }
                else
                  fl_line(jj,yy-ss+1, jj,yy+ss);
              } }
        } // circlemethod
      } } // draw simulation cell

    // usleep(speed.extradelay);

    fl_line_style(FL_SOLID,0,NULL);
    if (mymessage) {
      int atx=BORDER+6,aty=BORDER+MENUH+35;

      fl_color(OI_YELLOW);
      fl_rectf(BORDER,BORDER+MENUH,640,90);
      fl_font(FL_HELVETICA,32);
      fl_color(OI_RED);
      // see Fl::repeat_timeout() below
      if (mymessage==1) {
        fl_draw("Molecular Dynamics failed!",atx,aty); aty+=40;
        fl_draw("Switching temporarily to Monte Carlo…",atx,aty); }
      else {
        fl_draw("Too low density in NPT ensemble!",atx,aty); aty+=40;
        fl_draw("Switching to Metropolis MC…",atx,aty); }
      fl_font(FL_HELVETICA,16);
      fl_color(FL_BLACK); }
  }

  static void timer_callback(void *userdata) {
    Simulate *sim = (Simulate*)userdata;
    static double tcycle=1/speed.FPS;
    clock_t c0=clock();

    if (isMC(thermostat) && !setd->value())
      sliders.d->show();
    else
      sliders.d->hide();

    if (isMD(thermostat) && thermostat!=NVE)
      sliders.tau->show();
    else
      sliders.tau->hide();

    if (isMC(thermostat))
      setd->show();
    else
      setd->hide();

    if (isNPT(thermostat)) {
       sliders.P->show();
       sliders.rho->hide(); }
    else {
       sliders.P->hide();
       sliders.rho->show(); }

    if (thermostat==CREUTZ || thermostat==NVE)
      sliders.T->hide();
    else
      sliders.T->show();

    switch (bc) {
      case BOX:
        buttons.wallx->show();
        buttons.wallxL->show();
        buttons.wally->show();
        buttons.wallyL->show();
        invertwalls->show();
        sliders.g->show();
        break;
      case SLIT:
        buttons.wallx->hide();
        buttons.wallxL->hide();
        buttons.wally->show();
        buttons.wallyL->show();
        invertwalls->show();
        sliders.g->show();
        break;
      case PERIODIC:
        buttons.wallx->hide();
        buttons.wallxL->hide();
        buttons.wally->hide();
        buttons.wallyL->hide();
        invertwalls->hide();
        sliders.g->value(0);
        sliders.g->hide();
        break;
    }
    if (fabs(sliders.g->value())<0.0017) sliders.g->value(0);

    sim->redraw(); // → sim->draw() = all calculations and molecule drawing
    panel->infopanel->redraw();
    panel->resultpanel->redraw();
    menu->redraw();

    files.record=buttons.record->value();
    if (files.record && measure==NONE) measure=QUANTITIES;
    //    fprintf(stderr,"draw: head=%p %g record=%d\n",head,t,record);
    if (!files.record && head) {
      /* just turned off */
      files.record=showclearM();
      buttons.record->value(files.record); // continue (will be turned on again)
    }
    Fl::flush();
    //    Fl::check();

    if (debug) fprintf(stderr,"                                          timeout repeated at %.4f\n",mytime());

    /* info message red-on-yellow pops up for 2.5 s */
    if (mymessage) {
      mymessage=0;
      Fl::repeat_timeout(2.5, timer_callback, userdata); }
    else {
      double t=(clock()-c0)/(double)CLOCKS_PER_SEC;
      if (t>2*tcycle) tcycle*=2;
      else if (t<0.5*tcycle) tcycle*=0.5;
      else tcycle=sqrt(tcycle*sqrt(tcycle*t));
      /* tcycle is the estimated CPU time of 1 cycle, smoothed */
      //      fprintf(stderr,"tcycle=%g\n",tcycle);
      // calculate the timeout: the minimum is speed.MINTIMEOUT/speed.FPS
      t=speed.timerdelay-tcycle;
      if (t<speed.MINTIMEOUT/speed.FPS) t=speed.MINTIMEOUT/speed.FPS;
      Fl::repeat_timeout(t, timer_callback, userdata); }
  }
public:
  Simulate(int X,int Y,int W,int H) : Fl_Box(X,Y,W,H) {
    // I cannot find this method in the manual, but without
    // box(FL_FLAT_BOX) the window smashes everything together
    box(FL_FLAT_BOX);
    color(FL_WHITE);
    //    color(FL_BACKGROUND2_COLOR); // gray
    Fl::add_timeout(speed.timerdelay, timer_callback, (void*)this);
  }
};
