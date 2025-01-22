// The main simulation loop:
// Simulate.draw() calculates MC/MD and shows molecules
// Also, the timer control is here
class Simulate : public Fl_Box {
  void get_data () {
    int i;
    static double lastT;

    // get information from sliders; cf. cb_parms()
    slider2speed(sliders.speed->value());

    //    printf("%d %d %g\n",time(NULL),stride,timerdelay);
    block=pow(10,sliders.block->value())+0.5;

    T=exp(sliders.T->value());
    if (T!=lastT) {
      calculateB2();
      lastT=T; }

    d=exp(DSCALE*sliders.d->value());
    tau=exp(sliders.tau->value());

    if (bc==PERIODIC) gravity=0;
    else gravity=sliders.g->value();

    if (mauto && method<=MCNPT) {
      if (mcstart<=0) {
        // initial MC finished, switching to MD

        if (delayed_y_color) { // hack for DIFFUSION (color half AFTER MC)
          y_split();
          delayed_y_color=0; }

        MDreset(mauto);
        mauto=AUTO; /* = 0 */ }
      else
        mcstart-=speed.stride; }

    insdel((int)(pow(sliders.N->value(),NPOW)+0.5));

    if (isNPT(method)) {
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
        if (isMD(method) && (bc==BOX)) {
          if (debug) fprintf(stderr,"MDerror: r[%d] = [%g %g]\n",i,r[i].x,r[i].y);
          MDerror();
          r[i].x=rnd()*(L-1)+0.5; }
        while (r[i].x<0) r[i].x+=L;
        while (r[i].x>L) r[i].x-=L; }
      if ((r[i].y<=0) || (r[i].y>=L)) {
        if (isMD(method) && (bc!=PERIODIC)) {
          if (debug) fprintf(stderr,"MDerror: r[%d] = [%g %g]\n",i,r[i].x,r[i].y);
          MDerror();
          r[i].y=rnd()*(L-1)+0.5; }
        while (r[i].y<0) r[i].y+=L;
        while (r[i].y>L) r[i].y-=L; } }
  }

  void meas() {
    int i,j,k,l;
    unsigned ih;
    vector ri,rt;

    /* calculate virial pressure and RDF, by stride
       NB: in MD, Ekin is be known AFTER a MD step has finished
       NB: MSD calculated at the last block */
    En.Upot=0;
    En.P=zero; /* P=sum r*f = -virial */
    En.fwx=En.fwxL=En.fwy=En.fwyL=0;
    addP=measure==RDF?addPRDF:addPonly;
    loop (i,0,N) {
      ri=r[i];
      En.Upot+=gravity*r[i].y; // potential of gravity
      if (bc<PERIODIC) {
        double fwy=fwally(ri.y),fwyL=fwallyL(L-ri.y);
        En.P.y+=ri.y*fwy+(L-ri.y)*fwyL;
        En.fwy+=fwy;
        En.fwyL+=fwyL;
        En.Upot+=uwally(ri.y)+uwallyL(L-ri.y);
        if (bc==BOX) {
          double fwx=fwallx(ri.x),fwxL=fwallxL(L-ri.x);
          En.P.x+=ri.x*fwx+(L-ri.x)*fwxL;
          En.fwx+=fwx;
          En.fwxL+=fwxL;
          En.Upot+=uwallx(ri.x)+uwallxL(L-ri.x); } }
      switch (bc) {
        case BOX:
          loop (j,i+1,N)
            addP(ri.x-r[j].x,ri.y-r[j].y);
          break;
        case SLIT:
          loop (j,i+1,N)
            if (ss.C2<Lh)
              addP(ni(ri.x-r[j].x),ri.y-r[j].y);
            else
              loop (j,i+1,N)
                loopto (k,-nimg,nimg) {
                  rt.x=ri.x-r[j].x+L*k;
                  addP(rt.x,ri.y-r[j].y); }
          break;
        case PERIODIC:
          if (ss.C2<Lh)
            loop (j,i+1,N)
              addP(ni(ri.x-r[j].x),ni(ri.y-r[j].y));
          else
            loop (j,i+1,N)
              loopto (k,-nimg,nimg) {
                rt.x=ri.x-r[j].x+L*k;
                loopto (l,-nimg,nimg) {
                  rt.y=ri.y-r[j].y+L*l;
                  addP(rt.x,rt.y); } }
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

    //    if (method==METROPOLIS) Tk=T;
    //    En.P.x=(N*Tk+En.P.x)/Sqr(L);
    //    En.P.y=(N*Tk+En.P.y)/Sqr(L);
    // because of MTK, Ekin used instead of N*Tk (no "kinetic pressure correction")

    if (method==CREUTZ) En.Ekin=0; // this is not fully correct...
    if (method==METROPOLIS || method==MCNPT)
      En.Ek.x=En.Ek.y=En.Ekin=En.Nf*T/2;
    En.P.x=(En.Ek.x+En.P.x)/Sqr(L);
    En.P.y=(En.Ek.y+En.P.y)/Sqr(L);

    En.fwx/=L;
    En.fwxL/=L;
    En.fwy/=L;
    En.fwyL/=L;
  }

  void draw() { // in addition to drawing, all calculations are called from here
    int i;
    double V=0; // to suppress warning uninitilized

    get_data();
    if (run->value()) {
      loop (i,0,speed.stride) {

        if (i==speed.stride-1) {
          /* measurements BEFORE the last cycle so that Ekin is in sync
             with En.Upot AFTER this loop has finished */
          if (measure) {
            meas();
            sum.P.x+=En.P.x;
            sum.P.y+=En.P.y;
            sum.fwx+=En.fwx;
            sum.fwxL+=En.fwxL;
            sum.fwy+=En.fwy;
            sum.fwyL+=En.fwyL;
            V=N/L2rho(L);
            sum.V+=V; /* some corrections for not PERIODIC */
            sum.U+=En.Upot;
            if (method<=MCNPT) {
              /* the current bag + En.Upot is constant for CREUTZ */
              Econserved=En.Upot+bag;
              sum.Econserved+=Econserved;
              sum.Ekin+=bag; } } }

        if (isMD(method)) {
          MDstep();
          sum.Tk+=Tk/speed.stride; }
        else {
          accepted=Vaccepted=0;
          MCsweep();
          sum.accr+=(double)accepted/(N*speed.stride);
          sum.Vaccr+=(double)Vaccepted/speed.stride; }
      } // speed.stride

      /* last cycle finished, both En.Upot and Ekin are in sync */
      if (method==NVE) {
        Econserved=En.Upot+En.Ekin;
        sum.Econserved+=Econserved; }

      if (isNPT(method))
        En.H=En.Upot+En.Ekin+P*V;
      else
        En.H=En.Upot+En.Ekin+(En.P.x+En.P.y)*V/2;

      sum.H+=En.H;

      if (isMD(method))
        sum.Ekin+=En.Ekin;
      if (method==NOSE_HOOVER || method==MTK)
        sum.Econserved+=Econserved+En.Upot;

      iblock++;

      if (files.record) {
        double Ek=isMC(method)?bag:En.Ekin;
        if (!measureVar) {
          sum0.V=V;
          sum0.Upot=En.Upot;
          sum0.H=En.H;
          sum0.Ekin=Ek; }

        sums.V+=V-sum0.V; sumq.V+=Sqr(V-sum0.V);
        sums.H+=En.H-sum0.H; sumq.H+=Sqr(En.H-sum0.H);
        sums.Upot+=En.Upot-sum0.Upot; sumq.Upot+=Sqr(En.Upot-sum0.Upot);
        sums.Ekin+=Ek-sum0.Ekin; sumq.Ekin+=Sqr(Ek-sum0.Ekin);
        measureVar++; } }

    else
      // simulation not running due to button [run] in off state
      usleep(speed.extradelay+5000);

    norm();

    if (drawmode->value()<3) {
      // draw background
      enum colormode_e cmode=(colormode_e) colormode->value();
      static int boxsize=BOXSIZE; // size of the box in pix (will change)

      if (cmode==0) { molblack(); cmode=CM_KEEP; colormode->value(cmode); }

      if (drawmode->value()==0) Fl_Box::draw();

      shownbrs=0;
      switch (cmode) {
        case CM_NEIGHBORS:
          shownbrs=1;
          neighbors();
          break;
        case CM_YSPLIT:
          y_split();
          cmode=CM_KEEP;
          colormode->value(cmode);
          break;
        case CM_RANDOM:
          randomcolors();
          cmode=CM_KEEP;
          colormode->value(cmode);
          break;
        case CM_ONERED:
          molcol[0]=OI_BLACK;
          loop (i,1,N) molcol[i]=OI_ORANGE;
          cmode=CM_KEEP;
          colormode->value(cmode);
        default: // CM_KEEP, CM_BLACK
          break; }

      fl_color(FL_WHITE);
      if (drawmode->value()==1)
        loop (i,0,Sqr(boxsize)/64) fl_point(rnd()*boxsize+BORDER,rnd()*boxsize+MENUH+BORDER);

      boxsize=w()<h() ? w()-2*BORDER : h()-2*BORDER;
      double scale=boxsize/L;

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
      cmode=CM_BLACK; // the same color
      // I do not understand the following line :(
      loop (i,1,N) if (molcol[i]!=molcol[0]) { cmode=CM_YSPLIT; break; }
      if (cmode==CM_BLACK) fl_color(molcol[0]);

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
    wasMDerror=0;
    if (errmessage) {
      int atx=BORDER+6,aty=BORDER+MENUH+35;

      fl_color(OI_YELLOW);
      fl_rectf(BORDER,BORDER+MENUH,640,90);
      fl_font(FL_HELVETICA,32);
      fl_color(OI_RED);
      // see Fl::repeat_timeout() below
      if (errmessage==MDFAILED) {
        wasMDerror++;
        fl_draw("Molecular Dynamics failed!",atx,aty); aty+=40;
        fl_draw("Switching temporarily to Monte Carlo…",atx,aty); }
      else if (errmessage==MSDFAILED) {
        fl_draw("Cannot follow periodic b.c.!",atx,aty); aty+=40;
        fl_draw("MSD turned off for now…",atx,aty); }
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

    if (isMC(method) && !setd->value())
      sliders.d->show();
    else
      sliders.d->hide();

    if (isMD(method) && method!=NVE)
      sliders.tau->show();
    else
      sliders.tau->hide();

    if (isMC(method))
      setd->show();
    else
      setd->hide();

    if (isNPT(method)) {
       sliders.P->show();
       sliders.rho->hide(); }
    else {
       sliders.P->hide();
       sliders.rho->show(); }

    if (method==CREUTZ || method==NVE)
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

#if 0 // OLD
    /* this sets t=0 for the whole 1st block */
    if (files.record && head==NULL) t=0;
#endif

    if (files.record && !files.record0) {
      int i;

      MSDskip=0;

      // initialization of MSD
      loop (i,0,N) {
        r0[i].x=rpbc[i].x=r[i].x/L;
        r0[i].y=rpbc[i].y=r[i].y/L; }
      t=0; }

    files.record0=files.record;

    if (files.record && measure==NONE) measure=QUANTITIES;
    //    fprintf(stderr,"draw: head=%p %g record=%d\n",head,t,files.record);
    if (!files.record && head) {
      /* just turned off */
      files.record=showclearM();
      buttons.record->value(files.record); // continue (will be turned on again)
    }
    Fl::flush();
    //    Fl::check();

    if (debug) fprintf(stderr,"                                          timeout repeated at %.4f\n",mytime());

    /* info message red-on-yellow pops up for 2.5 s */
    if (errmessage) {
      errmessage=NOERROR;
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
