/*
  2D linked-cell list method
  cutoff<=Lh required
  simplified MACSIMUS code (version LINKCELL=3)
*/

linklist_t site_storage[MAXN];

/*
 * A linked-cell list is built according to the x,y coordinates of sites.
 * Memory needed is allocated in the 1st pass and re-used later.
 * Coordinates x,y are normalized to [0,L),
 * Type (linklist_t***) is interpreted as
   array[ncell][ncell] of type (linklist_t*)
 */
linklist_t ***linklist(void) /************************************* linklist */
{
  int i,ix,iy;
  static int oldncell;
  static linklist_t ***list=NULL; // 2D array of pointers to the linked lists
  linklist_t *l;

  if (ss.ncell<2) ss.ncell=2; // should never happen

  double rL=ss.ncell/L*(1-2e-16); // by 2 in the lsb less
  
  if (ss.ncell!=oldncell && list) {
    free(list);
    list=NULL; }

  if (!list) {
    /* 2D array of cells allocated at once, cf. MACSIMUS alloc2Darray() */
    list=(linklist_t***)malloc(sizeof(list[0])*ss.ncell);
    list[0]=(linklist_t**)malloc(sizeof(list[0][0])*Sqr(ss.ncell));
    loop (ix,0,ss.ncell) list[ix]=list[0]+ix*ss.ncell;
    oldncell=ss.ncell; }

  /* clear the linked lists */
  loop (ix,0,ss.ncell)
    loop (iy,0,ss.ncell) list[ix][iy]=NULL;

  l=site_storage;

  loop (i,0,N) {
    // cf. Simulate.norm()
    if (bc!=BOX) { // normalize periodic in x
      if (r[i].x<0) r[i].x+=L;
      if (r[i].x>=L) r[i].x-=L; }
    l->r.x=r[i].x;
    l->i=i;
    ix=(int)(r[i].x*rL);

    if (bc==PERIODIC) { // normalize periodic in y
      if (r[i].y<0) r[i].y+=L;
      if (r[i].y>=L) r[i].y-=L; }
    l->r.y=r[i].y;
    iy=(int)(r[i].y*rL);

    l->f=a+i;
    if (method==VICSEK) {
      l->v=v+i;
      l->vnew=vllast+i;
      l->wnbrs=wnbrs+i; }
    l->next=list[ix][iy];
    list[ix][iy]=l;
    l++; }

  return list;
} /* linklist */

static linklist_t lsloc;

/*** functions for the linked-list loop ***/

/* DEBUG
   To compare |xi-xj|, |yi-yj| obtained by the direct pair and linked list.
   Lines marked PC (Pair in Cutoff) and LC (Linked list in Cutoff) should match.
   PE and LE mark points calculated outside cutoff.
   Activated by option -D-1, exut after printout.
*/
void debuglc(linklist_t *l) /*************************************** debuglc */
{
  for (; l; l=l->next) {
    double rr=Sqr(lsloc.r.x-l->r.x)+Sqr(lsloc.r.y-l->r.y);
    printf("%16.11f %16.11f %d %d L%c\n",
           fabs(lsloc.r.x-l->r.x),fabs(lsloc.r.y-l->r.y),
           lsloc.i,l->i,
           "EC"[rr<ss.C2q]); }
  return;
}

/* for MD, pressure tensor calculated (Pneeded) */
void MDlcPneeded(linklist_t *l) /******************************* MDlcPneeded */
{
  vector rt;
  double ff;
  
  for (; l; l=l->next) {
    rt.x=lsloc.r.x-l->r.x;
    rt.y=lsloc.r.y-l->r.y;
    ff=ss.lcfunc(Sqr(rt.x)+Sqr(rt.y));
    if (ff) {
      En.P.x+=Sqr(rt.x)*ff;
      En.P.y+=Sqr(rt.y)*ff;
      lsloc.f->x+=rt.x*ff; lsloc.f->y+=rt.y*ff;
      l->f->x-=rt.x*ff; l->f->y-=rt.y*ff; } }
}

/* for MD, pressure tensor not calculated (!Pneeded) */
void MDlc(linklist_t *l) /********************************************* MDlc */
{
  vector rt;
  double ff;
  
  for (; l; l=l->next) {
    rt.x=lsloc.r.x-l->r.x;
    rt.y=lsloc.r.y-l->r.y;
    ff=ss.lcfunc(Sqr(rt.x)+Sqr(rt.y));
    if (ff) { // does this pay off?
      lsloc.f->x+=rt.x*ff; lsloc.f->y+=rt.y*ff;
      l->f->x-=rt.x*ff; l->f->y-=rt.y*ff; } }
}

/* for MD, pressure tensor not calculated (!Pneeded) */
void MDlcVicsek(linklist_t *l) /********************************* MDlcVicsek */
{
  vector rt;
  double ff,rr;
  
  for (; l; l=l->next) {
    rt.x=lsloc.r.x-l->r.x;
    rt.y=lsloc.r.y-l->r.y;
    rr=Sqr(rt.x)+Sqr(rt.y);
    ff=ss.lcfunc(rr);
    if (ff) {
      lsloc.f->x+=rt.x*ff; lsloc.f->y+=rt.y*ff;
      l->f->x-=rt.x*ff; l->f->y-=rt.y*ff;
    
      /* Vicsek: */
      rr=ss.C2q-rr;
      l->vnew->x+=rr*lsloc.v->x;
      l->vnew->y+=rr*lsloc.v->y;
      lsloc.vnew->x+=rr*l->vnew->x;
      lsloc.vnew->y+=rr*l->vnew->y;
      *lsloc.wnbrs+=rr;
      *l->wnbrs+=rr; } }
}

/* for MCNPT; Metropolis displacement not implemented via linklist (yet) */
void MCNPTlc(linklist_t *l) /*************************************** MCNPTlc */
{
  vector rt;
  double rr;
  
  for (; l; l=l->next) {
    rt.x=lsloc.r.x-l->r.x;
    rt.y=lsloc.r.y-l->r.y;
    rr=Sqr(rt.x)+Sqr(rt.y);
    ss.DeltaU+=ss.lcfunc(rr*ss.rrscale)-ss.lcfunc(rr); }
}

/*
more versions of linked-cell list would be more difficult to implement:
- for measure: Simulate::meas(), the problem is with a longer range for RDF
- for MC: difficult - other type of algorithm would be needed
x for MCNPT: this would be doable by using the maximum cutoff scaled by L,Lnew
*/

void lc(DO_f DO) /******************************************************* lc */
{
  linklist_t ***list=linklist(),**lxs,*ls;
  double rL=ss.ncell/L*(1-2e-16); // by 2 in the lsb less
  double cell=L/ss.ncell*(1-2e-16); // cellsize,cell=1/rL
  double DX,dy,dyq;
  int ix,iy,ixs,iys,ixm,iym,iy0,ix1,iy1;
  vector cellr;

  /* pair sum for testing */
  if (DO==debuglc) {
    printf("# box=%.13g rl=%.14g\n",L,rL);
    loop (ix,0,N) {
      printf("%16.11f %16.11f %3d cfg\n",r[ix].x,r[ix].y,ix);
      loop (ixm,ix+1,N) {
        double xx=0,yy=0;
        switch (bc) {
          case BOX:      xx=Sqr  (r[ix].x-r[ixm].x); yy=Sqr  (r[ix].y-r[ixm].y); break;
          case SLIT:     xx=Sqrni(r[ix].x-r[ixm].x); yy=Sqr  (r[ix].y-r[ixm].y); break;
          case PERIODIC: xx=Sqrni(r[ix].x-r[ixm].x); yy=Sqrni(r[ix].y-r[ixm].y); break; }
        printf("%16.11f %16.11f %d %d P%c\n",
               sqrt(xx), sqrt(yy),
               ix,ixm,
               "EC"[xx+yy<ss.C2q]); } } }

  loop (ixs,0,ss.ncell) {
    lxs=list[ixs];
    loop (iys,0,ss.ncell) {
      for (ls=lxs[iys]; ls; ls=ls->next) {
        /* central site *ls selected:
           Now draw a circle of radius C2 around it and find matching cells.
           `lsloc' is the local copy of the central cell, *ls;
           see MACSIMUS lc.C and lcforce-nolsloc.c for more info. */
        lsloc=*ls;
        cellr.x=rL*lsloc.r.x;
        cellr.y=rL*lsloc.r.y;
        ix1=cellr.x+ss.C2*rL;

        if (bc==BOX && ix1>=ss.ncell) ix1=ss.ncell-1; // not periodic in x

        loopto (ix,ixm=ixs,ix1) {
          if (ix==ixs) DX=0;
          else DX=ix-cellr.x;
          if (ixm==ss.ncell) { // cannot happen for BOX
            lsloc.r.x-=L; ixm-=ss.ncell; }
          dyq=ss.C2q-Sqr(DX*cell); /* OPT: dyq can be replaced by dy */

          if (dyq<0) {
            if (dyq/ss.C2q<-1e-12) fprintf(stderr,"linkcell problem: dyq=%g, set to 0\n",dyq);
            if (dyq/ss.C2q<-1e-12) printf("linkcell problem: dyq=%g, set to 0\n",dyq);
            dyq=0; }

          dy=sqrt(dyq)*rL;
          iy1 = (int)(cellr.y+dy);
          if (bc!=PERIODIC && iy1>=ss.ncell) iy1=ss.ncell-1; // not periodic in y

          if (ixm==ixs) {
            iym = iy0 = iys;
            loopto (iy,iy0,iy1) {
              if (iym==ss.ncell) { // cannot happen for BOX|SLIT
                lsloc.r.y-=L; iym-=ss.ncell; }
              if (iym==iys) /* the same cell */
                DO(lsloc.next);
              else
                DO(list[ixm][iym]);
              iym++; } }
          else {
            iym = iy0 = (int)(cellr.y+ss.ncell-dy)-ss.ncell;
            if (iy0<0) {
              if (bc!=PERIODIC)
                iym=0; // BOX|SLIT
              else {
                lsloc.r.y+=L; iym+=ss.ncell; } }

            loopto (iy,iy0,iy1) {
              if (iym==ss.ncell) { // cannot happen for BOX|SLIT
                lsloc.r.y-=L; iym-=ss.ncell; }
              DO(list[ixm][iym]);
              iym++; } }
          lsloc.r.y=ls->r.y; /* restore original */
          ixm++; } /*ix*/
        //        lsloc.r.x=ls->r.x; /* to omit this: lsloc.r no longer needed */
        } /*ls*/
      } /*iys*/
    } /*ixs*/
  if (DO==debuglc) exit(0);
}
