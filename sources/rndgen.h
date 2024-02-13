/*
  common header file for random number generators rndgen0.c, rndgen1.c, ...

  NB: rnd() and rndcos() may give, albeit with tiny probablity ~1e-16, the
  interval end values. See rndgen.c for more info! Particular version may be
  guaranteed to give rnd() in [0,1) or (0,1), etc., though.
  
  Similarly, rndsph and rnddisk may give points by a tiny value outside the
  sphere or disk.
*/


// extern char rndinfo[80];           /* info string */
unsigned rndinit(int tablen,int seed);/* returns max Irnd+1 */
unsigned Irnd(void);               /* function returning random integer upto max */
double rnd(void);                  /* random number in [0,1] */
unsigned irnd(unsigned n);         /* random number in {0,...,n-1} */
double rndcos(void);               /* random number in [-1,1] */

/* see rndetc.c: */
double rndsph(double *r);          /* random vector in 1-sphere */
double rnddisk(double *r);         /* random vector in 1-disk */
void rndonsph(double *v);          /* random vector on 1-sphere */
double rndgauss(void);             /* normalized Gausian distribution */
