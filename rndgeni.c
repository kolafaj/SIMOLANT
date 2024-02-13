/* 

to be directly #included

random number generator - four-tap register-shift generator
^^^^^^^^^^^^^^^^^^^^^^^   grid of rnd() is at least 2^-53 (as fine as possible)
                          high quality and fast

initialization: rndinit(tablen,seed)
where tablen = irrelevant
      seed = any integer (0 uses seed from time)

RNDMETHOD is a sum of:
  rnd():
  0 = rnd in [0,1] 
      if cast to IEEE double, can be exactly: 
        0 with probability 2^-64=5.42e-20
        1 with probability 2^-54=5.55e-17
      if kept in register (10 bytes long double), can be exactly:
        0 with probability 2^-64=5.42e-20
        1 never
  1 = rnd in [0,1), IEEE numbers assumed 
      (i.e., exact (double)1 is never reached)
      
  rndcos():
  0 = rndcos in [-1,1] (can be exactly +1 or -1 with prob. < 1e-16)
  1 = rndcos in [-1,1), IEEE numbers assumed
  1+2 = rndcos in (-1,1), IEEE numbers assumed
  
  irnd:
  0: based on 32 bits (simpler and faster, good enough for most cases)
  4: based on 64 bits (rarely important)

NOTE: using RNDMETHOD>=2 is deprecated because a portable code should NOT
depend on whether a particular real value is reached. Note that generally
code like "double x,y; x=y; if (x==y)" is NOT guaranteed to pass because the 
first x may be kept in a register (10 bytes long) and the second not.

irnd(N) uses integer aritmetic and IS guaranteed to give values 
in {0,..,N-1} and never N
*/

#ifndef RNDMETHOD 
#define RNDMETHOD 0
#endif

#include "rndgen.h"
#include "rndseed.c"

#define A 471  
#define B 1586 
#define C 6988 
#define D 9689 
#define M 16383
#define trnd (++nd,ra[nd&M]=ra[(nd-A)&M]^ra[(nd-B)&M]^ra[(nd-C)&M]^ra[(nd-D)&M])

static unsigned ra[M+1],nd;

/* now a higher-quality generator (rndgen1) used as initializer */
#define rndX2 (int)113L      /* any not too big prime */
#define rndM2 (int)19004269L /* max prime < 2^31/rndX2 */
#define rndX1 (int)78125L    /* so that the cycle mod 2^32 is max (2^29) */

struct irndgen_s { int seed1,seed2,fac,tab[1/*var len*/]; };

unsigned initIrnd(struct irndgen_s *rndgen)
{
 int *p,r;

 p = rndgen->tab + (rndgen->seed2 = (rndgen->seed2*rndX2)%rndM2) / rndgen->fac;
 r = *p;
 *p = rndgen->seed1 *= rndX1;

 return (unsigned)r;
}

unsigned rndinit(int tablen,int seed)
{
  struct irndgen_s *rndgen;
  int i;

  if (seed==0) seed=seed_from_time();

  if (tablen<=0) tablen=47;
  if (tablen>1023) tablen=1023;
#ifdef alloc  
  alloc(rndgen,sizeof(struct irndgen_s)+tablen*sizeof(int));
#else  
  rndgen=(struct irndgen_s*)malloc(sizeof(struct irndgen_s)+tablen*sizeof(int));
#endif

  rndgen->seed1 = seed | 1; /* to be odd */
  rndgen->seed2 = abs(seed)%(rndX2-1)+1; /* to be in [1,rndX2-1] */
  rndgen->fac = rndM2/tablen+1;
  
  for (i=0; i<7; i++) rndgen->seed1 *= rndX1;
  
  for (i=0; i<tablen; i++) rndgen->tab[i] = rndgen->seed1 *= rndX1;

  for (nd=0; nd<=M; nd++) ra[nd]=initIrnd(rndgen);
  /* to have all bits sufficiently independent on seed */
  for (nd=0; nd<=M; nd++) ra[nd]^=initIrnd(rndgen)>>7;
  for (nd=0; nd<=M; nd++) ra[nd]^=initIrnd(rndgen)>>13;
  for (nd=0; nd<=M; nd++) ra[nd]^=initIrnd(rndgen)>>19;

  free(rndgen);
  return 0xffffffffU;
}

#undef rndX2
#undef rndX1
#undef rndM2

unsigned Irnd(void)
{
  return trnd;
}

unsigned irnd(unsigned N)
{
#if RNDMETHOD&4
  unsigned t0=((long long unsigned)(trnd)*N)>>32;
  return ((long long unsigned)(trnd)*N+t0)>>32;
#else  
  return ((long long unsigned)(trnd)*N)>>32;
#endif  
}
  
double rnd(void)
{
#if RNDMETHOD&1
  /* random number in [0,1) */
  double r=(trnd&0xfffff800U)*(1./4294967296e0);
  return (r+(double)trnd)*(1./4294967296e0);
#else  
  /* random number in [0,1] */
  double r=trnd*(1./4294967296e0);
  return (r+(double)trnd)*(1./4294967296e0);
#endif
}

double rndcos(void)
{
#if RNDMETHOD&2
  /* random number in (-1,1) */
  double r=(trnd&0xfffff800U|0x3ff)*(1./4294967296e0);
  return (r+(double)trnd)*(1./2147483648e0)-1;
#elif RNDMETHOD&1
  /* random number in [-1,1) */
  double r=(trnd&0xfffffc00U)*(1./4294967296e0); 
  return (r+(double)trnd)*(1./2147483648e0)-1;
#else  
  /* random number in [-1,1] */
  double r=trnd*(1./4294967296e0);
  return (r+trnd)*(1./2147483648e0)-1;
#endif
}

#undef A
#undef B
#undef C
#undef D
#undef M

#include "rndetc.c"
