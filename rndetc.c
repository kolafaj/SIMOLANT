double rndgauss(void) /******************************************** rndgauss */
/* normal Gaussian distribution, exact version */
{
  static unsigned phase;
  static double gauss2;

  if (phase++&1)
    return gauss2;
  else {
    double x=rndcos()*PI,y=sqrt(-2*log(rnd()));
    gauss2=sin(x)*y;
    return cos(x)*y; }
}
