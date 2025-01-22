// Timing of animations support

/* range for sliders.speed is [MINSPEED,MAXSPEED+1.5]:
   [MINSPEED,1] ⇒ delays added to slow the motion
   [1,MAXSPEED] ⇒ stride (every sliders.speed-th configuration shown)
                  frame rate = FPS 
   >MAXSPEED+1 ⇒ stride=MAXSPEED and the maximum possible speed 
                 (may by shaky for small systems) */
#define MINSPEED -5
#define MAXSPEED 10

static struct speed_s {
  double FPS=60; // frames per second for showing (ref. value): see -F
  double extradelay=1e-5; // [s] delay added - is this needed? : see -U
  double timerdelay; // timer delay (tcycle to be subtracted, if possible)
  double init=3.5; // initial value of the speed slider: see -S
  double MINTIMEOUT=0.2; // minimum timeout is MINTIMEOUT*tcycle: see -T
  int stride; // every stride-th frame shown
} speed;

// calculate the stride and timer-delay from the speed slider
static void slider2speed(double val) /************************* slider2speed */
{
  if (val>=1) {
    speed.stride=(int)val;
    speed.timerdelay=1./speed.FPS;
    if (speed.stride==MAXSPEED+1) { // max speed
      speed.stride=MAXSPEED;
      speed.timerdelay=0; } }
  else {
    speed.stride=1;
    speed.timerdelay=exp((-val+1)/2)/speed.FPS; }
  //  fprintf(stderr,"%d %g\n",speed.stride,speed.timerdelay);
} // slider2speed()
