// from the standard JK interface

#define Sqr(X) ((X)*(X))
#define Cub(X) ((X)*(X)*(X))

/* C-style loop (last not included) */
#define loop(I,FROM,TO) for ((I)=(FROM); (I)<(TO); (I)++)

/* Pascal etc. style loop (last included) */
#define loopto(I,FROM,TO) for ((I)=(FROM); (I)<=(TO); (I)++)

/* loop over linked list.  Example:
   struct mylist_s {
     struct mylist_s *next; // next item; NULL terminates the list
     char text[8]; } 
   *head,*l;
   looplist (l,head) printf("%s\n",l->text);
*/
#define looplist(PTR,HEAD) for ((PTR)=(HEAD); (PTR); (PTR)=(PTR)->next)

#define put_(_X) printf("%11s=%-13.6g",#_X,(double)(_X));
#define put(_X) printf("%11s=%-13.6g\n",#_X,(double)(_X));

