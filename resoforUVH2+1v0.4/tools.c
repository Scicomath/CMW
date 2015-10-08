/* tools.c            functions for reso.c main program calculations 
*   Josef Sollfrank                  Nov. .98             */


#include	<string.h>
#include	<stdio.h>
#include	<math.h>

#include	"tools.h"

/* converts an integer to a string !*/

double	lin_int(x1,x2,f1,f2,x)
        double x1, x2, f1, f2, x;
        
        {
         double aa, bb;
         
	 if (x2 == x1) 
           aa = 0.0;
	 else
           aa =(f2-f1)/(x2-x1);
         bb = f1 - aa * x1;
         
         return aa*x + bb;
        }


int    poweri(xx,l)
	int xx;
	int l;
	{ 
	 int i;
	 double re ;

	 re = xx ;
	 if (l == 0) return 1;
	 else {
	   for(i=1;i<l;i++) re = re*xx ; 
	   return re ;
	   }
	 }


void    convei(zahl,p)
	int zahl;
	char p[SSL];
	{
	  int dig, j, ii = 0, sta, ani, z ;
          double dum;

	  if (zahl < 0) {p[0] = 45; ++ii;} 
	  z = abs(zahl);

	  dum = log10( (double) z);
          ani = (int) dum;
	  sta = poweri(10,ani);

	  for(j=ii; j <= ani+ii; j++) {
	     dig = z / sta ;
	     p[j]  = dig + 48 ;
	     z = z - dig*sta;
	     sta = sta / 10 ;
	  }
	  p[j] = '\0';
	}


