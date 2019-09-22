/*********************************************************
*  diff.c:  Differentiation with forward, central and    *    
*           Richardson extrapolation difference methods  *
*********************************************************/

#include <stdio.h>
#include <math.h>

#define h       1e-3                /* stepsize for all methods */
#define xmax    7                   /* range for calculation */
#define xmin    0
#define xstep   0.1          	    /* stepsize in x  */

void main()
{ 
   double dc, result, x;
   double f(double);			    /* function for differentiation */ 
   double fprime(double);		    /* exact derivative */ 
					               					
   FILE *output;			    /* save data in diff.dat */
   output = fopen("diff.dat","w");
   
   fprintf(output,"x         \t");                /* print header */
   fprintf(output,"Exact     \t");                    
   fprintf(output,"Forward   \t");   
   fprintf(output,"Central   \t");   
   fprintf(output,"Richardson Extrap \n\n");    /* Richardson Extrapolation, D2(h/4)  */

   for(x=xmin; x<=xmax; x+=xstep)
   {

      fprintf(output,"%f\t", x);
      
      result=fprime(x);                       /* exact result */
      fprintf(output, "%.10f\t", result);
 
      result=(f(x+h)-f(x))/h;                 /* forward difference */
      fprintf(output, "%.10f\t", result);
   
      result=(f(x+h/2)-f(x-h/2))/h;           /* central difference */
      fprintf(output, "%.10f\t", result);
   
      result=(8.0*(f(x+h/4.0)-f(x-h/4.0))-(f(x+h/2.0)-f(x-h/2.0)))/(3.0*h);  /* Richardson Extrap, D2(h/4) */
      fprintf(output, "%.10f\n", result);     /* extrapolated diff */
   }
   printf("data stored in diff.dat\n");
   fclose(output);
}

/*  end of main program */

/* the function for differentiation */

double f(double x)		
{
   return(x*exp(x));
}

double fprime(double x)		
{
   return(exp(x)+x*exp(x));
}