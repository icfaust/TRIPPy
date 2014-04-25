#include <math.h>

void intercept2d(double outval[], double pt1[3], double pt2[3], double outliner[], double outlinez[], int ix)
{ 
  /* initialize variables */
  int i;
  double delr,delz,del1,del2,denom;
  *outval = INFINITY;

  delr = pt1[0] - pt2[0];
  delz = pt1[2] - pt2[2];

  for(i= ix-1 ;i--;)
    {

      del1 = outliner[i+1]-outliner[i];
      del2 = outlinez[i+1]-outlinez[i];

      /* If the inverse is possible  */
      denom = delr*del2 - delz*del1;
      if(denom)
	{
	  del1 = (del2*(pt1[0] - outliner[i]) - del1*(pt1[2] - outlinez[i]))/denom;
	  del2 = (delr*(pt1[2] - outlinez[i]) - delz*(pt1[0] - outliner[i]))/denom;
	  /*Does the point intercept, and is it smaller than the previous smallest value?  */
	  /* del2 is the length between outline at point i+1 to point i where the intercept occurs */
	  /* del1 is the length between pt1 and pt2 where the intercept occurs */

	  if((del2 > 0) & (del2 < 1) & (del1 < *outval) & (del1 > 0))
	    { 
	      *outval = del1;
	    }
	}
    }
}

void interceptCyl(double *s, double pt0[3], double norm[3], double outliner[], double outlinez[], int ix)
{ 
  /* parabolic intercepts of lines using quadratic formula */

  /* initialize variables */
  int i;
  double A,A0,B,B0,C,C0,delr,delz,temp,s1,s2;
  *s = INFINITY;

  C0 = pow(pt0[0],2) + pow(pt0[1],2);
  B0 = 2*(pt0[0]*norm[0]+pt0[1]*norm[1]);
  A0 = pow(norm[0],2) + pow(norm[1],2);

  for(i= ix-1 ;i--;)
    {
      A = 0;
      B = 0;
      C = 0;
      
      delr = outliner[i+1]-outliner[i];
      delz = outlinez[i+1]-outlinez[i];
     
      if(delz)
	{
	  temp = outliner[i] + (delr/delz)*(pt0[2]-outlinez[i]);
	  
	  A = A0 - pow((delr/delz)*norm[2],2);
	  B = B0 - 2*temp*(delr/delz)*norm[2];
	  C = C0 - pow(temp,2);
	 
	}
      else if((norm[2] != 0) & (delr != 0)) /* if is not a purely radial line, and points actually a line */
	{ 
	  /*prevents rest of quadratic iteration from occuring */
	  s1 = (outlinez[i] - pt0[2])/norm[2];
	  temp = sqrt(pow(pt0[0]+norm[0]*s1,2) + pow(pt0[1]+norm[1]*s1,2))/delr;
	  if((temp > 0) & (temp <= 1) & (s1 < *s) & (s1 > 0))
	    { 
	      *s = s1;
	    } 
	}
      
      if(A) /*the quadratic form*/
	{
	  temp = B*B - 4*A*C; /*reuse a variable, not exactly a good idea */
	  if(temp >= 0) /* if there is an intercept */
	    {
	      temp = sqrt(temp);
	      s1 = -.5*(temp + B)/A;
	      s2 = .5*(temp - B)/A;
	      temp = (s1*norm[2] + (pt0[2] - outlinez[i]))/delz; /*length along cylinder parameterization of intercept */
	      if((temp > 0) & (temp <= 1) & (s1 < *s) & (s1 > 0))
		{ 
		  *s = s1;
		}
	      
	      temp = (s2*norm[2] + (pt0[2] - outlinez[i]))/delz; /*length along cylinder parameterization of intercept */
	      if((temp > 0) & (temp <= 1) & (s2 < *s) & (s2 > 0))
		{ 
		  *s = s2;
		}
	      
	    }
	}
      else if((B != 0) & (C != 0)) /*for directly radial views */
	{
	  s1 = -B/C;
	  temp = (s1*norm[2]+(pt0[2]-outlinez[0]))/delz; /*length along cylinder parameterization of intercept */
	  if((temp > 0) & (temp <= 1) & (s1 < *s) & (s1 > 0))
	    { 
	      *s = s1;
	    }
	}	  
    }
}


