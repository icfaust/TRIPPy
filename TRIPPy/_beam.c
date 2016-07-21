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

void interceptCyl(double s[], double pt0[][3], double norm[][3], double outliner[], double outlinez[], int ix, int jx)
{ 
  /* parabolic intercepts of lines using quadratic formula */

  /* initialize variables */
  int i,j;
  double A,A0,B,B0,C,C0,delr,delz,temp,s1,s2,stemp,res=1e-6;

 
  for(j=jx;j--;)
    {

      stemp = INFINITY;
      
      C0 = pow(pt0[j][0],2) + pow(pt0[j][1],2);
      B0 = 2*(pt0[j][0]*norm[j][0]+pt0[j][1]*norm[j][1]);
      A0 = pow(norm[j][0],2) + pow(norm[j][1],2);
      
      for(i= ix;i--;)
	{
	  A = 0;
	  B = 0;
	  C = 0;
	  
	  delr = outliner[i+1] - outliner[i];
	  delz = outlinez[i+1] - outlinez[i];
	  
	  if(delz)
	    {
	      temp = outliner[i] + (delr/delz)*(pt0[j][2] - outlinez[i]);
	      
	      A = A0 - pow((delr/delz)*norm[j][2],2);
	      B = B0 - 2*temp*(delr/delz)*norm[j][2];
	      C = C0 - pow(temp,2);
	      
	    }
	  else if((norm[j][2] != 0) & (delr != 0)) /* if is not a purely radial line, and points actually a line */
	    { 
	      /*prevents rest of quadratic iteration from occuring */
	      s1 = (outlinez[i] - pt0[j][2])/norm[j][2];
	      temp = (sqrt(pow(pt0[j][0] + norm[j][0]*s1,2) + pow(pt0[j][1] + norm[j][1]*s1,2)) - outliner[i])/delr;
	      
	      if((temp > 0) & (temp <= 1) & (s1 < stemp) & (s1 > res))
		{
		  stemp = s1;
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
		  temp = (s1*norm[j][2] + (pt0[j][2] - outlinez[i]))/delz; /*length along cylinder parameterization of intercept */
		  if((temp > 0) & (temp <= 1) & (s1 < stemp) & (s1 > res))
		    {
		      stemp = s1;
		    }
		  
		  temp = (s2*norm[j][2] + (pt0[j][2] - outlinez[i]))/delz; /*length along cylinder parameterization of intercept */
		  if((temp > 0) & (temp <= 1) & (s2 < stemp) & (s2 > res))
		    {
		      stemp = s2;
		    }
		  
		}
	    }
	  else if((B != 0) & (C != 0)) /*for directly radial views */
	    {
	      s1 = -B/C;
	      temp = (s1*norm[j][2] + (pt0[j][2] - outlinez[i]))/delz; /*length along cylinder parameterization of intercept */
	      if((temp > 0) & (temp <= 1) & (s1 < stemp) & (s1 > res))
		{
		  stemp = s1;
		}
	    }
	}

      s[j]= stemp;
    }

}

void lineCirc(double out[][5], double pt0[3], double norm[3], double r[], double z[],int ix)
{ int i;
  double eps,epst,f0r,f1r;
  
  f0r = pow(pt0[0],2) + pow(pt0[1],2);
  f1r = pow(norm[0],2) + pow(norm[1],2);
  eps = pt0[0]*norm[0] + pt0[1]*norm[1];

  for(i=0;i < ix; i++)
    { 
      epst = eps + (pt0[2] - z[i])*norm[2];
      out[i][0] = f1r;
      out[i][1] = 2*(eps + epst*f1r);
      out[i][2] = f0r + 4*eps*epst + f1r*(pow(epst,2) - pow(r[i],2)*f1r);
      out[i][3] = 2*(f0r*epst + eps*(pow(epst,2) - pow(r[i],2)*f1r));
      out[i][4] = f0r*pow(epst,2) - pow(r[i],2)*pow(eps,2);
    }
}

void bessel_fourier_kernel(double theta[], double m, double zero, double rho, double out[], int ix)
{ int i;

  for(i=ix;i--;)
    {
      out[i] = cos(m*theta[i])*sin(zero*(cos(theta[i]) - rho));
    }
}

void idx_add(double output[], int idx1[], int idx2[], double data[], double mult[], double ds, int lim2, int ix, int jx1, int jx2)
{ int i,j,lim1 = 0;

  for(i=0; i < ix; i++)
    {
      for(j=0; j < jx1; j++)
	{
	  output[lim2 + idx1[lim1 + j]] = output[lim2 + idx1[lim1 + j]] + mult[lim1 + j]*(ds - data[lim1 + j]);
	  output[lim2 + idx2[lim1 + j]] = output[lim2 + idx2[lim1 + j]] + mult[lim1 + j]*data[lim1 + j];
	}
      lim1 = lim1 + jx1;
      lim2 = lim2 + jx2;
    }
}


void idx_add2(double output[], int idx1[], int idx2[], double data[], double mult[], double ds, int ix, int jx1, int jx2)
{ int i,j,lim1 = 0,lim2 = 0;

  for(i=0; i < ix; i++)
    {
      for(j=0; j < jx1; j++)
	{
	  output[lim2 + idx1[lim1 + j]] = output[lim2 + idx1[lim1 + j]] + mult[lim1 + j]*data[lim1 + j];
	  output[lim2 + idx2[lim1 + j]] = output[lim2 + idx2[lim1 + j]] + mult[lim1 + j]*(ds - data[lim1 + j]);
	}
      lim1 = lim1 + jx1;
      lim2 = lim2 + jx2;
    }
}
