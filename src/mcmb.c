
/****************************************************************/
/*  Markov Chain Marginal Bootstrap Package for Quantile Regression   */ 
/*  Maria Kocherginsky <mkocherg@uchicago.edu>,Xuming He <he@uiuc.edu>*/ 
/*                                                                    */
/*  This code modified by R. Koenker (July 30, 2003, April 21 2010)   */ 
/*  to employ R random number generation conventions and to conform   */ 
/*  to the conventions of the quantreg package.  See the function     */ 
/*  boot.rq.mcmb for details on the interface to summary.rq           */ 
/**********************************************************************/
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>

#define MAXN 200000
#define MAXP 100
#define XACCF .001
#define XACCD .001

#define MAXIT 100
#define BMAXIT 100

#define UNUSED (-1.11e30)
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50


int allZero;

/****************************************************************/
/*  Error Handling      */ 
/****************************************************************/
void error(const char * format, ...);
/****************************************************************/
/*  Various Utilities off the Web (for use with sort2())        */ 
/****************************************************************/

unsigned int *lvector(int  nl, int  nh)
/* allocate a int vector with subscript range v[nl..nh] */
{
    unsigned int *v = (unsigned int *)malloc((size_t)((nh-nl+2) * sizeof(int)));
    if (v == NULL) error("allocation failure in lvector()");
    return (v-nl+1);
}

void free_lvector(unsigned int *v, int nl, int nh)
/* free a int int vector allocated with lvector() */
{
    free(v+nl-1);
}

/****************************************************************/
/*  Multiply two vectors                                        */ 
/****************************************************************/

double mprodx(double *x, double *c, int pp){
  int i;
  double sum;
  
  sum=0;
  for(i=0;i<pp;i++){
    sum=sum+x[i]*c[i];}
  
  return sum;
}

/****************************************************************/
/*  Sign Function                                               */ 
/****************************************************************/

double sign(double x){
  double sign;

  if(x>0.0)
    sign=1;
  else if (x<0.0)
    sign=-1;
  else 
    sign=0;	
  
  return sign;
}

/****************************************************************/
/*  Sort Function                                               */ 
/****************************************************************/

void sort2(unsigned int  n, double arr[], double brr[]){
  unsigned int  i,ir=n,j,k,l=1,*istack;
  int jstack=0;
  double a,b,temp;
  istack=lvector(1,NSTACK);

  for (;;) { 
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (!jstack) {
	free_lvector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack]; 
      l=istack[jstack-1];
      jstack -= 2;
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
      SWAP(brr[k],brr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	SWAP(brr[l],brr[ir])
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	SWAP(brr[l+1],brr[ir])
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	SWAP(brr[l],brr[l+1])
      }
      i=l+1; 
      j=ir;
      a=arr[l+1]; 
      b=brr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break; 
	SWAP(arr[i],arr[j]) 
	SWAP(brr[i],brr[j])
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      brr[l+1]=brr[j];
      brr[j]=b;
      jstack += 2;
      if (jstack > NSTACK) error("NSTACK too small in sort2.\n");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}



/****************************************************************/
/* This function finds a weighted taustar quantile, instead of  */
/* using bisection to solve the MCMB equations. See new paper.  */
/****************************************************************/

double func(double *x, double *y,  double tau, double *tTilda, double *A,  double 
sum_right, double sumxij, double sumabsxij, int j, int pp, int nn){
  int i,m;
  double *xj, *yj, *z, *wt, wtsum;
  unsigned int mm;
  double taustar, pwtsum, ans, large;

  xj=(double *) calloc(nn+1, sizeof(double));
  yj=(double *) calloc(nn+1, sizeof(double));  
  z=(double *) calloc(nn+2, sizeof(double));  
  wt=(double *) calloc(nn+2, sizeof(double));

  
  for(i=0;i<nn;i++){
    yj[i]=y[i];
    xj[i]=x[i*pp+j];
  }
  xj[nn]=-sum_right/tau;
  yj[nn]=10e16;
  wtsum=sumabsxij + fabs(xj[nn]);
  
  /*update the first n elements*/

  /*sort function sorts arr[1...n], so set the first elements of z and wt to 0, so that the actual values in z and wt start from position z[1] and wt[1]*/
  m=1;
  z[0]=0; wt[0]=0; 
  for(i=0;i<nn;i++){
    if(fabs(xj[i]) > 10e-16){
      z[m]=(y[i]-mprodx(&x[i*pp],tTilda,pp)+tTilda[j]*xj[i])/xj[i];
      wt[m]=fabs(xj[i])/wtsum;
      m=m+1;
    }
    else{error("fabs(xj[i])<10e-16\n");}
  }
  z[m]=10e16*sign(xj[nn]);
  large=z[m];
  wt[m]=fabs(xj[nn])/wtsum;
  
  /*calculate taustar*/
  taustar=(tau-0.5)*(sumxij+xj[nn])/(wtsum)+0.5;
  
  if(m==0){
    error("Error: one design variable contains all 0s.\n");
    allZero=1;}
  
  mm=m;
  sort2(m, z, wt);

  pwtsum=0; i=1;
  ans=z[1];
  while(pwtsum<=taustar & i<nn+1){
    pwtsum=pwtsum+wt[i];
    ans=z[i];
    i++;
  }

  /*resample if pick the largest z*/ 
  if(fabs(ans) > 10e15){
    error("Picked infinity; need to resample\n");
    return 1.0;
  }
  

  free(xj);
  free(yj);
  free(z);
  free(wt);

  return ans;

}

/****************************************************************/
/* RQMCMB                                                       */
/****************************************************************/


void bootnp(double *x, double *y,  double *tau, double *theta_tilda, double *A, double *zstar, 
	    double *sumxij, double *sumabsxij, int *n, int *p, int *success, double *theta, int *MAXK, int *seed){
  
  int i, j,  jj, k, nn, pp;
  double sum, s[MAXP],tau2, tTilda[MAXP];
  int rand_ind;
  extern int allZero;
  
  pp=(int) *p;
  nn=(int) *n;
  tau2=(double) *tau;
  allZero=0;
  
  
  for(i=0;i<pp;i++){
    theta[i]=theta_tilda[i];
    tTilda[i]=theta_tilda[i];
  }
  
  GetRNGstate();
  success[0]=1;  
  
  for(k=0;k<*MAXK; k++){
    
    /*bootstrap from x[i][j] individually for each j, return sum */
    for(j=0;j<pp;j++){
      sum=0; 
      for(i=0;i<nn;i++){
        rand_ind=(int)(nn*j+nn*unif_rand());
	sum=sum+zstar[rand_ind];
      }
      s[j]=sqrt(nn)/sqrt(nn-pp)* sum;
    }
  PutRNGstate();
        
    for(j=0;j<pp;j++){
      theta[(k+1)*pp+j]=func(x, y, tau2, tTilda, A, s[j], sumxij[j], sumabsxij[j], j, pp, nn);
      if(allZero==1){
	success[0]=0;
	return;}
      if(theta[(k+1)*pp+j]==1.0){
        for(jj=0;jj<pp;jj++){
          tTilda[jj]=theta[k*pp+jj];}
        k=k-1;
        break;
      }
      else { 
	tTilda[j]=theta[(k+1)*pp+j];
      }
    }
  }
  
}












