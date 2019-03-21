#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double *data;
double rslt[3];

extern void ferr(double *,int,double *);

void main(argc,argv)
int argc;
char *argv[];
{
  FILE *fptr;
  double tdata,xdata;
  int i,j,beg,end,ndata,blk,eflag;

  if(argc < 5)
    {
      printf("avg file blksize beg end xdata(opt)\n");
      exit(1);
    }

  sscanf(argv[2],"%d",&blk);
  sscanf(argv[3],"%d",&beg);
  sscanf(argv[4],"%d",&eflag);

  if( (fptr = fopen(argv[1],"r")) == NULL)
    {
      if(eflag != -2) printf("could not open the file %s to read\n",argv[1]);
      exit(1);
    }

  j = 0;
  while(fscanf(fptr,"%lf\n",&tdata) != EOF) j++;
  fclose(fptr);
  if(eflag == -1) end = j;
  else if(eflag == -2)
    {
      sscanf(argv[5],"%lf",&xdata);
      end = j;
    }
  else end = eflag;

  if(j < end || beg > end)
    {
      if(eflag != -2) printf("data file does not long enough or beg > end \n");
      exit(1);
    }
  ndata = end-beg+1;
  /*  printf("meas of %s \n total data = %d \n",argv[1],end-beg+1); */

  ndata /= blk;
  if(ndata == 0)
    {
      if(eflag != -2) printf("no data in the given range!\n");
      exit(1);
    }


  fptr = fopen(argv[1],"r");
  data = (double *)malloc(ndata*sizeof(double));

  for(i=1;i<beg;i++) fscanf(fptr,"%lf\n",&tdata);
  i=0;
  j=0;
  while(1)
    {
      if(j == 0) data[i]=0.;
      if(fscanf(fptr,"%lf\n",&tdata) == EOF)
        {
          if(eflag != -2) printf("premature end of file detected!\n");
          exit(1);
        }
      data[i] += tdata;
          j++;
      if( j == blk)
        {
          data[i] /= blk;
          j=0;
          i++;
        }
      if (i == ndata) break;
    }
  fclose(fptr);

  ferr(data,ndata,rslt);
  if(eflag == -2) printf("%24.16g \t %24.16g \t %24.16g\n",
                         xdata,rslt[0],rslt[1]);
  else printf("%24.16g \t %24.16g \t %24.16g blks: %d\n",
                         rslt[0],rslt[1],rslt[2],ndata);

  free(data);

}

float   cut_frct=0;     /* when non-zero, this is the fraction of the   */
                        /* length of the data sample to be used as a    */
                        /* cutoff for the relative time displacement in */
                        /* integration determining the autocorrelation  */
                        /* time.  When zero we continue to integrate    */
                        /* out to the place the correlation becomes     */
                        /* negative.                                    */

void ferr(ans,num,rslt)
/******************************************************************
 'ans' points to an array of 'num' floating point numbers to be averaged.
'rslt' points to a block of three floating point numbers where the results
are written: first the average, then the error, then the correlation length
******************************************************************/
int num;
double *ans,*rslt;

{
double cor=0,delta,sigma,avg=0,sqrt(),corlen;
int displ,n,
        max_displ;      /* the maximum displacement used in computing   */
                        /* auto correlation time.                       */
char prnt[3];  /* prnt != 'y' turns off printing*/

if(num == 0)
        {
        rslt[0]=0; rslt[1]=0; rslt[2]=0;
        return;
        }
for(n=0;n<num;n++) avg += ans[n];
avg = avg/num;
rslt[0]=avg;

if(prnt[0]=='y' && prnt[1] == 'e') printf("\n avg = %le",avg);

if(cut_frct==0) max_displ=num;
else
        {
        max_displ=(float)num*cut_frct;
        }
for(displ=0; displ<max_displ;displ++)
        {
        delta = 0;
        for(n=0;n<num-displ;n++) delta += (ans[n]-avg)*(ans[n+displ]-avg);
        if(displ) delta = 2*delta;
        else sigma = delta/num;
        if( (cut_frct==0)&&delta<=0 ) break;
        cor += delta/(num-displ);
        if(prnt[0]=='y' && prnt[1] == 'e')
                printf("\ndispl = %d\t cor = %le  delta = %le",
                        displ,sqrt(cor),sqrt(delta/(num-displ) ) );
        }
if(sigma) corlen=cor/sigma;
        else corlen=0;

if(prnt[0]=='y' && prnt[1] == 'e')
        printf("\nerror = %le   corr # = %d  corr len = %lf",
                sqrt(cor/num),displ,corlen);
rslt[1]=sqrt(cor/num);
rslt[2]=corlen;
}

void wrtofile(fname,index,data,ndata)
char *fname;
int index,ndata;
double *data;
{
  FILE *ofptr;
  char ofname[20];
  int i;

  sprintf(ofname,"%s_nt%d",fname,index);
  ofptr = fopen(ofname,"w");

  for (i=0;i<ndata;i++)
    fprintf(ofptr,"%13.10e\n",data[i]);
  fclose(ofptr);

  return;
}
