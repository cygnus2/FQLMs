/* construct gauge inequivalent states of the U(1) model  on triangular lattice */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include "define.h"
//#include <algorithm> //std::count
//#include <vector>    //std::vector
// For some notes see below

void conststates(){
  extern void add(char*);
  extern void addATpos(char*,int);
  extern int checkGL(bool*, bool*, int);
  int count,fluxcount;
  int flagFLUX;
  int base=6;
  long long int n,NST;
  char str[VOL2];
  int i,j,nPOINT;
  int ind_lin,ind_chk;
  int linkset[VOL];

  bool conf[2*VOL],isSET[2*VOL];
  int flagGI,flagJUMP;
  int p,pp,l,cc;
  int q,xf,yf,xb,yb;
  int q0x,q0y,qx,qy;
  int Q,e1,e2,e3,e4;
  char cG;

  count = 0; fluxcount=0; 
  NST = pow(6,VOL2);
  for(i=0;i<VOL2;i++) str[i]='0';
  str[VOL2]='\0';
  
  //for(n=0;n<NST;n++){
  n=0;
  while(n<NST){
     
    //if((n%2176782336)==0) printf("%lld\n",n);
    //if((n%1679616)==0) printf("%lld\n",n);
    for(l=0;l<VOL;l++) { linkset[l]=0; isSET[2*l]=false; isSET[2*l+1]=false; }
    flagGI=1; flagJUMP=0;
    //printf("state=%s\n",str);
    for(l=0;l<VOL2;l++){
      cG=str[l];
      q=chk2lin[l]; 
      xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
      yf=next[DIM+2][q]; yb=next[DIM-2][q];
      // convert to 1d co-ordinates on the lattice configuration
      q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
      // from Fig 1 (left to right) in https://arxiv.org/pdf/1311.2459.pdf
      switch(cG){
        case '0' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=true;  conf[qy]=false; 
                  isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '1' : conf[q0x]=true;  conf[q0y]=true;  conf[qx]=true;  conf[qy]=true;  
                  isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '2' : conf[q0x]=false; conf[q0y]=false; conf[qx]=false; conf[qy]=false;
                  isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '3' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=true;  conf[qy]=false;
                  isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '4' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=false; conf[qy]=true; 
                  isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '5' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=false; conf[qy]=true; 
                  isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
      } //close switch
      // check Gauss' Law at this point on the odd sublattices
      if(linkset[xf]==4){ 
         flagGI = checkGL(conf,isSET,xf);
         if(flagGI==0){  addATpos(str,l);
           //printf("jumping at pos %d\n",l); 
           n=n+pow(6,VOL2-1-l); flagJUMP=1; break; }
         //if(flagGI==0) break;
      }
      if(linkset[xb]==4){ 
         flagGI = checkGL(conf,isSET,xb);
         if(flagGI==0){ addATpos(str,l);
           n=n+pow(6,VOL2-1-l); flagJUMP=1; break; }
      }
      if(linkset[yf]==4){ 
         flagGI = checkGL(conf,isSET,yf);
         if(flagGI==0){ addATpos(str,l);
           n=n+pow(6,VOL2-1-l); flagJUMP=1; break; } 
      }
      if(linkset[yb]==4){ 
         flagGI = checkGL(conf,isSET,yb);
         if(flagGI==0){ addATpos(str,l);
           n=n+pow(6,VOL2-1-l); flagJUMP=1; break; }
      }
    }//close l
    if(flagGI){
      for(p=VOL2;p<VOL;p++){
         pp=chk2lin[p]; 
         flagGI = checkGL(conf,isSET,pp);
         if(flagGI==0) break;
      }
    }// if flagGI closure 
     
    /* check if all links have been set for gauge invariant configs */
    //if(flagGI){ 
    //    for(l=0;l<VOL;l++){
    //        if(linkset[l]!=4) printf("unsatisfied links for site=%d\n",l);
    //        e1=(conf[2*l])? 1:-1; e2=(conf[2*l+1])? 1:-1;
    //        e3=(conf[2*next[DIM-1][l]])? 1:-1; e4=(conf[2*next[DIM-2][l]+1])? 1:-1;
    //        if(e1+e2-e3-e4) printf("GI violated at site=%d\n",l);
    //    }}
    count += flagGI;
    if(flagJUMP==0){ add(str); n++; }
  }
  printf("No of gauge invariant states: %d\n",count);
}

int checkGL(bool *conf, bool *isSET, int site){
  int e1,e2,e3,e4;
  int pp,qx,qy;
  qx=2*next[DIM-1][site]; 
  qy=2*next[DIM-2][site]+1; 
  pp=2*site;
  // safety check; remove later
  if((isSET[pp]==false)||(isSET[pp+1]==false)||(isSET[qx]==false)||(isSET[qy]==false)){
    printf("Links unset message \n");
  }
  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
  // if Q=e1+e2-e3-e4=0 return 1 (Gauss' Law OK) else return 0
  if(e1+e2-e3-e4) return 0;
  else return 1;
}

// Function to convert a given decimal number to a base 'base'
void fromDeci(long long int inputNum, int base, char* str){
    int index;
    index=VOL2;  // state looks reversed to fit neighbor indexing
    // (number)_10 --> (number)_(base) by repeatedly dividing it 
    // by base and taking remainder
    str[VOL2] = '\0';
    while (inputNum > 0){
        str[--index] = (char)(inputNum%base + '0');
        inputNum /= base;
    }
    while(index > 0){
        str[--index] = (char)('0');
    }
}

void add(char *str){
  int k,carry;
  carry=0; k=VOL2-1;
  if(str[k]<'5') { str[k]=str[k]+1; return; }
  else { str[k]='0'; carry=1; } 
  while(k>0){
    k--;
    if((carry==1)&&(str[k]=='5')) str[k]='0';
    else if(str[k]<'5') { str[k]=str[k]+1; carry=0; return; }
  }
}

void addATpos(char *str, int pos){
  int k,carry;
  carry=0;k=pos;
  if(str[k]<'5') { str[k]=str[k]+1; return; }
  else { str[k]='0'; carry=1; } 
  while(k>0){
    k--;
    if((carry==1)&&(str[k]=='5')) str[k]='0';
    else if(str[k]<'5') { str[k]=str[k]+1; carry=0; return; }
  }
}

// Here are some notes about the implementation of the Gauss' Law. In what 
// I thought would help, I implemented the checking of the odd sites using
// binary arithmetic. I.E., if e1,e2,e3,and e4 are the four boolean values,
//                         |
//                        e2
//                         |
//                   --e3--o--e1--
//                         |
//                        e4
//                         |
// then one can construct a truth table, and from these see that the Gauss
// law of e1+e2-e3-e4=0 is implemented only when there are even number of
// true and false values and when (e1==e2)&&(e3==e4)&&(e1!=e2).
// It turns out that the integer arithmetic as implemented now is slightly faster.

