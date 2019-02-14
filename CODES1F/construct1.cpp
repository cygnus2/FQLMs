/* construct gauge inequivalent states of the U(1) model  on triangular lattice */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm> // for copy
#include<iterator>  // for ostream_iterator
#include<math.h>
#include<vector>
#include "define.h"
//#include <algorithm> //std::count
//#include <vector>    //std::vector
// For some notes see below


void conststates(){
  extern void add(std::vector<char>&);
  extern void setGL(bool*, bool*, char, int);
  std::vector<int>chkpts;
  std::vector<char>str;
  int count,fluxcount;
  int flagFLUX;
  int base=6;
  long long int n,NST;
  int i,j,nPOINT;
  int ind_lin,ind_chk;

  if(LX<LY){
   printf("Please switch the LX and the LY. \n");
   exit(0);
  }
  if((LX!=4)&&(LX!=6)){
   printf("Please use the other routine construct_old for fixing states. \n");
   exit(0);
  }

  /* count no of points where independent links can be assigned */
  nPOINT=0;
  for(j=0;j<LY;j=j+2){
  for(i=0;i<LX;i=i+2){
    ind_lin=j*LX+i;
    chkpts.push_back(ind_lin);
    nPOINT++;
  }}
  for(i=1;i<LY;i=i+2){
   ind_lin=i*LX+i;
   chkpts.push_back(ind_lin);
   nPOINT++;
  }
  printf("LX=%d, LY=%d, #-points=%d\n",LX,LY,nPOINT);
  for(i=0;i<nPOINT;i++){
   printf("(x,y)=(%d,%d)\n",chkpts[i]%LX,chkpts[i]/LX);
  }
  
  count = 0; fluxcount=0; 
  NST = pow(6,nPOINT);
  str.reserve(nPOINT);
  for(i=0;i<nPOINT;i++) str.push_back('0');
  str[nPOINT]='\0';
  std::cout << str.size() << std::endl;
  
  for(n=0;n<NST;n++){
    bool conf[2*VOL],isSET[2*VOL];
    int x,y,flagGI;
    int p,pp,l,cc;
    int x0,y0,q0,q0x,q0y,qx,qy;
    int Q,e1,e2,e3,e4;
    char cG;
    long long int k;
    
    /* print state */
    //for(int ii=0; ii<str.size(); ii++) std::cout << str[ii];
    // std::cout << std::endl; 
     
    /* set links using Gauss Law at vertices */
    for(l=0;l<nPOINT;l++){
      /* the current site, and the GL-type */
      q0=chkpts[l]; cG=str[l];
      setGL(conf,isSET,cG,q0);
      // next, recursively set GL on unsatisfied vertices
       
    }
    //count += flagGI;
    add(str);
  }
  printf("No of gauge invariant states: %d\n",count);
}

// Function to convert a given decimal number to a base 'base'
void fromDeci(long long int inputNum, int base, char* str){
    int index;
    index=0;  // state looks reversed to fit neighbor indexing
    // (number)_10 --> (number)_(base) by repeatedly dividing it 
    // by base and taking remainder
    while (inputNum > 0){
        str[index++] = (char)(inputNum%base + '0');
        inputNum /= base;
    }
    while(index < VOL2){
        str[index++] = (char)('0');
    }
    str[VOL2] = '\0';
}

void add(std::vector<char> &str){
  int k,carry;
  carry=0;
  if(str[0]<'5') { str[0]=str[0]+1; return; }
  else { str[0]='0'; carry=1; } 
  k=1;
  while(k<VOL2){
    if((carry==1)&&(str[k]=='5')) str[k]='0';
    else if(str[k]<'5') { str[k]=str[k]+1; carry=0; return; }
    k++;  
  }
}

void setGI(bool conf[2*VOL], bool isSET[2*VOL], char cG, int q0){
  int qx,qy;
  /* note that q0 has the point in linear lattice co-ordinates */
  qx=next[DIM-1][q0]; qy=next[DIM-2][q0];
  // convert to 1d co-ordinates on the lattice configuration
  q0x=2*q0; q0y=2*q0+1; qx=2*qx;  qy=2*qy+1;
  // from Fig 1 (left to right) in https://arxiv.org/pdf/1311.2459.pdf
  switch(cG){
     case '0' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=true;  conf[qy]=false; 
               isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;  break;
     case '1' : conf[q0x]=true;  conf[q0y]=true;  conf[qx]=true;  conf[qy]=true;  
               isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;  break;
     case '2' : conf[q0x]=false; conf[q0y]=false; conf[qx]=false; conf[qy]=false; 
               isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;  break;
     case '3' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=true;  conf[qy]=false; 
               isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;  break;
     case '4' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=false; conf[qy]=true;  
               isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;  break;
     case '5' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=false; conf[qy]=true;  
               isSET[q0x]=true; isSET[q0y]=true; isSET[qx]=true; isSET[qy]=true;  break;
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

