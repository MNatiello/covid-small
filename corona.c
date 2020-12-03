/*
 * Compile: Use Makefile
 * Version 11  2020-09-27
 */
#include "corona.h"

void GetPars(void){
  unsigned int j;
  double p,aux=0.;
  double auxU=0.;

 /*       vir and pr are given here
           These parameters are read in inputscript
           betaT,betaU                            scale of contagion ability
           detT, detU                               scale of detection/remotion ability 
  */

/* viremic level, stage j */
   vir[0]= 0.0; 
   vir[1]= 0.04061609517588571;
   vir[2]= 0.1201001729227572;
   vir[3]= 0.1695124040711986;
   vir[4]= 0.1712557902092357;
   vir[5]= 0.1413344669889077;
   vir[6]= 0.1018975776295615;
   vir[7]= 0.06670354441952196;
   vir[8]= 0.04061609517588571;

/* promotion rate to next stage */
  for (j=0;j<CONT;j++) {pr[j]=1.;};
  pr[1]=1./2.10;  pr[2]=1./1.86;  /* special case for stages 1 and 2 */

  /* remotion for T */
  if((delayT<1) || (delayT>(CONT-1))) {fprintf(stderr, "Wrong delay  %d \n",delayT); exit(4);};
    
  for (j=0;j<delayT;j++) {
	p=detT*vir[0]/(pr[j]+detT*vir[0]);
        aux = (1.-p)*aux+p;
        remT[j] = detT*vir[0];
        };
 for (j=delayT;j<CONT;j++) {
	p=detT*vir[j]/(pr[j]+detT*vir[j]);
        aux = (1.-p)*aux+p;
        remT[j] = detT*vir[j];
        };
 fprintf(stderr,"remotion rateT %f Prob %f Acum %f %d\n",remT[j-1],p,aux,j);
 eps=1.-aux;                                                                                                   /* aux are the removed, eps the undetected */

  /* remotion for U */
  if((delayU<1) || (delayU>(CONT-1))) {fprintf(stderr, "Wrong delay %d \n",delayU); exit(4);};
    
  for (j=0;j<delayU;j++) {
	p=detU*vir[0]/(pr[j]+detU*vir[0]);
        auxU = (1.-p)*auxU+p;
        remU[j] = detU*vir[0];
        };
 for (j=delayU;j<CONT;j++) {
	p=detU*vir[j]/(pr[j]+detU*vir[j]);
        auxU = (1.-p)*auxU+p;
        remU[j] = detU*vir[j];
        };
  fprintf(stderr,"remotion rateU %f Prob %f Acum %f %d\n",remU[j-1],p,auxU,j);
  slow=auxU;                                                                                                       /* aux are the removed, slow are the detected */
    
  return ;
}

void GetRates(unsigned int N){
  unsigned int j,count;
  count=0;
  W[0]=ext*((double)X[0])/((double)(N+1));                                                                 /* external contagion. N here is all but one */
  count++;
  for(j=0;j<CONT;j++) {                                                                                               /* Infection by a detected case */
    W[count+j]=betaT*vir[j]*((double)X[1+j])*((double)X[0])/((double)N);
  };
  count+=CONT;
  for(j=0;j<CONT;j++) {                                                                                               /* Infection by a non detected case */
    W[count+j]=betaU*vir[j]*((double)X[CONT+1+j])*((double)X[0])/((double)N);
  };

  count+=CONT;
   for(j=0;j<CONT;j++) {                                                                                              /* Promotion of infectedT to next infective stage */
     W[count+j]=pr[j]*X[1+j];
  };

  count+=CONT;
   for(j=0;j<CONT;j++) {                                                                                              /* Promotion of infectedU to next infective stage */
     W[count+j]=pr[j]*X[CONT+1+j];
  };
  count+=CONT;
  for(j=0;j<CONT;j++) {                                                                                              /* Remotion of detected infectedT */
     W[count+j]=remT[j]*X[1+j];
  };
 count+=CONT;
    for(j=0;j<CONT;j++) {                                                                                              /* Remotion of detected infectedU */
     W[count+j]=remU[j]*X[CONT+1+j];
  };

  return ;
 }

unsigned int BinS(double x){
  unsigned int i, first, last, middle;
  double ratesum[Events], aux;
  aux=0;
  for (i=0;i<Events;i++){
    aux+=W[i];
    ratesum[i]=aux;
  }
  first=0; last=Events-1; 
  while (first < last) {
    middle=(first+last)/2;
    if (ratesum[middle] < x) first = middle+1;    
    else last=middle;
  };
  return first;   
}

void UpdatePops(unsigned int ev){
  double rn;
  if (ev == (unsigned int)0) {
      rn=(double) (dsfmt_genrand_close_open(&dsfmt));
      X[0]--;
      if(rn>slow) {X[1+CONT]++;} else {X[1]++;};                        /* external contagion */
      }
  else if (ev <= (unsigned int)CONT) {
     rn=(double) (dsfmt_genrand_close_open(&dsfmt));
     X[0]--;
     if (rn>eps) {X[1]++;} else {X[1+CONT]++;};          
     /*fprintf(stderr,"D-infected\n");*/
    }
  else if (ev <= (unsigned int)(2*CONT)) {
     rn=(double) (dsfmt_genrand_close_open(&dsfmt));
     X[0]--;
     if (rn>slow) {X[1+CONT]++;} else {X[1]++;};
     /*fprintf(stderr,"N-infected\n");*/
  }
  else if (ev < (unsigned int)(3*CONT))    {X[ev-2*CONT]--; X[ev-2*CONT+1]++;}
  else if (ev == (unsigned int)(3*CONT)) {X[ev-2*CONT]--; X[2*CONT+1]++;} /*these should be recovered */
  else if (ev < (unsigned int)(4*CONT))    {X[ev-2*CONT]--; X[ev-2*CONT+1]++;}
  else if (ev == (unsigned int)(4*CONT)) {X[ev-2*CONT]--; X[2*CONT+2]++;}/*these should be recovered*/
  else if (ev <= (unsigned int)(5*CONT)) {X[ev-4*CONT]--; Cases++ ;} //update cases
  else if (ev <= (unsigned int)(6*CONT)) {X[ev-4*CONT]--; Cases++ ;} //update cases
    else {fprintf(stderr, "Wrong event index %d  \n",ev); exit(1);};
  return ;
   }    

void main(int argc, char *argv[])
{ 
  unsigned int N,Nrates,day;
  long int count;
  unsigned int Duration,Stat,event;
  unsigned int S,T,U;
  unsigned int i,j,k;
  double dt,dd;
  double R;
  char fn[20], name[20],number[4];
  FILE *file,*Ifile,*Ofile;

  /* Reading section */
  
if(argc<16){
 fprintf(stderr,"Usage: name, Stat, Duration,ext,betaT,betaU,detT,detU,delayT,delayU,S,T,U,seed,wipe,Tmin; \n");
 exit(2);
}

Ifile=fopen("IConditions","r"); // busca condiciones inciales
Ofile=fopen("FConditions","w");


 sscanf(argv[1],"%s",&name);       fprintf(stderr,"Name: %s\n",name); fflush(stderr);
 sscanf(argv[2],"%d",&Stat);          fprintf(stderr,"Stat=%i\n",Stat); fflush(stderr);
 sscanf(argv[3],"%d",&Duration);   fprintf(stderr,"Duration=%i\n",Duration); fflush(stderr);
 sscanf(argv[4],"%lf",&ext);         fprintf(stderr,"external contagion rate =%lf\n",ext); fflush(stderr);
 sscanf(argv[5],"%lf",&betaT);         fprintf(stderr,"contagion scale T =%lf\n",betaT); fflush(stderr);
 sscanf(argv[6],"%lf",&betaU);         fprintf(stderr,"contagion scale U=%lf\n",betaU); fflush(stderr);
 sscanf(argv[7],"%lf",&detT);         fprintf(stderr,"remotion scale T =%lf\n",detT); fflush(stderr);
 sscanf(argv[8],"%lf",&detU);         fprintf(stderr,"remotion scale U =%lf\n",detU); fflush(stderr);
 sscanf(argv[9],"%d",&delayT);         fprintf(stderr,"detection delay T =%d\n",delayT); fflush(stderr);
 sscanf(argv[10],"%d",&delayU);         fprintf(stderr,"detection delay  U =%d\n",delayU); fflush(stderr);
 sscanf(argv[11],"%d",&S);              fprintf(stderr,"Susceptibles =%i\n",S); fflush(stderr);
 sscanf(argv[12],"%d",&T);              fprintf(stderr,"Initial T=%i\n",T); fflush(stderr);
 sscanf(argv[13],"%d",&U);           fprintf(stderr,"Initial U =%i\n",U); fflush(stderr);
 sscanf(argv[14],"%d",&seed);       fprintf(stderr,"random seed =%i\n",seed); fflush(stderr);
 sscanf(argv[15],"%d",&wipe);       fprintf(stderr,"If non-zero delete short-lived epidemies: wipe=%i\n",wipe); fflush(stderr);
 sscanf(argv[16],"%d",&Tmin);       fprintf(stderr,"Proposed cutoff time, Tmin=%i\n",Tmin); fflush(stderr); 

 if((wipe>0) && (Tmin >= Duration)){
	Tmin= Duration-1;
	fprintf(stderr,"Tmin too large, new Tmin %d\n", Tmin);	 }

 if(Stat>MaxStat){ 
   fprintf(stderr,"Stat: %d reset to %d \n",Stat,MaxStat);
   Stat=MaxStat;
 }
 
 if((Duration>Maxdays) || (Duration<1)){
   fprintf(stderr,"Improper duration: %d  Maxdays: %d Retry \n",Duration,Maxdays);
   exit(3);
 }

/* initialise random generator */
 dsfmt_init_gen_rand(&dsfmt,seed);
  
 GetPars();
 
  /* Repetitions loop */
 for (i=0; i<Stat; i++){

 for(j=0;j<Events;j++){W[j]=0.;};
 for(j=0;j<SIZE;j++){X[j]=0;};
 
Nrates=S+T+U-1;
if(NULL == Ifile){
 /*load initial condition*/
 X[0]=S;
 X[1]=T;
 X[CONT+1]=U;
 Cases = 0; // inicializa Casos por corrida
 dt=0.;
 day=1;
 fprintf(stderr,"Corrida sin archivo IConditions\n");
}
else
{
 fprintf(stderr,"Corrida CON archivo IConditions\n");
        wipe = -1; //turn off wipping
	fscanf(Ifile,"%u",&day);
  for (j=0;j<SIZE;j++)
	  fscanf(Ifile,"%u ",X+j);
  fscanf(Ifile,"%u %lf",&Cases,&dt);
  if ( dt>=(double)day ) {dt=(double) (day-1) + (1./86400);} /* correccion por si dt excede el día */
};
   
 /*Run cycle */
 count=0;
  sprintf(number,"%d",i);
  strcpy(fn,name);
  strcat(fn,number);
  file=fopen(fn,"a");

 while (dt <= Duration){

  /* Compute number of infected. Starts in j=1 */
    N=0;
    for (j=1;j<SIZE-2;j++) {N+=X[j];};
    
    /* Check if we keep this simulation */
    if( (N<1) && (wipe>0) && (day<Tmin) ) {
      i--;
      fclose(file);
      file=fopen(fn,"w"); //close wiped file and reopen
      fprintf(stderr,"Wiped %s %i < %i Tmin. wipe= %d\n",fn,day,Tmin,wipe);
      dt=(double)Duration+(double)(1./86400.);  /*dt=Duration+1;*/
    }
    else {
    GetRates(Nrates);
    R=0.;
    for (j=0;j<Events;j++){R+=W[j];};
    if(R>(double)0) {
             count++;
             dd=(double) (R*dsfmt_genrand_close_open(&dsfmt));
             event=BinS(dd); 
             dd=(-log((double)(dsfmt_genrand_close_open(&dsfmt))))/R;         /* Time for present event */
             while ((dd>1) && (day<=Duration)) {
	       dt++;
	       dd--;
	       fprintf(file,"%d ",day); 
               for(j=0;j<SIZE;j++){fprintf(file,"%d ",X[j]);};
               fprintf(file,"%d ",Cases);  // Cases recorded
	       fprintf(file,"\n");
               fflush(file);
	       day++;
	     };
             dt+= dd;                                                          /* Elapsed time */
	     UpdatePops(event);
             if ((dt>day) && (day<=Duration)) {
	        fprintf(file,"%d ",day); 
                for(j=0;j<SIZE;j++){fprintf(file,"%d ",X[j]);};
                fprintf(file,"%d ",Cases);  // Cases recorded
		fprintf(file,"\n");
                fflush(file);
	       day++;
	     };
    }
    else{
		fprintf(stderr,"End of realisation because of nonpositive total rate: R= %f at day %d \n",R,day);
		for (k=day;k<=Duration;k++){
                     fprintf(file,"%d ",k); 
		     for(j=0;j<SIZE;j++){fprintf(file,"%d ",X[j]);};
		     fprintf(file,"%d ",Cases); 
		     fprintf(file,"\n");
		     fflush(file);
		     day++;
		};
		dt=(double)Duration+(double)(1./86400.); /*               dt=Duration+1; */
           }
    };
  };
  // output to Ofile
  //
if(day > Duration) // do not save wipped ends
{
fprintf(Ofile,"%u ",day);
  for (k=0;k<SIZE;k++) fprintf(Ofile,"%u ",X[k]);
  fprintf(Ofile,"%u %lf\n",Cases,dt);
};

  fprintf(stderr, "Events: %lu \n",count);	
  fclose(file);
 };
 if(NULL != Ifile){fflush(Ifile); fclose(Ifile);
   fprintf(stderr, "file %s closed\n","IConditions");
 };
 fflush(Ofile); fclose(Ofile);
 fprintf(stderr, "file %s closed \n","FConditions");
 fprintf(stderr, "exiting without errors...\n");
 exit(0);
};
