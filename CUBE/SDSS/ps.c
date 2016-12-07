#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

fftw_complex *DeltaKi;
int Nd,Nd21;
float boxsize;
double *K3d;

void cal_numps(int nst,char filename[256])
{
  int i,j,k,id;
  float kk,volume,hre,him;
  float pi2L,Nor,tt1,tt2;
  double pssim[10000][3];
  int  nps[10000],pid;
  FILE *fp;

  printf("Now output the power spectrum and Hamilton forces for %d initial condition\n",nst);
    
  for(i=0;i<Nd2;i++)
  {
  	pssim[i][0]=0.0;
  	pssim[i][1]=0.0;
  	nps[i]=0;
  }

  pi2L=2*PI/boxsize;  
  volume=pow(boxsize,3);
  for(i=0;i<Nd;i++)
  {
     for(j=0;j<Nd;j++)
     {
     	for(k=0;k<Nd21;k++)
     	{ 
     	   id=(i*Nd+j)*Nd21+k;
     	   kk=sqrt(K3d[i]*K3d[i]+K3d[j]*K3d[j]+K3d[k]*K3d[k]);
     	   
     	   pid=(int)(kk/pi2L);
     	   if(pid < Nd2)
     	   {
     	   	tt1=DeltaKi[id][0]*DeltaKi[id][0]+DeltaKi[id][1]*DeltaKi[id][1];

     	   	if(k == 0)
     	   	{
     	   		pssim[pid][0]=pssim[pid][0]+kk;
     	   		pssim[pid][1]=pssim[pid][1]+tt1;
     	   		nps[pid]++;
     	   	}else
     	   	{
     	   		pssim[pid][0]=pssim[pid][0]+2*kk;
     	   		pssim[pid][1]=pssim[pid][1]+2*tt1;
     	   		nps[pid]=nps[pid]+2;
     	   	}
     	   }
     	}
     }
  }

  sprintf(filename,"%s/powspe/PSHF_%.4d.dat",outdir,nst);
  
  if(!(fp = fopen(filename, "w")))
  {
    printf("can't open file `%s'\n", filename);
    exit(1);
  }

  for(i=1;i<Nd2;i++)
  {
  	pssim[i][0]=pssim[i][0]/nps[i];
  	pssim[i][1]=pssim[i][1]/nps[i];
//  	tt1=PowerSpe_Z0(pssim[i][0]);
  	
  	fprintf(fp,"%e %e %e %d\n",pssim[i][0],pssim[i][1]*volume,tt1,nps[i]);
  }  
  fclose(fp);

  printf("output power spectrum and Hamilton forces into file %s done\n", filename);
}

void read_deltaki(char filename[256])
{
   FILE *fp;
   int nst,i,j,ccid;
   double ipp;
   
   printf("Now read DeltaKi data from previously dumped file:\n%s\n", filename);
   
   DeltaKi=(fftw_complex *)malloc(sizeof(fftw_complex)*Nd*Nd*Nd21);
   
   if(!(fp = fopen(filename, "r")))
   {
    	printf("can't open file `%s'\n", filename);
    	exit(1);
   } 
   fread(&nst,sizeof(int),1,fp);
   fread(&ipp,sizeof(double),1,fp);
   fread(DeltaKi,sizeof(fftw_complex), Nd*Nd*Nd21, fp);
   fclose(fp);
   
   printf("Ip is %e, nstep is %d\n",ipp, nst);

}

void Initial_GFKXYZ()
{
   int i;
   double pi2L;

   K3d=(double *)malloc(sizeof(double)*Nd);
   pi2L=2*PI/boxsize;
   for(i=0;i<Nd;i++)
   {
      K3d[i]=i*pi2L;
      if(i > Nd2)
      {
     	K3d[i]=(i-Nd)*pi2L;
      }   	
   }
}


int main()
{
   int nst;
   char filename[256];

   Nd=500;
   Nd21=Nd/2+1;
   boxsize=500; // Mpc/h


   sprintf(filename,"???Here input the filename of the deltaK???");
   read_deltaki(filename);
   
   Initial_GFKXYZ();
   
   nst=1500;
   sprintf(filename,"??the output power spectra??");
   cal_numps(nst, filename);

   return 1;
}
