/* C functions for the test particle module */
/* Incase you forget how to compile the shared lib file
 * gcc -o functions.so -shared -fPIC functions.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void method(float A,float *B){
    printf("Hello World\n");
    printf("A = %.3f and B = %.3f",A,B[0]);
    B[0] = A;

}

int roll_idx(int i, int nx, int inc){
    if (i<nx-inc) return i+inc; else return i-nx+inc;
}

int rolled_1d_idx(int x, int y, int z, int inc, int ax, int nx, int ny, int nz){
// index in 1D map of a 3D array after rolling by inc along axis ax
   int j,l;
         if (ax==0){ j=roll_idx(x,nx,inc); l=j*ny*nz+y*nz+nz;}
    else if (ax==1){ j=roll_idx(y,ny,inc); l=x*ny*nz+j*nz+nz;}
    else if (ax==2){ j=roll_idx(z,nz,inc); l=x*ny*nz+y*nz+j ;}
   return l;
}

double mean(double *A, int nx, int ny, int nz){
//double mean(const void *inA, int nx, int ny, int nz){
//   const double *A = (double *) inA;
   int i,lenA;
   double mean;
      mean=0.;
      lenA=nx*ny*nz; //sizeof(A)/sizeof(A[0]);
      for (i=0;i<lenA;i++){
         mean=mean+A[i];
      };
   return mean/lenA;
}

double rms(double *A, int nx, int ny, int nz){
   int i,lenA;
   double rms;
      rms=0.;
      lenA=nx*ny*nz; //sizeof(A)/sizeof(A[0]);
      for (i=0;i<lenA;i++){
         rms=rms+pow(A[i],2);
      };
   return rms/lenA;
}

double kurtosis(double *A, int nx, int ny, int nz){
   int i,lenA;
   double sq,qd;
      sq=0.; qd=0.;
      lenA=nx*ny*nz; //sizeof(A)/sizeof(A[0]);
      for (i=0;i<lenA;i++){
         sq=sq+pow(A[i],2);
         qd=qd+pow(A[i],4);
      };
   return qd/sq;
}
void c_sdk(double *A, int ax, int nx, int ny, int nz, double *ddx, double *ark){
   int i,x,y,z,k,l;
   double rmsv;
   double B[nx*ny*nz];
   for(i=0;i<nx/2;i++){
      ark[i]=0.;
      // Find B=A-np.roll(A,inc)
      for(x=0;x<nx;x++){
         for(y=0;y<ny;y++){
            for(z=0;z<nz;z++){
            k=x*ny*nz+y*nz+nz;
            l=rolled_1d_idx(x,y,z,i,ax,nx,ny,nz);
            B[k]=A[k]-A[l];
            };
         };
      };
      printf("Rolled A and found the increment");
      // Find rms value of B
      rmsv=rms(B,nx,ny,nz);
      printf("RMS Value of increment of A is %.4f",rmsv);
      // If rms(B) != 0, divide B by rms(B)
      if(rmsv!=0){
         for(x=0;x<nx*ny*nz;x++){
            B[x]=B[x]/rmsv;
         };
      };
      // Find kurtosis of B, assign to ark[i]
      printf("Found the kurtosis of increment");
      ark[i]=kurtosis(B,nx,ny,nz);
      // Assign ddx
      ddx[i]=(double) i;
   };
}

void c_correlation(double *A, double *B, double dx, int nlen, int nx, int ny, int nz, int ax, double *r, double *corr){
   int i,j,k,l,x,y,z;
   double meana, meanb;
   meana=mean(A,nx,ny,nz);
   meanb=mean(B,nx,ny,nz);
   printf("%.4f, %.4f \n \n",meana,meanb);
      for (i=0;i<nlen;i++){
         corr[i]=0;
         for (x=0;x<nx;x++){
            for (y=0;y<ny;y++){
               for (z=0;z<nz;z++){
                k=x*ny*nz+y*nz+nz;
                l=rolled_1d_idx(x,y,z,i,ax,nx,ny,nz);
                corr[i]=corr[i]+A[k]*B[l];
               };
            };
         };
         corr[i]=corr[i]/(nx*ny*nz*meana*meanb);
         r[i] = i*dx;
      };
}

////void c_correlation(double *A, double *B, double dx, int nlen, int nx, int ny, int nz, int ax, double *r, double *corr){
//void c_correlation(const void *inA, const void *inB, double dx, int nlen, int nx, int ny, int nz, int ax, void *inr, void *incorr){
//   int i,j,k,l,x,y,z;
//   const double *A = (double *) inA;
//   const double *B = (double *) inB;
//   double * r = (double *) inr;
//   double * corr = (double *) incorr;
//   double AA[nx][ny][nz], BB[nx][ny][nz];
//   for(x=0;x<nx;x++){
//      for(y=0;y<ny;y++){
//         for(z=0;z<nz;z++){
//            i=x*ny*nz+y*nz+z;
//            AA[x][y][z]=A[i];
//            BB[x][y][z]=B[i];
//         };
//      };
//   };
//   double meana, meanb;
//   meana=mean(AA,nx,ny,nz);
//   meanb=mean(BB,nx,ny,nz);
//   printf("%.4f, %.4f \n \n",meana,meanb);
//      for (i=0;i<nlen;i++){
//         corr[i]=0;
//         for (x=0;x<nx;x++){
//            for (y=0;y<ny;y++){
//               for (z=0;z<nz;z++){
//                       if (ax==0){ j=roll_idx(x,nx,i); k=y; l=z; }
//                  else if (ax==1){ j=x; k=roll_idx(y,ny,i); l=z; }
//                  else if (ax==2){ j=x; k=y; l=roll_idx(z,nz,i); }
//                corr[i]=corr[i]+AA[x][y][z]*BB[j][k][l];
//               };
//            };
//         };
//         corr[i]=corr[i]/(nx*ny*nz*meana*meanb);
//         r[i] = i*dx;
//      };
//}
