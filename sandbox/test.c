#include <stdio.h>

int roll_idx(int i, int nx, int inc){
    if (i<nx-inc) return i+inc; else return i-nx+inc;
}
int main(){
   int a[10]={0,1,2,3,4,5,6,7,8,9};
   int inc=3;
   for (int i=0; i<10;i++){
    printf("%d,%d\n",a[i],a[roll_idx(i,10,3)]); 
   };
}
