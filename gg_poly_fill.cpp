#include<iostream>
#include<graphics.h>
using namespace std;
int main(){
int i,j,n,k,x[20],y[20],ymin=10000,ymax=0,dy,dx,in_x[100],temp; 
float slope[100];
int window1 =initwindow(800,800);
cout<<"Enter the number of vertices"<<endl;
cin>>n;
cout<<"Enter the coordinates of edges"<<endl;
for(i=0;i<n;i++){
cin>>x[i]>>y[i];
if(y[i]>ymax)
ymax=y[i];
if(y[i]<ymin)
ymin=y[i];
}
x[n]=x[0];y[n]=y[0];
for(i=0;i<n;i++)
line(x[i],y[i],x[i+1],y[i+1]);
delay(4000);
for(i=0;i<n;i++){
dy=y[i+1]-y[i]; 
dx=x[i+1]-x[i]; 
if(dy==0)
slope[i]=1.0;
if(dx==0)
slope[i]=0.0;
if(dy!=0 && dx!=0)
slope[i]=(float)dx/dy;
}
for(i=ymin;i<=ymax;i++){
k=0;
for(j=0;j<n;j++){
if((y[j]<=i && y[j+1]>i) || (y[j]>i && y[j+1]<=i)){
in_x[k]=(int)(x[j]+ slope[j]*(i-y[j]));
k++;
}
}
for(int m=0;m<k-1;m++){
for(int l=0;l<k-1;l++){
if(in_x[l]>in_x[l+1]){
temp=in_x[l];
in_x[l]=in_x[l+1];
in_x[l+1]=temp;
}
}
}
setcolor(2);
for(int p=0;p<k;p+=2)
{
line(in_x[p],i,in_x[p+1],i);
delay(100);
}
}
system("pause");
return 1;
}

