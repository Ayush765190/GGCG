#include<iostream>
#include<graphics.h>
#include<conio.h>
#include<math.h>
using namespace std;



void plotpoints(int x, int y,int *p){
    putpixel(x+p[0],getmaxy()-(y+p[1]),255);
    putpixel(y+p[0],getmaxy()-(x+p[1]),255);
    putpixel(y+p[0],getmaxy()-(-x+p[1]),255);
    putpixel(x+p[0],getmaxy()-(-y+p[1]),255);
    putpixel(-x+p[0],getmaxy()-(-y+p[1]),255);
    putpixel(-y+p[0],getmaxy()-(-x+p[1]),255);
    putpixel(-y+p[0],getmaxy()-(x+p[1]),255);
    putpixel(-x+p[0],getmaxy()-(y+p[1]),255);
    cout<<x<<", "<<-y<<endl;
}
void circle(int *p,int r){
    int x=0,y=r;
    int d=1-r;
    plotpoints(x,y,p);
    while(x<y){
        if(d<=0){
            d=d+2*x+3;
        }
        else{
            d=d+2*(x-y)+5;
            y--;
        }
        x++;
        plotpoints(x,y,p);
    }
}
int main(){
    int gd = DETECT, gm;
	char pathtodriver[] = "";
	initgraph(&gd, &gm, pathtodriver);
    int *p=new int(2);
    int r=50;
    p[0]=100;
    p[1]=100;
    circle(p,r);

   
    getch();
    closegraph();
    return 0;
}
