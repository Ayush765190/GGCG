///DDA
#include <iostream>
#include <graphics.h>
#include <cmath>

void drawLineDDA(int x1, int y1, int x2, int y2) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps;

    if (abs(dx) > abs(dy)) {
        steps = abs(dx);
    } else {
        steps = abs(dy);
    }

    float incX = dx / (float)steps;
    float incY = dy / (float)steps;
    float x = x1, y = y1;

    for (int i = 0; i < steps; i++) {
        putpixel(round(x), round(y), WHITE);
        x += incX;
        y += incY;
    }
}

int main() {
    int gd = DETECT, gm;
    initgraph(&gd, &gm, "");

    int x1 = 100, y1 = 100; // Starting point
    int x2 = 400, y2 = 300; // Ending point

    drawLineDDA(x1, y1, x2, y2);

    getch();
    closegraph();
    return 0;
}


////mid point ellipse
#include<iostream>
#include<graphics.h>
#include<conio.h>
#include<math.h>
using namespace std;

void plotpoints(int x, int y,int *p){
    putpixel(x+p[0],getmaxy()-(y+p[1]),255);
    putpixel(x+p[0],getmaxy()-(-y+p[1]),255);
    putpixel(-x+p[0],getmaxy()-(-y+p[1]),255);
    putpixel(-x+p[0],getmaxy()-(y+p[1]),255);
}

void Ellipse(int a,int b, int *p){
	int x=0,y=b;
	int sa=a*a;
	int sb=b*b;
	double d1=sb-sa*b+0.25*sa;
	plotpoints(x,y,p);
	while( sa*(y-0.5) > sb*(x+1)){ // Region 1
		if(d1<0) { //choose E   E= b^2 (2x + 3)
			d1+=sb*((2*x)+3);
		}
		else{	//choose SE    SE= b^2 (2x + 3) + a^2 (-2y + 2)
			d1+=sb*((2*x)+3) + sa*(-(2*y)+2);
			y--;
		}
		x++;
		plotpoints(x,y,p);
	}
	double d2 = sb*(x+0.5)*(x+0.5) + sa*(y-1)*(y-1) -sa*sb;
	while (y>0){ // Region 2
		if(d2<0){ // choose SE  SE= b^2 (2x + 2) + a^2 (-2y + 3)
			d2+= sb*((2*x)+2) + sa*(-(2*y)+3);
			x++;
		}
		else { // choose S    S= a^2 (-2y + 3)
			d2+= sa*(-(2*y)+3);
		}
		y--;
		plotpoints(x,y,p);
	}
}
int main(){
    int gd = DETECT, gm;
	char pathtodriver[] = "";
	initgraph(&gd, &gm, pathtodriver);
    int *p=new int(2);

	
    int a =50;
    int b =60;

    p[0]=200;
    p[1]=300;

	Ellipse(a,b,p);
    getch();
    closegraph();
    return 0;
}

////bresenhams circle
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <graphics.h>

using namespace std;

void drawCirclePixels(int c_x, int c_y, int x, int y, int val)
{
  putpixel(c_x + x, c_y + y, val);
  putpixel(c_x + y, c_y + x, val);
  putpixel(c_x + y, c_y - x, val);
  putpixel(c_x + x, c_y - y, val);
  putpixel(c_x - x, c_y - y, val);
  putpixel(c_x - y, c_y - x, val);
  putpixel(c_x - y, c_y + x, val);
  putpixel(c_x - x, c_y + y, val);
  return;
}

void BresenhamCircle(int x1, int y1, int r)
{
  // Setup
  int gd = DETECT, gm;
  initgraph(&gd, &gm, "");

  // Get middle of window + given value as centre
  int x_c = round(x1 + getmaxx() / 2);
  int y_c = round(-y1 + getmaxy() / 2);

  // Initial value of d
  int d = round(5.0 / 4.0 - r);

  // Draw initial pixel
  drawCirclePixels(x_c, y_c, 0, -r, WHITE);

  // Output to terminal
  cout << "\nIst OCTANT\n-------------" << endl;
  cout << "\ni\tPixel\td\tPlotted Values" << endl;
  cout << "0\t  \t  \t(" << x1 << "," << y1 + r << ")" << endl;

  int i = 0;
  string pixel = "";
  int x = 0;
  int y = r;

  while (y >= x)
  {
    i = i + 1;
    int d_temp = d;

    // Choose E pixel
    if (d < 0)
    {
      d += 2 * x + 3;
      x += 1;
      pixel = "E";
    }
    // Choose SE pixel
    else
    {
      d += 2 * (x - y) + 5;
      x += 1;
      y -= 1;
      pixel = "SE";
    }
    drawCirclePixels(x_c, y_c, x, -y, WHITE);

    // Output to terminal
    cout << i << "\t" << pixel << "\t" << d_temp << "\t(" << x << "," << y << ")" << endl;
  }
  cout << endl;

  // Clean up
  getch();
  closegraph();
}

int main(void)
{
  int x, y, r;
  cout << "Enter Centre (x y): ";
  cin >> x >> y;
  cout << "Enter Radius (r): ";
  cin >> r;

  BresenhamCircle(x, y, r);

  return 0;
}
