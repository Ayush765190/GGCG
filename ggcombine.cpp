// bresenhams line
#include <iostream>
#include <iomanip>
#include <cmath>
#include <graphics.h>

using namespace std;

void BresenhamLine(int x1, int y1, int x2, int y2) {
    // Setup
    int gd = DETECT, gm;
    initgraph(&gd, &gm, "");

    // Calculate dy, dx, a, b
    double dx = x2 - x1;
    double dy = y2 - y1;
    double a = 2 * dy;
    double b = -2 * dx;

    // Initial value of d
    double d = 2 * dy - dx;

    // Draw initial pixel
    putpixel(x1, y1, WHITE);

    // Output to terminal
    cout << "\ni\tPixel\td\tPlotted Values" << endl;
    cout << "0\t  \t  \t(" << round(x1) << "," << round(y1) << ")" << endl;

    double x = x1;
    double y = y1;

    string pixel = "";

    for (int i = 1;; i++) {
        double d_temp = d;

        // Choose NE pixel
        if (d >= 0) {
            d = d + a + b;
            x = x + 1;
            y = y + 1;
            pixel = "NE";
        }
        // Choose E pixel
        else {
            d = d + a;
            x = x + 1;
            pixel = "E";
        }
        // Exit condition
        if (x > x2 || y > y2)
            break;

        // draw pixel
        putpixel(x, y, WHITE);

        // Output to terminal
        cout << i << "\t" << pixel << "\t" << d_temp << "\t(" << round(x) << "," << round(y) << ")" << endl;
    }

    // Clean up
    getch();
    closegraph();

    cout << endl;
}

int main() {
    int x1 = -200;
    int y1 = -100;
    int x2 = 200;
    int y2 = 180;
    BresenhamLine(x1, y1, x2, y2);

    return 0;
}





//mid point circle
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


//line clip cohen sutherland
#include <iostream>
#include <graphics.h>
using namespace std;

// Function to compute region code for a point
int computeCode(int x, int y, int xmin, int xmax, int ymin, int ymax) {
    // Define region codes
    const int INSIDE = 0; // 0000
    const int LEFT = 1;   // 0001
    const int RIGHT = 2;  // 0010
    const int BOTTOM = 4; // 0100
    const int TOP = 8;    // 1000

    int code = INSIDE;
    if (x < xmin) // to the left of rectangle
        code |= LEFT;
    else if (x > xmax) // to the right of rectangle
        code |= RIGHT;
    if (y < ymin) // below the rectangle
        code |= BOTTOM;
    else if (y > ymax) // above the rectangle
        code |= TOP;
    return code;
}

// Function to clip a line using Cohen-Sutherland algorithm
void cohenSutherland(int x1, int y1, int x2, int y2, int xmin, int xmax, int ymin, int ymax) {
    int code1 = computeCode(x1, y1, xmin, xmax, ymin, ymax);
    int code2 = computeCode(x2, y2, xmin, xmax, ymin, ymax);
    bool accept = false;

    while (true) {
        if (!(code1 | code2)) { // Both endpoints lie inside the window
            accept = true;
            break;
        } else if (code1 & code2) { // Both endpoints lie outside the window
            break;
        } else {
            int codeOut = code1 ? code1 : code2;

            float x, y;
            if (codeOut & 8) {           // point is above the clip window
                x = x1 + (x2 - x1) * (ymax - y1) / (float)(y2 - y1);
                y = ymax;
            } else if (codeOut & 4) { // point is below the clip window
                x = x1 + (x2 - x1) * (ymin - y1) / (float)(y2 - y1);
                y = ymin;
            } else if (codeOut & 2) {  // point is to the right of clip window
                y = y1 + (y2 - y1) * (xmax - x1) / (float)(x2 - x1);
                x = xmax;
            } else if (codeOut & 1) {   // point is to the left of clip window
                y = y1 + (y2 - y1) * (xmin - x1) / (float)(x2 - x1);
                x = xmin;
            }

            if (codeOut == code1) {
                x1 = x;
                y1 = y;
                code1 = computeCode(x1, y1, xmin, xmax, ymin, ymax);
            } else {
                x2 = x;
                y2 = y;
                code2 = computeCode(x2, y2, xmin, xmax, ymin, ymax);
            }
        }
    }

    if (accept) {
        // Draw the clipped line
        line(x1, y1, x2, y2);
    } else {
        cout << "Line lies completely outside the clip window." << endl;
    }
}

int main() {
    int gd = DETECT, gm;
    initgraph(&gd, &gm, "");

    // Input the coordinates of the clipping window
    int xmin, xmax, ymin, ymax;
    cout << "Enter coordinates of the clipping window (xmin ymin xmax ymax): ";
    cin >> xmin >> ymin >> xmax >> ymax;

    // Draw the clipping window
    rectangle(xmin, ymin, xmax, ymax);

    // Input the coordinates of the line
    int x1, y1, x2, y2;
    cout << "Enter coordinates of the line (x1 y1 x2 y2): ";
    cin >> x1 >> y1 >> x2 >> y2;

    // Clip and draw the line
    cohenSutherland(x1, y1, x2, y2, xmin, xmax, ymin, ymax);

    delay(5000);
    closegraph();
    return 0;
}


///clip polygon sutherland hodgeman
#include<iostream>
#include <graphics.h>
using namespace std;
  
const int MAX_POINTS = 20;
  
// Returns x-value of point of intersection of two
// lines
int x_intersect(int x1, int y1, int x2, int y2,
                int x3, int y3, int x4, int y4)
{
    int num = (x1*y2 - y1*x2) * (x3-x4) -
              (x1-x2) * (x3*y4 - y3*x4);
    int den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4);
    return num/den;
}
  
// Returns y-value of point of intersection of
// two lines
int y_intersect(int x1, int y1, int x2, int y2,
                int x3, int y3, int x4, int y4)
{
    int num = (x1*y2 - y1*x2) * (y3-y4) -
              (y1-y2) * (x3*y4 - y3*x4);
    int den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4);
    return num/den;
}
  
// This functions clips all the edges w.r.t one clip
// edge of clipping area
void clip(int poly_points[][2], int &poly_size,
          int x1, int y1, int x2, int y2)
{
    int new_points[MAX_POINTS][2], new_poly_size = 0;
  
    // (ix,iy),(kx,ky) are the co-ordinate values of
    // the points
    for (int i = 0; i < poly_size; i++)
    {
        // i and k form a line in polygon
        int k = (i+1) % poly_size;
        int ix = poly_points[i][0], iy = poly_points[i][1];
        int kx = poly_points[k][0], ky = poly_points[k][1];
  
        // Calculating position of first point
        // w.r.t. clipper line
        int i_pos = (x2-x1) * (iy-y1) - (y2-y1) * (ix-x1);
  
        // Calculating position of second point
        // w.r.t. clipper line
        int k_pos = (x2-x1) * (ky-y1) - (y2-y1) * (kx-x1);
  
        // Case 1 : When both points are inside
        if (i_pos < 0  && k_pos < 0)
        {
            //Only second point is added
            new_points[new_poly_size][0] = kx;
            new_points[new_poly_size][1] = ky;
            new_poly_size++;
        }
  
        // Case 2: When only first point is outside
        else if (i_pos >= 0  && k_pos < 0)
        {
            // Point of intersection with edge
            // and the second point is added
            new_points[new_poly_size][0] = x_intersect(x1,
                              y1, x2, y2, ix, iy, kx, ky);
            new_points[new_poly_size][1] = y_intersect(x1,
                              y1, x2, y2, ix, iy, kx, ky);
            new_poly_size++;
  
            new_points[new_poly_size][0] = kx;
            new_points[new_poly_size][1] = ky;
            new_poly_size++;
        }
  
        // Case 3: When only second point is outside
        else if (i_pos < 0  && k_pos >= 0)
        {
            //Only point of intersection with edge is added
            new_points[new_poly_size][0] = x_intersect(x1,
                              y1, x2, y2, ix, iy, kx, ky);
            new_points[new_poly_size][1] = y_intersect(x1,
                              y1, x2, y2, ix, iy, kx, ky);
            new_poly_size++;
        }
  
        // Case 4: When both points are outside
        else
        {
            //No points are added
        }
    }
  
    // Copying new points into original array
    // and changing the no. of vertices
    poly_size = new_poly_size;
    for (int i = 0; i < poly_size; i++)
    {
        poly_points[i][0] = new_points[i][0];
        poly_points[i][1] = new_points[i][1];
    }
}
  
// Implements Sutherland–Hodgman algorithm
void suthHodgClip(int poly_points[][2], int poly_size,
                  int clipper_points[][2], int clipper_size)
{
    //i and k are two consecutive indexes
    for (int i=0; i<clipper_size; i++)
    {
        int k = (i+1) % clipper_size;
  
        // We pass the current array of vertices, it's size
        // and the end points of the selected clipper line
        clip(poly_points, poly_size, clipper_points[i][0],
             clipper_points[i][1], clipper_points[k][0],
             clipper_points[k][1]);
    }
  
    // Printing vertices of clipped polygon
    cout << "\nClipped Polygon : " << endl;
    for (int i=0; i < poly_size; i++)
        cout << '(' << poly_points[i][0] <<
                ", " << poly_points[i][1] << ") ";
    cout << endl << endl;

    // Drawing Clipped Polygon
    int poly_clipped[50];
    for (int q = 0; q < poly_size; q++)
    {
        for (int t = 0; t < 2; t++)
        {
            poly_clipped[q * 2 + t] = poly_points[q][t];
        }
    }
    setcolor(BLUE);
    poly_clipped[2 * poly_size] = poly_clipped[0];
    poly_clipped[2 * poly_size + 1] = poly_clipped[1];
    drawpoly(poly_size + 1, poly_clipped);

    getch();
}
  
//Driver code
int main()
{
    int gd = DETECT, gm, errorcode;
    initgraph(&gd, &gm, NULL);

    // Defining polygon vertices in clockwise order
    int poly_size = 3;
    int poly_points[20][2] = {{100,150}, {200,250},
                              {300,100}};
  
    // Defining clipper polygon vertices in clockwise order
    // 1st Example with square clipper
    int clipper_size = 4;
    int clipper_points[][2] = {{100,100}, {100,200},
                              {200,200}, {200,100} };
    setcolor(RED);
    rectangle(100, 100, 200, 200);

    setcolor(YELLOW);
    int poly[50];
    for (int q = 0; q < poly_size; q++)
    {
        for (int t = 0; t < 2; t++)
        {
            poly[q * 2 + t] = poly_points[q][t];
        }
    }
    poly[2 * poly_size] = poly[0];
    poly[2 * poly_size + 1] = poly[1];
    drawpoly(poly_size + 1, poly);
  
    //Calling the clipping function
    suthHodgClip(poly_points, poly_size, clipper_points,
                 clipper_size);

    getch();
  
    return 0;
}


////scan line polygon fill algo

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


////2d transformation
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <graphics.h>
#include <iostream>
#define COORD_SHIFT 100

using namespace std;

void clrscr()
{
#ifdef _WIN32
  system("cls");
#elif __unix__
  system("clear");
#endif
}

double **inputFigure(int n)
{
  cout << "Enter the matrix for the 2-D shape (homogeneous):\n";

  double **figure = NULL;
  figure = new double *[n];

  for (int i = 0; i < n; i++)
  {
    figure[i] = new double[3];
    for (int j = 0; j < 3; j++)
    {
      cin >> figure[i][j];
    }
  }

  return figure;
}

void drawFigure(double **points, int n)
{
  setcolor(WHITE);
  for (int i = 0; i < n; i++)
  {
    line(COORD_SHIFT + points[i][0],
         COORD_SHIFT + points[i][1],
         COORD_SHIFT + points[(i + 1) % n][0],
         COORD_SHIFT + points[(i + 1) % n][1]);
  }

  delay(5e3);
  cleardevice();
}

double **translate(double **figure, int dim, int m, int n)
{
  double **_figure = NULL;
  int T[dim][3] = {{1, 0, 0}, {0, 1, 0}, {m, n, 1}};

  _figure = new double *[dim];

  for (int i = 0; i < dim; i++)
  {
    _figure[i] = new double[3];
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < dim; k++)
      {
        _figure[i][j] += figure[i][k] * T[k][j];
      }
    }
  }

  return _figure;
}

double **rotate(double **figure, int dim, double theta)
{
  double **_figure = NULL;
  double T[dim][3] = {{cos(theta * M_PI / 180.0), sin(theta * M_PI / 180.0), 0},
                      {-sin(theta * M_PI / 180.0), cos(theta * M_PI / 180.0), 0},
                      {0, 0, 1}};

  _figure = new double *[dim];

  for (int i = 0; i < dim; i++)
  {
    _figure[i] = new double[3];
    for (int j = 0; j < 2; j++)
    {
      for (int k = 0; k < dim; k++)
      {
        _figure[i][j] += figure[i][k] * T[k][j];
      }
    }
  }

  return _figure;
}

double **scale(double **figure, int dim, int m, int n)
{
  double **_figure = NULL;
  int T[dim][3] = {{m, 0, 0}, {0, n, 0}, {0, 0, 1}};

  _figure = new double *[dim];

  for (int i = 0; i < dim; i++)
  {
    _figure[i] = new double[3];
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < dim; k++)
      {
        _figure[i][j] += figure[i][k] * T[k][j];
      }
    }
  }

  return _figure;
}

double **reflect(double **figure, int dim, int c)
{
  double **_figure = NULL;
  int T[dim][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  switch (c)
  {
  case 1:
    T[1][1] = -1;
    break;
  case 2:
    T[0][0] = -1;
    break;
  case 3:
    T[0][0] = 0;
    T[0][1] = 1;
    T[1][0] = 1;
    T[1][1] = 0;
    break;
  case 4:
    T[0][0] = -1;
    T[1][1] = -1;
    break;
  default:
    return NULL;
    break;
  }

  _figure = new double *[dim];

  for (int i = 0; i < dim; i++)
  {
    _figure[i] = new double[3];
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < dim; k++)
      {
        _figure[i][j] += figure[i][k] * T[k][j];
      }
    }
  }

  return _figure;
}

double **shear(double **figure, int dim, int m, int n)
{
  double **_figure = NULL;
  int T[dim][3] = {{1, n, 0}, {m, 1, 0}, {0, 0, 1}};

  _figure = new double *[dim];

  for (int i = 0; i < dim; i++)
  {
    _figure[i] = new double[3];
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < dim; k++)
      {
        _figure[i][j] += figure[i][k] * T[k][j];
      }
    }
  }

  return _figure;
}

void menu(double **figure, int dim)
{
  int ch = 0;
  double **_figure;

  do
  {
    clrscr();
    cout << "\nMenu\n-------\n(1) Translation\n(2) Rotation";
    cout << "\n(3) Scaling\n(4) Reflection\n(5) Shearing";
    cout << "\n(6) View Figure\n(7) Exit\n\nEnter Choice: ";
    cin >> ch;
    cout << endl;
    switch (ch)
    {
    case 1:
      int m, n;

      cout << "Enter translation in x-axis: ";
      cin >> m;
      cout << "Enter translation in y-axis: ";
      cin >> n;

      _figure = translate(figure, dim, m, n);

      cout << "Drawing Original Figure...\n";
      drawFigure(figure, dim);

      cout << "Drawing Transformed Figure...\n";
      drawFigure(_figure, dim);
      break;
    case 2:
      double theta;

      cout << "Enter rotation angle (degrees): ";
      cin >> theta;

      _figure = rotate(figure, dim, theta);

      cout << "Drawing Original Figure...\n";
      drawFigure(figure, dim);

      cout << "Drawing Transformed Figure...\n";
      drawFigure(_figure, dim);
      break;
    case 3:
      cout << "Enter scaling in x-axis: ";
      cin >> m;
      cout << "Enter scaling in y-axis: ";
      cin >> n;

      _figure = scale(figure, dim, m, n);

      cout << "Drawing Original Figure...\n";
      drawFigure(figure, dim);

      cout << "Drawing Transformed Figure...\n";
      drawFigure(_figure, dim);
      break;
    case 4:
      cout << "Reflect along\n(1) x-axis\n(2) y-axis\n(3) y = x\n(4) y = -x\n"
           << "\nEnter Choice: ";
      cin >> m;

      _figure = reflect(figure, dim, m);

      cout << "Drawing Original Figure...\n";
      drawFigure(figure, dim);

      cout << "Drawing Transformed Figure...\n";
      drawFigure(_figure, dim);
      break;
    case 5:
      cout << "Enter shearing in x-axis: ";
      cin >> m;
      cout << "Enter shearing in y-axis: ";
      cin >> n;

      _figure = shear(figure, dim, m, n);

      cout << "Drawing Original Figure...\n";
      drawFigure(figure, dim);

      cout << "Drawing Transformed Figure...\n";
      drawFigure(_figure, dim);
      break;
    case 6:
      cout << "Drawing Original Figure...\n";
      drawFigure(figure, dim);
      break;
    case 7:
    default:
      break;
    }

    delete _figure;

    cout << endl
         << "Finished..."
         << endl;

    if (ch != 7)
    {
      cout << "\nPress Enter to continue ...\n";
      cin.ignore();
      cin.get();
    }
  } while (ch != 7);
};

int main(void)
{
  int n;
  double **fig;
  int gd = DETECT, gm;

  initgraph(&gd, &gm, NULL);

  cout << "Enter number of points in the figure: ";
  cin >> n;

  fig = inputFigure(n);

  menu(fig, n);

  delete fig;
  closegraph();

  return 0;
}


////3d transformation
#include <iostream>
#include <direct.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <graphics.h>
#include <process.h>
using namespace std;

int gd = DETECT, gm;
double x1, x2, y2;

void draw_cube(double edge[20][3])
{
  double y1;
  initgraph(&gd, &gm, NULL);
  int i;
  clearviewport();
  for (i = 0; i < 19; i++)
  {
    x1 = edge[i][0] + edge[i][2] * (cos(2.3562));
    y1 = edge[i][1] - edge[i][2] * (sin(2.3562));
    x2 = edge[i + 1][0] + edge[i + 1][2] * (cos(2.3562));
    y2 = edge[i + 1][1] - edge[i + 1][2] * (sin(2.3562));
    line(x1 + 320, 240 - y1, x2 + 320, 240 - y2);
  }
  line(320, 240, 320, 25);
  line(320, 240, 550, 240);
  line(320, 240, 150, 410);
  getch();
  closegraph();
}

void scale(double edge[20][3])
{
  double a, b, c;
  int i;
  cout << "Enter The Scaling Factors: ";
  cin >> a >> b >> c;
  initgraph(&gd, &gm, NULL);
  clearviewport();
  for (i = 0; i < 20; i++)
  {
    edge[i][0] = edge[i][0] * a;
    edge[i][1] = edge[i][1] * b;
    edge[i][2] = edge[i][2] * c;
  }
  draw_cube(edge);
  closegraph();
}

void translate(double edge[20][3])
{
  int a, b, c;
  int i;
  cout << "Enter The Translation Factors: ";
  cin >> a >> b >> c;
  initgraph(&gd, &gm, NULL);
  clearviewport();
  for (i = 0; i < 20; i++)
  {
    edge[i][0] += a;
    edge[i][0] += b;
    edge[i][0] += c;
  }
  draw_cube(edge);
  closegraph();
}

void rotate(double edge[20][3])
{
  int ch;
  int i;
  double temp, theta, temp1;
  cout << "-=[ Rotation About ]=-" << endl;
  cout << "1:==> X-Axis " << endl;
  cout << "2:==> Y-Axis" << endl;
  cout << "3:==> Z-Axis " << endl;
  cout << "Enter Your Choice: ";
  cin >> ch;
  switch (ch)
  {
  case 1:
    cout << "Enter The Angle: ";
    cin >> theta;
    theta = (theta * 3.14) / 180;
    for (i = 0; i < 20; i++)
    {
      edge[i][0] = edge[i][0];
      temp = edge[i][1];
      temp1 = edge[i][2];
      edge[i][1] = temp * cos(theta) - temp1 * sin(theta);
      edge[i][2] = temp * sin(theta) + temp1 * cos(theta);
    }
    draw_cube(edge);
    break;

  case 2:
    cout << "Enter The Angle: ";
    cin >> theta;
    theta = (theta * 3.14) / 180;
    for (i = 0; i < 20; i++)
    {
      edge[i][1] = edge[i][1];
      temp = edge[i][0];
      temp1 = edge[i][2];
      edge[i][0] = temp * cos(theta) + temp1 * sin(theta);
      edge[i][2] = -temp * sin(theta) + temp1 * cos(theta);
    }
    draw_cube(edge);
    break;

  case 3:
    cout << "Enter The Angle: ";
    cin >> theta;
    theta = (theta * 3.14) / 180;
    for (i = 0; i < 20; i++)
    {
      edge[i][2] = edge[i][2];
      temp = edge[i][0];
      temp1 = edge[i][1];
      edge[i][0] = temp * cos(theta) - temp1 * sin(theta);
      edge[i][1] = temp * sin(theta) + temp1 * cos(theta);
    }
    draw_cube(edge);
    break;
  }
}

void reflect(double edge[20][3])
{
  int ch;
  int i;
  cout << "-=[ Reflection About ]=-" << endl;
  cout << "1:==> X-Axis" << endl;
  cout << "2:==> Y-Axis " << endl;
  cout << "3:==> Z-Axis " << endl;
  cout << "Enter Your Choice: ";
  cin >> ch;
  switch (ch)
  {
  case 1:
    for (i = 0; i < 20; i++)
    {
      edge[i][0] = edge[i][0];
      edge[i][1] = -edge[i][1];
      edge[i][2] = -edge[i][2];
    }
    draw_cube(edge);
    break;

  case 2:
    for (i = 0; i < 20; i++)
    {
      edge[i][1] = edge[i][1];
      edge[i][0] = -edge[i][0];
      edge[i][2] = -edge[i][2];
    }
    draw_cube(edge);
    break;

  case 3:
    for (i = 0; i < 20; i++)
    {
      edge[i][2] = edge[i][2];
      edge[i][0] = -edge[i][0];
      edge[i][1] = -edge[i][1];
    }
    draw_cube(edge);
    break;
  }
}

void perspect(double edge[20][3])
{
  int ch;
  int i;
  double p, q, r;
  cout << "-=[ Perspective Projection About ]=-" << endl;
  cout << "1:==> X-Axis " << endl;
  cout << "2:==> Y-Axis " << endl;
  cout << "3:==> Z-Axis" << endl;
  cout << "Enter Your Choice := ";
  cin >> ch;
  switch (ch)
  {
  case 1:
    cout << " Enter P := ";
    cin >> p;
    for (i = 0; i < 20; i++)
    {
      edge[i][0] = edge[i][0] / (p * edge[i][0] + 1);
      edge[i][1] = edge[i][1] / (p * edge[i][0] + 1);
      edge[i][2] = edge[i][2] / (p * edge[i][0] + 1);
    }
    draw_cube(edge);
    break;

  case 2:
    cout << " Enter Q := ";
    cin >> q;
    for (i = 0; i < 20; i++)
    {
      edge[i][1] = edge[i][1] / (edge[i][1] * q + 1);
      edge[i][0] = edge[i][0] / (edge[i][1] * q + 1);
      edge[i][2] = edge[i][2] / (edge[i][1] * q + 1);
    }
    draw_cube(edge);
    break;

  case 3:
    cout << " Enter R := ";
    cin >> r;
    for (i = 0; i < 20; i++)
    {
      edge[i][2] = edge[i][2] / (edge[i][2] * r + 1);
      edge[i][0] = edge[i][0] / (edge[i][2] * r + 1);
      edge[i][1] = edge[i][1] / (edge[i][2] * r + 1);
    }
    draw_cube(edge);
    break;
  }
  closegraph();
}

int main()
{
  int choice;
  double edge[20][3] = {
      100, 0, 0,
      100, 100, 0,
      0, 100, 0,
      0, 100, 100,
      0, 0, 100,
      0, 0, 0,
      100, 0, 0,
      100, 0, 100,
      100, 75, 100,
      75, 100, 100,
      100, 100, 75,
      100, 100, 0,
      100, 100, 75,
      100, 75, 100,
      75, 100, 100,
      0, 100, 100,
      0, 100, 0,
      0, 0, 0,
      0, 0, 100,
      100, 0, 100};
  while (1)
  {
    cout << "1:==> Draw Cube " << endl;
    cout << "2:==> Scaling " << endl;
    cout << "3:==> Rotation " << endl;
    cout << "4:==> Reflection " << endl;
    cout << "5:==> Translation " << endl;
    cout << "6:==> Perspective Projection " << endl;
    cout << "7:==> Exit " << endl;
    cout << "Enter Your Choice := ";
    cin >> choice;
    switch (choice)
    {
    case 1:
      draw_cube(edge);
      break;

    case 2:
      scale(edge);
      break;

    case 3:
      rotate(edge);
      break;

    case 4:
      reflect(edge);
      break;

    case 5:
      translate(edge);
      break;

    case 6:
      perspect(edge);
      break;

    case 7:
      exit(0);

    default:
      cout << "\nPress A Valid Key...!!! ";
      getch();
      break;
    }
    closegraph();
  }
  
  return 0;
}

////hermite curve
#include <iostream>
#include <graphics.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

struct point
{
  int x, y;
};

void hermite(point p1, point p4, double r1, double r4)
{
  float x, y, t;
  for (t = 0.0; t <= 1.0; t += 0.001)
  {
    x = (2 * t * t * t - 3 * t * t + 1) * p1.x + (-2 * t * t * t + 3 * t * t) * p4.x + (t * t * t - 2 * t * t + t) * r1 + (t * t * t - t * t) * r4;
    y = (2 * t * t * t - 3 * t * t + 1) * p1.y + (-2 * t * t * t + 3 * t * t) * p4.y + (t * t * t - 2 * t * t + 1) * r1 + (t * t * t - t * t) * r4;
    putpixel(x, y, YELLOW);
  }
}

int main()
{
  /* request auto detection */
  int gdriver = DETECT, gmode, errorcode;

  /* initialize graphics and local variables */
  initgraph(&gdriver, &gmode, NULL);

  /* read result of initialization */
  errorcode = graphresult();

  /* an error occurred */
  if (errorcode != grOk)
  {
    printf("Graphics error: %s\n", grapherrormsg(errorcode));
    printf("Press any key to halt:");
    getch();
    exit(1);
  }

  double r1, r4;
  point p1, p2;
  cout << "Enter 2 hermite points: " << endl;
  cin >> p1.x >> p1.y >> p2.x >> p2.y;
  cout << "Enter tangents at p1 and p4: " << endl;
  cin >> r1 >> r4;
  hermite(p1, p2, r1, r4);
  putpixel(p1.x, p1.y, WHITE);
  putpixel(p2.x, p2.y, WHITE);
  getch();
  closegraph();

  return 0;
}


////bezire curve
#include<graphics.h>
#include<math.h>
#include<conio.h>
#include<stdio.h>
int main()
{
int x[4],y[4],i;
double put_x,put_y,t;
int gr=DETECT,gm;
initgraph(&gr,&gm,NULL);
printf("\n****** Bezier Curve ***********");
printf("\n Please enter x and y coordinates ");
for(i=0;i<4;i++)                 
{
scanf("%d%d",&x[i],&y[i]);
putpixel(x[i],y[i],3);                // Control Points
}

for(t=0.0;t<=1.0;t=t+0.001)             // t always lies between 0 and 1
{
put_x = pow(1-t,3)*x[0] + 3*t*pow(1-t,2)*x[1] + 3*t*t*(1-t)*x[2] + pow(t,3)*x[3]; // Formula to draw curve
put_y =  pow(1-t,3)*y[0] + 3*t*pow(1-t,2)*y[1] + 3*t*t*(1-t)*y[2] + pow(t,3)*y[3];
putpixel(put_x,put_y, WHITE);            // putting pixel 
}
getch();
closegraph();

return 0;
}



