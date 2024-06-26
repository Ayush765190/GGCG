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
