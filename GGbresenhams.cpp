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
