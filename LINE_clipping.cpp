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
