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
