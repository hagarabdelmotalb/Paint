#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_4
#define M_PI_4 (M_PI / 4)
#endif
#pragma comment(linker, "/STACK:8000000")
#define WIDTH 800
#define HEIGHT 600


#include <tchar.h>
#include <windows.h>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <stack>
#include <queue>
#include <algorithm>
#include <codecvt>
#include <iostream>

#pragma comment(lib, "comdlg32.lib")

using namespace std;

bool useClipping = false;
bool showClippingWindow = false;
bool showRectangleClippingWindow = false;
bool showCircleClippingWindow = false;
vector<POINT> splinePoints;
POINT ellipseCenter, ellipseRadiusPoint;
bool ellipseSecondClick = false;
vector<POINT> ellipsePoints;


const int winLeft = 100, winTop = 100;
const int winRight = 400, winBottom = 300;
const int squareSide = min(winRight - winLeft, winBottom - winTop);
const int squareLeft = (winLeft + winRight - squareSide) / 2;
const int squareTop = (winTop + winBottom - squareSide) / 2;
const int squareRight = squareLeft + squareSide;
const int squareBottom = squareTop + squareSide;
const int circleCenterX = (winLeft + winRight) / 2;
const int circleCenterY = (winTop + winBottom) / 2;
const int circleRadius = min((winRight - winLeft), (winBottom - winTop)) / 2;

enum ShapeType
{
    LINE_DDA,
    LINE_MID,
    LINE_PARAMETRIC,
    CIRCLE_DIRECT,
    CIRCLE_POLAR,
    CIRCLE_ITERATIVE_POLAR,
    CIRCLE_MID,
    CIRCLE_MODIFIED_MID,
    POLYGON,
    CARDINAL_SPLINE,
    ELLIPSE_DIRECT,
    ELLIPSE_POLAR,
    ELLIPSE_MIDPOINT,
    POINT_CLIP_RECT,
    POINT_CLIP_SQUARE,
    LINE_CLIP_RECT,
    LINE_CLIP_SQUARE,
    POLYGON_CLIP_RECT,
    POINT_CLIP_CIRCLE,
    LINE_CLIP_CIRCLE,
    NONE
};

enum FillType
{
    FILL_RECURSIVE,
    FILL_NON_RECURSIVE,
    FILL_CONVEX,
    FILL_NON_CONVEX,
    FILL_NONE
};

struct Shape
{
    ShapeType type;
    POINT p1, p2;
    COLORREF color;
    vector<POINT> points;
    bool isPolygon;
    FillType fillType;
    COLORREF fillColor;

    Shape() : type(NONE), p1({ 0, 0 }), p2({ 0, 0 }), color(RGB(0, 0, 0)), isPolygon(false), fillType(FILL_NONE), fillColor(RGB(0, 0, 0)) {}
    Shape(ShapeType t, POINT pt1, POINT pt2, COLORREF c)
            : type(t), p1(pt1), p2(pt2), color(c), isPolygon(false), fillType(FILL_NONE), fillColor(RGB(0, 0, 0)) {}
    Shape(ShapeType t, const vector<POINT>& pts, COLORREF c, bool poly)
            : type(t), points(pts), color(c), isPolygon(poly), fillType(FILL_NONE), fillColor(RGB(0, 0, 0)) {}
};

vector<Shape> shapes;
vector<POINT> polygonPoints;
COLORREF currentColor = RGB(0, 0, 0);
ShapeType currentShape = NONE;
FillType currentFill = FILL_NONE;
POINT tempPoint;
bool isFirstClick = true;
int lastCircleX = -1, lastCircleY = -1, lastCircleR = -1;
bool hasLastCircle = false;
int selectedQuarter = 1;

TCHAR szClassName[] = _T("DrawingApp");

const int INSIDE = 0, LEFT = 1, RIGHT = 2, BOTTOM = 4, TOP = 8;

int ComputeCodeRectangle(int x, int y)
{
    int code = INSIDE;
    if (x < winLeft) code |= LEFT;
    else if (x > winRight) code |= RIGHT;
    if (y < winTop) code |= TOP;
    else if (y > winBottom) code |= BOTTOM;
    return code;
}

int ComputeCodeSquare(int x, int y)
{
    int code = INSIDE;
    if (x < squareLeft) code |= LEFT;
    else if (x > squareRight) code |= RIGHT;
    if (y < squareTop) code |= TOP;
    else if (y > squareBottom) code |= BOTTOM;
    return code;
}

bool IsPointInsideCircle(int x, int y)
{
    int dx = x - circleCenterX;
    int dy = y - circleCenterY;
    return (dx * dx + dy * dy) <= (circleRadius * circleRadius);
}

void DrawPixel(HDC hdc, int x, int y, COLORREF color)
{
    SetPixel(hdc, x, y, color);
}

void ClipPointRectangle(HDC hdc, int x, int y, COLORREF color)
{
    if (x >= winLeft && x <= winRight && y >= winTop && y <= winBottom)
    {
        HBRUSH brush = CreateSolidBrush(color);
        SelectObject(hdc, brush);
        Ellipse(hdc, x - 4, y - 4, x + 4, y + 4);
        DeleteObject(brush);
    }
}

void ClipPointSquare(HDC hdc, int x, int y, COLORREF color)
{
    if (x >= squareLeft && x <= squareRight && y >= squareTop && y <= squareBottom)
    {
        HBRUSH brush = CreateSolidBrush(color);
        SelectObject(hdc, brush);
        Ellipse(hdc, x - 4, y - 4, x + 4, y + 4);
        DeleteObject(brush);
    }
}

void ClipPointCircle(HDC hdc, int x, int y, COLORREF color)
{
    if (IsPointInsideCircle(x, y))
    {
        HBRUSH brush = CreateSolidBrush(color);
        SelectObject(hdc, brush);
        Ellipse(hdc, x - 4, y - 4, x + 4, y + 4);
        DeleteObject(brush);
    }
}

void ClipLineRectangle(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    int code1 = ComputeCodeRectangle(x1, y1);
    int code2 = ComputeCodeRectangle(x2, y2);
    bool accept = false;

    while (true)
    {
        if ((code1 | code2) == 0)
        {
            accept = true;
            break;
        }
        else if (code1 & code2)
        {
            break;
        }
        else
        {
            int codeOut = code1 ? code1 : code2;
            int x, y;

            if (codeOut & TOP)
            {
                x = x1 + (x2 - x1) * (winTop - y1) / (y2 - y1);
                y = winTop;
            }
            else if (codeOut & BOTTOM)
            {
                x = x1 + (x2 - x1) * (winBottom - y1) / (y2 - y1);
                y = winBottom;
            }
            else if (codeOut & RIGHT)
            {
                y = y1 + (y2 - y1) * (winRight - x1) / (x2 - x1);
                x = winRight;
            }
            else
            {
                y = y1 + (y2 - y1) * (winLeft - x1) / (x2 - x1);
                x = winLeft;
            }

            if (codeOut == code1)
            {
                x1 = x;
                y1 = y;
                code1 = ComputeCodeRectangle(x1, y1);
            }
            else
            {
                x2 = x;
                y2 = y;
                code2 = ComputeCodeRectangle(x2, y2);
            }
        }
    }

    if (accept)
    {
        HPEN pen = CreatePen(PS_SOLID, 2, color);
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        MoveToEx(hdc, x1, y1, NULL);
        LineTo(hdc, x2, y2);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);
    }
}

void ClipLineSquare(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    int code1 = ComputeCodeSquare(x1, y1);
    int code2 = ComputeCodeSquare(x2, y2);
    bool accept = false;

    while (true)
    {
        if ((code1 | code2) == 0)
        {
            accept = true;
            break;
        }
        else if (code1 & code2)
        {
            break;
        }
        else
        {
            int codeOut = code1 ? code1 : code2;
            int x, y;

            if (codeOut & TOP)
            {
                x = x1 + (x2 - x1) * (squareTop - y1) / (y2 - y1);
                y = squareTop;
            }
            else if (codeOut & BOTTOM)
            {
                x = x1 + (x2 - x1) * (squareBottom - y1) / (y2 - y1);
                y = squareBottom;
            }
            else if (codeOut & RIGHT)
            {
                y = y1 + (y2 - y1) * (squareRight - x1) / (x2 - x1);
                x = squareRight;
            }
            else
            {
                y = y1 + (y2 - y1) * (squareLeft - x1) / (x2 - x1);
                x = squareLeft;
            }

            if (codeOut == code1)
            {
                x1 = x;
                y1 = y;
                code1 = ComputeCodeSquare(x1, y1);
            }
            else
            {
                x2 = x;
                y2 = y;
                code2 = ComputeCodeSquare(x2, y2);
            }
        }
    }

    if (accept)
    {
        HPEN pen = CreatePen(PS_SOLID, 2, color);
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        MoveToEx(hdc, x1, y1, NULL);
        LineTo(hdc, x2, y2);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);
    }
}

void ClipLineCircle(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    bool inside1 = IsPointInsideCircle(x1, y1);
    bool inside2 = IsPointInsideCircle(x2, y2);

    if (inside1 && inside2)
    {
        HPEN pen = CreatePen(PS_SOLID, 2, color);
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        MoveToEx(hdc, x1, y1, NULL);
        LineTo(hdc, x2, y2);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);
        return;
    }

    double dx = x2 - x1;
    double dy = y2 - y1;
    double a = dx * dx + dy * dy;
    double b = 2 * (dx * (x1 - circleCenterX) + dy * (y1 - circleCenterY));
    double c = (x1 - circleCenterX) * (x1 - circleCenterX) + (y1 - circleCenterY) * (y1 - circleCenterY) - circleRadius * circleRadius;

    double disc = b * b - 4 * a * c;
    if (disc < 0)
        return;

    double t1 = (-b + sqrt(disc)) / (2 * a);
    double t2 = (-b - sqrt(disc)) / (2 * a);

    vector<POINT> intersections;
    if (t1 >= 0 && t1 <= 1)
    {
        int x = (int)(x1 + t1 * dx);
        int y = (int)(y1 + t1 * dy);
        intersections.push_back({ x, y });
    }
    if (t2 >= 0 && t2 <= 1 && t2 != t1)
    {
        int x = (int)(x1 + t2 * dx);
        int y = (int)(y1 + t2 * dy);
        intersections.push_back({ x, y });
    }

    if (inside1 && !inside2 && !intersections.empty())
    {
        HPEN pen = CreatePen(PS_SOLID, 2, color);
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        MoveToEx(hdc, x1, y1, NULL);
        LineTo(hdc, intersections[0].x, intersections[0].y);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);
    }
    else if (!inside1 && inside2 && !intersections.empty())
    {
        HPEN pen = CreatePen(PS_SOLID, 2, color);
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        MoveToEx(hdc, intersections[0].x, intersections[0].y, NULL);
        LineTo(hdc, x2, y2);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);
    }
    else if (!intersections.empty() && intersections.size() == 2)
    {
        HPEN pen = CreatePen(PS_SOLID, 2, color);
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        MoveToEx(hdc, intersections[0].x, intersections[0].y, NULL);
        LineTo(hdc, intersections[1].x, intersections[1].y);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);
    }
}

void ClipPolygonRectangle(HDC hdc, const vector<POINT>& points, COLORREF color)
{
    if (points.size() < 2) return;

    struct Edge { int x1, y1, x2, y2; };
    vector<Edge> edges;
    for (size_t i = 0; i < points.size(); i++)
    {
        size_t j = (i + 1) % points.size();
        edges.push_back({ points[i].x, points[i].y, points[j].x, points[j].y });
    }

    vector<Edge> newEdges;
    auto clipEdge = [&](int boundary, bool isX, int value, bool isMin) {
        newEdges.clear();
        for (const auto& e : edges)
        {
            int code1 = ComputeCodeRectangle(e.x1, e.y1);
            int code2 = ComputeCodeRectangle(e.x2, e.y2);
            int flag = isX ? (isMin ? LEFT : RIGHT) : (isMin ? TOP : BOTTOM);
            int limit = value;

            bool inside1 = !(code1 & flag);
            bool inside2 = !(code2 & flag);

            if (inside1 && inside2)
            {
                newEdges.push_back(e);
            }
            else if (inside1)
            {
                int x, y;
                if (isX)
                {
                    y = e.y1 + (e.y2 - e.y1) * (limit - e.x1) / (e.x2 - e.x1);
                    x = limit;
                }
                else
                {
                    x = e.x1 + (e.x2 - e.x1) * (limit - e.y1) / (e.y2 - e.y1);
                    y = limit;
                }
                newEdges.push_back({ e.x1, e.y1, x, y });
            }
            else if (inside2)
            {
                int x, y;
                if (isX)
                {
                    y = e.y1 + (e.y2 - e.y1) * (limit - e.x1) / (e.x2 - e.x1);
                    x = limit;
                }
                else
                {
                    x = e.x1 + (e.x2 - e.x1) * (limit - e.y1) / (e.y2 - e.y1);
                    y = limit;
                }
                newEdges.push_back({ x, y, e.x2, e.y2 });
            }
            else
            {
                int x, y;
                if (isX)
                {
                    if ((e.x1 < limit && e.x2 > limit) || (e.x1 > limit && e.x2 < limit))
                    {
                        y = e.y1 + (e.y2 - e.y1) * (limit - e.x1) / (e.x2 - e.x1);
                        if (y >= winTop && y <= winBottom)
                        {
                            newEdges.push_back({ limit, y, limit, y });
                        }
                    }
                }
                else
                {
                    if ((e.y1 < limit && e.y2 > limit) || (e.y1 > limit && e.y2 < limit))
                    {
                        x = e.x1 + (e.x2 - e.x1) * (limit - e.y1) / (e.y2 - e.y1);
                        if (x >= winLeft && x <= winRight)
                        {
                            newEdges.push_back({ x, limit, x, limit });
                        }
                    }
                }
            }
        }
        edges = newEdges;
    };

    clipEdge(LEFT, true, winLeft, true);
    clipEdge(RIGHT, true, winRight, false);
    clipEdge(TOP, false, winTop, true);
    clipEdge(BOTTOM, false, winBottom, false);

    for (const auto& e : edges)
    {
        ClipLineRectangle(hdc, e.x1, e.y1, e.x2, e.y2, color);
    }
}

void DrawLineDDA(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    float xinc = dx / (float)steps;
    float yinc = dy / (float)steps;
    float x = x1, y = y1;
    for (int i = 0; i <= steps; ++i)
    {
        DrawPixel(hdc, round(x), round(y), color);
        x += xinc;
        y += yinc;
    }
}
void Draw8Points(HDC hdc, int xc, int yc, int x, int y) {
    //std::cout << "Drawing at: (" << xc + x << ", " << yc + y << ")\n";


    SetPixel(hdc, xc + x, yc + y, RGB(0, 128, 0));
    SetPixel(hdc, xc - x, yc + y, RGB(0, 128, 0));
    SetPixel(hdc, xc + x, yc - y, RGB(0, 128, 0));
    SetPixel(hdc, xc - x, yc - y, RGB(0, 128, 0));
}

void EllipseMidpoint(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
    int x = 0;
    int y = b;
    int aSq = a * a;
    int bSq = b * b;
    int d = bSq - aSq * b + (aSq / 4);
    COLORREF originalColor = SetTextColor(hdc, color); // Set the drawing color

    // Region 1
    while (aSq * (y - 0.5) > bSq * (x + 1))
    {
        SetPixel(hdc, xc + x, yc + y, color);
        SetPixel(hdc, xc - x, yc + y, color);
        SetPixel(hdc, xc + x, yc - y, color);
        SetPixel(hdc, xc - x, yc - y, color);

        if (d < 0)
        {
            d += bSq * (2 * x + 3);
        }
        else
        {
            d += bSq * (2 * x + 3) + aSq * (-2 * y + 2);
            y--;
        }
        x++;
    }

    // Region 2
    d = bSq * (x + 0.5) * (x + 0.5) + aSq * (y - 1) * (y - 1) - aSq * bSq;
    while (y > 0)
    {
        SetPixel(hdc, xc + x, yc + y, color);
        SetPixel(hdc, xc - x, yc + y, color);
        SetPixel(hdc, xc + x, yc - y, color);
        SetPixel(hdc, xc - x, yc - y, color);

        if (d < 0)
        {
            d += bSq * (2 * x + 2) + aSq * (-2 * y + 3);
            x++;
        }
        else
        {
            d += aSq * (-2 * y + 3);
        }
        y--;
    }

    SetTextColor(hdc, originalColor); // Restore the original color
}

void DrawEllipse_Polar(HDC hdc, int xc, int yc, int a, int b, COLORREF color) {
    double theta = 0;
    double dtheta = 1.0 / max(a, b);

    while (theta < 2 * 3.14159) {
        int x = (int)(a * cos(theta));
        int y = (int)(b * sin(theta));
        SetPixel(hdc, xc + x, yc + y, color);
        theta += dtheta;
    }
}

void DrawLineMidpoint(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    int dx = abs(x2 - x1), dy = abs(y2 - y1);
    int sx = x1 < x2 ? 1 : -1;
    int sy = y1 < y2 ? 1 : -1;
    bool steep = dy > dx;
    if (steep)
    {
        swap(x1, y1);
        swap(x2, y2);
        swap(dx, dy);
        swap(sx, sy);
    }
    int d = 2 * dy - dx;
    int y = y1;
    for (int x = x1; sx == 1 ? x <= x2 : x >= x2; x += sx)
    {
        if (steep)
            DrawPixel(hdc, y, x, color);
        else
            DrawPixel(hdc, x, y, color);
        if (d > 0)
        {
            y += sy;
            d -= 2 * dx;
        }
        d += 2 * dy;
    }
}

void DrawLineParametric(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    float length = sqrt(dx * dx + dy * dy);
    int steps = (int)length;
    if (steps == 0) steps = 1;
    float t_step = 1.0f / steps;
    for (float t = 0; t <= 1.0f; t += t_step)
    {
        int x = (int)round(x1 + t * dx);
        int y = (int)round(y1 + t * dy);
        DrawPixel(hdc, x, y, color);
    }
}

void DrawCircleDirect(HDC hdc, int xc, int yc, int r, COLORREF color)
{
    for (int x = 0; x <= r / sqrt(2); ++x)
    {
        int y = (int)round(sqrt(r * r - x * x));
        DrawPixel(hdc, xc + x, yc + y, color);
        DrawPixel(hdc, xc - x, yc + y, color);
        DrawPixel(hdc, xc + x, yc - y, color);
        DrawPixel(hdc, xc - x, yc - y, color);
        DrawPixel(hdc, xc + y, yc + x, color);
        DrawPixel(hdc, xc - y, yc + x, color);
        DrawPixel(hdc, xc + y, yc - x, color);
        DrawPixel(hdc, xc - y, yc - x, color);
    }
}

void DrawCirclePolar(HDC hdc, int xc, int yc, int r, COLORREF color)
{
    double theta = 0;
    double step = 1.0 / r;
    while (theta <= M_PI_4)
    {
        int x = (int)round(r * cos(theta));
        int y = (int)round(r * sin(theta));
        DrawPixel(hdc, xc + x, yc + y, color);
        DrawPixel(hdc, xc - x, yc + y, color);
        DrawPixel(hdc, xc + x, yc - y, color);
        DrawPixel(hdc, xc - x, yc - y, color);
        DrawPixel(hdc, xc + y, yc + x, color);
        DrawPixel(hdc, xc - y, yc + x, color);
        DrawPixel(hdc, xc + y, yc - x, color);
        DrawPixel(hdc, xc - y, yc - x, color);
        theta += step;
    }
}

void DrawCircleIterativePolar(HDC hdc, int xc, int yc, int r, COLORREF color)
{
    double x = r;
    double y = 0;
    double theta = 1.0 / r;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    while (x >= y)
    {
        int ix = (int)round(x);
        int iy = (int)round(y);
        DrawPixel(hdc, xc + ix, yc + iy, color);
        DrawPixel(hdc, xc - ix, yc + iy, color);
        DrawPixel(hdc, xc + ix, yc - iy, color);
        DrawPixel(hdc, xc - ix, yc - iy, color);
        DrawPixel(hdc, xc + iy, yc + ix, color);
        DrawPixel(hdc, xc - iy, yc + ix, color);
        DrawPixel(hdc, xc + iy, yc - ix, color);
        DrawPixel(hdc, xc - iy, yc - ix, color);

        double x_new = x * cos_theta - y * sin_theta;
        double y_new = x * sin_theta + y * cos_theta;
        x = x_new;
        y = y_new;
    }
}

void DrawCircleMidpoint(HDC hdc, int xc, int yc, int r, COLORREF color)
{
    int x = 0, y = r;
    int d = 1 - r;
    while (x <= y)
    {
        DrawPixel(hdc, xc + x, yc + y, color);
        DrawPixel(hdc, xc - x, yc + y, color);
        DrawPixel(hdc, xc + x, yc - y, color);
        DrawPixel(hdc, xc - x, yc - y, color);
        DrawPixel(hdc, xc + y, yc + x, color);
        DrawPixel(hdc, xc - y, yc + x, color);
        DrawPixel(hdc, xc + y, yc - x, color);
        DrawPixel(hdc, xc - y, yc - x, color);

        if (d < 0)
            d += 2 * x + 3;
        else
        {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
    }
}

void DrawCircleModifiedMidpoint(HDC hdc, int xc, int yc, int r, COLORREF color)
{
    int x = 0, y = r;
    int d = 1 - r;
    int deltaE = 3;
    int deltaSE = 5 - 2 * r;

    while (x <= y)
    {
        DrawPixel(hdc, xc + x, yc + y, color);
        DrawPixel(hdc, xc - x, yc + y, color);
        DrawPixel(hdc, xc + x, yc - y, color);
        DrawPixel(hdc, xc - x, yc - y, color);
        DrawPixel(hdc, xc + y, yc + x, color);
        DrawPixel(hdc, xc - y, yc + x, color);
        DrawPixel(hdc, xc + y, yc - x, color);
        DrawPixel(hdc, xc - y, yc - x, color);

        if (d < 0)
        {
            d += deltaE;
            deltaE += 2;
            deltaSE += 2;
        }
        else
        {
            d += deltaSE;
            deltaE += 2;
            deltaSE += 4;
            y--;
        }
        x++;
    }
}


void DrawEllipse_Direct(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
    for (int x = 0; x <= a; ++x)
    {
        double y = b * sqrt(1.0 - (x * x) / (double)(a * a));
        DrawPixel(hdc, xc + x, yc + (int)y, color);
        DrawPixel(hdc, xc - x, yc + (int)y, color);
        DrawPixel(hdc, xc + x, yc - (int)y, color);
        DrawPixel(hdc, xc - x, yc - (int)y, color);
    }

    for (int y = 0; y <= b; ++y)
    {
        double x = a * sqrt(1.0 - (y * y) / (double)(b * b));
        DrawPixel(hdc, xc + (int)x, yc + y, color);
        DrawPixel(hdc, xc - (int)x, yc + y, color);
        DrawPixel(hdc, xc + (int)x, yc - y, color);
        DrawPixel(hdc, xc - (int)x, yc - y, color);
    }
}

void DrawPolygon(HDC hdc, const vector<POINT>& points, COLORREF color)
{
    if (points.size() < 2) return;

    for (size_t i = 0; i < points.size() - 1; i++)
    {
        DrawLineDDA(hdc, points[i].x, points[i].y, points[i + 1].x, points[i + 1].y, color);
    }
    DrawLineDDA(hdc, points.back().x, points.back().y, points[0].x, points[0].y, color);
}

void FillConvexPolygon(HDC hdc, const vector<POINT>& points, COLORREF color)
{
    if (points.size() < 3) return;

    long ymin = points[0].y;
    long ymax = points[0].y;
    for (const auto& pt : points)
    {
        ymin = min(ymin, pt.y);
        ymax = max(ymax, pt.y);
    }

    for (int y = ymin; y <= ymax; y++)
    {
        vector<int> intersections;
        for (size_t i = 0; i < points.size(); i++)
        {
            size_t j = (i + 1) % points.size();
            int x1 = points[i].x, y1 = points[i].y;
            int x2 = points[j].x, y2 = points[j].y;

            if ((y1 <= y && y2 > y) || (y2 <= y && y1 > y))
            {
                float x = x1 + (float)(y - y1) * (x2 - x1) / (y2 - y1);
                intersections.push_back((int)round(x));
            }
        }

        sort(intersections.begin(), intersections.end());
        for (size_t i = 0; i < intersections.size(); i += 2)
        {
            if (i + 1 >= intersections.size()) break;
            int x1 = intersections[i];
            int x2 = intersections[i + 1];
            for (int x = x1; x <= x2; x++)
            {
                DrawPixel(hdc, x, y, color);
            }
        }
    }
}

void FillNonConvexPolygon(HDC hdc, const vector<POINT>& points, COLORREF color)
{
    if (points.size() < 3) return;

    long ymin = points[0].y;
    long ymax = points[0].y;
    for (const auto& pt : points)
    {
        ymin = min(ymin, pt.y);
        ymax = max(ymax, pt.y);
    }

    for (int y = ymin; y <= ymax; y++)
    {
        vector<float> intersections;
        for (size_t i = 0; i < points.size(); i++)
        {
            size_t j = (i + 1) % points.size();
            int x1 = points[i].x, y1 = points[i].y;
            int x2 = points[j].x, y2 = points[j].y;

            if ((y1 <= y && y2 > y) || (y2 <= y && y1 > y))
            {
                if (y2 != y1)
                {
                    float x = x1 + (float)(y - y1) * (x2 - x1) / (y2 - y1);
                    intersections.push_back(x);
                }
            }
        }

        sort(intersections.begin(), intersections.end());
        for (size_t i = 0; i < intersections.size(); i += 2)
        {
            if (i + 1 >= intersections.size()) break;
            int x1 = (int)round(intersections[i]);
            int x2 = (int)round(intersections[i + 1]);
            for (int x = x1; x <= x2; x++)
            {
                DrawPixel(hdc, x, y, color);
            }
        }
    }
}

void FloodFillRecursive(HDC hdc, int x, int y, COLORREF fillColor, COLORREF borderColor, RECT* clientRect, std::vector<std::vector<bool>>* visited)
{
    if (x < clientRect->left || x >= clientRect->right || y < clientRect->top || y >= clientRect->bottom)
        return;

    int visitedX = x - clientRect->left;
    int visitedY = y - clientRect->top;

    if ((*visited)[visitedY][visitedX])
        return;

    COLORREF current = GetPixel(hdc, x, y);
    static COLORREF targetColor = current;

    if (x == clientRect->left && y == clientRect->top)
        targetColor = current;

    if (current != targetColor || current == fillColor || current == borderColor)
        return;

    (*visited)[visitedY][visitedX] = true;
    SetPixelV(hdc, x, y, fillColor);

    FloodFillRecursive(hdc, x + 1, y, fillColor, borderColor, clientRect, visited);
    FloodFillRecursive(hdc, x - 1, y, fillColor, borderColor, clientRect, visited);
    FloodFillRecursive(hdc, x, y + 1, fillColor, borderColor, clientRect, visited);
    FloodFillRecursive(hdc, x, y - 1, fillColor, borderColor, clientRect, visited);
}

void FloodFillRecursiveWrapper(HDC hdc, int x, int y, COLORREF fillColor, COLORREF borderColor)
{
    RECT clipRect;
    GetClipBox(hdc, &clipRect);

    std::vector<std::vector<bool>> visited(
            clipRect.bottom - clipRect.top,
            std::vector<bool>(clipRect.right - clipRect.left, false)
    );

    FloodFillRecursive(hdc, x, y, fillColor, borderColor, &clipRect, &visited);
}

void FloodFillNonRecursive(HDC hdc, int x, int y, COLORREF fillColor, COLORREF borderColor)
{
    RECT clipRect;
    GetClipBox(hdc, &clipRect);

    COLORREF targetColor = GetPixel(hdc, x, y);

    if (targetColor == fillColor || targetColor == borderColor)
        return;

    std::vector<std::vector<bool>> visited(
            clipRect.bottom - clipRect.top,
            std::vector<bool>(clipRect.right - clipRect.left, false)
    );

    std::stack<POINT> pixels;
    pixels.push({ x, y });

    while (!pixels.empty())
    {
        POINT pt = pixels.top();
        pixels.pop();

        if (pt.x < clipRect.left || pt.x >= clipRect.right ||
            pt.y < clipRect.top || pt.y >= clipRect.bottom)
            continue;

        int visitedX = pt.x - clipRect.left;
        int visitedY = pt.y - clipRect.top;

        if (visited[visitedY][visitedX])
            continue;

        visited[visitedY][visitedX] = true;

        if (GetPixel(hdc, pt.x, pt.y) == targetColor)
        {
            SetPixelV(hdc, pt.x, pt.y, fillColor);

            pixels.push({ pt.x - 1, pt.y });
            pixels.push({ pt.x + 1, pt.y });
            pixels.push({ pt.x, pt.y - 1 });
            pixels.push({ pt.x, pt.y + 1 });
        }
    }
}

void FillCircle_Lines(HDC hdc, int xc, int yc, int R, int Q, COLORREF C)
{
    for (int angle = 0; angle <= 90; angle++)
    {
        float angle_Rad = angle * M_PI / 180.0;
        int x = (int)(R * cos(angle_Rad));
        int y = (int)(R * sin(angle_Rad));
        int x2 = xc, y2 = yc;
        switch (Q)
        {
            case 1:
                x2 += x;
                y2 += y;
                break;
            case 2:
                x2 -= x;
                y2 += y;
                break;
            case 3:
                x2 -= x;
                y2 -= y;
                break;
            case 4:
                x2 += x;
                y2 -= y;
                break;
            default:
                return;
        }
        DrawLineDDA(hdc, xc, yc, x2, y2, C);
    }
}

void FillCircleInQuarter(HDC hdc, int xc, int yc, int R, int Q, COLORREF C)
{
    for (int r = 6; r <= R; r += 10)
    {
        for (int angle = 0; angle <= 90; angle += 3)
        {
            float angle_Rad = angle * M_PI / 180.0;
            int x = (int)(r * cos(angle_Rad));
            int y = (int)(r * sin(angle_Rad));
            int x2 = xc, y2 = yc;
            switch (Q)
            {
                case 1:
                    x2 += x;
                    y2 += y;
                    break;
                case 2:
                    x2 -= x;
                    y2 += y;
                    break;
                case 3:
                    x2 -= x;
                    y2 -= y;
                    break;
                case 4:
                    x2 += x;
                    y2 -= y;
                    break;
                default:
                    return;
            }
            DrawCircleMidpoint(hdc, x2, y2, 1, C);
        }
    }
}

void FillCircle_Circle(HDC hdc, int xc, int yc, int R, COLORREF C)
{
    for (int r = R; r >= 1; r -= 5)
    {
        DrawCircleMidpoint(hdc, xc, yc, r, C);
    }
}

void DrawHermiteCurve(HDC hdc, POINT p0, POINT t0, POINT p1, POINT t1, COLORREF color)
{
    for (double t = 0; t <= 1; t += 0.001)
    {
        double h1 = 2 * t * t * t - 3 * t * t + 1;
        double h2 = -2 * t * t * t + 3 * t * t;
        double h3 = t * t * t - 2 * t * t + t;
        double h4 = t * t * t - t * t;

        double x = h1 * p0.x + h2 * p1.x + h3 * t0.x + h4 * t1.x;
        double y = h1 * p0.y + h2 * p1.y + h3 * t0.y + h4 * t1.y;
        DrawPixel(hdc, round(x), round(y), color);
    }
}

void FillSquareWithHermite(HDC hdc, int left, int top, int size, COLORREF color)
{
    for (int x = 0; x <= size; x += 2)
    {
        POINT p0 = { left + x, top };
        POINT p1 = { left + x, top + size };
        POINT t0 = { 20, 0 };
        POINT t1 = { -20, 0 };
        DrawHermiteCurve(hdc, p0, t0, p1, t1, color);
    }
}

void DrawBezierCurve(HDC hdc, POINT p0, POINT p1, POINT p2, POINT p3, COLORREF color)
{
    for (double t = 0; t <= 1; t += 0.001)
    {
        double x = pow(1 - t, 3) * p0.x +
                   3 * t * pow(1 - t, 2) * p1.x +
                   3 * t * t * (1 - t) * p2.x +
                   t * t * t * p3.x;

        double y = pow(1 - t, 3) * p0.y +
                   3 * t * pow(1 - t, 2) * p1.y +
                   3 * t * t * (1 - t) * p2.y +
                   t * t * t * p3.y;

        DrawPixel(hdc, round(x), round(y), color);
    }
}

void FillRectangleWithBezier(HDC hdc, int left, int top, int width, int height, COLORREF color)
{
    for (int y = 0; y <= height; y += 2)
    {
        POINT p0 = { left, top + y };
        POINT p3 = { left + width, top + y };
        POINT p1 = { left + width / 3, top + y + 20 };
        POINT p2 = { left + 2 * width / 3, top + y - 20 };
        DrawBezierCurve(hdc, p0, p1, p2, p3, color);
    }
}

void DrawCardinalSpline(HDC hdc, const vector<POINT>& pts, COLORREF color, float tension = 0.5f)
{
    if (pts.size() < 4) return;

    for (size_t i = 1; i + 2 < pts.size(); ++i)
    {
        POINT p0 = pts[i - 1];
        POINT p1 = pts[i];
        POINT p2 = pts[i + 1];
        POINT p3 = pts[i + 2];

        for (float t = 0; t <= 1.0f; t += 0.01f)
        {
            float t2 = t * t;
            float t3 = t2 * t;

            float b1 = -tension * t3 + 2 * tension * t2 - tension * t;
            float b2 = (2 - tension) * t3 + (tension - 3) * t2 + 1;
            float b3 = (tension - 2) * t3 + (3 - 2 * tension) * t2 + tension * t;
            float b4 = tension * t3 - tension * t2;

            int x = (int)(b1 * p0.x + b2 * p1.x + b3 * p2.x + b4 * p3.x);
            int y = (int)(b1 * p0.y + b2 * p1.y + b3 * p2.y + b4 * p3.y);

            DrawPixel(hdc, x, y, color);
        }
    }
}

void ClearScreen(HWND hwnd)
{
    shapes.clear();
    polygonPoints.clear();
    splinePoints.clear();
    isFirstClick = true;
    ellipseSecondClick = false;
    InvalidateRect(hwnd, NULL, TRUE);
}

void SaveShapes(HWND hwnd, const wstring& file_path)
{
    std::wofstream out(file_path.c_str());
    if (!out)
    {
        TCHAR errorMsg[256];
        _stprintf_s(errorMsg, _T("Failed to open file for saving: %s"), file_path.c_str());
        MessageBox(hwnd, errorMsg, _T("Error"), MB_OK | MB_ICONERROR);
        return;
    }


    for (const auto& s : shapes)
    {
        out << static_cast<int>(s.type) << L' ';
        switch (s.type)
        {
            case POINT_CLIP_RECT:
            case POINT_CLIP_SQUARE:
            case POINT_CLIP_CIRCLE:
                out << 1 << L' ' << s.p1.x << L' ' << s.p1.y << L' ';
                break;
            case LINE_DDA:
            case LINE_MID:
            case LINE_PARAMETRIC:
            case LINE_CLIP_RECT:
            case LINE_CLIP_SQUARE:
            case LINE_CLIP_CIRCLE:
            case CIRCLE_DIRECT:
            case CIRCLE_POLAR:
            case CIRCLE_ITERATIVE_POLAR:
            case CIRCLE_MID:
            case CIRCLE_MODIFIED_MID:
            case ELLIPSE_DIRECT:
                out << 2 << L' ' << s.p1.x << L' ' << s.p1.y << L' ' << s.p2.x << L' ' << s.p2.y << L' ';
                break;
            case POLYGON:
            case POLYGON_CLIP_RECT:
            case CARDINAL_SPLINE:
                out << s.points.size() << L' ';
                for (const auto& pt : s.points)
                {
                    out << pt.x << L' ' << pt.y << L' ';
                }
                break;
            default:
                continue; // Skip invalid shape types
        }
        out << static_cast<int>(s.color) << L' ' << s.isPolygon << L' '
            << static_cast<int>(s.fillType) << L' ' << static_cast<int>(s.fillColor) << L'\n';
    }

    out.close();
}

void LoadShapes(HWND hwnd, const wstring& file_path)
{
    shapes.clear();
    polygonPoints.clear();
    splinePoints.clear();
    hasLastCircle = false;
    lastCircleX = -1;
    lastCircleY = -1;
    lastCircleR = -1;

    std::wifstream in(file_path.c_str());
    if (!in.is_open())
    {
        TCHAR errorMsg[256];
        _stprintf_s(errorMsg, _T("Failed to open file for loading: %s"), file_path.c_str());
        MessageBox(hwnd, errorMsg, _T("Error"), MB_OK | MB_ICONERROR);
        return;
    }


    in.imbue(std::locale(in.getloc(), new std::codecvt_utf8<wchar_t>));

    Shape s;
    while (in.good())
    {
        int type, pointCount, color, isPolygon, fillType, fillColor;
        in >> type >> pointCount;
        if (in.eof()) break;
        if (in.fail() || type < 0 || type >= static_cast<int>(NONE))
        {
            MessageBox(hwnd, _T("Invalid shape type in file!"), _T("Error"), MB_OK | MB_ICONERROR);
            in.close();
            return;
        }
        s.type = static_cast<ShapeType>(type);

        bool validPointCount = true;
        wstring errorMsg = L"Invalid point count for shape type " + to_wstring(type) + L" (expected: ";
        switch (s.type)
        {
            case POINT_CLIP_RECT:
            case POINT_CLIP_SQUARE:
            case POINT_CLIP_CIRCLE:
                if (pointCount != 1) validPointCount = false;
                errorMsg += L"1, got: " + to_wstring(pointCount) + L")";
                break;
            case LINE_DDA:
            case LINE_MID:
            case LINE_PARAMETRIC:
            case LINE_CLIP_RECT:
            case LINE_CLIP_SQUARE:
            case LINE_CLIP_CIRCLE:
            case CIRCLE_DIRECT:
            case CIRCLE_POLAR:
            case CIRCLE_ITERATIVE_POLAR:
            case CIRCLE_MID:
            case CIRCLE_MODIFIED_MID:
            case ELLIPSE_DIRECT:
                if (pointCount != 2) validPointCount = false;
                errorMsg += L"2, got: " + to_wstring(pointCount) + L")";
                break;
            case POLYGON:
            case POLYGON_CLIP_RECT:
            case CARDINAL_SPLINE:
                if (pointCount < 2) validPointCount = false;
                errorMsg += L"at least 2, got: " + to_wstring(pointCount) + L")";
                break;
            default:
                validPointCount = false;
                errorMsg += L"unknown type)";
        }

        if (!validPointCount)
        {
            TCHAR msg[256];
            _stprintf_s(msg, _T("%s"), errorMsg.c_str());
            MessageBox(hwnd, msg, _T("Error"), MB_OK | MB_ICONERROR);
            in.close();
            return;
        }

        s.points.clear();
        for (int i = 0; i < pointCount; ++i)
        {
            int x, y;
            in >> x >> y;
            s.points.push_back({ x, y });
        }
        in >> color >> isPolygon >> fillType >> fillColor;
        s.color = static_cast<COLORREF>(color);
        s.isPolygon = (isPolygon != 0);
        s.fillType = static_cast<FillType>(fillType);
        if (s.fillType > FILL_NONE) s.fillType = FILL_NONE; // Default to FILL_NONE if invalid
        s.fillColor = static_cast<COLORREF>(fillColor);

        // Assign p1 and p2 for two-point shapes
        if (pointCount >= 1) s.p1 = s.points[0];
        if (pointCount >= 2) s.p2 = s.points[1];

        shapes.push_back(s);
    }

    in.close();
    InvalidateRect(hwnd, NULL, TRUE);
}

COLORREF GetCurrentColor()
{
    CHOOSECOLOR cc = { sizeof(CHOOSECOLOR) };
    static COLORREF acrCustClr[16];
    cc.lpCustColors = acrCustClr;
    cc.Flags = CC_RGBINIT | CC_FULLOPEN;
    cc.rgbResult = currentColor;
    if (ChooseColor(&cc)) currentColor = cc.rgbResult;
    return currentColor;
}

void AddMenus(HWND hwnd)
{
    HMENU hMenubar = CreateMenu();
    HMENU hFile = CreateMenu();
    HMENU hDraw = CreateMenu();
    HMENU hFill = CreateMenu();
    HMENU hOption = CreateMenu();
    HMENU hClipping = CreateMenu();
    HMENU hFillLines = CreateMenu();
    HMENU hFillCircles = CreateMenu();
    HMENU hQuarter = CreateMenu();

    AppendMenu(hQuarter, MF_STRING, 11, _T("Quarter 1"));
    AppendMenu(hQuarter, MF_STRING, 12, _T("Quarter 2"));
    AppendMenu(hQuarter, MF_STRING, 13, _T("Quarter 3"));
    AppendMenu(hQuarter, MF_STRING, 14, _T("Quarter 4"));

    AppendMenu(hFillLines, MF_STRING, 28, _T("Apply Fill"));
    AppendMenu(hFillLines, MF_POPUP, (UINT_PTR)hQuarter, _T("Select Quarter"));

    AppendMenu(hFillCircles, MF_STRING, 29, _T("Apply Fill"));
    AppendMenu(hFillCircles, MF_POPUP, (UINT_PTR)hQuarter, _T("Select Quarter"));

    AppendMenu(hFile, MF_STRING, 1, _T("Save"));
    AppendMenu(hFile, MF_STRING, 2, _T("Load"));
    AppendMenu(hFile, MF_SEPARATOR, 0, NULL);
    AppendMenu(hFile, MF_STRING, 3, _T("Exit"));

    AppendMenu(hDraw, MF_STRING, 4, _T("Line (DDA)"));
    AppendMenu(hDraw, MF_STRING, 5, _T("Line (Midpoint)"));
    AppendMenu(hDraw, MF_STRING, 18, _T("Line (Parametric)"));
    AppendMenu(hDraw, MF_SEPARATOR, 0, NULL);
    AppendMenu(hDraw, MF_STRING, 6, _T("Circle (Midpoint)"));
    AppendMenu(hDraw, MF_STRING, 19, _T("Circle (Direct)"));
    AppendMenu(hDraw, MF_STRING, 20, _T("Circle (Polar)"));
    AppendMenu(hDraw, MF_STRING, 21, _T("Circle (Iterative Polar)"));
    AppendMenu(hDraw, MF_STRING, 22, _T("Circle (Modified Midpoint)"));
    AppendMenu(hDraw, MF_STRING, 32, _T("Ellipse (Direct)"));
    AppendMenu(hDraw, MF_STRING, 42, _T("Ellipse (Polar)"));
    AppendMenu(hDraw, MF_STRING, 43, _T("Ellipse (Midpoint)"));
    AppendMenu(hDraw, MF_STRING, 31, _T("Cardinal Spline"));
    AppendMenu(hDraw, MF_SEPARATOR, 0, NULL);
    AppendMenu(hDraw, MF_STRING, 25, _T("Polygon"));
    AppendMenu(hDraw, MF_STRING, 15, _T("Fill Square (Hermite)"));
    AppendMenu(hDraw, MF_STRING, 16, _T("Fill Rectangle (Bezier)"));
    AppendMenu(hDraw, MF_STRING, 17, _T("Fill Circle (Concentric Circles)"));
    AppendMenu(hDraw, MF_POPUP, (UINT_PTR)hFillLines, _T("Fill Circle with Lines"));
    AppendMenu(hDraw, MF_POPUP, (UINT_PTR)hFillCircles, _T("Fill Circle with Circles in Quarter"));

    AppendMenu(hClipping, MF_STRING, 33, _T("Point (Rectangle Clip)"));
    AppendMenu(hClipping, MF_STRING, 34, _T("Point (Square Clip)"));
    AppendMenu(hClipping, MF_STRING, 39, _T("Point (Circle Clip)"));
    AppendMenu(hClipping, MF_STRING, 30, _T("Line (Rectangle Clip)"));
    AppendMenu(hClipping, MF_STRING, 35, _T("Line (Square Clip)"));
    AppendMenu(hClipping, MF_STRING, 40, _T("Line (Circle Clip)"));
    AppendMenu(hClipping, MF_STRING, 36, _T("Polygon (Rectangle Clip)"));

    AppendMenu(hFill, MF_STRING, 23, _T("Flood Fill (Recursive)"));
    AppendMenu(hFill, MF_STRING, 24, _T("Flood Fill (Non-Recursive)"));
    AppendMenu(hFill, MF_STRING, 26, _T("Convex Fill"));
    AppendMenu(hFill, MF_STRING, 27, _T("Non-Convex Fill"));

    AppendMenu(hOption, MF_STRING, 7, _T("Choose Color"));
    AppendMenu(hOption, MF_STRING, 8, _T("Clear Screen"));
    AppendMenu(hOption, MF_STRING, 37, _T("Toggle Square Clipping Window"));
    AppendMenu(hOption, MF_STRING, 38, _T("Toggle Rectangle Clipping Window"));
    AppendMenu(hOption, MF_STRING, 41, _T("Toggle Circle Clipping Window"));

    AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hFile, _T("File"));
    AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hDraw, _T("Draw"));
    AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hClipping, _T("Clipping"));
    AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hFill, _T("Fill"));
    AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hOption, _T("Options"));

    SetMenu(hwnd, hMenubar);
}

LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
static HDC hdc;

switch (msg)
{
case WM_CREATE:
AddMenus(hwnd);
break;

case WM_COMMAND:
switch (LOWORD(wParam))
{
case 1: // Save
{
TCHAR szFile[260] = _T("shapes.txt");
OPENFILENAME ofn = { sizeof(OPENFILENAME) };
ofn.hwndOwner = hwnd;
ofn.lpstrFilter = _T("Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0");
ofn.lpstrFile = szFile;
ofn.nMaxFile = sizeof(szFile) / sizeof(TCHAR);
ofn.Flags = OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;
ofn.lpstrDefExt = _T("txt");

if (GetSaveFileName(&ofn))
{
SaveShapes(hwnd, szFile);
}
break;
}
case 2: // Load
{
TCHAR szFile[260] = _T("");
OPENFILENAME ofn = { sizeof(OPENFILENAME) };
ofn.hwndOwner = hwnd;
ofn.lpstrFilter = _T("Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0");
ofn.lpstrFile = szFile;
ofn.nMaxFile = sizeof(szFile) / sizeof(TCHAR);
ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
ofn.lpstrDefExt = _T("txt");

if (GetOpenFileName(&ofn))
{
    LoadShapes(hwnd, std::wstring(szFile));

}
break;
}
case 3:
PostQuitMessage(0);
break;
case 4:
currentShape = LINE_DDA;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 5:
currentShape = LINE_MID;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 6:
currentShape = CIRCLE_MID;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 7:
GetCurrentColor();
break;
case 8:
ClearScreen(hwnd);
break;
case 11:
selectedQuarter = 1;
break;
case 12:
selectedQuarter = 2;
break;
case 13:
selectedQuarter = 3;
break;
case 14:
selectedQuarter = 4;
break;
case 15:
hdc = GetDC(hwnd);
FillSquareWithHermite(hdc, 100, 100, 150, currentColor);
ReleaseDC(hwnd, hdc);
break;
case 16:
hdc = GetDC(hwnd);
FillRectangleWithBezier(hdc, 300, 100, 200, 100, currentColor);
ReleaseDC(hwnd, hdc);
break;
case 17:
hdc = GetDC(hwnd);
FillCircle_Circle(hdc, lastCircleX, lastCircleY, lastCircleR, currentColor);
ReleaseDC(hwnd, hdc);
break;
case 18:
currentShape = LINE_PARAMETRIC;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 19:
currentShape = CIRCLE_DIRECT;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 20:
currentShape = CIRCLE_POLAR;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 21:
currentShape = CIRCLE_ITERATIVE_POLAR;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 22:
currentShape = CIRCLE_MODIFIED_MID;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 23:
currentFill = FILL_RECURSIVE;
currentShape = NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 24:
currentFill = FILL_NON_RECURSIVE;
currentShape = NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 25:
currentShape = POLYGON;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 26:
currentFill = FILL_CONVEX;
currentShape = NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 27:
currentFill = FILL_NON_CONVEX;
currentShape = NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 28:
if (hasLastCircle)
{
hdc = GetDC(hwnd);
FillCircle_Lines(hdc, lastCircleX, lastCircleY, lastCircleR, selectedQuarter, currentColor);
ReleaseDC(hwnd, hdc);
}
break;
case 29:
if (hasLastCircle)
{
hdc = GetDC(hwnd);
FillCircleInQuarter(hdc, lastCircleX, lastCircleY, lastCircleR, selectedQuarter, currentColor);
ReleaseDC(hwnd, hdc);
}
break;
case 30:
currentShape = LINE_CLIP_RECT;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
break;
case 31:
currentShape = CARDINAL_SPLINE;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
splinePoints.clear();
break;
case 32:
currentShape = ELLIPSE_DIRECT;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = false;
ellipseSecondClick = false;
break;
case 33:
currentShape = POINT_CLIP_RECT;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
break;
case 34:
currentShape = POINT_CLIP_SQUARE;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
break;
case 35:
currentShape = LINE_CLIP_SQUARE;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
break;
case 36:
currentShape = POLYGON_CLIP_RECT;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
polygonPoints.clear();
break;
case 37:
showClippingWindow = !showClippingWindow;
InvalidateRect(hwnd, NULL, TRUE);
break;
case 38:
showRectangleClippingWindow = !showRectangleClippingWindow;
InvalidateRect(hwnd, NULL, TRUE);
break;
case 39:
currentShape = POINT_CLIP_CIRCLE;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
break;
case 40:
currentShape = LINE_CLIP_CIRCLE;
currentFill = FILL_NONE;
isFirstClick = true;
useClipping = true;
ellipseSecondClick = false;
break;
case 41:
showCircleClippingWindow = !showCircleClippingWindow;
InvalidateRect(hwnd, NULL, TRUE);
break;
case 42:
currentShape = ELLIPSE_POLAR;
currentFill = FILL_NONE;
ellipseSecondClick = false;
break;
case 43:
currentShape = ELLIPSE_MIDPOINT;
break;
}
break;

case WM_LBUTTONDOWN:
{
int x = LOWORD(lParam);
int y = HIWORD(lParam);
POINT pt = { LOWORD(lParam), HIWORD(lParam) };

if (currentShape == ELLIPSE_MIDPOINT) {
if (!ellipseSecondClick) {
ellipsePoints.clear();
ellipsePoints.push_back(pt);
ellipseSecondClick = true;
}
else {
ellipsePoints.push_back(pt);
Shape s;
s.type = ELLIPSE_MIDPOINT;
s.points = ellipsePoints;
s.color = currentColor;
shapes.push_back(s);
ellipseSecondClick = false;
InvalidateRect(hwnd, NULL, FALSE);
}
}

if (currentShape == ELLIPSE_POLAR) {
if (!ellipseSecondClick) {
ellipseCenter.x = LOWORD(lParam);
ellipseCenter.y = HIWORD(lParam);
ellipseSecondClick = true;
}
else {
ellipseRadiusPoint.x = LOWORD(lParam);
ellipseRadiusPoint.y = HIWORD(lParam);


int a = abs(ellipseRadiusPoint.x - ellipseCenter.x);
int b = abs(ellipseRadiusPoint.y - ellipseCenter.y);


Shape s;
s.type = ELLIPSE_POLAR;
s.points.push_back(ellipseCenter);
s.points.push_back(ellipseRadiusPoint);
s.color = currentColor;
s.isPolygon = false;
shapes.push_back(s);


ellipseSecondClick = false;
InvalidateRect(hwnd, NULL, TRUE);
}
break;
}

if (currentFill != FILL_NONE)
{
hdc = GetDC(hwnd);
COLORREF borderColor = RGB(0, 0, 0);

if (currentFill == FILL_RECURSIVE)
{
RECT clientRect;
GetClientRect(hwnd, &clientRect);
FloodFillRecursiveWrapper(hdc, x, y, currentColor, borderColor);
}
else if (currentFill == FILL_NON_RECURSIVE)
{
FloodFillNonRecursive(hdc, x, y, currentColor, borderColor);
}
else if (currentFill == FILL_CONVEX)
{
for (auto it = shapes.rbegin(); it != shapes.rend(); ++it)
{
if (it->isPolygon && it->points.size() >= 3)
{
FillConvexPolygon(hdc, it->points, currentColor);
it->fillType = FILL_CONVEX;
it->fillColor = currentColor;
break;
}
}
}
else if (currentFill == FILL_NON_CONVEX)
{
for (auto it = shapes.rbegin(); it != shapes.rend(); ++it)
{
if (it->isPolygon && it->points.size() >= 3)
{
FillNonConvexPolygon(hdc, it->points, currentColor);
it->fillType = FILL_NON_CONVEX;
it->fillColor = currentColor;
break;
}
}
}
ReleaseDC(hwnd, hdc);
break;
}

if (currentShape == POLYGON || currentShape == POLYGON_CLIP_RECT)
{
if (isFirstClick)
{
polygonPoints.clear();
polygonPoints.push_back({ x, y });
isFirstClick = false;
}
else
{
polygonPoints.push_back({ x, y });
}
InvalidateRect(hwnd, NULL, TRUE);
}
else if (currentShape == ELLIPSE_DIRECT)
{
if (!ellipseSecondClick)
{
ellipseCenter = { x, y };
ellipseSecondClick = true;
hdc = GetDC(hwnd);
DrawPixel(hdc, x, y, currentColor);
ReleaseDC(hwnd, hdc);
}
else
{
shapes.push_back({ ELLIPSE_DIRECT, ellipseCenter, { x, y }, currentColor });
ellipseSecondClick = false;
InvalidateRect(hwnd, NULL, TRUE);
}
}
else if (currentShape == CARDINAL_SPLINE)
{
splinePoints.push_back({ x, y });
InvalidateRect(hwnd, NULL, TRUE);
}
else if (currentShape == POINT_CLIP_RECT || currentShape == POINT_CLIP_SQUARE || currentShape == POINT_CLIP_CIRCLE)
{
hdc = GetDC(hwnd);
if (currentShape == POINT_CLIP_RECT)
{
ClipPointRectangle(hdc, x, y, currentColor);
shapes.push_back({ POINT_CLIP_RECT, { x, y }, { x, y }, currentColor });
}
else if (currentShape == POINT_CLIP_SQUARE)
{
ClipPointSquare(hdc, x, y, currentColor);
shapes.push_back({ POINT_CLIP_SQUARE, { x, y }, { x, y }, currentColor });
}
else
{
ClipPointCircle(hdc, x, y, currentColor);
shapes.push_back({ POINT_CLIP_CIRCLE, { x, y }, { x, y }, currentColor });
}
ReleaseDC(hwnd, hdc);
}
else
{
if (isFirstClick)
{
tempPoint = { x, y };
isFirstClick = false;
}
else
{
hdc = GetDC(hwnd);
if (currentShape == LINE_DDA)
{
if (useClipping)
{
ClipLineRectangle(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_CLIP_RECT, tempPoint, { x, y }, currentColor });
}
else
{
DrawLineDDA(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_DDA, tempPoint, { x, y }, currentColor });
}
}
else if (currentShape == LINE_MID)
{
DrawLineMidpoint(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_MID, tempPoint, { x, y }, currentColor });
}
else if (currentShape == LINE_PARAMETRIC)
{
DrawLineParametric(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_PARAMETRIC, tempPoint, { x, y }, currentColor });
}
else if (currentShape == LINE_CLIP_RECT)
{
ClipLineRectangle(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_CLIP_RECT, tempPoint, { x, y }, currentColor });
}
else if (currentShape == LINE_CLIP_SQUARE)
{
ClipLineSquare(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_CLIP_SQUARE, tempPoint, { x, y }, currentColor });
}
else if (currentShape == LINE_CLIP_CIRCLE)
{
ClipLineCircle(hdc, tempPoint.x, tempPoint.y, x, y, currentColor);
shapes.push_back({ LINE_CLIP_CIRCLE, tempPoint, { x, y }, currentColor });
}
else
{
int r = (int)sqrt(pow(x - tempPoint.x, 2) + pow(y - tempPoint.y, 2));
if (currentShape == CIRCLE_MID)
DrawCircleMidpoint(hdc, tempPoint.x, tempPoint.y, r, currentColor);
else if (currentShape == CIRCLE_DIRECT)
DrawCircleDirect(hdc, tempPoint.x, tempPoint.y, r, currentColor);
else if (currentShape == CIRCLE_POLAR)
DrawCirclePolar(hdc, tempPoint.x, tempPoint.y, r, currentColor);
else if (currentShape == CIRCLE_ITERATIVE_POLAR)
DrawCircleIterativePolar(hdc, tempPoint.x, tempPoint.y, r, currentColor);
else if (currentShape == CIRCLE_MODIFIED_MID)
DrawCircleModifiedMidpoint(hdc, tempPoint.x, tempPoint.y, r, currentColor);
shapes.push_back({ currentShape, tempPoint, { x, y }, currentColor });
lastCircleX = tempPoint.x;
lastCircleY = tempPoint.y;
lastCircleR = r;
hasLastCircle = true;
}
isFirstClick = true;
useClipping = false;
ReleaseDC(hwnd, hdc);
}
}
break;
}

case WM_RBUTTONDOWN:
{
if (currentShape == POLYGON || currentShape == POLYGON_CLIP_RECT)
{
if (polygonPoints.size() >= 3 && !isFirstClick)
{
shapes.push_back({ currentShape, polygonPoints, currentColor, true });
polygonPoints.clear();
isFirstClick = true;
InvalidateRect(hwnd, NULL, TRUE);
}
}
else if (currentShape == CARDINAL_SPLINE && splinePoints.size() >= 4)
{
shapes.push_back({ CARDINAL_SPLINE, splinePoints, currentColor, true });
splinePoints.clear();
InvalidateRect(hwnd, NULL, TRUE);
}
break;
}

case WM_MOUSEMOVE:
{
if ((currentShape == POLYGON || currentShape == POLYGON_CLIP_RECT) && !polygonPoints.empty())
{
InvalidateRect(hwnd, NULL, TRUE);
}
break;
}
case WM_PAINT:
{
PAINTSTRUCT ps;
HDC hdc = BeginPaint(hwnd, &ps);

HBRUSH hWhiteBrush = (HBRUSH)GetStockObject(WHITE_BRUSH);
RECT clientRect;
GetClientRect(hwnd, &clientRect);
FillRect(hdc, &clientRect, hWhiteBrush);

HBRUSH nullBrush = (HBRUSH)GetStockObject(NULL_BRUSH);
HGDIOBJ oldBrush = SelectObject(hdc, nullBrush);



if (showRectangleClippingWindow) {
HPEN rectPen = CreatePen(PS_SOLID, 2, RGB(0, 0, 255));
HGDIOBJ oldPen = SelectObject(hdc, rectPen);
Rectangle(hdc, winLeft, winTop, winRight, winBottom);
SelectObject(hdc, oldPen);
DeleteObject(rectPen);
}

if (showClippingWindow) {
HPEN squarePen = CreatePen(PS_SOLID, 2, RGB(255, 0, 0));
HGDIOBJ oldPen = SelectObject(hdc, squarePen);
Rectangle(hdc, squareLeft, squareTop, squareRight, squareBottom);
SelectObject(hdc, oldPen);
DeleteObject(squarePen);
}

if (showCircleClippingWindow) {
HPEN circlePen = CreatePen(PS_SOLID, 2, RGB(0, 255, 0));
HGDIOBJ oldPen = SelectObject(hdc, circlePen);
Ellipse(hdc, circleCenterX - circleRadius, circleCenterY - circleRadius,
circleCenterX + circleRadius, circleCenterY + circleRadius);
SelectObject(hdc, oldPen);
DeleteObject(circlePen);
}

SelectObject(hdc, oldBrush);

for (const auto& s : shapes) {
if (s.fillType == FILL_CONVEX || s.fillType == FILL_NON_CONVEX) {
if (s.fillType == FILL_CONVEX) {
FillConvexPolygon(hdc, s.points, s.fillColor);
}
else {
FillNonConvexPolygon(hdc, s.points, s.fillColor);
}
}
}

for (const auto& s : shapes) {
if (currentShape == ELLIPSE_POLAR && ellipseSecondClick) {
POINT cursorPos;
GetCursorPos(&cursorPos);
ScreenToClient(hwnd, &cursorPos);
int a = abs(cursorPos.x - ellipseCenter.x);
int b = abs(cursorPos.y - ellipseCenter.y);
DrawEllipse_Polar(hdc, ellipseCenter.x, ellipseCenter.y, a, b, currentColor);
}
else if (s.type == CARDINAL_SPLINE && s.points.size() >= 4) {
DrawCardinalSpline(hdc, s.points, s.color);
}
else if (s.isPolygon) {
if (s.type == POLYGON_CLIP_RECT)
ClipPolygonRectangle(hdc, s.points, s.color);
else
DrawPolygon(hdc, s.points, s.color);
}
else {
switch (s.type) {
case POINT_CLIP_RECT:
ClipPointRectangle(hdc, s.p1.x, s.p1.y, s.color);
break;
case POINT_CLIP_SQUARE:
ClipPointSquare(hdc, s.p1.x, s.p1.y, s.color);
break;
case POINT_CLIP_CIRCLE:
ClipPointCircle(hdc, s.p1.x, s.p1.y, s.color);
break;
case LINE_CLIP_RECT:
ClipLineRectangle(hdc, s.p1.x, s.p1.y, s.p2.x, s.p2.y, s.color);
break;
case LINE_CLIP_SQUARE:
ClipLineSquare(hdc, s.p1.x, s.p1.y, s.p2.x, s.p2.y, s.color);
break;
case LINE_CLIP_CIRCLE:
ClipLineCircle(hdc, s.p1.x, s.p1.y, s.p2.x, s.p2.y, s.color);
break;
case LINE_DDA:
DrawLineDDA(hdc, s.p1.x, s.p1.y, s.p2.x, s.p2.y, s.color);
break;
case LINE_MID:
DrawLineMidpoint(hdc, s.p1.x, s.p1.y, s.p2.x, s.p2.y, s.color);
break;
case LINE_PARAMETRIC:
DrawLineParametric(hdc, s.p1.x, s.p1.y, s.p2.x, s.p2.y, s.color);
break;
case CIRCLE_MID:
case CIRCLE_DIRECT:
case CIRCLE_POLAR:
case CIRCLE_ITERATIVE_POLAR:
case CIRCLE_MODIFIED_MID: {
int r = (int)sqrt(pow(s.p2.x - s.p1.x, 2) + pow(s.p2.y - s.p1.y, 2));
switch (s.type) {
case CIRCLE_MID: DrawCircleMidpoint(hdc, s.p1.x, s.p1.y, r, s.color); break;
case CIRCLE_DIRECT: DrawCircleDirect(hdc, s.p1.x, s.p1.y, r, s.color); break;
case CIRCLE_POLAR: DrawCirclePolar(hdc, s.p1.x, s.p1.y, r, s.color); break;
case CIRCLE_ITERATIVE_POLAR: DrawCircleIterativePolar(hdc, s.p1.x, s.p1.y, r, s.color); break;
case CIRCLE_MODIFIED_MID: DrawCircleModifiedMidpoint(hdc, s.p1.x, s.p1.y, r, s.color); break;
}
break;
}
case ELLIPSE_DIRECT: {
int a = abs(s.p2.x - s.p1.x);
int b = abs(s.p2.y - s.p1.y);
DrawEllipse_Direct(hdc, s.p1.x, s.p1.y, a, b, s.color);
break;
}
case ELLIPSE_MIDPOINT: {
if (s.points.size() >= 2) {
int a = abs(s.points[1].x - s.points[0].x);
int b = abs(s.points[1].y - s.points[0].y);
EllipseMidpoint(hdc, s.points[0].x, s.points[0].y, a, b, s.color);
}
break; }
case ELLIPSE_POLAR:
{
int a = abs(s.points[1].x - s.points[0].x);
int b = abs(s.points[1].y - s.points[0].y);
DrawEllipse_Polar(hdc, s.points[0].x, s.points[0].y, a, b, s.color);
break;

}
}
}
}

if (!polygonPoints.empty()) {
HPEN pen = CreatePen(PS_SOLID, 2, currentColor);
HGDIOBJ oldPen = SelectObject(hdc, pen);

for (size_t i = 0; i < polygonPoints.size(); i++) {
Ellipse(hdc, polygonPoints[i].x - 3, polygonPoints[i].y - 3,
polygonPoints[i].x + 3, polygonPoints[i].y + 3);

if (i > 0) {
MoveToEx(hdc, polygonPoints[i - 1].x, polygonPoints[i - 1].y, NULL);
LineTo(hdc, polygonPoints[i].x, polygonPoints[i].y);
}
}

POINT pt;
GetCursorPos(&pt);
ScreenToClient(hwnd, &pt);
if (polygonPoints.size() > 0) {
MoveToEx(hdc, polygonPoints.back().x, polygonPoints.back().y, NULL);
LineTo(hdc, pt.x, pt.y);
}

SelectObject(hdc, oldPen);
DeleteObject(pen);
}
else if (!splinePoints.empty() && splinePoints.size() >= 4) {
DrawCardinalSpline(hdc, splinePoints, currentColor);
}

EndPaint(hwnd, &ps);
break;
}

case WM_DESTROY:
PostQuitMessage(0);
break;

default:
return DefWindowProc(hwnd, msg, wParam, lParam);
}
return 0;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nCmdShow)
{
WNDCLASSEX wc = { sizeof(WNDCLASSEX) };
wc.style = CS_HREDRAW | CS_VREDRAW;
wc.lpfnWndProc = WindowProcedure;
wc.hInstance = hInstance;
wc.hCursor = LoadCursor(NULL, IDC_CROSS);
wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
wc.lpszClassName = szClassName;

if (!RegisterClassEx(&wc))
{
MessageBox(NULL, _T("Failed to register window class!"), _T("Error"), MB_OK | MB_ICONERROR);
return 0;
}

HWND hwnd = CreateWindow(szClassName, _T("2D Drawing App"),
                         WS_OVERLAPPEDWINDOW, 100, 100, 800, 600, NULL, NULL, hInstance, NULL);

if (hwnd == NULL)
{
MessageBox(NULL, _T("Failed to create window!"), _T("Error"), MB_OK | MB_ICONERROR);
return 0;
}

ShowWindow(hwnd, SW_SHOWNORMAL);
UpdateWindow(hwnd);

MSG msg = { 0 };
while (GetMessage(&msg, NULL, 0, 0))
{
TranslateMessage(&msg);
DispatchMessage(&msg);
}

return (int)msg.wParam;
}
