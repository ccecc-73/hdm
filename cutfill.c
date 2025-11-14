#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define EPSILON 1e-12

// 点结构体
typedef struct {
    double x;
    double y;
} Point;

// 交点结构体
typedef struct {
    double x;
    double y;
    int sji; // 设计线索引
    int dmj; // 地面线索引
} IntersectionPoint;

// 扩展线段函数
void extendLine(double x1, double y1, double x2, double y2, double distance, double* result) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double length = sqrt(dx * dx + dy * dy);
    
    if (length < EPSILON) {
        result[0] = x1;
        result[1] = y1;
        return;
    }
    
    double unitX = dx / length;
    double unitY = dy / length;
    
    result[0] = x1 + unitX * distance;
    result[1] = y1 + unitY * distance;
}

// 线段交点检测函数
bool intersectionsSegment(double x1, double y1, double x2, double y2,
                         double x3, double y3, double x4, double y4, double* result) {
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    
    if (fabs(denom) < EPSILON) {
        return false; // 平行或重合
    }
    
    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
    double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;
    
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        result[0] = x1 + t * (x2 - x1);
        result[1] = y1 + t * (y2 - y1);
        return true;
    }
    
    return false;
}

// 计算多边形面积（鞋带公式）
double calculatePolygonArea(Point* points, int pointCount) {
    double area = 0.0;
    
    for (int i = 0; i < pointCount - 1; i++) {
        area += points[i].x * points[i + 1].y - points[i].y * points[i + 1].x;
    }
    
    return fabs(area) / 2.0;
}

// 计算带符号的多边形面积
double calculateSignedArea(Point* points, int pointCount) {
    double area = 0.0;
    
    for (int i = 0; i < pointCount - 1; i++) {
        area += points[i].x * points[i + 1].y - points[i].y * points[i + 1].x;
    }
    
    return area / 2.0;
}

// 主函数：计算填挖面积
__declspec(dllexport) void cutAndFillArea(double* xsA, double* ysA, int lengthA,
                   double* xsB, double* ysB, int lengthB,
                   double epsilon, double* result) {
    
    // 扩展线段端点
    double temp[2];
    
    extendLine(xsA[0], ysA[0], xsA[1], ysA[1], -epsilon, temp);
    xsA[0] = temp[0];
    ysA[0] = temp[1];
    
    extendLine(xsB[0], ysB[0], xsB[1], ysB[1], -epsilon, temp);
    xsB[0] = temp[0];
    ysB[0] = temp[1];
    
    extendLine(xsA[lengthA - 1], ysA[lengthA - 1], 
               xsA[lengthA - 2], ysA[lengthA - 2], -epsilon, temp);
    xsA[lengthA - 1] = temp[0];
    ysA[lengthA - 1] = temp[1];
    
    extendLine(xsB[lengthB - 1], ysB[lengthB - 1], 
               xsB[lengthB - 2], ysB[lengthB - 2], -epsilon, temp);
    xsB[lengthB - 1] = temp[0];
    ysB[lengthB - 1] = temp[1];
    
    // 存储交点信息
    IntersectionPoint* xys = malloc(100 * sizeof(IntersectionPoint));
    int intersectionCount = 0;
    double intersectionResult[2];
    
    double fill = 0;
    double cut = 0;
    
    // 查找所有交点
    for (int i = 0; i < lengthA - 1; i++) {
        double x1 = xsA[i], y1 = ysA[i];
        double x2 = xsA[i + 1], y2 = ysA[i + 1];
        
        for (int j = 0; j < lengthB - 1; j++) {
            double x3 = xsB[j], y3 = ysB[j];
            double x4 = xsB[j + 1], y4 = ysB[j + 1];
            
            if (intersectionsSegment(x1, y1, x2, y2, x3, y3, x4, y4, intersectionResult)) {
                if (intersectionCount < 100) {
                    xys[intersectionCount].x = intersectionResult[0];
                    xys[intersectionCount].y = intersectionResult[1];
                    xys[intersectionCount].sji = i;
                    xys[intersectionCount].dmj = j;
                    intersectionCount++;
                }
            }
        }
    }
    
    // 处理交点序列，计算填挖面积
    for (int i = 0; i < intersectionCount - 1; i++) {
        IntersectionPoint xy0 = xys[i];
        IntersectionPoint xy1 = xys[i + 1];
        int i0 = xy0.sji, i1 = xy1.sji;
        int j0 = xy0.dmj, j1 = xy1.dmj;
        
        // 估计最大点数并分配内存
        int maxPoints = (i1 - i0 + 1) + (j1 - j0 + 1) + 2;
        Point* pts = malloc(maxPoints * sizeof(Point));
        int pointCount = 0;
        
        // 添加交点
        pts[pointCount++] = (Point){xy0.x, xy0.y};
        
        // 添加设计线（sj）上的点
        for (int j = i0 + 1; j <= i1; j++) {
            if (pointCount < maxPoints) {
                pts[pointCount++] = (Point){xsA[j], ysA[j]};
            }
        }
        
        // 添加第二个交点
        if (pointCount < maxPoints) {
            pts[pointCount++] = (Point){xy1.x, xy1.y};
        }
        
        // 添加地面线（dm）上的点（逆序）
        for (int j = j1; j > j0; j--) {
            if (pointCount < maxPoints) {
                pts[pointCount++] = (Point){xsB[j], ysB[j]};
            }
        }
        
        // 闭合多边形
        if (pointCount < maxPoints) {
            pts[pointCount++] = (Point){xy0.x, xy0.y};
        }
        
        // 计算面积
        double area = calculatePolygonArea(pts, pointCount);
        double signedArea = calculateSignedArea(pts, pointCount);
        
        if (signedArea < 0) {
            fill += area; // 负面积为填方
        } else {
            cut += area; // 正面积为挖方
        }
        
        free(pts);
    }
    
    // 四舍五入到4位小数
    result[0] = round(fill * 10000) / 10000.0;
    result[1] = round(cut * 10000) / 10000.0;
    
    free(xys);
}
