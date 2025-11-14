#include <math.h>
#include <stdio.h>
#include <stdbool.h>


double hua_DmsToRadians(double dms)
{
    int degrees = (int)dms; // 提取度
    double fractional = dms - degrees; // 获取小数部分
    int minutes = (int)(fractional * 100); // 提取分
    double seconds = (fractional * 100 - minutes) * 100; // 提取秒

    double totalDegrees = degrees + minutes / 60.0 + seconds / 3600.0; // 转换为十进制度
    return totalDegrees * M_PI / 180.0; // 转换为弧度
}
__declspec(dllexport) int Add(int a, int b) {
    return a + b;
}

__declspec(dllexport) bool zs(double xyk, double xyx, double xyy, double xyhd, double xycd, 
            double xyqdr, double xyzdr, double xyzy, double jsk, 
            double jsb, double jd, double result[3]) 
{
    const double epsilon = 0.01; 
jd=hua_DmsToRadians(jd);
    // 第一个条件分支: 圆弧
    if (fabs(xyqdr - xyzdr) < epsilon && fabs(xyqdr) > 0.01) 
    {
        double centerX = xyx + xyqdr * cos(xyhd + xyzy * M_PI / 2.0);
        double centerY = xyy + xyqdr * sin(xyhd + xyzy * M_PI / 2.0);
        double deltaArcLength = jsk - xyk;
        double deltaAngle = deltaArcLength / xyqdr * xyzy;
        double targetAzimuth = xyhd + deltaAngle;

        // 标准化方位角到 [0, 2*PI) 范围
        while (targetAzimuth < 0) targetAzimuth += 2.0 * M_PI;
        while (targetAzimuth >= 2.0 * M_PI) targetAzimuth -= 2.0 * M_PI;

        double angle = xyhd + xyzy * M_PI / 2.0 + deltaAngle;
        double targetX = centerX - xyqdr * cos(angle) + jsb * cos(angle - xyzy * M_PI / 2.0 + jd);
        double targetY = centerY - xyqdr * sin(angle) + jsb * sin(angle - xyzy * M_PI / 2.0 + jd);

        result[0] = targetX;
        result[1] = targetY;
        result[2] = targetAzimuth;
        return true;
    }

    // 第二个条件分支: 直线
    if (xyqdr < epsilon && xyzdr < epsilon && fabs(xyzy) < epsilon) 
    {
        double targetX = xyx + (jsk - xyk) * cos(xyhd) + jsb * cos(xyhd + jd);
        double targetY = xyy + (jsk - xyk) * sin(xyhd) + jsb * sin(xyhd + jd);

        result[0] = targetX;
        result[1] = targetY;
        result[2] = xyhd; // 直线情况下，方位角保持不变
        return true;
    }

    // 第三个条件分支: 缓和曲线 (默认情况)
    double effectiveQdr = (xyqdr < 0.001) ? 99999999.0 : xyqdr;
    double effectiveZdr = (xyzdr < 0.001) ? 99999999.0 : xyzdr;

    double f0 = xyhd;
    double q = xyzy;
    double c = 1.0 / effectiveQdr;
    double d = (effectiveQdr - effectiveZdr) / (2.0 * xycd * effectiveQdr * effectiveZdr);

    // 高斯-勒让德积分系数 (索引 0-3)
    double rr[4] = {0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226};
    double vv[4] = {0.0694318442, 0.3300094782, 0.6699905218, 0.9305681558};

    double w = jsk - xyk;
    double xs = 0.0;
    double ys = 0.0;

    for (int i = 0; i < 4; i++) 
    {
        double ff = f0 + q * vv[i] * w * (c + vv[i] * w * d);
        xs += rr[i] * cos(ff);
        ys += rr[i] * sin(ff);
    }

    double fhz3 = f0 + q * w * (c + w * d);
    // 标准化最终方位角 fhz3 到 [0, 2*PI) 范围
    while (fhz3 < 0) fhz3 += 2.0 * M_PI;
    while (fhz3 >= 2.0 * M_PI) fhz3 -= 2.0 * M_PI;

    double fhz1 = xyx + w * xs + jsb * cos(fhz3 + jd);
    double fhz2 = xyy + w * ys + jsb * sin(fhz3 + jd);

    result[0] = fhz1;
    result[1] = fhz2;
    result[2] = fhz3;
    return true;
}