import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class LL {
/**
 * 
 * @param lon
 * @param lat
 * @return
 */
    public static double[] web_proj(double lon, double lat) {
        double x = lon * 20037508.34 / 180;
        double y = Math.log(Math.tan((90 + lat) * Math.PI / 360)) / (Math.PI / 180);
        y = y * 20037508.34 / 180;
        return new double[] { x, y };
    }

    public static double[] web_unproj(double x, double y) {
        double lon = x / 20037508.34 * 180;
        double lat = 180 / Math.PI * (2 * Math.atan(Math.exp(y / 20037508.34 * Math.PI * 180)) - Math.PI / 2);
        return new double[] { lon, lat };
    }

    public static int[] web_getTileFromLatLonZ(double lon, double lat, int zoom) {

        int n = (int) Math.pow(2, zoom);
        int x = (int) Math.floor((lon + 180) / 360 * n);
        int y = (int) Math.floor(
                (1 - Math.log(Math.tan(lat * Math.PI / 180) + 1 / Math.cos(lat * Math.PI / 180)) / Math.PI) / 2 * n);
        return new int[] { zoom, x, y };
    }

    public static double[] getAbsolutePixelFromLonLat(double lon, double lat, int zoom, int delx, int dely) {
        final int tileSize = 256;
        int numTilesAtZoom = (int) Math.pow(2, zoom);

        // 将经纬度转换为瓦片坐标(行 列 编号)
        int tileX = (int) Math.floor((lon + 180) / 360 * numTilesAtZoom);
        double latRad = lat * Math.PI / 180;
        int tileY = (int) Math
                .floor((1 - Math.log(Math.tan(latRad) + 1 / Math.cos(latRad)) / Math.PI) / 2 * numTilesAtZoom);

        // 计算像素坐标在瓦片内的偏移量
        double pixelXInTile = ((lon + 180) / 360 * numTilesAtZoom - tileX) * tileSize + delx * 256;
        double pixelYInTile = ((1 - Math.log(Math.tan(latRad) + 1 / Math.cos(latRad)) / Math.PI) / 2 * numTilesAtZoom
                - tileY) * tileSize + dely * 256;

        // 计算绝对像素坐标
        // double absolutePixelX = tileX * tileSize + pixelXInTile;
        // double absolutePixelY = (tileY + dely) * tileSize + pixelYInTile;

        return new double[] { pixelXInTile, pixelYInTile };
    }

    public static String num2K(double meters) {
        int km = (int) Math.floor(meters / 1000);
        double m = meters - km * 1000;
        m = Math.round(m * 100) / 100.0; // 四舍五入到两位小数

        if (m == (int) m) {
            // 整数米，如 2.00 → K1+002
            return "K" + km + "+" + String.format("%03d", (int) m);
        } else {
            // 小数米，如 2.35 → K0+002.35
            // 手动格式化小数部分
            int integerPart = (int) m;
            double fractionalPart = m - integerPart;
            // 将小数部分转换为整数形式（如0.35变成35），并格式化为两位整数
            String fractionalStr = String.format("%02d", (int) Math.round(fractionalPart * 100));
            return "K" + km + "+" + String.format("%03d", integerPart) + "." + fractionalStr;
        }
    }

    public static double dmsToRadians(double dms) {
        int degrees = (int) dms;
        double fractional = dms - degrees;
        int minutes = (int) (fractional * 100);
        double seconds = (fractional * 100 - minutes) * 100;

        double totalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;
        return totalDegrees * Math.PI / 180;
    }

    public static double radiansToDMS(double radians) {
        double degrees = radians * (180 / Math.PI);
        int d = (int) Math.floor(degrees);
        double remaining = (degrees - d) * 60;
        int m = (int) Math.floor(remaining);
        double s = (remaining - m) * 60;

        // 手动格式化分和秒，确保两位数格式
        String mm = m < 10 ? "0" + m : String.valueOf(m);
        // 秒保留一位小数，并格式化为两位数
        String ss;
        if (s < 10) {
            // 使用String.format确保一位小数，并用0填充
            ss = "0" + String.format("%.1f", s);
        } else {
            ss = String.format("%.1f", s);
        }

        String dfm = d + mm + ss;
        // 将组合字符串解析为double，并除以10000，最后四舍五入到5位小数
        double result = Double.parseDouble(dfm) / 10000.0;

        return result;
    }

    public static String radiansToDMS_度分秒(double radians) {
        double degrees = radians * (180 / Math.PI);
        int d = (int) Math.floor(degrees);
        double remaining = (degrees - d) * 60;
        int m = (int) Math.floor(remaining);
        double s = (remaining - m) * 60;

        // 格式化分，确保两位数
        String mm = m < 10 ? "0" + m : String.valueOf(m);
        // 格式化秒，保留一位小数并确保两位数
        String ss;
        if (s < 10) {
            // 使用String.format控制小数位数，并手动补零
            ss = "0" + String.format("%.1f", s);
        } else {
            ss = String.format("%.1f", s);
        }

        return d + "°" + mm + "′" + ss + "″";
    }

    public static double[] fwj(double x0, double y0, double x1, double y1) {
        double x = x1 - x0;
        double y = y1 - y0;
        double cd = Math.sqrt(x * x + y * y); // 距离
        double hd = Math.atan2(y, x); // 角度（弧度）

        if (hd < 0) {
            hd += 2 * Math.PI; // 确保角度在 [0, 2π) 范围内
        }
        return new double[] { cd, hd };
    }

    public static double[] zs(double xyk, double xyx, double xyy, double xyhd, double xycd, double xyqdr, double xyzdr,
            double xyzy, double jsk, double jsb, double jd) {
        jd = dmsToRadians(jd); // 假设已存在此方法

        // 第一个条件分支: 圆弧 (起点半径和终点半径近似相等且不为零)
        if (Math.abs(xyqdr - xyzdr) < 0.01 && xyqdr != 0) {
            double centerX = xyx + xyqdr * Math.cos(xyhd + xyzy * Math.PI / 2);
            double centerY = xyy + xyqdr * Math.sin(xyhd + xyzy * Math.PI / 2);
            double deltaArcLength = jsk - xyk;
            double deltaAngle = deltaArcLength / xyqdr * xyzy;
            double targetAzimuth = xyhd + deltaAngle;
            if (targetAzimuth < 0) {
                targetAzimuth += 2 * Math.PI;
            }
            double angle = xyhd + xyzy * Math.PI / 2 + deltaAngle;
            double targetX = centerX - xyqdr * Math.cos(angle) + jsb * Math.cos(angle - xyzy * Math.PI / 2 + jd);
            double targetY = centerY - xyqdr * Math.sin(angle) + jsb * Math.sin(angle - xyzy * Math.PI / 2 + jd);

            return new double[] { targetX, targetY, targetAzimuth };
        }

        // 第二个条件分支: 直线 (起点半径、终点半径和转向都近似为零)
        if (xyqdr < 0.01 && xyzdr < 0.01 && Math.abs(xyzy) < 0.01) {
            double targetX = xyx + (jsk - xyk) * Math.cos(xyhd) + jsb * Math.cos(xyhd + jd);
            double targetY = xyy + (jsk - xyk) * Math.sin(xyhd) + jsb * Math.sin(xyhd + jd);
            return new double[] { targetX, targetY, xyhd };
        }

        // 第三个条件分支: 缓和曲线 (默认情况)
        if (xyqdr < 0.001)
            xyqdr = 99999999;
        if (xyzdr < 0.001)
            xyzdr = 99999999;

        double f0 = xyhd;
        double q = xyzy;
        double c = 1 / xyqdr;
        double d = (xyqdr - xyzdr) / (2 * xycd * xyqdr * xyzdr);

        // 初始化 rr 和 vv 数组，Java数组索引从0开始
        double[] rr = new double[4]; // 大小为4，索引0-3
        double[] vv = new double[4];

        rr[0] = 0.1739274226;
        rr[1] = 0.3260725774;
        rr[2] = rr[1]; // 0.3260725774
        rr[3] = rr[0]; // 0.1739274226

        vv[0] = 0.0694318442;
        vv[1] = 0.3300094782;
        vv[2] = 1 - vv[1]; // 0.6699905218
        vv[3] = 1 - vv[0]; // 0.9305681558

        double w = jsk - xyk;
        double xs = 0;
        double ys = 0;

        // 循环索引改为从0到3
        for (int i = 0; i < 4; i++) {
            double ff = f0 + q * vv[i] * w * (c + vv[i] * w * d);
            xs += rr[i] * Math.cos(ff);
            ys += rr[i] * Math.sin(ff);
        }

        double fhz3 = f0 + q * w * (c + w * d);
        // 规范化角度到 [0, 2π) 范围
        fhz3 = fhz3 % (2 * Math.PI);
        if (fhz3 < 0) {
            fhz3 += 2 * Math.PI;
        }

        double fhz1 = xyx + w * xs + jsb * Math.cos(fhz3 + jd);
        double fhz2 = xyy + w * ys + jsb * Math.sin(fhz3 + jd);
        return new double[] { fhz1, fhz2, fhz3 };
    }

    public static double[] dantiaoxianludange(double[][] pqx, double k, double b, double z) {
        for (int i = 0; i < pqx.length; i++) {
            double dtk = pqx[i][0];
            double dtx = pqx[i][1];
            double dty = pqx[i][2];
            double dtfwj = pqx[i][3];
            double dtcd = pqx[i][4];
            double dtr1 = pqx[i][5];
            double dtr2 = pqx[i][6];
            double dtzy = pqx[i][7];

            if (k >= dtk && k <= dtk + dtcd) {
                double hudu = dmsToRadians(dtfwj); // 假设已存在此方法
                return zs(dtk, dtx, dty, hudu, dtcd, dtr1, dtr2, dtzy, k, b, z); // 假设已存在此方法
            }
        }
        return new double[] { 0, 0, 0 };
    }

    public static double[] fs(double fsx, double fsy, double[][] pqx) {
        // 调用fwj方法，返回double[]，不再使用dynamic
        double[] jljd = fwj(pqx[0][1], pqx[0][2], fsx, fsy);
        double k = pqx[0][0];
        double hudu = dmsToRadians(pqx[0][3]); // 调用度分秒转弧度方法
        double cz = jljd[0] * Math.cos(jljd[1] - hudu);
        double pj = jljd[0] * Math.sin(jljd[1] - hudu);

        int hang = pqx.length - 1;
        double qdlc = pqx[0][0];
        double zdlc = pqx[hang][0] + pqx[hang][4];
        int jisuancishu = 0;

        while (Math.abs(cz) > 0.01) {
            k += cz;
            jisuancishu += 1;

            if (k < qdlc) {
                return new double[] { -1, -1 };
            }
            if (k > zdlc) {
                return new double[] { -2, -2 };
            }
            if (jisuancishu > 15) {
                return new double[] { -3, -3 };
            }

            double[] xy = dantiaoxianludange(pqx, k, 0, 0); // 调用线元计算坐标方法
            jljd = fwj(xy[0], xy[1], fsx, fsy);
            cz = jljd[0] * Math.cos(jljd[1] - xy[2]);
            pj = jljd[0] * Math.sin(jljd[1] - xy[2]);
        }

        // 使用乘除运算实现四舍五入到3位小数，避免使用BigDecimal
        double roundedK = Math.round(k * 1000.0) / 1000.0;
        double roundedPj = Math.round(pj * 1000.0) / 1000.0;

        return new double[] { roundedK, roundedPj };
    }

    private static double gaocheng(double bpdlc, double bpdgc, double r, double qp, double hp, double t, double k) {
        double f = qp - hp;
        r = r * Math.abs(f) / f; // 使用Math.abs取绝对值
        double x;
        if (k <= bpdlc - t) {
            x = 0;
        } else if (k >= bpdlc + t) {
            x = 0;
            qp = hp;
        } else {
            x = k - bpdlc + t;
        }
        // 使用Math.pow进行幂运算
        return (bpdgc - (bpdlc - k) * qp - Math.pow(x, 2) / 2 / r);
    }

    public static double h(double k, double[][] sqxb) {
        double hp = 0;
        int length = sqxb.length; // 使用.length获取数组长度

        for (int i = 1; i < length - 1; i++) {
            double r = sqxb[i][2]; // Java使用[i][j]访问二维数组元素
            if (r < 0.001) {
                r = 0.001;
            }

            double qp = (sqxb[i][1] - sqxb[i - 1][1]) / (sqxb[i][0] - sqxb[i - 1][0]);
            hp = (sqxb[i + 1][1] - sqxb[i][1]) / (sqxb[i + 1][0] - sqxb[i][0]);
            double f = qp - hp;
            double t = r * Math.abs(f) / 2; // 使用Math.abs取绝对值

            if (k <= sqxb[i][0] + t) {
                double result = gaocheng(sqxb[i][0], sqxb[i][1], r, qp, hp, t, k); // 调用gaocheng方法
                // 使用乘除运算实现四舍五入到3位小数
                return Math.round(result * 1000.0) / 1000.0;
            }
        }

        // 处理最后一个元素
        if (k <= sqxb[length - 1][0]) {
            double linearResult = sqxb[length - 1][1] + (k - sqxb[length - 1][0] * hp);
            return Math.round(linearResult * 1000.0) / 1000.0;
        }

        return -1;
    }

    public static double[] gauss_proj(double L, double B, double lonCenter) {
        double pi = 3.141592653589793238463;
        double p0 = 206264.8062470963551564;
        double e = 0.00669438002290;
        double e1 = 0.00673949677548;
        double b = 6356752.3141;
        double a = 6378137.0;

        // 将度转换为弧度
        B = B * pi / 180;
        L = L * pi / 180;

        double L_num;
        double L_center;

        // 确定中央子午线经度
        if (lonCenter >= 359) {
            L_num = Math.floor(L * 180 / pi / 3.0 + 0.5);
            L_center = 3 * L_num;
        } else {
            L_center = lonCenter;
        }

        // 计算经度与中央子午线的差值（以秒为单位）
        double l = (L * 180 / pi - L_center) * 3600;

        // 计算子午线弧长公式系数
        double M0 = a * (1 - e);
        double M2 = 3.0 / 2.0 * e * M0;
        double M4 = 5.0 / 4.0 * e * M2;
        double M6 = 7.0 / 6.0 * e * M4;
        double M8 = 9.0 / 8.0 * e * M6;

        double a0 = M0 + M2 / 2.0 + 3.0 / 8.0 * M4 + 5.0 / 16.0 * M6 + 35.0 / 128.0 * M8;
        double a2 = M2 / 2.0 + M4 / 2 + 15.0 / 32.0 * M6 + 7.0 / 16.0 * M8;
        double a4 = M4 / 8.0 + 3.0 / 16.0 * M6 + 7.0 / 32.0 * M8;
        double a6 = M6 / 32.0 + M8 / 16.0;
        double a8 = M8 / 128.0;

        // 计算子午线弧长Xz
        double Xz = a0 * B - a2 / 2.0 * Math.sin(2 * B) + a4 / 4.0 * Math.sin(4 * B) - a6 / 6.0 * Math.sin(6 * B)
                + a8 / 8.0 * Math.sin(8 * B);

        // 计算高斯投影公式中的各个参数
        double c = a * a / b;
        double V = Math.sqrt(1 + e1 * Math.cos(B) * Math.cos(B));
        double N = c / V;
        double t = Math.tan(B);
        double n = e1 * Math.cos(B) * Math.cos(B);

        // 计算高斯投影系数m1-m6
        double m1 = N * Math.cos(B);
        double m2 = N / 2.0 * Math.sin(B) * Math.cos(B);
        double m3 = N / 6.0 * Math.pow(Math.cos(B), 3) * (1 - t * t + n);
        double m4 = N / 24.0 * Math.sin(B) * Math.pow(Math.cos(B), 3) * (5 - t * t + 9 * n);
        double m5 = N / 120.0 * Math.pow(Math.cos(B), 5) * (5 - 18 * t * t + Math.pow(t, 4) + 14 * n - 58 * n * t * t);
        double m6 = N / 720.0 * Math.sin(B) * Math.pow(Math.cos(B), 5) * (61 - 58 * t * t + Math.pow(t, 4));

        // 计算高斯平面坐标x和y
        double x = Xz + m2 * l * l / Math.pow(p0, 2) + m4 * Math.pow(l, 4) / Math.pow(p0, 4)
                + m6 * Math.pow(l, 6) / Math.pow(p0, 6);
        double y0 = m1 * l / p0 + m3 * Math.pow(l, 3) / Math.pow(p0, 3) + m5 * Math.pow(l, 5) / Math.pow(p0, 5);
        double y = y0 + 500000; // 添加500km假东偏移

        return new double[] { x, y, L_center };
    }

    public static double[] gauss_unproj(double x, double y, double l0) {
        double pi = 3.141592653589793238463;
        double e = 0.00669438002290;
        double e1 = 0.00673949677548;
        double b = 6356752.3141;
        double a = 6378137.0;

        // 移除500000米的假东偏移
        double y1 = y - 500000;

        // 计算子午线弧长公式系数
        double M0 = a * (1 - e);
        double M2 = 3.0 / 2.0 * e * M0;
        double M4 = 5.0 / 4.0 * e * M2;
        double M6 = 7.0 / 6.0 * e * M4;
        double M8 = 9.0 / 8.0 * e * M6;

        double a0 = M0 + M2 / 2.0 + 3.0 / 8.0 * M4 + 5.0 / 16.0 * M6 + 35.0 / 128.0 * M8;
        double a2 = M2 / 2.0 + M4 / 2 + 15.0 / 32.0 * M6 + 7.0 / 16.0 * M8;
        double a4 = M4 / 8.0 + 3.0 / 16.0 * M6 + 7.0 / 32.0 * M8;
        double a6 = M6 / 32.0 + M8 / 16.0;

        // 迭代计算底点纬度Bf
        double Bf = x / a0;
        double B0 = Bf;

        // 使用do-while循环确保至少执行一次迭代
        do {
            B0 = Bf;
            double FBf = -a2 / 2.0 * Math.sin(2 * B0) + a4 / 4.0 * Math.sin(4 * B0) - a6 / 6.0 * Math.sin(6 * B0);
            Bf = (x - FBf) / a0;
        } while (Math.abs(Bf - B0) > 0.0000001);

        // 计算高斯投影反算公式中的各个参数
        double t = Math.tan(Bf);
        double c = a * a / b;
        double V = Math.sqrt(1 + e1 * Math.cos(Bf) * Math.cos(Bf));
        double N = c / V;
        double M = c / Math.pow(V, 3);
        double n = e1 * Math.cos(Bf) * Math.cos(Bf);

        // 计算高斯投影反算系数n1-n6
        double n1 = 1 / (N * Math.cos(Bf));
        double n2 = -t / (2.0 * M * N);
        double n3 = -(1 + 2 * t * t + n) / (6.0 * Math.pow(N, 3) * Math.cos(Bf));
        double n4 = t * (5 + 3 * t * t + n - 9 * n * t * t) / (24.0 * M * Math.pow(N, 3));
        double n5 = (5 + 28 * t * t + 24 * Math.pow(t, 4) + 6 * n + 8 * n * t * t)
                / (120.0 * Math.pow(N, 5) * Math.cos(Bf));
        double n6 = -t * (61 + 90 * t * t + 45 * Math.pow(t, 4)) / (720.0 * M * Math.pow(N, 5));

        // 计算纬度B（度）和经度L（度）
        double B = (Bf + n2 * y1 * y1 + n4 * Math.pow(y1, 4) + n6 * Math.pow(y1, 6)) / pi * 180;
        double L0 = l0;
        double l = n1 * y1 + n3 * Math.pow(y1, 3) + n5 * Math.pow(y1, 5);
        double L = L0 + l / pi * 180;

        return new double[] { L, B };
    }

    public static double[] utm_proj(double longitude, double latitude) {
        double EQUATORIAL_RADIUS = 6378137.0;
        double FLATTENING = 1 / 298.257223563;
        double ECC_SQUARED = 2 * FLATTENING - Math.pow(FLATTENING, 2);
        double ECC_PRIME_SQUARED = ECC_SQUARED / (1 - ECC_SQUARED);
        double SCALE_FACTOR = 0.9996;
        double FALSE_EASTING = 500000.0;
        double FALSE_NORTHING_S = 10000000.0;

        double centralMeridian = 0;
        int zoneNumber = (int) Math.floor((longitude + 180) / 6) + 1;
        centralMeridian = (zoneNumber - 1) * 6 - 180 + 3;

        double latRad = latitude * Math.PI / 180.0;
        double lonRad = longitude * Math.PI / 180.0;
        double lonCenterRad = centralMeridian * Math.PI / 180.0;

        double N = EQUATORIAL_RADIUS / Math.sqrt(1 - ECC_SQUARED * Math.pow(Math.sin(latRad), 2));
        double T = Math.pow(Math.tan(latRad), 2);
        double C = ECC_PRIME_SQUARED * Math.pow(Math.cos(latRad), 2);
        double A = (lonRad - lonCenterRad) * Math.cos(latRad);

        double M = EQUATORIAL_RADIUS
                * ((1 - ECC_SQUARED / 4 - 3 * Math.pow(ECC_SQUARED, 2) / 64 - 5 * Math.pow(ECC_SQUARED, 3) / 256)
                        * latRad -
                        (3 * ECC_SQUARED / 8 + 3 * Math.pow(ECC_SQUARED, 2) / 32 + 45 * Math.pow(ECC_SQUARED, 3) / 1024)
                                * Math.sin(2 * latRad)
                        +
                        (15 * Math.pow(ECC_SQUARED, 2) / 256 + 45 * Math.pow(ECC_SQUARED, 3) / 1024)
                                * Math.sin(4 * latRad)
                        -
                        (35 * Math.pow(ECC_SQUARED, 3) / 3072) * Math.sin(6 * latRad));

        double easting = SCALE_FACTOR * N * (A +
                (1 - T + C) * Math.pow(A, 3) / 6 +
                (5 - 18 * T + Math.pow(T, 2) + 72 * C - 58 * ECC_PRIME_SQUARED) * Math.pow(A, 5) / 120) + FALSE_EASTING;

        double northing = SCALE_FACTOR * (M +
                N * Math.tan(latRad) * (Math.pow(A, 2) / 2 +
                        (5 - T + 9 * C + 4 * Math.pow(C, 2)) * Math.pow(A, 4) / 24 +
                        (61 - 58 * T + Math.pow(T, 2) + 600 * C - 330 * ECC_PRIME_SQUARED) * Math.pow(A, 6) / 720));

        if (latitude < 0) {
            northing += FALSE_NORTHING_S;
        }

        return new double[] { northing, easting, zoneNumber };
    }

    public static double[] utm_unproj(double northing, double easting, boolean isNorthern, int zoneNumber) {
        double EQUATORIAL_RADIUS = 6378137.0;
        double FLATTENING = 1 / 298.257223563;
        double ECC_SQUARED = 2 * FLATTENING - Math.pow(FLATTENING, 2);
        double ECC_PRIME_SQUARED = ECC_SQUARED / (1 - ECC_SQUARED);
        double SCALE_FACTOR = 0.9996;
        double FALSE_EASTING = 500000.0;
        double FALSE_NORTHING_S = 10000000.0;

        // 移除假偏移
        double x = easting - FALSE_EASTING;
        double y = isNorthern ? northing : northing - FALSE_NORTHING_S;

        // 计算中央子午线经度
        double centralMeridian = (zoneNumber - 1) * 6 - 180 + 3;
        double lonCenterRad = centralMeridian * Math.PI / 180.0;

        // 计算底点纬度
        double M = y / SCALE_FACTOR;
        double mu = M / (EQUATORIAL_RADIUS
                * (1 - ECC_SQUARED / 4 - 3 * Math.pow(ECC_SQUARED, 2) / 64.0 - 5 * Math.pow(ECC_SQUARED, 3) / 256.0));

        double e1 = (1 - Math.sqrt(1 - ECC_SQUARED)) / (1 + Math.sqrt(1 - ECC_SQUARED));

        // 使用迭代法计算底点纬度phi1Rad
        double phi1Rad = mu + (3 * e1 / 2 - 27 * Math.pow(e1, 3) / 32) * Math.sin(2 * mu) +
                (21 * Math.pow(e1, 2) / 16 - 55 * Math.pow(e1, 4) / 32) * Math.sin(4 * mu) +
                (151 * Math.pow(e1, 3) / 96) * Math.sin(6 * mu);

        // 计算辅助参数
        double N1 = EQUATORIAL_RADIUS / Math.sqrt(1 - ECC_SQUARED * Math.pow(Math.sin(phi1Rad), 2));
        double T1 = Math.pow(Math.tan(phi1Rad), 2);
        double C1 = ECC_PRIME_SQUARED * Math.pow(Math.cos(phi1Rad), 2);
        double R1 = EQUATORIAL_RADIUS * (1 - ECC_SQUARED)
                / Math.pow(1 - ECC_SQUARED * Math.pow(Math.sin(phi1Rad), 2), 1.5);

        double D = x / (N1 * SCALE_FACTOR);

        // 计算纬度和经度（弧度）
        double latRad = phi1Rad - (N1 * Math.tan(phi1Rad) / R1) *
                (Math.pow(D, 2) / 2 -
                        (5 + 3 * T1 + 10 * C1 - 4 * Math.pow(C1, 2) - 9 * ECC_PRIME_SQUARED) * Math.pow(D, 4) / 24 +
                        (61 + 90 * T1 + 298 * C1 + 45 * Math.pow(T1, 2) - 252 * ECC_PRIME_SQUARED - 3 * Math.pow(C1, 2))
                                * Math.pow(D, 6) / 720);

        double lonRad = lonCenterRad + (D - (1 + 2 * T1 + C1) * Math.pow(D, 3) / 6 +
                (5 - 2 * C1 + 28 * T1 - 3 * Math.pow(C1, 2) + 8 * ECC_PRIME_SQUARED + 24 * Math.pow(T1, 2))
                        * Math.pow(D, 5) / 120)
                / Math.cos(phi1Rad);

        // 将弧度转换为度
        double longitude = lonRad * 180 / Math.PI;
        double latitude = latRad * 180 / Math.PI;

        return new double[] { longitude, latitude };
    }

    public static int Utm_zone(double longitude) {
        return (int) ((longitude + 180) / 6) + 1;
    }

    public static double[] cs4(double[] source, double[] target) {
        // 参数校验
        if (source == null || target == null || source.length != target.length) {
            System.out.println("坐标数组长度必须相等");
            return new double[] { 0, 0, 0, 1 };
        }
        if (source.length < 4 || source.length % 2 != 0) {
            System.out.println("至少需要2个点且坐标为偶数");
            return new double[] { 0, 0, 0, 1 };
        }

        int pointCount = source.length / 2;

        // 计算重心坐标
        double sumX1 = 0, sumY1 = 0, sumX2 = 0, sumY2 = 0;
        for (int i = 0; i < pointCount; i++) {
            sumX1 += source[2 * i];
            sumY1 += source[2 * i + 1];
            sumX2 += target[2 * i];
            sumY2 += target[2 * i + 1];
        }

        double meanX1 = sumX1 / pointCount;
        double meanY1 = sumY1 / pointCount;
        double meanX2 = sumX2 / pointCount;
        double meanY2 = sumY2 / pointCount;

        // 中心化坐标
        double[] centeredSource = new double[source.length];
        double[] centeredTarget = new double[target.length];

        for (int i = 0; i < pointCount; i++) {
            centeredSource[2 * i] = source[2 * i] - meanX1;
            centeredSource[2 * i + 1] = source[2 * i + 1] - meanY1;
            centeredTarget[2 * i] = target[2 * i] - meanX2;
            centeredTarget[2 * i + 1] = target[2 * i + 1] - meanY2;
        }

        // 构建矩阵H和向量B
        double H11 = 0, H12 = 0, H21 = 0, H22 = 0;
        double B1 = 0, B2 = 0;

        for (int i = 0; i < pointCount; i++) {
            double x1 = centeredSource[2 * i];
            double y1 = centeredSource[2 * i + 1];
            double x2 = centeredTarget[2 * i];
            double y2 = centeredTarget[2 * i + 1];

            H11 += x1 * x1 + y1 * y1;
            H12 += 0;
            H21 += 0;
            H22 += x1 * x1 + y1 * y1;

            B1 += x1 * x2 + y1 * y2;
            B2 += x1 * y2 - y1 * x2;
        }

        // 求解参数
        double det = H11 * H22 - H12 * H21;
        if (Math.abs(det) < 1e-15) {
            System.out.println("矩阵奇异，无法求解参数");
            return new double[] { 0, 0, 0, 1 };
        }

        double a = (H22 * B1 - H12 * B2) / det;
        double b = (-H21 * B1 + H11 * B2) / det;

        double scale = Math.sqrt(a * a + b * b);
        double rotation = Math.atan2(b, a);

        double deltaX = meanX2 - (a * meanX1 - b * meanY1);
        double deltaY = meanY2 - (b * meanX1 + a * meanY1);

        System.out.println("四参数计算:");
        System.out.println("deltaX: " + deltaX);
        System.out.println("deltaY: " + deltaY);
        System.out.println("rotation (radians): " + rotation);
        System.out.println("scale: " + scale);

        return new double[] { deltaX, deltaY, rotation, scale };
    }

    public static double[] fourParameterTransform(double x, double y, double deltaX, double deltaY, double rotation,
            double scale) {
        // 计算旋转和缩放后的坐标
        double convertedX = scale * (x * Math.cos(rotation) - y * Math.sin(rotation)) + deltaX;
        double convertedY = scale * (x * Math.sin(rotation) + y * Math.cos(rotation)) + deltaY;
        return new double[] { convertedX, convertedY };
    }

    // 交点转线元
    public static double[][] jd2pqx(double[][] data) {
        List<Double> listXY1 = new ArrayList<>();
        double[][] lr18 = data;
        double sk = lr18[0][2];
        double JD2Mileage = 0.0;
        double hzx = 0, hzy = 0;

        // JD点坐标结构
        double JD1_X, JD1_Y;
        double JD2_X, JD2_Y;
        double JD3_X, JD3_Y;

        for (int i = 1; i < lr18.length - 1; i++) {
            // 获取前一个JD点坐标
            if (i == 1) {
                JD1_X = lr18[i - 1][0];
                JD1_Y = lr18[i - 1][1];
            } else {
                JD1_X = hzx;
                JD1_Y = hzy;
            }

            JD2_X = lr18[i][0];
            JD2_Y = lr18[i][1];
            JD3_X = lr18[i + 1][0];
            JD3_Y = lr18[i + 1][1];

            double R = lr18[i][4];
            double Ls1 = lr18[i][2];
            double Ls2 = lr18[i][3];
            int xyzy = 1;

            // 计算方位角
            double[] fwj12 = fwj(JD1_X, JD1_Y, JD2_X, JD2_Y);
            double[] fwj23 = fwj(JD2_X, JD2_Y, JD3_X, JD3_Y);
            double azimuth12 = fwj12[1];
            double azimuth23 = fwj23[1];

            // 计算转向角
            double alpha = azimuth23 - azimuth12;
            if (alpha < 0) {
                alpha = -alpha;
            }
            if (alpha > Math.PI) {
                alpha = Math.PI * 2 - alpha;
            }

            // 判断转向方向（左转-1，右转1）
            double area = calculateTriangleArea(JD1_X, JD1_Y, JD2_X, JD2_Y, JD3_X, JD3_Y);
            if (area < 0) {
                xyzy = -1;
            }

            // 计算曲线要素（需要实现calculateCurveElements方法）
            double[] curveElements = calculateCurveElements(Ls1, Ls2, R, alpha);
            double T1 = curveElements[0];
            double T2 = curveElements[1];
            double Ly = curveElements[2];
            double L = curveElements[3];

            // 计算主点桩号
            double dist1 = fwj12[0];
            JD2Mileage = sk + dist1;
            double ZH = JD2Mileage - T1;
            double HY = ZH + Ls1;
            // double QZ = ZH + Ls1 + (L - Ls2 - Ls1) / 2;
            double YH = ZH + L - Ls2;
            // double HZ = ZH + L;

            // 前直线段
            if (dist1 > T1 && dist1 - T1 > 0.01) {
                double xycd = dist1 - T1;
                Collections.addAll(listXY1, sk, JD1_X, JD1_Y, azimuth12, xycd, 0.0, 0.0, 0.0);
            }

            // 计算ZH点坐标
            double zhx = JD2_X - T1 * Math.cos(azimuth12);
            double zhy = JD2_Y - T1 * Math.sin(azimuth12);
            double zhk = sk + dist1 - T1;

            // 第一缓和曲线
            if (Ls1 != 0) {
                Collections.addAll(listXY1, zhk, zhx, zhy, azimuth12, Ls1, 0.0, R, (double) xyzy);
            }

            // 圆曲线段
            double[] xy = new double[] { 0, 0, 0 };
            if (Ly != 0 && Ls1 != 0) {
                xy = zs(ZH, zhx, zhy, azimuth12, Ls1, 0, R, xyzy, ZH + Ls1, 0, 90);
                Collections.addAll(listXY1, HY, xy[0], xy[1], xy[2], Ly, R, R, (double) xyzy);
            } else {
                Collections.addAll(listXY1, HY, zhx, zhy, azimuth12, Ly, R, R, (double) xyzy);
            }

            // 第二缓和曲线
            if (Ls2 != 0) {
                xy = zs(HY, xy[0], xy[1], xy[2], Ly, R, R, xyzy, HY + Ly, 0, 90);
                Collections.addAll(listXY1, YH, xy[0], xy[1], xy[2], Ls2, R, 0.0, (double) xyzy);
            }

            // 计算下一个交点坐标
            hzx = JD2_X + T2 * Math.cos(azimuth23);
            hzy = JD2_Y + T2 * Math.sin(azimuth23);
            sk = zhk + L;

            // 终点直线段（最后一个交点）
            if (i == lr18.length - 2) {
                double dist2 = fwj23[0];
                if (dist2 - T2 > 0.01) {
                    Collections.addAll(listXY1, sk, hzx, hzy, azimuth23, dist2 - T2, 0.0, 0.0, 0.0);
                }
            }
        }

        // 将List转换为二维数组
        int rows = listXY1.size() / 8;
        int cols = 8;
        double[][] pqx = new double[rows][cols];

        for (int i = 0; i < listXY1.size(); i++) {
            int row = i / cols;
            int col = i % cols;
            pqx[row][col] = listXY1.get(i);
        }

        // 数据后处理（四舍五入和角度转换）
        for (int i = 0; i < pqx.length; i++) {
            for (int j = 0; j < pqx[i].length; j++) {
                if (j != 3) {
                    pqx[i][j] = Math.round(pqx[i][j] * 1000.0) / 1000.0;
                } else {
                    pqx[i][j] = radiansToDMS(pqx[i][j]); // 需要实现此方法
                }
            }
        }

        return pqx;
    }

    private static double calculateTriangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {

        return 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    }

    private static double[] calculateCurveElements(double Ls1, double Ls2, double R, double alpha) {
        double T1 = 0;
        double T2 = 0;
        // double E;
        double L;
        // double J;
        // 计算缓和曲线参数 m1, m2, p1, p2
        double m1 = Ls1 / 2 - Math.pow(Ls1, 3) / (240 * Math.pow(R, 2)) - Math.pow(Ls1, 5) / (34560 * Math.pow(R, 4));
        double m2 = Ls2 / 2 - Math.pow(Ls2, 3) / (240 * Math.pow(R, 2)) - Math.pow(Ls2, 5) / (34560 * Math.pow(R, 4));
        double p1 = Math.pow(Ls1, 2) / (24 * R) - Math.pow(Ls1, 4) / (2688 * Math.pow(R, 3));
        double p2 = Math.pow(Ls2, 2) / (24 * R) - Math.pow(Ls2, 4) / (2688 * Math.pow(R, 3));
        // 计算切线长 T1 和 T2
        T1 = m1 + (R + p2 - (R + p1) * Math.cos(alpha)) / Math.sin(alpha);
        T2 = m2 + (R + p1 - (R + p2) * Math.cos(alpha)) / Math.sin(alpha);
        // 计算外矢距 E
        // E = (R + (p1 + p2) / 2) / Math.cos(alpha / 2) - R;
        // 计算缓和曲线角
        double beta01 = Ls1 / (2 * R);
        double beta02 = Ls2 / (2 * R);
        // 计算圆曲线长 Ly 和总曲线长 L
        double Ly = R * (alpha - beta01 - beta02);
        L = Ls1 + Ls2 + Ly;
        // 计算切曲差 J
        // J = T1 + T2 - L;
        // 返回所有计算结果的数组
        return new double[] { T1, T2, Ly, L, Ly, p1, p2, beta01, beta02 };
    }

    public static List<Point> offsetPolyline(double[] xCoords, double[] yCoords, double d) {
        int n = xCoords.length;
        if (n < 2 || yCoords.length != n)
            throw new IllegalArgumentException("输入坐标数组长度必须相等且至少包含2个点");

        List<Point> result = new ArrayList<>();

        for (int i = 0; i < n - 1; i++) {
            Point p1 = new Point(xCoords[i], yCoords[i]);
            Point p2 = new Point(xCoords[i + 1], yCoords[i + 1]);

            double dx = p2.x - p1.x;
            double dy = p2.y - p1.y;
            double len = Math.sqrt(dx * dx + dy * dy);

            if (len < 1e-10) {
                if (i == 0)
                    result.add(p1);
                if (i == n - 2)
                    result.add(p2);
                continue;
            }

            double ux = dx / len;
            double uy = dy / len;
            double nx = -uy;
            double ny = ux;

            Point offsetStart = new Point(p1.x + nx * d, p1.y + ny * d);
            Point offsetEnd = new Point(p2.x + nx * d, p2.y + ny * d);

            if (i == 0) {
                result.add(offsetStart);
            }

            if (i < n - 2) {
                Point p3 = new Point(xCoords[i + 2], yCoords[i + 2]);
                double dx2 = p3.x - p2.x;
                double dy2 = p3.y - p2.y;
                double len2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);

                if (len2 < 1e-10) {
                    result.add(offsetEnd);
                    continue;
                }
                double ux2 = dx2 / len2;
                double uy2 = dy2 / len2;
                double nx2 = -uy2;
                double ny2 = ux2;
                Point offset2Start = new Point(p2.x + nx2 * d, p2.y + ny2 * d);
                Point offset2End = new Point(p3.x + nx2 * d, p3.y + ny2 * d);
                double x1 = offsetStart.x, y1 = offsetStart.y;
                double x2 = offsetEnd.x, y2 = offsetEnd.y;
                double x3 = offset2Start.x, y3 = offset2Start.y;
                double x4 = offset2End.x, y4 = offset2End.y;
                double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
                if (Math.abs(denom) > 1e-10) {
                    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
                    double x = x1 + t * (x2 - x1);
                    double y = y1 + t * (y2 - y1);
                    result.add(new Point(x, y));
                } else {
                    result.add(offsetEnd);
                }
            } else {
                result.add(offsetEnd);
            }
        }
        return result;
    }

    // 线性插值
    public static double fromXgetY(double[] x, double[] y, double targetX) {
        double tolerance = 1e-6;

        // 使用二分查找确定 targetX 在数组 x 中的插入点
        int index = Arrays.binarySearch(x, targetX);

        // 如果恰好找到目标值，直接返回对应的 y 值
        if (index >= 0) {
            return y[index];
        }

        // 计算目标值应插入的位置（取补码）
        int rightIndex = ~index; // 等价于 C# 中的 ~index
        int leftIndex = rightIndex - 1;

        // 处理边界情况：如果目标值超出数组范围
        if (rightIndex >= x.length) {
            return y[y.length - 1]; // 返回最后一个 y 值
        }
        if (leftIndex < 0) {
            return y[0]; // 返回第一个 y 值
        }

        double xLeft = x[leftIndex];
        double xRight = x[rightIndex];

        // 如果左右边界 x 值过于接近，避免除以零，返回较大的 y 值
        if (Math.abs(xRight - xLeft) < tolerance) {
            return Math.max(y[leftIndex], y[rightIndex]);
        }

        double yLeft = y[leftIndex];
        double yRight = y[rightIndex];

        // 线性插值公式: y = yLeft + (yRight - yLeft) * (targetX - xLeft) / (xRight - xLeft)
        return yLeft + (yRight - yLeft) * (targetX - xLeft) / (xRight - xLeft);
    }

    public static double[] cutAndFillArea(double[] xsA, double[] ysA, double[] xsB, double[] ysB, double epsilon,
            double jgx) {
        // 扩展线段端点
        double[] left = extendLine(xsA[0], ysA[0], xsA[1], ysA[1], -epsilon);
        xsA[0] = left[0];
        ysA[0] = left[1];
        left = extendLine(xsB[0], ysB[0], xsB[1], ysB[1], -epsilon);
        xsB[0] = left[0];
        ysB[0] = left[1];
        left = extendLine(xsA[xsA.length - 1], ysA[ysA.length - 1], xsA[xsA.length - 2], ysA[ysA.length - 2], -epsilon);
        xsA[xsA.length - 1] = left[0];
        ysA[ysA.length - 1] = left[1];
        left = extendLine(xsB[xsB.length - 1], ysB[ysB.length - 1], xsB[xsB.length - 2], ysB[ysB.length - 2], -epsilon);
        xsB[xsB.length - 1] = left[0];
        ysB[ysB.length - 1] = left[1];

        // 存储交点信息
        List<IntersectionPoint> xys = new ArrayList<>();
        double fill = 0;
        double cut = 0;

        // 查找所有交点
        for (int i = 0; i < xsA.length - 1; i++) {
            double x1 = xsA[i], y1 = ysA[i];
            double x2 = xsA[i + 1], y2 = ysA[i + 1];

            for (int j = 0; j < xsB.length - 1; j++) {
                double x3 = xsB[j], y3 = ysB[j];
                double x4 = xsB[j + 1], y4 = ysB[j + 1];

                double[] jd = intersectionsSegment(x1, y1, x2, y2, x3, y3, x4, y4);
                if (jd.length > 0) {
                    xys.add(new IntersectionPoint(jd[0], jd[1], i, j));
                }
            }
        }

        // 处理交点序列，计算填挖面积
        for (int i = 0; i < xys.size() - 1; i++) {
            IntersectionPoint xy0 = xys.get(i);
            IntersectionPoint xy1 = xys.get(i + 1);
            int i0 = xy0.sji, i1 = xy1.sji;
            int j0 = xy0.dmj, j1 = xy1.dmj;

            List<Point> pts = new ArrayList<>();
            pts.add(new Point(xy0.x, xy0.y));

            // 添加设计线（sj）上的点
            for (int j = i0 + 1; j <= i1; j++) {
                pts.add(new Point(xsA[j], ysA[j]));
            }
            pts.add(new Point(xy1.x, xy1.y));

            // 添加地面线（dm）上的点（逆序）
            for (int j = j1; j > j0; j--) {
                pts.add(new Point(xsB[j], ysB[j]));
            }

            // 闭合多边形
            pts.add(new Point(xy0.x, xy0.y));

            // 计算多边形面积（鞋带公式）
            double area = 0;
            for (int k = 0; k < pts.size() - 1; k++) {
                Point p1 = pts.get(k);
                Point p2 = pts.get(k + 1);
                area += p1.x * p2.y - p1.y * p2.x;
            }
            area = Math.abs(area) / 2.0;

            // 根据面积正负判断是填方还是挖方
            double signedArea = 0;
            for (int k = 0; k < pts.size() - 1; k++) {
                Point p1 = pts.get(k);
                Point p2 = pts.get(k + 1);
                signedArea += p1.x * p2.y - p1.y * p2.x;
            }

            if (signedArea < 0) {
                fill += area; // 负面积为填方
            } else {
                cut += area; // 正面积为挖方
            }
        }

        return new double[] { Math.round(fill * 10000) / 10000.0, Math.round(cut * 10000) / 10000.0 };
    }

    public static double[] extendLine(double x1, double y1, double x2, double y2, double distance) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        double length = Math.sqrt(dx * dx + dy * dy);

        if (length == 0) {
            return new double[] { x1, y1 };
        }

        double ux = dx / length;
        double uy = dy / length;
        double newX, newY;

        if (distance < 0) {
            // 从起点反向延长
            newX = x1 + distance * ux;
            newY = y1 + distance * uy;
        } else {
            // 从终点正向延长
            newX = x2 + distance * ux;
            newY = y2 + distance * uy;
        }

        return new double[] { newX, newY };
    }

    public static double[] intersectionsSegment(double x1, double y1, double x2, double y2,
            double x3, double y3, double x4, double y4) {
        final double epsilon = 1e-10; // 添加精度容差

        // 计算向量叉积判断相对位置
        double d1 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
        double d2 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);
        double d3 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
        double d4 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);

        // 使用精度容差判断是否严格相交
        boolean isStrictlyIntersecting = (d1 * d2 <= epsilon) && (d3 * d4 <= epsilon);

        if (!isStrictlyIntersecting) {
            return new double[0]; // 无交点返回空数组
        }

        // 计算直线交点（线段无限延长）
        double a1 = y2 - y1;
        double b1 = x1 - x2;
        double c1 = x2 * y1 - x1 * y2;
        double a2 = y4 - y3;
        double b2 = x3 - x4;
        double c2 = x4 * y3 - x3 * y4;

        double denominator = a1 * b2 - a2 * b1;

        // 检查分母是否为零（平行或共线）
        if (Math.abs(denominator) < epsilon) {
            return new double[0]; // 平行或共线，无唯一交点
        }

        double px = (b1 * c2 - b2 * c1) / denominator;
        double py = (a2 * c1 - a1 * c2) / denominator;

        // 验证交点是否在两条线段范围内
        if (isPointOnSegment(px, py, x1, y1, x2, y2, epsilon) &&
                isPointOnSegment(px, py, x3, y3, x4, y4, epsilon)) {
            return new double[] { px, py };
        }

        return new double[0];
    }

    // 辅助方法：判断点是否在线段上
    private static boolean isPointOnSegment(double px, double py,
            double x1, double y1, double x2, double y2,
            double epsilon) {
        // 检查点是否在线段的边界框内
        boolean inXRange = (px >= Math.min(x1, x2) - epsilon) &&
                (px <= Math.max(x1, x2) + epsilon);
        boolean inYRange = (py >= Math.min(y1, y2) - epsilon) &&
                (py <= Math.max(y1, y2) + epsilon);

        if (!inXRange || !inYRange) {
            return false;
        }

        // 检查点是否共线（叉积接近零）
        double crossProduct = (px - x1) * (y2 - y1) - (py - y1) * (x2 - x1);
        return Math.abs(crossProduct) < epsilon;
    }

    public static double slope(double k, double[][] slope) {
        for (int i = 0; i < slope.length - 1; i++) {
            if (k >= slope[i][0] && k <= slope[i + 1][0]) {
                return slope[i][1]
                        + (slope[i + 1][1] - slope[i][1]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
            }
        }
        return 0;

    }

}

class Point {
    public double x;
    public double y;

    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    @Override
    public String toString() {
        return String.format("(%.3f, %.3f)", x, y);
    }
}

class IntersectionPoint {
    double x, y;
    int sji, dmj;

    IntersectionPoint(double x, double y, int sji, int dmj) {
        this.x = x;
        this.y = y;
        this.sji = sji;
        this.dmj = dmj;
    }
}