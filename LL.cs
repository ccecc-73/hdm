using System.Text;

namespace hengduanmian;
public class LL
{
    //1
    #region   数字转桩号格式 
    /// <summary>
    /// 将数字转换为 K 格式，例如：12345 -> K12+345
    /// </summary>
    /// <param name="meters">输入的数字</param>
    /// <returns>格式化后的字符串</returns>
    /// <example>LL.Num2K(12345)</example>
    public string hua_Num2K(double meters)
    {
        int km = (int)Math.Floor(meters / 1000);
        double m = meters - km * 1000;
        m = Math.Round(m, 2);

        if (m == (int)m)
        {
            // 整数米，如 2.00 → k1+002
            return $"K{km}+{(int)m:D3}";
        }
        else
        {
            // 小数米，如 2.35 → k0+002.35
            return $"K{km}+{m:000.00}";
        }
    }
    #endregion

    //2
    #region 格式为 DD.MMSS转弧度
    /// <summary>
    /// 将度分秒（DMS）格式转换为弧度
    /// 假设输入的 dms 是一个浮点数，格式为 DD.MMSS
    /// </summary>
    /// <param name="dms">度分秒格式的数值</param>
    /// <returns>对应的弧度值</returns>
    public static double hua_DmsToRadians(double dms)
    {
        int degrees = (int)dms;
        double fractional = dms - degrees;
        int minutes = (int)(fractional * 100);
        double seconds = (fractional * 100 - minutes) * 100;

        double totalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;
        return totalDegrees * Math.PI / 180;
    }
    #endregion
 

    public static double hua_radiansToDMS(double radians)
    {
        double degrees = radians * (180 / Math.PI);
        int d = (int)Math.Floor(degrees);
        double remaining = (degrees - d) * 60;
        int m = (int)Math.Floor(remaining);
        double s = (remaining - m) * 60;
        string mm = m < 10 ? $"0{m}" : m.ToString();
        string ss = s < 10 ? $"0{s:F1}" : $"{s:F1}";
        string dfm = $"{d}{mm}{ss}";

        return Math.Round(double.Parse(dfm) / 10000.0, 5);
    }

    public static string hua_radiansToDMS_度分秒(double radians)
    {
        double degrees = radians * (180 / Math.PI);
        int d = (int)Math.Floor(degrees);
        double remaining = (degrees - d) * 60;
        int m = (int)Math.Floor(remaining);
        double s = (remaining - m) * 60;
        string mm = m < 10 ? $"0{m}" : m.ToString();
        string ss = s < 10 ? $"0{s:F1}" : $"{s:F1}";
        return $"{d}°{mm}′{ss}″";
    }
    //3
    #region  距离和方位角
    /// <summary>
    /// 计算两点之间的极坐标距离和角度（弧度）
    /// </summary>
    /// <param name="x0">起点 X 坐标</param>
    /// <param name="y0">起点 Y 坐标</param>
    /// <param name="x1">终点 X 坐标</param>
    /// <param name="y1">终点 Y 坐标</param>
    /// <returns>包含两个元素的数组：[距离(cd), 角度(hd]) ]</returns>
    public static double[] hua_Fwj(double x0, double y0, double x1, double y1)
    {
        double x = x1 - x0;
        double y = y1 - y0;
        double cd = Math.Sqrt(x * x + y * y); // 距离
        double hd = Math.Atan2(y, x);         // 角度（弧度）

        if (hd < 0)
        {
            hd += 2 * Math.PI; // 确保角度在 [0, 2π) 范围内
        }
        return [cd, hd];
    }
    #endregion



    [DllImport("LL.dll", CallingConvention = CallingConvention.Cdecl)]
    static extern bool zs(double xyk, double xyx, double xyy, double xyhd, double xycd,
                               double xyqdr, double xyzdr, double xyzy, double jsk,
                               double jsb, double jd, double[] result); // 注意这里使用 double[]

    //4
    #region 由桩号计算xy
    /// <summary>
    /// 
    /// </summary>
    /// <param name="xyk">线元起点桩号</param>
    /// <param name="xyx">线元起点 北坐标</param>
    /// <param name="xyy">线元起点 东坐标</param>
    /// <param name="xyhd">线元起点 方位角(弧度)</param>
    /// <param name="xycd">线元长度</param>
    /// <param name="xyqdr">线元起点半径 直线为0 其他正常输入</param>
    /// <param name="xyzdr">规则相同起点</param>
    /// <param name="xyzy">线元方向 左转 -1  右转 1  直线 0</param>
    /// <param name="jsk">计算的里程</param>
    /// <param name="jsb">计算的宽度  左负右正 中桩=0</param>
    /// <param name="jd">右夹角  整数正常输入  若90°03分07.8秒 则输入 90.03078</param>
    /// <returns></returns>
    public static double[] hua_Zs(double xyk, double xyx, double xyy, double xyhd, double xycd, double xyqdr, double xyzdr, double xyzy, double jsk, double jsb, double jd)
    {
        double[] resultArray = new double[3];

        // 2. 调用 C 函数，并传入所有参数
        bool success = zs(xyk, xyx, xyy, xyhd, xycd, xyqdr, xyzdr, xyzy, jsk, jsb, jd, resultArray);

        // 3. 检查调用是否成功并处理结果
           return resultArray;
       
        /**
        jd = hua_DmsToRadians(jd);
        // 第一个条件分支  圆弧
        if (Math.Abs(xyqdr - xyzdr) < 0.01 && xyqdr > 0)
        {
            double centerX = xyx + xyqdr * Math.Cos(xyhd + xyzy * Math.PI / 2);
            double centerY = xyy + xyqdr * Math.Sin(xyhd + xyzy * Math.PI / 2);
            double deltaArcLength = jsk - xyk;
            double deltaAngle = deltaArcLength / xyqdr * xyzy;
            double targetAzimuth = xyhd + deltaAngle;
            if (targetAzimuth < 0)
                targetAzimuth += 2 * Math.PI;
            double angle = xyhd + xyzy * Math.PI / 2 + deltaAngle;
            double targetX = centerX - xyqdr * Math.Cos(angle) + jsb * Math.Cos(angle - xyzy * Math.PI / 2 + jd );
            double targetY = centerY - xyqdr * Math.Sin(angle) + jsb * Math.Sin(angle - xyzy * Math.PI / 2 + jd );

            return [targetX, targetY, targetAzimuth];
        }

        // 第二个条件分支 直线
        if (xyqdr < 0.01 && xyzdr < 0.01 && xyzy < 0.01)
        {
            double targetX = xyx + (jsk - xyk) * Math.Cos(xyhd) + jsb * Math.Cos(xyhd + jd );
            double targetY = xyy + (jsk - xyk) * Math.Sin(xyhd) + jsb * Math.Sin(xyhd + jd );

            return [targetX, targetY, xyhd];
        }

        // 第三个条件分支（默认情况） 缓和曲线
        if (xyqdr < 0.001) xyqdr = 99999999;
        if (xyzdr < 0.001) xyzdr = 99999999;

        double f0 = xyhd;
        double q = xyzy;
        double c = 1 / xyqdr;
        double d = (xyqdr - xyzdr) / (2 * xycd * xyqdr * xyzdr);

        // 初始化 rr 和 vv 数组，索引从1开始，因此大小为5（索引1-4）
        double[] rr = new double[5];
        double[] vv = new double[5];

        rr[1] = 0.1739274226;
        rr[2] = 0.3260725774;
        rr[3] = rr[2];
        rr[4] = rr[1];

        vv[1] = 0.0694318442;
        vv[2] = 0.3300094782;
        vv[3] = 1 - vv[2];
        vv[4] = 1 - vv[1];

        double w = jsk - xyk;
        double xs = 0;
        double ys = 0;

        for (int i = 1; i < 5; i++)
        {
            double ff = f0 + q * vv[i] * w * (c + vv[i] * w * d);
            xs += rr[i] * Math.Cos(ff);
            ys += rr[i] * Math.Sin(ff);
        }

        double fhz3 = f0 + q * w * (c + w * d);
        if (fhz3 < 0)
            fhz3 += 2 * Math.PI;
        if (fhz3 >= 2 * Math.PI)
            fhz3 -= 2 * Math.PI;

        double fhz1 = xyx + w * xs + jsb * Math.Cos(fhz3 + jd );
        double fhz2 = xyy + w * ys + jsb * Math.Sin(fhz3 + jd);
        return [fhz1, fhz2, fhz3];
        **/
    }
    #endregion
    
    //5
    #region 单条路线单个坐标计算
    /// <summary>
    /// 计算指定里程处的坐标和方位角
    /// </summary>
    /// <param name="pqx">路径点数组，每个元素是一个包含多个属性的数组</param>
    /// <param name="k">目标里程</param>
    /// <param name="b">宽度参数</param>
    /// <param name="z">角度参数</param>
    /// <returns>包含目标X坐标、Y坐标和方位角的数组</returns>
    public static double[] hua_Dantiaoxianludange(double[,] pqx, double k, double b, double z)
    {
        for (int i = 0; i < pqx.GetLength(0); i++)
        {
            double dtk = pqx[i, 0];
            double dtx = pqx[i, 1];
            double dty = pqx[i, 2];
            double dtfwj = pqx[i, 3];
            double dtcd = pqx[i, 4];
            double dtr1 = pqx[i, 5];
            double dtr2 = pqx[i, 6];
            double dtzy = pqx[i, 7];
            if (k >= dtk && k <= dtk + dtcd)
            {
                double hudu = hua_DmsToRadians(dtfwj);
                double[] jsxy1 = hua_Zs(dtk, dtx, dty, hudu, dtcd, dtr1, dtr2, dtzy, k, b, z);
                return [Math.Round(jsxy1[0], 3), Math.Round(jsxy1[1], 3), jsxy1[2]];
            }
        }
        return [0, 0, 0];
    }

    public static double[] hua_Dantiaoxianludange(double[][] pqx, double k, double b, double z)
    {
        for (int i = 0; i < pqx.GetLength(0); i++)
        {
            double dtk = pqx[i][ 0];
            double dtx = pqx[i][1];
            double dty = pqx[i][2];
            double dtfwj = pqx[i][3];
            double dtcd = pqx[i][4];
            double dtr1 = pqx[i][5];
            double dtr2 = pqx[i][6];
            double dtzy = pqx[i][7];
            if (k >= dtk && k <= dtk + dtcd)
            {
                double hudu = hua_DmsToRadians(dtfwj);
                double[] jsxy1 = hua_Zs(dtk, dtx, dty, hudu, dtcd, dtr1, dtr2, dtzy, k, b, z);
                return [Math.Round(jsxy1[0], 3), Math.Round(jsxy1[1], 3), jsxy1[2]];
            }
        }
        return [0, 0, 0];
    }
    #endregion

    //6
    #region xy=>k
    /// <summary>
    /// 
    /// </summary>
    /// <param name="pqx"></param>
    /// <param name="fsx"></param>
    /// <param name="fsy"></param>
    /// <returns></returns>
    public static double[] hua_Fs(double fsx, double fsy, double[,] pqx)
    {
        dynamic jljd = hua_Fwj(pqx[0, 1], pqx[0, 2], fsx, fsy);
        double k = pqx[0, 0];
        double hudu = hua_DmsToRadians(pqx[0, 3]);
        double cz = jljd[0] * Math.Cos(jljd[1] - hudu);
        double pj = jljd[0] * Math.Sin(jljd[1] - hudu);
        int hang = pqx.GetLength(0) - 1;
        double qdlc = pqx[0, 0];
        double zdlc = pqx[hang, 0] + pqx[hang, 4];
        int jisuancishu = 0;
        while (Math.Abs(cz) > 0.01)
        {
            k += cz;
            jisuancishu += 1;
            if (k < qdlc)
            {
                return [-1, -1];
            }
            if (k > zdlc)
            {
                return [-2, -2];
            }
            if (jisuancishu > 15)
            {
                return [-3, -3];
            }
            double[] xy = hua_Dantiaoxianludange(pqx, k, 0, 0);
            jljd = hua_Fwj((double)xy[0], (double)xy[1], fsx, fsy);
            cz = jljd[0] * Math.Cos(jljd[1] - xy[2]);
            pj = jljd[0] * Math.Sin(jljd[1] - xy[2]);
        }
        return [Math.Round(k, 3), Math.Round(pj, 3)];
    }
    #endregion
    //6
    #region xy=>k
    /// <summary>
    /// 
    /// </summary>
    /// <param name="pqx"></param>
    /// <param name="fsx"></param>
    /// <param name="fsy"></param>
    /// <returns></returns>
    public static double[] hua_Fs(double fsx, double fsy, double[][] pqx)
    {
        dynamic jljd = hua_Fwj(pqx[0][ 1], pqx[0][2], fsx, fsy);
        double k = pqx[0][0];
        double hudu = hua_DmsToRadians(pqx[0][3]);
        double cz = jljd[0] * Math.Cos(jljd[1] - hudu);
        double pj = jljd[0] * Math.Sin(jljd[1] - hudu);
        int hang = pqx.GetLength(0) - 1;
        double qdlc = pqx[0][0];
        double zdlc = pqx[hang][0] + pqx[hang][4];
        int jisuancishu = 0;
        while (Math.Abs(cz) > 0.01)
        {
            k += cz;
            jisuancishu += 1;
            if (k < qdlc)
            {
                return [-1, -1];
            }
            if (k > zdlc)
            {
                return [-2, -2];
            }
            if (jisuancishu > 15)
            {
                return [-3, -3];
            }
            double[] xy = hua_Dantiaoxianludange(pqx, k, 0, 0);
            jljd = hua_Fwj((double)xy[0], (double)xy[1], fsx, fsy);
            cz = jljd[0] * Math.Cos(jljd[1] - xy[2]);
            pj = jljd[0] * Math.Sin(jljd[1] - xy[2]);
        }
        return [Math.Round(k, 3), Math.Round(pj, 3)];
    }
    #endregion

    #region 
    private static double Gaocheng(double bpdlc, double bpdgc, double r, double qp, double hp, double t, double k)
    {
        double f = qp - hp;
        r = r * Math.Abs(f) / f;
        double x;
        if (k <= bpdlc - t)
        {
            x = 0;
        }
        else if (k >= bpdlc + t)
        {
            x = 0;
            qp = hp;
        }
        else
        {
            x = k - bpdlc + t;
        }

        return (bpdgc - (bpdlc - k) * qp - Math.Pow(x, 2) / 2 / r);
    }
    #endregion

    //7高程计算
    #region h
    /// <summary>
    /// 
    /// </summary>
    /// <param name="k"></param>
    /// <param name="sqxb"></param>
    /// <returns></returns>
    public static double hua_H(double k, double[,] sqxb)
    {
        double hp = 0;
        int length = sqxb.GetLength(0);
        for (int i = 1; i < length - 1; i++)
        {
            double r = sqxb[i, 2];
            if (r < 0.001)
                r = 0.001;
            double qp = (sqxb[i, 1] - sqxb[i - 1, 1]) / (sqxb[i, 0] - sqxb[i - 1, 0]);
            hp = (sqxb[i + 1, 1] - sqxb[i, 1]) / (sqxb[i + 1, 0] - sqxb[i, 0]);
            double f = qp - hp;
            double t = r * Math.Abs(f) / 2;
            if (k <= sqxb[i, 0] + t)
                return Math.Round(Gaocheng(sqxb[i, 0], sqxb[i, 1], r, qp, hp, t, k), 3);
        }
        //the last
        if (k <= sqxb[length - 1, 0])
        {
            return Math.Round(sqxb[length - 1, 1] + (k - sqxb[length - 1, 0]) * hp, 3);
        }
        return -1;
    }
    #endregion
    //7高程计算
    #region h
    /// <summary>
    /// 
    /// </summary>
    /// <param name="k"></param>
    /// <param name="sqxb"></param>
    /// <returns></returns>
    public static double H(double k, double[][] sqxb)
    {
        double hp = 0;
        int length = sqxb.GetLength(0);
        for (int i = 1; i < length - 1; i++)
        {
            double r = sqxb[i][2];
            if (r < 0.001)
                r = 0.001;
            double qp = (sqxb[i][1] - sqxb[i - 1][1]) / (sqxb[i][0] - sqxb[i - 1][0]);
            hp = (sqxb[i + 1][1] - sqxb[i][1]) / (sqxb[i + 1][0] - sqxb[i][0]);
            double f = qp - hp;
            double t = r * Math.Abs(f) / 2;
            if (k <= sqxb[i][0] + t)
                return Math.Round(Gaocheng(sqxb[i][0], sqxb[i][1], r, qp, hp, t, k), 3);
        }
        //the last
        if (k <= sqxb[length - 1][0])
        {
            return Math.Round(sqxb[length - 1][1] + (k - sqxb[length - 1][0]) * hp, 3);
        }
        return -1;
    }
    #endregion
    //高斯投影
    /// <summary>
    /// 
    /// </summary>
    /// <param name="L">经度</param>
    /// <param name="B">纬度</param>
    /// <param name="lonCenter">中心经度,若>360则自动计算</param>
    /// <returns>[北坐标,东坐标,中心经度]</returns>
    public static double[] hua_Gauss_proj(double L, double B, double lonCenter = 360.0)
    {
        double pi = 3.141592653589793238463;
        double p0 = 206264.8062470963551564;
        double e = 0.00669438002290;
        double e1 = 0.00673949677548;
        double b = 6356752.3141;
        double a = 6378137.0;
        B = B * pi / 180;
        L = L * pi / 180;
        double L_num;
        double L_center;
        if (lonCenter >=359)
        {
            L_num = Math.Floor(L * 180 / pi / 3.0 + 0.5);
            L_center = 3 * L_num;
        }
        else
        {
            L_center = lonCenter;
        }
        double l = (L / pi * 180 - L_center) * 3600;
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
        double Xz = a0 * B - a2 / 2.0 * Math.Sin(2 * B) + a4 / 4.0 * Math.Sin(4 * B) - a6 / 6.0 * Math.Sin(6 * B) + a8 / 8.0 * Math.Sin(8 * B);
        double c = a * a / b;
        double V = Math.Sqrt(1 + e1 * Math.Cos(B) * Math.Cos(B));
        double N = c / V;
        double t = Math.Tan(B);
        double n = e1 * Math.Cos(B) * Math.Cos(B);
        double m1 = N * Math.Cos(B);
        double m2 = N / 2.0 * Math.Sin(B) * Math.Cos(B);
        double m3 = N / 6.0 * Math.Pow(Math.Cos(B), 3) * (1 - t * t + n);
        double m4 = N / 24.0 * Math.Sin(B) * Math.Pow(Math.Cos(B), 3) * (5 - t * t + 9 * n);
        double m5 = N / 120.0 * Math.Pow(Math.Cos(B), 5) * (5 - 18 * t * t + Math.Pow(t, 4) + 14 * n - 58 * n * t * t);
        double m6 = N / 720.0 * Math.Sin(B) * Math.Pow(Math.Cos(B), 5) * (61 - 58 * t * t + Math.Pow(t, 4));
        double x = Xz + m2 * l * l / Math.Pow(p0, 2) + m4 * Math.Pow(l, 4) / Math.Pow(p0, 4) + m6 * Math.Pow(l, 6) / Math.Pow(p0, 6);
        double y0 = m1 * l / p0 + m3 * Math.Pow(l, 3) / Math.Pow(p0, 3) + m5 * Math.Pow(l, 5) / Math.Pow(p0, 5);
        double y = y0 + 500000;
        return [x, y, L_center];
    }
    //高斯 反投影
    /// <summary>
    /// 
    /// </summary>
    /// <param name="x">北坐标</param>
    /// <param name="y">东坐标</param>
    /// <param name="l0">中心经度</param>
    /// <returns>[经度,纬度]</returns>
    public static double[] hua_Gauss_unproj(double x, double y, double l0)
    {
        double pi = 3.141592653589793238463;
        double e = 0.00669438002290;
        double e1 = 0.00673949677548;
        double b = 6356752.3141;
        double a = 6378137.0;
        double y1 = y - 500000;
        double M0 = a * (1 - e);
        double M2 = 3.0 / 2.0 * e * M0;
        double M4 = 5.0 / 4.0 * e * M2;
        double M6 = 7.0 / 6.0 * e * M4;
        double M8 = 9.0 / 8.0 * e * M6;
        double a0 = M0 + M2 / 2.0 + 3.0 / 8.0 * M4 + 5.0 / 16.0 * M6 + 35.0 / 128.0 * M8;
        double a2 = M2 / 2.0 + M4 / 2 + 15.0 / 32.0 * M6 + 7.0 / 16.0 * M8;
        double a4 = M4 / 8.0 + 3.0 / 16.0 * M6 + 7.0 / 32.0 * M8;
        double a6 = M6 / 32.0 + M8 / 16.0;
        double Bf = x / a0;
        double B0 = Bf;
        while (Math.Abs(Bf - B0) > 0.0000001 || B0 == Bf)
        {
            B0 = Bf;
            double FBf = -a2 / 2.0 * Math.Sin(2 * B0) + a4 / 4.0 * Math.Sin(4 * B0) - a6 / 6.0 * Math.Sin(6 * B0);
            Bf = (x - FBf) / a0;
        }
        double t = Math.Tan(Bf);
        double c = a * a / b;
        double V = Math.Sqrt(1 + e1 * Math.Cos(Bf) * Math.Cos(Bf));
        double N = c / V;
        double M = c / Math.Pow(V, 3);
        double n = e1 * Math.Cos(Bf) * Math.Cos(Bf);
        double n1 = 1 / (N * Math.Cos(Bf));
        double n2 = -t / (2.0 * M * N);
        double n3 = -(1 + 2 * t * t + n) / (6.0 * Math.Pow(N, 3) * Math.Cos(Bf));
        double n4 = t * (5 + 3 * t * t + n - 9 * n * t * t) / (24.0 * M * Math.Pow(N, 3));
        double n5 = (5 + 28 * t * t + 24 * Math.Pow(t, 4) + 6 * n + 8 * n * t * t) / (120.0 * Math.Pow(N, 5) * Math.Cos(Bf));
        double n6 = -t * (61 + 90 * t * t + 45 * Math.Pow(t, 4)) / (720.0 * M * Math.Pow(N, 5));
        double B = (Bf + n2 * y1 * y1 + n4 * Math.Pow(y1, 4) + n6 * Math.Pow(y1, 6)) / pi * 180;
        double L0 = l0;
        double l = n1 * y1 + n3 * Math.Pow(y1, 3) + n5 * Math.Pow(y1, 5);
        double L = L0 + l / pi * 180;
        return [L, B];
    }

    // 经纬度转UTM坐标
/// <summary>
/// 
/// </summary>
/// <param name="longitude"></param>
/// <param name="latitude"></param>
/// <returns></returns>
    public static double[] hua_Utm_proj(double longitude, double latitude)
    {
        double EQUATORIAL_RADIUS = 6378137.0;
        double FLATTENING = 1 / 298.257223563;
        double ECC_SQUARED = 2 * FLATTENING - Math.Pow(FLATTENING, 2);
        double ECC_PRIME_SQUARED = ECC_SQUARED / (1 - ECC_SQUARED);
        double SCALE_FACTOR = 0.9996;
        double FALSE_EASTING = 500000.0;
        double FALSE_NORTHING_S = 10000000.0;
        double centralMeridian = 0;
        int zoneNumber = (int)Math.Floor((longitude + 180) / 6) + 1;
        centralMeridian = (zoneNumber - 1) * 6 - 180 + 3;
        double latRad = (latitude) * Math.PI / 180.0;
        double lonRad = (longitude) * Math.PI / 180.0;
        double lonCenterRad = (centralMeridian) * Math.PI / 180.0;
        double N = EQUATORIAL_RADIUS / Math.Sqrt(1 - ECC_SQUARED * Math.Pow(Math.Sin(latRad), 2));
        double T = Math.Pow(Math.Tan(latRad), 2);
        double C = ECC_PRIME_SQUARED * Math.Pow(Math.Cos(latRad), 2);
        double A = (lonRad - lonCenterRad) * Math.Cos(latRad);
        double M = EQUATORIAL_RADIUS * ((1 - ECC_SQUARED / 4 - 3 * Math.Pow(ECC_SQUARED, 2) / 64 - 5 * Math.Pow(ECC_SQUARED, 3) / 256) * latRad - (3 * ECC_SQUARED / 8 + 3 * Math.Pow(ECC_SQUARED, 2) / 32 + 45 * Math.Pow(ECC_SQUARED, 3) / 1024) * Math.Sin(2 * latRad) + (15 * Math.Pow(ECC_SQUARED, 2) / 256 + 45 * Math.Pow(ECC_SQUARED, 3) / 1024) * Math.Sin(4 * latRad) - (35 * Math.Pow(ECC_SQUARED, 3) / 3072) * Math.Sin(6 * latRad));
        double easting = SCALE_FACTOR * N * (A + (1 - T + C) * Math.Pow(A, 3) / 6 + (5 - 18 * T + Math.Pow(T, 2) + 72 * C - 58 * ECC_PRIME_SQUARED) * Math.Pow(A, 5) / 120) + FALSE_EASTING;
        double northing = SCALE_FACTOR * (M + N * Math.Tan(latRad) * (Math.Pow(A, 2) / 2 + (5 - T + 9 * C + 4 * Math.Pow(C, 2)) * Math.Pow(A, 4) / 24 + (61 - 58 * T + Math.Pow(T, 2) + 600 * C - 330 * ECC_PRIME_SQUARED) * Math.Pow(A, 6) / 720));
        if (latitude < 0) northing += FALSE_NORTHING_S;
        return [northing, easting, zoneNumber];
    }

    // UTM坐标转经纬度
    /// <summary>
    /// 
    /// </summary>
    /// <param name="northing"></param>
    /// <param name="easting"></param>
    /// <param name="isNorthern"></param>
    /// <param name="zoneNumber">(经度 + 180) / 6) + 1) 舍去小数部分</param>
    /// <returns></returns>
    public static double[] hua_Utm_unproj(double northing, double easting, bool isNorthern,int zoneNumber)
    {
        double EQUATORIAL_RADIUS = 6378137.0;
        double FLATTENING = 1 / 298.257223563;
        double ECC_SQUARED = 2 * FLATTENING - Math.Pow(FLATTENING, 2);
        double ECC_PRIME_SQUARED = ECC_SQUARED / (1 - ECC_SQUARED);
        double SCALE_FACTOR = 0.9996;
        double FALSE_EASTING = 500000.0;
        double FALSE_NORTHING_S = 10000000.0;
        double x = easting - FALSE_EASTING;
        double y = isNorthern ? northing : northing - FALSE_NORTHING_S;

        double centralMeridian = (zoneNumber - 1) * 6 - 180 + 3;
        
        double lonCenterRad = (centralMeridian) * Math.PI / 180.0;
        double M = y / SCALE_FACTOR;
        double mu = M / (EQUATORIAL_RADIUS * (1 - ECC_SQUARED / 4 - 3 * Math.Pow(ECC_SQUARED, 2) / 64.0 - 5 * Math.Pow(ECC_SQUARED, 3) / 256.0));
        double e1 = (1 - Math.Sqrt(1 - ECC_SQUARED)) / (1 + Math.Sqrt(1 - ECC_SQUARED));
        double phi1Rad = mu + (3 * e1 / 2 - 27 * Math.Pow(e1, 3) / 32) * Math.Sin(2 * mu)+ (21 * Math.Pow(e1, 2) / 16 - 55 * Math.Pow(e1, 4) / 32) * Math.Sin(4 * mu) + (151 * Math.Pow(e1, 3) / 96) * Math.Sin(6 * mu);
        double N1 = EQUATORIAL_RADIUS / Math.Sqrt(1 - ECC_SQUARED * Math.Pow(Math.Sin(phi1Rad), 2));
        double T1 = Math.Pow(Math.Tan(phi1Rad), 2);
        double C1 = ECC_PRIME_SQUARED * Math.Pow(Math.Cos(phi1Rad), 2);
        double R1 = EQUATORIAL_RADIUS * (1 - ECC_SQUARED) / Math.Pow(1 - ECC_SQUARED * Math.Pow(Math.Sin(phi1Rad), 2), 1.5);
        double D = x / (N1 * SCALE_FACTOR);
        double latRad = phi1Rad - (N1 * Math.Tan(phi1Rad) / R1)* (Math.Pow(D, 2) / 2 - (5 + 3 * T1 + 10 * C1 - 4 * Math.Pow(C1, 2) - 9 * ECC_PRIME_SQUARED)* Math.Pow(D, 4) / 24 + (61 + 90 * T1 + 298 * C1 + 45 * Math.Pow(T1, 2)- 252 * ECC_PRIME_SQUARED - 3 * Math.Pow(C1, 2)) * Math.Pow(D, 6) / 720);
        double lonRad = lonCenterRad + (D - (1 + 2 * T1 + C1) * Math.Pow(D, 3) / 6 + (5 - 2 * C1 + 28 * T1 - 3 * Math.Pow(C1, 2) + 8 * ECC_PRIME_SQUARED + 24 * Math.Pow(T1, 2)) * Math.Pow(D, 5) / 120) / Math.Cos(phi1Rad);
        return [ (lonRad) * 180 / Math.PI, (latRad) * 180 / Math.PI];
    }
    //utm投影带号
    public static int utm_zone(double longitude)
    {
        return (int)Math.Floor((longitude + 180) / 6) + 1;
    }
    //四参数
    // 应用转换公式
    //double convertedX = params.scale* (x* Math.cos(params.rotation) - y* Math.sin(params.rotation)) + params.deltaX;
    //        double convertedY = params.scale* (x* Math.sin(params.rotation) + y* Math.cos(params.rotation)) + params.deltaY;
    /// <summary>
    /// 计算两组坐标系下的四参数转换参数
    /// </summary>
    /// <param name="source">[x1,y1,x2,y2...]</param>
    /// <param name="target">[x1,y1,x2,y2...]</param>
    /// <returns>[deltaX, deltaY, rotation(弧度), scale]</returns>
    /// <exception cref="ArithmeticException"></exception>
    public static double[] hua_Cs4(double[] source, double[] target)
    {
        if (source == null || target == null || source.Length != target.Length)
        {
            Console.WriteLine("坐标数组长度必须相等");
            return [0, 0, 0, 1];
        }
        if (source.Length < 4 || source.Length % 2 != 0)
        {
            Console.WriteLine("至少需要2个点且坐标为偶数");
            return [0, 0, 0, 1];
        }
        int pointCount = source.Length / 2;
        double sumX1 = 0, sumY1 = 0, sumX2 = 0, sumY2 = 0;
        for (int i = 0; i < pointCount; i++)
        {
            sumX1 += source[2 * i];
            sumY1 += source[2 * i + 1];
            sumX2 += target[2 * i];
            sumY2 += target[2 * i + 1];
        }
        double meanX1 = sumX1 / pointCount;
        double meanY1 = sumY1 / pointCount;
        double meanX2 = sumX2 / pointCount;
        double meanY2 = sumY2 / pointCount;
        double[] centeredSource = new double[source.Length];
        double[] centeredTarget = new double[target.Length];
        for (int i = 0; i < pointCount; i++)
        {
            centeredSource[2 * i] = source[2 * i] - meanX1;
            centeredSource[2 * i + 1] = source[2 * i + 1] - meanY1;
            centeredTarget[2 * i] = target[2 * i] - meanX2;
            centeredTarget[2 * i + 1] = target[2 * i + 1] - meanY2;
        }
        double H11 = 0, H12 = 0, H21 = 0, H22 = 0;
        double B1 = 0, B2 = 0;
        for (int i = 0; i < pointCount; i++)
        {
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
        double det = H11 * H22 - H12 * H21;
        if (Math.Abs(det) < 1e-15)
        {
            Console.WriteLine("矩阵奇异，无法求解参数");
        }
        double a = (H22 * B1 - H12 * B2) / det;
        double b = (-H21 * B1 + H11 * B2) / det;
        double scale = Math.Sqrt(a * a + b * b);
        double rotation = Math.Atan2(b, a);
        double deltaX = meanX2 - (a * meanX1 - b * meanY1);
        double deltaY = meanY2 - (b * meanX1 + a * meanY1);
        Console.WriteLine($"四参数计算:\ndeltaX: {deltaX}\ndeltaY:{deltaY}\nrotation (radians): {rotation}\nscale: {scale}");
        return [deltaX, deltaY, rotation, scale];
    }
    /// <summary>
    /// 已知四参数进行坐标转换
    /// </summary>
    /// <param name="x">原坐标x</param>
    /// <param name="y">原坐标y</param>
    /// <param name="deltaX">平移x</param>
    /// <param name="deltaY">平移y</param>
    /// <param name="rotation">旋转弧度</param>
    /// <param name="scale">缩放</param>
    /// <returns>[x,y]</returns>
    public static double[] hua_FourParameterTransform(double x, double y, double deltaX, double deltaY, double rotation, double scale)
    {
        double convertedX = scale * (x * Math.Cos(rotation) - y * Math.Sin(rotation)) + deltaX;
        double convertedY = scale * (x * Math.Sin(rotation) + y * Math.Cos(rotation)) + deltaY;
        //Console.WriteLine($"转换后坐标:\nx: {convertedX}\ny: {convertedY}");
        return [convertedX, convertedY];
    }




    //以下为交点转线元
    /// <summary>
    /// 
    /// </summary>
    /// <param name="data"></param>
    /// <returns></returns>
    public static double[,] hua_JD2PQX(double[,] data)
    {
        List<double> listXY1 = [];
        var lr18 = data;
        var sk = lr18[0, 2];
        var JD2Mileage = 0.0;
        double hzx = 0, hzy = 0;
        (double X, double Y) JD1;
        (double X, double Y) JD2;
        (double X, double Y) JD3;
        for (int i = 1; i < lr18.GetLength(0) - 1; i++)
        {
            JD1 = (lr18[i - 1, 0], lr18[i - 1, 1]);

            if (i > 1)
            {
                JD1 = (hzx, hzy);
            }
            JD2 = (lr18[i, 0], lr18[i, 1]);

            JD3 = (lr18[i + 1, 0], lr18[i + 1, 1]);

            var R = lr18[i, 4];
            var Ls1 = lr18[i, 2];
            var Ls2 = lr18[i, 3];
            var xyzy = 1;
            var azimuth12 = calculateAzimuth(JD1, JD2);//方位角
            var azimuth23 = calculateAzimuth(JD2, JD3);
            var alpha = azimuth23 - azimuth12;
            if (alpha < 0)
            {
                alpha = -alpha;
            }
            if (alpha > Math.PI)
            {
                alpha = Math.PI * 2 - alpha;
            }
            var area = calculateTriangleArea(JD1.X, JD1.Y, JD2.X, JD2.Y, JD3.X, JD3.Y);
            if (area < 0)
                xyzy = -1;
            var (T1, T2, Ly, L, a2, a3, a4, a5, a6) = calculateCurveElements(Ls1, Ls2, R, alpha);//切线
            double ZH, HY, QZ, YH, HZ;
            double zhx, zhy;
            var dist1 = calculateDistance(JD1.X, JD1.Y, JD2.X, JD2.Y);
            JD2Mileage = sk + dist1;
            ZH = JD2Mileage - T1;
            HY = ZH + Ls1;
            QZ = ZH + Ls1 + (L - Ls2 - Ls1) / 2;
            YH = ZH + L - Ls2;
            HZ = ZH + L;

            //"jd" + i + "前直线"
            if (dist1 > T1 && dist1 - T1 > 0.01)
            {
                double xycd = dist1 - T1;

                listXY1.AddRange([sk, JD1.X, JD1.Y, azimuth12, xycd, 0, 0, 0]);
            }
            zhx = JD2.X - T1 * Math.Cos(azimuth12);
            zhy = JD2.Y - T1 * Math.Sin(azimuth12);
            var zhk = sk + dist1 - T1;
            // "jd" + i + "一缓", 
            if (Ls1 != 0)
            {
                listXY1.AddRange([zhk, zhx, zhy, azimuth12, Ls1, 0, R, xyzy]);
            }
            //var xy = {
            //    x: 0,
            //    y: 0,
            //    "fwj": 0,
            //    rad: 0
            //}
        ;
            //"jd" + i + "圆弧", 
            double[] xy = [0, 0, 0];
            if (Ly != 0 && Ls1 != 0)
            {
                xy = LL.hua_Zs(ZH, zhx, zhy, azimuth12, Ls1, 0, R, xyzy, ZH + Ls1, 0, 90);
                //[targetX, targetY, targetAzimuth];
                listXY1.AddRange([HY, xy[0], xy[1], xy[2], Ly, R, R, xyzy]);
            }
            else
            {
                //"圆弧jd" + i,
                //xy = {
                //x: zhx,
                //        y: zhy,
                //        "fwj": 0,
                //        rad: azimuth12
                //        }
                ;
                listXY1.AddRange([HY, zhx, zhy, azimuth12, Ly, R, R, xyzy]);
            }
            //"jd" + i + "二缓", 
            if (Ls2 != 0)
            {
                xy = hua_Zs(HY, xy[0], xy[1], xy[2], Ly, R, R, xyzy, HY + Ly, 0, 90);
                listXY1.AddRange([YH, xy[0], xy[1], xy[2], Ls2, R, 0, xyzy]);
            }
            hzx = JD2.X + T2 * Math.Cos(azimuth23);
            hzy = JD2.Y + T2 * Math.Sin(azimuth23);
            sk = zhk + L;
            //"jd" + i + "终点的直线", 
            if (i == lr18.GetLength(0) - 2)
            {
                var dist2 = calculateDistance(JD2.X, JD2.Y, JD3.X, JD3.Y);
                if (dist2 - T2 > 0.01)
                {
                    listXY1.AddRange([sk, hzx, hzy, azimuth23, dist2 - T2, 0, 0, 0]);
                }
            }
        }
        int rows = listXY1.Count / 8;
        int cols = 8;
        double[,] pqx = new double[rows, 8];
        for (int i = 0; i < listXY1.Count; i++)
        {
            int row = i / cols;  // 计算行索引： 
            int col = i % cols;  // 计算列索引： 
            pqx[row, col] = listXY1[i];
        }
        for (int i = 0; i < pqx.GetLength(0); i++)
        {
            for (int j = 0; j < pqx.GetLength(1); j++)
            {
                if (j != 3)
                    pqx[i, j] = Math.Round(pqx[i, j], 3);
                else
                    pqx[i, 3] = hua_radiansToDMS(pqx[i, 3]);
            }
        }
        return pqx;
    }



    private static double calculateAzimuth((double X, double Y) point1, (double X, double Y) point2)
    {
        var dx = point2.X - point1.X;
        var dy = point2.Y - point1.Y;
        var azimuth = Math.Atan2(dy, dx);
        if (azimuth < 0)
        {
            azimuth += 2 * Math.PI;
        }
        return azimuth;
    }

    private static double calculateTriangleArea(double x1, double y1, double x2, double y2, double x3, double y3)
    {

        return 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    }

    private static (double, double, double, double, double, double, double, double, double) calculateCurveElements(double Ls1, double Ls2, double R, double alpha)
    {
        double T1 = 0;
        double T2 = 0;
        double E, L, J;
        var m1 = Ls1 / 2 - Math.Pow(Ls1, 3) / (240 * Math.Pow(R, 2)) - Math.Pow(Ls1, 5) / (34560 * Math.Pow(R, 4));
        var m2 = Ls2 / 2 - Math.Pow(Ls2, 3) / (240 * Math.Pow(R, 2)) - Math.Pow(Ls2, 5) / (34560 * Math.Pow(R, 4));
        var p1 = Math.Pow(Ls1, 2) / (24 * R) - Math.Pow(Ls1, 4) / (2688 * R * R * R);
        var p2 = Math.Pow(Ls2, 2) / (24 * R) - Math.Pow(Ls2, 4) / (2688 * R * R * R);
        T1 = m1 + (R + p2 - (R + p1) * Math.Cos(alpha)) / Math.Sin(alpha);
        T2 = m2 + (R + p1 - (R + p2) * Math.Cos(alpha)) / Math.Sin(alpha);
        E = (R + (p1 + p2) / 2) / Math.Cos(alpha / 2) - R;
        var beta01 = Ls1 / (2 * R);
        var beta02 = Ls2 / (2 * R);
        double Ly = R * (alpha - beta01 - beta02);
        L = Ls1 + Ls2 + Ly;
        J = T1 + T2 - L;
        return (T1, T2, Ly, L, Ly, p1, p2, beta01, beta02);
    }
    private static double calculateDistance(double x1, double y1, double x2, double y2)
    {
        var dx = Math.Pow(x2 - x1, 2);
        var dy = Math.Pow(y2 - y1, 2);
        return Math.Sqrt(dx + dy);
    }
    //以上为交点转线元


    //polyline offset
    #region

    public static List<(double X, double Y)> hua_OffsetPolyline(double[] xCoords, double[] yCoords, double d)
    {
        int n = xCoords.Length;
        if (n < 2 || yCoords.Length != n)
            throw new ArgumentException("输入坐标数组长度必须相等且至少包含2个点");

        var result = new List<(double X, double Y)>();

        for (int i = 0; i < n - 1; i++)
        {
            (double X, double Y) p1 = (xCoords[i], yCoords[i]);
            (double X, double Y) p2 = (xCoords[i + 1], yCoords[i + 1]);

            double dx = p2.X - p1.X;
            double dy = p2.Y - p1.Y;
            double len = Math.Sqrt(dx * dx + dy * dy);

            if (len < 1e-10)
            {
                if (i == 0) result.Add(p1);
                if (i == n - 2) result.Add(p2);
                continue;
            }

            double ux = dx / len;
            double uy = dy / len;
            double nx = -uy;
            double ny = ux;

            (double X, double Y) offsetStart = (p1.X + nx * d, p1.Y + ny * d);
            (double X, double Y) offsetEnd = (p2.X + nx * d, p2.Y + ny * d);

            if (i == 0)
            {
                result.Add(offsetStart);
            }

            if (i < n - 2)
            {
                (double X, double Y) p3 = (xCoords[i + 2], yCoords[i + 2]);
                double dx2 = p3.X - p2.X;
                double dy2 = p3.Y - p2.Y;
                double len2 = Math.Sqrt(dx2 * dx2 + dy2 * dy2);

                if (len2 < 1e-10)
                {
                    result.Add(offsetEnd);
                    continue;
                }

                double ux2 = dx2 / len2;
                double uy2 = dy2 / len2;
                double nx2 = -uy2;
                double ny2 = ux2;

                (double X, double Y) offset2Start = (p2.X + nx2 * d, p2.Y + ny2 * d);
                (double X, double Y) offset2End = (p3.X + nx2 * d, p3.Y + ny2 * d);

                double x1 = offsetStart.X, y1 = offsetStart.Y;
                double x2 = offsetEnd.X, y2 = offsetEnd.Y;
                double x3 = offset2Start.X, y3 = offset2Start.Y;
                double x4 = offset2End.X, y4 = offset2End.Y;

                double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

                if (Math.Abs(denom) > 1e-10)
                {
                    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
                    double x = x1 + t * (x2 - x1);
                    double y = y1 + t * (y2 - y1);
                    result.Add((x, y));
                }
                else
                {
                    result.Add(offsetEnd);
                }
            }
            else
            {
                result.Add(offsetEnd);
            }
        }

        return result;
    }

    #endregion


    //fromXgetY
    #region
    /// <summary>
    /// 根据给定的x坐标数组和y坐标数组，计算指定x值对应的y值。
    /// 当x值相等时，返回较大的y值。
    /// </summary>
    /// <param name="x">已知点的x坐标数组（必须升序排列）</param>
    /// <param name="y">已知点的y坐标数组</param>
    /// <param name="targetX">要计算的目标x值</param>
    /// <param name="tolerance">判断x值相等的容差（默认1e-6）</param>
    /// <returns>目标x值对应的y值</returns>
    internal static double hua_fromXgetY(double[] x, double[] y, double targetX)
    {
        double tolerance = 1e-6;
        //  使用二分查找确定targetX所在区间[6](@ref)
        int index = Array.BinarySearch(x, targetX);
        // 如果恰好找到目标x值，直接返回对应的y值
        if (index >= 0)
        {
            return y[index];
        }
        int rightIndex = ~index;
        int leftIndex = rightIndex - 1;
        double xLeft = x[leftIndex];
        double xRight = x[rightIndex];
        if (Math.Abs(xRight - xLeft) < tolerance)
        {
            return Math.Max(y[leftIndex], y[rightIndex]);
        }
        double yLeft = y[leftIndex];
        double yRight = y[rightIndex];
        return yLeft + (yRight - yLeft) * (targetX - xLeft) / (xRight - xLeft);
    }
    #endregion


    #region
    [DllImport("cutfill.dll", CallingConvention = CallingConvention.Cdecl)]
    static extern void cutAndFillArea(double[] xsA, double[] ysA, int lengthA, double[] xsB, double[] ysB, int lengthB, double epsilon, double[] result);
    /// <summary>
    /// 
    /// </summary>
    /// <param name="xsA">一维数组</param>
    /// <param name="ysA"></param>
    /// <param name="xsB"></param>
    /// <param name="ysB"></param>
    /// <param name="epsilon">误差</param>
    /// <returns></returns>
    internal static double[] hua_cutAndFillArea(double[] xsA, double[] ysA, double[] xsB, double[] ysB, double epsilon)
    {
        double[] result = new double[2];
        if(xsA.Length != ysA.Length )
        {
            return [0.0, 0.0];
        }
        if(xsB.Length != ysB.Length)
        {
            return [0.0, 0.0];
        }

        cutAndFillArea(xsA, ysA, xsA.Length, xsB, ysB, xsB.Length, epsilon, result);
        return [Math.Round(result[0], 4), Math.Round(result[1], 4)];
       
    }
    #endregion


    //线段
    #region
    /// <summary>
    /// 两线段交点,如果没有则为[]空
    /// </summary>
    /// <param name="x1"></param>
    /// <param name="y1"></param>
    /// <param name="x2"></param>
    /// <param name="y2"></param>
    /// <param name="x3"></param>
    /// <param name="y3"></param>
    /// <param name="x4"></param>
    /// <param name="y4"></param>
    /// <returns></returns>
    [ExcelFunction(Description = "线段延长", Category = "横断面相关")]
    public static double[] hua_Intersections_Segment(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
    {
        double d1 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3); // 点(x1,y1)相对于线段(x3,y3)-(x4,y4)
        double d2 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3); // 点(x2,y2)相对于线段(x3,y3)-(x4,y4)
        double d3 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1); // 点(x3,y3)相对于线段(x1,y1)-(x2,y2)
        double d4 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1); // 点(x4,y4)相对于线段(x1,y1)-(x2,y2)

        bool isStrictlyIntersecting = (d1 * d2 <= 0) && (d3 * d4 <= 0);

        if (!isStrictlyIntersecting)
            return [];


        double a1 = y2 - y1;
        double b1 = x1 - x2;
        double c1 = x2 * y1 - x1 * y2;
        double a2 = y4 - y3;
        double b2 = x3 - x4;
        double c2 = x4 * y3 - x3 * y4;
        double denominator = a1 * b2 - a2 * b1;
        double px = (b1 * c2 - b2 * c1) / denominator;
        double py = (a2 * c1 - a1 * c2) / denominator;
        return [px, py];
    }
    #endregion


    //延长
    #region
    [ExcelFunction(Description = "线段延长", Category = "横断面相关")]
    public static double[] hua_ExtendLine(double x1, double y1, double x2, double y2, double distance)
    {
        double dx = x2 - x1;
        double dy = y2 - y1;
        double length = Math.Sqrt(dx * dx + dy * dy);
        if (length == 0)
        {
            return [x1, y1];
        }
        double ux = dx / length;
        double uy = dy / length;
        double newX; double newY;
        if (distance < 0)
        {
            newX = x1 + distance * ux;
            newY = y1 + distance * uy;
        }
        else
        {
            newX = x2 + distance * ux;
            newY = y2 + distance * uy;
        }
        return [newX, newY];
    }
    #endregion

    //slope
    #region
    public static double[] slope(double k, double[][] slope)
    {
        for (var i = 0; i < slope.Length - 1; i++)
        {
            if (k >= slope[i][0] && k <= slope[i + 1][0])
            {
                double L = slope[i][1]
                        + (slope[i + 1][1] - slope[i][1]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                double R = slope[i][2]
                        + (slope[i + 1][2] - slope[i][2]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                return [Math.Round(L, 3), Math.Round(R, 3)];
            }
        }
        return [0, 0];
    }
    #endregion
    #region
    public static double hua_slope(double k, double[,] slope)
    {
        for (var i = 0; i < slope.GetLength(0) - 1; i++)
        {
            if (k >= slope[i,0] && k <= slope[i + 1,0])
            {
                return slope[i,1]+ (slope[i + 1,1] - slope[i,1]) / (slope[i + 1,0] - slope[i,0]) * (k - slope[i,0]);                        
            }
        }
        return 0;
    }
    #endregion

    //根据给定的x坐标数组和y坐标数组，计算指定x值对应的y值。
    #region
    /// <summary>
    /// 根据给定的x坐标数组和y坐标数组，计算指定x值对应的y值。
    /// 当x值相等时，返回较大的y值。
    /// </summary>
    /// <param name="x">已知点的x坐标数组（必须升序排列）</param>
    /// <param name="y">已知点的y坐标数组</param>
    /// <param name="targetX">要计算的目标x值</param>
    /// <param name="tolerance">判断x值相等的容差（默认1e-6）</param>
    /// <returns>目标x值对应的y值</returns>
    internal static double Interpolate(double[] x, double[] y, double targetX)
    {
        double tolerance = 1e-6;
        //  使用二分查找确定targetX所在区间[6](@ref)
        int index = Array.BinarySearch(x, targetX);
        // 如果恰好找到目标x值，直接返回对应的y值
        if (index >= 0)
        {
            return y[index];
        }
        int rightIndex = ~index;
        int leftIndex = rightIndex - 1;
        double xLeft = x[leftIndex];
        double xRight = x[rightIndex];
        if (Math.Abs(xRight - xLeft) < tolerance)
        {
            return Math.Max(y[leftIndex], y[rightIndex]);
        }
        double yLeft = y[leftIndex];
        double yRight = y[rightIndex];
        return yLeft + (yRight - yLeft) * (targetX - xLeft) / (xRight - xLeft);
    }
    #endregion


    // 折线交点
    #region
    [ExcelFunction(Description = "折线交点", Category = "横断面相关")]
    public static double[,] hua_Intersections_polyline(double[] xsA, double[] ysA, double[] xsB, double[] ysB, double epsilon)
    {
        #region
        //extend sj
        if (Math.Abs(epsilon) > 0.001)
        {
            double[] left = hua_ExtendLine(xsA[0], ysA[0], xsA[1], ysA[1], -epsilon);
            xsA[0] = left[0]; ysA[0] = left[1];
            //extend dm
            left = hua_ExtendLine(xsB[0], ysB[0], xsB[1], ysB[1], -epsilon);
            xsB[0] = left[0]; ysB[0] = left[1];
            // right extend sj
            left = hua_ExtendLine(xsA[xsA.Length - 1], ysA[xsA.Length - 1], xsA[xsA.Length - 2], ysA[xsA.Length - 2], -epsilon);
            xsA[xsA.Length - 1] = left[0]; ysA[xsA.Length - 1] = left[1];
            //right extend dm
            left = hua_ExtendLine(xsB[xsA.Length - 1], ysB[xsA.Length - 1], xsB[xsA.Length - 2], ysB[xsA.Length - 2], -epsilon);
            xsB[xsA.Length - 1] = left[0]; ysB[xsA.Length - 1] = left[1];
        }
        #endregion
        //jd list
        List<(double x, double y, int sji, int dmj)> xys = [];
        for (int i = 0; i < xsA.Length - 1; i++)
        {
            double x1 = xsA[i];
            double y1 = ysA[i];
            double x2 = xsA[i + 1];
            double y2 = ysA[i + 1];
            for (int j = 0; j < xsB.Length - 1; j++)
            {
                double x3 = xsB[j];
                double y3 = ysB[j];
                double x4 = xsB[j + 1];
                double y4 = ysB[j + 1];
                double[] jd;
                jd = hua_Intersections_Segment(x1, y1, x2, y2, x3, y3, x4, y4);
                if (jd.Length > 0)
                {
                    xys.Add((jd[0], jd[1], i, j));
                }
            }
        }
        double[,] xy = new double[xys.Count, 4];
        for (int i = 0; i < xys.Count; i++)
        {
            xy[i, 0] = xys[i].x;
            xy[i, 1] = xys[i].y;
            xy[i, 2] = xys[i].sji;
            xy[i, 3] = xys[i].dmj;
        }
        return xy;
    }
    #endregion


    //直线交点
    [ExcelFunction(Description = "直线交点", Category = "横断面相关")]
    public static double[] hua_Intersections_line(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
    {
        double a1 = y2 - y1;
        double b1 = x1 - x2;
        double c1 = x2 * y1 - x1 * y2;
        double a2 = y4 - y3;
        double b2 = x3 - x4;
        double c2 = x4 * y3 - x3 * y4;
        double denominator = a1 * b2 - a2 * b1;
        if (Math.Abs(denominator) < 1e-10)
            return [0, 0];
        double px = (b1 * c2 - b2 * c1) / denominator;
        double py = (a2 * c1 - a1 * c2) / denominator;
        return [px, py];
    }


    //read txt
    /// <summary>
    /// 读取文本文件并解析为字典（键：文件名，值：二维double数组）
    /// </summary>
    /// <param name="filePath">文本文件路径</param>
    /// <param name="encoding">文件编码（默认UTF-8，需显式传递其他编码）</param>
    /// <returns>包含文件名和解析结果的字典</returns>
    /// <exception cref="FileNotFoundException">文件未找到</exception>
    /// <exception cref="FormatException">数据解析失败</exception>
    public static double[][] ReadTextFileToDictionary(
        string filePath,
        Encoding? encoding = null) // 默认参数改为 null（编译时常量）
    {
        if (!File.Exists(filePath))
        {
            return [];
        }

        string fileName = Path.GetFileName(filePath);
        Encoding actualEncoding = encoding ?? Encoding.UTF8;
        string[] allLines = File.ReadAllLines(filePath, actualEncoding);
        List<double[]> validRows = new List<double[]>();
        for (int lineIndex = 0; lineIndex < allLines.Length; lineIndex++)
        {
            string line = allLines[lineIndex];
            string trimmedLine = line.Trim();
            if (string.IsNullOrWhiteSpace(trimmedLine))
                continue;
            string[] tokens = trimmedLine.Split(
                new[] { ',', ' ' },
                StringSplitOptions.RemoveEmptyEntries
            );
            double[] doubleRow = new double[tokens.Length];
            for (int tokenIndex = 0; tokenIndex < tokens.Length; tokenIndex++)
            {
                if (!double.TryParse(tokens[tokenIndex], out double value))
                    throw new FormatException(
                        $"解析失败 → 文件：{fileName} → 行号：{lineIndex + 1} → 列号：{tokenIndex + 1} → 内容：{tokens[tokenIndex]}"
                    );
                doubleRow[tokenIndex] = value;
            }

            validRows.Add(doubleRow);
        }
        return validRows.ToArray();

    }
}


