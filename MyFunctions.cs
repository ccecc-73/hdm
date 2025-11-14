using hdm;
using Microsoft.Office.Interop.Excel;
using Microsoft.VisualBasic.Logging;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace hengduanmian;
public class MyFunctions : IExcelAddIn
{
    public static string xllPath;
    #region 填挖面积
    [ExcelFunction(Description = "填挖面积", Category = "横断面相关")]
    public static double[] hdm_CutAndFill(

       [ExcelArgument(Name = "sjx", Description = "设计点x，横向或纵向连续单元格")] double[] xsA,
       [ExcelArgument(Name = "sjy", Description = "设计点y，一维Range")] double[] ysA,
       [ExcelArgument(Name = "dmx", Description = "地面点x，一维Range")] double[] xsB,
       [ExcelArgument(Name = "dmy", Description = "地面点y，一维Range")] double[] ysB,
       [ExcelArgument(Name = "延长长度", Description = "比如1")] double epsilon)
    {
        //dynamic app = Marshal.GetActiveObject("AutoCAD.Application");
        //dynamic doc = app.ActiveDocument;
        return LL.hua_cutAndFillArea(xsA, ysA, xsB, ysB, epsilon);
    }
    #endregion

    #region k2xy
    [ExcelFunction(Description = "已知里程，计算对应的x y", Category = "横断面相关")]
    public static double[] hdm_k2xy(
         [ExcelArgument(Name = "name", Description = "\n第二段\n卡科拉")] string name,
         [ExcelArgument(Name = "k", Description = "桩号")] double k,
         [ExcelArgument(Name = "b", Description = "宽度，左负右正")] double b,
         [ExcelArgument(Name = "z", Description = "右夹角  整数正常输入  \n若90°03分07.8秒 则输入 90.03078")] double z
         )
    {
        if (Data.pqx.TryGetValue(name, out double[][]? pqx))
        {
            return LL.hua_Dantiaoxianludange(pqx, k, b, z);
        }
        else
        {
            pqx = LL.ReadTextFileToDictionary(Data.path + "\\" + name + "\\pqx.txt");
            if (pqx.Length > 0)
            {
                Data.pqx[name] = pqx;
                return LL.hua_Dantiaoxianludange(pqx, k, b, z);
            }
            return [0, 0, 0];
        }
    }
    #endregion

    #region xy2k
    [ExcelFunction(Description = "已知x,y计算对应的k b", Category = "横断面相关")]
    public static double[] hdm_xy2k(
        [ExcelArgument(Name = "name", Description = "路线名称")] string name,
        [ExcelArgument(Name = "x", Description = "x")] double x,
        [ExcelArgument(Name = "y", Description = "y")] double y)
    {
        if (Data.pqx.TryGetValue(name, out double[][] pqx))
        {
            return LL.hua_Fs( x, y, pqx);
        }
        else
        {
            pqx = LL.ReadTextFileToDictionary(Data.path + "\\" + name + "\\pqx.txt");
            if (pqx.Length > 0)
            {
                Data.pqx[name] = pqx;
                return LL.hua_Fs( x, y, pqx);
            }
            return [0, 0];
        }
    }
    #endregion

    #region 中桩高程
    [ExcelFunction(Description = "已知里程，计算对应的h", Category = "横断面相关")]
    public static double hdm_center_h(
         [ExcelArgument(Name = "name", Description = "路线名称")] string name,
        [ExcelArgument(Name = "k", Description = "桩号")] double k)
    {

        if (Data.sqx.TryGetValue(name, out double[][] sqx))
        {
            return LL.H(k, sqx);
        }
        else
        {
            sqx = LL.ReadTextFileToDictionary(Data.path + "\\" + name + "\\sqx.txt");
            if (sqx.Length > 0)
            {
                Data.sqx[name] = sqx;
                return LL.H(k, sqx);
            }
            return 0.0;
        }
    }
    #endregion



    #region slope
    [ExcelFunction(Description = "已知里程，计算对应的 横坡 ", Category = "横断面相关")]
    public static double[] hdm_slope(
          [ExcelArgument(Name = "name", Description = "路线名称")] string name,
        [ExcelArgument(Name = "k", Description = "桩号")] double k)
    {
        if (Data.slope.TryGetValue(name, out double[][] slope))
        { }
        else
        {
            slope = LL.ReadTextFileToDictionary(Data.path + "\\" + name + "\\hp.txt");
            if (slope.Length > 0)
            {
                Data.slope[name] = slope;
            }
        }

        for (var i = 0; i < slope.Length - 1; i++)
        {
            if (k >= slope[i][0] && k <= slope[i + 1][0])
            {
                double L = slope[i][1] + (slope[i + 1][1] - slope[i][1]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                double R = slope[i][2] + (slope[i + 1][2] - slope[i][2]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                return [L, R];
            }
        }
        return [0.0, 0.0];
    }

    #endregion

    #region 结构层宽度
    [ExcelFunction(Description = "结构层宽度", Category = "横断面相关")]
    public static double[] hdm_结构层宽度(
        [ExcelArgument(Name = "name", Description = "路线名称")] string name,
       [ExcelArgument(Name = "k", Description = "里程")] double k,
       [ExcelArgument(Name = "kd", Description = "顶层宽度(左负右正)")] double kd,
       [ExcelArgument(Name = "hd", Description = "距离顶层厚度")] double hd)
    {
        var hp = hdm_slope(name, k);
        var hpL = hp[0];
        var hpR = hp[1];

        var h = hdm_center_h(name, k);
        if (kd < 0)
        {
            var y1 = h - kd * hpL;
            var y2 = h - hd;
            var jd1 = LL.hua_Intersections_line(kd, y1, kd - 2, y1 - 1, 0, y2, -1, y2 + hpL);
            return jd1;
        }
        else
        {
            var y1 = h + kd * hpR;
            var y2 = h - hd;
            var jd1 = LL.hua_Intersections_line(kd, y1, kd + 2, y1 - 1, 0, y2, 1, y2 + hpR);
            return jd1;
        }
    }
    #endregion

    #region 已知折线,点坐标x,计算对应的y
    [ExcelFunction(Description = "填挖判断", Category = "横断面相关")]
    public static double hdm_fromeXgetY(

       [ExcelArgument(Name = "地面数据", Description = "选择连续的三列数据 <桩号,偏距,高程>")] double[,] dmxy,
       [ExcelArgument(Name = "设计K", Description = "计算的设计桩号")] double k,
       [ExcelArgument(Name = "设计宽度", Description = "左负右正")] double b)
    {
        List<double> x = new List<double>();
        List<double> y = new List<double>();
        for (int i = 0; i < dmxy.GetLength(0); i++)
        {
            if (Math.Abs(k - dmxy[i, 0]) <= 0.01)
            {
                x.Add(dmxy[i, 1]);
                y.Add(dmxy[i, 2]);
            }
        }
        return Math.Round(LL.Interpolate(x.ToArray(), y.ToArray(), b), 4);
    }
    #endregion



    #region width
    [ExcelFunction(Description = "已知里程，计算对应的 路面宽度 ", Category = "横断面相关")]
    public static double hdm_width(
          [ExcelArgument(Name = "name", Description = "路线名称")] string name,
        [ExcelArgument(Name = "k", Description = "桩号")] double k)
    {
        if (Data.width.TryGetValue(name, out double[][] width))
        { }
        else
        {
            width = LL.ReadTextFileToDictionary(Data.path + "\\" + name + "\\width.txt");
            if (width.Length > 0)
            {
                Data.width[name] = width;
            }
        }

        for (var i = 0; i < width.Length; i++)
        {
            if (k >= width[i][0] && k <= width[i][1])
            {

                return width[i][2];
            }
        }
        return 0;
    }

    #endregion


    //灵活选择数据

    #region 交点数据转为线元
    [ExcelFunction(Description = "交点数据转为线元", Category = "灵活选择数据")]
    public static double[,] hdm_交点转线元([ExcelArgument(Name = "jds", Description = "[第一行为:jdx,jdy,起点里程,0,0]\n剩余为[jdx,jdy,ls1,ls2,r]")] double[,] jds)
    {
        var pqx = LL.hua_JD2PQX(jds);
        return pqx;

    }
    #endregion



    #region k2xy 自己选择数据
    /// <summary>
    /// 
    /// </summary>
    /// <param name="k"></param>
    /// <param name="b"></param>
    /// <param name="z"></param>
    /// <param name="data"></param>
    /// <returns></returns>
    [ExcelFunction(
        Name = "hdm_k2xy_data",
        HelpTopic = "http://127.0.0.1:8080/help.html#hdm_k2xy_data", // Link to external help
                                                  //HelpTopic = "MyAddFunction.chm!1001"
        Description = "已知里程，计算对应的x y,灵活选择数据", 
        Category = "灵活选择数据"
        )]
    public static double[] hdm_k2xy_data(
         [ExcelArgument(Name = "k", Description = "桩号")] double k,
         [ExcelArgument(Name = "b", Description = "宽度，左负右正")] double b,
         [ExcelArgument(Name = "z", Description = "右夹角  整数正常输入  \n若90°03分07.8秒 则输入 90.03078")] double z,
         [ExcelArgument(Name = "data", Description = "线元数据 可以选择连续单元格\n每行格式为\n1:起点里程\n2:起点x\n3:起点y\n4:十进制起点方位角123°08′02.3″则输入123.08023\n5:线元长度\n6:起点半径正数,直线为0\n7:终点半径\n8:转向:直线=0,左=-1,右=1")] double[,] data
         )
    {
        if (data.GetLength(0) > 0)
        {
            return LL.hua_Dantiaoxianludange(data, k, b, z);
        }
        return [0, 0, 0];
    }
    #endregion

    #region xy2k 自己选择数据
    [ExcelFunction(Description = "已知x,y计算对应的k b,灵活选择数据", Category = "灵活选择数据")]
    public static double[] hdm_xy2k_data(
        [ExcelArgument(Name = "x", Description = "x")] double x,
        [ExcelArgument(Name = "y", Description = "y")] double y,
        [ExcelArgument(Name = "data", Description = "线元数据 可以选择连续单元格\n每行格式为\n1:起点里程\n2:起点x\n3:起点y\n4:十进制起点方位角123°08′02.3″则输入123.08023\n5:线元长度\n6:起点半径正数,直线为0\n7:终点半径\n8:转向:直线=0,左=-1,右=1")] double[,] data
         )
    {
        if (data.GetLength(0) > 0)
        {
            return LL.hua_Fs( x, y,data);
        }
        return [0, 0];
    }
    #endregion

    #region 中桩高程 自己选择数据
    [ExcelFunction(Description = "已知里程，计算对应的h,灵活选择数据", Category = "灵活选择数据")]
    public static double hdm_center_h_data(
    [ExcelArgument(Name = "k", Description = "桩号")] double k,
    [ExcelArgument(Name = "data", Description = "变坡点数据 可以选择连续单元格\n第一和最后一行格式为[变坡点里程,高程,半径输入0\n1:剩余为:[变坡点里程,高程,半径(都是正数)]:")] double[,] data)
    {
        if (data.GetLength(0) > 0)
        {
            return LL.hua_H(k, data);
        }
        return 0.0;
    }

    #endregion




    //地图相关


    #region 在地图中绘制点
    [ExcelFunction(Description = "在地图中绘制点", Category = "地图相关")]
    public static double map_point(
           [ExcelArgument(Name = "lon", Description = "经度十进制")] double lon,
          [ExcelArgument(Name = "lat", Description = "纬度十进制")] double lat)
    {

        if (UserControl1.Instance != null)
        {
            // 调用 JavaScript 中的 mainJavascriptFunction 函数
            UserControl1.Instance.CallJavaScriptFunction("mainJavascriptFunction", lon, lat, "Hello from Excel Ribbon!", DateTime.Now.ToString());
        }
        else
        {
            MessageBox.Show("WebView2 用户控件未初始化。请先确保控件已加载。");
        }
        return 0.0;

    }
    #endregion


    #region proj 投影  自己选择数据
    [ExcelFunction(Description = "已知经纬度，计算对应的lon,lat,灵活选择数据", Category = "地图相关")]
    public static double[] map_proj(
         [ExcelArgument(Name = "proj", Description = "投影方式\n 输入 \"gaosi\" 代表高斯投影\n 输入 'utm' 代码utm投影")] string proj,
         [ExcelArgument(Name = "lon", Description = "十进制 经度")] double lon,
         [ExcelArgument(Name = "lat", Description = "十进制 纬度")] double lat,
         [ExcelArgument(Name = "center", Description = "\n高斯输入中心经度\nutm时=(经度 + 180) / 6) + 1) 舍去小数部分\n一个范围内通常一致")] double center
       )
    {
        if (proj.Contains("gao"))
        {
            return LL.hua_Gauss_proj(lon, lat, center);
        }
        if (proj.Contains("utm"))
        {
            return LL.hua_Utm_proj(lon, lat);
        }
        return [0, 0, 0];
    }
    #endregion
    //double northing, double easting, double longitude, bool isNorthern int zoneNumber
    #region  unproj 自己选择数据
    [ExcelFunction(Description = "已知x,y 计算经纬度", Category = "地图相关")]
    public static double[] map_unproj(
         [ExcelArgument(Name = "proj", Description = "投影方式\n 输入 \"gaosi\" 代表高斯投影\n 输入 'utm' 代码utm投影")] string proj,
         [ExcelArgument(Name = "x", Description = "北坐标")] double x,
         [ExcelArgument(Name = "y", Description = "东坐标")] double y,
         [ExcelArgument(Name = "center", Description = "中心经度,utm时任意")] double center,
         [ExcelArgument(Name = "north", Description = "\nutm才需要\n北半球 输入 true \n 南半球输入 false")] bool north,
         [ExcelArgument(Name = "zoneNumber", Description = "高斯任意,utm时=(经度 + 180) / 6) + 1) 舍去小数部分\n一个范围内通常一致")] int zoneNumber
       )
    {
        if (proj.Contains("gao"))
        {
            return LL.hua_Gauss_unproj(x,y,center);
        }
        if (proj.Contains("utm"))
        {
            return LL.hua_Utm_unproj(x,y,north,zoneNumber);
        }
        return [0, 0, 0];
    }
    #endregion







    #region 四参数计算 自己选择数据  source, double[] target
    [ExcelFunction(Description = "仿射变化,返回 <ΔX,ΔY,θ,scale> ", Category = "地图相关")]
    public static double[] map_4FourParameter_四参数(
         [ExcelArgument(Name = "source", Description = "原坐标 每行为 x,y")] double[,] source,
         [ExcelArgument(Name = "target", Description = "目标坐标 每行为 x,y")] double[,] target
       )
    {
        double[] s=new double[source.Length];
        double[] t = new double[target.Length ];
        if (source.Length == target.Length)
        {
            for (int i = 0; i < source.GetLength(0); i++)
            {
                s[i * 2] = source[i, 0];
                s[i * 2 + 1] = source[i, 1];
                t[i * 2] = target[i, 0];
                t[i * 2 + 1] = target[i, 1];
            }
            return LL.hua_Cs4(s, t);
        }
        return [0, 0, 0];
    }
    #endregion

    #region proj反投影,一步到位
    [ExcelFunction(Description = "已知经纬度，计算对应的xy", Category = "地图相关")]
    public static double[] map_xy2gps(
        [ExcelArgument(Name = "x", Description = "北坐标")] double x,
        [ExcelArgument(Name = "y", Description = "东坐标")] double y,
         [ExcelArgument(Name = "proj", Description = "投影方式\n 输入 \"gaosi\" 代表高斯投影\n 输入 'utm' 代码utm投影")] string proj,
         [ExcelArgument(Name = "lonlat", Description = "十进制经纬度\n每行为[经度,纬度]")] double[,] lonlat,
         [ExcelArgument(Name = "xy", Description = "控制点xy\n每行为[x,y]")] double[,] xy
         //[ExcelArgument(Name = "center", Description = "中心经度\n高斯输入中心经度\nutm时=(经度 + 180) / 6) + 1) 舍去小数部分\n一个范围内通常一致")] double center,
         //[ExcelArgument(Name = "north", Description = "北半球 输入 true \n 南半球输入 false")] bool north
       )
    {
        double[] sorce_xy = new double[lonlat.Length];
        double[] target_xy = new double[lonlat.Length];
        double center = 0;
        //测区平均经度

        for (int i = 0; i < lonlat.GetLength(0); i++)
        {
            center += lonlat[i, 0];
        }
        center = center/ lonlat.GetLength(0);
        if (proj.Contains("gao"))
        {

            for (int i = 0; i < lonlat.GetLength(0); i++)
            {

                var ls = LL.hua_Gauss_proj(lonlat[i, 0], lonlat[i,1], center);
                sorce_xy[i*2] = ls[0];
                sorce_xy[i*2+1]=ls[1];
                target_xy[i * 2] = xy[i,0];
                target_xy[i*2+1]=xy[i,1];
            }
            var cs4 = LL.hua_Cs4( target_xy, sorce_xy);
            var lsxy = LL.hua_FourParameterTransform(x, y, cs4[0], cs4[1], cs4[2], cs4[3]);
            return LL.hua_Gauss_unproj(lsxy[0], lsxy[1],center);
        }
        if (proj.Contains("utm"))
        {
            bool north = true;
            if (lonlat[0, 1] < 0)
            {
                north = false;
            }
            for (int i = 0; i < lonlat.GetLength(0); i++)
            {
                var ls = LL.hua_Utm_proj(lonlat[i, 0], lonlat[i, 1]);
                sorce_xy[i * 2] = ls[0];
                sorce_xy[i * 2 + 1] = ls[1];
                target_xy[i * 2] = xy[i, 0];
                target_xy[i * 2 + 1] = xy[i, 1];
            }
            var cs4 = LL.hua_Cs4(target_xy, sorce_xy);
            var lsxy = LL.hua_FourParameterTransform(x, y, cs4[0], cs4[1], cs4[2], cs4[3]);
            //LL.Utm_unproj(x,y,center,north,zoneNumber);
            int zone = (int)Math.Floor((center + 180) / 6) + 1;
            return [..LL.hua_Utm_unproj(lsxy[0], lsxy[1],north,zone),center];
        }
        return [0, 0, 0];
    }
    #endregion

    #region proj投影,一步到位
    [ExcelFunction(Description = "已知经纬度，计算对应的xy", Category = "地图相关")]
    public static double[] map_gps2xy(
        [ExcelArgument(Name = "lon", Description = "十进制经度")] double lon,
        [ExcelArgument(Name = "lat", Description = "十进制纬度")] double lat,
         [ExcelArgument(Name = "proj", Description = "投影方式\n 输入 \"gaosi\" 代表高斯投影\n 输入 'utm' 代码utm投影")] string proj,
         [ExcelArgument(Name = "lonlat", Description = "控制点对应的经纬度\n每行为[经度,纬度]")] double[,] lonlat,
         [ExcelArgument(Name = "xy", Description = "控制点xy\n每行为[x,y]")] double[,] xy,
         [ExcelArgument(Name = "center", Description ="中心经度\n高斯输入中心经度\nutm时=(经度 + 180) / 6) + 1) 舍去小数部分\n一个范围内通常一致")] double center
       )
    {
        if (proj.Contains("gao"))
        {
            double[] sorce_xy = new double[lonlat.Length];
            double[] target_xy = new double[lonlat.Length];
            for (int i = 0; i < lonlat.GetLength(0); i++)
            {
                var ls = LL.hua_Gauss_proj(lonlat[i, 0], lonlat[i, 1], center);
                sorce_xy[i * 2] = ls[0];
                sorce_xy[i * 2 + 1] = ls[1];
                target_xy[i * 2] = xy[i, 0];
                target_xy[i * 2 + 1] = xy[i, 1];
            }
            var cs4 = LL.hua_Cs4(sorce_xy, target_xy);
            var lsxy = LL.hua_Gauss_proj(lon, lat, center);
            return LL.hua_FourParameterTransform(lsxy[0], lsxy[1], cs4[0], cs4[1], cs4[2], cs4[3]);
        }
        if (proj.Contains("utm"))
        {

            double[] sorce_xy = new double[lonlat.Length];
            double[] target_xy = new double[lonlat.Length];
            for (int i = 0; i < lonlat.GetLength(0); i++)
            {
                var ls = LL.hua_Utm_proj(lonlat[i, 0], lonlat[i, 1]);
                sorce_xy[i * 2] = ls[0];
                sorce_xy[i * 2 + 1] = ls[1];
                target_xy[i * 2] = xy[i, 0];
                target_xy[i * 2 + 1] = xy[i, 1];
            }
            var cs4 = LL.hua_Cs4(sorce_xy, target_xy);
            var lsxy = LL.hua_Utm_proj(lon, lat);
            //
            return LL.hua_FourParameterTransform(lsxy[0], lsxy[1], cs4[0], cs4[1], cs4[2], cs4[3]);
        }
        return [0, 0, 0];
    }
    #endregion

    #region start end
    public void AutoOpen()
    {
        xllPath = ExcelDnaUtil.XllPath;  // 获取.xll完整路径
        Data.path = System.IO.Path.GetDirectoryName(xllPath);
        // MessageBox.Show(Data.path);
        IntelliSenseServer.Install();

    }
    public void AutoClose()
    {
        IntelliSenseServer.Uninstall();
    }
    #endregion




    public static string hdm_mulu()
    {
        return directory.mulu();
    }


}
