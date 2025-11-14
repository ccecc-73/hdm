using hengduanmian;
using Microsoft.Office.Interop.Excel;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
namespace hdm
{
    public class Autocad
    {
        static dynamic acad = null;
       static Autocad()
        {
            acad= AcadHelper.GetActiveAutoCAD();
            acad.Visible = true;
        }
        public static dynamic AddLine(double x,double y,double x1,double y1)
        {
            acad??=acad= AcadHelper.GetActiveAutoCAD();
           return acad.ActiveDocument.ModelSpace.AddLine(new double[] {x,y,0},new double[] {x1,y1,0});
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="textString"></param>
        /// <param name="insertionPoint"></param>
        /// <param name="height"></param>
        /// <param name="Rotate"></param>
        /// <param name="AcAlignment">文字插入点  0:水平居中,垂直靠下  !=0 水平左下</param>
        /// <returns></returns>
        public static dynamic AddText(string textString,double insertX, double insertY,double height,double Rotate, bool alignBottomCenter)
        {
            double[] insertionPoint = [insertX, insertY, 0];
            acad ??= acad = AcadHelper.GetActiveAutoCAD();
            dynamic textObj= acad.ActiveDocument.ModelSpace.AddText(textString, insertionPoint, height);
            if (alignBottomCenter)
            {           
                textObj.Alignment =13;
                textObj.TextAlignmentPoint = insertionPoint;
            }
            textObj.Rotate(insertionPoint, Rotate*3.1415926/180);
            return textObj;
        }

        //Set plineObj = ThisDrawing.ModelSpace.AddLightWeightPolyline(points)
        /// <summary>
        /// 
        /// </summary>
        /// <param name="points">[x,y,x1,y1,x2,y2...]</param>
        /// <returns></returns>
        public static dynamic AddLightWeightPolyline(double[,] points)
        {
            acad ??= acad = AcadHelper.GetActiveAutoCAD();
            List<double> list = new List<double>();
            foreach (double p in points)
                list.Add(p);
            double[]p1=list.ToArray();
            return acad.ActiveDocument.ModelSpace.AddLightWeightPolyline(p1);
        }


        //ThisDrawing.ModelSpace.AddCircle(centerPoint, radius)
        public static dynamic AddCircle(double x, double y, double radius)
        {
            acad ??= acad = AcadHelper.GetActiveAutoCAD();
            return acad.ActiveDocument.ModelSpace.AddCircle(new double[] { x, y, 0 }, radius);
        }

        //ADDLAYER
        public static string AddLayer(string name)
        {
            acad ??= acad = AcadHelper.GetActiveAutoCAD();
            dynamic layer = acad.ActiveDocument.Layers.Add(name);
            acad.ActiveDocument.ActiveLayer = layer;
            return  layer.Name;
        }
        //SendCommand 
        public static string SendCommand(string c)
        {
            acad ??= acad = AcadHelper.GetActiveAutoCAD();
            acad.ActiveDocument.SendCommand(c);
            return c;

        }


    }
}
