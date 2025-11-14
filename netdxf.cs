 
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.VisualStyles;


namespace WindowsFormsApp1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
           // InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            // 创建 DXF 文档
            //DxfDocument doc = new DxfDocument();
            //var myLayer = new Layer("MySpecialLayer");
            //myLayer.Color = AciColor.Blue; // Set layer color to blue
            //doc.Layers.Add(myLayer);
            //// 文本参数
            //string textContent = "旋转30° 底部对齐";
            //double height = 2.5;                     // 文本高度
            //double rotationDegrees = 30;             // 旋转角度（30°）
            //double baseY = 0;                       // 底部对齐的基准线 Y 坐标（如 Y=0）

            //// 插入点：X 坐标自由设置，Y 坐标与基准线对齐
            //Vector2 insertionPoint = new Vector2(5, baseY);

            //for (int i = 0; i < 100000; i++)
            //{
            //    // 创建单行文本
            //    Text text = new Text(
            //    text: textContent,
            //    position: insertionPoint,
            //    height: height
            //);
            //    text.Rotation = 30;// rotationDegrees * Math.PI / 180.0;
            //    text.Alignment = TextAlignment.BottomCenter;
            //    text.Layer = myLayer;

            //    doc.Entities.Add(text);
            //}


            //// （可选）添加基准线（Y=0）用于可视化验证
            //Line baseLine = new Line(new Vector2(-10, baseY), new Vector2(10, baseY));
            //doc.Entities.Add(baseLine);
            //for (int i = 0; i < 100000; i++)
            //{
            //    Polyline2D poly = new Polyline2D();
            //    poly.Vertexes.Add(new Polyline2DVertex(5, 10));
            //    poly.Vertexes.Add(new Polyline2DVertex(15, 12));
            //    poly.Vertexes.Add(new Polyline2DVertex(7, 18));
            //    poly.Vertexes.Add(new Polyline2DVertex(5, 62));

            //    doc.Entities.Add(poly);
            //}


            //// 保存 DXF 文件
            //doc.Save("RotatedBottomAlignedText.dxf");
        }
    }
}
