using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using hengduanmian;
namespace hdm
{
    public interface ITCP;//net8  need
    [ComDefaultInterface(typeof(ITCP))]
    public partial class Road : UserControl, ITCP
    {
        public Road()
        {
            InitializeComponent();
        }

        private void Road_Load(object sender, EventArgs e)
        {

        }

        private async void button1_Click(object sender, EventArgs e)
        {
            dynamic acad = null;
            try
            {
                acad = AcadHelper.GetActiveAutoCAD();
                acad.Visible = true;
            }

            catch { }
            double[] center = [];
            double r = 0.2;
            for (int i = 0; i < 100; i++)
            {
                //center = [i, 0, 0];
                //acad.ActiveDocument.ModelSpace.AddCircle(center, r);
                //dynamic line = Autocad.addLine(i, 0, i, i + 5);
                //acad.ActiveDocument.Utility.prompt("Press any key..." + line.Length + "\n");
                //Thread.Sleep(100);

            }
            double[] x = [1.396, 9.635, 25.528, 41.176, 41.176, 45.394, 45.394, 85.513, 145.533];
            double[] y = [1.023, 6.69, 7.249, 4.867, -6.936, -6.936, 4.425, 45.869, 11.634];

            int j = 0;
            for (int k = 0; k < 100; k++)
            {
                j = k * 160;
                for (int i = 0; i < x.Length - 1; i++)
                {
                    double x1 = x[i] + j;
                    double y1 = y[i];
                    double x2 = x[i + 1] + j;
                    double y2 = y[i + 1];
                    dynamic line = Autocad.AddLine(x1, y1, x2, y2);
                  //  dynamic txt = Autocad.AddText("第" + i + "个", new double[] { (x1 + x2) / 2.0, (y1 + y2) / 2.0, 0 }, 1, line.Angle, 1);
                     //Autocad.AddText("第" + i + "个", new double[] { (x1 + x2) / 2.0, (y1 + y2) / 2.0, 0 }, 1, 3.1415/2, 0);


                }
                //Thread.Sleep(1000);
                acad.ActiveDocument.Utility.prompt($"第{k}组\n");
            }



        }

        private void Road_Load_1(object sender, EventArgs e)
        {

        }
    }
}
