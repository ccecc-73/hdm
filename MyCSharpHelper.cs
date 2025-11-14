using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace hdm
{
    [System.Runtime.InteropServices.ComVisible(true)]
    [ClassInterface(ClassInterfaceType.AutoDual)]
    //public interface ITCP2;//net8  need
    //[ComDefaultInterface(typeof(ITCP2))]
    public class MyCSharpHelper
    {
 
        // JS 可以同步调用此方法并直接获取返回值
        public void sendLonLat(double lon,double lat,long hang,long lie)
        {
            Application app = (Application)ExcelDnaUtil.Application;
           // var sheet = app.ActiveSheet;
            app.Cells[hang,lie] =lon;
            app.Cells[hang, lie+1] = lat;
        }

        // 一个属性，JS 可以像读取字段一样读取它
        public string Version => "1.0.0";

        // 也可以执行一些操作，比如显示消息框
        public void ShowMessage(string msg)
        {
            MessageBox.Show($"来自网页的消息: {msg}");
        }
    }
}
