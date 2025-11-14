
using hengduanmian;
using System.Threading.Tasks;
namespace hdm
{
    [ComVisible(true)]
    public class MyRibbon : ExcelRibbon
    {
        Road tcproad { get; set; }
        CustomTaskPane tcp { get; set; }
        UserControl1 tcpweb { get; set; }
        CustomTaskPane tcpwebpane { get; set; }

        public override string GetCustomUI(string RibbonID)
        {
            return RibbonResources.Ribbon;
        }

        public override object? LoadImage(string imageId)
        {
            // This will return the image resource with the name specified in the image='xxxx' tag
            return RibbonResources.ResourceManager.GetObject(imageId);
        }


        public async Task darwCAD(IRibbonControl control)
        {
            string buttonID = control.Id;
            dynamic acad = null;
            try
            {
                acad = AcadHelper.GetActiveAutoCAD();
                acad.Visible = true;
            }
            catch { MessageBox.Show("可能未启动autocad"); }
            switch (buttonID)
            {
                case "drawLine":
                    
                   // acad.ActiveDocument.ModelSpace.AddCircle(center, r);
                    break;
            }
        }
        public async void  OnButtonPressed(IRibbonControl control)
        {
            string buttonID = control.Id;
            dynamic acad = null;
            try
            {
                acad = AcadHelper.GetActiveAutoCAD();
                acad.Visible = true;
            }
            catch {  }

            switch (buttonID)
            {
                #region 画圆
                case "addcircle":
                    if (acad != null)
                    {
                        double[] center = [];
                        double r = 0.2;
                        for (int i = 0; i < 100; i++)
                        {
                            center = [i, 0, 0];
                            acad.ActiveDocument.ModelSpace.AddCircle(center, r);
                            //Thread.Sleep(500);
                        }

                        System.Windows.Forms.MessageBox.Show("ok!");
                    }
                    else
                    {
                        System.Windows.Forms.MessageBox.Show("可能未启动autocad");
                    }
                  
                    break;
                #endregion

                #region 获取点坐标
                case "getpoint":
                    if (acad != null)
                    {
                        //acad.ActiveDocument.ModelSpace.AddCircle(center, r);
                        Application app = (Application)ExcelDnaUtil.Application;
                        var cell = app.ActiveCell;
                        var hang = cell.Row;
                        var lie = cell.Column;
                        object basePoint = new double[] { 0, 0, 0 }; // 可以是 new double[] { 0, 0, 0 } 或 null
                        string prompt = "\n请选择一个点: ";
                        object result = acad.ActiveDocument.Utility.GetPoint(basePoint, prompt);
                        try
                        {
                            while (result != null)
                            {
                                double[] point = (double[])result;
                                app.Cells[hang, lie] = Math.Round(point[0], 3);
                                app.Cells[hang, lie + 1] = Math.Round(point[1], 3);
                                hang += 1;
                                result = acad.ActiveDocument.Utility.GetPoint(basePoint, prompt);
                            }
                        }
                        catch (Exception ex) { }
                    }
                    else
                    {
                        System.Windows.Forms.MessageBox.Show("可能未启动autocad");
                    }

                    break;
                #endregion


                case "road1":
                    //tcproad ??= new Road();
                    //tcp ??= CustomTaskPaneFactory.CreateCustomTaskPane(tcproad, "路线列表");
                    //tcp.Width = 200;
                    //tcp.DockPosition = MsoCTPDockPosition.msoCTPDockPositionLeft;
                    //tcp.Visible = !tcp.Visible;

                    //UserControl1 tcpweb { get; set; }
                    //CustomTaskPane tcpwebpane { get; set; }

                    tcpweb ??= new UserControl1();
                    if (tcpweb.lx)
                    {
                        tcpwebpane ??= CustomTaskPaneFactory.CreateCustomTaskPane(tcpweb, System.IO.Path.GetDirectoryName(MyFunctions.xllPath) + "/2d/index.html");
                    }
                    else
                    {
                        tcpwebpane ??= CustomTaskPaneFactory.CreateCustomTaskPane(tcpweb, "在线");
                    }


                    tcpwebpane.Width = 400;
                    tcpwebpane.DockPosition = MsoCTPDockPosition.msoCTPDockPositionLeft;
                    tcpwebpane.Visible = !tcpwebpane.Visible;

                    break;

                case "btnCallJS":

                    if (UserControl1.Instance != null)
                    {
                        await UserControl1.Instance.CallJavaScriptFunction("mainJavascriptFunction", "Hello from Excel Ribbon!", DateTime.Now.ToString());
                    }
                    else
                    {
                        MessageBox.Show("WebView2 用户控件未初始化。请先确保控件已加载。");
                    }

                    break;

            }

        }



        // Ribbon 按钮回调方法
        public async Task OnCallJSClick(IRibbonControl control)
        {
            if (UserControl1.Instance != null)
            {
                // 调用 JavaScript 中的 mainJavascriptFunction 函数
                await UserControl1.Instance.CallJavaScriptFunction("mainJavascriptFunction", "Hello from Excel Ribbon!", DateTime.Now.ToString());
            }
            else
            {
                MessageBox.Show("WebView2 用户控件未初始化。请先确保控件已加载。");
            }
        }

        public async void OnShowAlertClick(IRibbonControl control)
        {
            await UserControl1.Instance?.CallJavaScriptFunction("showAlert", "这是通过 Excel Ribbon 按钮触发的提示！");
        }

        public async void OnUpdateContentClick(IRibbonControl control)
        {
            string excelData = "从 Excel 中获取的数据: " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss");
            await UserControl1.Instance?.CallJavaScriptFunction("updateContent", excelData);
        }

        public async void OnSendDataClick(IRibbonControl control)
        {
            // 发送复杂数据到 JavaScript
            var sampleData = new
            {
                source = "Excel",
                value = 42,
                timestamp = DateTime.Now,
                items = new[] { "Item1", "Item2", "Item3" }
            };

            await UserControl1.Instance?.CallJavaScriptFunction("processDataFromExcel", sampleData);
        }
    }
}