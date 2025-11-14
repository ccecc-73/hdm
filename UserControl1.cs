using hengduanmian;
using Microsoft.Office.Interop.Excel;
using Microsoft.Web.WebView2.WinForms;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;
using System.Windows.Forms;
namespace hdm
{
    public interface ITCP1;//net8  need
    [ComDefaultInterface(typeof(ITCP1))]
    public partial class UserControl1 : UserControl, ITCP1
    {
        // 静态实例用于全局访问（简化示例，生产环境需注意生命周期管理）
        public static UserControl1 Instance { get; private set; }
        private WebView2 webView2;
        public bool lx = false;
        private HttpStaticFileServer? _httpServer; // HTTP 服务器实例
        public UserControl1()
        {
            InitializeComponent();
            Instance = this; // 设置静态实例
            string indexHtml = System.IO.Path.GetDirectoryName(MyFunctions.xllPath) + "/2d/index.html";

            if (File.Exists(indexHtml))
            {
                lx = true;
            }
            else
            {
                lx = false;
            }
            InitializeAsync();
            this.Disposed += WebViewDockableControl_Disposed;

        }

        private async void InitializeAsync()
        {
            //// 等待WebView2核心环境准备就绪
            //await webView21.EnsureCoreWebView2Async(null);
            //// 初始化完成后，导航到网页
            //webView21.CoreWebView2.Navigate("https://ccecc-73.github.io/index.html");
            try
            {
                // 初始化 WebView2 核心组件
                await webView21.EnsureCoreWebView2Async(null);

                // 关键步骤：将 C# 对象注入 JS 上下文，并命名为 "myCSharpObj"
                webView21.CoreWebView2.AddHostObjectToScript("myCSharpObj", new MyCSharpHelper());



                string webRoot = System.IO.Path.GetDirectoryName(MyFunctions.xllPath) + "/2d";// index.html";

                if (Directory.Exists(webRoot))
                {
                    _httpServer = new HttpStaticFileServer(webRoot);
                    _httpServer.Start();
                    webView21.CoreWebView2.Navigate(" http://127.0.0.1/index.html");
                }
                else
                {

                    webView21.CoreWebView2.Navigate("https://ccecc-73.github.io/index.html");
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"WebView2 初始化失败：{ex.Message}", "错误",
                    MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        /// <summary>
        /// 窗体关闭时释放 HTTP 服务器
        /// </summary>
        private void WebViewDockableControl_Disposed(object? sender, EventArgs e)
        {
            _httpServer?.Dispose(); // 停止服务器
        }

        private void webView21_Click(object sender, EventArgs e)
        {

        }


        /// <summary>
        /// 供 Ribbon 按钮调用的公共方法，用于执行 JavaScript 函数。
        /// </summary>
        /// <param name="functionName">要调用的 JavaScript 函数名</param>
        /// <param name="parameters">传递给函数的参数</param>
        public async Task CallJavaScriptFunction(string functionName, params object[] parameters)
        {
            if (webView21?.CoreWebView2 == null)
            {
                MessageBox.Show("WebView2 未就绪，无法调用 JavaScript。");
                return;
            }

            try
            {
                string jsCall;

                if (parameters != null && parameters.Length > 0)
                {
                    // 使用 System.Text.Json 序列化参数
                    string serializedArgs = JsonSerializer.Serialize(parameters);
                    // 使用扩展运算符 ... 将数组展开为参数列表
                    jsCall = $"{functionName}(...{serializedArgs});";
                }
                else
                {
                    jsCall = $"{functionName}();";
                }

                // 执行 JavaScript 代码
                string result = await   webView21.CoreWebView2.ExecuteScriptAsync(jsCall);

                // 处理返回值（result 是 JSON 字符串）
                if (!string.IsNullOrEmpty(result) && result != "null")
                {
                    System.Diagnostics.Debug.WriteLine($"JavaScript 返回值: {result}");
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"调用 JavaScript 函数 '{functionName}' 时出错: {ex.Message}");
            }
        }
        /// <summary>
        /// 导航到指定 URL
        /// </summary>
        public void Navigate(string url)
        {
            if (webView21?.CoreWebView2 != null)
            {
                webView21.CoreWebView2.Navigate(url);
            }
        }
        /// <summary>
        /// 执行任意 JavaScript 代码
        /// </summary>
        public async Task<string> ExecuteScriptAsync(string script)
        {
            if (webView21?.CoreWebView2 != null)
            {
                return await webView21.CoreWebView2.ExecuteScriptAsync(script);
            }
            return null;
        }
    }
}
