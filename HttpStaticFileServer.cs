using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Text;
using System.Threading.Tasks;

namespace hdm
{
    /// <summary>
    /// 轻量级静态文件 HTTP 服务器（基于 HttpListener）
    /// </summary>
    public class HttpStaticFileServer : IDisposable
    {
        private readonly HttpListener _listener;
        private readonly string _webRoot; // Web 内容根目录（如 WebContent）
        private readonly int _port;       // 监听端口
        private Task _requestLoopTask;    // 请求处理循环任务
        private bool _isDisposed;

        /// <summary>
        /// 初始化服务器
        /// </summary>
        /// <param name="webRoot">Web 内容根目录（绝对路径）</param>
        /// <param name="port">监听端口（若为 0 则自动选可用端口）</param>
        public HttpStaticFileServer(string webRoot, int port = 80)
        {
            _webRoot = Path.GetFullPath(webRoot); // 转为绝对路径
            _port = port;

            // 初始化 HttpListener 并注册前缀
            _listener = new HttpListener();
            string prefix = @"http://127.0.0.1/";
            _listener.Prefixes.Add(prefix);
            _listener.Start();

            // 若端口为 0，更新为实际分配的端口
            // if (_port == 0) _port = ((IPEndPoint)_listener.LocalEndpoint).Port;
        }

        /// <summary>
        /// 启动请求处理循环（后台线程）
        /// </summary>
        public void Start()
        {
            _requestLoopTask = Task.Run(ListenForRequestsAsync);
        }

        /// <summary>
        /// 监听请求并处理（异步循环）
        /// </summary>
        private async Task ListenForRequestsAsync()
        {
            try
            {
                while (_listener.IsListening)
                {
                    // 等待客户端请求
                    var context = await _listener.GetContextAsync();
                    // 异步处理请求（不阻塞循环）
                    _ = HandleRequestAsync(context);
                }
            }
            catch (ObjectDisposedException)
            {
                // 正常关闭监听器时触发
            }
        }

        /// <summary>
        /// 处理单个 HTTP 请求
        /// </summary>
        private async Task HandleRequestAsync(HttpListenerContext context)
        {
            var request = context.Request;
            var response = context.Response;

            try
            {
                // 1. 解析请求路径对应的物理文件路径
                string requestPath = request.Url.AbsolutePath.TrimStart('/');
                string filePath = Path.Combine(_webRoot, requestPath);
                filePath = Path.GetFullPath(filePath); // 防止路径遍历攻击（如 ../../../）

                // 校验文件是否在 Web 根目录下（防止越权访问）
                if (!filePath.StartsWith(_webRoot, StringComparison.OrdinalIgnoreCase))
                {
                    response.StatusCode = 403; // 禁止访问
                    await WriteResponseAsync(response, "Forbidden: Access outside web root.");
                    return;
                }

                // 2. 校验文件是否存在
                if (!File.Exists(filePath))
                {
                    response.StatusCode = 404; // 未找到
                    await WriteResponseAsync(response, "Not Found");
                    return;
                }

                // 3. 获取 MIME 类型（根据文件扩展名）
                string mimeType = GetMimeType(filePath);
                response.ContentType = mimeType;

                // 4. 返回文件内容
                using var fileStream = File.OpenRead(filePath);
                await fileStream.CopyToAsync(response.OutputStream);
            }
            catch (Exception ex)
            {
                response.StatusCode = 500; // 服务器内部错误
                await WriteResponseAsync(response, $"Internal Server Error: {ex.Message}");
            }
            finally
            {
                response.Close(); // 必须关闭响应
            }
        }

        /// <summary>
        /// 根据文件扩展名获取 MIME 类型
        /// </summary>
        private string GetMimeType(string filePath)
        {
            var extension = Path.GetExtension(filePath).ToLowerInvariant();
            return extension switch
            {
                ".html" => "text/html",
                ".css" => "text/css",
                ".js" => "application/javascript",
                ".png" => "image/png",
                ".jpg" or ".jpeg" => "image/jpeg",
                ".gif" => "image/gif",
                ".svg" => "image/svg+xml",
                ".json" => "application/json",
                ".ico" => "image/x-icon",
                _ => "application/octet-stream" // 默认二进制流
            };
        }

        /// <summary>
        /// 辅助方法：写入响应内容
        /// </summary>
        private async Task WriteResponseAsync(HttpListenerResponse response, string content)
        {
            byte[] buffer = Encoding.UTF8.GetBytes(content);
            response.ContentLength64 = buffer.Length;
            await response.OutputStream.WriteAsync(buffer, 0, buffer.Length);
        }

        /// <summary>
        /// 释放服务器资源（停止监听）
        /// </summary>
        public void Dispose()
        {
            if (_isDisposed) return;
            _isDisposed = true;

            _listener.Stop(); // 停止监听
            _listener.Close(); // 关闭监听器

            // 等待请求循环结束（可选）
            if (_requestLoopTask != null)
                _requestLoopTask.Wait(TimeSpan.FromSeconds(1));

            GC.SuppressFinalize(this);
        }
    }
}
