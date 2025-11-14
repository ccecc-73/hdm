using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;

namespace hdm
{
    public class CSharpHostObject
    {
        // 同步方法示例：获取消息
        public string GetMessage()
        {
            return "Hello from C#!";
        }

        // 异步方法示例：获取异步数据
        public async Task<string> GetAsyncData()
        {
            await Task.Delay(100); // 模拟异步工作
            return "Async data from C#";
        }

        // 带参数的方法：显示对话框
        public void ShowDialog(string message)
        {
            MessageBox.Show($"JS 说: {message}", "来自网页的提示");
        }

        // 处理复杂数据：相加两个数字
        public int AddNumbers(int a, int b)
        {
            return a + b;
        }

        // 处理 JSON 数据（使用 System.Text.Json）
        public string ProcessJsonData(string jsonData)
        {
            try
            {
                // 使用 System.Text.Json 反序列化
                var data = JsonSerializer.Deserialize<JsonElement>(jsonData);
                string name = data.GetProperty("name").GetString() ?? "Unknown";
                int score = data.GetProperty("score").GetInt32();

                return $"玩家 {name} 的得分是 {score}";
            }
            catch (Exception ex)
            {
                return $"处理 JSON 数据时出错: {ex.Message}";
            }
        }
    }
}
