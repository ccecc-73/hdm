using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace hdm
{
    public class directory
    {
        public static string mulu()
        {
            string xllPath = ExcelDnaUtil.XllPath;  // 获取.xll完整路径
            xllPath = System.IO.Path.GetDirectoryName(xllPath);  // 提取目录
            string currentDirectory = xllPath;// Directory.GetCurrentDirectory();

            // 创建 DirectoryInfo 对象
            DirectoryInfo dirInfo = new DirectoryInfo(currentDirectory);

            // 获取当前目录下的所有子目录（不包括子目录的子目录）
            DirectoryInfo[] subDirectories = dirInfo.GetDirectories();
            StringBuilder sb = new StringBuilder();
            foreach (DirectoryInfo dir in subDirectories)
            {
              sb.AppendLine (dir.Name);
            }
            return sb.ToString();
        }

    }
}
