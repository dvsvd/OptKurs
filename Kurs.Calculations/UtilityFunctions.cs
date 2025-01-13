using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    static public partial class UtilityFunctions
    {
        static public double[] ParseString(string s)
        {
            string[] nums = s.Replace(',', '.').Split(' ');
            double[] ret = new double[nums.Length];
            for (int i = 0; i < ret.Length; i++)
            {
                ret[i] = double.Parse(nums[i], System.Globalization.CultureInfo.InvariantCulture);
            }
            return ret;
        }
        static public double ParseNumber(string s) {
            DataTable dt = new DataTable();
                return Convert.ToDouble(dt.Compute
                    (s.Replace("Pi", Math.PI.ToString(),
                    StringComparison.CurrentCultureIgnoreCase).Replace(',', '.'), ""));
        }
    }
}
