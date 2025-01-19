using Kurs.Calculations;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Data;

namespace Kurs
{
    internal class DoubleRoundConverter : IValueConverter
    {

        public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if(!value.IsNumber())
                throw new InvalidOperationException();
            return Math.Round(System.Convert.ToDouble(value), 3).ToString();
        }

        public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            return value;
        }
    }
}
