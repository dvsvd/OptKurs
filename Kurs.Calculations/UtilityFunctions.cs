using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    static public class UtilityFunctions
    {
        static public double[] ParseString(string s)
        {
            string[] nums = s.Replace(',', '.').Split(' ');
            double[] ret = new double[nums.Length];
            for(int i = 0; i < ret.Length; i++)
            {
                ret[i] = double.Parse(nums[i], System.Globalization.CultureInfo.InvariantCulture);
            }
            return ret;
        }
        static public class DoubleParser
        {
            enum State
            {
                START,
                DIGIT,
                LETTER_E,
                LETTER_P,
                SIGN,
                ERROR
            };
            static public double Parse(string s)
            {
            State state = State.START;
            int sign = 1;
            double ret = 0.0;
                for(int i = 0; i < s.Length;)
                {
                    switch(state)
                    {
                        case State.START:
                            switch(s[i])
                            {
                                case '-':
                                    sign = -1;
                                    goto case '+';
                                case '+':
                                    i++;
                                    break;
                                case 'P':
                                case 'p':
                                case 'E':
                                case 'e':
                                case '0':
                                case '1':
                                case '2':
                                case '3':
                                case '4':
                                case '5':
                                case '6':
                                case '7':
                                case '8':
                                case '9':

                                    break;
                            }
                            break;
                    }
                }
                return ret;
            }
        }
    }
}
