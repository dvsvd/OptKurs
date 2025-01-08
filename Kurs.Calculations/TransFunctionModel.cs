using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    public class TransFunctionModel(double[] A, double[] B)
    {
        public ComplexPolynomial A { get; private set; } = A;
        public ComplexPolynomial B { get; private set; } = B;
        public Complex W(Complex s)
        {
            return A.Evaluate(s) / B.Evaluate(s);
        }
        public Complex W(double omega)
        {
            return W(new Complex(0.0, omega));
        }
    }
}
