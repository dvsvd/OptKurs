using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    public class TransFunction(double[] A, double[] B)
    {
        public ComplexPolynomial A { get; private set; } = A;
        public ComplexPolynomial B { get; private set; } = B;
        public Complex W(Complex s)
        {
            return B.Evaluate(s) / A.Evaluate(s);
        }
        public Complex W(double omega)
        {
            return W(new Complex(0.0, omega));
        }
        public Complex[] W(IEnumerable<Complex> ss)
        {
            return ss.Select(W).ToArray();
        }
        public Complex[] W(IEnumerable<double> omegas)
        {
            return W(omegas.Select(omega => new Complex(0.0, omega)).ToArray());
        }
    }
}
