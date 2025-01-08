using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    /// <summary>
    /// Класс, представляющий многочлен комплексного аргумента
    /// </summary>
    /// <param name="coefficients">Массив коэффициентов, где младший индекс - коэф. старшей степени</param>
    public class ComplexPolynomial(double[] coefficients)
    {
        public double[] C { get; private set; } = coefficients;
        /// <summary>
        /// Функция вычисления значения многочлена
        /// </summary>
        /// <param name="s">Комплексный аргумент</param>
        /// <returns>Комплексный результат</returns>
        public Complex Evaluate(Complex s)
        {
            int j;
            Complex p = new Complex(C[j = C.Length - 1], 0.0);
            while (j > 0)
                p = p * s + C[--j];
            return p;
        }
        public static implicit operator ComplexPolynomial(double[] coefficients) => new(coefficients);
        public static explicit operator double[](ComplexPolynomial obj) => obj.C;
    }
}
