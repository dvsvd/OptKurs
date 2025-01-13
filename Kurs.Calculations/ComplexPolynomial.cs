using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using static Kurs.Calculations.Calculator;

namespace Kurs.Calculations
{
    /// <summary>
    /// Класс, представляющий многочлен комплексного аргумента
    /// </summary>
    /// <param name="coefficients">Массив коэффициентов, где младший индекс - коэф. старшей степени</param>
    public class ComplexPolynomial(double[] coefficients)
    {
        public enum PolynomialStability
        {
            Unstable,
            BorderStable,
            AperiodicBorderStable,
            OscillatingBorderStable,
            Stable
        }
        public ComplexPolynomial(Func<Complex, Complex> f)
        {

        }
        public double[] C { get; private set; } = coefficients;
        public int Length { get { return C.Length; } }
        /// <summary>
        /// Эта
        /// </summary>
        public double Eta { get { return -((Polynomial)this).Roots().Select(r => r.Real).Max(); } }
        /// <summary>
        /// Устойчивость многочлена
        /// </summary>
        public PolynomialStability Stability { get { return RouthHurwitz(); } }
        /// <summary>
        /// Производная многочлена
        /// </summary>
        public ComplexPolynomial Derivative { get
            {
                double[] derivedCoefficients = new double[C.Length - 1];
                for (int i = C.Length - 1; i > 0; i--)
                {
                    derivedCoefficients[i - 1] = C[i] * (C.Length - 1 - i); // производная: cx^n => n*cx^(n-1)
                }
                return new ComplexPolynomial(derivedCoefficients);
            }
        }
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
        public PolynomialStability RouthHurwitz()
        {
            int n = C.Length - 1; // степень полинома
            double[,] M = new double[n, n];
            for (int i = 0; i < M.GetLength(0); i++)
            {
                for (int j = 0; j < M.GetLength(1); j++)
                {
                    int index = n - j * 2 + i - 1;
                    if (!(index < 0 || index > n))
                        M[i, j] = C[index];
                }
            }
            // Подсчет определителей матрицы Гурвица
            // Определители и миноры уже подсчитаны алгоритмом Барейса
            int k = 0;
            double[][] Mjag = M.Cast<double>().ToArray()
                .GroupBy(x => k++ / M.GetLength(1))
                .Select(y => y.ToArray()).ToArray();
            double det = DET(ref Mjag);
            //Определение устойчивости
            double[] dets = new double[Mjag.GetLength(0)];
            for (int i = 0; i < Mjag.GetLength(0); i++)
            {
                dets[i] = Mjag[i][i];
            }
            if (dets.Any(e => e < 0))
                return PolynomialStability.Unstable;                // Полином неустойчив
            if (dets[n - 1] == 0 && dets[n - 2] == 0)
                return PolynomialStability.BorderStable;            // Граница устойчивости
            if (dets[n - 1] == 0)
                return PolynomialStability.AperiodicBorderStable;   // Апериодическая граница устойчивости
            if (dets[n - 2] == 0)
                return PolynomialStability.OscillatingBorderStable; // Колебательная граница устойчивости
            return PolynomialStability.Stable;                      // Полином устойчив
        }            
        
        public static implicit operator ComplexPolynomial(double[] coefficients) => new(coefficients);
        public static explicit operator double[](ComplexPolynomial obj) => obj.C;
        public static implicit operator MathNet.Numerics.Polynomial(ComplexPolynomial obj) => new(obj.C.Reverse());
    }
}
