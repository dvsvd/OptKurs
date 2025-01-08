using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    public class Calculator
    {
        public static void MainAlgoritm(double[] A, double[] B,
            double omegaStart, double omegaEnd, int numSteps)
        {

        }
        public class PolynomialPair<T>(T[] result, T[]? remainder = null)
        {
            public T[] Result = result;
            public T[]? Remainder = remainder;
        }
        /// <summary>
        /// Метод деления многочлена на многочлен.
        /// Алгоритм полагает, что порядок элемента - это его индекс
        /// </summary>
        /// <typeparam name="T">Тип элемента</typeparam>
        /// <param name="u">Делимое</param>
        /// <param name="v">Делитель</param>
        /// <returns>Пара многочлена-результата и многочлена-остатка</returns>
        /// <exception cref="ArgumentException"></exception>
        public static PolynomialPair<T> PolynomialDivision<T>(T[] u, T[] v) where T : INumber<T>
        {
            T[] q = new T[u.Length], r = new T[u.Length];
            int k, j, n = u.Length -1, nv = v.Length - 1;
            while(nv >= 0 && v[nv] == T.Zero)
            if (nv < 0)
                throw new ArgumentException("Деление на нулевой многочлен");
            Array.Copy(u, r, u.Length);
            for (k = n - nv; k >= 0; k--)
            {
                q[k] = r[nv + k] / v[nv];
                for (j = nv + k - 1; j >= k; j--)
                    r[j] -= q[k] * v[j - k];
            }
            for (j = nv; j <= n; j++)
                r[j] = T.Zero;
            return new PolynomialPair<T>(u, r);
        }
        /// <summary>
        /// Метод деления многочлена на многочлен.
        /// Алгоритм полагает, что порядок элемента - это его индекс
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="numerator">Делимое</param>
        /// <param name="denominator">Делитель</param>
        /// <returns>Пара многочлена-результата и многочлена-остатка</returns>
        /// <exception cref="ArgumentException">Генерируется, если степень делимого больше степени делителя</exception>
        public static PolynomialPair<T> MyPolynomialDivision<T>(T[] numerator, T[] denominator) where T : INumber<T>
        {
            List<T> ret = [];
            List<T> rem = [];
            List<T> subt = [];
            if (numerator.Length < denominator.Length)
                throw new ArgumentException("Numerator's degree must be greater than or equal to denominator's degree");
            //Полагаем разность в порядках как разность в длинах
            int diff = numerator.Length - denominator.Length;
            for(int i = numerator.Length - 1; i > 0; i--)
            {
                T res = numerator[i] / denominator[i - diff];
                for(int j = 0; j < denominator.Length; j++)
                {
                    subt.Add(denominator[j] * res);
                }

                ret.Add(res);
                subt.Clear();
            }
            return new PolynomialPair<T>([..ret], [..rem]);
        }
    }
}
