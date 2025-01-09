using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Kurs.Calculations
{
    using System;
    using System.Numerics;

    public class LaplaceInversion
    {
        private const double PI = 3.141592653589793238462643; // pi
        private const int MAX_LAPORDER = 60;
        private static int test_case = 1;

        private static void Upperswap(ref double u, double v)
        {
            if (v > u)
            {
                u = v;
            }
        }

        public static double LaplaceInvert(TransFunction F, double t, double tolerance)
        {
            // Variable declaration
            int i, n, m, r; // counters & intermediate array indices
            int M = 40; // order of Taylor Expansion (must be less than MAX_LAPORDER)
            double DeHoogFactor = 4.0; // DeHoog time factor
            double T; // Period of DeHoog Inversion formula
            double gamma; // Integration limit parameter
            Complex h2M, R2M, z, dz, s; // Temporary variables

            Complex[] Fctrl = new Complex[2 * MAX_LAPORDER + 1];
            Complex[,] e = new Complex[2 * MAX_LAPORDER, MAX_LAPORDER];
            Complex[,] q = new Complex[2 * MAX_LAPORDER, MAX_LAPORDER];
            Complex[] d = new Complex[2 * MAX_LAPORDER + 1];
            Complex[] A = new Complex[2 * MAX_LAPORDER + 2];
            Complex[] B = new Complex[2 * MAX_LAPORDER + 2];

            // Calculate period and integration limits
            T = DeHoogFactor * t;
            gamma = -0.5 * Math.Log(tolerance) / T;

            // Calculate F(s) at evaluation points gamma + i * PI / T for i = 0 to 2 * M - 1
            Fctrl[0] = 0.5 * F.W(gamma);
            for (i = 1; i <= 2 * M; i++)
            {
                s = new Complex(gamma, i * PI / T);
                Fctrl[i] = F.W(s);
            }

            // Evaluate e and q
            for (i = 0; i < 2 * M; i++)
            {
                e[i, 0] = 0.0;
                q[i, 1] = Fctrl[i + 1] / Fctrl[i];
            }
            e[2 * M, 0] = 0.0;

            for (r = 1; r <= M - 1; r++) // one minor correction
            {
                for (i = 2 * (M - r); i >= 0; i--)
                {
                    if ((i < 2 * (M - r)) && (r > 1))
                    {
                        q[i, r] = q[i + 1, r - 1] * e[i + 1, r - 1] / e[i, r - 1];
                    }
                    e[i, r] = q[i + 1, r] - q[i, r] + e[i + 1, r - 1];
                }
            }

            // Populate d vector
            d[0] = Fctrl[0];
            for (m = 1; m <= M; m++)
            {
                d[2 * m - 1] = -q[0, m];
                d[2 * m] = -e[0, m];
            }

            // Evaluate A, B
            z = new Complex(Math.Cos(PI * t / T), Math.Sin(PI * t / T));

            A[0] = 0.0; B[0] = 1.0; // A_{-1}, B_{-1} in De Hoog
            A[1] = d[0]; B[1] = 1.0;
            for (n = 2; n <= 2 * M + 1; n++)
            {
                dz = d[n - 1] * z;
                A[n] = A[n - 1] + dz * A[n - 2];
                B[n] = B[n - 1] + dz * B[n - 2];
            }

            // Eqn. 23 in De Hoog
            h2M = 0.5 * (1.0 + z * (d[2 * M - 1] - d[2 * M]));
            R2M = -h2M * (1.0 - Complex.Sqrt(1.0 + (z * d[2 * M] / h2M / h2M)));

            // Eqn. 24 in De Hoog
            A[2 * M + 1] = A[2 * M] + R2M * A[2 * M - 1];
            B[2 * M + 1] = B[2 * M] + R2M * B[2 * M - 1];

            // Final result
            return (A[2 * M] / B[2 * M]).Real; // or use .Magnitude / .Phase depending on need
        }
    }
    public class Calculator
    {
        // Function to compute CorPhi
        internal static double[] CorPhi(double[] phi)
        {

            // Normalize the angles in phi
            for (int i = 1; i < phi.Length; i++)
            {
                int nValue = (int)Math.Round((phi[i] - phi[i - 1]) / (2 * Math.PI));
                if (nValue != 0)
                {
                    phi[i] -= 2 * Math.PI * nValue;
                }
            }
            return phi;
        }

        public static void MainAlgorithm(double[] A, double[] B,
            double omegaStart, double omegaEnd, int numSteps,
            double velocityErrorCoef, double amplitudeAmortText, double phaseAmortText,
            double phaseBorderLevelText)
        {
            Debug.Print("Main Algorithm start");
            // Step 1. Values
            TransFunction WInit = new TransFunction(A, B);
            if (omegaEnd < omegaStart)
            {
                (omegaEnd, omegaStart) = (omegaStart, omegaEnd);
            }
            double step = (omegaEnd - omegaStart) / numSteps;
            List<Complex> WResult = new(numSteps);
            double[] P, Q, Ampl, Phase;
            IEnumerable<double> GetOmega()
            {
                double omega = omegaStart;
                for (int i = 0; i < numSteps; i++)
                {
                    if (omega > omegaEnd)
                        yield break;
                    yield return omega += step;
                }
            }
            foreach (double omega_i in GetOmega())
            {
                WResult.Add(WInit.W(omega_i));
            }
            P = WResult.Select(wi => wi.Real).ToArray();
            Q = WResult.Select(wi => wi.Imaginary).ToArray();
            Ampl = WResult.Select(wi => wi.Magnitude).ToArray();
            Phase = CorPhi(WResult.Select(wi => wi.Phase).ToArray());
            Complex K = WInit.W(0),
                Ki = 1.0 / (K * velocityErrorCoef);
            double[] L = Ampl.Select(Ai => 20 * Math.Log10(Ai)).ToArray();
            List<Complex> htArg = new(numSteps);
            foreach (double omega_i in GetOmega())
            {
                htArg.Add(WInit.W(omega_i) / omega_i);
            }
            // Step 2. Plots
            Debug.Print("Main Algorithm end");
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
