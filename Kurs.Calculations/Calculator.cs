using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace Kurs.Calculations
{
    using System;
    using System.Numerics;
    using System.Reflection;
    using System.Runtime.InteropServices;
    using static System.Runtime.InteropServices.JavaScript.JSType;

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

        public static double LaplaceInvert(Func<Complex, Complex> F, double t, double tolerance = 0.1)
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
            if (!((t > 0.0) || (t < 0.0)))
                t = 1e-6;
            T = DeHoogFactor * t;
            gamma = -0.5 * Math.Log(tolerance) / T;

            // Calculate F(s) at evaluation points gamma + i * PI / T for i = 0 to 2 * M - 1
            Fctrl[0] = 0.5 * F(gamma);
            for (i = 1; i <= 2 * M; i++)
            {
                s = new Complex(gamma, i * PI / T);
                Fctrl[i] = F(s);
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
            
            return 1.0 / T * Math.Exp(gamma * t) * (A[2 * M] / B[2 * M]).Real; // or use .Magnitude / .Phase depending on need
        }
    }
    public class Calculator (double[] A, double[] B,
            double omegaStart, double omegaEnd, int numSteps,
            double velocityErrorCoef, double amplitudeAmortText, double phaseAmortText,
            double phaseBorderLevelText)
    {
        public event EventHandler Draw1Event;
        public event EventHandler Draw2Event;
        public TransFunction WInit { get; private set; }
        public double[] Omegas { get; private set; }
        public double[] P { get; private set; }
        public double[] Q { get; private set; }
        public double[] Ampl { get; private set; }
        public double[] Phase { get; private set; }
        public double[] L { get; private set; }
        public double[] ts { get; private set; }
        public double[] h_t { get; private set; }
        public double[] w_t { get; private set; }
        public double t_u { get; private set; }
        public double Tset { get; private set; }
        public double T3set { get; private set; }
        //public double Ttop  { get; private set; }
        public double Sigma  { get; private set; }
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
        public static IEnumerable<double> DoubleRange(double min, double max, double step)
        {
            for (double value = min; value <= max; value += step)
            {
                yield return value;
            }
        }

        public static double DET(ref double[][] M)
        {
            // Bareiss Det Algorithm
            if (M.GetLength(0) != M[0].GetLength(0))
                throw new InvalidOperationException("Matrix is not square");
            int dim = M.GetLength(0);

            if (dim <= 0)
            {
                return 0;
            }

            double sign = 1;

            for (int k = 0; k < dim - 1; k++)
            {
                //Pivot - row swap needed
                if (M[k][k] == 0)
                {
                    int m = 0;
                    for (m = k + 1; m < dim; m++)
                    {
                        if (M[m][k] != 0)
                        {
                            (M[m], M[k]) = (M[k], M[m]);
                            sign = -sign;
                            break;
                        }
                    }

                    //No entries != 0 found in column k -> det = 0
                    if (m == dim)
                    {
                        return 0;
                    }
                }

                //Apply formula
                for (int i = k + 1; i < dim; i++)
                {
                    for (int j = k; j < dim; j++)
                    {
                        M[i][j] = M[k][k] * M[i][j] - M[i][k] * M[k][j];
                        if (k != 0)
                        {
                            M[i][j] /= M[k - 1][k - 1];
                        }
                    }
                }
            }

            return sign * M[dim - 1][dim - 1];
        }

        public static double Overshoot(Func<double, double> h, double W_0, ComplexPolynomial P)
        {
            double h_0 = h(0);
            if (P.Stability == ComplexPolynomial.PolynomialStability.Unstable)
                throw new UnstablePolynomialException("Неустойчивый многочлен в рассчёте перерегулирования");
            double T_top = Math.Abs((3 / P.Eta) * 3);
            double[] m1 = DoubleRange(0, T_top, T_top / 1000).Select(h).ToArray();
            double h_max = m1.Max(), h_min = m1.Min();
            if ((h_0 < W_0) && (W_0 < h_max))
                return (h_max - W_0) / (W_0 - h_0);
            if ((h_0 > W_0) && (W_0 > h_min))
                return (h_min - W_0) / (W_0 - h_0);
            return 0.0;
        }
        private static double SettingTime(Func<double, double> h, double W0, double Winf, ComplexPolynomial P)
        {
            if (P.Stability != ComplexPolynomial.PolynomialStability.Stable)
                throw new UnstablePolynomialException("Неустойчивый многочлен в рассчёте установки времени");

            double T_top = (3 / P.Eta) * 10;              // Верхняя оценка времени
            return TransientPeriod(h, W0, Winf, T_top);
        }

        // Расчет времени установления
        private static double TransientPeriod(Func<double, double> h, double W0, double Winf, double T_top)
        {
            int nOfPoints = Convert.ToInt32(Math.Ceiling(T_top));
            // Определяем функцию для генерации значений
            IEnumerable<double> calculate_t_k(int n, double ttop) {
                    for(int i = 0; i < n; i++) yield return i * ttop / n;
            }
            // Расчет значений t_k
            double[] t_k = calculate_t_k(nOfPoints, T_top).ToArray();
            // Расчет значений h_k
            double[] h_k = t_k.Select(h).ToArray();
            double h_k_max = h_k.Max(),
                h_k_min = h_k.Min(),
                max_v = Math.Max(Math.Abs(Winf - h_k_max), Math.Abs(Winf - h_k_min)),
                delta = 0.05 * max_v;
            for (int i = 0; i < nOfPoints; i++)
            {
                // Возможно здесь надо будет перевернуть логику
                double[] hv = h_k.Skip(i).ToArray().Select(e => e - W0).ToArray();
                if (Math.Max(hv.Max(), -hv.Min()) <= delta)
                    return t_k[i];
            }
            return T_top;
        }
        public void MainAlgorithm()
        {
            Debug.Print("Main Algorithm start");
            // Step 1. Values
            WInit = new TransFunction(A, B); // W(s)
            if (omegaEnd < omegaStart)
            {
                (omegaEnd, omegaStart) = (omegaStart, omegaEnd);
            }
            double step = (omegaEnd - omegaStart) / numSteps;
            Complex[] WResult; // W(omega)
            Omegas = new double[numSteps];
            double omega = omegaStart;
            for (int i = 0; i < numSteps; i++, omega += step)
            {
                if (omega > omegaEnd)
                    break;
                Omegas[i] = omega;
            }
            WResult = WInit.W(Omegas);
            P = WResult.Select(wi => wi.Real).ToArray();
            Q = WResult.Select(wi => wi.Imaginary).ToArray();
            Ampl = WResult.Select(wi => wi.Magnitude).ToArray();
            Phase = CorPhi(WResult.Select(wi => wi.Phase).ToArray());
            Complex K = WInit.W(0),
                Ki = 1.0 / (K * velocityErrorCoef);
            L = Ampl.Select(Ai => 20 * Math.Log10(Ai)).ToArray();
            double h(double t) => LaplaceInversion.LaplaceInvert((s) => WInit.W(s) / s, t);
            double w(double t) => LaplaceInversion.LaplaceInvert((s) => WInit.W(s), t);
            double WInf = h(0);
            try
            {
                Tset = SettingTime(h, WInit.W(0).Real, WInf, A);
                T3set = Tset * 3;
            }
            catch (UnstablePolynomialException)
            {
                T3set = 1000;
            }
            ts = DoubleRange(0, T3set, T3set / 1000).ToArray();
            h_t = ts.Select(h).ToArray();
            w_t = ts.Select(w).ToArray();
            try
            {
                Sigma = Overshoot(h, WInit.W(0).Real, A);
            }
            catch(UnstablePolynomialException)
            {
                Sigma = double.NaN;
            }
            // Step 2. Plots
            Draw1Event?.Invoke(this, EventArgs.Empty);
            // Step 3. Controller parameters
            
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
