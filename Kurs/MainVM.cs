using Kurs.Calculations;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;
using System.ComponentModel;
using System.Runtime.CompilerServices;
using OxyPlot.Legends;
using MathNet.Numerics;

namespace Kurs
{

    public class MainVM : INotifyPropertyChanged
    {
        // Данные
        private Calculator calc = new(null, null, 0, 0, 0, 0, 0, 0, 0);

        public event PropertyChangedEventHandler? PropertyChanged;

        public PlotModel? PModel { get; private set; }
        public PlotModel? QModel { get; private set; }
        public PlotModel? AmplModel { get; private set; }
        public PlotModel? PhaseModel { get; private set; }
        public PlotModel? LModel { get; private set; }
        public PlotModel? PhaseLogModel { get; private set; }
        public PlotModel? htModel { get; private set; }
        public PlotModel? ht2Model { get; private set; }
        public PlotModel? ht3Model { get; private set; }
        public PlotModel? wtModel { get; private set; }
        public PlotModel? KpModel { get; private set; }
        public PlotModel? KdModel { get; private set; }
        public PlotModel? OsModel { get; private set; }
        public PlotModel? TSetModel { get; private set; }
        public string? RFormula => $"\\frac{{{Ki} + {KpToPrint}s + {KdToPrint}s^{{2}}}}{{s}}";
        public string? WZFormula
        {
            get
            {
                if (calc.WInit == null)
                    return "";
                Polynomial RP = new([Ki, KpToPrint, KdToPrint]);
                Polynomial q, r, WZ;
                (q, r) = (RP * (Polynomial)calc.WInit.B).DivideRemainder(calc.WInit.A);
                WZ = new((q + r).Coefficients.Select(e => Math.Round(e, 3)));
                return $"\\frac{{{WZ}}}{{{WZ + 1}}}";
            }
        }
        public string? RFormula2 => $"\\frac{{{Ki} + {KpToPrint2}s + {KdToPrint2}s^{{2}}}}{{s}}";
        public string? WZFormula2
        {
            get
            {
                if (calc.WInit == null)
                    return "";
                Polynomial RP = new([Ki, KpToPrint2, KdToPrint2]);
                Polynomial q, r, WZ;
                (q, r) = (RP * (Polynomial)calc.WInit.B).DivideRemainder(calc.WInit.A);
                WZ = new((q + r).Coefficients.Select(e => Math.Round(e, 3)));
                return $"\\frac{{{WZ}}}{{{WZ + 1}}}";
            }
        }
        public double Tset { get { return calc.Tset; } }
        public double Sigma { get { return calc.Sigma; } }
        public double Ki { get { return calc.Ki; } }
        public double KpToPrint { get { return calc.KpToPrint; } }
        public double KdToPrint { get { return calc.KdToPrint; } }
        public double KpToPrint2 { get { return calc.KpToPrint2; } }
        public double KdToPrint2 { get { return calc.KdToPrint2; } }
        public double TMinEval { get { return calc.TMinEval; } }
        public double OmegaMinEval { get { return calc.OmegaMinEval; } }
        public double SigmaEnd { get { return calc.SigmaEnd; } }
        public double TSettingEndNon3 { get { return calc.TSettingEndNon3; } }
        public double OmegaMinOvershoot { get { return calc.OmegaMinOvershoot; } }
        public double TSetOvershoot { get { return calc.TSetOvershoot; } }
        public double SigmaOfMinOvershoot { get { return calc.SigmaOfMinOvershoot; } }
        public double[]? TEndRange { get { return calc.TEndRange; } }
        public double[]? TEndRange2 { get { return calc.TEndRange2; } }
        public double[]? h_endless_vals { get { return calc.h_endless_vals; } }
        public double[]? h_endless_omega_min_vals { get { return calc.h_endless_omega_min_vals; } }
        public int ExecutionTime { get { return (calc.EndTime - calc.StartTime) / 1000; } }
        public double Memory => Math.Round((double)GC.GetGCMemoryInfo().TotalCommittedBytes / (1024 * 1024), 1);
        public System.Timers.Timer MemoryTimer { get; private set; } = new(3 * 1000);
        public string TimeText
        {
            get
            {
                return ExecutionTime < 0 ? "Выполняется..." : $"Выполнение заняло {ExecutionTime} с";
            }
        }
        public MainVM()
        {
            MemoryTimer.Elapsed += MemoryTimer_Elapsed;
            MemoryTimer.AutoReset = true;
            MemoryTimer.Enabled = true;
        }

        private void MemoryTimer_Elapsed(object? sender, System.Timers.ElapsedEventArgs e)
        {
            OnPropertyChanged(nameof(Memory));
        }

        public async void ExecuteMainAlgorithmAsync(string AText, string BText,
            string InitialValueText, string EndValueText, string NumStepsText,
            string VelocityErrorCoefText, string AmplitudeAmortText, string PhaseAmortText,
            string PhaseBorderLevelText)
        {
            try
            {
                var v = UtilityFunctions.ParseNumber(InitialValueText);
                calc = new Calculator(UtilityFunctions.ParseString(AText),
                    UtilityFunctions.ParseString(BText),
                    UtilityFunctions.ParseNumber(InitialValueText),
                    UtilityFunctions.ParseNumber(EndValueText),
                    int.Parse(NumStepsText),
                    UtilityFunctions.ParseNumber(VelocityErrorCoefText),
                    UtilityFunctions.ParseNumber(AmplitudeAmortText),
                    UtilityFunctions.ParseNumber(PhaseAmortText),
                    UtilityFunctions.ParseNumber(PhaseBorderLevelText)
                    );
                calc.Draw1Event += Draw1;
                calc.Draw2Event += Draw2;

                await Task.Run(calc.MainAlgorithm);
            }
            catch(Exception ex)
            {
                MessageBox.Show("При работе алгоритма возникла ошибка: " + ex.ToString());
            }
        }
        public void Draw1(object? sender, EventArgs e)
        {
            PModel = new PlotModel { Title = "Вещественная частотная характеристика P(ω)" };
            PModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            PModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            PModel.Series.Add(CreateLine(calc.Omegas, calc.P, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(PModel));
            QModel = new PlotModel { Title = "Мнимая частотная характеристика Q(ω)" };
            QModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            QModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            QModel.Series.Add(CreateLine(calc.Omegas, calc.Q, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(QModel));
            AmplModel = new PlotModel { Title = "Амплитудно-частотная характеристика A(ω)" };
            AmplModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            AmplModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            AmplModel.Series.Add(CreateLine(calc.Omegas, calc.Ampl, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(AmplModel));
            PhaseModel = new PlotModel { Title = "Фазо-частотная характеристика ϕ(ω)" };
            PhaseModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            PhaseModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            PhaseModel.Series.Add(CreateLine(calc.Omegas, calc.Phase, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(PhaseModel));
            var omega_log = Calculator.DoubleRange(0.01, 1.5, (1.5 - 0.01) / 500.0).ToArray();
            LModel = new PlotModel { Title = "Логарифмическая амплитудно-частотная характеристика L(ω)", ClipTitle=false };
            LModel.Axes.Add(new LogarithmicAxis()
            { Position = AxisPosition.Left, Minimum = calc.L.Min(), Maximum = calc.L.Max(), Title = "" });
            LModel.Axes.Add(new LogarithmicAxis()
            { Position = AxisPosition.Bottom, Minimum = omega_log.Min(), Maximum = omega_log.Max(), Title = "" });
            LModel.Series.Add(CreateLine(omega_log, calc.L, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(LModel));
            PhaseLogModel = new PlotModel { Title = "Логарифмическая фазо-частотная характеристика Φ(ω)" };
            PhaseLogModel.Axes.Add(new LogarithmicAxis()
            { Position = AxisPosition.Left, Minimum = calc.Phase.Min(), Maximum = calc.Phase.Max(), Title = "" });
            PhaseLogModel.Axes.Add(new LogarithmicAxis()
            { Position = AxisPosition.Bottom, Minimum = omega_log.Min(), Maximum = omega_log.Max(), Title = "" });
            //PhaseLogModel.Series.Add(CreateLine(omega_log, calc.Phase, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            PhaseLogModel.Series.Add(CreateLine(omega_log, calc.WInit.W(calc.Omegas.Select(Math.Log10)).Select(wi => wi.Phase), 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(PhaseLogModel));
            htModel = new PlotModel { Title = "Переходная характеристика h(t)" };
            htModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            htModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            htModel.Series.Add(CreateLine(calc.ts, calc.h_t, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(htModel));
            wtModel = new PlotModel { Title = "Импульсная характеристика w(t)" };
            wtModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            wtModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            wtModel.Series.Add(CreateLine(calc.ts, calc.w_t, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(wtModel));
            OnPropertyChanged(nameof(Tset));
            OnPropertyChanged(nameof(Sigma));
        }
        public void Draw2(object? sender, EventArgs e)
        {
            KpModel = new PlotModel { Title = "Коэффициент усиления пропорционального регулятора Кп(ω)" };
            KpModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            KpModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            KpModel.Series.Add(CreateLine(calc.Omegas, calc.Kps, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(KpModel));

            KdModel = new PlotModel { Title = "Коэффициент усиления дифференциального регулятора Кд(ω)" };
            KdModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            KdModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            KdModel.Series.Add(CreateLine(calc.Omegas, calc.Kds, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(KdModel));

            TSetModel = new PlotModel { Title = "Оценка времени установления" };
            TSetModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            TSetModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            TSetModel.Series.Add(CreateLine(calc.Omegas, calc.USTs, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0), "Верхняя граница регулирования"));
            TSetModel.Series.Add(CreateLine(calc.Omegas, calc.LSTs, 2.0, OxyColor.FromRgb(0x0, 0x0, 0xD1), "Нижняя граница регулирования"));
            TSetModel.Legends.Add(new Legend() { LegendPlacement = LegendPlacement.Outside,
                LegendPosition = LegendPosition.BottomLeft,
                LegendOrientation = LegendOrientation.Horizontal
            });
            OnPropertyChanged(nameof(TSetModel));

            OsModel = new PlotModel { Title = "Оценка перерегулирования" };
            OsModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            OsModel.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            OsModel.Series.Add(CreateLine(calc.Omegas, calc.OSs, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(OsModel));
            OnPropertyChanged(nameof(Ki));

            ht2Model = new PlotModel { Title = "Переходная характеристика h(t)" };
            ht2Model.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            ht2Model.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            ht2Model.Series.Add(CreateLine(calc.TEndRange, calc.h_endless_vals, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(ht2Model));

            ht3Model = new PlotModel { Title = "Переходная характеристика h(t)" };
            ht3Model.Axes.Add(new LinearAxis() { Position = AxisPosition.Left, Minimum = 0, Title = "" });
            ht3Model.Axes.Add(new LinearAxis() { Position = AxisPosition.Bottom, Minimum = 0, Title = "" });
            ht3Model.Series.Add(CreateLine(calc.TEndRange2, calc.h_endless_omega_min_vals, 2.0, OxyColor.FromRgb(0xD1, 0x0, 0x0)));
            OnPropertyChanged(nameof(ht3Model));

            OnPropertyChanged(nameof(RFormula));
            OnPropertyChanged(nameof(WZFormula));
            OnPropertyChanged(nameof(KdToPrint));
            OnPropertyChanged(nameof(KpToPrint));
            OnPropertyChanged(nameof(SigmaEnd));
            OnPropertyChanged(nameof(TSettingEndNon3));
            OnPropertyChanged(nameof(KdToPrint2));
            OnPropertyChanged(nameof(KpToPrint2));
            OnPropertyChanged(nameof(RFormula2));
            OnPropertyChanged(nameof(WZFormula2));
            OnPropertyChanged(nameof(TSetOvershoot));
            OnPropertyChanged(nameof(SigmaOfMinOvershoot));
            OnPropertyChanged(nameof(TimeText));
        }
        static private LineSeries CreateLine(IEnumerable<double> x, IEnumerable<double> y, double lineThick, OxyColor colorValue, string title = "")
        {
            var lineSeries = new LineSeries { StrokeThickness = lineThick, Color = colorValue, Title = title };
            int trim = Math.Abs(x.Count() - y.Count());
            if (x.Count() > y.Count())
                x = x.SkipLast(trim);
            else if (y.Count() > x.Count())
                y = y.SkipLast(trim);
            for (int i = 0; i < x.Count(); i++)
            {
                lineSeries.Points.Add(new DataPoint(x.ElementAt(i), y.ElementAt(i)));
            }
            return lineSeries;
        }
        private void OnPropertyChanged([CallerMemberName] string propertyName = "")
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
