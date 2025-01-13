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

namespace Kurs
{

    public class MainVM : INotifyPropertyChanged
    {
        // Данные
        private Calculator calc = new Calculator(null, null, 0, 0, 0, 0, 0, 0, 0);

        public event PropertyChangedEventHandler? PropertyChanged;

        public PlotModel PModel { get; private set; }
        public PlotModel QModel { get; private set; }
        public PlotModel AmplModel { get; private set; }
        public PlotModel PhaseModel { get; private set; }
        public PlotModel LModel { get; private set; }
        public PlotModel PhaseLogModel { get; private set; }
        public PlotModel htModel { get; private set; }
        public PlotModel wtModel { get; private set; }
        public double Tset { get { return calc.Tset; } }
        public double Sigma { get { return calc.Sigma; } }
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
        static private LineSeries CreateLine(IEnumerable<double> x, IEnumerable<double> y, double lineThick, OxyColor colorValue)
        {
            var lineSeries = new LineSeries { StrokeThickness = lineThick, Color = colorValue };
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
