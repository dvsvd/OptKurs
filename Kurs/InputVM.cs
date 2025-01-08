using Kurs.Calculations;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;

namespace Kurs
{
    public class InputVM
    {
        public async void ExecuteMainAlgorithmAsync(string AText, string BText,
            string InitialValueText, string EndValueText, string NumStepsText,
            string VelocityErrorCoefText, string AmplitudeAmortText, string PhaseAmortText,
            string PhaseBorderLevelText)
        {
            await Task.Run(() => Calculator.MainAlgoritm(
            UtilityFunctions.ParseString(AText),
                UtilityFunctions.ParseString(BText),
                double.Parse(InitialValueText, System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(EndValueText, System.Globalization.CultureInfo.InvariantCulture),
                int.Parse(NumStepsText, System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(VelocityErrorCoefText, System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(AmplitudeAmortText, System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(PhaseAmortText, System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(PhaseBorderLevelText, System.Globalization.CultureInfo.InvariantCulture)
                ));
        }
    }
}
