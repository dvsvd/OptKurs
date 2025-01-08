using Kurs.Calculations;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Kurs
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private async void CalculateButton_Click(object sender, RoutedEventArgs e)
        {
            await Task.Run(() => Calculator.MainAlgoritm(
                UtilityFunctions.ParseString(ATextBox.Text),
                UtilityFunctions.ParseString(BTextBox.Text),
                double.Parse(InitialValueTextBox.Text, System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(EndValueTextBox.Text, System.Globalization.CultureInfo.InvariantCulture),
                int.Parse(NumStepsTextBox.Text, System.Globalization.CultureInfo.InvariantCulture),
                ));
        }
    }
}