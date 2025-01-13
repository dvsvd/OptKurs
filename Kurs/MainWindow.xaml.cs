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

        private void CalculateButton_Click(object sender, RoutedEventArgs e)
        {
            mainVM.ExecuteMainAlgorithmAsync(ATextBox.Text, BTextBox.Text,
                InitialValueTextBox.Text, EndValueTextBox.Text, NumStepsTextBox.Text,
                VelocityErrorCoefTextBox.Text, AmplitudeAmortTextBox.Text, PhaseAmortTextBox.Text,
                PhaseBorderLevelTextBox.Text);
            MainTab.SelectedItem = Plots;
        }
    }
}