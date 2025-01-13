using System.Configuration;
using System.Data;
using System.Runtime.CompilerServices;
using System.Windows;
using System.Windows.Threading;

namespace Kurs
{
    /// <summary>
    /// Interaction logic for App.xaml
    /// </summary>
    public partial class App : Application
    {
        public App()
        {
            DispatcherUnhandledException += (object sender,
                DispatcherUnhandledExceptionEventArgs e) =>
            {
                MessageBox.Show("При работе алгоритма возникла ошибка: " + e.Exception, "Error", MessageBoxButton.OK, MessageBoxImage.Error);
                e.Handled = true;
            };
        }
    }

}
