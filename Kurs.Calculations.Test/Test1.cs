namespace Kurs.Calculations.Test
{
    [TestClass]
    public sealed class Test1
    {
        [ClassInitialize]
        public static void ClassInit(TestContext context)
        {
            // This method is called once for the test class, before any tests of the class are run.
        }

        [ClassCleanup]
        public static void ClassCleanup()
        {
            // This method is called once for the test class, after all tests of the class are run.
        }

        [TestMethod]
        public void TestMethod1()
        {
        }
    }
}
