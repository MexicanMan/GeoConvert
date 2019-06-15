using GeoConvert;
using GeoConvert.Matrix;
using NUnit.Framework;

namespace Tests
{
    public class Tests
    {
        private const double eps = 1e-8;

        private GeoConvertor geo;

        [SetUp]
        public void TestSetUp()
        {
            geo = new GeoConvertor();
        }

        [Test]
        public void GeoECEFTest()
        {
            double[] llaBefore = new double[] { 55.754066, 37.621734, 153 };

            double[] ecef = geo.ECEFfromLLA(llaBefore[0], llaBefore[1], llaBefore[2]);
            double[] llaAfter = geo.LLAfromECEF(ecef[0], ecef[1], ecef[2]);

            for (int i = 0; i < 3; i++)
                Assert.AreEqual(llaBefore[i], llaAfter[i], eps);
        }

        [Test]
        public void GeoENUTest()
        {
            double[] llaRef = new double[] { 55.753708, 37.620034, 154 };
            double[] llaBefore = new double[] { 55.754066, 37.621734, 153 };

            double[] enu = geo.ENUfromLLA(llaBefore[0], llaBefore[1], llaBefore[2],
                                          llaRef[0], llaRef[1], llaRef[2]);
            double[] llaAfter = geo.LLAfromENU(enu[0], enu[1], enu[2],
                                               llaRef[0], llaRef[1], llaRef[2]);

            for (int i = 0; i < 3; i++)
                Assert.AreEqual(llaBefore[i], llaAfter[i], eps);
        }

        [Test]
        public void MatrixInverseTest()
        {
            double[,] matrix = new double[,] { { 1, 2, 3 }, 
                                               { -3, 2, 1 }, 
                                               { 4, -1, 1 } };
            double[,] realInverted = new double[,] { { 1.5, -2.5, -2 }, 
                                                     { 3.5, -5.5, -5 }, 
                                                     { -2.5, 4.5, 4 } };

            double[,] inverted = MatrixWorker.Invert(matrix);

            int n = realInverted.GetLength(0);
            Assert.AreEqual(n, inverted.GetLength(0));
            Assert.AreEqual(n, inverted.GetLength(1));

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    Assert.AreEqual(realInverted[i, j], inverted[i, j]);
        }
    }
}