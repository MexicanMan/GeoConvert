using GeoConvert.Matrix;
using System;

namespace GeoConvert
{
    /// <summary>
    /// Convert from LLA to ENU and backwards using WGS84
    /// </summary>
    public class GeoConvertor
    {
        private const double WGS84a = 6378137.0;  // Equatorial radius
        private const double f = 1.0 / 298.257224;  // Inverse flattering
        private const double WGS84b = WGS84a - WGS84a * f;  // Polar radius
        private const double WGS84a2 = WGS84a * WGS84a;
        private const double WGS84b2 = WGS84b * WGS84b;

        /// <summary>
        /// Convert LLA to ECEF
        /// </summary>
        /// <param name="lat">Latitude</param>
        /// <param name="lon">Longitude</param>
        /// <param name="alt">Altitude</param>
        /// <returns>ECEF coordinates (x, y, z)</returns>
        public double[] ECEFfromLLA(double lat, double lon, double alt)
        {
            lat = Math.PI * lat / 180.0;
            lon = Math.PI * lon / 180.0;

            double cosLat = Math.Cos(lat);
            double sinLat = Math.Sin(lat);

            double cosLong = Math.Cos(lon);
            double sinLong = Math.Sin(lon);

            double c = 1 / Math.Sqrt(cosLat * cosLat + (1 - f) * (1 - f) * sinLat * sinLat);
            double s = (1 - f) * (1 - f) * c;

            double x = (WGS84a * c + alt) * cosLat * cosLong;
            double y = (WGS84a * c + alt) * cosLat * sinLong;
            double z = (WGS84a * s + alt) * sinLat;

            return new double[3] { x, y, z };
        }

        /// <summary>
        /// Convert ECEF to LLA
        /// </summary>
        /// <param name="x">X ECEF</param>
        /// <param name="y">Y ECEF</param>
        /// <param name="z">Z ECEF</param>
        /// <returns>LLA coordinates (lat, lon, alt)</returns>
        public double[] LLAfromECEF(double x, double y, double z)
        {
            double ea = Math.Sqrt((WGS84a * WGS84a - WGS84b * WGS84b) / (WGS84a * WGS84a));
            double eb = Math.Sqrt((WGS84a * WGS84a - WGS84b * WGS84b) / (WGS84b * WGS84b));
            double p = Math.Sqrt(x * x + y * y);

            double theta = Math.Atan2(z * WGS84a, p * WGS84b);
            double lon = Math.Atan2(y, x);
            double lat = Math.Atan2(z + eb * eb * WGS84b * Math.Pow(Math.Sin(theta), 3),
                                    p - ea * ea * WGS84a * Math.Pow(Math.Cos(theta), 3));
            double N = WGS84a / Math.Sqrt(1 - ea * ea * Math.Sin(lat) * Math.Sin(lat));
            double alt = p / Math.Cos(lat) - N;

            return new double[3] { lat * (180.0 / Math.PI), lon * (180.0 / Math.PI), alt };
        }

        /// <summary>
        /// Create a transformation matrix used for convert from ENU to ECEF
        /// </summary>
        /// <param name="lat">Latitude of ENU center</param>
        /// <param name="lon">Longitude of ENU center</param>
        /// <param name="alt">Altitude of ENU center</param>
        /// <returns>Transformation matrix</returns>
        private double[,] ECEFfromENUTransform(double lat, double lon, double alt)
        {
            double[] p = ECEFfromLLA(lat, lon, alt);

            lat = Math.PI * lat / 180.0;
            lon = Math.PI * lon / 180.0;
            double sa = Math.Sin(lat);
            double ca = Math.Cos(lat);
            double so = Math.Sin(lon);
            double co = Math.Cos(lon);

            return new double[,] {{ -so, -sa * co, ca * co, p[0] },
                                  { co, -sa * so, ca * so, p[1] },
                                  { 0, ca, sa, p[2] },
                                  { 0, 0, 0, 1 }};
        }

        /// <summary>
        /// Convert ENU to LLA
        /// </summary>
        /// <param name="x">East in ENU</param>
        /// <param name="y">North in ENU</param>
        /// <param name="z">Up in ENU</param>
        /// <param name="reflat">Latitude of ENU coordinate system center</param>
        /// <param name="reflon">Longitude of ENU coordinate system center</param>
        /// <param name="refalt">Altitude of ENU coordinate system center</param>
        /// <returns>LLA coordinates of given ENU point</returns>
        public double[] LLAfromENU(double x, double y, double z, double reflat, double reflon, double refalt)
        {
            double[,] T = ECEFfromENUTransform(reflat, reflon, refalt);
            double ex = T[0, 0] * x + T[0, 1] * y + T[0, 2] * z + T[0, 3];
            double ey = T[1, 0] * x + T[1, 1] * y + T[1, 2] * z + T[1, 3];
            double ez = T[2, 0] * x + T[2, 1] * y + T[2, 2] * z + T[2, 3];

            return LLAfromECEF(ex, ey, ez);
        }

        /// <summary>
        /// Convert LLA to ENU
        /// </summary>
        /// <param name="lat">Latitude</param>
        /// <param name="lon">Longitude</param>
        /// <param name="alt">Altitude</param>
        /// <param name="reflat">Latitude of ENU coordinate system center</param>
        /// <param name="reflon">Longitude of ENU coordinate system center</param>
        /// <param name="refalt">Altitude of ENU coordinate system center</param>
        /// <returns>ENU coordinates of given LLA point</returns>
        public double[] ENUfromLLA(double lat, double lon, double alt, double reflat, double reflon, double refalt)
        {
            double[,] T = MatrixWorker.Invert(ECEFfromENUTransform(reflat, reflon, refalt));
            double[] p = ECEFfromLLA(lat, lon, alt);

            double tx = T[0, 0] * p[0] + T[0, 1] * p[1] + T[0, 2] * p[2] + T[0, 3];
            double ty = T[1, 0] * p[0] + T[1, 1] * p[1] + T[1, 2] * p[2] + T[1, 3];
            double tz = T[2, 0] * p[0] + T[2, 1] * p[1] + T[2, 2] * p[2] + T[2, 3];

            return new double[3] { tx, ty, tz };
        }
    }
}
