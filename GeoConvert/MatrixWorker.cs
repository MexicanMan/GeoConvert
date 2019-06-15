using System;
using System.Collections.Generic;
using System.Text;

namespace GeoConvert.Matrix
{
    public static class MatrixWorker
    {
        /// <summary>
        /// Calculate the inverse of a matrix using Gauss-Jordan elimination
        /// </summary>
        /// <param name="matrix">Matrix to invert</param>
        /// <returns>Inverted matrix</returns>
        public static double[,] Invert(double[,] matrix)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int n, r;
            double scale;

            // Validate the matrix size
            if (rows != cols)
                throw new ArgumentException("Given matrix is not a square!", "matrix");

            n = rows;
            double[,] inverse = new double[rows, cols];

            // Initialize the inverse to the identity
            for (r = 0; r < n; ++r)
                for (int c = 0; c < n; ++c)
                    inverse[r, c] = (r == c) ? 1.0 : 0.0;

            // Process the matrix one column at a time
            for (int c = 0; c < n; ++c)
            {
                // Scale the current row to start with 1
                // Swap rows if the current value is too close to 0.0
                if (Math.Abs(matrix[c, c]) <= 2.0 * double.Epsilon)
                {
                    for (r = c + 1; r < n; ++r)
                        if (Math.Abs(matrix[r, c]) <= 2.0 * double.Epsilon)
                        {
                            RowSwap(matrix, n, c, r);
                            RowSwap(inverse, n, c, r);
                            break;
                        }
                    if (r >= n)
                        throw new Exception("Given matrix is singular!");
                }
                scale = 1.0 / matrix[c, c];
                RowScale(matrix, n, scale, c);
                RowScale(inverse, n, scale, c);

                // Zero out the rest of the column
                for (r = 0; r < n; ++r)
                {
                    if (r != c)
                    {
                        scale = -matrix[r, c];
                        RowScaleAdd(matrix, n, scale, c, r);
                        RowScaleAdd(inverse, n, scale, c, r);
                    }
                }
            }

            return inverse;
        }

        /// <summary>
        /// Swap 2 rows in a matrix
        /// </summary>
        /// <param name="matrix">The matrix to operate on</param>
        /// <param name="n">The size of the matrix</param>
        /// <param name="r0">First row to swap</param>
        /// <param name="r1">Second row to swap</param>
        private static void RowSwap(double[,] matrix, int n, int r0, int r1)
        {
            double tmp;

            for (int i = 0; i < n; ++i)
            {
                tmp = matrix[r0, i];
                matrix[r0, i] = matrix[r1, i];
                matrix[r1, i] = tmp;
            }
        }

        /// <summary>
        /// Perform scale and add a row in a matrix to another row
        /// </summary>
        /// <param name="matrix">The matrix to operate on</param>
        /// <param name="n">The size of the matrix</param>
        /// <param name="a">The scale factor to apply to row <paramref name="r0"/></param>
        /// <param name="r0">The row to scale</param>
        /// <param name="r1">The row to add and store to</param>
        private static void RowScaleAdd(double[,] matrix, int n, double a, int r0, int r1)
        {
            for (int i = 0; i < n; ++i)
                matrix[r1, i] += a * matrix[r0, i];
        }

        /// <summary>
        /// Scale a row in a matrix by a constant factor
        /// </summary>
        /// <param name="matrix">The matrix to operate on</param>
        /// <param name="n">The size of the matrix</param>
        /// <param name="a">The factor to scale row <paramref name="r"/></param>
        /// <param name="r">The row to scale</param>
        private static void RowScale(double[,] matrix, int n, double a, int r)
        {
            for (int i = 0; i < n; ++i)
                matrix[r, i] *= a;
        }

    }
}
