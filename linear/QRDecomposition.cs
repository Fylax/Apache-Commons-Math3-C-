// Licensed to the Apache Software Foundation (ASF) under one or more
// contributor license agreements.  See the NOTICE file distributed with
// this work for additional information regarding copyright ownership.
// The ASF licenses this file to You under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with
// the License.  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
using Math3.exception;
using Math3.util;
using System;

namespace Math3.linear
{
    /// <summary>
    /// Calculates the QR-decomposition of a matrix.
    /// <para>The QR-decomposition of a matrix A consists of two matrices Q and R
    /// that satisfy: A = QR, Q is orthogonal (Q^T Q = I), and R is
    /// upper triangular. If A is m&times;n, Q is m&times;m and R m&times;n.</para>
    /// <para>This class compute the decomposition using Householder reflectors.</para>
    /// <para>For efficiency purposes, the decomposition in packed form is transposed.
    /// This allows inner loop to iterate inside rows, which is much more cache-efficient
    /// in Java.</para>
    /// <para>This class is based on the class with similar name from the
    /// <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
    /// following changes:</para>
    /// <list type="bullet">
    /// <item>a <see cref="getQT()">getQT</see> method has been added,</item>
    /// <item>the <c>solve</c> and <c>isFullRank</c> methods have been replaced
    /// by a <see cref="getSolver()">getSolver</see> method and the equivalent methods
    /// provided by the returned <see cref="DecompositionSolver"/>.</item>
    /// </list>
    /// </summary>
    /// <remarks>
    /// See <a href="http://mathworld.wolfram.com/QRDecomposition.html">MathWorld</a><para/>
    /// See <a href="http://en.wikipedia.org/wiki/QR_decomposition">Wikipedia</a>
    /// </remarks>
    public class QRDecomposition
    {
        /// <summary>
        /// A packed TRANSPOSED representation of the QR decomposition.
        /// <para>The elements BELOW the diagonal are the elements of the UPPER triangular
        /// matrix R, and the rows ABOVE the diagonal are the Householder reflector vectors
        /// from which an explicit form of Q can be recomputed if desired.</para>
        /// </summary>
        private double[][] qrt;
        
        /// <summary>
        /// The diagonal elements of R.
        /// </summary>
        private double[] rDiag;
        
        /// <summary>
        /// Cached value of Q.
        /// </summary>
        private RealMatrix cachedQ;
        
        /// <summary>
        /// Cached value of QT.
        /// </summary>
        private RealMatrix cachedQT;
        
        /// <summary>
        /// Cached value of R.
        /// </summary>
        private RealMatrix cachedR;
        
        /// <summary>
        /// Cached value of H.
        /// </summary>
        private RealMatrix cachedH;
        
        /// <summary>
        /// Singularity threshold.
        /// </summary>
        private readonly double threshold;

        /// <summary>
        /// Calculates the QR-decomposition of the given matrix.
        /// The singularity threshold defaults to zero.
        /// </summary>
        /// <param name="matrix">The matrix to decompose.</param>
        /// <remarks>
        /// See <see cref="QRDecomposition(RealMatrix,double)"/>
        /// </remarks>
        public QRDecomposition(RealMatrix matrix) : this(matrix, 0d) { }

        /// <summary>
        /// Calculates the QR-decomposition of the given matrix.
        /// </summary>
        /// <param name="matrix">The matrix to decompose.</param>
        /// <param name="threshold">Singularity threshold.</param>
        public QRDecomposition(RealMatrix matrix,
                               double threshold)
        {
            this.threshold = threshold;

            int m = matrix.getRowDimension();
            int n = matrix.getColumnDimension();
            qrt = matrix.transpose().getData();
            rDiag = new double[FastMath.min(m, n)];
            cachedQ = null;
            cachedQT = null;
            cachedR = null;
            cachedH = null;

            decompose(qrt);

        }

        /// <summary>
        /// Decompose matrix.
        /// </summary>
        /// <param name="matrix">transposed matrix</param>
        protected void decompose(double[][] matrix)
        {
            for (int minor = 0; minor < FastMath.min(qrt.Length, qrt[0].Length); minor++)
            {
                performHouseholderReflection(minor, qrt);
            }
        }

        /// <summary>
        /// Perform Householder reflection for a minor A(minor, minor) of A.
        /// </summary>
        /// <param name="minor">minor index</param>
        /// <param name="matrix">transposed matrix</param>
        protected void performHouseholderReflection(int minor, double[][] matrix)
        {

            double[] qrtMinor = qrt[minor];

            /*
             * Let x be the first column of the minor, and a^2 = |x|^2.
             * x will be in the positions qr[minor][minor] through qr[m][minor].
             * The first column of the transformed minor will be (a,0,0,..)'
             * The sign of a is chosen to be opposite to the sign of the first
             * component of x. Let's find a:
             */
            double xNormSqr = 0;
            for (int row = minor; row < qrtMinor.Length; row++)
            {
                double c = qrtMinor[row];
                xNormSqr += c * c;
            }
            double a = (qrtMinor[minor] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
            rDiag[minor] = a;

            if (a != 0.0)
            {

                /*
                 * Calculate the normalized reflection vector v and transform
                 * the first column. We know the norm of v beforehand: v = x-ae
                 * so |v|^2 = <x-ae,x-ae> = <x,x>-2a<x,e>+a^2<e,e> =
                 * a^2+a^2-2a<x,e> = 2a*(a - <x,e>).
                 * Here <x, e> is now qr[minor][minor].
                 * v = x-ae is stored in the column at qr:
                 */
                qrtMinor[minor] -= a; // now |v|^2 = -2a*(qr[minor][minor])

                /*
                 * Transform the rest of the columns of the minor:
                 * They will be transformed by the matrix H = I-2vv'/|v|^2.
                 * If x is a column vector of the minor, then
                 * Hx = (I-2vv'/|v|^2)x = x-2vv'x/|v|^2 = x - 2<x,v>/|v|^2 v.
                 * Therefore the transformation is easily calculated by
                 * subtracting the column vector (2<x,v>/|v|^2)v from x.
                 *
                 * Let 2<x,v>/|v|^2 = alpha. From above we have
                 * |v|^2 = -2a*(qr[minor][minor]), so
                 * alpha = -<x,v>/(a*qr[minor][minor])
                 */
                for (int col = minor + 1; col < qrt.Length; col++)
                {
                    double[] qrtCol = qrt[col];
                    double alpha = 0;
                    for (int row = minor; row < qrtCol.Length; row++)
                    {
                        alpha -= qrtCol[row] * qrtMinor[row];
                    }
                    alpha /= a * qrtMinor[minor];

                    // Subtract the column vector alpha*v from x.
                    for (int row = minor; row < qrtCol.Length; row++)
                    {
                        qrtCol[row] -= alpha * qrtMinor[row];
                    }
                }
            }
        }


        /// <summary>
        /// Returns the matrix R of the decomposition.
        /// <para>R is an upper-triangular matrix</para>
        /// </summary>
        /// <returns>the R matrix</returns>
        public RealMatrix getR()
        {

            if (cachedR == null)
            {

                // R is supposed to be m x n
                int n = qrt.Length;
                int m = qrt[0].Length;
                double[][] ra = new double[m][];
                // copy the diagonal from rDiag and the upper triangle of qr
                for (int row = FastMath.min(m, n) - 1; row >= 0; row--)
                {
                    ra[row][row] = rDiag[row];
                    for (int col = row + 1; col < n; col++)
                    {
                        ra[row][col] = qrt[col][row];
                    }
                }
                cachedR = MatrixUtils.createRealMatrix(ra);
            }

            // return the cached matrix
            return cachedR;
        }

        /// <summary>
        /// Returns the matrix Q of the decomposition.
        /// <para>Q is an orthogonal matrix</para>
        /// </summary>
        /// <returns>the Q matrix</returns>
        public RealMatrix getQ()
        {
            if (cachedQ == null)
            {
                cachedQ = getQT().transpose();
            }
            return cachedQ;
        }

        /// <summary>
        /// Returns the transpose of the matrix Q of the decomposition.
        /// <para>Q is an orthogonal matrix</para> 
        /// </summary>
        /// <returns>the transpose of the Q matrix, Q^T</returns>
        public RealMatrix getQT()
        {
            if (cachedQT == null)
            {

                // QT is supposed to be m x m
                int n = qrt.Length;
                int m = qrt[0].Length;
                double[][] qta = new double[m][];

                /*
                 * Q = Q1 Q2 ... Q_m, so Q is formed by first constructing Q_m and then
                 * applying the Householder transformations Q_(m-1),Q_(m-2),...,Q1 in
                 * succession to the result
                 */
                for (int minor = m - 1; minor >= FastMath.min(m, n); minor--)
                {
                    qta[minor][minor] = 1.0d;
                }

                for (int minor = FastMath.min(m, n) - 1; minor >= 0; minor--)
                {
                    double[] qrtMinor = qrt[minor];
                    qta[minor][minor] = 1.0d;
                    if (qrtMinor[minor] != 0.0)
                    {
                        for (int col = minor; col < m; col++)
                        {
                            double alpha = 0;
                            for (int row = minor; row < m; row++)
                            {
                                alpha -= qta[col][row] * qrtMinor[row];
                            }
                            alpha /= rDiag[minor] * qrtMinor[minor];

                            for (int row = minor; row < m; row++)
                            {
                                qta[col][row] += -alpha * qrtMinor[row];
                            }
                        }
                    }
                }
                cachedQT = MatrixUtils.createRealMatrix(qta);
            }

            // return the cached matrix
            return cachedQT;
        }

        /// <summary>
        /// Returns the Householder reflector vectors.
        /// <para>H is a lower trapezoidal matrix whose columns represent
        /// each successive Householder reflector vector. This matrix is used
        /// to compute Q.</para>
        /// </summary>
        /// <returns>a matrix containing the Householder reflector vectors</returns>
        public RealMatrix getH()
        {
            if (cachedH == null)
            {

                int n = qrt.Length;
                int m = qrt[0].Length;
                double[][] ha = new double[m][];
                for (int i = 0; i < m; ++i)
                {
                    for (int j = 0; j < FastMath.min(i + 1, n); ++j)
                    {
                        ha[i][j] = qrt[j][i] / -rDiag[j];
                    }
                }
                cachedH = MatrixUtils.createRealMatrix(ha);
            }

            // return the cached matrix
            return cachedH;
        }

        /// <summary>
        /// Get a solver for finding the A &times; X = B solution in least square sense.
        /// <para>
        /// Least Square sense means a solver can be computed for an overdetermined system,
        /// (i.e. a system with more equations than unknowns, which corresponds to a tall A
        /// matrix with more rows than columns). In any case, if the matrix is singular
        /// within the tolerance set at <see 
        /// cref="QRDecomposition.QRDecomposition(RealMatrix, double)">construction</see>, 
        /// an error will be triggered when the 
        /// <see cref="DecompositionSolver#solve(RealVector)">solve</see> method will be called.
        /// </para>
        /// </summary>
        /// <returns>a solver</returns>
        public DecompositionSolver getSolver()
        {
            return new Solver(qrt, rDiag, threshold);
        }

        /// <summary>
        /// Specialized solver.
        /// </summary>
        private class Solver : DecompositionSolver
        {
            /// <summary>
            /// A packed TRANSPOSED representation of the QR decomposition.
            /// <para>The elements BELOW the diagonal are the elements of the UPPER triangular
            /// matrix R, and the rows ABOVE the diagonal are the Householder reflector vectors
            /// from which an explicit form of Q can be recomputed if desired.</para>
            /// </summary>
            private readonly double[][] qrt;
            
            /// <summary>
            /// The diagonal elements of R.
            /// </summary>
            private readonly double[] rDiag;
            
            /// <summary>
            /// Singularity threshold.
            /// </summary>
            private readonly double threshold;

            /// <summary>
            /// Build a solver from decomposed matrix.
            /// </summary>
            /// <param name="qrt">Packed TRANSPOSED representation of the QR decomposition.</param>
            /// <param name="rDiag">Diagonal elements of R.</param>
            /// <param name="threshold">Singularity threshold.</param>
            public Solver(double[][] qrt, double[] rDiag, double threshold)
            {
                this.qrt = qrt;
                this.rDiag = rDiag;
                this.threshold = threshold;
            }

            /// <inheritdoc/>
            public Boolean isNonSingular()
            {
                foreach (double diag in rDiag)
                {
                    if (FastMath.abs(diag) <= threshold)
                    {
                        return false;
                    }
                }
                return true;
            }

            /// <inheritdoc/>
            public RealVector solve(RealVector b)
            {
                int n = qrt.Length;
                int m = qrt[0].Length;
                if (b.getDimension() != m)
                {
                    throw new DimensionMismatchException(b.getDimension(), m);
                }
                if (!isNonSingular())
                {
                    throw new SingularMatrixException();
                }

                double[] x = new double[n];
                double[] y = b.toArray();

                // apply Householder transforms to solve Q.y = b
                for (int minor = 0; minor < FastMath.min(m, n); minor++)
                {

                    double[] qrtMinor = qrt[minor];
                    double dotProduct = 0;
                    for (int row = minor; row < m; row++)
                    {
                        dotProduct += y[row] * qrtMinor[row];
                    }
                    dotProduct /= rDiag[minor] * qrtMinor[minor];

                    for (int row = minor; row < m; row++)
                    {
                        y[row] += dotProduct * qrtMinor[row];
                    }
                }

                // solve triangular system R.x = y
                for (int row = rDiag.Length - 1; row >= 0; --row)
                {
                    y[row] /= rDiag[row];
                    double yRow = y[row];
                    double[] qrtRow = qrt[row];
                    x[row] = yRow;
                    for (int i = 0; i < row; i++)
                    {
                        y[i] -= yRow * qrtRow[i];
                    }
                }

                return new ArrayRealVector(x, false);
            }

            /// <inheritdoc/>
            public RealMatrix solve(RealMatrix b)
            {
                int n = qrt.Length;
                int m = qrt[0].Length;
                if (b.getRowDimension() != m)
                {
                    throw new DimensionMismatchException(b.getRowDimension(), m);
                }
                if (!isNonSingular())
                {
                    throw new SingularMatrixException();
                }

                int columns = b.getColumnDimension();
                int blockSize = BlockRealMatrix.BLOCK_SIZE;
                int cBlocks = (columns + blockSize - 1) / blockSize;
                double[][] xBlocks = BlockRealMatrix.createBlocksLayout(n, columns);
                double[][] y = new double[b.getRowDimension()][];
                double[] alpha = new double[blockSize];

                for (int kBlock = 0; kBlock < cBlocks; ++kBlock)
                {
                    int kStart = kBlock * blockSize;
                    int kEnd = FastMath.min(kStart + blockSize, columns);
                    int kWidth = kEnd - kStart;

                    // get the right hand side vector
                    b.copySubMatrix(0, m - 1, kStart, kEnd - 1, y);

                    // apply Householder transforms to solve Q.y = b
                    for (int minor = 0; minor < FastMath.min(m, n); minor++)
                    {
                        double[] qrtMinor = qrt[minor];
                        double factor = 1.0 / (rDiag[minor] * qrtMinor[minor]);

                        for (int i = 0; i < kWidth; ++i)
                        {
                            alpha[i] = 0.0d;
                        }
                        for (int row = minor; row < m; ++row)
                        {
                            double d = qrtMinor[row];
                            double[] yRow = y[row];
                            for (int k = 0; k < kWidth; ++k)
                            {
                                alpha[k] += d * yRow[k];
                            }
                        }
                        for (int k = 0; k < kWidth; ++k)
                        {
                            alpha[k] *= factor;
                        }

                        for (int row = minor; row < m; ++row)
                        {
                            double d = qrtMinor[row];
                            double[] yRow = y[row];
                            for (int k = 0; k < kWidth; ++k)
                            {
                                yRow[k] += alpha[k] * d;
                            }
                        }
                    }

                    // solve triangular system R.x = y
                    for (int j = rDiag.Length - 1; j >= 0; --j)
                    {
                        int jBlock = j / blockSize;
                        int jStart = jBlock * blockSize;
                        double factor = 1.0 / rDiag[j];
                        double[] yJ = y[j];
                        double[] xBlock = xBlocks[jBlock * cBlocks + kBlock];
                        int index = (j - jStart) * kWidth;
                        for (int k = 0; k < kWidth; ++k)
                        {
                            yJ[k] *= factor;
                            xBlock[index++] = yJ[k];
                        }

                        double[] qrtJ = qrt[j];
                        for (int i = 0; i < j; ++i)
                        {
                            double rIJ = qrtJ[i];
                            double[] yI = y[i];
                            for (int k = 0; k < kWidth; ++k)
                            {
                                yI[k] -= yJ[k] * rIJ;
                            }
                        }
                    }
                }

                return new BlockRealMatrix(n, columns, xBlocks, false);
            }

            /// <inheritdoc/>
            /// <exception cref="SingularMatrixException"> if the decomposed matrix is singular.
            /// </exception>
            public RealMatrix getInverse()
            {
                return solve(MatrixUtils.createRealIdentityMatrix(qrt[0].Length));
            }
        }
    }
}