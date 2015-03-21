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
    /// Calculates the Cholesky decomposition of a matrix.
    /// <para>The Cholesky decomposition of a real symmetric positive-definite
    /// matrix A consists of a lower triangular matrix L with same size such
    /// that: A = LL^T. In a sense, this is the square root of A.</para>
    /// <para>This class is based on the class with similar name from the
    /// <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
    /// following changes:</para>
    /// <list type="bullet">
    /// <item>a <see cref="getLT()"/> method has been added,</item>
    /// <item>the <c>isspd</c> method has been removed, since the constructor of
    /// this class throws a <see cref="NonPositiveDefiniteMatrixException"/> when a
    /// matrix cannot be decomposed,</item>
    /// <item>a <see cref="getDeterminant()"/> method has been added,</item>
    /// <item>the <c>solve</c> method has been replaced by a <see cref="getSolver()"/>
    /// method and the equivalent method provided by the returned
    /// <see cref="DecompositionSolver"/>.</item>
    /// </list>
    /// </summary>
    /// <remarks>
    /// See <a href="http://mathworld.wolfram.com/CholeskyDecomposition.html">MathWorld</a>
    /// <para/>
    /// See <a href="http://en.wikipedia.org/wiki/Cholesky_decomposition">Wikipedia</a>
    /// </remarks>
    public class CholeskyDecomposition
    {
        /// <summary>
        /// Default threshold above which off-diagonal elements are considered too different
        /// and matrix not symmetric.
        /// </summary>
        public const double DEFAULT_RELATIVE_SYMMETRY_THRESHOLD = 1.0e-15;
        
        /// <summary>
        /// Default threshold below which diagonal elements are considered null
        /// and matrix not positive definite.
        /// </summary>
        public const double DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD = 1.0e-10;
        
        /// <summary>
        /// Row-oriented storage for L^T matrix data.
        /// </summary>
        private double[][] lTData;
        
        /// <summary>
        /// Cached value of L.
        /// </summary>
        private RealMatrix cachedL;
        
        /// <summary>
        /// Cached value of LT.
        /// </summary>
        private RealMatrix cachedLT;

        /// <summary>
        /// Calculates the Cholesky decomposition of the given matrix.
        /// <para>
        /// Calling this constructor is equivalent to call
        /// <see cref="CholeskyDecomposition(RealMatrix, double, double)"/> with the
        /// thresholds set to the default values
        /// <see cref="DEFAULT_RELATIVE_SYMMETRY_THRESHOLD"/> and 
        /// <see cref="DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD"/>
        /// </para>
        /// </summary>
        /// <param name="matrix">the matrix to decompose</param>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        /// <exception cref="NonSymmetricMatrixException"> if the matrix is not symmetric.
        /// </exception>
        /// <exception cref="NonPositiveDefiniteMatrixException"> if the matrix is not
        /// strictly positive definite.</exception>
        /// <remarks>
        /// See <see cref="CholeskyDecomposition(RealMatrix, double, double)"/><para/>
        /// See <see cref="DEFAULT_RELATIVE_SYMMETRY_THRESHOLD"/><para/>
        /// See <see cref="DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD"/>
        /// </remarks>
        public CholeskyDecomposition(RealMatrix matrix) : this(matrix, DEFAULT_RELATIVE_SYMMETRY_THRESHOLD, DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD) { }

        /// <summary>
        /// Calculates the Cholesky decomposition of the given matrix.
        /// </summary>
        /// <param name="matrix">the matrix to decompose</param>
        /// <param name="relativeSymmetryThreshold">threshold above which off-diagonal
        /// elements are considered too different and matrix not symmetric</param>
        /// <param name="absolutePositivityThreshold">threshold below which diagonal
        /// elements are considered null and matrix not positive definite</param>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        /// <exception cref="NonSymmetricMatrixException"> if the matrix is not symmetric.
        /// </exception>
        /// <exception cref="NonPositiveDefiniteMatrixException"> if the matrix is not
        /// strictly positive definite.</exception>
        /// <remarks>
        /// See <see cref="CholeskyDecomposition(RealMatrix)"/><para/>
        /// See <see cref="DEFAULT_RELATIVE_SYMMETRY_THRESHOLD"/><para/>
        /// See <see cref="DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD"/>
        /// </remarks>
        public CholeskyDecomposition(RealMatrix matrix, double relativeSymmetryThreshold, double absolutePositivityThreshold)
        {
            if (!matrix.isSquare())
            {
                throw new NonSquareMatrixException(matrix.getRowDimension(),
                                                   matrix.getColumnDimension());
            }

            int order = matrix.getRowDimension();
            lTData = matrix.getData();
            cachedL = null;
            cachedLT = null;

            // check the matrix before transformation
            for (int i = 0; i < order; ++i)
            {
                double[] lI = lTData[i];

                // check off-diagonal elements (and reset them to 0)
                for (int j = i + 1; j < order; ++j)
                {
                    double[] lJ = lTData[j];
                    double lIJ = lI[j];
                    double lJI = lJ[i];
                    double maxDelta =
                        relativeSymmetryThreshold * FastMath.max(FastMath.abs(lIJ), FastMath.abs(lJI));
                    if (FastMath.abs(lIJ - lJI) > maxDelta)
                    {
                        throw new NonSymmetricMatrixException(i, j, relativeSymmetryThreshold);
                    }
                    lJ[i] = 0;
                }
            }

            // transform the matrix
            for (int i = 0; i < order; ++i)
            {

                double[] ltI = lTData[i];

                // check diagonal element
                if (ltI[i] <= absolutePositivityThreshold)
                {
                    throw new NonPositiveDefiniteMatrixException(ltI[i], i, absolutePositivityThreshold);
                }

                ltI[i] = FastMath.sqrt(ltI[i]);
                double inverse = 1.0 / ltI[i];

                for (int q = order - 1; q > i; --q)
                {
                    ltI[q] *= inverse;
                    double[] ltQ = lTData[q];
                    for (int p = q; p < order; ++p)
                    {
                        ltQ[p] -= ltI[q] * ltI[p];
                    }
                }
            }
        }

        /// <summary>
        /// Returns the matrix L of the decomposition.
        /// <para>L is an lower-triangular matrix</para>
        /// </summary>
        /// <returns>the L matrix</returns>
        public RealMatrix getL()
        {
            if (cachedL == null)
            {
                cachedL = getLT().transpose();
            }
            return cachedL;
        }

        /// <summary>
        /// Returns the transpose of the matrix L of the decomposition.
        /// <para>L^T is an upper-triangular matrix</para>
        /// </summary>
        /// <returns>the transpose of the matrix L of the decomposition</returns>
        public RealMatrix getLT()
        {

            if (cachedLT == null)
            {
                cachedLT = MatrixUtils.createRealMatrix(lTData);
            }

            // return the cached matrix
            return cachedLT;
        }

        /// <summary>
        /// Return the determinant of the matrix
        /// </summary>
        /// <returns>determinant of the matrix</returns>
        public double getDeterminant()
        {
            double determinant = 1.0;
            for (int i = 0; i < lTData.Length; ++i)
            {
                double lTii = lTData[i][i];
                determinant *= lTii * lTii;
            }
            return determinant;
        }

        /// <summary>
        /// Get a solver for finding the A &times; X = B solution in least square sense.
        /// </summary>
        /// <returns>a solver</returns>
        public DecompositionSolver getSolver()
        {
            return new Solver(lTData);
        }

        /// <summary>
        /// Specialized solver.
        /// </summary>
        private class Solver : DecompositionSolver
        {
            /// <summary>
            /// Row-oriented storage for L^T matrix data.
            /// </summary>
            private readonly double[][] lTData;

            /// <summary>
            /// Build a solver from decomposed matrix.
            /// </summary>
            /// <param name="lTData">row-oriented storage for L^T matrix data</param>
            internal Solver(double[][] lTData)
            {
                this.lTData = lTData;
            }

            /// <inheritdoc/>
            public Boolean isNonSingular()
            {
                // if we get this far, the matrix was positive definite, hence non-singular
                return true;
            }

            /// <inheritdoc/>
            public RealVector solve(RealVector b)
            {
                int m = lTData.Length;
                if (b.getDimension() != m)
                {
                    throw new DimensionMismatchException(b.getDimension(), m);
                }

                double[] x = b.toArray();

                // Solve LY = b
                for (int j = 0; j < m; j++)
                {
                    double[] lJ = lTData[j];
                    x[j] /= lJ[j];
                    double xJ = x[j];
                    for (int i = j + 1; i < m; i++)
                    {
                        x[i] -= xJ * lJ[i];
                    }
                }

                // Solve LTX = Y
                for (int j = m - 1; j >= 0; j--)
                {
                    x[j] /= lTData[j][j];
                    double xJ = x[j];
                    for (int i = 0; i < j; i++)
                    {
                        x[i] -= xJ * lTData[i][j];
                    }
                }

                return new ArrayRealVector(x, false);
            }

            /// <inheritdoc/>
            public RealMatrix solve(RealMatrix b)
            {
                int m = lTData.Length;
                if (b.getRowDimension() != m)
                {
                    throw new DimensionMismatchException(b.getRowDimension(), m);
                }

                int nColB = b.getColumnDimension();
                double[][] x = b.getData();

                // Solve LY = b
                for (int j = 0; j < m; j++)
                {
                    double[] lJ = lTData[j];
                    double lJJ = lJ[j];
                    double[] xJ = x[j];
                    for (int k = 0; k < nColB; ++k)
                    {
                        xJ[k] /= lJJ;
                    }
                    for (int i = j + 1; i < m; i++)
                    {
                        double[] xI = x[i];
                        double lJI = lJ[i];
                        for (int k = 0; k < nColB; ++k)
                        {
                            xI[k] -= xJ[k] * lJI;
                        }
                    }
                }

                // Solve LTX = Y
                for (int j = m - 1; j >= 0; j--)
                {
                    double lJJ = lTData[j][j];
                    double[] xJ = x[j];
                    for (int k = 0; k < nColB; ++k)
                    {
                        xJ[k] /= lJJ;
                    }
                    for (int i = 0; i < j; i++)
                    {
                        double[] xI = x[i];
                        double lIJ = lTData[i][j];
                        for (int k = 0; k < nColB; ++k)
                        {
                            xI[k] -= xJ[k] * lIJ;
                        }
                    }
                }
                return new Array2DRowRealMatrix(x);
            }

            /// <summary>
            /// Get the inverse of the decomposed matrix.
            /// </summary>
            /// <returns>the inverse matrix.</returns>
            public RealMatrix getInverse()
            {
                return solve(MatrixUtils.createRealIdentityMatrix(lTData.Length));
            }
        }
    }
}