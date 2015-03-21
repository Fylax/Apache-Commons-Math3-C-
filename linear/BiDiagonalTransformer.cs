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
using Math3.util;
using System;

namespace Math3.linear
{
    /// <summary>
    /// Class transforming any matrix to bi-diagonal shape.
    /// <para>Any m &times; n matrix A can be written as the product of three matrices:
    /// A = U &times; B &times; V^T with U an m &times; m orthogonal matrix,
    /// B an m &times; n bi-diagonal matrix (lower diagonal if m &lt; n, upper diagonal
    /// otherwise), and V an n &times; n orthogonal matrix.</para>
    /// <para>Transformation to bi-diagonal shape is often not a goal by itself, but it is
    /// an intermediate step in more general decomposition algorithms like 
    /// <see cref="SingularValueDecomposition">Singular Value Decomposition</see>. This class
    /// is therefore intended for internal use by the library and is not public. As a consequence
    /// of  this explicitly limited scope, many methods directly returns references to
    /// internal arrays, not copies.</para>
    /// </summary>
    internal class BiDiagonalTransformer
    {
        /// <summary>
        /// Householder vectors.
        /// </summary>
        private readonly double[][] householderVectors;

        /// <summary>
        /// Main diagonal.
        /// </summary>
        private readonly double[] main;

        /// <summary>
        /// Secondary diagonal.
        /// </summary>
        private readonly double[] secondary;

        /// <summary>
        /// Cached value of U.
        /// </summary>
        private RealMatrix cachedU;

        /// <summary>
        /// Cached value of B.
        /// </summary>
        private RealMatrix cachedB;

        /// <summary>
        /// Cached value of V.
        /// </summary>
        private RealMatrix cachedV;

        /// <summary>
        /// Build the transformation to bi-diagonal shape of a matrix.
        /// </summary>
        /// <param name="matrix">the matrix to transform.</param>
        public BiDiagonalTransformer(RealMatrix matrix)
        {
            int m = matrix.getRowDimension();
            int n = matrix.getColumnDimension();
            int p = FastMath.min(m, n);
            householderVectors = matrix.getData();
            main = new double[p];
            secondary = new double[p - 1];
            cachedU = null;
            cachedB = null;
            cachedV = null;

            // transform matrix
            if (m >= n)
            {
                transformToUpperBiDiagonal();
            }
            else
            {
                transformToLowerBiDiagonal();
            }

        }

        /// <summary>
        /// Returns the matrix U of the transform.
        /// <para>U is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the U matrix</returns>
        public RealMatrix getU()
        {
            if (cachedU == null)
            {

                int m = householderVectors.Length;
                int n = householderVectors[0].Length;
                int p = main.Length;
                int diagOffset = (m >= n) ? 0 : 1;
                double[] diagonal = (m >= n) ? main : secondary;
                double[][] ua = new double[m][];

                // fill up the part of the matrix not affected by Householder transforms
                for (int k = m - 1; k >= p; --k)
                {
                    ua[k][k] = 1;
                }

                // build up first part of the matrix by applying Householder transforms
                for (int k = p - 1; k >= diagOffset; --k)
                {
                    double[] hK = householderVectors[k];
                    ua[k][k] = 1;
                    if (hK[k - diagOffset] != 0.0)
                    {
                        for (int j = k; j < m; ++j)
                        {
                            double alpha = 0;
                            for (int i = k; i < m; ++i)
                            {
                                alpha -= ua[i][j] * householderVectors[i][k - diagOffset];
                            }
                            alpha /= diagonal[k - diagOffset] * hK[k - diagOffset];

                            for (int i = k; i < m; ++i)
                            {
                                ua[i][j] += -alpha * householderVectors[i][k - diagOffset];
                            }
                        }
                    }
                }
                if (diagOffset > 0)
                {
                    ua[0][0] = 1;
                }
                cachedU = MatrixUtils.createRealMatrix(ua);
            }

            // return the cached matrix
            return cachedU;

        }

        /// <summary>
        /// Returns the bi-diagonal matrix B of the transform.
        /// </summary>
        /// <returns>the B matrix</returns>
        public RealMatrix getB()
        {
            if (cachedB == null)
            {

                int m = householderVectors.Length;
                int n = householderVectors[0].Length;
                double[][] ba = new double[m][];
                for (int i = 0; i < main.Length; ++i)
                {
                    ba[i][i] = main[i];
                    if (m < n)
                    {
                        if (i > 0)
                        {
                            ba[i][i - 1] = secondary[i - 1];
                        }
                    }
                    else
                    {
                        if (i < main.Length - 1)
                        {
                            ba[i][i + 1] = secondary[i];
                        }
                    }
                }
                cachedB = MatrixUtils.createRealMatrix(ba);
            }

            // return the cached matrix
            return cachedB;

        }

        /// <summary>
        /// Returns the matrix V of the transform.
        /// <para>V is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the V matrix</returns>
        public RealMatrix getV()
        {
            if (cachedV == null)
            {
                int m = householderVectors.Length;
                int n = householderVectors[0].Length;
                int p = main.Length;
                int diagOffset = (m >= n) ? 1 : 0;
                double[] diagonal = (m >= n) ? secondary : main;
                double[][] va = new double[n][];

                // fill up the part of the matrix not affected by Householder transforms
                for (int k = n - 1; k >= p; --k)
                {
                    va[k][k] = 1;
                }

                // build up first part of the matrix by applying Householder transforms
                for (int k = p - 1; k >= diagOffset; --k)
                {
                    double[] hK = householderVectors[k - diagOffset];
                    va[k][k] = 1;
                    if (hK[k] != 0.0)
                    {
                        for (int j = k; j < n; ++j)
                        {
                            double beta = 0;
                            for (int i = k; i < n; ++i)
                            {
                                beta -= va[i][j] * hK[i];
                            }
                            beta /= diagonal[k - diagOffset] * hK[k];

                            for (int i = k; i < n; ++i)
                            {
                                va[i][j] += -beta * hK[i];
                            }
                        }
                    }
                }
                if (diagOffset > 0)
                {
                    va[0][0] = 1;
                }
                cachedV = MatrixUtils.createRealMatrix(va);
            }

            // return the cached matrix
            return cachedV;

        }

        /// <summary>
        /// Get the Householder vectors of the transform.
        /// <para>Note that since this class is only intended for internal use,
        /// it returns directly a reference to its internal arrays, not a copy.</para> 
        /// </summary>
        /// <returns>the main diagonal elements of the B matrix</returns>
        internal double[][] getHouseholderVectorsRef()
        {
            return householderVectors;
        }

        /// <summary>
        /// Get the main diagonal elements of the matrix B of the transform.
        /// <para>Note that since this class is only intended for internal use,
        /// it returns directly a reference to its internal arrays, not a copy.</para>
        /// </summary>
        /// <returns>the main diagonal elements of the B matrix</returns>
        internal double[] getMainDiagonalRef()
        {
            return main;
        }

        /// <summary>
        /// Get the secondary diagonal elements of the matrix B of the transform.
        /// <para>Note that since this class is only intended for internal use,
        /// it returns directly a reference to its internal arrays, not a copy.</para>
        /// </summary>
        /// <returns>the secondary diagonal elements of the B matrix</returns>
        internal double[] getSecondaryDiagonalRef()
        {
            return secondary;
        }

        /// <summary>
        /// Check if the matrix is transformed to upper bi-diagonal.
        /// </summary>
        /// <returns>true if the matrix is transformed to upper bi-diagonal</returns>
        internal Boolean isUpperBiDiagonal()
        {
            return householderVectors.Length >= householderVectors[0].Length;
        }

        /// <summary>
        /// Transform original matrix to upper bi-diagonal form.
        /// <para>Transformation is done using alternate Householder transforms
        /// on columns and rows.</para>
        /// </summary>
        private void transformToUpperBiDiagonal()
        {

            int m = householderVectors.Length;
            int n = householderVectors[0].Length;
            for (int k = 0; k < n; k++)
            {

                //zero-out a column
                double xNormSqr = 0;
                for (int i = k; i < m; ++i)
                {
                    double c = householderVectors[i][k];
                    xNormSqr += c * c;
                }
                double[] hK = householderVectors[k];
                double a = (hK[k] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
                main[k] = a;
                if (a != 0.0)
                {
                    hK[k] -= a;
                    for (int j = k + 1; j < n; ++j)
                    {
                        double alpha = 0;
                        for (int i = k; i < m; ++i)
                        {
                            double[] hI = householderVectors[i];
                            alpha -= hI[j] * hI[k];
                        }
                        alpha /= a * householderVectors[k][k];
                        for (int i = k; i < m; ++i)
                        {
                            double[] hI = householderVectors[i];
                            hI[j] -= alpha * hI[k];
                        }
                    }
                }

                if (k < n - 1)
                {
                    //zero-out a row
                    xNormSqr = 0;
                    for (int j = k + 1; j < n; ++j)
                    {
                        double c = hK[j];
                        xNormSqr += c * c;
                    }
                    double b = (hK[k + 1] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
                    secondary[k] = b;
                    if (b != 0.0)
                    {
                        hK[k + 1] -= b;
                        for (int i = k + 1; i < m; ++i)
                        {
                            double[] hI = householderVectors[i];
                            double beta = 0;
                            for (int j = k + 1; j < n; ++j)
                            {
                                beta -= hI[j] * hK[j];
                            }
                            beta /= b * hK[k + 1];
                            for (int j = k + 1; j < n; ++j)
                            {
                                hI[j] -= beta * hK[j];
                            }
                        }
                    }
                }

            }
        }

        /// <summary>
        /// Transform original matrix to lower bi-diagonal form.
        /// <para>Transformation is done using alternate Householder transforms
        /// on rows and columns.</para>
        /// </summary>
        private void transformToLowerBiDiagonal()
        {

            int m = householderVectors.Length;
            int n = householderVectors[0].Length;
            for (int k = 0; k < m; k++)
            {

                //zero-out a row
                double[] hK = householderVectors[k];
                double xNormSqr = 0;
                for (int j = k; j < n; ++j)
                {
                    double c = hK[j];
                    xNormSqr += c * c;
                }
                double a = (hK[k] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
                main[k] = a;
                if (a != 0.0)
                {
                    hK[k] -= a;
                    for (int i = k + 1; i < m; ++i)
                    {
                        double[] hI = householderVectors[i];
                        double alpha = 0;
                        for (int j = k; j < n; ++j)
                        {
                            alpha -= hI[j] * hK[j];
                        }
                        alpha /= a * householderVectors[k][k];
                        for (int j = k; j < n; ++j)
                        {
                            hI[j] -= alpha * hK[j];
                        }
                    }
                }

                if (k < m - 1)
                {
                    //zero-out a column
                    double[] hKp1 = householderVectors[k + 1];
                    xNormSqr = 0;
                    for (int i = k + 1; i < m; ++i)
                    {
                        double c = householderVectors[i][k];
                        xNormSqr += c * c;
                    }
                    double b = (hKp1[k] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
                    secondary[k] = b;
                    if (b != 0.0)
                    {
                        hKp1[k] -= b;
                        for (int j = k + 1; j < n; ++j)
                        {
                            double beta = 0;
                            for (int i = k + 1; i < m; ++i)
                            {
                                double[] hI = householderVectors[i];
                                beta -= hI[j] * hI[k];
                            }
                            beta /= b * hKp1[k];
                            for (int i = k + 1; i < m; ++i)
                            {
                                double[] hI = householderVectors[i];
                                hI[j] -= beta * hI[k];
                            }
                        }
                    }
                }

            }
        }
    }
}