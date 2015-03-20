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
using Math3.exception.util;
using Math3.util;
using System;

namespace Math3.linear
{
    /// <summary>
    /// Calculates the compact Singular Value Decomposition of a matrix.
    /// <para>
    /// The Singular Value Decomposition of matrix A is a set of three matrices: U,
    /// &Sigma; and V such that A = U &times; &Sigma; &times; V^T. Let A be
    /// a m &times; n matrix, then U is a m &times; p orthogonal matrix, &Sigma; is a
    /// p &times; p diagonal matrix with positive or null elements, V is a p &times;
    /// n orthogonal matrix (hence V<sup>T</sup> is also orthogonal) where
    /// p=min(m,n).
    /// </para>
    /// <para>This class is similar to the class with similar name from the
    /// <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
    /// following changes:</para>
    /// <list type="bullet">
    /// <item>the <c>norm2</c> method which has been renamed as <see cref="getNorm()">
    /// getNorm</see>,</item>
    /// <item>the <c>cond</c> method which has been renamed as
    /// <see cref="getConditionNumber()">getConditionNumber</see>,</item>
    /// <item>the <c>rank</c> method which has been renamed as <see cref="getRank()">
    /// getRank</see>,</item>
    /// <item>a <see cref="getUT()">getUT></see> method has been added,</item>
    /// <item>a <see cref="getVT()">getVT"</see> method has been added,</item>
    /// <item>a <see cref="getSolver()">getSolver</see> method has been added,</item>
    /// <item>a <see cref="getCovariance(double)">getCovariance</see> method has been added.</item>
    /// </list>
    /// See <a href="http://mathworld.wolfram.com/SingularValueDecomposition.html">MathWorld</a>
    /// <para/>
    /// See <a href="http://en.wikipedia.org/wiki/Singular_value_decomposition">Wikipedia</a>
    /// </summary>
    public class SingularValueDecomposition
    {
        /// <summary>
        /// Relative threshold for small singular values.
        /// </summary>
        private static readonly double EPS = (1.0 * Math.Pow(2, -52));
        
        /// <summary>
        /// Absolute threshold for small singular values.
        /// </summary>
        private static readonly double TINY = (1.0 * Math.Pow(2, -966));
        
        /// <summary>
        /// Computed singular values.
        /// </summary>
        private readonly double[] singularValues;
        
        /// <summary>
        /// max(row dimension, column dimension).
        /// </summary>
        private readonly int m;
        
        /// <summary>
        /// min(row dimension, column dimension).
        /// </summary>
        private readonly int n;
        
        /// <summary>
        /// Indicator for transposed matrix.
        /// </summary>
        private readonly Boolean transposed;
        
        /// <summary>
        /// Cached value of U matrix.
        /// </summary>
        private readonly RealMatrix cachedU;
        
        /// <summary>
        /// Cached value of transposed U matrix.
        /// </summary>
        private RealMatrix cachedUt;
        
        /// <summary>
        /// Cached value of S (diagonal) matrix.
        /// </summary>
        private RealMatrix cachedS;
        
        /// <summary>
        /// Cached value of V matrix.
        /// </summary>
        private readonly RealMatrix cachedV;
        
        /// <summary>
        /// Cached value of transposed V matrix.
        /// </summary>
        private RealMatrix cachedVt;
        
        /// <summary>
        /// Tolerance value for small singular values, calculated once we have
        /// populated "singularValues".
        /// </summary>
        private readonly double tol;

        /// <summary>
        /// Calculates the compact Singular Value Decomposition of the given matrix.
        /// </summary>
        /// <param name="matrix">Matrix to decompose.</param>
        public SingularValueDecomposition(RealMatrix matrix)
        {
            double[][] A;

            // "m" is always the largest dimension.
            if (matrix.getRowDimension() < matrix.getColumnDimension())
            {
                transposed = true;
                A = matrix.transpose().getData();
                m = matrix.getColumnDimension();
                n = matrix.getRowDimension();
            }
            else
            {
                transposed = false;
                A = matrix.getData();
                m = matrix.getRowDimension();
                n = matrix.getColumnDimension();
            }

            singularValues = new double[n];
            double[][] U = new double[m][];
            double[][] V = new double[n][];
            double[] e = new double[n];
            double[] work = new double[m];
            // Reduce A to bidiagonal form, storing the diagonal elements
            // in s and the super-diagonal elements in e.
            int nct = FastMath.min(m - 1, n);
            int nrt = FastMath.max(0, n - 2);
            for (int k = 0; k < FastMath.max(nct, nrt); k++)
            {
                if (k < nct)
                {
                    // Compute the transformation for the k-th column and
                    // place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    singularValues[k] = 0;
                    for (int i = k; i < m; i++)
                    {
                        singularValues[k] = FastMath.hypot(singularValues[k], A[i][k]);
                    }
                    if (singularValues[k] != 0)
                    {
                        if (A[k][k] < 0)
                        {
                            singularValues[k] = -singularValues[k];
                        }
                        for (int i = k; i < m; i++)
                        {
                            A[i][k] /= singularValues[k];
                        }
                        A[k][k] += 1;
                    }
                    singularValues[k] = -singularValues[k];
                }
                for (int j = k + 1; j < n; j++)
                {
                    if (k < nct &&
                        singularValues[k] != 0)
                    {
                        // Apply the transformation.
                        double t = 0;
                        for (int i = k; i < m; i++)
                        {
                            t += A[i][k] * A[i][j];
                        }
                        t = -t / A[k][k];
                        for (int i = k; i < m; i++)
                        {
                            A[i][j] += t * A[i][k];
                        }
                    }
                    // Place the k-th row of A into e for the
                    // subsequent calculation of the row transformation.
                    e[j] = A[k][j];
                }
                if (k < nct)
                {
                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (int i = k; i < m; i++)
                    {
                        U[i][k] = A[i][k];
                    }
                }
                if (k < nrt)
                {
                    // Compute the k-th row transformation and place the
                    // k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0;
                    for (int i = k + 1; i < n; i++)
                    {
                        e[k] = FastMath.hypot(e[k], e[i]);
                    }
                    if (e[k] != 0)
                    {
                        if (e[k + 1] < 0)
                        {
                            e[k] = -e[k];
                        }
                        for (int i = k + 1; i < n; i++)
                        {
                            e[i] /= e[k];
                        }
                        e[k + 1] += 1;
                    }
                    e[k] = -e[k];
                    if (k + 1 < m &&
                        e[k] != 0)
                    {
                        // Apply the transformation.
                        for (int i = k + 1; i < m; i++)
                        {
                            work[i] = 0;
                        }
                        for (int j = k + 1; j < n; j++)
                        {
                            for (int i = k + 1; i < m; i++)
                            {
                                work[i] += e[j] * A[i][j];
                            }
                        }
                        for (int j = k + 1; j < n; j++)
                        {
                            double t = -e[j] / e[k + 1];
                            for (int i = k + 1; i < m; i++)
                            {
                                A[i][j] += t * work[i];
                            }
                        }
                    }

                    // Place the transformation in V for subsequent
                    // back multiplication.
                    for (int i = k + 1; i < n; i++)
                    {
                        V[i][k] = e[i];
                    }
                }
            }
            // Set up the final bidiagonal matrix or order p.
            int p = n;
            if (nct < n)
            {
                singularValues[nct] = A[nct][nct];
            }
            if (m < p)
            {
                singularValues[p - 1] = 0;
            }
            if (nrt + 1 < p)
            {
                e[nrt] = A[nrt][p - 1];
            }
            e[p - 1] = 0;

            // Generate U.
            for (int j = nct; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    U[i][j] = 0;
                }
                U[j][j] = 1;
            }
            for (int k = nct - 1; k >= 0; k--)
            {
                if (singularValues[k] != 0)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        double t = 0;
                        for (int i = k; i < m; i++)
                        {
                            t += U[i][k] * U[i][j];
                        }
                        t = -t / U[k][k];
                        for (int i = k; i < m; i++)
                        {
                            U[i][j] += t * U[i][k];
                        }
                    }
                    for (int i = k; i < m; i++)
                    {
                        U[i][k] = -U[i][k];
                    }
                    U[k][k] = 1 + U[k][k];
                    for (int i = 0; i < k - 1; i++)
                    {
                        U[i][k] = 0;
                    }
                }
                else
                {
                    for (int i = 0; i < m; i++)
                    {
                        U[i][k] = 0;
                    }
                    U[k][k] = 1;
                }
            }

            // Generate V.
            for (int k = n - 1; k >= 0; k--)
            {
                if (k < nrt &&
                    e[k] != 0)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        double t = 0;
                        for (int i = k + 1; i < n; i++)
                        {
                            t += V[i][k] * V[i][j];
                        }
                        t = -t / V[k + 1][k];
                        for (int i = k + 1; i < n; i++)
                        {
                            V[i][j] += t * V[i][k];
                        }
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    V[i][k] = 0;
                }
                V[k][k] = 1;
            }

            // Main iteration loop for the singular values.
            int pp = p - 1;
            while (p > 0)
            {
                int k;
                int kase;
                // Here is where a test for too many iterations would go.
                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.
                // kase = 1     if s(p) and e[k-1] are negligible and k<p
                // kase = 2     if s(k) is negligible and k<p
                // kase = 3     if e[k-1] is negligible, k<p, and
                //              s(k), ..., s(p) are not negligible (qr step).
                // kase = 4     if e(p-1) is negligible (convergence).
                for (k = p - 2; k >= 0; k--)
                {
                    double threshold
                        = TINY + EPS * (FastMath.abs(singularValues[k]) +
                                        FastMath.abs(singularValues[k + 1]));

                    // the following condition is written this way in order
                    // to break out of the loop when NaN occurs, writing it
                    // as "if (FastMath.abs(e[k]) <= threshold)" would loop
                    // indefinitely in case of NaNs because comparison on NaNs
                    // always return false, regardless of what is checked
                    // see issue MATH-947
                    if (!(FastMath.abs(e[k]) > threshold))
                    {
                        e[k] = 0;
                        break;
                    }

                }

                if (k == p - 2)
                {
                    kase = 4;
                }
                else
                {
                    int ks;
                    for (ks = p - 1; ks >= k; ks--)
                    {
                        if (ks == k)
                        {
                            break;
                        }
                        double t = (ks != p ? FastMath.abs(e[ks]) : 0) +
                            (ks != k + 1 ? FastMath.abs(e[ks - 1]) : 0);
                        if (FastMath.abs(singularValues[ks]) <= TINY + EPS * t)
                        {
                            singularValues[ks] = 0;
                            break;
                        }
                    }
                    if (ks == k)
                    {
                        kase = 3;
                    }
                    else if (ks == p - 1)
                    {
                        kase = 1;
                    }
                    else
                    {
                        kase = 2;
                        k = ks;
                    }
                }
                k++;
                // Perform the task indicated by kase.
                switch (kase)
                {
                    // Deflate negligible s(p).
                    case 1:
                        {
                            double f = e[p - 2];
                            e[p - 2] = 0;
                            for (int j = p - 2; j >= k; j--)
                            {
                                double t = FastMath.hypot(singularValues[j], f);
                                double cs = singularValues[j] / t;
                                double sn = f / t;
                                singularValues[j] = t;
                                if (j != k)
                                {
                                    f = -sn * e[j - 1];
                                    e[j - 1] = cs * e[j - 1];
                                }

                                for (int i = 0; i < n; i++)
                                {
                                    t = cs * V[i][j] + sn * V[i][p - 1];
                                    V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                                    V[i][j] = t;
                                }
                            }
                        }
                        break;
                    // Split at negligible s(k).
                    case 2:
                        {
                            double f = e[k - 1];
                            e[k - 1] = 0;
                            for (int j = k; j < p; j++)
                            {
                                double t = FastMath.hypot(singularValues[j], f);
                                double cs = singularValues[j] / t;
                                double sn = f / t;
                                singularValues[j] = t;
                                f = -sn * e[j];
                                e[j] = cs * e[j];

                                for (int i = 0; i < m; i++)
                                {
                                    t = cs * U[i][j] + sn * U[i][k - 1];
                                    U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                                    U[i][j] = t;
                                }
                            }
                        }
                        break;
                    // Perform one qr step.
                    case 3:
                        {
                            // Calculate the shift.
                            double maxPm1Pm2 = FastMath.max(FastMath.abs(singularValues[p - 1]),
                                                                  FastMath.abs(singularValues[p - 2]));
                            double scale = FastMath.max(FastMath.max(FastMath.max(maxPm1Pm2,
                                                                                        FastMath.abs(e[p - 2])),
                                                                           FastMath.abs(singularValues[k])),
                                                              FastMath.abs(e[k]));
                            double sp = singularValues[p - 1] / scale;
                            double spm1 = singularValues[p - 2] / scale;
                            double epm1 = e[p - 2] / scale;
                            double sk = singularValues[k] / scale;
                            double ek = e[k] / scale;
                            double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                            double c = (sp * epm1) * (sp * epm1);
                            double shift = 0;
                            if (b != 0 ||
                                c != 0)
                            {
                                shift = FastMath.sqrt(b * b + c);
                                if (b < 0)
                                {
                                    shift = -shift;
                                }
                                shift = c / (b + shift);
                            }
                            double f = (sk + sp) * (sk - sp) + shift;
                            double g = sk * ek;
                            // Chase zeros.
                            for (int j = k; j < p - 1; j++)
                            {
                                double t = FastMath.hypot(f, g);
                                double cs = f / t;
                                double sn = g / t;
                                if (j != k)
                                {
                                    e[j - 1] = t;
                                }
                                f = cs * singularValues[j] + sn * e[j];
                                e[j] = cs * e[j] - sn * singularValues[j];
                                g = sn * singularValues[j + 1];
                                singularValues[j + 1] = cs * singularValues[j + 1];

                                for (int i = 0; i < n; i++)
                                {
                                    t = cs * V[i][j] + sn * V[i][j + 1];
                                    V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                                    V[i][j] = t;
                                }
                                t = FastMath.hypot(f, g);
                                cs = f / t;
                                sn = g / t;
                                singularValues[j] = t;
                                f = cs * e[j] + sn * singularValues[j + 1];
                                singularValues[j + 1] = -sn * e[j] + cs * singularValues[j + 1];
                                g = sn * e[j + 1];
                                e[j + 1] = cs * e[j + 1];
                                if (j < m - 1)
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = cs * U[i][j] + sn * U[i][j + 1];
                                        U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                                        U[i][j] = t;
                                    }
                                }
                            }
                            e[p - 2] = f;
                        }
                        break;
                    // Convergence.
                    default:
                        {
                            // Make the singular values positive.
                            if (singularValues[k] <= 0)
                            {
                                singularValues[k] = singularValues[k] < 0 ? -singularValues[k] : 0;

                                for (int i = 0; i <= pp; i++)
                                {
                                    V[i][k] = -V[i][k];
                                }
                            }
                            // Order the singular values.
                            while (k < pp)
                            {
                                if (singularValues[k] >= singularValues[k + 1])
                                {
                                    break;
                                }
                                double t = singularValues[k];
                                singularValues[k] = singularValues[k + 1];
                                singularValues[k + 1] = t;
                                if (k < n - 1)
                                {
                                    for (int i = 0; i < n; i++)
                                    {
                                        t = V[i][k + 1];
                                        V[i][k + 1] = V[i][k];
                                        V[i][k] = t;
                                    }
                                }
                                if (k < m - 1)
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = U[i][k + 1];
                                        U[i][k + 1] = U[i][k];
                                        U[i][k] = t;
                                    }
                                }
                                k++;
                            }
                            p--;
                        }
                        break;
                }
            }

            // Set the small value tolerance used to calculate rank and pseudo-inverse
            tol = FastMath.max(m * singularValues[0] * EPS,
                               FastMath.sqrt(Precision.SAFE_MIN));

            if (!transposed)
            {
                cachedU = MatrixUtils.createRealMatrix(U);
                cachedV = MatrixUtils.createRealMatrix(V);
            }
            else
            {
                cachedU = MatrixUtils.createRealMatrix(V);
                cachedV = MatrixUtils.createRealMatrix(U);
            }
        }

        /// <summary>
        /// Returns the matrix U of the decomposition.
        /// <para>U is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the U matrix</returns>
        /// <remarks>
        /// See <see cref="getUT()"/>
        /// </remarks>
        public RealMatrix getU()
        {
            // return the cached matrix
            return cachedU;

        }

        /// <summary>
        /// Returns the transpose of the matrix U of the decomposition.
        /// <para>U is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the U matrix (or null if decomposed matrix is singular)</returns>
        /// <remarks>
        /// See <see cref="getU()"/>
        /// </remarks>
        public RealMatrix getUT()
        {
            if (cachedUt == null)
            {
                cachedUt = getU().transpose();
            }
            // return the cached matrix
            return cachedUt;
        }

        /// <summary>
        /// Returns the diagonal matrix &Sigma; of the decomposition.
        /// <para>&Sigma; is a diagonal matrix. The singular values are provided in
        /// non-increasing order, for compatibility with Jama.</para>
        /// </summary>
        /// <returns>the &Sigma; matrix</returns>
        public RealMatrix getS()
        {
            if (cachedS == null)
            {
                // cache the matrix for subsequent calls
                cachedS = MatrixUtils.createRealDiagonalMatrix(singularValues);
            }
            return cachedS;
        }

        /// <summary>
        /// Returns the diagonal elements of the matrix &Sigma; of the decomposition.
        /// <para>The singular values are provided in non-increasing order, for
        /// compatibility with Jama.</para>
        /// </summary>
        /// <returns>the diagonal elements of the &Sigma; matrix</returns>
        public double[] getSingularValues()
        {
            return (Double[])singularValues.Clone();
        }

        /// <summary>
        /// Returns the matrix V of the decomposition.
        /// <para>V is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the V matrix (or null if decomposed matrix is singular)</returns>
        /// <remarks>
        /// See <see cref="getVT()"/>
        /// </remarks>
        public RealMatrix getV()
        {
            // return the cached matrix
            return cachedV;
        }

        /// <summary>
        /// Returns the transpose of the matrix V of the decomposition.
        /// <para>V is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the V matrix (or null if decomposed matrix is singular)</returns>
        /// <remarks>
        /// See <see cref="getV()"/>
        /// </remarks>
        public RealMatrix getVT()
        {
            if (cachedVt == null)
            {
                cachedVt = getV().transpose();
            }
            // return the cached matrix
            return cachedVt;
        }

        /// <summary>
        /// Returns the n &times; n covariance matrix.
        /// <para>The covariance matrix is V &times; J &times; V^T
        /// where J is the diagonal matrix of the inverse of the squares of
        /// the singular values.</para>
        /// </summary>
        /// <param name="minSingularValue">value below which singular values are ignored
        /// (a 0 or negative value implies all singular value will be used)</param>
        /// <returns>covariance matrix</returns>
        /// <exception cref="IllegalArgumentException"> if minSingularValue is larger than
        /// the largest singular value, meaning all singular values are ignored</exception>
        public RealMatrix getCovariance(double minSingularValue)
        {
            // get the number of singular values to consider
            int p = singularValues.Length;
            int dimension = 0;
            while (dimension < p &&
                   singularValues[dimension] >= minSingularValue)
            {
                ++dimension;
            }

            if (dimension == 0)
            {
                throw new NumberIsTooLargeException<Double, Double>(new LocalizedFormats("TOO_LARGE_CUTOFF_SINGULAR_VALUE"), minSingularValue, singularValues[0], true);
            }

            double[][] data = new double[dimension][];
            getVT().walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitorAnonymous(data, singularValues), 0, dimension - 1, 0, p - 1);

            RealMatrix jv = new Array2DRowRealMatrix(data, false);
            return jv.transpose().multiply(jv);
        }

        private class DefaultRealMatrixPreservingVisitorAnonymous : DefaultRealMatrixPreservingVisitor
        {
            private double[][] data;
            private readonly double[] singularValues;

            public DefaultRealMatrixPreservingVisitorAnonymous(double[][] data, double[] singularValues)
            {
                this.data = data;
                this.singularValues = singularValues;
            }

            /// <inheritdoc/>
            public new void visit(int row, int column, double value)
            {
                data[row][column] = value / singularValues[row];
            }
        }

        /// <summary>
        /// Returns the L_2 norm of the matrix.
        /// <para>The L_2 norm is max(|A &times; u|_2 /
        /// |u|_2), where |.|_2 denotes the vectorial 2-norm
        /// (i.e. the traditional euclidian norm).</para>
        /// </summary>
        /// <returns>norm</returns>
        public double getNorm()
        {
            return singularValues[0];
        }

        /// <summary>
        /// Return the condition number of the matrix.
        /// </summary>
        /// <returns>condition number of the matrix</returns>
        public double getConditionNumber()
        {
            return singularValues[0] / singularValues[n - 1];
        }

        /// <summary>
        /// Computes the inverse of the condition number.
        /// In cases of rank deficiency, the <see cref="#getConditionNumber()">
        /// number</see> will become undefined.
        /// </summary>
        /// <returns>the inverse of the condition number.</returns>
        public double getInverseConditionNumber()
        {
            return singularValues[n - 1] / singularValues[0];
        }

        /// <summary>
        /// Return the effective numerical matrix rank.
        /// <para>The effective numerical rank is the number of non-negligible
        /// singular values. The threshold used to identify non-negligible
        /// terms is max(m,n) &times; ulp(s_1) where ulp(s_1)
        /// is the least significant bit of the largest singular value.</para>
        /// </summary>
        /// <returns>effective numerical matrix rank</returns>
        public int getRank()
        {
            int r = 0;
            for (int i = 0; i < singularValues.Length; i++)
            {
                if (singularValues[i] > tol)
                {
                    r++;
                }
            }
            return r;
        }

        /// <summary>
        /// Get a solver for finding the A &times; X = B solution in least square sense.
        /// </summary>
        /// <returns>a solver</returns>
        public DecompositionSolver getSolver()
        {
            return new Solver(singularValues, getUT(), getV(), getRank() == m, tol);
        }

        /// <summary>
        /// Specialized solver.
        /// </summary>
        private class Solver : DecompositionSolver
        {
            /// <summary>
            /// Pseudo-inverse of the initial matrix.
            /// </summary>
            private readonly RealMatrix pseudoInverse;
            
            /// <summary>
            /// Singularity indicator.
            /// </summary>
            private Boolean nonSingular;

            /// <summary>
            /// Build a solver from decomposed matrix.
            /// </summary>
            /// <param name="singularValues">Singular values.</param>
            /// <param name="uT">U^T matrix of the decomposition.</param>
            /// <param name="v">V matrix of the decomposition.</param>
            /// <param name="nonSingular">Singularity indicator.</param>
            /// <param name="tol">tolerance for singular values</param>
            public Solver(double[] singularValues, RealMatrix uT, RealMatrix v, Boolean nonSingular, double tol)
            {
                double[][] suT = uT.getData();
                for (int i = 0; i < singularValues.Length; ++i)
                {
                    double a;
                    if (singularValues[i] > tol)
                    {
                        a = 1 / singularValues[i];
                    }
                    else
                    {
                        a = 0;
                    }
                    double[] suTi = suT[i];
                    for (int j = 0; j < suTi.Length; ++j)
                    {
                        suTi[j] *= a;
                    }
                }
                pseudoInverse = v.multiply(new Array2DRowRealMatrix(suT, false));
                this.nonSingular = nonSingular;
            }

            /// <summary>
            /// Solve the linear equation A &times; X = B in least square sense.
            /// <para>
            /// The m&times;n matrix A may not be square, the solution X is such that
            /// ||A &times; X - B|| is minimal.
            /// </para>
            /// </summary>
            /// <param name="b">Right-hand side of the equation A &times; X = B</param>
            /// <returns>a vector X that minimizes the two norm of A &times; X - B</returns>
            /// <exception cref="DimensionMismatchException">
            /// if the matrices dimensions do not match.</exception>
            public RealVector solve(RealVector b)
            {
                return pseudoInverse.operate(b);
            }

            /// <summary>
            /// Solve the linear equation A &times; X = B in least square sense.
            /// <para>
            /// The m&times;n matrix A may not be square, the solution X is such that
            /// ||A &times; X - B|| is minimal.
            /// </para>
            /// </summary>
            /// <param name="b">Right-hand side of the equation A &times; X = B</param>
            /// <returns>a matrix X that minimizes the two norm of A &times; X - B</returns>
            /// <exception cref="DimensionMismatchException">
            /// if the matrices dimensions do not match.</exception>
            public RealMatrix solve(RealMatrix b)
            {
                return pseudoInverse.multiply(b);
            }

            /// <summary>
            /// Check if the decomposed matrix is non-singular.
            /// </summary>
            /// <returns><c>true</c> if the decomposed matrix is non-singular.</returns>
            public Boolean isNonSingular()
            {
                return nonSingular;
            }

            /// <summary>
            /// Get the pseudo-inverse of the decomposed matrix.
            /// </summary>
            /// <returns>the inverse matrix.</returns>
            public RealMatrix getInverse()
            {
                return pseudoInverse;
            }
        }
    }
}