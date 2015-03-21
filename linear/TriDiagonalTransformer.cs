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

namespace Math3.linear
{
    /// <summary>
    /// Class transforming a symmetrical matrix to tridiagonal shape.
    /// <para>A symmetrical m &times; m matrix A can be written as the product of three matrices:
    /// A = Q &times; T &times; Q^T with Q an orthogonal matrix and T a symmetrical
    /// tridiagonal matrix. Both Q and T are m &times; m matrices.</para>
    /// <para>This implementation only uses the upper part of the matrix, the part below the
    /// diagonal is not accessed at all.</para>
    /// <para>Transformation to tridiagonal shape is often not a goal by itself, but it is
    /// an intermediate step in more general decomposition algorithms like 
    /// <see cref="EigenDecomposition">eigen decomposition</see>. This class is therefore 
    /// intended for internal use by the library and is not public. As a consequence of this 
    /// explicitly limited scope, many methods directly returns references to internal arrays,
    /// not copies.</p>
    /// </summary>
    internal class TriDiagonalTransformer
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
        /// Cached value of Q.
        /// </summary>
        private RealMatrix cachedQ;
        
        /// <summary>
        /// Cached value of Qt.
        /// </summary>
        private RealMatrix cachedQt;
        
        /// <summary>
        /// Cached value of T.
        /// </summary>
        private RealMatrix cachedT;

        /// <summary>
        /// Build the transformation to tridiagonal shape of a symmetrical matrix.
        /// <para>The specified matrix is assumed to be symmetrical without any check.
        /// Only the upper triangular part of the matrix is used.</para>
        /// </summary>
        /// <param name="matrix">Symmetrical matrix to transform.</param>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        public TriDiagonalTransformer(RealMatrix matrix)
        {
            if (!matrix.isSquare())
            {
                throw new NonSquareMatrixException(matrix.getRowDimension(),
                                                   matrix.getColumnDimension());
            }

            int m = matrix.getRowDimension();
            householderVectors = matrix.getData();
            main = new double[m];
            secondary = new double[m - 1];
            cachedQ = null;
            cachedQt = null;
            cachedT = null;

            // transform matrix
            transform();
        }

        /// <summary>
        /// Returns the matrix Q of the transform.
        /// <para>Q is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
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
        /// Returns the transpose of the matrix Q of the transform.
        /// <para>Q is an orthogonal matrix, i.e. its transpose is also its inverse.</para>
        /// </summary>
        /// <returns>the Q matrix</returns>
        public RealMatrix getQT()
        {
            if (cachedQt == null)
            {
                int m = householderVectors.Length;
                double[][] qta = new double[m][];

                // build up first part of the matrix by applying Householder transforms
                for (int k = m - 1; k >= 1; --k)
                {
                    double[] hK = householderVectors[k - 1];
                    qta[k][k] = 1;
                    if (hK[k] != 0.0)
                    {
                        double inv = 1.0 / (secondary[k - 1] * hK[k]);
                        double beta = 1.0 / secondary[k - 1];
                        qta[k][k] = 1 + beta * hK[k];
                        for (int i = k + 1; i < m; ++i)
                        {
                            qta[k][i] = beta * hK[i];
                        }
                        for (int j = k + 1; j < m; ++j)
                        {
                            beta = 0;
                            for (int i = k + 1; i < m; ++i)
                            {
                                beta += qta[j][i] * hK[i];
                            }
                            beta *= inv;
                            qta[j][k] = beta * hK[k];
                            for (int i = k + 1; i < m; ++i)
                            {
                                qta[j][i] += beta * hK[i];
                            }
                        }
                    }
                }
                qta[0][0] = 1;
                cachedQt = MatrixUtils.createRealMatrix(qta);
            }

            // return the cached matrix
            return cachedQt;
        }

        /// <summary>
        /// Returns the tridiagonal matrix T of the transform.
        /// </summary>
        /// <returns>the T matrix</returns>
        public RealMatrix getT()
        {
            if (cachedT == null)
            {
                int m = main.Length;
                double[][] ta = new double[m][];
                for (int i = 0; i < m; ++i)
                {
                    ta[i][i] = main[i];
                    if (i > 0)
                    {
                        ta[i][i - 1] = secondary[i - 1];
                    }
                    if (i < main.Length - 1)
                    {
                        ta[i][i + 1] = secondary[i];
                    }
                }
                cachedT = MatrixUtils.createRealMatrix(ta);
            }

            // return the cached matrix
            return cachedT;
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
        /// Get the main diagonal elements of the matrix T of the transform.
        /// <para>Note that since this class is only intended for internal use,
        /// it returns directly a reference to its internal arrays, not a copy.</para>
        /// </summary>
        /// <returns>the main diagonal elements of the T matrix</returns>
        internal double[] getMainDiagonalRef()
        {
            return main;
        }

        /// <summary>
        /// Get the secondary diagonal elements of the matrix T of the transform.
        /// <para>Note that since this class is only intended for internal use,
        /// it returns directly a reference to its internal arrays, not a copy.</para>
        /// </summary>
        /// <returns>the secondary diagonal elements of the T matrix</returns>
        internal double[] getSecondaryDiagonalRef()
        {
            return secondary;
        }

        /// <summary>
        /// Transform original matrix to tridiagonal form.
        /// <para>Transformation is done using Householder transforms.</para>
        /// </summary>
        private void transform()
        {
            int m = householderVectors.Length;
            double[] z = new double[m];
            for (int k = 0; k < m - 1; k++)
            {

                //zero-out a row and a column simultaneously
                double[] hK = householderVectors[k];
                main[k] = hK[k];
                double xNormSqr = 0;
                for (int j = k + 1; j < m; ++j)
                {
                    double c = hK[j];
                    xNormSqr += c * c;
                }
                double a = (hK[k + 1] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
                secondary[k] = a;
                if (a != 0.0)
                {
                    // apply Householder transform from left and right simultaneously

                    hK[k + 1] -= a;
                    double beta = -1 / (a * hK[k + 1]);

                    // compute a = beta A v, where v is the Householder vector
                    // this loop is written in such a way
                    //   1) only the upper triangular part of the matrix is accessed
                    //   2) access is cache-friendly for a matrix stored in rows
                    for (int i = k + 1; i < m; ++i)
                    {
                        z[i] = 0;
                    }
                    for (int i = k + 1; i < m; ++i)
                    {
                        double[] hI = householderVectors[i];
                        double hKI = hK[i];
                        double zI = hI[i] * hKI;
                        for (int j = i + 1; j < m; ++j)
                        {
                            double hIJ = hI[j];
                            zI += hIJ * hK[j];
                            z[j] += hIJ * hKI;
                        }
                        z[i] = beta * (z[i] + zI);
                    }

                    // compute gamma = beta vT z / 2
                    double gamma = 0;
                    for (int i = k + 1; i < m; ++i)
                    {
                        gamma += z[i] * hK[i];
                    }
                    gamma *= beta / 2;

                    // compute z = z - gamma v
                    for (int i = k + 1; i < m; ++i)
                    {
                        z[i] -= gamma * hK[i];
                    }

                    // update matrix: A = A - v zT - z vT
                    // only the upper triangular part of the matrix is updated
                    for (int i = k + 1; i < m; ++i)
                    {
                        double[] hI = householderVectors[i];
                        for (int j = i; j < m; ++j)
                        {
                            hI[j] -= hK[i] * z[j] + z[i] * hK[j];
                        }
                    }
                }
            }
            main[m - 1] = householderVectors[m - 1][m - 1];
        }
    }
}