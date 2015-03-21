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
    /// Class transforming a general real matrix to Hessenberg form.
    /// <para>A m &times; m matrix A can be written as the product of three matrices: A = P
    /// &times; H &times; P^T with P an orthogonal matrix and H a Hessenberg
    /// matrix. Both P and H are m &times; m matrices.</para>
    /// <para>Transformation to Hessenberg form is often not a goal by itself, but it is an
    /// intermediate step in more general decomposition algorithms like
    /// <see cref="EigenDecomposition">eigen decomposition</see>. This class is therefore
    /// intended for internal use by the library and is not public. As a consequence
    /// of this explicitly limited scope, many methods directly returns references to
    /// internal arrays, not copies.</para>
    /// <para>This class is based on the method orthes in class EigenvalueDecomposition
    /// from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</para>
    /// </summary>
    /// <remarks>
    /// See <a href="http://mathworld.wolfram.com/HessenbergDecomposition.html">MathWorld</a>
    /// <para/>
    /// See <a href="http://en.wikipedia.org/wiki/Householder_transformation">Householder 
    /// Transformations</a>
    /// </remarks>
    internal class HessenbergTransformer
    {
        /// <summary>
        /// Householder vectors.
        /// </summary>
        private readonly double[][] householderVectors;
        
        /// <summary>
        /// Temporary storage vector.
        /// </summary>
        private readonly double[] ort;
        
        /// <summary>
        /// Cached value of P.
        /// </summary>
        private RealMatrix cachedP;
        
        /// <summary>
        /// Cached value of Pt.
        /// </summary>
        private RealMatrix cachedPt;
        
        /// <summary>
        /// Cached value of H.
        /// </summary>
        private RealMatrix cachedH;

        /// <summary>
        /// Build the transformation to Hessenberg form of a general matrix.
        /// </summary>
        /// <param name="matrix">matrix to transform</param>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square</exception>
        public HessenbergTransformer(RealMatrix matrix)
        {
            if (!matrix.isSquare())
            {
                throw new NonSquareMatrixException(matrix.getRowDimension(),
                        matrix.getColumnDimension());
            }

            int m = matrix.getRowDimension();
            householderVectors = matrix.getData();
            ort = new double[m];
            cachedP = null;
            cachedPt = null;
            cachedH = null;

            // transform matrix
            transform();
        }

        /// <summary>
        /// Returns the matrix P of the transform.
        /// <para>P is an orthogonal matrix, i.e. its inverse is also its transpose.</para>
        /// </summary>
        /// <returns>the P matrix</returns>
        public RealMatrix getP()
        {
            if (cachedP == null)
            {
                int n = householderVectors.Length;
                int high = n - 1;
                double[][] pa = new double[n][];

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        pa[i][j] = (i == j) ? 1 : 0;
                    }
                }

                for (int m = high - 1; m >= 1; m--)
                {
                    if (householderVectors[m][m - 1] != 0.0)
                    {
                        for (int i = m + 1; i <= high; i++)
                        {
                            ort[i] = householderVectors[i][m - 1];
                        }

                        for (int j = m; j <= high; j++)
                        {
                            double g = 0.0;

                            for (int i = m; i <= high; i++)
                            {
                                g += ort[i] * pa[i][j];
                            }

                            // Double division avoids possible underflow
                            g = (g / ort[m]) / householderVectors[m][m - 1];

                            for (int i = m; i <= high; i++)
                            {
                                pa[i][j] += g * ort[i];
                            }
                        }
                    }
                }

                cachedP = MatrixUtils.createRealMatrix(pa);
            }
            return cachedP;
        }

        /// <summary>
        /// Returns the transpose of the matrix P of the transform.
        /// <para>P is an orthogonal matrix, i.e. its inverse is also its transpose.</para>
        /// </summary>
        /// <returns>the transpose of the P matrix</returns>
        public RealMatrix getPT()
        {
            if (cachedPt == null)
            {
                cachedPt = getP().transpose();
            }

            // return the cached matrix
            return cachedPt;
        }

        /// <summary>
        /// Returns the Hessenberg matrix H of the transform.
        /// </summary>
        /// <returns>the H matrix</returns>
        public RealMatrix getH()
        {
            if (cachedH == null)
            {
                int m = householderVectors.Length;
                double[][] h = new double[m][];
                for (int i = 0; i < m; ++i)
                {
                    if (i > 0)
                    {
                        // copy the entry of the lower sub-diagonal
                        h[i][i - 1] = householderVectors[i][i - 1];
                    }

                    // copy upper triangular part of the matrix
                    for (int j = i; j < m; ++j)
                    {
                        h[i][j] = householderVectors[i][j];
                    }
                }
                cachedH = MatrixUtils.createRealMatrix(h);
            }

            // return the cached matrix
            return cachedH;
        }

        /// <summary>
        /// Get the Householder vectors of the transform.
        /// <para>Note that since this class is only intended for internal use, it returns
        /// directly a reference to its internal arrays, not a copy.</para>
        /// </summary>
        /// <returns>the main diagonal elements of the B matrix</returns>
        double[][] getHouseholderVectorsRef()
        {
            return householderVectors;
        }

        /// <summary>
        /// Transform original matrix to Hessenberg form.
        /// <para>Transformation is done using Householder transforms.</para>
        /// </summary>
        private void transform()
        {
            int n = householderVectors.Length;
            int high = n - 1;

            for (int m = 1; m <= high - 1; m++)
            {
                // Scale column.
                double scale = 0;
                for (int i = m; i <= high; i++)
                {
                    scale += FastMath.abs(householderVectors[i][m - 1]);
                }

                if (!Precision.equals(scale, 0))
                {
                    // Compute Householder transformation.
                    double h = 0;
                    for (int i = high; i >= m; i--)
                    {
                        ort[i] = householderVectors[i][m - 1] / scale;
                        h += ort[i] * ort[i];
                    }
                    double g = (ort[m] > 0) ? -FastMath.sqrt(h) : FastMath.sqrt(h);

                    h -= ort[m] * g;
                    ort[m] -= g;

                    // Apply Householder similarity transformation
                    // H = (I - u*u' / h) * H * (I - u*u' / h)

                    for (int j = m; j < n; j++)
                    {
                        double f = 0;
                        for (int i = high; i >= m; i--)
                        {
                            f += ort[i] * householderVectors[i][j];
                        }
                        f /= h;
                        for (int i = m; i <= high; i++)
                        {
                            householderVectors[i][j] -= f * ort[i];
                        }
                    }

                    for (int i = 0; i <= high; i++)
                    {
                        double f = 0;
                        for (int j = high; j >= m; j--)
                        {
                            f += ort[j] * householderVectors[i][j];
                        }
                        f /= h;
                        for (int j = m; j <= high; j++)
                        {
                            householderVectors[i][j] -= f * ort[j];
                        }
                    }

                    ort[m] = scale * ort[m];
                    householderVectors[m][m - 1] = scale * g;
                }
            }
        }
    }
}