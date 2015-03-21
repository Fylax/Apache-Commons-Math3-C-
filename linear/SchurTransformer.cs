using Math3.exception;
using Math3.exception.util;
using Math3.util;
using System;

namespace Math3.linear
{
    /// <summary>
    /// Class transforming a general real matrix to Schur form.
    /// <para>A m &times; m matrix A can be written as the product of three matrices: A = P
    /// &times; T &times; P^T with P an orthogonal matrix and T an quasi-triangular
    /// matrix. Both P and T are m &times; m matrices.</para>
    /// <para>Transformation to Schur form is often not a goal by itself, but it is an
    /// intermediate step in more general decomposition algorithms like
    /// <see cref="EigenDecomposition">eigen decomposition</see>. This class is therefore
    /// intended for internal use by the library and is not public. As a consequence
    /// of this explicitly limited scope, many methods directly returns references to
    /// internal arrays, not copies.</para>
    /// <para>This class is based on the method hqr2 in class EigenvalueDecomposition
    /// from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</para>
    /// </summary>
    /// <remarks>
    /// See <a href="http://mathworld.wolfram.com/SchurDecomposition.html">Schur Decomposition -
    /// MathWorld</a><para/>
    /// See <a href="http://en.wikipedia.org/wiki/Schur_decomposition">Schur Decomposition - 
    /// Wikipedia</a><para/>
    /// See <a href="http://en.wikipedia.org/wiki/Householder_transformation">Householder 
    /// Transformations</a>
    /// </remarks>
    internal class SchurTransformer
    {
        /// <summary>
        /// Maximum allowed iterations for convergence of the transformation.
        /// </summary>
        private const int MAX_ITERATIONS = 100;

        /// <summary>
        /// P matrix.
        /// </summary>
        private readonly double[][] matrixP;
        
        /// <summary>
        /// T matrix.
        /// </summary>
        private readonly double[][] matrixT;
        
        /// <summary>
        /// Cached value of P.
        /// </summary>
        private RealMatrix cachedP;
        
        /// <summary>
        /// Cached value of T.
        /// </summary>
        private RealMatrix cachedT;
        
        /// <summary>
        /// Cached value of PT.
        /// </summary>
        private RealMatrix cachedPt;

        /// <summary>
        /// Epsilon criteria taken from JAMA code (originally was 2^-52).
        /// </summary>
        private readonly double epsilon = Precision.EPSILON;

        /// <summary>
        /// Build the transformation to Schur form of a general real matrix.
        /// </summary>
        /// <param name="matrix">matrix to transform</param>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square</exception>
        public SchurTransformer(RealMatrix matrix)
        {
            if (!matrix.isSquare())
            {
                throw new NonSquareMatrixException(matrix.getRowDimension(),
                                                   matrix.getColumnDimension());
            }

            HessenbergTransformer transformer = new HessenbergTransformer(matrix);
            matrixT = transformer.getH().getData();
            matrixP = transformer.getP().getData();
            cachedT = null;
            cachedP = null;
            cachedPt = null;

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
                cachedP = MatrixUtils.createRealMatrix(matrixP);
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
        /// Returns the quasi-triangular Schur matrix T of the transform.
        /// </summary>
        /// <returns>the T matrix</returns>
        public RealMatrix getT()
        {
            if (cachedT == null)
            {
                cachedT = MatrixUtils.createRealMatrix(matrixT);
            }

            // return the cached matrix
            return cachedT;
        }

        /// <summary>
        /// Transform original matrix to Schur form.
        /// </summary>
        /// <exception cref="MaxCountExceededException"> if the transformation does not converge
        /// </exception>
        private void transform()
        {
            int n = matrixT.Length;

            // compute matrix norm
            double norm = getNorm();

            // shift information
            ShiftInfo shift = new ShiftInfo();

            // Outer loop over eigenvalue index
            int iteration = 0;
            int iu = n - 1;
            while (iu >= 0)
            {

                // Look for single small sub-diagonal element
                int il = findSmallSubDiagonalElement(iu, norm);

                // Check for convergence
                if (il == iu)
                {
                    // One root found
                    matrixT[iu][iu] += shift.exShift;
                    iu--;
                    iteration = 0;
                }
                else if (il == iu - 1)
                {
                    // Two roots found
                    double p = (matrixT[iu - 1][iu - 1] - matrixT[iu][iu]) / 2.0;
                    double q = p * p + matrixT[iu][iu - 1] * matrixT[iu - 1][iu];
                    matrixT[iu][iu] += shift.exShift;
                    matrixT[iu - 1][iu - 1] += shift.exShift;

                    if (q >= 0)
                    {
                        double z = FastMath.sqrt(FastMath.abs(q));
                        if (p >= 0)
                        {
                            z = p + z;
                        }
                        else
                        {
                            z = p - z;
                        }
                        double x = matrixT[iu][iu - 1];
                        double s = FastMath.abs(x) + FastMath.abs(z);
                        p = x / s;
                        q = z / s;
                        double r = FastMath.sqrt(p * p + q * q);
                        p /= r;
                        q /= r;

                        // Row modification
                        for (int j = iu - 1; j < n; j++)
                        {
                            z = matrixT[iu - 1][j];
                            matrixT[iu - 1][j] = q * z + p * matrixT[iu][j];
                            matrixT[iu][j] = q * matrixT[iu][j] - p * z;
                        }

                        // Column modification
                        for (int i = 0; i <= iu; i++)
                        {
                            z = matrixT[i][iu - 1];
                            matrixT[i][iu - 1] = q * z + p * matrixT[i][iu];
                            matrixT[i][iu] = q * matrixT[i][iu] - p * z;
                        }

                        // Accumulate transformations
                        for (int i = 0; i <= n - 1; i++)
                        {
                            z = matrixP[i][iu - 1];
                            matrixP[i][iu - 1] = q * z + p * matrixP[i][iu];
                            matrixP[i][iu] = q * matrixP[i][iu] - p * z;
                        }
                    }
                    iu -= 2;
                    iteration = 0;
                }
                else
                {
                    // No convergence yet
                    computeShift(il, iu, iteration, shift);

                    // stop transformation after too many iterations
                    if (++iteration > MAX_ITERATIONS)
                    {
                        throw new MaxCountExceededException<Int32>(new LocalizedFormats("CONVERGENCE_FAILED"), MAX_ITERATIONS);
                    }

                    // the initial houseHolder vector for the QR step
                    double[] hVec = new double[3];

                    int im = initQRStep(il, iu, shift, hVec);
                    performDoubleQRStep(il, im, iu, shift, hVec);
                }
            }
        }

        /// <summary>
        /// Computes the L1 norm of the (quasi-)triangular matrix T.
        /// </summary>
        /// <returns>the L1 norm of matrix T</returns>
        private double getNorm()
        {
            double norm = 0.0;
            for (int i = 0; i < matrixT.Length; i++)
            {
                // as matrix T is (quasi-)triangular, also take the sub-diagonal element into account
                for (int j = FastMath.max(i - 1, 0); j < matrixT.Length; j++)
                {
                    norm += FastMath.abs(matrixT[i][j]);
                }
            }
            return norm;
        }

        /// <summary>
        /// Find the first small sub-diagonal element and returns its index.
        /// </summary>
        /// <param name="startIdx">the starting index for the search</param>
        /// <param name="norm">the L1 norm of the matrix</param>
        /// <returns>the index of the first small sub-diagonal element</returns>
        private int findSmallSubDiagonalElement(int startIdx, double norm)
        {
            int l = startIdx;
            while (l > 0)
            {
                double s = FastMath.abs(matrixT[l - 1][l - 1]) + FastMath.abs(matrixT[l][l]);
                if (s == 0.0)
                {
                    s = norm;
                }
                if (FastMath.abs(matrixT[l][l - 1]) < epsilon * s)
                {
                    break;
                }
                l--;
            }
            return l;
        }

        /// <summary>
        /// Compute the shift for the current iteration.
        /// </summary>
        /// <param name="l">the index of the small sub-diagonal element</param>
        /// <param name="idx">the current eigenvalue index</param>
        /// <param name="iteration">the current iteration</param>
        /// <param name="shift">holder for shift information</param>
        private void computeShift(int l, int idx, int iteration, ShiftInfo shift)
        {
            // Form shift
            shift.x = matrixT[idx][idx];
            shift.y = shift.w = 0.0;
            if (l < idx)
            {
                shift.y = matrixT[idx - 1][idx - 1];
                shift.w = matrixT[idx][idx - 1] * matrixT[idx - 1][idx];
            }

            // Wilkinson's original ad hoc shift
            if (iteration == 10)
            {
                shift.exShift += shift.x;
                for (int i = 0; i <= idx; i++)
                {
                    matrixT[i][i] -= shift.x;
                }
                double s = FastMath.abs(matrixT[idx][idx - 1]) + FastMath.abs(matrixT[idx - 1][idx - 2]);
                shift.x = 0.75 * s;
                shift.y = 0.75 * s;
                shift.w = -0.4375 * s * s;
            }

            // MATLAB's new ad hoc shift
            if (iteration == 30)
            {
                double s = (shift.y - shift.x) / 2.0;
                s = s * s + shift.w;
                if (s > 0.0)
                {
                    s = FastMath.sqrt(s);
                    if (shift.y < shift.x)
                    {
                        s = -s;
                    }
                    s = shift.x - shift.w / ((shift.y - shift.x) / 2.0 + s);
                    for (int i = 0; i <= idx; i++)
                    {
                        matrixT[i][i] -= s;
                    }
                    shift.exShift += s;
                    shift.x = shift.y = shift.w = 0.964;
                }
            }
        }

        /// <summary>
        /// Initialize the householder vectors for the QR step.
        /// </summary>
        /// <param name="il">the index of the small sub-diagonal element</param>
        /// <param name="iu">the current eigenvalue index</param>
        /// <param name="shift">shift information holder</param>
        /// <param name="hVec">the initial houseHolder vector</param>
        /// <returns>the start index for the QR step</returns>
        private int initQRStep(int il, int iu, ShiftInfo shift, double[] hVec)
        {
            // Look for two consecutive small sub-diagonal elements
            int im = iu - 2;
            while (im >= il)
            {
                double z = matrixT[im][im];
                double r = shift.x - z;
                double s = shift.y - z;
                hVec[0] = (r * s - shift.w) / matrixT[im + 1][im] + matrixT[im][im + 1];
                hVec[1] = matrixT[im + 1][im + 1] - z - r - s;
                hVec[2] = matrixT[im + 2][im + 1];

                if (im == il)
                {
                    break;
                }

                double lhs = FastMath.abs(matrixT[im][im - 1]) * (FastMath.abs(hVec[1]) + FastMath.abs(hVec[2]));
                double rhs = FastMath.abs(hVec[0]) * (FastMath.abs(matrixT[im - 1][im - 1]) +
                                                            FastMath.abs(z) +
                                                            FastMath.abs(matrixT[im + 1][im + 1]));

                if (lhs < epsilon * rhs)
                {
                    break;
                }
                im--;
            }

            return im;
        }

        /// <summary>
        /// Perform a double QR step involving rows l:idx and columns m:n
        /// </summary>
        /// <param name="il">the index of the small sub-diagonal element</param>
        /// <param name="im">the start index for the QR step</param>
        /// <param name="iu">the current eigenvalue index</param>
        /// <param name="shift">shift information holder</param>
        /// <param name="hVec">the initial houseHolder vector</param>
        private void performDoubleQRStep(int il, int im, int iu, ShiftInfo shift, double[] hVec)
        {
            int n = matrixT.Length;
            double p = hVec[0];
            double q = hVec[1];
            double r = hVec[2];

            for (int k = im; k <= iu - 1; k++)
            {
                Boolean notlast = k != (iu - 1);
                if (k != im)
                {
                    p = matrixT[k][k - 1];
                    q = matrixT[k + 1][k - 1];
                    r = notlast ? matrixT[k + 2][k - 1] : 0.0;
                    shift.x = FastMath.abs(p) + FastMath.abs(q) + FastMath.abs(r);
                    if (Precision.equals(shift.x, 0.0, epsilon))
                    {
                        continue;
                    }
                    p /= shift.x;
                    q /= shift.x;
                    r /= shift.x;
                }
                double s = FastMath.sqrt(p * p + q * q + r * r);
                if (p < 0.0)
                {
                    s = -s;
                }
                if (s != 0.0)
                {
                    if (k != im)
                    {
                        matrixT[k][k - 1] = -s * shift.x;
                    }
                    else if (il != im)
                    {
                        matrixT[k][k - 1] = -matrixT[k][k - 1];
                    }
                    p += s;
                    shift.x = p / s;
                    shift.y = q / s;
                    double z = r / s;
                    q /= p;
                    r /= p;

                    // Row modification
                    for (int j = k; j < n; j++)
                    {
                        p = matrixT[k][j] + q * matrixT[k + 1][j];
                        if (notlast)
                        {
                            p += r * matrixT[k + 2][j];
                            matrixT[k + 2][j] -= p * z;
                        }
                        matrixT[k][j] -= p * shift.x;
                        matrixT[k + 1][j] -= p * shift.y;
                    }

                    // Column modification
                    for (int i = 0; i <= FastMath.min(iu, k + 3); i++)
                    {
                        p = shift.x * matrixT[i][k] + shift.y * matrixT[i][k + 1];
                        if (notlast)
                        {
                            p += z * matrixT[i][k + 2];
                            matrixT[i][k + 2] -= p * r;
                        }
                        matrixT[i][k] -= p;
                        matrixT[i][k + 1] -= p * q;
                    }

                    // Accumulate transformations
                    int high = matrixT.Length - 1;
                    for (int i = 0; i <= high; i++)
                    {
                        p = shift.x * matrixP[i][k] + shift.y * matrixP[i][k + 1];
                        if (notlast)
                        {
                            p += z * matrixP[i][k + 2];
                            matrixP[i][k + 2] -= p * r;
                        }
                        matrixP[i][k] -= p;
                        matrixP[i][k + 1] -= p * q;
                    }
                }  // (s != 0)
            }  // k loop

            // clean up pollution due to round-off errors
            for (int i = im + 2; i <= iu; i++)
            {
                matrixT[i][i - 2] = 0.0;
                if (i > im + 2)
                {
                    matrixT[i][i - 3] = 0.0;
                }
            }
        }

        /// <summary>
        /// Internal data structure holding the current shift information.
        /// Contains variable names as present in the original JAMA code.
        /// </summary>
        private class ShiftInfo
        {
            // CHECKSTYLE: stop all

            /// <summary>
            /// x shift info
            /// </summary>
            internal double x;
            
            /// <summary>
            /// y shift info
            /// </summary>
            internal double y;
            
            /// <summary>
            /// w shift info
            /// </summary>
            internal double w;
            
            /// <summary>
            /// Indicates an exceptional shift.
            /// </summary>
            internal double exShift;

            // CHECKSTYLE: resume all
        }
    }
}