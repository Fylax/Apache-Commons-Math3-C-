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
    /// Implementation of a diagonal matrix.
    /// </summary>
    public class DiagonalMatrix : AbstractRealMatrix
    {
        /// <summary>
        /// Entries of the diagonal.
        /// </summary>
        private readonly double[] data;

        /// <summary>
        /// Creates a matrix with the supplied dimension.
        /// </summary>
        /// <param name="dimension">Number of rows and columns in the new matrix.</param>
        /// <exception cref="NotStrictlyPositiveException"> if the dimension is
        /// not positive.</exception>
        public DiagonalMatrix(int dimension)
            : base(dimension, dimension)
        {
            data = new double[dimension];
        }

        /// <summary>
        /// Creates a matrix using the input array as the underlying data.
        /// <para/>
        /// The input array is copied, not referenced.
        /// </summary>
        /// <param name="d">Data for the new matrix.</param>
        public DiagonalMatrix(double[] d) : this(d, true) { }

        /// <summary>
        /// Creates a matrix using the input array as the underlying data.
        /// <para/>
        /// If an array is created specially in order to be embedded in a
        /// this instance and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>.
        /// This will prevent the copying and improve performance as no new
        /// array will be built and no data will be copied.
        /// </summary>
        /// <param name="d">Data for new matrix.</param>
        /// <param name="copyArray">if <c>true</c>, the input array will be copied,
        /// otherwise it will be referenced.</param>
        public DiagonalMatrix(double[] d, Boolean copyArray)
        {
            MathUtils.checkNotNull(d);
            data = copyArray ? (double[])d.Clone() : d;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if the requested dimensions are not equal.</exception>
        public override RealMatrix createMatrix(int rowDimension, int columnDimension)
        {
            if (rowDimension != columnDimension)
            {
                throw new DimensionMismatchException(rowDimension, columnDimension);
            }

            return new DiagonalMatrix(rowDimension);
        }

        /// <inheritdoc/>
        public override RealMatrix copy()
        {
            return new DiagonalMatrix(data);
        }

        /// <summary>
        /// Compute the sum of <c>this</c> and <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to be added.</param>
        /// <returns><c>this + m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        public DiagonalMatrix add(DiagonalMatrix m)
        {
            // Safety check.
            MatrixUtils.checkAdditionCompatible(this, m);

            int dim = getRowDimension();
            double[] outData = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                outData[i] = data[i] + m.data[i];
            }

            return new DiagonalMatrix(outData, false);
        }

        /// <summary>
        /// Returns <c>this</c> minus <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to be subtracted.</param>
        /// <returns><c>this - m</c></returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        public DiagonalMatrix subtract(DiagonalMatrix m)
        {
            MatrixUtils.checkSubtractionCompatible(this, m);

            int dim = getRowDimension();
            double[] outData = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                outData[i] = data[i] - m.data[i];
            }

            return new DiagonalMatrix(outData, false);
        }

        /// <summary>
        /// Returns the result of postmultiplying <c>this</c> by <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to postmultiply by</param>
        /// <returns><c>this * m</c></returns>
        /// <exception cref="DimensionMismatchException"> if
        /// <c>columnDimension(this) != rowDimension(m)</c></exception>
        public DiagonalMatrix multiply(DiagonalMatrix m)
        {
            MatrixUtils.checkMultiplicationCompatible(this, m);

            int dim = getRowDimension();
            double[] outData = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                outData[i] = data[i] * m.data[i];
            }

            return new DiagonalMatrix(outData, false);
        }

        /// <summary>
        /// Returns the result of postmultiplying <c>this</c> by <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to postmultiply by</param>
        /// <returns><c>this * m</c></returns>
        /// <exception cref="DimensionMismatchException"> if
        /// <c>columnDimension(this) != rowDimension(m)</c></exception>
        public new RealMatrix multiply(RealMatrix m)
        {
            if (m is DiagonalMatrix)
            {
                return multiply((DiagonalMatrix)m);
            }
            else
            {
                MatrixUtils.checkMultiplicationCompatible(this, m);
                int nRows = m.getRowDimension();
                int nCols = m.getColumnDimension();
                double[][] product = new double[nRows][];
                for (int r = 0; r < nRows; r++)
                {
                    for (int c = 0; c < nCols; c++)
                    {
                        product[r][c] = data[r] * m.getEntry(r, c);
                    }
                }
                return new Array2DRowRealMatrix(product, false);
            }
        }

        /// <inheritdoc/>
        public new double[][] getData()
        {
            int dim = getRowDimension();
            double[][] outp = new double[dim][];

            for (int i = 0; i < dim; i++)
            {
                outp[i][i] = data[i];
            }

            return outp;
        }

        /// <summary>
        /// Gets a reference to the underlying data array.
        /// </summary>
        /// <returns>1-dimensional array of entries.</returns>
        public double[] getDataRef()
        {
            return data;
        }

        /// <inheritdoc/>
        public override double getEntry(int row, int column)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            return row == column ? data[row] : 0;
        }

        /// <inheritdoc/>
        /// <exception cref="NumberIsTooLargeException"> if <c>row != column</c> and value is
        /// non-zero.</exception>
        public override void setEntry(int row, int column, double value)
        {
            if (row == column)
            {
                MatrixUtils.checkRowIndex(this, row);
                data[row] = value;
            }
            else
            {
                ensureZero(value);
            }
        }

        /// <inheritdoc/>
        /// <exception cref="NumberIsTooLargeException"> if <c>row != column</c> and increment is
        /// non-zero.</exception>
        public new void addToEntry(int row, int column, double increment)
        {
            if (row == column)
            {
                MatrixUtils.checkRowIndex(this, row);
                data[row] += increment;
            }
            else
            {
                ensureZero(increment);
            }
        }

        /// <inheritdoc/>
        public new void multiplyEntry(int row, int column, double factor)
        {
            // we don't care about non-diagonal elements for multiplication
            if (row == column)
            {
                MatrixUtils.checkRowIndex(this, row);
                data[row] *= factor;
            }
        }

        /// <inheritdoc/>
        public override int getRowDimension()
        {
            return data.Length;
        }

        /// <inheritdoc/>
        public override int getColumnDimension()
        {
            return data.Length;
        }

        /// <inheritdoc/>
        public new double[] operate(double[] v)
        {
            return multiply(new DiagonalMatrix(v, false)).getDataRef();
        }

        /// <inheritdoc/>
        public new double[] preMultiply(double[] v)
        {
            return operate(v);
        }

        /// <inheritdoc/>
        public new RealVector preMultiply(RealVector v)
        {
            double[] vectorData;
            if (v is ArrayRealVector)
            {
                vectorData = ((ArrayRealVector)v).getDataRef();
            }
            else
            {
                vectorData = v.toArray();
            }
            return MatrixUtils.createRealVector(preMultiply(vectorData));
        }

        /// <summary>
        /// Ensure a value is zero.
        /// </summary>
        /// <param name="value">value to check</param>
        /// <exception cref="NumberIsTooLargeException"> if value is not zero</exception>
        private void ensureZero(double value)
        {
            if (!Precision.equals(0.0, value, 1))
            {
                throw new NumberIsTooLargeException<Double, Int32>(FastMath.abs(value), 0, true);
            }
        }

        /// <summary>
        /// Computes the inverse of this diagonal matrix.
        /// <para>
        /// Note: this method will use a singularity threshold of 0,
        /// use <see cref="inverse(double)"/> if a different threshold is needed.
        /// </para>
        /// </summary>
        /// <returns>the inverse of <c>m</c></returns>
        /// <exception cref="SingularMatrixException"> if the matrix is singular</exception>
        public DiagonalMatrix inverse()
        {
            return inverse(0);
        }

        /// <summary>
        /// Computes the inverse of this diagonal matrix.
        /// </summary>
        /// <param name="threshold">Singularity threshold.</param>
        /// <returns>the inverse of <c>m</c></returns>
        /// <exception cref="SingularMatrixException"> if the matrix is singular</exception>
        public DiagonalMatrix inverse(double threshold)
        {
            if (isSingular(threshold))
            {
                throw new SingularMatrixException();
            }

            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; i++)
            {
                result[i] = 1.0 / data[i];
            }
            return new DiagonalMatrix(result, false);
        }

        /// <summary>
        /// Returns whether this diagonal matrix is singular, i.e. any diagonal entry
        /// is equal to <c>0</c> within the given threshold.
        /// </summary>
        /// <param name="threshold">Singularity threshold.</param>
        /// <returns><c>true</c> if the matrix is singular, <c>false</c> otherwise</returns>
        public Boolean isSingular(double threshold)
        {
            for (int i = 0; i < data.Length; i++)
            {
                if (Precision.equals(data[i], 0.0, threshold))
                {
                    return true;
                }
            }
            return false;
        }
    }
}