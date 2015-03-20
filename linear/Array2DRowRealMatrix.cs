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
    /// Implementation of <see cref="RealMatrix"/> using a <c>double[][]</c> array to
    /// store entries.
    /// </summary>
    public class Array2DRowRealMatrix : AbstractRealMatrix
    {
        /// <summary>
        /// Entries of the matrix.
        /// </summary>
        private double[][] data;

        /// <summary>
        /// Creates a matrix with no data
        /// </summary>
        public Array2DRowRealMatrix() { }

        /// <summary>
        /// Create a new RealMatrix with the supplied row and column dimensions.
        /// </summary>
        /// <param name="rowDimension">Number of rows in the new matrix.</param>
        /// <param name="columnDimension">Number of columns in the new matrix.</param>
        /// <exception cref="NotStrictlyPositiveException"> if the row or column dimension is
        /// not positive.</exception>
        public Array2DRowRealMatrix(int rowDimension, int columnDimension)
            : base(rowDimension, columnDimension)
        {
            data = new double[rowDimension][];
        }

        /// <summary>
        /// Create a new <c>RealMatrix</c> using the input array as the underlying
        /// data array.
        /// <para>The input array is copied, not referenced. This constructor has
        /// the same effect as calling <see cref="Array2DRowRealMatrix(double[][], boolean)"/>
        /// with the second argument set to <c>true</c>.</para>
        /// </summary>
        /// <param name="d">Data for the new matrix.</param>
        /// <exception cref="DimensionMismatchException"> if <c>d</c> is not rectangular.
        /// </exception>
        /// <exception cref="NoDataException"> if <c>d</c> row or column dimension is zero.
        /// </exception>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <remarks>
        /// See <see cref="Array2DRowRealMatrix(double[][], boolean)"/>
        /// </remarks>
        public Array2DRowRealMatrix(double[][] d)
        {
            copyIn(d);
        }

        /// <summary>
        /// Create a new RealMatrix using the input array as the underlying
        /// data array.
        /// If an array is built specially in order to be embedded in a
        /// RealMatrix and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>. This will prevent the copying and improve
        /// performance as no new array will be built and no data will be copied.
        /// </summary>
        /// <param name="d">Data for new matrix.</param>
        /// <param name="copyArray">if <c>true</c>, the input array will be copied,
        /// otherwise it will be referenced.</param>
        /// <exception cref="DimensionMismatchException"> if <c>d</c> is not rectangular.
        /// </exception>
        /// <exception cref="NoDataException"> if <c>d</c> row or column dimension is zero.
        /// </exception>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <remarks>
        /// See <see cref="Array2DRowRealMatrix(double[][])"/>
        /// </remarks>
        public Array2DRowRealMatrix(double[][] d, Boolean copyArray)
        {
            if (copyArray)
            {
                copyIn(d);
            }
            else
            {
                if (d == null)
                {
                    throw new NullArgumentException();
                }
                int nRows = d.Length;
                if (nRows == 0)
                {
                    throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_ROW"));
                }
                int nCols = d[0].Length;
                if (nCols == 0)
                {
                    throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_COLUMN"));
                }
                for (int r = 1; r < nRows; r++)
                {
                    if (d[r].Length != nCols)
                    {
                        throw new DimensionMismatchException(d[r].Length, nCols);
                    }
                }
                data = d;
            }
        }

        /// <summary>
        /// Create a new (column) RealMatrix using <c>v</c> as the
        /// data for the unique column of the created matrix.
        /// The input array is copied.
        /// </summary>
        /// <param name="v">Column vector holding data for new matrix.</param>
        public Array2DRowRealMatrix(double[] v)
        {
            int nRows = v.Length;
            data = new double[nRows][];
            for (int row = 0; row < nRows; row++)
            {
                data[row][0] = v[row];
            }
        }

        /// <inheritdoc/>
        public override RealMatrix createMatrix(int rowDimension, int columnDimension)
        {
            return new Array2DRowRealMatrix(rowDimension, columnDimension);
        }

        /// <inheritdoc/>
        public override RealMatrix copy()
        {
            return new Array2DRowRealMatrix(copyOut(), false);
        }

        /// <summary>
        /// Compute the sum of <c>this</c> and <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to be added.</param>
        /// <returns><c>this + m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        public Array2DRowRealMatrix add(Array2DRowRealMatrix m)
        {
            // Safety check.
            MatrixUtils.checkAdditionCompatible(this, m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            double[][] outData = new double[rowCount][];
            for (int row = 0; row < rowCount; row++)
            {
                double[] dataRow = data[row];
                double[] mRow = m.data[row];
                double[] outDataRow = outData[row];
                for (int col = 0; col < columnCount; col++)
                {
                    outDataRow[col] = dataRow[col] + mRow[col];
                }
            }

            return new Array2DRowRealMatrix(outData, false);
        }

        /// <summary>
        /// Returns <c>this</c> minus <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to be subtracted.</param>
        /// <returns><c>this - m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        public Array2DRowRealMatrix subtract(Array2DRowRealMatrix m)
        {
            MatrixUtils.checkSubtractionCompatible(this, m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            double[][] outData = new double[rowCount][];
            for (int row = 0; row < rowCount; row++)
            {
                double[] dataRow = data[row];
                double[] mRow = m.data[row];
                double[] outDataRow = outData[row];
                for (int col = 0; col < columnCount; col++)
                {
                    outDataRow[col] = dataRow[col] - mRow[col];
                }
            }

            return new Array2DRowRealMatrix(outData, false);
        }

        /// <summary>
        /// Returns the result of postmultiplying <c>this</c> by <c>m</c>. 
        /// </summary>
        /// <param name="m">matrix to postmultiply by</param>
        /// <returns><c>this * m</c></returns>
        /// <exception cref="DimensionMismatchException"> if
        /// <c>columnDimension(this) != rowDimension(m)</c></exception>
        public Array2DRowRealMatrix multiply(Array2DRowRealMatrix m)
        {
            MatrixUtils.checkMultiplicationCompatible(this, m);

            int nRows = this.getRowDimension();
            int nCols = m.getColumnDimension();
            int nSum = this.getColumnDimension();

            double[][] outData = new double[nRows][];
            // Will hold a column of "m".
            double[] mCol = new double[nSum];
            double[][] mData = m.data;

            // Multiply.
            for (int col = 0; col < nCols; col++)
            {
                // Copy all elements of column "col" of "m" so that
                // will be in contiguous memory.
                for (int mRow = 0; mRow < nSum; mRow++)
                {
                    mCol[mRow] = mData[mRow][col];
                }

                for (int row = 0; row < nRows; row++)
                {
                    double[] dataRow = data[row];
                    double sum = 0;
                    for (int i = 0; i < nSum; i++)
                    {
                        sum += dataRow[i] * mCol[i];
                    }
                    outData[row][col] = sum;
                }
            }

            return new Array2DRowRealMatrix(outData, false);
        }

        /// <inheritdoc/>
        public new double[][] getData()
        {
            return copyOut();
        }

        /// <summary>
        /// Get a reference to the underlying data array.
        /// </summary>
        /// <returns>2-dimensional array of entries.</returns>
        public double[][] getDataRef()
        {
            return data;
        }

        /// <inheritdoc/>
        public new void setSubMatrix(double[][] subMatrix, int row, int column)
        {
            if (data == null)
            {
                if (row > 0)
                {
                    throw new MathIllegalStateException(new LocalizedFormats("FIRST_ROWS_NOT_INITIALIZED_YET"), row);
                }
                if (column > 0)
                {
                    throw new MathIllegalStateException(new LocalizedFormats("FIRST_COLUMNS_NOT_INITIALIZED_YET"), column);
                }
                MathUtils.checkNotNull(subMatrix);
                int nRows = subMatrix.Length;
                if (nRows == 0)
                {
                    throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_ROW"));
                }

                int nCols = subMatrix[0].Length;
                if (nCols == 0)
                {
                    throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_COLUMN"));
                }
                data = new double[subMatrix.Length][];
                for (int i = 0; i < data.Length; ++i)
                {
                    if (subMatrix[i].Length != nCols)
                    {
                        throw new DimensionMismatchException(subMatrix[i].Length, nCols);
                    }
                    Array.Copy(subMatrix[i], 0, data[i + row], column, nCols);
                }
            }
            else
            {
                base.setSubMatrix(subMatrix, row, column);
            }

        }

        /// <inheritdoc/>
        public override double getEntry(int row, int column)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            return data[row][column];
        }

        /// <inheritdoc/>
        public override void setEntry(int row, int column, double value)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            data[row][column] = value;
        }

        /// <inheritdoc/>
        public new void addToEntry(int row, int column, double increment)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            data[row][column] += increment;
        }

        /// <inheritdoc/>
        public new void multiplyEntry(int row, int column, double factor)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            data[row][column] *= factor;
        }

        /// <inheritdoc/>
        public override int getRowDimension()
        {
            return (data == null) ? 0 : data.Length;
        }

        /// <inheritdoc/>
        public override int getColumnDimension()
        {
            return ((data == null) || (data[0] == null)) ? 0 : data[0].Length;
        }

        /// <inheritdoc/>
        public new double[] operate(double[] v)
        {
            int nRows = this.getRowDimension();
            int nCols = this.getColumnDimension();
            if (v.Length != nCols)
            {
                throw new DimensionMismatchException(v.Length, nCols);
            }
            double[] outp = new double[nRows];
            for (int row = 0; row < nRows; row++)
            {
                double[] dataRow = data[row];
                double sum = 0;
                for (int i = 0; i < nCols; i++)
                {
                    sum += dataRow[i] * v[i];
                }
                outp[row] = sum;
            }
            return outp;
        }

        /// <inheritdoc/>
        public new double[] preMultiply(double[] v)
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (v.Length != nRows)
            {
                throw new DimensionMismatchException(v.Length, nRows);
            }

            double[] outp = new double[nCols];
            for (int col = 0; col < nCols; ++col)
            {
                double sum = 0;
                for (int i = 0; i < nRows; ++i)
                {
                    sum += data[i][col] * v[i];
                }
                outp[col] = sum;
            }

            return outp;

        }

        /// <inheritdoc/>
        public new double walkInRowOrder(RealMatrixChangingVisitor visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int i = 0; i < rows; ++i)
            {
                double[] rowI = data[i];
                for (int j = 0; j < columns; ++j)
                {
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInRowOrder(RealMatrixPreservingVisitor visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int i = 0; i < rows; ++i)
            {
                double[] rowI = data[i];
                for (int j = 0; j < columns; ++j)
                {
                    visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInRowOrder(RealMatrixChangingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int i = startRow; i <= endRow; ++i)
            {
                double[] rowI = data[i];
                for (int j = startColumn; j <= endColumn; ++j)
                {
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInRowOrder(RealMatrixPreservingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int i = startRow; i <= endRow; ++i)
            {
                double[] rowI = data[i];
                for (int j = startColumn; j <= endColumn; ++j)
                {
                    visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInColumnOrder(RealMatrixChangingVisitor visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int j = 0; j < columns; ++j)
            {
                for (int i = 0; i < rows; ++i)
                {
                    double[] rowI = data[i];
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInColumnOrder(RealMatrixPreservingVisitor visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int j = 0; j < columns; ++j)
            {
                for (int i = 0; i < rows; ++i)
                {
                    visitor.visit(i, j, data[i][j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInColumnOrder(RealMatrixChangingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int j = startColumn; j <= endColumn; ++j)
            {
                for (int i = startRow; i <= endRow; ++i)
                {
                    double[] rowI = data[i];
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInColumnOrder(RealMatrixPreservingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int j = startColumn; j <= endColumn; ++j)
            {
                for (int i = startRow; i <= endRow; ++i)
                {
                    visitor.visit(i, j, data[i][j]);
                }
            }
            return visitor.end();
        }

        /// <summary>
        /// Get a fresh copy of the underlying data array.
        /// </summary>
        /// <returns>a copy of the underlying data array.</returns>
        private double[][] copyOut()
        {
            int nRows = this.getRowDimension();
            double[][] outp = new double[nRows][];
            // can't copy 2-d array in one shot, otherwise get row references
            for (int i = 0; i < nRows; i++)
            {
                Array.Copy(data[i], 0, outp[i], 0, data[i].Length);
            }
            return outp;
        }

        /// <summary>
        /// Replace data with a fresh copy of the input array.
        /// </summary>
        /// <param name="inp">Data to copy.</param>
        /// <exception cref="DimensionMismatchException"> if the input array not rectangular.
        /// </exception>
        /// <exception cref="NullArgumentException"> if the input array is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if the input array is empty.</exception>
        private void copyIn(double[][] inp)
        {
            setSubMatrix(inp, 0, 0);
        }
    }
}
