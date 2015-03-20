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
    /// Implementation of FieldMatrix<T> using a <see cref="FieldElement"/>[][] array to
    /// store entries.
    /// <para>
    /// As specified in the <see cref="FieldMatrix"/> interface, matrix element indexing
    /// is 0-based -- e.g., <c>getEntry(0, 0)</c>
    /// returns the element in the first row, first column of the matrix.
    /// </para>
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public class Array2DRowFieldMatrix<T> : AbstractFieldMatrix<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Entries of the matrix
        /// </summary>
        private T[][] data;

        /// <summary>
        /// Creates a matrix with no data
        /// </summary>
        /// <param name="field">field to which the elements belong</param>
        public Array2DRowFieldMatrix(Field<T> field) : base(field) { }

        /// <summary>
        /// Create a new <c>FieldMatrix(T)</c> with the supplied row and column dimensions.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="rowDimension">Number of rows in the new matrix.</param>
        /// <param name="columnDimension">Number of columns in the new matrix.</param>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension
        /// is not positive.</exception>
        public Array2DRowFieldMatrix(Field<T> field, int rowDimension, int columnDimension)
            : base(field, rowDimension, columnDimension)
        {
            data = MathArrays.buildArray(field, rowDimension, columnDimension);
        }

        /// <summary>
        /// Create a new <c>FieldMatrix(T)</c> using the input array as the underlying
        /// data array.
        /// <para>The input array is copied, not referenced. This constructor has
        /// the same effect as calling 
        /// <see cref="Array2DRowFieldMatrix(FieldElement[][], boolean)"/>
        /// with the second argument set to <c>true</c>.</para>
        /// </summary>
        /// <param name="d">Data for the new matrix.</param>
        /// <exception cref="DimensionMismatchException"> if <c>d</c> is not rectangular.
        /// </exception>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if there are not at least one row and one column.
        /// </exception>
        /// <remarks>
        /// See <see cref="Array2DRowFieldMatrix(FieldElement[][], boolean)"/>
        /// </remarks>
        public Array2DRowFieldMatrix(T[][] d) : this(extractField(d), d) { }

        /// <summary>
        /// Create a new <c>FieldMatrix<T></c> using the input array as the underlying
        /// data array.
        /// <para>The input array is copied, not referenced. This constructor has
        /// the same effect as calling 
        /// <see cref="Array2DRowFieldMatrix(FieldElement[][], boolean)"/>
        /// with the second argument set to <c>true</c>.</para>
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="d">Data for the new matrix.</param>
        /// <exception cref="DimensionMismatchException"> if <c>d</c> is not rectangular.
        /// </exception>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if there are not at least one row and one column.
        /// </exception>
        /// <remarks>
        /// See <see cref="Array2DRowFieldMatrix(FieldElement[][], boolean)"/>
        /// </remarks>
        public Array2DRowFieldMatrix(Field<T> field, T[][] d)
            : base(field)
        {
            copyIn(d);
        }

        /// <summary>
        /// Create a new <c>FieldMatrix<T></c> using the input array as the underlying
        /// data array.
        /// <para>If an array is built specially in order to be embedded in a
        /// <c>FieldMatrix<T></c> and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>. This will prevent the copying and improve
        /// performance as no new array will be built and no data will be copied.</para>
        /// </summary>
        /// <param name="d">Data for the new matrix.</param>
        /// <param name="copyArray">Whether to copy or reference the input array.</param>
        /// <exception cref="DimensionMismatchException"> if <c>d</c> is not rectangular.
        /// </exception>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if there are not at least one row and one column.
        /// </exception>
        /// <remarks>
        /// See <see cref="Array2DRowFieldMatrix(FieldElement[][], boolean)"/>
        /// </remarks>
        public Array2DRowFieldMatrix(T[][] d, Boolean copyArray) : this(extractField(d), d, copyArray) { }

        /// <summary>
        /// Create a new <c>FieldMatrix<T></c> using the input array as the underlying
        /// data array.
        /// <para>If an array is built specially in order to be embedded in a
        /// <c>FieldMatrix<T></c> and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>. This will prevent the copying and improve
        /// performance as no new array will be built and no data will be copied.</para>
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="d">Data for the new matrix.</param>
        /// <param name="copyArray">Whether to copy or reference the input array.</param>
        /// <exception cref="DimensionMismatchException"> if <c>d</c> is not rectangular.
        /// </exception>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if there are not at least one row and one column.
        /// </exception>
        /// <remarks>
        /// See <see cref="Array2DRowFieldMatrix(FieldElement[][], boolean)"/>
        /// </remarks>
        public Array2DRowFieldMatrix(Field<T> field, T[][] d, Boolean copyArray)
            : base(field)
        {
            if (copyArray)
            {
                copyIn(d);
            }
            else
            {
                MathUtils.checkNotNull(d);
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
                        throw new DimensionMismatchException(nCols, d[r].Length);
                    }
                }
                data = d;
            }
        }

        /// <summary>
        /// Create a new (column) <c>FieldMatrix(T)</c> using <c>v</c> as the
        /// data for the unique column of the created matrix.
        /// The input array is copied.
        /// </summary>
        /// <param name="v">Column vector holding data for new matrix.</param>
        /// <exception cref="NoDataException"> if v is empty</exception>
        public Array2DRowFieldMatrix(T[] v) : this(extractField(v), v) { }

        /// <summary>
        /// Create a new (column) <c>FieldMatrix<T></c> using <c>v</c> as the
        /// data for the unique column of the created matrix.
        /// The input array is copied.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="v">Column vector holding data for new matrix.</param>
        public Array2DRowFieldMatrix(Field<T> field, T[] v)
            : base(field)
        {
            int nRows = v.Length;
            data = MathArrays.buildArray(getField(), nRows, 1);
            for (int row = 0; row < nRows; row++)
            {
                data[row][0] = v[row];
            }
        }

        /// <inheritdoc/>
        public override FieldMatrix<T> createMatrix(int rowDimension, int columnDimension)
        {
            return new Array2DRowFieldMatrix<T>(getField(), rowDimension, columnDimension);
        }

        /// <inheritdoc/>
        public override FieldMatrix<T> copy()
        {
            return new Array2DRowFieldMatrix<T>(getField(), copyOut(), false);
        }

        /// <summary>
        /// Add <c>m</c> to this matrix.
        /// </summary>
        /// <param name="m">Matrix to be added.</param>
        /// <returns><c>this</c> + m.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as this matrix.</exception>
        public Array2DRowFieldMatrix<T> add(Array2DRowFieldMatrix<T> m)
        {
            // safety check
            checkAdditionCompatible(m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            T[][] outData = MathArrays.buildArray(getField(), rowCount, columnCount);
            for (int row = 0; row < rowCount; row++)
            {
                T[] dataRow = data[row];
                T[] mRow = m.data[row];
                T[] outDataRow = outData[row];
                for (int col = 0; col < columnCount; col++)
                {
                    outDataRow[col] = dataRow[col].add(mRow[col]);
                }
            }

            return new Array2DRowFieldMatrix<T>(getField(), outData, false);
        }

        /// <summary>
        /// Subtract <c>m</c> from this matrix.
        /// </summary>
        /// <param name="m">Matrix to be subtracted.</param>
        /// <returns><c>this</c> + m.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as this matrix.</exception>
        public Array2DRowFieldMatrix<T> subtract(Array2DRowFieldMatrix<T> m)
        {
            // safety check
            checkSubtractionCompatible(m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            T[][] outData = MathArrays.buildArray(getField(), rowCount, columnCount);
            for (int row = 0; row < rowCount; row++)
            {
                T[] dataRow = data[row];
                T[] mRow = m.data[row];
                T[] outDataRow = outData[row];
                for (int col = 0; col < columnCount; col++)
                {
                    outDataRow[col] = dataRow[col].subtract(mRow[col]);
                }
            }

            return new Array2DRowFieldMatrix<T>(getField(), outData, false);

        }

        /// <summary>
        /// Postmultiplying this matrix by <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to postmultiply by.</param>
        /// <returns><c>this</c> * m.</returns>
        /// <exception cref="DimensionMismatchException"> if the number of columns of this
        /// matrix is not equal to the number of rows of <c>m</c>.</exception>
        public Array2DRowFieldMatrix<T> multiply(Array2DRowFieldMatrix<T> m)
        {
            // safety check
            checkMultiplicationCompatible(m);

            int nRows = this.getRowDimension();
            int nCols = m.getColumnDimension();
            int nSum = this.getColumnDimension();
            T[][] outData = MathArrays.buildArray(getField(), nRows, nCols);
            for (int row = 0; row < nRows; row++)
            {
                T[] dataRow = data[row];
                T[] outDataRow = outData[row];
                for (int col = 0; col < nCols; col++)
                {
                    T sum = getField().getZero();
                    for (int i = 0; i < nSum; i++)
                    {
                        sum = sum.add(dataRow[i].multiply(m.data[i][col]));
                    }
                    outDataRow[col] = sum;
                }
            }

            return new Array2DRowFieldMatrix<T>(getField(), outData, false);

        }

        /// <inheritdoc/>
        public new T[][] getData()
        {
            return copyOut();
        }

        /// <summary>
        /// Get a reference to the underlying data array.
        /// This methods returns internal data, not fresh copy of it.
        /// </summary>
        /// <returns>the 2-dimensional array of entries.</returns>
        public T[][] getDataRef()
        {
            return data;
        }

        /// <inheritdoc/>
        public new void setSubMatrix(T[][] subMatrix, int row, int column)
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
                data = MathArrays.buildArray(getField(), subMatrix.Length, nCols);
                for (int i = 0; i < data.Length; ++i)
                {
                    if (subMatrix[i].Length != nCols)
                    {
                        throw new DimensionMismatchException(nCols, subMatrix[i].Length);
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
        public override T getEntry(int row, int column)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            return data[row][column];
        }

        /// <inheritdoc/>
        public override void setEntry(int row, int column, T value)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            data[row][column] = value;
        }

        /// <inheritdoc/>
        public override void addToEntry(int row, int column, T increment)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            data[row][column] = data[row][column].add(increment);
        }

        /// <inheritdoc/>
        public override void multiplyEntry(int row, int column, T factor)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            data[row][column] = data[row][column].multiply(factor);
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
        public new T[] operate(T[] v)
        {
            int nRows = this.getRowDimension();
            int nCols = this.getColumnDimension();
            if (v.Length != nCols)
            {
                throw new DimensionMismatchException(v.Length, nCols);
            }
            T[] outp = MathArrays.buildArray(getField(), nRows);
            for (int row = 0; row < nRows; row++)
            {
                T[] dataRow = data[row];
                T sum = getField().getZero();
                for (int i = 0; i < nCols; i++)
                {
                    sum = sum.add(dataRow[i].multiply(v[i]));
                }
                outp[row] = sum;
            }
            return outp;
        }

        /// <inheritdoc/>
        public new T[] preMultiply(T[] v)
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (v.Length != nRows)
            {
                throw new DimensionMismatchException(v.Length, nRows);
            }

            T[] outp = MathArrays.buildArray(getField(), nCols);
            for (int col = 0; col < nCols; ++col)
            {
                T sum = getField().getZero();
                for (int i = 0; i < nRows; ++i)
                {
                    sum = sum.add(data[i][col].multiply(v[i]));
                }
                outp[col] = sum;
            }

            return outp;
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int i = 0; i < rows; ++i)
            {
                T[] rowI = data[i];
                for (int j = 0; j < columns; ++j)
                {
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int i = 0; i < rows; ++i)
            {
                T[] rowI = data[i];
                for (int j = 0; j < columns; ++j)
                {
                    visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int i = startRow; i <= endRow; ++i)
            {
                T[] rowI = data[i];
                for (int j = startColumn; j <= endColumn; ++j)
                {
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int i = startRow; i <= endRow; ++i)
            {
                T[] rowI = data[i];
                for (int j = startColumn; j <= endColumn; ++j)
                {
                    visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInColumnOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int j = 0; j < columns; ++j)
            {
                for (int i = 0; i < rows; ++i)
                {
                    T[] rowI = data[i];
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInColumnOrder(FieldMatrixPreservingVisitor<T> visitor)
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
        public new T walkInColumnOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int j = startColumn; j <= endColumn; ++j)
            {
                for (int i = startRow; i <= endRow; ++i)
                {
                    T[] rowI = data[i];
                    rowI[j] = visitor.visit(i, j, rowI[j]);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInColumnOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
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
        private T[][] copyOut()
        {
            int nRows = this.getRowDimension();
            T[][] outp = MathArrays.buildArray(getField(), nRows, getColumnDimension());
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
        private void copyIn(T[][] inp)
        {
            setSubMatrix(inp, 0, 0);
        }
    }
}
