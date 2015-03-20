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
using System.Collections.Generic;
using System.Text;

namespace Math3.linear
{
    /// <summary>
    /// Basic implementation of <see cref="FieldMatrix"/> methods regardless of the underlying
    /// storage.
    /// <para>All the methods implemented here use <see cref="getEntry(int, int)"/> to access
    /// matrix elements. Derived class can provide faster implementations. </para>
    /// </summary>
    /// <typeparam name="T">Type of the field elements.</typeparam>
    public abstract class AbstractFieldMatrix<T> : FieldMatrix<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Field to which the elements belong.
        /// </summary>
        private readonly Field<T> field;

        /// <summary>
        /// Constructor for use with Serializable
        /// </summary>
        protected AbstractFieldMatrix()
        {
            field = null;
        }

        /// <summary>
        /// Creates a matrix with no data
        /// </summary>
        /// <param name="field">field to which the elements belong</param>
        protected AbstractFieldMatrix(Field<T> field)
        {
            this.field = field;
        }

        /// <summary>
        /// Create a new FieldMatrix(T) with the supplied row and column dimensions.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="rowDimension">Number of rows in the new matrix.</param>
        /// <param name="columnDimension">Number of columns in the new matrix.</param>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension is not
        /// positive.</exception>
        protected AbstractFieldMatrix(Field<T> field, int rowDimension, int columnDimension)
        {
            if (rowDimension <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("DIMENSION"), rowDimension);
            }
            if (columnDimension <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("DIMENSION"), columnDimension);
            }
            this.field = field;
        }

        /// <summary>
        /// Get the elements type from an array.
        /// </summary>
        /// <typeparam name="U">Type of the field elements.</typeparam>
        /// <param name="d">Data array.</param>
        /// <returns>the field to which the array elements belong.</returns>
        /// <exception cref="NullArgumentException"> if the array is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if the array is empty.</exception>
        protected static Field<U> extractField<U>(U[][] d) where U : FieldElement<U>
        {
            if (d == null)
            {
                throw new NullArgumentException();
            }
            if (d.Length == 0)
            {
                throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_ROW"));
            }
            if (d[0].Length == 0)
            {
                throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_COLUMN"));
            }
            return d[0][0].getField();
        }

        /// <summary>
        /// Get the elements type from an array.
        /// </summary>
        /// <typeparam name="U">Type of the field elements.</typeparam>
        /// <param name="d">Data array.</param>
        /// <returns>the field to which the array elements belong.</returns>
        /// <exception cref="NoDataException"> if array is empty.</exception>
        protected static Field<U> extractField<U>(U[] d) where U : FieldElement<U>
        {
            if (d.Length == 0)
            {
                throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_ROW"));
            }
            return d[0].getField();
        }

        /// <summary>
        /// Build an array of elements.
        /// <para>
        /// Complete arrays are filled with field.getZero()
        /// </para>
        /// </summary>
        /// <typeparam name="U">Type of the field elements</typeparam>
        /// <param name="field">field to which array elements belong</param>
        /// <param name="rows">number of rows</param>
        /// <param name="columns">number of columns (may be negative to build partial
        /// arrays in the same way <c>new Field[rows][]</c> works)</param>
        /// <returns>a new array</returns>
        [Obsolete("replaced by MathArrays.buildArray(Field, int int")]
        protected static U[][] buildArray<U>(Field<U> field, int rows, int columns) where U : FieldElement<U>
        {
            return MathArrays.buildArray(field, rows, columns);
        }

        /// <summary>
        /// Build an array of elements.
        /// <para>
        /// Arrays are filled with field.getZero()
        /// </para>
        /// </summary>
        /// <typeparam name="U">the type of the field elements</typeparam>
        /// <param name="field">field to which array elements belong</param>
        /// <param name="length">of the array</param>
        /// <returns>a new array</returns>
        [Obsolete("replaced by MathArrays.buildArray(Field, int)")]
        protected static U[] buildArray<U>(Field<U> field, int length) where U : FieldElement<U>
        {
            return MathArrays.buildArray(field, length);
        }

        /// <inheritdoc/>
        public Field<T> getField()
        {
            return field;
        }

        /// <inheritdoc/>
        public abstract FieldMatrix<T> createMatrix(int rowDimension, int columnDimension);

        /// <inheritdoc/>
        public abstract FieldMatrix<T> copy();

        /// <inheritdoc/>
        public FieldMatrix<T> add(FieldMatrix<T> m)
        {
            // safety check
            checkAdditionCompatible(m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col).add(m.getEntry(row, col)));
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public FieldMatrix<T> subtract(FieldMatrix<T> m)
        {
            // safety check
            checkSubtractionCompatible(m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col).subtract(m.getEntry(row, col)));
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public FieldMatrix<T> scalarAdd(T d)
        {

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col).add(d));
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public FieldMatrix<T> scalarMultiply(T d)
        {
            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col).multiply(d));
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public FieldMatrix<T> multiply(FieldMatrix<T> m)
        {
            // safety check
            checkMultiplicationCompatible(m);

            int nRows = getRowDimension();
            int nCols = m.getColumnDimension();
            int nSum = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(nRows, nCols);
            for (int row = 0; row < nRows; ++row)
            {
                for (int col = 0; col < nCols; ++col)
                {
                    T sum = field.getZero();
                    for (int i = 0; i < nSum; ++i)
                    {
                        sum = sum.add(getEntry(row, i).multiply(m.getEntry(i, col)));
                    }
                    outp.setEntry(row, col, sum);
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public FieldMatrix<T> preMultiply(FieldMatrix<T> m)
        {
            return m.multiply(this);
        }

        /// <inheritdoc/>
        public FieldMatrix<T> power(int p)
        {
            if (p < 0)
            {
                throw new NotPositiveException<Int32>(p);
            }

            if (!isSquare())
            {
                throw new NonSquareMatrixException(getRowDimension(), getColumnDimension());
            }

            if (p == 0)
            {
                return MatrixUtils.createFieldIdentityMatrix(this.getField(), this.getRowDimension());
            }

            if (p == 1)
            {
                return this.copy();
            }

            int power = p - 1;

            /*
             * Only log_2(p) operations is used by doing as follows:
             * 5^214 = 5^128 * 5^64 * 5^16 * 5^4 * 5^2
             *
             * In general, the same approach is used for A^p.
             */

            char[] binaryRepresentation = Convert.ToString(power, 2).ToCharArray();
            List<Int32> nonZeroPositions = new List<Int32>();

            for (int i = 0; i < binaryRepresentation.Length; ++i)
            {
                if (binaryRepresentation[i] == '1')
                {
                    int pos = binaryRepresentation.Length - i - 1;
                    nonZeroPositions.Add(pos);
                }
            }

            List<FieldMatrix<T>> results = new List<FieldMatrix<T>>(binaryRepresentation.Length);

            results[0] = this.copy();

            for (int i = 1; i < binaryRepresentation.Length; ++i)
            {
                FieldMatrix<T> s = results[i - 1];
                FieldMatrix<T> r = s.multiply(s);
                results[i] = r;
            }

            FieldMatrix<T> result = this.copy();

            foreach (int i in nonZeroPositions)
            {
                result = result.multiply(results[i]);
            }

            return result;
        }

        /// <inheritdoc/>
        public T[][] getData()
        {
            T[][] data = MathArrays.buildArray(field, getRowDimension(), getColumnDimension());

            for (int i = 0; i < data.Length; ++i)
            {
                T[] dataI = data[i];
                for (int j = 0; j < dataI.Length; ++j)
                {
                    dataI[j] = getEntry(i, j);
                }
            }

            return data;
        }

        /// <inheritdoc/>
        public FieldMatrix<T> getSubMatrix(int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);

            FieldMatrix<T> subMatrix = createMatrix(endRow - startRow + 1, endColumn - startColumn + 1);
            for (int i = startRow; i <= endRow; ++i)
            {
                for (int j = startColumn; j <= endColumn; ++j)
                {
                    subMatrix.setEntry(i - startRow, j - startColumn, getEntry(i, j));
                }
            }

            return subMatrix;

        }

        /// <inheritdoc/>
        public FieldMatrix<T> getSubMatrix(int[] selectedRows, int[] selectedColumns)
        {
            // safety checks
            checkSubMatrixIndex(selectedRows, selectedColumns);

            // copy entries
            FieldMatrix<T> subMatrix =
                createMatrix(selectedRows.Length, selectedColumns.Length);
            subMatrix.walkInOptimizedOrder(new DefaultFieldMatrixChangingVisitorAnonymous1<T>(field.getZero(), selectedRows, selectedColumns, this));

            return subMatrix;

        }

        private class DefaultFieldMatrixChangingVisitorAnonymous1<U> : DefaultFieldMatrixChangingVisitor<U> where U : FieldElement<U>
        {
            private readonly int[] selectedRows;
            private readonly int[] selectedColumns;
            private readonly AbstractFieldMatrix<U> m;
            public DefaultFieldMatrixChangingVisitorAnonymous1(U zero, int[] selectedRows, int[] selectedColumns, AbstractFieldMatrix<U> m)
                : base(zero)
            {
                this.selectedRows = selectedRows;
                this.selectedColumns = selectedColumns;
                this.m = m;
            }
            /// <inheritdoc/>
            public new U visit(int row, int column, U value)
            {
                return m.getEntry(selectedRows[row], selectedColumns[column]);
            }
        }

        /// <inheritdoc/>
        public void copySubMatrix(int startRow, int endRow, int startColumn, int endColumn, T[][] destination)
        {
            // safety checks
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            int rowsCount = endRow + 1 - startRow;
            int columnsCount = endColumn + 1 - startColumn;
            if ((destination.Length < rowsCount) || (destination[0].Length < columnsCount))
            {
                throw new MatrixDimensionMismatchException(destination.Length,
                                                           destination[0].Length,
                                                           rowsCount,
                                                           columnsCount);
            }

            // copy entries
            walkInOptimizedOrder(new DefaultFieldMatrixPreservingVisitor<T>(field.getZero())
            {



            }, startRow, endRow, startColumn, endColumn);

        }

        private class DefaultFieldMatrixPreservingVisitorAnonymous1<U> : DefaultFieldMatrixPreservingVisitor<U> where U : FieldElement<U>
        {
            /// <summary>
            /// Initial row index.
            /// </summary>
            private int startRow;

            /// <summary>
            /// Initial column index.
            /// </summary>
            private int startColumn;

            private U[][] destination;

            public DefaultFieldMatrixPreservingVisitorAnonymous1(U zero, U[][] destination)
                : base(zero)
            {
                this.destination = destination;
            }

            /// <inheritdoc/>
            public new void start(int rows, int columns, int startRow, int endRow, int startColumn, int endColumn)
            {
                this.startRow = startRow;
                this.startColumn = startColumn;
            }

            /// <inheritdoc/>
            public new void visit(int row, int column, U value)
            {
                destination[row - startRow][column - startColumn] = value;
            }
        }

        /// <inheritdoc/>
        public void copySubMatrix(int[] selectedRows, int[] selectedColumns, T[][] destination)
        {
            // safety checks
            checkSubMatrixIndex(selectedRows, selectedColumns);
            if ((destination.Length < selectedRows.Length) ||
                (destination[0].Length < selectedColumns.Length))
            {
                throw new MatrixDimensionMismatchException(destination.Length,
                                                           destination[0].Length,
                                                           selectedRows.Length,
                                                           selectedColumns.Length);
            }

            // copy entries
            for (int i = 0; i < selectedRows.Length; i++)
            {
                T[] destinationI = destination[i];
                for (int j = 0; j < selectedColumns.Length; j++)
                {
                    destinationI[j] = getEntry(selectedRows[i], selectedColumns[j]);
                }
            }

        }

        /// <inheritdoc/>
        public void setSubMatrix(T[][] subMatrix, int row, int column)
        {
            if (subMatrix == null)
            {
                throw new NullArgumentException();
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

            for (int r = 1; r < nRows; ++r)
            {
                if (subMatrix[r].Length != nCols)
                {
                    throw new DimensionMismatchException(nCols, subMatrix[r].Length);
                }
            }

            checkRowIndex(row);
            checkColumnIndex(column);
            checkRowIndex(nRows + row - 1);
            checkColumnIndex(nCols + column - 1);

            for (int i = 0; i < nRows; ++i)
            {
                for (int j = 0; j < nCols; ++j)
                {
                    setEntry(row + i, column + j, subMatrix[i][j]);
                }
            }
        }

        /// <inheritdoc/>
        public FieldMatrix<T> getRowMatrix(int row)
        {
            checkRowIndex(row);
            int nCols = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(1, nCols);
            for (int i = 0; i < nCols; ++i)
            {
                outp.setEntry(0, i, getEntry(row, i));
            }

            return outp;

        }

        /// <inheritdoc/>
        public void setRowMatrix(int row, FieldMatrix<T> matrix)
        {
            checkRowIndex(row);
            int nCols = getColumnDimension();
            if ((matrix.getRowDimension() != 1) ||
                (matrix.getColumnDimension() != nCols))
            {
                throw new MatrixDimensionMismatchException(matrix.getRowDimension(),
                                                           matrix.getColumnDimension(),
                                                           1, nCols);
            }
            for (int i = 0; i < nCols; ++i)
            {
                setEntry(row, i, matrix.getEntry(0, i));
            }

        }

        /// <inheritdoc/>
        public FieldMatrix<T> getColumnMatrix(int column)
        {
            checkColumnIndex(column);
            int nRows = getRowDimension();
            FieldMatrix<T> outp = createMatrix(nRows, 1);
            for (int i = 0; i < nRows; ++i)
            {
                outp.setEntry(i, 0, getEntry(i, column));
            }

            return outp;

        }

        /// <inheritdoc/>
        public void setColumnMatrix(int column, FieldMatrix<T> matrix)
        {
            checkColumnIndex(column);
            int nRows = getRowDimension();
            if ((matrix.getRowDimension() != nRows) ||
                (matrix.getColumnDimension() != 1))
            {
                throw new MatrixDimensionMismatchException(matrix.getRowDimension(),
                                                           matrix.getColumnDimension(),
                                                           nRows, 1);
            }
            for (int i = 0; i < nRows; ++i)
            {
                setEntry(i, column, matrix.getEntry(i, 0));
            }

        }

        /// <inheritdoc/>
        public FieldVector<T> getRowVector(int row)
        {
            return new ArrayFieldVector<T>(field, getRow(row), false);
        }

        /// <inheritdoc/>
        public void setRowVector(int row, FieldVector<T> vector)
        {
            checkRowIndex(row);
            int nCols = getColumnDimension();
            if (vector.getDimension() != nCols)
            {
                throw new MatrixDimensionMismatchException(1, vector.getDimension(),
                                                           1, nCols);
            }
            for (int i = 0; i < nCols; ++i)
            {
                setEntry(row, i, vector.getEntry(i));
            }

        }

        /// <inheritdoc/>
        public FieldVector<T> getColumnVector(int column)
        {
            return new ArrayFieldVector<T>(field, getColumn(column), false);
        }

        /// <inheritdoc/>
        public void setColumnVector(int column, FieldVector<T> vector)
        {
            checkColumnIndex(column);
            int nRows = getRowDimension();
            if (vector.getDimension() != nRows)
            {
                throw new MatrixDimensionMismatchException(vector.getDimension(), 1,
                                                           nRows, 1);
            }
            for (int i = 0; i < nRows; ++i)
            {
                setEntry(i, column, vector.getEntry(i));
            }

        }

        /// <inheritdoc/>
        public T[] getRow(int row)
        {
            checkRowIndex(row);
            int nCols = getColumnDimension();
            T[] outp = MathArrays.buildArray(field, nCols);
            for (int i = 0; i < nCols; ++i)
            {
                outp[i] = getEntry(row, i);
            }

            return outp;

        }

        /// <inheritdoc/>
        public void setRow(int row, T[] array)
        {
            checkRowIndex(row);
            int nCols = getColumnDimension();
            if (array.Length != nCols)
            {
                throw new MatrixDimensionMismatchException(1, array.Length, 1, nCols);
            }
            for (int i = 0; i < nCols; ++i)
            {
                setEntry(row, i, array[i]);
            }

        }

        /// <inheritdoc/>
        public T[] getColumn(int column)
        {
            checkColumnIndex(column);
            int nRows = getRowDimension();
            T[] outp = MathArrays.buildArray(field, nRows);
            for (int i = 0; i < nRows; ++i)
            {
                outp[i] = getEntry(i, column);
            }
            return outp;
        }

        /// <inheritdoc/>
        public void setColumn(int column, T[] array)
        {
            checkColumnIndex(column);
            int nRows = getRowDimension();
            if (array.Length != nRows)
            {
                throw new MatrixDimensionMismatchException(array.Length, 1, nRows, 1);
            }
            for (int i = 0; i < nRows; ++i)
            {
                setEntry(i, column, array[i]);
            }
        }

        /// <inheritdoc/>
        public abstract T getEntry(int row, int column);

        /// <inheritdoc/>
        public abstract void setEntry(int row, int column, T value);

        /// <inheritdoc/>
        public abstract void addToEntry(int row, int column, T increment);

        /// <inheritdoc/>
        public abstract void multiplyEntry(int row, int column, T factor);

        /// <inheritdoc/>
        public FieldMatrix<T> transpose()
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            FieldMatrix<T> outp = createMatrix(nCols, nRows);
            walkInOptimizedOrder(new DefaultFieldMatrixPreservingVisitorAnonymous2<T>(field.getZero(), outp));

            return outp;
        }

        private class DefaultFieldMatrixPreservingVisitorAnonymous2<U> : DefaultFieldMatrixPreservingVisitor<U> where U : FieldElement<U>
        {
            private readonly FieldMatrix<U> outp;
            public DefaultFieldMatrixPreservingVisitorAnonymous2(U zero, FieldMatrix<U> outp)
                : base(zero)
            {
                this.outp = outp;
            }
            /// <inheritdoc/>
            public new void visit(int row, int column, U value)
            {
                outp.setEntry(column, row, value);
            }
        }

        /// <inheritdoc/>
        public Boolean isSquare()
        {
            return getColumnDimension() == getRowDimension();
        }

        /// <inheritdoc/>
        public abstract int getRowDimension();

        /// <inheritdoc/>
        public abstract int getColumnDimension();

        /// <inheritdoc/>
        public T getTrace()
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (nRows != nCols)
            {
                throw new NonSquareMatrixException(nRows, nCols);
            }
            T trace = field.getZero();
            for (int i = 0; i < nRows; ++i)
            {
                trace = trace.add(getEntry(i, i));
            }
            return trace;
        }

        /// <inheritdoc/>
        public T[] operate(T[] v)
        {

            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (v.Length != nCols)
            {
                throw new DimensionMismatchException(v.Length, nCols);
            }

            T[] outp = MathArrays.buildArray(field, nRows);
            for (int row = 0; row < nRows; ++row)
            {
                T sum = field.getZero();
                for (int i = 0; i < nCols; ++i)
                {
                    sum = sum.add(getEntry(row, i).multiply(v[i]));
                }
                outp[row] = sum;
            }

            return outp;
        }

        /// <inheritdoc/>
        public FieldVector<T> operate(FieldVector<T> v)
        {
            try
            {
                return new ArrayFieldVector<T>(field, operate(((ArrayFieldVector<T>)v).getDataRef()), false);
            }
            catch (InvalidCastException)
            {
                int nRows = getRowDimension();
                int nCols = getColumnDimension();
                if (v.getDimension() != nCols)
                {
                    throw new DimensionMismatchException(v.getDimension(), nCols);
                }

                T[] outp = MathArrays.buildArray(field, nRows);
                for (int row = 0; row < nRows; ++row)
                {
                    T sum = field.getZero();
                    for (int i = 0; i < nCols; ++i)
                    {
                        sum = sum.add(getEntry(row, i).multiply(v.getEntry(i)));
                    }
                    outp[row] = sum;
                }

                return new ArrayFieldVector<T>(field, outp, false);
            }
        }

        /// <inheritdoc/>
        public T[] preMultiply(T[] v)
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (v.Length != nRows)
            {
                throw new DimensionMismatchException(v.Length, nRows);
            }

            T[] outp = MathArrays.buildArray(field, nCols);
            for (int col = 0; col < nCols; ++col)
            {
                T sum = field.getZero();
                for (int i = 0; i < nRows; ++i)
                {
                    sum = sum.add(getEntry(i, col).multiply(v[i]));
                }
                outp[col] = sum;
            }

            return outp;
        }

        /// <inheritdoc/>
        public FieldVector<T> preMultiply(FieldVector<T> v)
        {
            try
            {
                return new ArrayFieldVector<T>(field, preMultiply(((ArrayFieldVector<T>)v).getDataRef()), false);
            }
            catch (InvalidCastException)
            {
                int nRows = getRowDimension();
                int nCols = getColumnDimension();
                if (v.getDimension() != nRows)
                {
                    throw new DimensionMismatchException(v.getDimension(), nRows);
                }

                T[] outp = MathArrays.buildArray(field, nCols);
                for (int col = 0; col < nCols; ++col)
                {
                    T sum = field.getZero();
                    for (int i = 0; i < nRows; ++i)
                    {
                        sum = sum.add(getEntry(i, col).multiply(v.getEntry(i)));
                    }
                    outp[col] = sum;
                }

                return new ArrayFieldVector<T>(field, outp, false);
            }
        }

        /// <inheritdoc/>
        public T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int row = 0; row < rows; ++row)
            {
                for (int column = 0; column < columns; ++column)
                {
                    T oldValue = getEntry(row, column);
                    T newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int row = 0; row < rows; ++row)
            {
                for (int column = 0; column < columns; ++column)
                {
                    visitor.visit(row, column, getEntry(row, column));
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int row = startRow; row <= endRow; ++row)
            {
                for (int column = startColumn; column <= endColumn; ++column)
                {
                    T oldValue = getEntry(row, column);
                    T newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int row = startRow; row <= endRow; ++row)
            {
                for (int column = startColumn; column <= endColumn; ++column)
                {
                    visitor.visit(row, column, getEntry(row, column));
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInColumnOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int column = 0; column < columns; ++column)
            {
                for (int row = 0; row < rows; ++row)
                {
                    T oldValue = getEntry(row, column);
                    T newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInColumnOrder(FieldMatrixPreservingVisitor<T> visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int column = 0; column < columns; ++column)
            {
                for (int row = 0; row < rows; ++row)
                {
                    visitor.visit(row, column, getEntry(row, column));
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInColumnOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int column = startColumn; column <= endColumn; ++column)
            {
                for (int row = startRow; row <= endRow; ++row)
                {
                    T oldValue = getEntry(row, column);
                    T newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInColumnOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int column = startColumn; column <= endColumn; ++column)
            {
                for (int row = startRow; row <= endRow; ++row)
                {
                    visitor.visit(row, column, getEntry(row, column));
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public T walkInOptimizedOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            return walkInRowOrder(visitor);
        }

        /// <inheritdoc/>
        public T walkInOptimizedOrder(FieldMatrixPreservingVisitor<T> visitor)
        {
            return walkInRowOrder(visitor);
        }

        /// <inheritdoc/>
        public T walkInOptimizedOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            return walkInRowOrder(visitor, startRow, endRow, startColumn, endColumn);
        }

        /// <inheritdoc/>
        public T walkInOptimizedOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            return walkInRowOrder(visitor, startRow, endRow, startColumn, endColumn);
        }

        /// <summary>
        /// Get a string representation for this matrix.
        /// </summary>
        /// <returns>a string representation for this matrix</returns>
        public override String ToString()
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            StringBuilder res = new StringBuilder();
            String fullClassName = GetType().Name;
            String shortClassName = fullClassName.Substring(fullClassName.LastIndexOf('.') + 1);
            res.Append(shortClassName).Append("{");

            for (int i = 0; i < nRows; ++i)
            {
                if (i > 0)
                {
                    res.Append(",");
                }
                res.Append("{");
                for (int j = 0; j < nCols; ++j)
                {
                    if (j > 0)
                    {
                        res.Append(",");
                    }
                    res.Append(getEntry(i, j));
                }
                res.Append("}");
            }

            res.Append("}");
            return res.ToString();
        }

        /// <summary>
        /// Returns true iff <c>object</c> is a
        /// <c>FieldMatrix</c> instance with the same dimensions as this
        /// and all corresponding matrix entries are equal.
        /// </summary>
        /// <param name="obj">the object to test equality against.</param>
        /// <returns>true if object equals this</returns>
        public override Boolean Equals(Object obj)
        {
            if (obj == this)
            {
                return true;
            }
            if (obj is FieldMatrix<T>)
            {
                return false;
            }
            FieldMatrix<T> m = (FieldMatrix<T>)obj;
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (m.getColumnDimension() != nCols || m.getRowDimension() != nRows)
            {
                return false;
            }
            for (int row = 0; row < nRows; ++row)
            {
                for (int col = 0; col < nCols; ++col)
                {
                    if (!getEntry(row, col).Equals(m.getEntry(row, col)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// Computes a hashcode for the matrix.
        /// </summary>
        /// <returns>hashcode for matrix</returns>
        public override int GetHashCode()
        {
            int ret = 322562;
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            ret = ret * 31 + nRows;
            ret = ret * 31 + nCols;
            for (int row = 0; row < nRows; ++row)
            {
                for (int col = 0; col < nCols; ++col)
                {
                    ret = ret * 31 + (11 * (row + 1) + 17 * (col + 1)) * getEntry(row, col).GetHashCode();
                }
            }
            return ret;
        }

        /// <summary>
        /// Check if a row index is valid.
        /// </summary>
        /// <param name="row">Row index to check.</param>
        /// <exception cref="OutOfRangeException"> if <c>index</c> is not valid.</exception>
        protected void checkRowIndex(int row)
        {
            if (row < 0 || row >= getRowDimension())
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("ROW_INDEX"),
                                              row, 0, getRowDimension() - 1);
            }
        }

        /// <summary>
        /// Check if a columns index is valid.
        /// </summary>
        /// <param name="columns">Colums index to check.</param>
        /// <exception cref="OutOfRangeException"> if <c>index</c> is not valid.</exception>
        protected void checkColumnIndex(int column)
        {
            if (column < 0 || column >= getColumnDimension())
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("COLUMN_INDEX"), column, 0, getColumnDimension() - 1);
            }
        }

        /// <summary>
        /// Check if submatrix ranges indices are valid.
        /// Rows and columns are indicated counting from 0 to n-1.
        /// </summary>
        /// <param name="startRow">Initial row index.</param>
        /// <param name="endRow">Final row index.</param>
        /// <param name="startColumn">Initial column index.</param>
        /// <param name="endColumn">Final column index.</param>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        protected void checkSubMatrixIndex(int startRow, int endRow, int startColumn, int endColumn)
        {
            checkRowIndex(startRow);
            checkRowIndex(endRow);
            if (endRow < startRow)
            {
                throw new NumberIsTooSmallException<Int32, Int32>(new LocalizedFormats("INITIAL_ROW_AFTER_FINAL_ROW"), endRow, startRow, true);
            }

            checkColumnIndex(startColumn);
            checkColumnIndex(endColumn);
            if (endColumn < startColumn)
            {
                throw new NumberIsTooSmallException<Int32, Int32>(new LocalizedFormats("INITIAL_COLUMN_AFTER_FINAL_COLUMN"), endColumn, startColumn, true);
            }
        }

        /// <summary>
        /// Check if submatrix ranges indices are valid.
        /// Rows and columns are indicated counting from 0 to n-1.
        /// </summary>
        /// <param name="selectedRows">Array of row indices.</param>
        /// <param name="selectedColumns">Array of column indices.</param>
        /// <exception cref="NullArgumentException"> if the arrays are <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if the arrays have zero length.</exception>
        /// <exception cref="OutOfRangeException"> if row or column selections are not valid.
        /// </exception>
        protected void checkSubMatrixIndex(int[] selectedRows, int[] selectedColumns)
        {
            if (selectedRows == null ||
                selectedColumns == null)
            {
                throw new NullArgumentException();
            }
            if (selectedRows.Length == 0 ||
                selectedColumns.Length == 0)
            {
                throw new NoDataException();
            }

            foreach (int row in selectedRows)
            {
                checkRowIndex(row);
            }
            foreach (int column in selectedColumns)
            {
                checkColumnIndex(column);
            }
        }

        /// <summary>
        /// Check if a matrix is addition compatible with the instance.
        /// </summary>
        /// <param name="m">Matrix to check.</param>
        /// <exception cref="MatrixDimensionMismatchException"> if the matrix is not
        /// addition-compatible with instance.</exception>
        protected void checkAdditionCompatible(FieldMatrix<T> m)
        {
            if ((getRowDimension() != m.getRowDimension()) ||
                (getColumnDimension() != m.getColumnDimension()))
            {
                throw new MatrixDimensionMismatchException(m.getRowDimension(), m.getColumnDimension(),
                                                           getRowDimension(), getColumnDimension());
            }
        }

        /// <summary>
        /// Check if a matrix is subtraction compatible with the instance.
        /// </summary>
        /// <param name="m">Matrix to check.</param>
        /// <exception cref="MatrixDimensionMismatchExceptio">n if the matrix is not
        /// subtraction-compatible with instance.</exception>
        protected void checkSubtractionCompatible(FieldMatrix<T> m)
        {
            if ((getRowDimension() != m.getRowDimension()) ||
                (getColumnDimension() != m.getColumnDimension()))
            {
                throw new MatrixDimensionMismatchException(m.getRowDimension(), m.getColumnDimension(),
                                                           getRowDimension(), getColumnDimension());
            }
        }

        /// <summary>
        /// Check if a matrix is multiplication compatible with the instance.
        /// </summary>
        /// <param name="m">Matrix to check.</param>
        /// <exception cref="DimensionMismatchException"> if the matrix is not
        /// multiplication-compatible with instance.</exception>
        protected void checkMultiplicationCompatible(FieldMatrix<T> m)
        {
            if (getColumnDimension() != m.getRowDimension())
            {
                throw new DimensionMismatchException(m.getRowDimension(), getColumnDimension());
            }
        }
    }
}
