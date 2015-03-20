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
using System.Globalization;
using System.Text;

namespace Math3.linear
{
    /// <summary>
    /// Basic implementation of RealMatrix methods regardless of the underlying storage.
    /// <para>All the methods implemented here use <see cref="getEntry(int, int)"/> to access
    /// matrix elements. Derived class can provide faster implementations.</para>
    /// </summary>
    public abstract class AbstractRealMatrix : RealLinearOperator, RealMatrix
    {
        /// <summary>
        /// Default format.
        /// </summary>
        private static readonly RealMatrixFormat DEFAULT_FORMAT = RealMatrixFormat.getInstance(CultureInfo.InvariantCulture);
        static AbstractRealMatrix()
        {
            // set the minimum fraction digits to 1 to keep compatibility
            DEFAULT_FORMAT.getFormat().NumberDecimalDigits = 1;
        }

        /// <summary>
        /// Creates a matrix with no data
        /// </summary>
        protected AbstractRealMatrix() { }

        /// <summary>
        /// Create a new RealMatrix with the supplied row and column dimensions.
        /// </summary>
        /// <param name="rowDimension">the number of rows in the new matrix</param>
        /// <param name="columnDimension">the number of columns in the new matrix</param>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension is
        /// not positive</exception>
        protected AbstractRealMatrix(int rowDimension, int columnDimension)
        {
            if (rowDimension < 1)
            {
                throw new NotStrictlyPositiveException<Int32>(rowDimension);
            }
            if (columnDimension < 1)
            {
                throw new NotStrictlyPositiveException<Int32>(columnDimension);
            }
        }

        /// <inheritdoc/>
        public RealMatrix add(RealMatrix m)
        {
            MatrixUtils.checkAdditionCompatible(this, m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            RealMatrix outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col) + m.getEntry(row, col));
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public RealMatrix subtract(RealMatrix m)
        {
            MatrixUtils.checkSubtractionCompatible(this, m);

            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            RealMatrix outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col) - m.getEntry(row, col));
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public RealMatrix scalarAdd(double d)
        {
            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            RealMatrix outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col) + d);
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public RealMatrix scalarMultiply(double d)
        {
            int rowCount = getRowDimension();
            int columnCount = getColumnDimension();
            RealMatrix outp = createMatrix(rowCount, columnCount);
            for (int row = 0; row < rowCount; ++row)
            {
                for (int col = 0; col < columnCount; ++col)
                {
                    outp.setEntry(row, col, getEntry(row, col) * d);
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public RealMatrix multiply(RealMatrix m)
        {
            MatrixUtils.checkMultiplicationCompatible(this, m);

            int nRows = getRowDimension();
            int nCols = m.getColumnDimension();
            int nSum = getColumnDimension();
            RealMatrix outp = createMatrix(nRows, nCols);
            for (int row = 0; row < nRows; ++row)
            {
                for (int col = 0; col < nCols; ++col)
                {
                    double sum = 0;
                    for (int i = 0; i < nSum; ++i)
                    {
                        sum += getEntry(row, i) * m.getEntry(i, col);
                    }
                    outp.setEntry(row, col, sum);
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public RealMatrix preMultiply(RealMatrix m)
        {
            return m.multiply(this);
        }

        /// <inheritdoc/>
        public RealMatrix power(int p)
        {
            if (p < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("NOT_POSITIVE_EXPONENT"), p);
            }

            if (!isSquare())
            {
                throw new NonSquareMatrixException(getRowDimension(), getColumnDimension());
            }

            if (p == 0)
            {
                return MatrixUtils.createRealIdentityMatrix(this.getRowDimension());
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
            int maxI = -1;

            for (int i = 0; i < binaryRepresentation.Length; ++i)
            {
                if (binaryRepresentation[i] == '1')
                {
                    int pos = binaryRepresentation.Length - i - 1;
                    nonZeroPositions.Add(pos);

                    // The positions are taken in turn, so maxI is only changed once
                    if (maxI == -1)
                    {
                        maxI = pos;
                    }
                }
            }

            RealMatrix[] results = new RealMatrix[maxI + 1];
            results[0] = this.copy();

            for (int i = 1; i <= maxI; ++i)
            {
                results[i] = results[i - 1].multiply(results[i - 1]);
            }

            RealMatrix result = this.copy();

            foreach (Int32 i in nonZeroPositions)
            {
                result = result.multiply(results[i]);
            }

            return result;
        }

        /// <inheritdoc/>
        public double[][] getData()
        {
            double[][] data = new double[getRowDimension()][];

            for (int i = 0; i < data.Length; ++i)
            {
                double[] dataI = data[i];
                for (int j = 0; j < dataI.Length; ++j)
                {
                    dataI[j] = getEntry(i, j);
                }
            }

            return data;
        }

        /// <inheritdoc/>
        public double getNorm()
        {
            return walkInColumnOrder(new RealMatrixPreservingVisitorAnonymous1());
        }

        private class RealMatrixPreservingVisitorAnonymous1 : RealMatrixPreservingVisitor
        {
            /// <summary>
            /// Last row index.
            /// </summary>
            private double endRow;

            /// <summary>
            /// Sum of absolute values on one column.
            /// </summary>
            private double columnSum;

            /// <summary>
            /// Maximal sum across all columns.
            /// </summary>
            private double maxColSum;

            /// <inheritdoc/>
            public void start(int rows, int columns, int startRow, int endRow, int startColumn, int endColumn)
            {
                this.endRow = endRow;
                columnSum = 0;
                maxColSum = 0;
            }

            /// <inheritdoc/>
            public void visit(int row, int column, double value)
            {
                columnSum += FastMath.abs(value);
                if (row == endRow)
                {
                    maxColSum = FastMath.max(maxColSum, columnSum);
                    columnSum = 0;
                }
            }

            /// <inheritdoc/>
            public double end()
            {
                return maxColSum;
            }
        }

        /// <inheritdoc/>
        public double getFrobeniusNorm()
        {
            return walkInOptimizedOrder(new RealMatrixPreservingVisitorAnonymous2());
        }

        private class RealMatrixPreservingVisitorAnonymous2 : RealMatrixPreservingVisitor
        {
            /// <summary>
            /// Sum of squared entries.
            /// </summary>
            private double sum;

            /// <inheritdoc/>
            public void start(int rows, int columns, int startRow, int endRow, int startColumn, int endColumn)
            {
                sum = 0;
            }

            /// <inheritdoc/>
            public void visit(int row, int column, double value)
            {
                sum += value * value;
            }

            /// <inheritdoc/>
            public double end()
            {
                return FastMath.sqrt(sum);
            }
        }

        /// <inheritdoc/>
        public RealMatrix getSubMatrix(int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);

            RealMatrix subMatrix =
                createMatrix(endRow - startRow + 1, endColumn - startColumn + 1);
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
        public RealMatrix getSubMatrix(int[] selectedRows, int[] selectedColumns)
        {
            MatrixUtils.checkSubMatrixIndex(this, selectedRows, selectedColumns);

            RealMatrix subMatrix = createMatrix(selectedRows.Length, selectedColumns.Length);
            subMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitorAnonymous1(ref selectedRows, ref selectedColumns, this));

            return subMatrix;
        }

        private class DefaultRealMatrixChangingVisitorAnonymous1 : DefaultRealMatrixChangingVisitor
        {
            private readonly int[] selectedRows;
            private readonly int[] selectedColumns;
            private readonly AbstractRealMatrix m;
            public DefaultRealMatrixChangingVisitorAnonymous1(ref int[] selectedRows, ref int[] selectedColumns, AbstractRealMatrix m)
            {
                this.selectedRows = selectedRows;
                this.selectedColumns = selectedColumns;
                this.m = m;
            }
            /// <inheritdoc/>
            public new double visit(int row, int column, double value)
            {
                return m.getEntry(selectedRows[row], selectedColumns[column]);
            }
        }

        /// <inheritdoc/>
        public void copySubMatrix(int startRow, int endRow, int startColumn, int endColumn, double[][] destination)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            int rowsCount = endRow + 1 - startRow;
            int columnsCount = endColumn + 1 - startColumn;
            if ((destination.Length < rowsCount) || (destination[0].Length < columnsCount))
            {
                throw new MatrixDimensionMismatchException(destination.Length, destination[0].Length,
                                                           rowsCount, columnsCount);
            }

            for (int i = 1; i < rowsCount; i++)
            {
                if (destination[i].Length < columnsCount)
                {
                    throw new MatrixDimensionMismatchException(destination.Length, destination[i].Length,
                                                               rowsCount, columnsCount);
                }
            }

            walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitorAnonymous1(ref destination), startRow, endRow, startColumn, endColumn);
        }

        private class DefaultRealMatrixPreservingVisitorAnonymous1 : DefaultRealMatrixPreservingVisitor
        {
            /// <summary>
            /// Initial row index.
            /// </summary>
            private int startRow;

            /// <summary>
            /// Initial column index.
            /// </summary>
            private int startColumn;

            private Double[][] destination;

            public DefaultRealMatrixPreservingVisitorAnonymous1(ref Double[][] destination)
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
            public new void visit(int row, int column, double value)
            {
                destination[row - startRow][column - startColumn] = value;
            }
        }

        /// <inheritdoc/>
        public void copySubMatrix(int[] selectedRows, int[] selectedColumns, double[][] destination)
        {
            MatrixUtils.checkSubMatrixIndex(this, selectedRows, selectedColumns);
            int nCols = selectedColumns.Length;
            if ((destination.Length < selectedRows.Length) ||
                (destination[0].Length < nCols))
            {
                throw new MatrixDimensionMismatchException(destination.Length, destination[0].Length,
                                                           selectedRows.Length, selectedColumns.Length);
            }

            for (int i = 0; i < selectedRows.Length; i++)
            {
                double[] destinationI = destination[i];
                if (destinationI.Length < nCols)
                {
                    throw new MatrixDimensionMismatchException(destination.Length, destinationI.Length,
                                                               selectedRows.Length, selectedColumns.Length);
                }
                for (int j = 0; j < selectedColumns.Length; j++)
                {
                    destinationI[j] = getEntry(selectedRows[i], selectedColumns[j]);
                }
            }
        }

        /// <inheritdoc/>
        public void setSubMatrix(double[][] subMatrix, int row, int column)
        {
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

            for (int r = 1; r < nRows; ++r)
            {
                if (subMatrix[r].Length != nCols)
                {
                    throw new DimensionMismatchException(nCols, subMatrix[r].Length);
                }
            }

            MatrixUtils.checkRowIndex(this, row);
            MatrixUtils.checkColumnIndex(this, column);
            MatrixUtils.checkRowIndex(this, nRows + row - 1);
            MatrixUtils.checkColumnIndex(this, nCols + column - 1);

            for (int i = 0; i < nRows; ++i)
            {
                for (int j = 0; j < nCols; ++j)
                {
                    setEntry(row + i, column + j, subMatrix[i][j]);
                }
            }
        }

        /// <inheritdoc/>
        public RealMatrix getRowMatrix(int row)
        {
            MatrixUtils.checkRowIndex(this, row);
            int nCols = getColumnDimension();
            RealMatrix outp = createMatrix(1, nCols);
            for (int i = 0; i < nCols; ++i)
            {
                outp.setEntry(0, i, getEntry(row, i));
            }

            return outp;
        }

        /// <inheritdoc/>
        public void setRowMatrix(int row, RealMatrix matrix)
        {
            MatrixUtils.checkRowIndex(this, row);
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
        public RealMatrix getColumnMatrix(int column)
        {
            MatrixUtils.checkColumnIndex(this, column);
            int nRows = getRowDimension();
            RealMatrix outp = createMatrix(nRows, 1);
            for (int i = 0; i < nRows; ++i)
            {
                outp.setEntry(i, 0, getEntry(i, column));
            }

            return outp;
        }

        /// <inheritdoc/>
        public void setColumnMatrix(int column, RealMatrix matrix)
        {
            MatrixUtils.checkColumnIndex(this, column);
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
        public RealVector getRowVector(int row)
        {
            return new ArrayRealVector(getRow(row), false);
        }

        /// <inheritdoc/>
        public void setRowVector(int row, RealVector vector)
        {
            MatrixUtils.checkRowIndex(this, row);
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
        public RealVector getColumnVector(int column)
        {
            return new ArrayRealVector(getColumn(column), false);
        }

        /// <inheritdoc/>
        public void setColumnVector(int column, RealVector vector)
        {
            MatrixUtils.checkColumnIndex(this, column);
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
        public double[] getRow(int row)
        {
            MatrixUtils.checkRowIndex(this, row);
            int nCols = getColumnDimension();
            double[] outp = new double[nCols];
            for (int i = 0; i < nCols; ++i)
            {
                outp[i] = getEntry(row, i);
            }

            return outp;
        }

        /// <inheritdoc/>
        public void setRow(int row, double[] array)
        {
            MatrixUtils.checkRowIndex(this, row);
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
        public double[] getColumn(int column)
        {
            MatrixUtils.checkColumnIndex(this, column);
            int nRows = getRowDimension();
            double[] outp = new double[nRows];
            for (int i = 0; i < nRows; ++i)
            {
                outp[i] = getEntry(i, column);
            }

            return outp;
        }

        /// <inheritdoc/>
        public void setColumn(int column, double[] array)
        {
            MatrixUtils.checkColumnIndex(this, column);
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
        public void addToEntry(int row, int column, double increment)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            setEntry(row, column, getEntry(row, column) + increment);
        }

        /// <inheritdoc/>
        public void multiplyEntry(int row, int column, double factor)
        {
            MatrixUtils.checkMatrixIndex(this, row, column);
            setEntry(row, column, getEntry(row, column) * factor);
        }

        /// <inheritdoc/>
        public RealMatrix transpose()
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            RealMatrix outp = createMatrix(nCols, nRows);
            walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitorAnonymous2(outp));
            return outp;
        }

        private class DefaultRealMatrixPreservingVisitorAnonymous2 : DefaultRealMatrixPreservingVisitor
        {
            private readonly RealMatrix outp;
            public DefaultRealMatrixPreservingVisitorAnonymous2(RealMatrix outp)
            {
                this.outp = outp;
            }
            /// <inheritdoc/>
            public new void visit(int row, int column, double value)
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
        public double getTrace()
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (nRows != nCols)
            {
                throw new NonSquareMatrixException(nRows, nCols);
            }
            double trace = 0;
            for (int i = 0; i < nRows; ++i)
            {
                trace += getEntry(i, i);
            }
            return trace;
        }

        /// <inheritdoc/>
        public double[] operate(double[] v)
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            if (v.Length != nCols)
            {
                throw new DimensionMismatchException(v.Length, nCols);
            }

            double[] outp = new double[nRows];
            for (int row = 0; row < nRows; ++row)
            {
                double sum = 0;
                for (int i = 0; i < nCols; ++i)
                {
                    sum += getEntry(row, i) * v[i];
                }
                outp[row] = sum;
            }

            return outp;
        }

        /// <inheritdoc/>
        public override RealVector operate(RealVector v)
        {
            try
            {
                return new ArrayRealVector(operate(((ArrayRealVector)v).getDataRef()), false);
            }
            catch (InvalidCastException)
            {
                int nRows = getRowDimension();
                int nCols = getColumnDimension();
                if (v.getDimension() != nCols)
                {
                    throw new DimensionMismatchException(v.getDimension(), nCols);
                }

                double[] outp = new double[nRows];
                for (int row = 0; row < nRows; ++row)
                {
                    double sum = 0;
                    for (int i = 0; i < nCols; ++i)
                    {
                        sum += getEntry(row, i) * v.getEntry(i);
                    }
                    outp[row] = sum;
                }

                return new ArrayRealVector(outp, false);
            }
        }

        /// <inheritdoc/>
        public double[] preMultiply(double[] v)
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
                    sum += getEntry(i, col) * v[i];
                }
                outp[col] = sum;
            }

            return outp;
        }

        /// <inheritdoc/>
        public RealVector preMultiply(RealVector v)
        {
            try
            {
                return new ArrayRealVector(preMultiply(((ArrayRealVector)v).getDataRef()), false);
            }
            catch (InvalidCastException)
            {
                int nRows = getRowDimension();
                int nCols = getColumnDimension();
                if (v.getDimension() != nRows)
                {
                    throw new DimensionMismatchException(v.getDimension(), nRows);
                }

                double[] outp = new double[nCols];
                for (int col = 0; col < nCols; ++col)
                {
                    double sum = 0;
                    for (int i = 0; i < nRows; ++i)
                    {
                        sum += getEntry(i, col) * v.getEntry(i);
                    }
                    outp[col] = sum;
                }

                return new ArrayRealVector(outp, false);
            }
        }

        /// <inheritdoc/>
        public double walkInRowOrder(RealMatrixChangingVisitor visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int row = 0; row < rows; ++row)
            {
                for (int column = 0; column < columns; ++column)
                {
                    double oldValue = getEntry(row, column);
                    double newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public double walkInRowOrder(RealMatrixPreservingVisitor visitor)
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
        public double walkInRowOrder(RealMatrixChangingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int row = startRow; row <= endRow; ++row)
            {
                for (int column = startColumn; column <= endColumn; ++column)
                {
                    double oldValue = getEntry(row, column);
                    double newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public double walkInRowOrder(RealMatrixPreservingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
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
        public double walkInColumnOrder(RealMatrixChangingVisitor visitor)
        {
            int rows = getRowDimension();
            int columns = getColumnDimension();
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int column = 0; column < columns; ++column)
            {
                for (int row = 0; row < rows; ++row)
                {
                    double oldValue = getEntry(row, column);
                    double newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public double walkInColumnOrder(RealMatrixPreservingVisitor visitor)
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
        public double walkInColumnOrder(RealMatrixChangingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
            visitor.start(getRowDimension(), getColumnDimension(),
                          startRow, endRow, startColumn, endColumn);
            for (int column = startColumn; column <= endColumn; ++column)
            {
                for (int row = startRow; row <= endRow; ++row)
                {
                    double oldValue = getEntry(row, column);
                    double newValue = visitor.visit(row, column, oldValue);
                    setEntry(row, column, newValue);
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public double walkInColumnOrder(RealMatrixPreservingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
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
        public double walkInOptimizedOrder(RealMatrixChangingVisitor visitor)
        {
            return walkInRowOrder(visitor);
        }

        /// <inheritdoc/>
        public double walkInOptimizedOrder(RealMatrixPreservingVisitor visitor)
        {
            return walkInRowOrder(visitor);
        }

        /// <inheritdoc/>
        public double walkInOptimizedOrder(RealMatrixChangingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            return walkInRowOrder(visitor, startRow, endRow, startColumn, endColumn);
        }

        /// <inheritdoc/>
        public double walkInOptimizedOrder(RealMatrixPreservingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            return walkInRowOrder(visitor, startRow, endRow, startColumn, endColumn);
        }

        /// <summary>
        /// Get a string representation for this matrix.
        /// </summary>
        /// <returns>a string representation for this matrix</returns>
        public override String ToString()
        {
            StringBuilder res = new StringBuilder();
            String fullClassName = GetType().Name;
            String shortClassName = fullClassName.Substring(fullClassName.LastIndexOf('.') + 1);
            res.Append(shortClassName);
            res.Append(DEFAULT_FORMAT.format(this));
            return res.ToString();
        }

        /// <summary>
        /// Returns true iff <c>object</c> is a
        /// <c>RealMatrix</c> instance with the same dimensions as this
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
            if (!(obj is RealMatrix))
            {
                return false;
            }
            RealMatrix m = (RealMatrix)obj;
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
                    if (getEntry(row, col) != m.getEntry(row, col))
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
            int ret = 7;
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            ret = ret * 31 + nRows;
            ret = ret * 31 + nCols;
            for (int row = 0; row < nRows; ++row)
            {
                for (int col = 0; col < nCols; ++col)
                {
                    ret = ret * 31 + (11 * (row + 1) + 17 * (col + 1)) *
                        MathUtils.hash(getEntry(row, col));
                }
            }
            return ret;
        }

        /// <inheritdoc/>
        public abstract RealMatrix createMatrix(int rowDimension, int columnDimension);

        /// <inheritdoc/>
        public abstract RealMatrix copy();

        /// <inheritdoc/>
        public abstract double getEntry(int row, int column);

        /// <inheritdoc/>
        public abstract void setEntry(int row, int column, double value);
    }
}