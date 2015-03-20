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
    /// Sparse matrix implementation based on an open addressed map.
    /// <para>
    /// Caveat: This implementation assumes that, for any <c>x</c>,
    /// the equality <c>x * 0d == 0d</c> holds. But it is is not true for
    /// <c>NaN</c>. Moreover, zero entries will lose their sign.
    /// Some operations (that involve <c>NaN</c> and/or infinities) may
    /// thus give incorrect results.
    /// </para>
    /// </summary>
    public class OpenMapRealMatrix : AbstractRealMatrix, SparseRealMatrix
    {
        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        private readonly int rows;

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        private readonly int columns;

        /// <summary>
        /// Storage for (sparse) matrix elements.
        /// </summary>
        private readonly OpenIntToDoubleHashMap entries;

        /// <summary>
        /// Build a sparse matrix with the supplied row and column dimensions.
        /// </summary>
        /// <param name="rowDimension">Number of rows of the matrix.</param>
        /// <param name="columnDimension">Number of columns of the matrix.</param>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension is not
        /// positive.</exception>
        /// <exception cref="NumberIsTooLargeException"> if the total number of entries of the
        /// matrix is larger than <c>Int32.MaxValue</c>.</exception>
        public OpenMapRealMatrix(int rowDimension, int columnDimension)
            : base(rowDimension, columnDimension)
        {
            long lRow = rowDimension;
            long lCol = columnDimension;
            if (lRow * lCol >= Int32.MaxValue)
            {
                throw new NumberIsTooLargeException<Int64, Int32>(lRow * lCol, Int32.MaxValue, false);
            }
            this.rows = rowDimension;
            this.columns = columnDimension;
            this.entries = new OpenIntToDoubleHashMap(0.0);
        }

        /// <summary>
        /// Build a matrix by copying another one.
        /// </summary>
        /// <param name="matrix">matrix to copy.</param>
        public OpenMapRealMatrix(OpenMapRealMatrix matrix)
        {
            this.rows = matrix.rows;
            this.columns = matrix.columns;
            this.entries = new OpenIntToDoubleHashMap(matrix.entries);
        }

        /// <inheritdoc/>
        public override RealMatrix copy()
        {
            return new OpenMapRealMatrix(this);
        }

        /// <inheritdoc/>
        /// <exception cref="NumberIsTooLargeException"> if the total number of entries of the
        /// matrix is larger than <c>Int32.MaxValue</c>.</exception>
        public override RealMatrix createMatrix(int rowDimension, int columnDimension)
        {
            return new OpenMapRealMatrix(rowDimension, columnDimension);
        }

        /// <inheritdoc/>
        public override int getColumnDimension()
        {
            return columns;
        }

        /// <summary>
        /// Compute the sum of this matrix and <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to be added.</param>
        /// <returns><c>this</c> + <c>m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        public OpenMapRealMatrix add(OpenMapRealMatrix m)
        {

            MatrixUtils.checkAdditionCompatible(this, m);

            OpenMapRealMatrix outp = new OpenMapRealMatrix(this);
            for (OpenIntToDoubleHashMap.Iterator iterator = m.entries.iterator(); iterator.MoveNext(); )
            {
                int row = iterator.Current.Key / columns;
                int col = iterator.Current.Key - row * columns;
                outp.setEntry(row, col, getEntry(row, col) + iterator.Current.Value);
            }

            return outp;

        }

        /// <inheritdoc/>
        public new OpenMapRealMatrix subtract(RealMatrix m)
        {
            try
            {
                return subtract((OpenMapRealMatrix)m);
            }
            catch (InvalidCastException)
            {
                return (OpenMapRealMatrix)base.subtract(m);
            }
        }

        /// <summary>
        /// Subtract <c>m</c> from this matrix.
        /// </summary>
        /// <param name="m">Matrix to be subtracted.</param>
        /// <returns><c>this</c> - <c>m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        public OpenMapRealMatrix subtract(OpenMapRealMatrix m)
        {
            MatrixUtils.checkAdditionCompatible(this, m);

            OpenMapRealMatrix outp = new OpenMapRealMatrix(this);
            for (OpenIntToDoubleHashMap.Iterator iterator = m.entries.iterator(); iterator.MoveNext(); )
            {
                int row = iterator.Current.Key / columns;
                int col = iterator.Current.Key - row * columns;
                outp.setEntry(row, col, getEntry(row, col) - iterator.Current.Value);
            }

            return outp;
        }

        /// <inheritdoc/>
        /// <exception cref="NumberIsTooLargeException"> if <c>m</c> is an
        /// <c>OpenMapRealMatrix</c>, and the total number of entries of the product
        /// is larger than <c>Int32.MaxValue</c>.</exception>
        public new RealMatrix multiply(RealMatrix m)
        {
            try
            {
                return multiply((OpenMapRealMatrix)m);
            }
            catch (InvalidCastException)
            {

                MatrixUtils.checkMultiplicationCompatible(this, m);

                int outCols = m.getColumnDimension();
                BlockRealMatrix outp = new BlockRealMatrix(rows, outCols);
                for (OpenIntToDoubleHashMap.Iterator iterator = entries.iterator(); iterator.MoveNext(); )
                {
                    double value = iterator.Current.Value;
                    int key = iterator.Current.Key;
                    int i = key / columns;
                    int k = key % columns;
                    for (int j = 0; j < outCols; ++j)
                    {
                        outp.addToEntry(i, j, value * m.getEntry(k, j));
                    }
                }
                return outp;
            }
        }

        /// <summary>
        /// Postmultiply this matrix by <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to postmultiply by.</param>
        /// <returns><c>this</c> * <c>m</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the number of rows of <c>m</c>
        /// differ from the number of columns of <c>this</c> matrix.</exception>
        /// <exception cref="NumberIsTooLargeException"> if the total number of entries of the
        /// product is larger than <c>Int32.MaxValue</c>.</exception>
        public OpenMapRealMatrix multiply(OpenMapRealMatrix m)
        {
            // Safety check.
            MatrixUtils.checkMultiplicationCompatible(this, m);

            int outCols = m.getColumnDimension();
            OpenMapRealMatrix outp = new OpenMapRealMatrix(rows, outCols);
            for (OpenIntToDoubleHashMap.Iterator iterator = entries.iterator(); iterator.MoveNext(); )
            {
                double value = iterator.Current.Value;
                int key = iterator.Current.Key;
                int i = key / columns;
                int k = key % columns;
                for (int j = 0; j < outCols; ++j)
                {
                    int rightKey = m.computeKey(k, j);
                    if (m.entries.containsKey(rightKey))
                    {
                        int outKey = outp.computeKey(i, j);
                        double outValue =
                            outp.entries.get(outKey) + value * m.entries.get(rightKey);
                        if (outValue == 0.0)
                        {
                            outp.entries.remove(outKey);
                        }
                        else
                        {
                            outp.entries.put(outKey, outValue);
                        }
                    }
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public override double getEntry(int row, int column)
        {
            MatrixUtils.checkRowIndex(this, row);
            MatrixUtils.checkColumnIndex(this, column);
            return entries.get(computeKey(row, column));
        }

        /// <inheritdoc/>
        public override int getRowDimension()
        {
            return rows;
        }

        /// <inheritdoc/>
        public override void setEntry(int row, int column, double value)
        {
            MatrixUtils.checkRowIndex(this, row);
            MatrixUtils.checkColumnIndex(this, column);
            if (value == 0.0)
            {
                entries.remove(computeKey(row, column));
            }
            else
            {
                entries.put(computeKey(row, column), value);
            }
        }

        /// <inheritdoc/>
        public new void addToEntry(int row, int column, double increment)
        {
            MatrixUtils.checkRowIndex(this, row);
            MatrixUtils.checkColumnIndex(this, column);
            int key = computeKey(row, column);
            double value = entries.get(key) + increment;
            if (value == 0.0)
            {
                entries.remove(key);
            }
            else
            {
                entries.put(key, value);
            }
        }

        /// <inheritdoc/>
        public new void multiplyEntry(int row, int column, double factor)
        {
            MatrixUtils.checkRowIndex(this, row);
            MatrixUtils.checkColumnIndex(this, column);
            int key = computeKey(row, column);
            double value = entries.get(key) * factor;
            if (value == 0.0)
            {
                entries.remove(key);
            }
            else
            {
                entries.put(key, value);
            }
        }

        /// <summary>
        /// Compute the key to access a matrix element
        /// </summary>
        /// <param name="row">row index of the matrix element</param>
        /// <param name="column">column index of the matrix element</param>
        /// <returns>key within the map to access the matrix element</returns>
        private int computeKey(int row, int column)
        {
            return row * columns + column;
        }
    }
}