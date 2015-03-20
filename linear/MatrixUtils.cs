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
using Math3.fraction;
using Math3.util;
using System;

namespace Math3.linear
{
    /// <summary>
    /// A collection of static methods that operate on or return matrices.
    /// </summary>
    public class MatrixUtils
    {
        /// <summary>
        /// The default format for <see cref="RealMatrix"/> objects.
        /// </summary>
        public static readonly RealMatrixFormat DEFAULT_FORMAT = RealMatrixFormat.getInstance();

        /// <summary>
        /// A format for <see cref="RealMatrix"/> objects compatible with octave.
        /// </summary>
        public static readonly RealMatrixFormat OCTAVE_FORMAT = new RealMatrixFormat("[", "]", "", "", "; ", ", ");

        /// <summary>
        /// Private constructor.
        /// </summary>
        private MatrixUtils() : base() { }

        /// <summary>
        /// Returns a <see cref="RealMatrix"/> with specified dimensions.
        /// <para>The type of matrix returned depends on the dimension. Below
        /// 2^12 elements (i.e. 4096 elements or 64&times;64 for a
        /// square matrix) which can be stored in a 32kB array, a
        /// <see cref="Array2DRowRealMatrix"/> instance is built. Above this threshold a
        /// <see cref="BlockRealMatrix"/> instance is built.</para>
        /// <para>The matrix elements are all set to 0.0.</para>
        /// </summary>
        /// <param name="rows">number of rows of the matrix</param>
        /// <param name="columns">number of columns of the matrix</param>
        /// <returns>RealMatrix with specified dimensions</returns>
        /// <remarks>
        /// See <see cref="createRealMatrix(double[][])"/>
        /// </remarks>
        public static RealMatrix createRealMatrix(int rows, int columns)
        {
            if (rows * columns <= 4096)
            {
                return new Array2DRowRealMatrix(rows, columns);
            }
            else
            {
                return new BlockRealMatrix(rows, columns);
            }
        }

        /// <summary>
        /// Returns a <see cref="FieldMatrix"/> with specified dimensions.
        /// <para>The type of matrix returned depends on the dimension. Below
        /// 2^12 elements (i.e. 4096 elements or 64&times;64 for a
        /// square matrix), a <see cref="FieldMatrix"/> instance is built. Above
        /// this threshold a <see cref="BlockFieldMatrix"/> instance is built.</para>
        /// <para>The matrix elements are all set to field.getZero().</para>
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="field">field to which the matrix elements belong</param>
        /// <param name="rows">number of rows of the matrix</param>
        /// <param name="columns">number of columns of the matrix</param>
        /// <returns>FieldMatrix with specified dimensions</returns>
        /// <remarks>
        /// See <see cref="createFieldMatrix(FieldElement[][])"/>
        /// </remarks>
        public static FieldMatrix<T> createFieldMatrix<T>(Field<T> field, int rows, int columns) where T : FieldElement<T>
        {
            if (rows * columns <= 4096)
            {
                return new Array2DRowFieldMatrix<T>(field, rows, columns);
            }
            else
            {
                return new BlockFieldMatrix<T>(field, rows, columns);
            }
        }

        /// <summary>
        /// Returns a <see cref="RealMatrix"/> whose entries are the the values in the
        /// the input array.
        /// <para>The type of matrix returned depends on the dimension. Below
        /// 2^12 elements (i.e. 4096 elements or 64&times;64 for a
        /// square matrix) which can be stored in a 32kB array, a 
        /// <see cref="Array2DRowRealMatrix"/> instance is built. Above this threshold a 
        /// <see cref="BlockRealMatrix"/> instance is built.</para>
        /// <para>The input array is copied, not referenced.</para>
        /// </summary>
        /// <param name="data">input array</param>
        /// <returns>RealMatrix containing the values of the array</returns>
        /// <exception cref="DimensionMismatchException">
        /// if <c>data</c> is not rectangular (not all rows have the same length).</exception>
        /// <exception cref="NoDataException"> if a row or column is empty.</exception>
        /// <exception cref="NullArgumentException"> if either <c>data</c> or <c>data[0]</c>
        /// is <c>null</c>.</exception>
        /// <exception cref="DimensionMismatchException"> if <c>data</c> is not rectangular.
        /// </exception>
        /// <remarks>
        /// See <see cref="createRealMatrix(int, int)"/>
        /// </remarks>
        public static RealMatrix createRealMatrix(double[][] data)
        {
            if (data == null ||
                data[0] == null)
            {
                throw new NullArgumentException();
            }
            if (data.Length * data[0].Length <= 4096)
            {
                return new Array2DRowRealMatrix(data);
            }
            else
            {
                return new BlockRealMatrix(data);
            }
        }

        /// <summary>
        /// Returns a <see cref="FieldMatrix"/> whose entries are the the values in the
        /// the input array.
        /// <para>The type of matrix returned depends on the dimension. Below
        /// 2^12 elements (i.e. 4096 elements or 64&times;64 for a
        /// square matrix), a <see cref="FieldMatrix"/> instance is built. Above
        /// this threshold a <see cref="BlockFieldMatrix"/> instance is built.</para>
        /// <para>The input array is copied, not referenced.</para>
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="data">input array</param>
        /// <returns>a matrix containing the values of the array.</returns>
        /// <exception cref="DimensionMismatchException">
        /// if <c>data</c> is not rectangular (not all rows have the same length).</exception>
        /// <exception cref="NoDataException"> if a row or column is empty.</exception>
        /// <exception cref="NullArgumentException"> if either <c>data</c> or <c>data[0]</c>
        /// is <c>null</c>.</exception>
        /// <remarks>
        /// See <see cref="createFieldMatrix(Field, int, int)"/>
        /// </remarks>
        public static FieldMatrix<T> createFieldMatrix<T>(T[][] data) where T : FieldElement<T>
        {
            if (data == null ||
                data[0] == null)
            {
                throw new NullArgumentException();
            }
            if (data.Length * data[0].Length <= 4096)
            {
                return new Array2DRowFieldMatrix<T>(data);
            }
            else
            {
                return new BlockFieldMatrix<T>(data);
            }
        }

        /// <summary>
        /// Returns <c>dimension x dimension</c> identity matrix.
        /// </summary>
        /// <param name="dimension">dimension of identity matrix to generate</param>
        /// <returns>identity matrix</returns>
        /// <exception cref="IllegalArgumentException"> if dimension is not positive</exception>
        public static RealMatrix createRealIdentityMatrix(int dimension)
        {
            RealMatrix m = createRealMatrix(dimension, dimension);
            for (int i = 0; i < dimension; ++i)
            {
                m.setEntry(i, i, 1.0);
            }
            return m;
        }

        /// <summary>
        /// Returns <c>dimension x dimension</c> identity matrix.
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="field">field to which the elements belong</param>
        /// <param name="dimension">dimension of identity matrix to generate</param>
        /// <returns>identity matrix</returns>
        /// <exception cref="IllegalArgumentException"> if dimension is not positive</exception>
        public static FieldMatrix<T> createFieldIdentityMatrix<T>(Field<T> field, int dimension) where T : FieldElement<T>
        {
            T zero = field.getZero();
            T one = field.getOne();
            T[][] d = MathArrays.buildArray(field, dimension, dimension);
            for (int row = 0; row < dimension; row++)
            {
                T[] dRow = d[row];
                for (int i = 0; i < row; ++i)
                {
                    dRow[i] = zero;
                }
                dRow[row] = one;
            }
            return new Array2DRowFieldMatrix<T>(field, d, false);
        }

        /// <summary>
        /// Returns a diagonal matrix with specified elements.
        /// </summary>
        /// <param name="diagonal">diagonal elements of the matrix (the array elements
        /// will be copied)</param>
        /// <returns>diagonal matrix</returns>
        public static RealMatrix createRealDiagonalMatrix(double[] diagonal)
        {
            RealMatrix m = createRealMatrix(diagonal.Length, diagonal.Length);
            for (int i = 0; i < diagonal.Length; ++i)
            {
                m.setEntry(i, i, diagonal[i]);
            }
            return m;
        }

        /// <summary>
        /// Returns a diagonal matrix with specified elements.
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="diagonal">diagonal elements of the matrix (the array elements
        /// will be copied)</param>
        /// <returns>diagonal matrix</returns>
        public static FieldMatrix<T> createFieldDiagonalMatrix<T>(T[] diagonal) where T : FieldElement<T>
        {
            FieldMatrix<T> m =
                createFieldMatrix(diagonal[0].getField(), diagonal.Length, diagonal.Length);
            for (int i = 0; i < diagonal.Length; ++i)
            {
                m.setEntry(i, i, diagonal[i]);
            }
            return m;
        }

        /// <summary>
        /// Creates a <see cref="RealVector"/> using the data from the input array.
        /// </summary>
        /// <param name="data">the input data</param>
        /// <returns>a data.length RealVector</returns>
        /// <exception cref="NoDataException"> if <c>data</c> is empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>data</c> is <c>null</c>.</exception>
        public static RealVector createRealVector(double[] data)
        {
            if (data == null)
            {
                throw new NullArgumentException();
            }
            return new ArrayRealVector(data, true);
        }

        /// <summary>
        /// Creates a <see cref="FieldVector"/> using the data from the input array.
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="data">the input data</param>
        /// <returns>a data.length FieldVector</returns>
        /// <exception cref="NoDataException"> if <c>data</c> is empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>data</c> is <c>null</c>.</exception>
        /// <exception cref="ZeroException"> if <c>data</c> has 0 elements</exception>
        public static FieldVector<T> createFieldVector<T>(T[] data) where T : FieldElement<T>
        {
            if (data == null)
            {
                throw new NullArgumentException();
            }
            if (data.Length == 0)
            {
                throw new ZeroException(new LocalizedFormats("VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT"));
            }
            return new ArrayFieldVector<T>(data[0].getField(), data, true);
        }

        /// <summary>
        /// Create a row <see cref="RealMatrix"/> using the data from the input
        /// array.
        /// </summary>
        /// <param name="rowData">the input row data</param>
        /// <returns>a 1 x rowData.length RealMatrix</returns>
        /// <exception cref="NoDataException"> if <c>rowData</c> is empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>rowData</c> is <c>null</c>.</exception>
        public static RealMatrix createRowRealMatrix(double[] rowData)
        {
            if (rowData == null)
            {
                throw new NullArgumentException();
            }
            int nCols = rowData.Length;
            RealMatrix m = createRealMatrix(1, nCols);
            for (int i = 0; i < nCols; ++i)
            {
                m.setEntry(0, i, rowData[i]);
            }
            return m;
        }

        /// <summary>
        /// Create a row <see cref="FieldMatrix"/> using the data from the input
        /// array.
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="rowData">the input row data</param>
        /// <returns>a 1 x rowData.length FieldMatrix</returns>
        /// <exception cref="NoDataException"> if <c>rowData</c> is empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>rowData</c> is <c>null</c>.</exception>
        public static FieldMatrix<T> createRowFieldMatrix<T>(T[] rowData) where T : FieldElement<T>
        {
            if (rowData == null)
            {
                throw new NullArgumentException();
            }
            int nCols = rowData.Length;
            if (nCols == 0)
            {
                throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_COLUMN"));
            }
            FieldMatrix<T> m = createFieldMatrix(rowData[0].getField(), 1, nCols);
            for (int i = 0; i < nCols; ++i)
            {
                m.setEntry(0, i, rowData[i]);
            }
            return m;
        }

        /// <summary>
        /// Creates a column <see cref="RealMatrix"/> using the data from the input
        /// array.
        /// </summary>
        /// <param name="columnData">the input column data</param>
        /// <returns>a columnData x 1 RealMatrix</returns>
        /// <exception cref="NoDataException"> if <c>columnData</c> is empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>columnData</c> is <c>null</c>.
        /// </exception>
        public static RealMatrix createColumnRealMatrix(double[] columnData)
        {
            if (columnData == null)
            {
                throw new NullArgumentException();
            }
            int nRows = columnData.Length;
            RealMatrix m = createRealMatrix(nRows, 1);
            for (int i = 0; i < nRows; ++i)
            {
                m.setEntry(i, 0, columnData[i]);
            }
            return m;
        }

        /// <summary>
        /// Creates a column <see cref="FieldMatrix"/> using the data from the input
        /// array.
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="columnData">the input column data</param>
        /// <returns>a columnData x 1 FieldMatrix</returns>
        /// <exception cref="NoDataException"> if <c>data</c> is empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>columnData</c> is <c>null</c>.
        /// </exception>
        public static FieldMatrix<T> createColumnFieldMatrix<T>(T[] columnData) where T : FieldElement<T>
        {
            if (columnData == null)
            {
                throw new NullArgumentException();
            }
            int nRows = columnData.Length;
            if (nRows == 0)
            {
                throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_ROW"));
            }
            FieldMatrix<T> m = createFieldMatrix(columnData[0].getField(), nRows, 1);
            for (int i = 0; i < nRows; ++i)
            {
                m.setEntry(i, 0, columnData[i]);
            }
            return m;
        }

        /// <summary>
        /// Checks whether a matrix is symmetric, within a given relative tolerance.
        /// </summary>
        /// <param name="matrix">Matrix to check.</param>
        /// <param name="relativeTolerance">Tolerance of the symmetry check.</param>
        /// <param name="raiseException">If <c>true</c>, an exception will be raised if
        /// the matrix is not symmetric.</param>
        /// <returns><c>true</c> if <c>matrix</c> is symmetric.</returns>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        /// <exception cref="NonSymmetricMatrixException"> if the matrix is not symmetric.
        /// </exception>
        private static Boolean isSymmetricInternal(RealMatrix matrix, double relativeTolerance, Boolean raiseException)
        {
            int rows = matrix.getRowDimension();
            if (rows != matrix.getColumnDimension())
            {
                if (raiseException)
                {
                    throw new NonSquareMatrixException(rows, matrix.getColumnDimension());
                }
                else
                {
                    return false;
                }
            }
            for (int i = 0; i < rows; i++)
            {
                for (int j = i + 1; j < rows; j++)
                {
                    double mij = matrix.getEntry(i, j);
                    double mji = matrix.getEntry(j, i);
                    if (FastMath.abs(mij - mji) >
                        FastMath.max(FastMath.abs(mij), FastMath.abs(mji)) * relativeTolerance)
                    {
                        if (raiseException)
                        {
                            throw new NonSymmetricMatrixException(i, j, relativeTolerance);
                        }
                        else
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// Checks whether a matrix is symmetric.
        /// </summary>
        /// <param name="matrix">Matrix to check.</param>
        /// <param name="eps">Relative tolerance.</param>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        /// <exception cref="NonSymmetricMatrixException"> if the matrix is not symmetric.
        /// </exception>
        public static void checkSymmetric(RealMatrix matrix, double eps)
        {
            isSymmetricInternal(matrix, eps, true);
        }

        /// <summary>
        /// Checks whether a matrix is symmetric.
        /// </summary>
        /// <param name="matrix">Matrix to check.</param>
        /// <param name="eps">Relative tolerance.</param>
        /// <returns><c>true</c> if <c>matrix</c> is symmetric.</returns>
        public static Boolean isSymmetric(RealMatrix matrix, double eps)
        {
            return isSymmetricInternal(matrix, eps, false);
        }

        /// <summary>
        /// Check if matrix indices are valid.
        /// </summary>
        /// <param name="m">Matrix.</param>
        /// <param name="row">Row index to check.</param>
        /// <param name="column"Column index to check.></param>
        /// <exception cref="OutOfRangeException"> if <c>row</c> or <c>column</c> is not
        /// a valid index.</exception>
        public static void checkMatrixIndex(AnyMatrix m, int row, int column)
        {
            checkRowIndex(m, row);
            checkColumnIndex(m, column);
        }

        /// <summary>
        /// Check if a row index is valid.
        /// </summary>
        /// <param name="m">Matrix.</param>
        /// <param name="row">Row index to check.</param>
        /// <exception cref="OutOfRangeException"> if <c>row</c> is not a valid index.</exception>
        public static void checkRowIndex(AnyMatrix m, int row)
        {
            if (row < 0 ||
                row >= m.getRowDimension())
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("ROW_INDEX"), row, 0, m.getRowDimension() - 1);
            }
        }

        /// <summary>
        /// Check if a column index is valid.
        /// </summary>
        /// <param name="m">Matrix.</param>
        /// <param name="column">Column index to check.</param>
        /// <exception cref="OutOfRangeException"> if <c>column</c> is not a valid index.
        /// </exception>
        public static void checkColumnIndex(AnyMatrix m, int column)
        {
            if (column < 0 || column >= m.getColumnDimension())
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("COLUMN_INDEX"), column, 0, m.getColumnDimension() - 1);
            }
        }

        /// <summary>
        /// Check if submatrix ranges indices are valid.
        /// Rows and columns are indicated counting from 0 to <c>n - 1</c>.
        /// </summary>
        /// <param name="m">Matrix.</param>
        /// <param name="startRow">Initial row index.</param>
        /// <param name="endRow">Final row index.</param>
        /// <param name="startColumn">Initial column index.</param>
        /// <param name="endColumn">Final column index.</param>
        /// <exception cref="OutOfRangeException"> if the indices are invalid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        public static void checkSubMatrixIndex(AnyMatrix m, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkRowIndex(m, startRow);
            checkRowIndex(m, endRow);
            if (endRow < startRow)
            {
                throw new NumberIsTooSmallException<Int32, Int32>(new LocalizedFormats("INITIAL_ROW_AFTER_FINAL_ROW"), endRow, startRow, false);
            }

            checkColumnIndex(m, startColumn);
            checkColumnIndex(m, endColumn);
            if (endColumn < startColumn)
            {
                throw new NumberIsTooSmallException<Int32, Int32>(new LocalizedFormats("INITIAL_COLUMN_AFTER_FINAL_COLUMN"), endColumn, startColumn, false);
            }
        }

        /// <summary>
        /// Check if submatrix ranges indices are valid.
        /// Rows and columns are indicated counting from 0 to n-1.
        /// </summary>
        /// <param name="m">Matrix.</param>
        /// <param name="selectedRows">Array of row indices.</param>
        /// <param name="selectedColumns">Array of column indices.</param>
        /// <exception cref="NullArgumentException"> if <c>selectedRows</c> or
        /// <c>selectedColumns</c> are <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if the row or column selections are empty (zero
        /// length).</exception>
        /// <exception cref="OutOfRangeException"> if row or column selections are not valid.
        /// </exception>
        public static void checkSubMatrixIndex(AnyMatrix m, int[] selectedRows, int[] selectedColumns)
        {
            if (selectedRows == null)
            {
                throw new NullArgumentException();
            }
            if (selectedColumns == null)
            {
                throw new NullArgumentException();
            }
            if (selectedRows.Length == 0)
            {
                throw new NoDataException(new LocalizedFormats("EMPTY_SELECTED_ROW_INDEX_ARRAY"));
            }
            if (selectedColumns.Length == 0)
            {
                throw new NoDataException(new LocalizedFormats("EMPTY_SELECTED_COLUMN_INDEX_ARRAY"));
            }

            foreach (int row in selectedRows)
            {
                checkRowIndex(m, row);
            }
            foreach (int column in selectedColumns)
            {
                checkColumnIndex(m, column);
            }
        }

        /// <summary>
        /// Check if matrices are addition compatible.
        /// </summary>
        /// <param name="left">Left hand side matrix.</param>
        /// <param name="right">Right hand side matrix.</param>
        /// <exception cref="MatrixDimensionMismatchException"> if the matrices are not addition
        /// compatible.</exception>
        public static void checkAdditionCompatible(AnyMatrix left, AnyMatrix right)
        {
            if ((left.getRowDimension() != right.getRowDimension()) ||
                (left.getColumnDimension() != right.getColumnDimension()))
            {
                throw new MatrixDimensionMismatchException(left.getRowDimension(), left.getColumnDimension(),
                                                           right.getRowDimension(), right.getColumnDimension());
            }
        }

        /// <summary>
        /// Check if matrices are subtraction compatible
        /// </summary>
        /// <param name="left">Left hand side matrix.</param>
        /// <param name="right">Right hand side matrix.</param>
        /// <exception cref="MatrixDimensionMismatchException"> if the matrices are not addition
        /// compatible.</exception>
        public static void checkSubtractionCompatible(AnyMatrix left, AnyMatrix right)
        {
            if ((left.getRowDimension() != right.getRowDimension()) ||
                (left.getColumnDimension() != right.getColumnDimension()))
            {
                throw new MatrixDimensionMismatchException(left.getRowDimension(), left.getColumnDimension(),
                                                           right.getRowDimension(), right.getColumnDimension());
            }
        }

        /// <summary>
        /// Check if matrices are multiplication compatible
        /// </summary>
        /// <param name="left">Left hand side matrix.</param>
        /// <param name="right">Right hand side matrix.</param>
        /// <exception cref="DimensionMismatchException"> if matrices are not multiplication
        /// compatible.</exception>
        public static void checkMultiplicationCompatible(AnyMatrix left, AnyMatrix right)
        {

            if (left.getColumnDimension() != right.getRowDimension())
            {
                throw new DimensionMismatchException(left.getColumnDimension(),
                                                     right.getRowDimension());
            }
        }

        /// <summary>
        /// Convert a <see cref="FieldMatrix"/>/<see cref="Fraction"/> matrix to a 
        /// <see cref="RealMatrix"/>.
        /// </summary>
        /// <param name="m">Matrix to convert.</param>
        /// <returns>the converted matrix.</returns>
        public static Array2DRowRealMatrix fractionMatrixToRealMatrix(FieldMatrix<Fraction> m)
        {
            FractionMatrixConverter converter = new FractionMatrixConverter();
            m.walkInOptimizedOrder(converter);
            return converter.getConvertedMatrix();
        }

        /// <summary>
        /// Converter for <see cref="FieldMatrix"/>/<see cref="Fraction"/>.
        /// </summary>
        private class FractionMatrixConverter : DefaultFieldMatrixPreservingVisitor<Fraction>
        {
            /// <summary>
            /// Converted array.
            /// </summary>
            private double[][] data;

            /// <summary>
            /// Simple constructor.
            /// </summary>
            public FractionMatrixConverter() : base(Fraction.ZERO) { }

            /// <inheritdoc/>
            public new void start(int rows, int columns,
                              int startRow, int endRow, int startColumn, int endColumn)
            {
                data = new double[rows][];
            }

            /// <inheritdoc/>
            public new void visit(int row, int column, Fraction value)
            {
                data[row][column] = value.doubleValue();
            }

            /// <summary>
            /// Get the converted matrix.
            /// </summary>
            /// <returns>the converted matrix.</returns>
            internal Array2DRowRealMatrix getConvertedMatrix()
            {
                return new Array2DRowRealMatrix(data, false);
            }

        }

        /// <summary>
        /// Convert a <see cref="FieldMatrix"/>/<see cref="BigFraction"/> matrix to a 
        /// <see cref="RealMatrix"/>.
        /// </summary>
        /// <param name="m">Matrix to convert.</param>
        /// <returns>the converted matrix.</returns>
        public static Array2DRowRealMatrix bigFractionMatrixToRealMatrix(FieldMatrix<BigFraction> m)
        {
            BigFractionMatrixConverter converter = new BigFractionMatrixConverter();
            m.walkInOptimizedOrder(converter);
            return converter.getConvertedMatrix();
        }

        /// <summary>
        /// Converter for <see cref="FieldMatrix"/>/<see cref="BigFraction"/>.
        /// </summary>
        private class BigFractionMatrixConverter : DefaultFieldMatrixPreservingVisitor<BigFraction>
        {
            /// <summary>
            /// Converted array.
            /// </summary>
            private double[][] data;

            /// <summary>
            /// Simple constructor.
            /// </summary>
            public BigFractionMatrixConverter() : base(BigFraction.ZERO) { }

            /// <inheritdoc/>
            public new void start(int rows, int columns,
                              int startRow, int endRow, int startColumn, int endColumn)
            {
                data = new double[rows][];
            }

            /// <inheritdoc/>
            public new void visit(int row, int column, BigFraction value)
            {
                data[row][column] = value.doubleValue();
            }

            /// <summary>
            /// Get the converted matrix.
            /// </summary>
            /// <returns>the converted matrix.</returns>
            internal Array2DRowRealMatrix getConvertedMatrix()
            {
                return new Array2DRowRealMatrix(data, false);
            }
        }

        /// <summary>
        /// Solve  a  system of composed of a Lower Triangular Matrix
        /// <see cref="RealMatrix"/>.
        /// <para>
        /// This method is called to solve systems of equations which are
        /// of the lower triangular form. The matrix <see cref="RealMatrix"/>
        /// is assumed, though not checked, to be in lower triangular form.
        /// The vector <see cref="RealVector"/> is overwritten with the solution.
        /// The matrix is checked that it is square and its dimensions match
        /// the length of the vector.
        /// </para>
        /// </summary>
        /// <param name="rm">RealMatrix which is lower triangular</param>
        /// <param name="b">RealVector this is overwritten</param>
        /// <exception cref="DimensionMismatchException"> if the matrix and vector are not
        /// conformable</exception>
        /// <exception cref="NonSquareMatrixException"> if the matrix <c>rm</c> is not square
        /// </exception>
        /// <exception cref="MathArithmeticException"> if the absolute value of one of the diagonal
        /// coefficient of <c>rm</c> is lower than <see cref="Precision.SAFE_MIN"/></exception>
        public static void solveLowerTriangularSystem(RealMatrix rm, RealVector b)
        {
            if ((rm == null) || (b == null) || (rm.getRowDimension() != b.getDimension()))
            {
                throw new DimensionMismatchException(
                        (rm == null) ? 0 : rm.getRowDimension(),
                        (b == null) ? 0 : b.getDimension());
            }
            if (rm.getColumnDimension() != rm.getRowDimension())
            {
                throw new NonSquareMatrixException(rm.getRowDimension(),
                                                   rm.getColumnDimension());
            }
            int rows = rm.getRowDimension();
            for (int i = 0; i < rows; i++)
            {
                double diag = rm.getEntry(i, i);
                if (FastMath.abs(diag) < Precision.SAFE_MIN)
                {
                    throw new MathArithmeticException(new LocalizedFormats("ZERO_DENOMINATOR"));
                }
                double bi = b.getEntry(i) / diag;
                b.setEntry(i, bi);
                for (int j = i + 1; j < rows; j++)
                {
                    b.setEntry(j, b.getEntry(j) - bi * rm.getEntry(j, i));
                }
            }
        }

        /// <summary>
        /// Solver a  system composed  of an Upper Triangular Matrix
        /// <see cref="RealMatrix"/>.
        /// <para>
        /// This method is called to solve systems of equations which are
        /// of the lower triangular form. The matrix <see cref="RealMatrix"/>
        /// is assumed, though not checked, to be in upper triangular form.
        /// The vector <see cref="RealVector"/> is overwritten with the solution.
        /// The matrix is checked that it is square and its dimensions match
        /// the length of the vector.
        /// </para>
        /// </summary>
        /// <param name="rm">RealMatrix which is upper triangular</param>
        /// <param name="b">RealVector this is overwritten</param>
        /// <exception cref="DimensionMismatchException"> if the matrix and vector are not
        /// conformable</exception>
        /// <exception cref="NonSquareMatrixException"> if the matrix <c>rm</c> is not
        /// square</exception>
        /// <exception cref="MathArithmeticException"> if the absolute value of one of the diagonal
        /// coefficient of <c>rm</c> is lower than <see cref="Precision.SAFE_MIN"/></exception>
        public static void solveUpperTriangularSystem(RealMatrix rm, RealVector b)
        {
            if ((rm == null) || (b == null) || (rm.getRowDimension() != b.getDimension()))
            {
                throw new DimensionMismatchException(
                        (rm == null) ? 0 : rm.getRowDimension(),
                        (b == null) ? 0 : b.getDimension());
            }
            if (rm.getColumnDimension() != rm.getRowDimension())
            {
                throw new NonSquareMatrixException(rm.getRowDimension(),
                                                   rm.getColumnDimension());
            }
            int rows = rm.getRowDimension();
            for (int i = rows - 1; i > -1; i--)
            {
                double diag = rm.getEntry(i, i);
                if (FastMath.abs(diag) < Precision.SAFE_MIN)
                {
                    throw new MathArithmeticException(new LocalizedFormats("ZERO_DENOMINATOR"));
                }
                double bi = b.getEntry(i) / diag;
                b.setEntry(i, bi);
                for (int j = i - 1; j > -1; j--)
                {
                    b.setEntry(j, b.getEntry(j) - bi * rm.getEntry(j, i));
                }
            }
        }

        /// <summary>
        /// Computes the inverse of the given matrix by splitting it into
        /// 4 sub-matrices.
        /// </summary>
        /// <param name="m">Matrix whose inverse must be computed.</param>
        /// <param name="splitIndex">Index that determines the "split" line and
        /// column.
        /// The element corresponding to this index will part of the
        /// upper-left sub-matrix.</param>
        /// <returns>the inverse of <c>m</c>.</returns>
        /// <exception cref="NonSquareMatrixException"> if <c>m</c> is not square.</exception>
        public static RealMatrix blockInverse(RealMatrix m, int splitIndex)
        {
            int n = m.getRowDimension();
            if (m.getColumnDimension() != n)
            {
                throw new NonSquareMatrixException(m.getRowDimension(),
                                                   m.getColumnDimension());
            }

            int splitIndex1 = splitIndex + 1;

            RealMatrix a = m.getSubMatrix(0, splitIndex, 0, splitIndex);
            RealMatrix b = m.getSubMatrix(0, splitIndex, splitIndex1, n - 1);
            RealMatrix c = m.getSubMatrix(splitIndex1, n - 1, 0, splitIndex);
            RealMatrix d = m.getSubMatrix(splitIndex1, n - 1, splitIndex1, n - 1);

            SingularValueDecomposition aDec = new SingularValueDecomposition(a);
            DecompositionSolver aSolver = aDec.getSolver();
            if (!aSolver.isNonSingular())
            {
                throw new SingularMatrixException();
            }
            RealMatrix aInv = aSolver.getInverse();

            SingularValueDecomposition dDec = new SingularValueDecomposition(d);
            DecompositionSolver dSolver = dDec.getSolver();
            if (!dSolver.isNonSingular())
            {
                throw new SingularMatrixException();
            }
            RealMatrix dInv = dSolver.getInverse();

            RealMatrix tmp1 = a.subtract(b.multiply(dInv).multiply(c));
            SingularValueDecomposition tmp1Dec = new SingularValueDecomposition(tmp1);
            DecompositionSolver tmp1Solver = tmp1Dec.getSolver();
            if (!tmp1Solver.isNonSingular())
            {
                throw new SingularMatrixException();
            }
            RealMatrix result00 = tmp1Solver.getInverse();

            RealMatrix tmp2 = d.subtract(c.multiply(aInv).multiply(b));
            SingularValueDecomposition tmp2Dec = new SingularValueDecomposition(tmp2);
            DecompositionSolver tmp2Solver = tmp2Dec.getSolver();
            if (!tmp2Solver.isNonSingular())
            {
                throw new SingularMatrixException();
            }
            RealMatrix result11 = tmp2Solver.getInverse();

            RealMatrix result01 = aInv.multiply(b).multiply(result11).scalarMultiply(-1);
            RealMatrix result10 = dInv.multiply(c).multiply(result00).scalarMultiply(-1);

            RealMatrix result = new Array2DRowRealMatrix(n, n);
            result.setSubMatrix(result00.getData(), 0, 0);
            result.setSubMatrix(result01.getData(), 0, splitIndex1);
            result.setSubMatrix(result10.getData(), splitIndex1, 0);
            result.setSubMatrix(result11.getData(), splitIndex1, splitIndex1);

            return result;
        }

        /// <summary>
        /// Computes the inverse of the given matrix.
        /// <para>
        /// By default, the inverse of the matrix is computed using the QR-decomposition,
        /// unless a more efficient method can be determined for the input matrix.
        /// </para>
        /// <para>
        /// Note: this method will use a singularity threshold of 0,
        /// use <see cref="inverse(RealMatrix, double)"/> if a different threshold is needed.
        /// </para>
        /// </summary>
        /// <param name="matrix">Matrix whose inverse shall be computed</param>
        /// <returns>the inverse of <c>matrix</c></returns>
        /// <exception cref="NullArgumentException"> if <c>matrix</c> is <c>null</c></exception>
        /// <exception cref="SingularMatrixException"> if m is singular</exception>
        /// <exception cref="NonSquareMatrixException"> if matrix is not square</exception>
        public static RealMatrix inverse(RealMatrix matrix)
        {
            return inverse(matrix, 0);
        }

        /// <summary>
        /// Computes the inverse of the given matrix.
        /// <para>
        /// By default, the inverse of the matrix is computed using the QR-decomposition,
        /// unless a more efficient method can be determined for the input matrix.
        /// </para>
        /// </summary>
        /// <param name="matrix">Matrix whose inverse shall be computed</param>
        /// <param name="threshold">Singularity threshold</param>
        /// <returns>the inverse of <c>m</c></returns>
        /// <exception cref="NullArgumentException"> if <c>matrix</c> is <c>null</c></exception>
        /// <exception cref="SingularMatrixException"> if matrix is singular</exception>
        /// <exception cref="NonSquareMatrixException"> if matrix is not square</exception>
        public static RealMatrix inverse(RealMatrix matrix, double threshold)
        {
            MathUtils.checkNotNull(matrix);

            if (!matrix.isSquare())
            {
                throw new NonSquareMatrixException(matrix.getRowDimension(),
                                                   matrix.getColumnDimension());
            }

            if (matrix is DiagonalMatrix)
            {
                return ((DiagonalMatrix)matrix).inverse(threshold);
            }
            else
            {
                QRDecomposition decomposition = new QRDecomposition(matrix, threshold);
                return decomposition.getSolver().getInverse();
            }
        }
    }
}