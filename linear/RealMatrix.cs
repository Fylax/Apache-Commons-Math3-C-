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
namespace Math3.linear
{
    /// <summary>
    /// Interface defining a real-valued matrix with basic algebraic operations.
    /// <para>
    /// Matrix element indexing is 0-based -- e.g., <c>getEntry(0, 0)</c>
    /// returns the element in the first row, first column of the matrix.</para>
    /// </summary>
    public interface RealMatrix : AnyMatrix
    {
        /// <summary>
        /// Create a new RealMatrix of the same type as the instance with the
        /// supplied row and column dimensions.
        /// </summary>
        /// <param name="rowDimension">the number of rows in the new matrix</param>
        /// <param name="columnDimension">the number of columns in the new matrix</param>
        /// <returns>a new matrix of the same type as the instance</returns>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension is not
        /// positive.</exception>
        RealMatrix createMatrix(int rowDimension, int columnDimension);

        /// <summary>
        /// Returns a (deep) copy of this.
        /// </summary>
        /// <returns>matrix copy</returns>
        RealMatrix copy();

        /// <summary>
        /// Returns the sum of <c>this</c> and <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to be added</param>
        /// <returns><c>this + m</c></returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        RealMatrix add(RealMatrix m);

        /// <summary>
        /// Returns <c>this</c> minus <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to be subtracted</param>
        /// <returns><c>this - m</c></returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c>.</exception>
        RealMatrix subtract(RealMatrix m);

        /// <summary>
        /// Returns the result of adding <c>d</c> to each entry of <c>this</c>.
        /// </summary>
        /// <param name="d">value to be added to each entry</param>
        /// <returns><c>d + this</c></returns>
        RealMatrix scalarAdd(double d);

        /// <summary>
        /// Returns the result of multiplying each entry of <c>this</c> by
        /// <c>d</c>.
        /// </summary>
        /// <param name="d">value to multiply all entries by</param>
        /// <returns><c>d * this</c></returns>
        RealMatrix scalarMultiply(double d);

        /// <summary>
        /// Returns the result of postmultiplying <c>this</c> by <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to postmultiply by</param>
        /// <returns><c>this * m</c></returns>
        /// <exception cref="DimensionMismatchException"> if
        /// <c>columnDimension(this) != rowDimension(m)</c></exception>
        RealMatrix multiply(RealMatrix m);

        /// <summary>
        /// Returns the result of premultiplying <c>this</c> by <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to premultiply by</param>
        /// <returns><c>m * this</c></returns>
        /// <exception cref="DimensionMismatchException"> if
        /// <c>rowDimension(this) != columnDimension(m)</c></exception>
        RealMatrix preMultiply(RealMatrix m);

        /// <summary>
        /// Returns the result of multiplying <c>this</c> with itself <c>p</c>
        /// times. Depending on the underlying storage, instability for high powers
        /// might occur.
        /// </summary>
        /// <param name="p">raise <c>this</c> to power <c>p</c></param>
        /// <returns><c>this^p</c></returns>
        /// <exception cref="NotPositiveException"> if <c>p < 0</c></exception>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square</exception>
        RealMatrix power(int p);

        /// <summary>
        /// Returns matrix entries as a two-dimensional array.
        /// </summary>
        /// <returns>2-dimensional array of entries</returns>
        double[][] getData();

        /// <summary>
        /// Returns the <a href="http://mathworld.wolfram.com/MaximumAbsoluteRowSumNorm.html">
        /// maximum absolute row sum norm</a> of the matrix.
        /// </summary>
        /// <returns>norm</returns>
        double getNorm();

        /// <summary>
        /// Returns the <a href="http://mathworld.wolfram.com/FrobeniusNorm.html">
        /// Frobenius norm</a> of the matrix
        /// </summary>
        /// <returns>norm</returns>
        double getFrobeniusNorm();

        /// <summary>
        /// Gets a submatrix. Rows and columns are indicated
        /// counting from 0 to n-1.
        /// </summary>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index (inclusive)</param>
        /// <returns>The subMatrix containing the data of the
        /// specified rows and columns.</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        RealMatrix getSubMatrix(int startRow, int endRow, int startColumn,
                                int endColumn);

        /// <summary>
        /// Gets a submatrix. Rows and columns are indicated counting from 0 to n-1.
        /// </summary>
        /// <param name="selectedRows">Array of row indices.</param>
        /// <param name="selectedColumns">Array of column indices.</param>
        /// <returns>The subMatrix containing the data in the specified rows and
        /// columns</returns>
        /// <exception cref="NullArgumentException"> if the row or column selections are
        /// <c>null</c></exception>
        /// <exception cref="NoDataException"> if the row or column selections are empty (zero
        /// length).</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        RealMatrix getSubMatrix(int[] selectedRows, int[] selectedColumns);

        /// <summary>
        /// Copy a submatrix. Rows and columns are indicated counting from 0 to n-1.
        /// </summary>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index (inclusive)</param>
        /// <param name="destination">The arrays where the submatrix data should be copied
        /// (if larger than rows/columns counts, only the upper-left part will be
        /// used)</param>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the destination array is too
        /// small.</exception>
        void copySubMatrix(int startRow, int endRow, int startColumn, int endColumn, double[][] destination);

        /// <summary>
        /// Copy a submatrix. Rows and columns are indicated counting from 0 to n-1.
        /// </summary>
        /// <param name="selectedRows">Array of row indices.</param>
        /// <param name="selectedColumns">Array of column indices.</param>
        /// <param name="destination">The arrays where the submatrix data should be copied
        /// (if larger than rows/columns counts, only the upper-left part will be
        /// used)</param>
        /// <exception cref="NullArgumentException"> if the row or column selections are
        /// <c>null</c></exception>
        /// <exception cref="NoDataException"> if the row or column selections are empty (zero
        /// length).</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the destination array is too
        /// small.</exception>
        void copySubMatrix(int[] selectedRows, int[] selectedColumns, double[][] destination);

        /// <summary>
        /// Replace the submatrix starting at <c>row, column</c> using data in the
        /// input <c>subMatrix</c> array. Indexes are 0-based.
        /// <para>
        /// Example:<para/>
        /// Starting with <code>
        /// 1  2  3  4
        /// 5  6  7  8
        /// 9  0  1  2
        /// </code>
        /// and <c>subMatrix = {{3, 4} {5,6}}</c>, invoking
        /// <c>setSubMatrix(subMatrix,1,1))</c> will result in 
        /// <code>
        /// 1  2  3  4
        /// 5  3  4  8
        /// 9  5  6  2
        /// </code></para>
        /// </summary>
        /// <param name="subMatrix">array containing the submatrix replacement data</param>
        /// <param name="row">row coordinate of the top, left element to be replaced</param>
        /// <param name="column">column coordinate of the top, left element to be replaced</param>
        /// <exception cref="NoDataException"> if <c>subMatrix</c> is empty.</exception>
        /// <exception cref="OutOfRangeException"> if <c>subMatrix</c> does not fit into
        /// this matrix from element in <c>(row, column)</c>.</exception>
        /// <exception cref="DimensionMismatchException"> if <c>subMatrix</c> is not rectangular
        /// (not all rows have the same length) or empty.</exception>
        /// <exception cref="NullArgumentException"> if <c>subMatrix</c> is <c>null</c>.
        /// </exception>
        void setSubMatrix(double[][] subMatrix, int row, int column);

        /// <summary>
        /// Get the entries at the given row index as a row matrix.  Row indices start
        /// at 0.
        /// </summary>
        /// <param name="row">Row to be fetched.</param>
        /// <returns>Matrix.</returns>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.</exception>
        RealMatrix getRowMatrix(int row);

        /// <summary>
        /// Sets the specified <c>row</c> of <c>this</c> matrix to the entries of
        /// the specified row <c>matrix</c>. Row indices start at 0.
        /// </summary>
        /// <param name="row">Row to be set.</param>
        /// <param name="matrix">Row matrix to be copied (must have one row and the same
        /// number of columns as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the row dimension of the
        /// <c>matrix</c> is not <c>1</c>, or the column dimensions of <c>this</c>
        /// and <c>matrix</c> do not match.</exception>
        void setRowMatrix(int row, RealMatrix matrix);

        /// <summary>
        /// Get the entries at the given column index as a column matrix. Column
        /// indices start at 0.
        /// </summary>
        /// <param name="column">Column to be fetched.</param>
        /// <returns>Matrix.</returns>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        RealMatrix getColumnMatrix(int column);

        /// <summary>
        /// Sets the specified <c>column</c> of <c>this</c> matrix to the entries
        /// of the specified column <c>matrix</c>. Column indices start at 0.
        /// </summary>
        /// <param name="column">Column to be set.</param>
        /// <param name="matrix">Column matrix to be copied (must have one column and the
        /// same number of rows as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the column dimension of the
        /// <c>matrix</c> is not <c>1</c>, or the row dimensions of <c>this</c>
        /// and <c>matrix</c> do not match.</exception>
        void setColumnMatrix(int column, RealMatrix matrix);

        /// <summary>
        /// Returns the entries in row number <c>row</c> as a vector. Row indices
        /// start at 0.
        /// </summary>
        /// <param name="row">Row to be fetched.</param>
        /// <returns>a row vector.</returns>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        RealVector getRowVector(int row);

        /// <summary>
        /// Sets the specified <c>row</c> of <c>this</c> matrix to the entries of
        /// the specified <c>vector</c>. Row indices start at 0.
        /// </summary>
        /// <param name="row">Row to be set.</param>
        /// <param name="vector">row vector to be copied (must have the same number of
        /// column as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the <c>vector</c> dimension
        /// does not match the column dimension of <c>this</c> matrix.</exception>
        void setRowVector(int row, RealVector vector);

        /// <summary>
        /// Get the entries at the given column index as a vector. Column indices
        /// start at 0.
        /// </summary>
        /// <param name="column">Column to be fetched.</param>
        /// <returns>a column vector.</returns>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid
        /// </exception>
        RealVector getColumnVector(int column);

        /// <summary>
        /// Sets the specified <c>column</c> of <c>this</c> matrix to the entries
        /// of the specified <c>vector</c>. Column indices start at 0.
        /// </summary>
        /// <param name="column">Column to be set.</param>
        /// <param name="vector">column vector to be copied (must have the same number of
        /// rows as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the <c>vector</c> dimension
        /// does not match the row dimension of <c>this</c> matrix.</exception>
        void setColumnVector(int column, RealVector vector);

        /// <summary>
        /// Get the entries at the given row index. Row indices start at 0.
        /// </summary>
        /// <param name="row">Row to be fetched.</param>
        /// <returns>the array of entries in the row.</returns>
        /// <exception cref="OutOfRangeException"> if the specified row index is not valid.
        /// </exception>
        double[] getRow(int row);

        /// <summary>
        /// Sets the specified <c>row</c> of <c>this</c> matrix to the entries
        /// of the specified <c>array</c>. Row indices start at 0.
        /// </summary>
        /// <param name="row">Row to be set.</param>
        /// <param name="array">Row matrix to be copied (must have the same number of
        /// columns as the instance)</param>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the <c>array</c> length does
        /// not match the column dimension of <c>this</c> matrix.</exception>
        void setRow(int row, double[] array);

        /// <summary>
        /// Get the entries at the given column index as an array. Column indices
        /// start at 0.
        /// </summary>
        /// <param name="column">Column to be fetched.</param>
        /// <returns>the array of entries in the column.</returns>
        /// <exception cref="OutOfRangeException"> if the specified column index is not valid.
        /// </exception>
        double[] getColumn(int column);

        /// <summary>
        /// Sets the specified <c>column</c> of <c>this</c> matrix to the entries
        /// of the specified <c>array</c>. Column indices start at 0.
        /// </summary>
        /// <param name="column">Column to be set.</param>
        /// <param name="array">Column array to be copied (must have the same number of
        /// rows as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the <c>array</c> length does
        /// not match the row dimension of <c>this</c> matrix.</exception>
        void setColumn(int column, double[] array);

        /// <summary>
        /// Get the entry in the specified row and column. Row and column indices
        /// start at 0.
        /// </summary>
        /// <param name="row">Row index of entry to be fetched.</param>
        /// <param name="column">Column index of entry to be fetched.</param>
        /// <returns>the matrix entry at <c>(row, column)</c>.</returns>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.
        /// </exception>
        double getEntry(int row, int column);

        /// <summary>
        /// Set the entry in the specified row and column. Row and column indices
        /// start at 0.
        /// </summary>
        /// <param name="row">Row index of entry to be set.</param>
        /// <param name="column">Column index of entry to be set.</param>
        /// <param name="value">the new value of the entry.</param>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid
        /// </exception>
        void setEntry(int row, int column, double value);

        /// <summary>
        /// Adds (in place) the specified value to the specified entry of
        /// <c>this</c> matrix. Row and column indices start at 0.
        /// </summary>
        /// <param name="row">Row index of the entry to be modified.</param>
        /// <param name="column">Column index of the entry to be modified.</param>
        /// <param name="increment">value to add to the matrix entry.</param>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.
        /// </exception>
        void addToEntry(int row, int column, double increment);

        /// <summary>
        /// Multiplies (in place) the specified entry of <c>this</c> matrix by the
        /// specified value. Row and column indices start at 0.
        /// </summary>
        /// <param name="row">Row index of the entry to be modified.</param>
        /// <param name="column">Column index of the entry to be modified.</param>
        /// <param name="factor">Multiplication factor for the matrix entry.</param>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.
        /// </exception>
        void multiplyEntry(int row, int column, double factor);

        /// <summary>
        /// Returns the transpose of this matrix. 
        /// </summary>
        /// <returns>transpose matrix</returns>
        RealMatrix transpose();

        /// <summary>
        /// Returns the <a href="http://mathworld.wolfram.com/MatrixTrace.html">
        /// trace</a> of the matrix (the sum of the elements on the main diagonal).
        /// </summary>
        /// <returns>the trace.</returns>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        double getTrace();

        /// <summary>
        /// Returns the result of multiplying this by the vector <c>v</c>.
        /// </summary>
        /// <param name="v">the vector to operate on</param>
        /// <returns><c>this * v</c></returns>
        /// <exception cref="DimensionMismatchException"> if the length of <c>v</c> does not
        /// match the column dimension of <c>this</c>.</exception>
        double[] operate(double[] v);

        /// <summary>
        /// Returns the result of multiplying this by the vector <c>v</c>.
        /// </summary>
        /// <param name="v">the vector to operate on</param>
        /// <returns><c>this * v</c></returns>
        /// <exception cref="DimensionMismatchException"> if the dimension of <c>v</c> does not
        /// match the column dimension of <c>this</c>.</exception>
        RealVector operate(RealVector v);

        /// <summary>
        /// Returns the (row) vector result of premultiplying this by the vector <c>v</c>.
        /// </summary>
        /// <param name="v">the row vector to premultiply by</param>
        /// <returns><c>v * this</c></returns>
        /// <exception cref="DimensionMismatchException"> if the length of <c>v</c> does not
        /// match the row dimension of <c>this</c>.</exception>
        double[] preMultiply(double[] v);

        /// <summary>
        /// Returns the (row) vector result of premultiplying this by the vector <c>v</c>.
        /// </summary>
        /// <param name="v">the row vector to premultiply by</param>
        /// <returns><c>v * this</c></returns>
        /// <exception cref="DimensionMismatchException"> if the dimension of <c>v</c> does
        /// not match the row dimension of <c>this</c>.</exception>
        RealVector preMultiply(RealVector v);

        /// <summary>
        /// Visit (and possibly change) all matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor)"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/>
        /// <para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// <para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/>
        /// <para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// <para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/>
        /// <para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <returns>the value returned by <see cref="RealMatrixChangingVisitor.end()"/> at the end
        /// of the walk</returns>
        double walkInRowOrder(RealMatrixChangingVisitor visitor);

        /// <summary>
        /// Visit (but don't change) all matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <returns>the value returned by <see cref="RealMatrixPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        double walkInRowOrder(RealMatrixPreservingVisitor visitor);

        /// <summary>
        /// Visit (and possibly change) some matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index</param>
        /// <returns>the value returned by <see cref="RealMatrixChangingVisitor.end()"/> at the end
        /// of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        double walkInRowOrder(RealMatrixChangingVisitor visitor, int startRow,
            int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (but don't change) some matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index</param>
        /// <returns>the value returned by <see cref="RealMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        double walkInRowOrder(RealMatrixPreservingVisitor visitor, int startRow,
            int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (and possibly change) all matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <returns>the value returned by <see cref="RealMatrixChangingVisitor.end()"/> at the end
        /// of the walk</returns>
        double walkInColumnOrder(RealMatrixChangingVisitor visitor);

        /// <summary>
        /// Visit (but don't change) all matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <returns>the value returned by <see cref="RealMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        double walkInColumnOrder(RealMatrixPreservingVisitor visitor);

        /// <summary>
        /// Visit (and possibly change) some matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index</param>
        /// <returns>the value returned by <see cref="RealMatrixChangingVisitor.end()"/> at the 
        /// end of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        double walkInColumnOrder(RealMatrixChangingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (but don't change) some matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index</param>
        /// <returns>the value returned by <see cref="RealMatrixPreservingVisitor.end()"/> at 
        /// the end of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        double walkInColumnOrder(RealMatrixPreservingVisitor visitor, int startRow,
            int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (and possibly change) all matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <returns>the value returned by <see cref="RealMatrixChangingVisitor.end()"/> at the
        /// end of the walk</returns>
        double walkInOptimizedOrder(RealMatrixChangingVisitor visitor);

        /// <summary>
        /// Visit (but don't change) all matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <returns>the value returned by <see cref="RealMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        double walkInOptimizedOrder(RealMatrixPreservingVisitor visitor);

        /// <summary>
        /// Visit (and possibly change) some matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor, int, int, int, int"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index (inclusive)</param>
        /// <returns>the value returned by <see cref="RealMatrixChangingVisitor.end()"/> at the
        /// end of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        double walkInOptimizedOrder(RealMatrixChangingVisitor visitor,
            int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (but don't change) some matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInRowOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixChangingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInColumnOrder(RealMatrixPreservingVisitor, int, int, int, int"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixPreservingVisitor"/><para/>
        /// <see cref="walkInOptimizedOrder(RealMatrixChangingVisitor, int, int, int, int"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index (inclusive)</param>
        /// <returns>the value returned by <see cref="RealMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        double walkInOptimizedOrder(RealMatrixPreservingVisitor visitor, int startRow, int endRow, int startColumn, int endColumn);
    }
}