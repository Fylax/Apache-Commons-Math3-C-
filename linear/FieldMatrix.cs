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
    /// Interface defining field-valued matrix with basic algebraic operations.
    /// <para>
    /// Matrix element indexing is 0-based -- e.g., <c>getEntry(0, 0)</c>
    /// returns the element in the first row, first column of the matrix.</para>
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public interface FieldMatrix<T> : AnyMatrix where T : FieldElement<T>
    {
        /// <summary>
        /// Get the type of field elements of the matrix.
        /// </summary>
        /// <returns>the type of field elements of the matrix.</returns>
        Field<T> getField();

        /// <summary>
        /// Create a new FieldMatrix(T) of the same type as the instance with
        /// the supplied row and column dimensions.
        /// </summary>
        /// <param name="rowDimension">the number of rows in the new matrix</param>
        /// <param name="columnDimension">the number of columns in the new matrix</param>
        /// <returns>a new matrix of the same type as the instance</returns>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension 
        /// is not positive.</exception>
        FieldMatrix<T> createMatrix(int rowDimension, int columnDimension);

        /// <summary>
        /// Make a (deep) copy of this.
        /// </summary>
        /// <returns>a copy of this matrix.</returns>
        FieldMatrix<T> copy();

        /// <summary>
        /// Compute the sum of this and m.
        /// </summary>
        /// <param name="m">Matrix to be added.</param>
        /// <returns><c>this</c> + <c>m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not 
        /// the same size as <c>this</c> matrix.</exception>
        FieldMatrix<T> add(FieldMatrix<T> m);

        /// <summary>
        /// Subtract <c>m</c> from this matrix.
        /// </summary>
        /// <param name="m">Matrix to be subtracted.</param>
        /// <returns><c>this</c> - <c>m</c>.</returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not
        /// the same size as <c>this</c> matrix.</exception>
        FieldMatrix<T> subtract(FieldMatrix<T> m);

        /// <summary>
        /// Increment each entry of this matrix.
        /// </summary>
        /// <param name="d">Value to be added to each entry.</param>
        /// <returns><c>d</c> + <c>this</c>.</returns>
        FieldMatrix<T> scalarAdd(T d);

        /// <summary>
        /// Multiply each entry by <c>d</c>.
        /// </summary>
        /// <param name="d">Value to multiply all entries by.</param>
        /// <returns><c>d</c> * <c>this</c>.</returns>
        FieldMatrix<T> scalarMultiply(T d);

        /// <summary>
        /// Postmultiply this matrix by <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to postmultiply by.</param>
        /// <returns><c>this</c> * <c>m</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the number of columns of
        /// <c>this</c> matrix is not equal to the number of rows of matrix
        /// <c>m</c>.</exception>
        FieldMatrix<T> multiply(FieldMatrix<T> m);

        /// <summary>
        /// Premultiply this matrix by <c>m</c>.
        /// </summary>
        /// <param name="m">Matrix to premultiply by.</param>
        /// <returns><c>m</c> * <c>this</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the number of columns of <c>m</c>
        /// differs from the number of rows of <c>this</c> matrix.</exception>
        FieldMatrix<T> preMultiply(FieldMatrix<T> m);

        /// <summary>
        /// Returns the result multiplying this with itself <c>p</c> times.
        /// Depending on the type of the field elements, T, instability for high
        /// powers might occur.
        /// </summary>
        /// <param name="p">raise this to power p</param>
        /// <returns>this^p</returns>
        /// <exception cref="NotPositiveException"> if <c>p < 0</c></exception>
        /// <exception cref="NonSquareMatrixException"> if <c>this matrix</c> is not square</exception>
        FieldMatrix<T> power(int p);

        /// <summary>
        /// Returns matrix entries as a two-dimensional array.
        /// </summary>
        /// <returns>a 2-dimensional array of entries.</returns>
        T[][] getData();

        /// <summary>
        /// Get a submatrix. Rows and columns are indicated
        /// counting from 0 to n - 1.
        /// </summary>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index (inclusive)</param>
        /// <returns>the matrix containing the data of the specified rows and columns.</returns>
        /// <exception cref="NumberIsTooSmallException"> is <c>endRow < startRow</c> of
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        FieldMatrix<T> getSubMatrix(int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Get a submatrix. Rows and columns are indicated
        /// counting from 0 to n - 1.
        /// </summary>
        /// <param name="selectedRows">Array of row indices.</param>
        /// <param name="selectedColumns">Array of column indices.</param>
        /// <returns>the matrix containing the data in the
        /// rows and columns.</returns>
        /// <exception cref="NoDataException"> if <c>selectedRows</c> or
        /// <c>selectedColumns</c> is empty</exception>
        /// <exception cref="NullArgumentException"> if <c>selectedRows</c> or
        /// <c>selectedColumns</c> is <c>null</c>.</exception>
        /// <exception cref="OutOfRangeException"> if row or column selections are not valid.</exception>
        FieldMatrix<T> getSubMatrix(int[] selectedRows, int[] selectedColumns);

        /// <summary>
        /// Copy a submatrix. Rows and columns are indicated
        /// counting from 0 to n-1.
        /// </summary>
        /// <param name="startRow">Initial row index.</param>
        /// <param name="endRow">row index (inclusive).</param>
        /// <param name="startColumn">Initial column index.</param>
        /// <param name="endColumn">column index (inclusive).</param>
        /// <param name="destination">The arrays where the submatrix data should be copied
        /// (if larger than rows/columns counts, only the upper-left part will be used).</param>
        /// <exception cref="">MatrixDimensionMismatchException if the dimensions of
        /// <c>destination</c> do not match those of <c>this</c>.</exception>
        /// <exception cref="NumberIsTooSmallException"> is <c>endRow < startRow</c> of
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="IllegalArgumentException"> if the destination array is too small.</exception>
        void copySubMatrix(int startRow, int endRow, int startColumn, int endColumn, T[][] destination);

        /// <summary>
        /// Copy a submatrix. Rows and columns are indicated
        /// counting from 0 to n - 1.
        /// </summary>
        /// <param name="selectedRows">Array of row indices.</param>
        /// <param name="selectedColumns">Array of column indices.</param>
        /// <param name="destination">Arrays where the submatrix data should be copied
        /// (if larger than rows/columns counts, only the upper-left part will be used)</param>
        /// <exception cref="MatrixDimensionMismatchException"> if the dimensions of
        /// <c>destination</c> do not match those of <c>this</c>.</exception>
        /// <exception cref="NoDataException"> if <c>selectedRows</c> or
        /// <c>selectedColumns</c> is empty</exception>
        /// <exception cref="NullArgumentException"> if <c>selectedRows</c> or
        /// <c>selectedColumns</c> is <c>null</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        void copySubMatrix(int[] selectedRows, int[] selectedColumns, T[][] destination);

        /// <summary>
        /// Replace the submatrix starting at <c>(row, column)</c> using data in the
        /// input <c>subMatrix</c> array. Indexes are 0-based.
        /// <para>
        /// Example:<para/>
        /// Starting with
        /// <code>
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
        /// </code>
        /// </para>
        /// </summary>
        /// <param name="subMatrix">Array containing the submatrix replacement data.</param>
        /// <param name="row">Row coordinate of the top-left element to be replaced.</param>
        /// <param name="column">Column coordinate of the top-left element to be replaced.</param>
        /// <exception cref="OutOfRangeException"> if <c>subMatrix</c> does not fit into this
        /// matrix from element in <c>(row, column)</c>.</exception>
        /// <exception cref="NoDataException"> if a row or column of <c>subMatrix</c> is empty.
        /// </exception>
        /// <exception cref="DimensionMismatchException"> if <c>subMatrix</c> is not
        /// rectangular (not all rows have the same length).</exception>
        /// <exception cref="NullArgumentException"> if <c>subMatrix</c> is <c>null</c>.
        /// </exception>
        void setSubMatrix(T[][] subMatrix, int row, int column);

        /// <summary>
        /// Get the entries in row number <c>row</c>
        /// as a row matrix.
        /// </summary>
        /// <param name="row">Row to be fetched.</param>
        /// <returns>a row matrix.</returns>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        FieldMatrix<T> getRowMatrix(int row);

        /// <summary>
        /// Set the entries in row number <c>row</c> as a row matrix. 
        /// </summary>
        /// <param name="row">Row to be set.</param>
        /// <param name="matrix">Row matrix (must have one row and the same number
        /// of columns as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.</exception>
        /// <exception cref="MatrixDimensionMismatchException">
        /// if the matrix dimensions do not match one instance row.</exception>
        void setRowMatrix(int row, FieldMatrix<T> matrix);

        /// <summary>
        /// Get the entries in column number <c>column</c>
        /// as a column matrix.
        /// </summary>
        /// <param name="column">Column to be fetched.</param>
        /// <returns>a column matrix.</returns>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.</exception>
        FieldMatrix<T> getColumnMatrix(int column);

        /// <summary>
        /// Set the entries in column number <c>column</c>
        /// as a column matrix.
        /// </summary>
        /// <param name="column">Column to be set.</param>
        /// <param name="matrix">column matrix (must have one column and the same
        /// number of rows as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.</exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the matrix dimensions do
        /// not match one instance column.</exception>
        void setColumnMatrix(int column, FieldMatrix<T> matrix);

        /// <summary>
        /// Get the entries in row number <c>row</c>
        /// as a vector.
        /// </summary>
        /// <param name="row">Row to be fetched</param>
        /// <returns>a row vector.</returns>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.</exception>
        FieldVector<T> getRowVector(int row);

        /// <summary>
        /// Set the entries in row number <c>row</c>
        /// as a vector.
        /// </summary>
        /// <param name="row">Row to be set.</param>
        /// <param name="vector">row vector (must have the same number of columns
        /// as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.</exception>
        /// <exception cref="MatrixDimensionMismatchException">
        /// if the vector dimensions do not match one instance row.</exception>
        void setRowVector(int row, FieldVector<T> vector);

        /// <summary>
        /// Returns the entries in column number <c>column</c>
        /// as a vector.
        /// </summary>
        /// <param name="column">Column to be fetched.</param>
        /// <returns>a column vector.</returns>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        FieldVector<T> getColumnVector(int column);

        /// <summary>
        /// Set the entries in column number <c>column</c>
        /// as a vector.
        /// </summary>
        /// <param name="column">Column to be set.</param>
        /// <param name="vector">Column vector (must have the same number of rows
        /// as the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the vector dimension does not
        /// match one instance column.</exception>
        void setColumnVector(int column, FieldVector<T> vector);

        /// <summary>
        /// Get the entries in row number <c>row</c> as an array.
        /// </summary>
        /// <param name="row">Row to be fetched.</param>
        /// <returns>array of entries in the row.</returns>
        /// <exception cref="OutOfRangeException"> if the specified row index is not valid.
        /// </exception>
        T[] getRow(int row);

        /// <summary>
        /// Set the entries in row number <c>row</c>
        /// as a row matrix.
        /// </summary>
        /// <param name="row">Row to be set.</param>
        /// <param name="array">Row matrix (must have the same number of columns as
        /// the instance).</param>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the array size does not match
        /// one instance row.</exception>
        void setRow(int row, T[] array);

        /// <summary>
        /// Get the entries in column number <c>col</c> as an array.
        /// </summary>
        /// <param name="column">the column to be fetched</param>
        /// <returns>array of entries in the column</returns>
        /// <exception cref="OutOfRangeException"> if the specified column index is not valid.
        /// </exception>
        T[] getColumn(int column);

        /// <summary>
        /// Set the entries in column number <c>column</c>
        /// as a column matrix.
        /// </summary>
        /// <param name="column">the column to be set</param>
        /// <param name="array">column array (must have the same number of rows as the instance)
        /// </param>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid.
        /// </exception>
        /// <exception cref="MatrixDimensionMismatchException"> if the array size does not match
        /// one instance column.</exception>
        void setColumn(int column, T[] array);

        /// <summary>
        /// Returns the entry in the specified row and column.
        /// </summary>
        /// <param name="row">row location of entry to be fetched</param>
        /// <param name="column">column location of entry to be fetched</param>
        /// <returns>matrix entry in row,column</returns>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.
        /// </exception>
        T getEntry(int row, int column);

        /// <summary>
        /// Set the entry in the specified row and column.
        /// </summary>
        /// <param name="row">row location of entry to be set</param>
        /// <param name="column">column location of entry to be set</param>
        /// <param name="value">matrix entry to be set in row,column</param>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.
        /// </exception>
        void setEntry(int row, int column, T value);

        /// <summary>
        /// Change an entry in the specified row and column.
        /// </summary>
        /// <param name="row">Row location of entry to be set.</param>
        /// <param name="column">Column location of entry to be set.</param>
        /// <param name="increment">Value to add to the current matrix entry in
        /// <c>(row, column)</c>.</param>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.</exception>
        void addToEntry(int row, int column, T increment);

        /// <summary>
        /// Change an entry in the specified row and column.
        /// </summary>
        /// <param name="row">Row location of entry to be set.</param>
        /// <param name="column">Column location of entry to be set.</param>
        /// <param name="factor">Multiplication factor for the current matrix entry
        /// in <c>(row,column)</c></param>
        /// <exception cref="OutOfRangeException"> if the row or column index is not valid.</exception>
        void multiplyEntry(int row, int column, T factor);

        /// <summary>
        /// Returns the transpose of this matrix.
        /// </summary>
        /// <returns>transpose matrix</returns>
        FieldMatrix<T> transpose();

        /// <summary>
        /// Returns the <a href="http://mathworld.wolfram.com/MatrixTrace.html">
        /// trace</a> of the matrix (the sum of the elements on the main diagonal).
        /// </summary>
        /// <returns>trace</returns>
        /// <exception cref="NonSquareMatrixException"> if the matrix is not square.</exception>
        T getTrace();

        /// <summary>
        /// Returns the result of multiplying this by the vector <c>v</c>.
        /// </summary>
        /// <param name="v">the vector to operate on</param>
        /// <returns><c>this * v</c></returns>
        /// <exception cref="DimensionMismatchException"> if the number of columns of
        /// <c>this</c> matrix is not equal to the size of the vector <c>v</c>.</exception>
        T[] operate(T[] v);

        /// <summary>
        /// Returns the result of multiplying this by the vector <c>v</c>.
        /// </summary>
        /// <param name="v">the vector to operate on</param>
        /// <returns><c>this * v</c></returns>
        /// <exception cref="DimensionMismatchException"> if the number of columns of
        /// <c>this</c> matrix is not equal to the size of the vector <c>v</c>.</exception>
        FieldVector<T> operate(FieldVector<T> v);

        /// <summary>
        /// Returns the (row) vector result of premultiplying this by the vector
        /// <c>v</c>.
        /// </summary>
        /// <param name="v">the row vector to premultiply by</param>
        /// <returns><c>v * this</c></returns>
        /// <exception cref="DimensionMismatchException"> if the number of rows of <c>this</c>
        /// matrix is not equal to the size of the vector <c>v</c></exception>
        T[] preMultiply(T[] v);

        /// <summary>
        /// Returns the (row) vector result of premultiplying this by the vector
        /// <c>v</c>.
        /// </summary>
        /// <param name="v">the row vector to premultiply by</param>
        /// <returns><c>v * this</c></returns>
        /// <exception cref="DimensionMismatchException"> if the number of rows of <c>this</c>
        /// matrix is not equal to the size of the vector <c>v</c></exception>
        FieldVector<T> preMultiply(FieldVector<T> v);

        /// <summary>
        ///  Visit (and possibly change) all matrix entries in row order.
        ///  <para>Row order starts at upper left and iterating through all elements
        ///  of a row from left to right before going to the leftmost element
        ///  of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor visitor used to process all matrix entries
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/></param>
        /// <returns>the value returned by <see cref="FieldMatrixChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor);

        /// <summary>
        /// Visit (but don't change) all matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/></param>
        /// <returns>the value returned by <see cref="FieldMatrixPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor);

        /// <summary>
        /// Visit (and possibly change) some matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index</param>
        /// <returns>the value returned by <see cref="FieldMatrixChangingVisitor.end()"/>
        /// at the end  of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (but don't change) some matrix entries in row order.
        /// <para>Row order starts at upper left and iterating through all elements
        /// of a row from left to right before going to the leftmost element
        /// of the next row.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index</param>
        /// <returns>the value returned by <see cref="FieldMatrixPreservingVisitor.end()"/> 
        /// at the end of the walk</returns>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (and possibly change) all matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <returns>the value returned by <see cref="FieldMatrixChangingVisitor.end()"/> at 
        /// the end of the walk</returns>
        T walkInColumnOrder(FieldMatrixChangingVisitor<T> visitor);

        /// <summary>
        /// Visit (but don't change) all matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <returns>the value returned by <see cref="FieldMatrixPreservingVisitor.end()"/> at 
        /// the end of the walk</returns>
        T walkInColumnOrder(FieldMatrixPreservingVisitor<T> visitor);

        /// <summary>
        /// Visit (and possibly change) some matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index</param>
        /// <returns>the value returned by <see cref="FieldMatrixChangingVisitor.end()"/> at 
        /// the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        T walkInColumnOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (but don't change) some matrix entries in column order.
        /// <para>Column order starts at upper left and iterating through all elements
        /// of a column from top to bottom before going to the topmost element
        /// of the next column.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index</param>
        /// <returns>the value returned by <see cref="FieldMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        T walkInColumnOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (and possibly change) all matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <returns>the value returned by <see cref="FieldMatrixChangingVisitor.end()"/> at 
        /// the end of the walk</returns>
        T walkInOptimizedOrder(FieldMatrixChangingVisitor<T> visitor);

        /// <summary>
        /// Visit (but don't change) all matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <returns>the value returned by <see cref="FieldMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        T walkInOptimizedOrder(FieldMatrixPreservingVisitor<T> visitor);

        /// <summary>
        /// Visit (and possibly change) some matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index (inclusive)</param>
        /// <returns>the value returned by <see cref="FieldMatrixChangingVisitor.end()"/> at 
        /// the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        T walkInOptimizedOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit (but don't change) some matrix entries using the fastest possible order.
        /// <para>The fastest walking order depends on the exact matrix class. It may be
        /// different from traditional row or column orders.</para>
        /// </summary>
        /// <param name="visitor">visitor used to process all matrix entries<para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInRowOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInRowOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInColumnOrder(FieldMatrixPreservingVisitor, int, int, int, int)"/>
        /// <para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixPreservingVisitor)"/><para/>
        /// See <see cref="walkInOptimizedOrder(FieldMatrixChangingVisitor, int, int, int, int)"/>
        /// </param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">column index (inclusive)</param>
        /// <returns>the value returned by <see cref="FieldMatrixPreservingVisitor.end()"/> at
        /// the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>endRow < startRow</c> or
        /// <c>endColumn < startColumn</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        T walkInOptimizedOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn);
    }
}
