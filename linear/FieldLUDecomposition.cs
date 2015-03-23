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
    /// Calculates the LUP-decomposition of a square matrix.
    ///  <para>The LUP-decomposition of a matrix A consists of three matrices
    ///  L, U and P that satisfy: PA = LU, L is lower triangular, and U is
    ///  upper triangular and P is a permutation matrix. All matrices are
    ///  m&times;m.</para>
    ///  <para>Since <see cref="FieldElement">field elements</see> do not provide an ordering
    ///  operator, the permutation matrix is computed here only in order to avoid
    ///  a zero pivot element, no attempt is done to get the largest pivot
    ///  element.</para>
    ///  <para>This class is based on the class with similar name from the
    ///  <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</para>
    ///  <list type="bullet">
    ///    <item>a <see cref="getP()"/> method has been added,</item>
    ///    <item>the <c>det</c> method has been renamed as <see cref="getDeterminant()"/>
    ///    ,</item>
    ///    <item>the <c>getDoublePivot</c> method has been removed (but the int based
    ///    <see cref="getPivot()"/> method has been kept),</item>
    ///    <item>the <c>solve</c> and <c>isNonSingular</c> methods have been replaced
    ///    by a <see cref="getSolver()"/> method and the equivalent methods
    ///    provided by the returned <see cref="DecompositionSolver"/>.</item>
    ///  </list>
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    /// <remarks>
    /// See <a href="http://mathworld.wolfram.com/LUDecomposition.html">MathWorld</a><para/>
    /// See <a href="http://en.wikipedia.org/wiki/LU_decomposition">Wikipedia</a>
    /// </remarks>
    public class FieldLUDecomposition<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Field to which the elements belong.
        /// </summary>
        private readonly Field<T> field;

        /// <summary>
        /// Entries of LU decomposition.
        /// </summary>
        private T[][] lu;

        /// <summary>
        /// Pivot permutation associated with LU decomposition.
        /// </summary>
        private int[] pivot;

        /// <summary>
        /// Parity of the permutation associated with the LU decomposition.
        /// </summary>
        private Boolean even;

        /// <summary>
        /// Singularity indicator.
        /// </summary>
        private Boolean singular;

        /// <summary>
        /// Cached value of L.
        /// </summary>
        private FieldMatrix<T> cachedL;

        /// <summary>
        /// Cached value of U.
        /// </summary>
        private FieldMatrix<T> cachedU;

        /// <summary>
        /// Cached value of P.
        /// </summary>
        private FieldMatrix<T> cachedP;

        /// <summary>
        /// Calculates the LU-decomposition of the given matrix.
        /// </summary>
        /// <param name="matrix">The matrix to decompose.</param>
        /// <exception cref="NonSquareMatrixException"> if matrix is not square</exception>
        public FieldLUDecomposition(FieldMatrix<T> matrix)
        {
            if (!matrix.isSquare())
            {
                throw new NonSquareMatrixException(matrix.getRowDimension(), matrix.getColumnDimension());
            }

            int m = matrix.getColumnDimension();
            field = matrix.getField();
            lu = matrix.getData();
            pivot = new int[m];
            cachedL = null;
            cachedU = null;
            cachedP = null;

            // Initialize permutation array and parity
            for (int row = 0; row < m; row++)
            {
                pivot[row] = row;
            }
            even = true;
            singular = false;

            // Loop over columns
            for (int col = 0; col < m; col++)
            {

                T sum = field.getZero();

                // upper
                for (int row = 0; row < col; row++)
                {
                    T[] luRow = lu[row];
                    sum = luRow[col];
                    for (int i = 0; i < row; i++)
                    {
                        sum = sum.subtract(luRow[i].multiply(lu[i][col]));
                    }
                    luRow[col] = sum;
                }

                // lower
                int nonZero = col; // permutation row
                for (int row = col; row < m; row++)
                {
                    T[] luRow = lu[row];
                    sum = luRow[col];
                    for (int i = 0; i < col; i++)
                    {
                        sum = sum.subtract(luRow[i].multiply(lu[i][col]));
                    }
                    luRow[col] = sum;

                    if (lu[nonZero][col].Equals(field.getZero()))
                    {
                        // try to select a better permutation choice
                        ++nonZero;
                    }
                }

                // Singularity check
                if (nonZero >= m)
                {
                    singular = true;
                    return;
                }

                // Pivot if necessary
                if (nonZero != col)
                {
                    T tmp = field.getZero();
                    for (int i = 0; i < m; i++)
                    {
                        tmp = lu[nonZero][i];
                        lu[nonZero][i] = lu[col][i];
                        lu[col][i] = tmp;
                    }
                    int temp = pivot[nonZero];
                    pivot[nonZero] = pivot[col];
                    pivot[col] = temp;
                    even = !even;
                }

                // Divide the lower elements by the "winning" diagonal elt.
                T luDiag = lu[col][col];
                for (int row = col + 1; row < m; row++)
                {
                    T[] luRow = lu[row];
                    luRow[col] = luRow[col].divide(luDiag);
                }
            }

        }

        /// <summary>
        /// Returns the matrix L of the decomposition.
        /// <para>L is a lower-triangular matrix</para>
        /// </summary>
        /// <returns>the L matrix (or null if decomposed matrix is singular)</returns>
        public FieldMatrix<T> getL()
        {
            if ((cachedL == null) && !singular)
            {
                int m = pivot.Length;
                cachedL = new Array2DRowFieldMatrix<T>(field, m, m);
                for (int i = 0; i < m; ++i)
                {
                    T[] luI = lu[i];
                    for (int j = 0; j < i; ++j)
                    {
                        cachedL.setEntry(i, j, luI[j]);
                    }
                    cachedL.setEntry(i, i, field.getOne());
                }
            }
            return cachedL;
        }

        /// <summary>
        /// Returns the matrix U of the decomposition.
        /// <para>U is an upper-triangular matrix</para>
        /// </summary>
        /// <returns>the U matrix (or null if decomposed matrix is singular)</returns>
        public FieldMatrix<T> getU()
        {
            if ((cachedU == null) && !singular)
            {
                int m = pivot.Length;
                cachedU = new Array2DRowFieldMatrix<T>(field, m, m);
                for (int i = 0; i < m; ++i)
                {
                    T[] luI = lu[i];
                    for (int j = i; j < m; ++j)
                    {
                        cachedU.setEntry(i, j, luI[j]);
                    }
                }
            }
            return cachedU;
        }

        /// <summary>
        /// Returns the P rows permutation matrix.
        /// <para>P is a sparse matrix with exactly one element set to 1.0 in
        /// each row and each column, all other elements being set to 0.0.</para>
        /// <para>The positions of the 1 elements are given by the <see cref="getPivot()">
        /// pivot permutation vector</see>.</para>
        /// </summary>
        /// <returns>the P rows permutation matrix (or null if decomposed matrix is singular)
        /// </returns>
        /// <remarks>
        /// See <see cref="getPivot()"/>
        /// </remarks>
        public FieldMatrix<T> getP()
        {
            if ((cachedP == null) && !singular)
            {
                int m = pivot.Length;
                cachedP = new Array2DRowFieldMatrix<T>(field, m, m);
                for (int i = 0; i < m; ++i)
                {
                    cachedP.setEntry(i, pivot[i], field.getOne());
                }
            }
            return cachedP;
        }

        /// <summary>
        /// Returns the pivot permutation vector.
        /// </summary>
        /// <returns>the pivot permutation vector</returns>
        /// <remarks>
        /// See <see cref="getP()"/>
        /// </remarks>
        public int[] getPivot()
        {
            return (Int32[])pivot.Clone();
        }

        /// <summary>
        /// Return the determinant of the matrix.
        /// </summary>
        /// <returns>determinant of the matrix</returns>
        public T getDeterminant()
        {
            if (singular)
            {
                return field.getZero();
            }
            else
            {
                int m = pivot.Length;
                T determinant = even ? field.getOne() : field.getZero().subtract(field.getOne());
                for (int i = 0; i < m; i++)
                {
                    determinant = determinant.multiply(lu[i][i]);
                }
                return determinant;
            }
        }

        /// <summary>
        /// Get a solver for finding the A &times; X = B solution in exact linear sense.
        /// </summary>
        /// <returns>a solver</returns>
        public FieldDecompositionSolver<T> getSolver()
        {
            return new Solver<T>(field, lu, pivot, singular);
        }

        /// <summary>
        /// Specialized solver.
        /// </summary>
        /// <typeparam name="U"></typeparam>
        private class Solver<U> : FieldDecompositionSolver<U> where U : FieldElement<U>
        {

            /// <summary>
            /// Field to which the elements belong.
            /// </summary>
            private readonly Field<U> field;

            /// <summary>
            /// Entries of LU decomposition.
            /// </summary>
            private readonly U[][] lu;

            /// <summary>
            /// Pivot permutation associated with LU decomposition.
            /// </summary>
            private readonly int[] pivot;

            /// <summary>
            /// Singularity indicator.
            /// </summary>
            private readonly Boolean singular;

            /// <summary>
            /// Build a solver from decomposed matrix.
            /// </summary>
            /// <param name="field">field to which the matrix elements belong</param>
            /// <param name="lu">entries of LU decomposition</param>
            /// <param name="pivot">pivot permutation associated with LU decomposition</param>
            /// <param name="singular">singularity indicator</param>
            internal Solver(Field<U> field, U[][] lu, int[] pivot, Boolean singular)
            {
                this.field = field;
                this.lu = lu;
                this.pivot = pivot;
                this.singular = singular;
            }

            /// <inheritdoc/>
            public Boolean isNonSingular()
            {
                return !singular;
            }

            /// <inheritdoc/>
            public FieldVector<U> solve(FieldVector<U> b)
            {
                try
                {
                    return solve((ArrayFieldVector<U>)b);
                }
                catch (InvalidCastException)
                {

                    int m = pivot.Length;
                    if (b.getDimension() != m)
                    {
                        throw new DimensionMismatchException(b.getDimension(), m);
                    }
                    if (singular)
                    {
                        throw new SingularMatrixException();
                    }

                    // Apply permutations to b
                    U[] bp = MathArrays.buildArray(field, m);
                    for (int row = 0; row < m; row++)
                    {
                        bp[row] = b.getEntry(pivot[row]);
                    }

                    // Solve LY = b
                    for (int col = 0; col < m; col++)
                    {
                        U bpCol = bp[col];
                        for (int i = col + 1; i < m; i++)
                        {
                            bp[i] = bp[i].subtract(bpCol.multiply(lu[i][col]));
                        }
                    }

                    // Solve UX = Y
                    for (int col = m - 1; col >= 0; col--)
                    {
                        bp[col] = bp[col].divide(lu[col][col]);
                        U bpCol = bp[col];
                        for (int i = 0; i < col; i++)
                        {
                            bp[i] = bp[i].subtract(bpCol.multiply(lu[i][col]));
                        }
                    }

                    return new ArrayFieldVector<U>(field, bp, false);

                }
            }

            /// <summary>
            /// Solve the linear equation A &times; X = B.
            /// <para>The A matrix is implicit here. It is </para>
            /// </summary>
            /// <param name="b">right-hand side of the equation A &times; X = B</param>
            /// <returns>a vector X such that A &times; X = B</returns>
            /// <exception cref="DimensionMismatchException"> if the matrices dimensions do
            /// not match.</exception>
            /// <exception cref="SingularMatrixException"> if the decomposed matrix is singular.
            /// </exception>
            public ArrayFieldVector<U> solve(ArrayFieldVector<U> b)
            {
                int m = pivot.Length;
                int length = b.getDimension();
                if (length != m)
                {
                    throw new DimensionMismatchException(length, m);
                }
                if (singular)
                {
                    throw new SingularMatrixException();
                }

                // Apply permutations to b
                U[] bp = MathArrays.buildArray(field, m);
                for (int row = 0; row < m; row++)
                {
                    bp[row] = b.getEntry(pivot[row]);
                }

                // Solve LY = b
                for (int col = 0; col < m; col++)
                {
                    U bpCol = bp[col];
                    for (int i = col + 1; i < m; i++)
                    {
                        bp[i] = bp[i].subtract(bpCol.multiply(lu[i][col]));
                    }
                }

                // Solve UX = Y
                for (int col = m - 1; col >= 0; col--)
                {
                    bp[col] = bp[col].divide(lu[col][col]);
                    U bpCol = bp[col];
                    for (int i = 0; i < col; i++)
                    {
                        bp[i] = bp[i].subtract(bpCol.multiply(lu[i][col]));
                    }
                }

                return new ArrayFieldVector<U>(bp, false);
            }

            /// <inheritdoc/>
            public FieldMatrix<U> solve(FieldMatrix<U> b)
            {
                int m = pivot.Length;
                if (b.getRowDimension() != m)
                {
                    throw new DimensionMismatchException(b.getRowDimension(), m);
                }
                if (singular)
                {
                    throw new SingularMatrixException();
                }

                int nColB = b.getColumnDimension();

                // Apply permutations to b
                U[][] bp = MathArrays.buildArray(field, m, nColB);
                for (int row = 0; row < m; row++)
                {
                    U[] bpRow = bp[row];
                    int pRow = pivot[row];
                    for (int col = 0; col < nColB; col++)
                    {
                        bpRow[col] = b.getEntry(pRow, col);
                    }
                }

                // Solve LY = b
                for (int col = 0; col < m; col++)
                {
                    U[] bpCol = bp[col];
                    for (int i = col + 1; i < m; i++)
                    {
                        U[] bpI = bp[i];
                        U luICol = lu[i][col];
                        for (int j = 0; j < nColB; j++)
                        {
                            bpI[j] = bpI[j].subtract(bpCol[j].multiply(luICol));
                        }
                    }
                }

                // Solve UX = Y
                for (int col = m - 1; col >= 0; col--)
                {
                    U[] bpCol = bp[col];
                    U luDiag = lu[col][col];
                    for (int j = 0; j < nColB; j++)
                    {
                        bpCol[j] = bpCol[j].divide(luDiag);
                    }
                    for (int i = 0; i < col; i++)
                    {
                        U[] bpI = bp[i];
                        U luICol = lu[i][col];
                        for (int j = 0; j < nColB; j++)
                        {
                            bpI[j] = bpI[j].subtract(bpCol[j].multiply(luICol));
                        }
                    }
                }

                return new Array2DRowFieldMatrix<U>(field, bp, false);

            }

            /// <inheritdoc/>
            public FieldMatrix<U> getInverse()
            {
                int m = pivot.Length;
                U one = field.getOne();
                FieldMatrix<U> identity = new Array2DRowFieldMatrix<U>(field, m, m);
                for (int i = 0; i < m; ++i)
                {
                    identity.setEntry(i, i, one);
                }
                return solve(identity);
            }
        }
    }
}