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
    /// Cache-friendly implementation of FieldMatrix using a flat arrays to store
    /// square blocks of the matrix.
    /// <para>
    /// This implementation is specially designed to be cache-friendly. Square blocks are
    /// stored as small arrays and allow efficient traversal of data both in row major direction
    /// and columns major direction, one block at a time. This greatly increases performances
    /// for algorithms that use crossed directions loops like multiplication or transposition.
    /// </para>
    /// <para>
    /// The size of square blocks is a static parameter. It may be tuned according to the cache
    /// size of the target computer processor. As a rule of thumbs, it should be the largest
    /// value that allows three blocks to be simultaneously cached (this is necessary for example
    /// for matrix multiplication). The default value is to use 36x36 blocks.
    /// </para>
    /// <para>
    /// The regular blocks represent <see cref="BLOCK_SIZE"/> x <see cref="BLOCK_SIZE"/>
    /// squares. Blocks at right hand side and bottom side which may be smaller to fit matrix 
    /// dimensions. The square blocks are flattened in row major order in single dimension arrays
    /// which are therefore <see cref="BLOCK_SIZE"/><sup>2</sup> elements long for regular blocks.
    /// The blocks are themselves organized in row major order.
    /// </para>
    /// <para>
    /// As an example, for a block size of 36x36, a 100x60 matrix would be stored in 6 blocks.
    /// Block 0 would be a Field[1296] array holding the upper left 36x36 square, block 1 would be
    /// a Field[1296] array holding the upper center 36x36 square, block 2 would be a Field[1008]
    /// array holding the upper right 36x28 rectangle, block 3 would be a Field[864] array holding
    /// the lower left 24x36 rectangle, block 4 would be a Field[864] array holding the lower
    /// center 24x36 rectangle and block 5 would be a Field[672] array holding the lower right 
    /// 24x28 rectangle.
    /// </para>
    /// <para>
    /// The layout complexity overhead versus simple mapping of matrices to java
    /// arrays is negligible for small matrices (about 1%). The gain from cache efficiency leads
    /// to up to 3-fold improvements for matrices of moderate to large size.
    /// </para>
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public class BlockFieldMatrix<T> : AbstractFieldMatrix<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Block size.
        /// </summary>
        public const int BLOCK_SIZE = 36;

        /// <summary>
        /// Blocks of matrix entries.
        /// </summary>
        private readonly T[][] blocks;

        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        private readonly int rows;

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        private readonly int columns;

        /// <summary>
        /// Number of block rows of the matrix.
        /// </summary>
        private readonly int blockRows;

        /// <summary>
        /// Number of block columns of the matrix.
        /// </summary>
        private readonly int blockColumns;

        /// <summary>
        /// Create a new matrix with the supplied row and column dimensions.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="rows">Number of rows in the new matrix.</param>
        /// <param name="columns">Number of columns in the new matrix.</param>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension is not
        /// positive.</exception>
        public BlockFieldMatrix(Field<T> field, int rows, int columns)
            : base(field, rows, columns)
        {
            this.rows = rows;
            this.columns = columns;

            // number of blocks
            blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
            blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

            // allocate storage blocks, taking care of smaller ones at right and bottom
            blocks = createBlocksLayout(field, rows, columns);
        }

        /// <summary>
        /// Create a new dense matrix copying entries from raw layout data.
        /// <para>The input array <em>must</em> already be in raw layout.</para>
        /// <para>Calling this constructor is equivalent to call:
        /// <code>matrix = new BlockFieldMatrix(T)(getField(), rawData.Length, rawData[0].Length,
        ///                                   toBlocksLayout(rawData), false);</code>
        /// </para>
        /// </summary>
        /// <param name="rawData">Data for the new matrix, in raw layout.</param>
        /// <exception cref="DimensionMismatchException"> if the <c>blockData</c> shape is
        /// inconsistent with block layout.</exception>
        /// <remarks>
        /// See <see cref="BlockFieldMatrix(int, int, FieldElement[][], boolean)"/>
        /// </remarks>
        public BlockFieldMatrix(T[][] rawData) : this(rawData.Length, rawData[0].Length, toBlocksLayout(rawData), false) { }

        /// <summary>
        /// Create a new dense matrix copying entries from block layout data.
        /// <para>The input array must already be in blocks layout.</para>
        /// </summary>
        /// <param name="rows">the number of rows in the new matrix</param>
        /// <param name="columns">the number of columns in the new matrix</param>
        /// <param name="blockData">data for new matrix</param>
        /// <param name="copyArray">if true, the input array will be copied, otherwise
        /// it will be referenced</param>
        /// <exception cref="DimensionMismatchException"> if the <c>blockData</c> shape is
        /// inconsistent with block layout.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if row or column dimension is not
        /// positive.</exception>
        /// <remarks>
        /// See <see cref="createBlocksLayout(Field, int, int)"/><para/>
        /// See <see cref="toBlocksLayout(FieldElement[][])"/><para/>
        /// See <see cref="BlockFieldMatrix(FieldElement[][])"/>
        /// </remarks>
        public BlockFieldMatrix(int rows, int columns, T[][] blockData, Boolean copyArray)
            : base(extractField(blockData), rows, columns)
        {
            this.rows = rows;
            this.columns = columns;

            // number of blocks
            blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
            blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

            if (copyArray)
            {
                // allocate storage blocks, taking care of smaller ones at right and bottom
                blocks = MathArrays.buildArray(getField(), blockRows * blockColumns, -1);
            }
            else
            {
                // reference existing array
                blocks = blockData;
            }

            int index = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                for (int jBlock = 0; jBlock < blockColumns; ++jBlock, ++index)
                {
                    if (blockData[index].Length != iHeight * blockWidth(jBlock))
                    {
                        throw new DimensionMismatchException(blockData[index].Length,
                                                             iHeight * blockWidth(jBlock));
                    }
                    if (copyArray)
                    {
                        blocks[index] = (T[])blockData[index].Clone();
                    }
                }
            }
        }

        /// <summary>
        /// Convert a data array from raw layout to blocks layout.
        /// <para>
        /// Raw layout is the straightforward layout where element at row i and
        /// column j is in array element <code>rawData[i][j]</code>. Blocks layout
        /// is the layout used in <see cref="BlockFieldMatrix"/> instances, where the matrix
        /// is split in square blocks (except at right and bottom side where blocks may
        /// be rectangular to fit matrix size) and each block is stored in a flattened
        /// one-dimensional array.
        /// </para>
        /// <para>
        /// This method creates an array in blocks layout from an input array in raw layout.
        /// It can be used to provide the array argument of the
        /// <see cref="BlockFieldMatrix(int, int, FieldElement[][], boolean)"/>
        /// constructor.
        /// </para>
        /// </summary>
        /// <typeparam name="U">Type of the field elements.</typeparam>
        /// <param name="rawData">Data array in raw layout.</param>
        /// <returns>a new data array containing the same entries but in blocks layout</returns>
        /// <exception cref="DimensionMismatchException"> if <c>rawData</c> is not rectangular
        /// (not all rows have the same length).</exception>
        public static U[][] toBlocksLayout<U>(U[][] rawData) where U : FieldElement<U>
        {
            int rows = rawData.Length;
            int columns = rawData[0].Length;
            int blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
            int blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

            // safety checks
            for (int i = 0; i < rawData.Length; ++i)
            {
                int length = rawData[i].Length;
                if (length != columns)
                {
                    throw new DimensionMismatchException(columns, length);
                }
            }

            // convert array
            Field<U> field = extractField(rawData);
            U[][] blocks = MathArrays.buildArray(field, blockRows * blockColumns, -1);
            int blockIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                int iHeight = pEnd - pStart;
                for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                {
                    int qStart = jBlock * BLOCK_SIZE;
                    int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                    int jWidth = qEnd - qStart;

                    // allocate new block
                    U[] block = MathArrays.buildArray(field, iHeight * jWidth);
                    blocks[blockIndex] = block;

                    // copy data
                    int index = 0;
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        Array.Copy(rawData[p], qStart, block, index, jWidth);
                        index += jWidth;
                    }

                    ++blockIndex;
                }
            }

            return blocks;
        }

        /// <summary>
        /// Create a data array in blocks layout.
        /// <para>
        /// This method can be used to create the array argument of the
        /// <see cref="BlockFieldMatrix(int, int, FieldElement[][], boolean)"/>
        /// constructor.
        /// </para>
        /// </summary>
        /// <typeparam name="U">Type of the field elements.</typeparam>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="rows">Number of rows in the new matrix.</param>
        /// <param name="columns">Number of columns in the new matrix.</param>
        /// <returns>a new data array in blocks layout.</returns>
        /// <remarks>
        /// See <see cref="toBlocksLayout(FieldElement[][])"/>
        /// </remarks>
        public static U[][] createBlocksLayout<U>(Field<U> field, int rows, int columns) where U : FieldElement<U>
        {
            int blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
            int blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

            U[][] blocks = MathArrays.buildArray(field, blockRows * blockColumns, -1);
            int blockIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                int iHeight = pEnd - pStart;
                for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                {
                    int qStart = jBlock * BLOCK_SIZE;
                    int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                    int jWidth = qEnd - qStart;
                    blocks[blockIndex] = MathArrays.buildArray(field, iHeight * jWidth);
                    ++blockIndex;
                }
            }

            return blocks;
        }

        /// <inheritdoc/>
        public override FieldMatrix<T> createMatrix(int rowDimension, int columnDimension)
        {
            return new BlockFieldMatrix<T>(getField(), rowDimension, columnDimension);
        }

        /// <inheritdoc/>
        public override FieldMatrix<T> copy()
        {

            // create an empty matrix
            BlockFieldMatrix<T> copied = new BlockFieldMatrix<T>(getField(), rows, columns);

            // copy the blocks
            for (int i = 0; i < blocks.Length; ++i)
            {
                Array.Copy(blocks[i], 0, copied.blocks[i], 0, blocks[i].Length);
            }

            return copied;
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> add(FieldMatrix<T> m)
        {
            try
            {
                return add((BlockFieldMatrix<T>)m);
            }
            catch (InvalidCastException)
            {

                // safety check
                checkAdditionCompatible(m);

                BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, columns);

                // perform addition block-wise, to ensure good cache behavior
                int blockIndex = 0;
                for (int iBlock = 0; iBlock < outp.blockRows; ++iBlock)
                {
                    for (int jBlock = 0; jBlock < outp.blockColumns; ++jBlock)
                    {

                        // perform addition on the current block
                        T[] outBlock = outp.blocks[blockIndex];
                        T[] tBlock = blocks[blockIndex];
                        int pStart = iBlock * BLOCK_SIZE;
                        int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                        int qStart = jBlock * BLOCK_SIZE;
                        int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                        int k = 0;
                        for (int p = pStart; p < pEnd; ++p)
                        {
                            for (int q = qStart; q < qEnd; ++q)
                            {
                                outBlock[k] = tBlock[k].add(m.getEntry(p, q));
                                ++k;
                            }
                        }

                        // go to next block
                        ++blockIndex;
                    }
                }

                return outp;
            }
        }

        /// <summary>
        /// Compute the sum of <c>this</c> and <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to be added</param>
        /// <returns><c>this + m</c></returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c></exception>
        public BlockFieldMatrix<T> add(BlockFieldMatrix<T> m)
        {

            // safety check
            checkAdditionCompatible(m);

            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, columns);

            // perform addition block-wise, to ensure good cache behavior
            for (int blockIndex = 0; blockIndex < outp.blocks.Length; ++blockIndex)
            {
                T[] outBlock = outp.blocks[blockIndex];
                T[] tBlock = blocks[blockIndex];
                T[] mBlock = m.blocks[blockIndex];
                for (int k = 0; k < outBlock.Length; ++k)
                {
                    outBlock[k] = tBlock[k].add(mBlock[k]);
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> subtract(FieldMatrix<T> m)
        {
            try
            {
                return subtract((BlockFieldMatrix<T>)m);
            }
            catch (InvalidCastException)
            {

                // safety check
                checkSubtractionCompatible(m);

                BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, columns);

                // perform subtraction block-wise, to ensure good cache behavior
                int blockIndex = 0;
                for (int iBlock = 0; iBlock < outp.blockRows; ++iBlock)
                {
                    for (int jBlock = 0; jBlock < outp.blockColumns; ++jBlock)
                    {

                        // perform subtraction on the current block
                        T[] outBlock = outp.blocks[blockIndex];
                        T[] tBlock = blocks[blockIndex];
                        int pStart = iBlock * BLOCK_SIZE;
                        int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                        int qStart = jBlock * BLOCK_SIZE;
                        int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                        int k = 0;
                        for (int p = pStart; p < pEnd; ++p)
                        {
                            for (int q = qStart; q < qEnd; ++q)
                            {
                                outBlock[k] = tBlock[k].subtract(m.getEntry(p, q));
                                ++k;
                            }
                        }

                        // go to next block
                        ++blockIndex;
                    }
                }
                return outp;
            }
        }

        /// <summary>
        /// Compute <c>this - m</c>.
        /// </summary>
        /// <param name="m">matrix to be subtracted</param>
        /// <returns><c>this - m</c></returns>
        /// <exception cref="MatrixDimensionMismatchException"> if <c>m</c> is not the same
        /// size as <c>this</c></exception>
        public BlockFieldMatrix<T> subtract(BlockFieldMatrix<T> m)
        {
            // safety check
            checkSubtractionCompatible(m);

            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, columns);

            // perform subtraction block-wise, to ensure good cache behavior
            for (int blockIndex = 0; blockIndex < outp.blocks.Length; ++blockIndex)
            {
                T[] outBlock = outp.blocks[blockIndex];
                T[] tBlock = blocks[blockIndex];
                T[] mBlock = m.blocks[blockIndex];
                for (int k = 0; k < outBlock.Length; ++k)
                {
                    outBlock[k] = tBlock[k].subtract(mBlock[k]);
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> scalarAdd(T d)
        {
            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, columns);

            // perform subtraction block-wise, to ensure good cache behavior
            for (int blockIndex = 0; blockIndex < outp.blocks.Length; ++blockIndex)
            {
                T[] outBlock = outp.blocks[blockIndex];
                T[] tBlock = blocks[blockIndex];
                for (int k = 0; k < outBlock.Length; ++k)
                {
                    outBlock[k] = tBlock[k].add(d);
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> scalarMultiply(T d)
        {

            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, columns);

            // perform subtraction block-wise, to ensure good cache behavior
            for (int blockIndex = 0; blockIndex < outp.blocks.Length; ++blockIndex)
            {
                T[] outBlock = outp.blocks[blockIndex];
                T[] tBlock = blocks[blockIndex];
                for (int k = 0; k < outBlock.Length; ++k)
                {
                    outBlock[k] = tBlock[k].multiply(d);
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> multiply(FieldMatrix<T> m)
        {
            try
            {
                return multiply((BlockFieldMatrix<T>)m);
            }
            catch (InvalidCastException)
            {

                // safety check
                checkMultiplicationCompatible(m);

                BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, m.getColumnDimension());
                T zero = getField().getZero();

                // perform multiplication block-wise, to ensure good cache behavior
                int blockIndex = 0;
                for (int iBlock = 0; iBlock < outp.blockRows; ++iBlock)
                {

                    int pStart = iBlock * BLOCK_SIZE;
                    int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);

                    for (int jBlock = 0; jBlock < outp.blockColumns; ++jBlock)
                    {

                        int qStart = jBlock * BLOCK_SIZE;
                        int qEnd = FastMath.min(qStart + BLOCK_SIZE, m.getColumnDimension());

                        // select current block
                        T[] outBlock = outp.blocks[blockIndex];

                        // perform multiplication on current block
                        for (int kBlock = 0; kBlock < blockColumns; ++kBlock)
                        {
                            int kWidth = blockWidth(kBlock);
                            T[] tBlock = blocks[iBlock * blockColumns + kBlock];
                            int rStart = kBlock * BLOCK_SIZE;
                            int k = 0;
                            for (int p = pStart; p < pEnd; ++p)
                            {
                                int lStart = (p - pStart) * kWidth;
                                int lEnd = lStart + kWidth;
                                for (int q = qStart; q < qEnd; ++q)
                                {
                                    T sum = zero;
                                    int r = rStart;
                                    for (int l = lStart; l < lEnd; ++l)
                                    {
                                        sum = sum.add(tBlock[l].multiply(m.getEntry(r, q)));
                                        ++r;
                                    }
                                    outBlock[k] = outBlock[k].add(sum);
                                    ++k;
                                }
                            }
                        }

                        // go to next block
                        ++blockIndex;

                    }
                }
                return outp;
            }
        }

        /// <summary>
        /// Returns the result of postmultiplying <c>this</c> by <c>m</c>.
        /// </summary>
        /// <param name="m">matrix to postmultiply by</param>
        /// <returns><c>this * m</c></returns>
        /// <exception cref="DimensionMismatchException"> if the matrices are not compatible.
        /// </exception>
        public BlockFieldMatrix<T> multiply(BlockFieldMatrix<T> m)
        {
            // safety check
            checkMultiplicationCompatible(m);

            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, m.columns);
            T zero = getField().getZero();

            // perform multiplication block-wise, to ensure good cache behavior
            int blockIndex = 0;
            for (int iBlock = 0; iBlock < outp.blockRows; ++iBlock)
            {

                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);

                for (int jBlock = 0; jBlock < outp.blockColumns; ++jBlock)
                {
                    int jWidth = outp.blockWidth(jBlock);
                    int jWidth2 = jWidth + jWidth;
                    int jWidth3 = jWidth2 + jWidth;
                    int jWidth4 = jWidth3 + jWidth;

                    // select current block
                    T[] outBlock = outp.blocks[blockIndex];

                    // perform multiplication on current block
                    for (int kBlock = 0; kBlock < blockColumns; ++kBlock)
                    {
                        int kWidth = blockWidth(kBlock);
                        T[] tBlock = blocks[iBlock * blockColumns + kBlock];
                        T[] mBlock = m.blocks[kBlock * m.blockColumns + jBlock];
                        int k = 0;
                        for (int p = pStart; p < pEnd; ++p)
                        {
                            int lStart = (p - pStart) * kWidth;
                            int lEnd = lStart + kWidth;
                            for (int nStart = 0; nStart < jWidth; ++nStart)
                            {
                                T sum = zero;
                                int l = lStart;
                                int n = nStart;
                                while (l < lEnd - 3)
                                {
                                    sum = sum.
                                          add(tBlock[l].multiply(mBlock[n])).
                                          add(tBlock[l + 1].multiply(mBlock[n + jWidth])).
                                          add(tBlock[l + 2].multiply(mBlock[n + jWidth2])).
                                          add(tBlock[l + 3].multiply(mBlock[n + jWidth3]));
                                    l += 4;
                                    n += jWidth4;
                                }
                                while (l < lEnd)
                                {
                                    sum = sum.add(tBlock[l++].multiply(mBlock[n]));
                                    n += jWidth;
                                }
                                outBlock[k] = outBlock[k].add(sum);
                                ++k;
                            }
                        }
                    }

                    // go to next block
                    ++blockIndex;
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public new T[][] getData()
        {

            T[][] data = MathArrays.buildArray(getField(), getRowDimension(), getColumnDimension());
            int lastColumns = columns - (blockColumns - 1) * BLOCK_SIZE;

            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                int regularPos = 0;
                int lastPos = 0;
                for (int p = pStart; p < pEnd; ++p)
                {
                    T[] dataP = data[p];
                    int blockIndex = iBlock * blockColumns;
                    int dataPos = 0;
                    for (int jBlock = 0; jBlock < blockColumns - 1; ++jBlock)
                    {
                        Array.Copy(blocks[blockIndex++], regularPos, dataP, dataPos, BLOCK_SIZE);
                        dataPos += BLOCK_SIZE;
                    }
                    Array.Copy(blocks[blockIndex], lastPos, dataP, dataPos, lastColumns);
                    regularPos += BLOCK_SIZE;
                    lastPos += lastColumns;
                }
            }

            return data;
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> getSubMatrix(int startRow, int endRow, int startColumn, int endColumn)
        {
            // safety checks
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);

            // create the output matrix
            BlockFieldMatrix<T> outp =
                new BlockFieldMatrix<T>(getField(), endRow - startRow + 1, endColumn - startColumn + 1);

            // compute blocks shifts
            int blockStartRow = startRow / BLOCK_SIZE;
            int rowsShift = startRow % BLOCK_SIZE;
            int blockStartColumn = startColumn / BLOCK_SIZE;
            int columnsShift = startColumn % BLOCK_SIZE;

            // perform extraction block-wise, to ensure good cache behavior
            int pBlock = blockStartRow;
            for (int iBlock = 0; iBlock < outp.blockRows; ++iBlock)
            {
                int iHeight = outp.blockHeight(iBlock);
                int qBlock = blockStartColumn;
                for (int jBlock = 0; jBlock < outp.blockColumns; ++jBlock)
                {
                    int jWidth = outp.blockWidth(jBlock);

                    // handle one block of the output matrix
                    int outIndex = iBlock * outp.blockColumns + jBlock;
                    T[] outBlock = outp.blocks[outIndex];
                    int index = pBlock * blockColumns + qBlock;
                    int width = blockWidth(qBlock);

                    int heightExcess = iHeight + rowsShift - BLOCK_SIZE;
                    int widthExcess = jWidth + columnsShift - BLOCK_SIZE;
                    if (heightExcess > 0)
                    {
                        // the submatrix block spans on two blocks rows from the original matrix
                        if (widthExcess > 0)
                        {
                            // the submatrix block spans on two blocks columns from the original matrix
                            int width2 = blockWidth(qBlock + 1);
                            copyBlockPart(blocks[index], width,
                                          rowsShift, BLOCK_SIZE,
                                          columnsShift, BLOCK_SIZE,
                                          outBlock, jWidth, 0, 0);
                            copyBlockPart(blocks[index + 1], width2,
                                          rowsShift, BLOCK_SIZE,
                                          0, widthExcess,
                                          outBlock, jWidth, 0, jWidth - widthExcess);
                            copyBlockPart(blocks[index + blockColumns], width,
                                          0, heightExcess,
                                          columnsShift, BLOCK_SIZE,
                                          outBlock, jWidth, iHeight - heightExcess, 0);
                            copyBlockPart(blocks[index + blockColumns + 1], width2,
                                          0, heightExcess,
                                          0, widthExcess,
                                          outBlock, jWidth, iHeight - heightExcess, jWidth - widthExcess);
                        }
                        else
                        {
                            // the submatrix block spans on one block column from the original matrix
                            copyBlockPart(blocks[index], width,
                                          rowsShift, BLOCK_SIZE,
                                          columnsShift, jWidth + columnsShift,
                                          outBlock, jWidth, 0, 0);
                            copyBlockPart(blocks[index + blockColumns], width,
                                          0, heightExcess,
                                          columnsShift, jWidth + columnsShift,
                                          outBlock, jWidth, iHeight - heightExcess, 0);
                        }
                    }
                    else
                    {
                        // the submatrix block spans on one block row from the original matrix
                        if (widthExcess > 0)
                        {
                            // the submatrix block spans on two blocks columns from the original matrix
                            int width2 = blockWidth(qBlock + 1);
                            copyBlockPart(blocks[index], width,
                                          rowsShift, iHeight + rowsShift,
                                          columnsShift, BLOCK_SIZE,
                                          outBlock, jWidth, 0, 0);
                            copyBlockPart(blocks[index + 1], width2,
                                          rowsShift, iHeight + rowsShift,
                                          0, widthExcess,
                                          outBlock, jWidth, 0, jWidth - widthExcess);
                        }
                        else
                        {
                            // the submatrix block spans on one block column from the original matrix
                            copyBlockPart(blocks[index], width,
                                          rowsShift, iHeight + rowsShift,
                                          columnsShift, jWidth + columnsShift,
                                          outBlock, jWidth, 0, 0);
                        }
                    }
                    ++qBlock;
                }
                ++pBlock;
            }

            return outp;
        }

        /// <summary>
        /// Copy a part of a block into another one
        /// <para>This method can be called only when the specified part fits in both
        /// blocks, no verification is done here.</para> 
        /// </summary>
        /// <param name="srcBlock">source block</param>
        /// <param name="srcWidth">source block width (<see cref="BLOCK_SIZE"/> or smaller)</param>
        /// <param name="srcStartRow">start row in the source block</param>
        /// <param name="srcEndRow">end row (exclusive) in the source block</param>
        /// <param name="srcStartColumn">start column in the source block</param>
        /// <param name="srcEndColumn">end column (exclusive) in the source block</param>
        /// <param name="dstBlock">destination block</param>
        /// <param name="dstWidth">destination block width (<see cref="BLOCK_SIZE"/> or smaller)</param>
        /// <param name="dstStartRow">start row in the destination block</param>
        /// <param name="dstStartColumn">start column in the destination block</param>
        private void copyBlockPart(T[] srcBlock, int srcWidth, int srcStartRow, int srcEndRow, int srcStartColumn, int srcEndColumn, T[] dstBlock, int dstWidth, int dstStartRow, int dstStartColumn)
        {
            int length = srcEndColumn - srcStartColumn;
            int srcPos = srcStartRow * srcWidth + srcStartColumn;
            int dstPos = dstStartRow * dstWidth + dstStartColumn;
            for (int srcRow = srcStartRow; srcRow < srcEndRow; ++srcRow)
            {
                Array.Copy(srcBlock, srcPos, dstBlock, dstPos, length);
                srcPos += srcWidth;
                dstPos += dstWidth;
            }
        }

        /// <inheritdoc/>
        public new void setSubMatrix(T[][] subMatrix, int row, int column)
        {
            // safety checks
            MathUtils.checkNotNull(subMatrix);
            int refLength = subMatrix[0].Length;
            if (refLength == 0)
            {
                throw new NoDataException(new LocalizedFormats("AT_LEAST_ONE_COLUMN"));
            }
            int endRow = row + subMatrix.Length - 1;
            int endColumn = column + refLength - 1;
            checkSubMatrixIndex(row, endRow, column, endColumn);
            foreach (T[] subRow in subMatrix)
            {
                if (subRow.Length != refLength)
                {
                    throw new DimensionMismatchException(refLength, subRow.Length);
                }
            }

            // compute blocks bounds
            int blockStartRow = row / BLOCK_SIZE;
            int blockEndRow = (endRow + BLOCK_SIZE) / BLOCK_SIZE;
            int blockStartColumn = column / BLOCK_SIZE;
            int blockEndColumn = (endColumn + BLOCK_SIZE) / BLOCK_SIZE;

            // perform copy block-wise, to ensure good cache behavior
            for (int iBlock = blockStartRow; iBlock < blockEndRow; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                int firstRow = iBlock * BLOCK_SIZE;
                int iStart = FastMath.max(row, firstRow);
                int iEnd = FastMath.min(endRow + 1, firstRow + iHeight);

                for (int jBlock = blockStartColumn; jBlock < blockEndColumn; ++jBlock)
                {
                    int jWidth = blockWidth(jBlock);
                    int firstColumn = jBlock * BLOCK_SIZE;
                    int jStart = FastMath.max(column, firstColumn);
                    int jEnd = FastMath.min(endColumn + 1, firstColumn + jWidth);
                    int jLength = jEnd - jStart;

                    // handle one block, row by row
                    T[] block = blocks[iBlock * blockColumns + jBlock];
                    for (int i = iStart; i < iEnd; ++i)
                    {
                        Array.Copy(subMatrix[i - row], jStart - column,
                                         block, (i - firstRow) * jWidth + (jStart - firstColumn),
                                         jLength);
                    }

                }
            }
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> getRowMatrix(int row)
        {
            checkRowIndex(row);
            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), 1, columns);

            // perform copy block-wise, to ensure good cache behavior
            int iBlock = row / BLOCK_SIZE;
            int iRow = row - iBlock * BLOCK_SIZE;
            int outBlockIndex = 0;
            int outIndex = 0;
            T[] outBlock = outp.blocks[outBlockIndex];
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
            {
                int jWidth = blockWidth(jBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                int available = outBlock.Length - outIndex;
                if (jWidth > available)
                {
                    Array.Copy(block, iRow * jWidth, outBlock, outIndex, available);
                    outBlock = outp.blocks[++outBlockIndex];
                    Array.Copy(block, iRow * jWidth, outBlock, 0, jWidth - available);
                    outIndex = jWidth - available;
                }
                else
                {
                    Array.Copy(block, iRow * jWidth, outBlock, outIndex, jWidth);
                    outIndex += jWidth;
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public new void setRowMatrix(int row, FieldMatrix<T> matrix)
        {
            try
            {
                setRowMatrix(row, (BlockFieldMatrix<T>)matrix);
            }
            catch (InvalidCastException)
            {
                base.setRowMatrix(row, matrix);
            }
        }

        /// <summary>
        /// Sets the entries in row number <c>row</c>
        /// as a row matrix.  Row indices start at 0.
        /// </summary>
        /// <param name="row">the row to be set</param>
        /// <param name="matrix">row matrix (must have one row and the same number of columns
        /// as the instance)</param>
        /// <exception cref="MatrixDimensionMismatchException"> if the matrix dimensions do
        /// not match one instance row.</exception>
        /// <exception cref="OutOfRangeException"> if the specified row index is invalid.
        /// </exception>
        public void setRowMatrix(int row, BlockFieldMatrix<T> matrix)
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

            // perform copy block-wise, to ensure good cache behavior
            int iBlock = row / BLOCK_SIZE;
            int iRow = row - iBlock * BLOCK_SIZE;
            int mBlockIndex = 0;
            int mIndex = 0;
            T[] mBlock = matrix.blocks[mBlockIndex];
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
            {
                int jWidth = blockWidth(jBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                int available = mBlock.Length - mIndex;
                if (jWidth > available)
                {
                    Array.Copy(mBlock, mIndex, block, iRow * jWidth, available);
                    mBlock = matrix.blocks[++mBlockIndex];
                    Array.Copy(mBlock, 0, block, iRow * jWidth, jWidth - available);
                    mIndex = jWidth - available;
                }
                else
                {
                    Array.Copy(mBlock, mIndex, block, iRow * jWidth, jWidth);
                    mIndex += jWidth;
                }
            }
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> getColumnMatrix(int column)
        {
            checkColumnIndex(column);
            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), rows, 1);

            // perform copy block-wise, to ensure good cache behavior
            int jBlock = column / BLOCK_SIZE;
            int jColumn = column - jBlock * BLOCK_SIZE;
            int jWidth = blockWidth(jBlock);
            int outBlockIndex = 0;
            int outIndex = 0;
            T[] outBlock = outp.blocks[outBlockIndex];
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                for (int i = 0; i < iHeight; ++i)
                {
                    if (outIndex >= outBlock.Length)
                    {
                        outBlock = outp.blocks[++outBlockIndex];
                        outIndex = 0;
                    }
                    outBlock[outIndex++] = block[i * jWidth + jColumn];
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public new void setColumnMatrix(int column, FieldMatrix<T> matrix)
        {
            try
            {
                setColumnMatrix(column, (BlockFieldMatrix<T>)matrix);
            }
            catch (InvalidCastException)
            {
                base.setColumnMatrix(column, matrix);
            }
        }

        /// <summary>
        /// Sets the entries in column number <c>column</c>
        /// as a column matrix.  Column indices start at 0.
        /// </summary>
        /// <param name="column">Column to be set.</param>
        /// <param name="matrix">Column matrix (must have one column and the same number of rows
        /// as the instance).</param>
        /// <exception cref="MatrixDimensionMismatchException"> if the matrix dimensions do
        /// not match one instance column.</exception>
        /// <exception cref="OutOfRangeException"> if the specified column index is invalid
        /// </exception>
        void setColumnMatrix(int column, BlockFieldMatrix<T> matrix)
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

            // perform copy block-wise, to ensure good cache behavior
            int jBlock = column / BLOCK_SIZE;
            int jColumn = column - jBlock * BLOCK_SIZE;
            int jWidth = blockWidth(jBlock);
            int mBlockIndex = 0;
            int mIndex = 0;
            T[] mBlock = matrix.blocks[mBlockIndex];
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                for (int i = 0; i < iHeight; ++i)
                {
                    if (mIndex >= mBlock.Length)
                    {
                        mBlock = matrix.blocks[++mBlockIndex];
                        mIndex = 0;
                    }
                    block[i * jWidth + jColumn] = mBlock[mIndex++];
                }
            }
        }

        /// <inheritdoc/>
        public new FieldVector<T> getRowVector(int row)
        {
            checkRowIndex(row);
            T[] outData = MathArrays.buildArray(getField(), columns);

            // perform copy block-wise, to ensure good cache behavior
            int iBlock = row / BLOCK_SIZE;
            int iRow = row - iBlock * BLOCK_SIZE;
            int outIndex = 0;
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
            {
                int jWidth = blockWidth(jBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                Array.Copy(block, iRow * jWidth, outData, outIndex, jWidth);
                outIndex += jWidth;
            }

            return new ArrayFieldVector<T>(getField(), outData, false);
        }

        /// <inheritdoc/>
        public new void setRowVector(int row, FieldVector<T> vector)
        {
            try
            {
                setRow(row, ((ArrayFieldVector<T>)vector).getDataRef());
            }
            catch (InvalidCastException)
            {
                base.setRowVector(row, vector);
            }
        }

        /// <inheritdoc/>
        public new FieldVector<T> getColumnVector(int column)
        {
            checkColumnIndex(column);
            T[] outData = MathArrays.buildArray(getField(), rows);

            // perform copy block-wise, to ensure good cache behavior
            int jBlock = column / BLOCK_SIZE;
            int jColumn = column - jBlock * BLOCK_SIZE;
            int jWidth = blockWidth(jBlock);
            int outIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                for (int i = 0; i < iHeight; ++i)
                {
                    outData[outIndex++] = block[i * jWidth + jColumn];
                }
            }

            return new ArrayFieldVector<T>(getField(), outData, false);
        }

        /// <inheritdoc/>
        public new void setColumnVector(int column, FieldVector<T> vector)
        {
            try
            {
                setColumn(column, ((ArrayFieldVector<T>)vector).getDataRef());
            }
            catch (InvalidCastException)
            {
                base.setColumnVector(column, vector);
            }
        }

        /// <inheritdoc/>
        public new T[] getRow(int row)
        {
            checkRowIndex(row);
            T[] outp = MathArrays.buildArray(getField(), columns);

            // perform copy block-wise, to ensure good cache behavior
            int iBlock = row / BLOCK_SIZE;
            int iRow = row - iBlock * BLOCK_SIZE;
            int outIndex = 0;
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
            {
                int jWidth = blockWidth(jBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                Array.Copy(block, iRow * jWidth, outp, outIndex, jWidth);
                outIndex += jWidth;
            }

            return outp;
        }

        /// <inheritdoc/>
        public new void setRow(int row, T[] array)
        {
            checkRowIndex(row);
            int nCols = getColumnDimension();
            if (array.Length != nCols)
            {
                throw new MatrixDimensionMismatchException(1, array.Length, 1, nCols);
            }

            // perform copy block-wise, to ensure good cache behavior
            int iBlock = row / BLOCK_SIZE;
            int iRow = row - iBlock * BLOCK_SIZE;
            int outIndex = 0;
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
            {
                int jWidth = blockWidth(jBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                Array.Copy(array, outIndex, block, iRow * jWidth, jWidth);
                outIndex += jWidth;
            }
        }

        /// <inheritdoc/>
        public new T[] getColumn(int column)
        {
            checkColumnIndex(column);
            T[] outp = MathArrays.buildArray(getField(), rows);

            // perform copy block-wise, to ensure good cache behavior
            int jBlock = column / BLOCK_SIZE;
            int jColumn = column - jBlock * BLOCK_SIZE;
            int jWidth = blockWidth(jBlock);
            int outIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                for (int i = 0; i < iHeight; ++i)
                {
                    outp[outIndex++] = block[i * jWidth + jColumn];
                }
            }

            return outp;
        }

        /// <inheritdoc/>
        public new void setColumn(int column, T[] array)
        {
            checkColumnIndex(column);
            int nRows = getRowDimension();
            if (array.Length != nRows)
            {
                throw new MatrixDimensionMismatchException(array.Length, 1, nRows, 1);
            }

            // perform copy block-wise, to ensure good cache behavior
            int jBlock = column / BLOCK_SIZE;
            int jColumn = column - jBlock * BLOCK_SIZE;
            int jWidth = blockWidth(jBlock);
            int outIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int iHeight = blockHeight(iBlock);
                T[] block = blocks[iBlock * blockColumns + jBlock];
                for (int i = 0; i < iHeight; ++i)
                {
                    block[i * jWidth + jColumn] = array[outIndex++];
                }
            }
        }

        /// <inheritdoc/>
        public override T getEntry(int row, int column)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            int iBlock = row / BLOCK_SIZE;
            int jBlock = column / BLOCK_SIZE;
            int k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE);

            return blocks[iBlock * blockColumns + jBlock][k];
        }

        /// <inheritdoc/>
        public override void setEntry(int row, int column, T value)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            int iBlock = row / BLOCK_SIZE;
            int jBlock = column / BLOCK_SIZE;
            int k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE);

            blocks[iBlock * blockColumns + jBlock][k] = value;
        }

        /// <inheritdoc/>
        public override void addToEntry(int row, int column, T increment)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            int iBlock = row / BLOCK_SIZE;
            int jBlock = column / BLOCK_SIZE;
            int k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE);
            T[] blockIJ = blocks[iBlock * blockColumns + jBlock];

            blockIJ[k] = blockIJ[k].add(increment);
        }

        /// <inheritdoc/>
        public override void multiplyEntry(int row, int column, T factor)
        {
            checkRowIndex(row);
            checkColumnIndex(column);

            int iBlock = row / BLOCK_SIZE;
            int jBlock = column / BLOCK_SIZE;
            int k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE);
            T[] blockIJ = blocks[iBlock * blockColumns + jBlock];

            blockIJ[k] = blockIJ[k].multiply(factor);
        }

        /// <inheritdoc/>
        public new FieldMatrix<T> transpose()
        {
            int nRows = getRowDimension();
            int nCols = getColumnDimension();
            BlockFieldMatrix<T> outp = new BlockFieldMatrix<T>(getField(), nCols, nRows);

            // perform transpose block-wise, to ensure good cache behavior
            int blockIndex = 0;
            for (int iBlock = 0; iBlock < blockColumns; ++iBlock)
            {
                for (int jBlock = 0; jBlock < blockRows; ++jBlock)
                {

                    // transpose current block
                    T[] outBlock = outp.blocks[blockIndex];
                    T[] tBlock = blocks[jBlock * blockColumns + iBlock];
                    int pStart = iBlock * BLOCK_SIZE;
                    int pEnd = FastMath.min(pStart + BLOCK_SIZE, columns);
                    int qStart = jBlock * BLOCK_SIZE;
                    int qEnd = FastMath.min(qStart + BLOCK_SIZE, rows);
                    int k = 0;
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        int lInc = pEnd - pStart;
                        int l = p - pStart;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            outBlock[k] = tBlock[l];
                            ++k;
                            l += lInc;
                        }
                    }

                    // go to next block
                    ++blockIndex;

                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public override int getRowDimension()
        {
            return rows;
        }

        /// <inheritdoc/>
        public override int getColumnDimension()
        {
            return columns;
        }

        /// <inheritdoc/>
        public new T[] operate(T[] v)
        {
            if (v.Length != columns)
            {
                throw new DimensionMismatchException(v.Length, columns);
            }
            T[] outp = MathArrays.buildArray(getField(), rows);
            T zero = getField().getZero();

            // perform multiplication block-wise, to ensure good cache behavior
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                {
                    T[] block = blocks[iBlock * blockColumns + jBlock];
                    int qStart = jBlock * BLOCK_SIZE;
                    int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                    int k = 0;
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        T sum = zero;
                        int q = qStart;
                        while (q < qEnd - 3)
                        {
                            sum = sum.
                                  add(block[k].multiply(v[q])).
                                  add(block[k + 1].multiply(v[q + 1])).
                                  add(block[k + 2].multiply(v[q + 2])).
                                  add(block[k + 3].multiply(v[q + 3]));
                            k += 4;
                            q += 4;
                        }
                        while (q < qEnd)
                        {
                            sum = sum.add(block[k++].multiply(v[q++]));
                        }
                        outp[p] = outp[p].add(sum);
                    }
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public new T[] preMultiply(T[] v)
        {

            if (v.Length != rows)
            {
                throw new DimensionMismatchException(v.Length, rows);
            }
            T[] outp = MathArrays.buildArray(getField(), columns);
            T zero = getField().getZero();

            // perform multiplication block-wise, to ensure good cache behavior
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
            {
                int jWidth = blockWidth(jBlock);
                int jWidth2 = jWidth + jWidth;
                int jWidth3 = jWidth2 + jWidth;
                int jWidth4 = jWidth3 + jWidth;
                int qStart = jBlock * BLOCK_SIZE;
                int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                for (int iBlock = 0; iBlock < blockRows; ++iBlock)
                {
                    T[] block = blocks[iBlock * blockColumns + jBlock];
                    int pStart = iBlock * BLOCK_SIZE;
                    int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                    for (int q = qStart; q < qEnd; ++q)
                    {
                        int k = q - qStart;
                        T sum = zero;
                        int p = pStart;
                        while (p < pEnd - 3)
                        {
                            sum = sum.
                                  add(block[k].multiply(v[p])).
                                  add(block[k + jWidth].multiply(v[p + 1])).
                                  add(block[k + jWidth2].multiply(v[p + 2])).
                                  add(block[k + jWidth3].multiply(v[p + 3]));
                            k += jWidth4;
                            p += 4;
                        }
                        while (p < pEnd)
                        {
                            sum = sum.add(block[k].multiply(v[p++]));
                            k += jWidth;
                        }
                        outp[q] = outp[q].add(sum);
                    }
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                for (int p = pStart; p < pEnd; ++p)
                {
                    for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                    {
                        int jWidth = blockWidth(jBlock);
                        int qStart = jBlock * BLOCK_SIZE;
                        int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                        T[] block = blocks[iBlock * blockColumns + jBlock];
                        int k = (p - pStart) * jWidth;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            block[k] = visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor)
        {
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                for (int p = pStart; p < pEnd; ++p)
                {
                    for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                    {
                        int jWidth = blockWidth(jBlock);
                        int qStart = jBlock * BLOCK_SIZE;
                        int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                        T[] block = blocks[iBlock * blockColumns + jBlock];
                        int k = (p - pStart) * jWidth;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow,
                                int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(rows, columns, startRow, endRow, startColumn, endColumn);
            for (int iBlock = startRow / BLOCK_SIZE; iBlock < 1 + endRow / BLOCK_SIZE; ++iBlock)
            {
                int p0 = iBlock * BLOCK_SIZE;
                int pStart = FastMath.max(startRow, p0);
                int pEnd = FastMath.min((iBlock + 1) * BLOCK_SIZE, 1 + endRow);
                for (int p = pStart; p < pEnd; ++p)
                {
                    for (int jBlock = startColumn / BLOCK_SIZE; jBlock < 1 + endColumn / BLOCK_SIZE; ++jBlock)
                    {
                        int jWidth = blockWidth(jBlock);
                        int q0 = jBlock * BLOCK_SIZE;
                        int qStart = FastMath.max(startColumn, q0);
                        int qEnd = FastMath.min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn);
                        T[] block = blocks[iBlock * blockColumns + jBlock];
                        int k = (p - p0) * jWidth + qStart - q0;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            block[k] = visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInRowOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(rows, columns, startRow, endRow, startColumn, endColumn);
            for (int iBlock = startRow / BLOCK_SIZE; iBlock < 1 + endRow / BLOCK_SIZE; ++iBlock)
            {
                int p0 = iBlock * BLOCK_SIZE;
                int pStart = FastMath.max(startRow, p0);
                int pEnd = FastMath.min((iBlock + 1) * BLOCK_SIZE, 1 + endRow);
                for (int p = pStart; p < pEnd; ++p)
                {
                    for (int jBlock = startColumn / BLOCK_SIZE; jBlock < 1 + endColumn / BLOCK_SIZE; ++jBlock)
                    {
                        int jWidth = blockWidth(jBlock);
                        int q0 = jBlock * BLOCK_SIZE;
                        int qStart = FastMath.max(startColumn, q0);
                        int qEnd = FastMath.min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn);
                        T[] block = blocks[iBlock * blockColumns + jBlock];
                        int k = (p - p0) * jWidth + qStart - q0;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInOptimizedOrder(FieldMatrixChangingVisitor<T> visitor)
        {
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            int blockIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                {
                    int qStart = jBlock * BLOCK_SIZE;
                    int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                    T[] block = blocks[blockIndex];
                    int k = 0;
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            block[k] = visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                    ++blockIndex;
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInOptimizedOrder(FieldMatrixPreservingVisitor<T> visitor)
        {
            visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
            int blockIndex = 0;
            for (int iBlock = 0; iBlock < blockRows; ++iBlock)
            {
                int pStart = iBlock * BLOCK_SIZE;
                int pEnd = FastMath.min(pStart + BLOCK_SIZE, rows);
                for (int jBlock = 0; jBlock < blockColumns; ++jBlock)
                {
                    int qStart = jBlock * BLOCK_SIZE;
                    int qEnd = FastMath.min(qStart + BLOCK_SIZE, columns);
                    T[] block = blocks[blockIndex];
                    int k = 0;
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                    ++blockIndex;
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInOptimizedOrder(FieldMatrixChangingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(rows, columns, startRow, endRow, startColumn, endColumn);
            for (int iBlock = startRow / BLOCK_SIZE; iBlock < 1 + endRow / BLOCK_SIZE; ++iBlock)
            {
                int p0 = iBlock * BLOCK_SIZE;
                int pStart = FastMath.max(startRow, p0);
                int pEnd = FastMath.min((iBlock + 1) * BLOCK_SIZE, 1 + endRow);
                for (int jBlock = startColumn / BLOCK_SIZE; jBlock < 1 + endColumn / BLOCK_SIZE; ++jBlock)
                {
                    int jWidth = blockWidth(jBlock);
                    int q0 = jBlock * BLOCK_SIZE;
                    int qStart = FastMath.max(startColumn, q0);
                    int qEnd = FastMath.min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn);
                    T[] block = blocks[iBlock * blockColumns + jBlock];
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        int k = (p - p0) * jWidth + qStart - q0;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            block[k] = visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                }
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new T walkInOptimizedOrder(FieldMatrixPreservingVisitor<T> visitor, int startRow, int endRow, int startColumn, int endColumn)
        {
            checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
            visitor.start(rows, columns, startRow, endRow, startColumn, endColumn);
            for (int iBlock = startRow / BLOCK_SIZE; iBlock < 1 + endRow / BLOCK_SIZE; ++iBlock)
            {
                int p0 = iBlock * BLOCK_SIZE;
                int pStart = FastMath.max(startRow, p0);
                int pEnd = FastMath.min((iBlock + 1) * BLOCK_SIZE, 1 + endRow);
                for (int jBlock = startColumn / BLOCK_SIZE; jBlock < 1 + endColumn / BLOCK_SIZE; ++jBlock)
                {
                    int jWidth = blockWidth(jBlock);
                    int q0 = jBlock * BLOCK_SIZE;
                    int qStart = FastMath.max(startColumn, q0);
                    int qEnd = FastMath.min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn);
                    T[] block = blocks[iBlock * blockColumns + jBlock];
                    for (int p = pStart; p < pEnd; ++p)
                    {
                        int k = (p - p0) * jWidth + qStart - q0;
                        for (int q = qStart; q < qEnd; ++q)
                        {
                            visitor.visit(p, q, block[k]);
                            ++k;
                        }
                    }
                }
            }
            return visitor.end();
        }

        /// <summary>
        /// Get the height of a block.
        /// </summary>
        /// <param name="blockRow">row index (in block sense) of the block</param>
        /// <returns>height (number of rows) of the block</returns>
        private int blockHeight(int blockRow)
        {
            return (blockRow == blockRows - 1) ? rows - blockRow * BLOCK_SIZE : BLOCK_SIZE;
        }

        /// <summary>
        /// Get the width of a block.
        /// </summary>
        /// <param name="blockColumn">column index (in block sense) of the block</param>
        /// <returns>width (number of columns) of the block</returns>
        private int blockWidth(int blockColumn)
        {
            return (blockColumn == blockColumns - 1) ? columns - blockColumn * BLOCK_SIZE : BLOCK_SIZE;
        }
    }
}