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
    /// Interface defining a visitor for matrix entries.
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public interface FieldMatrixPreservingVisitor<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Start visiting a matrix.
        /// <para>This method is called once before any entry of the matrix is visited.</para>
        /// </summary>
        /// <param name="rows">number of rows of the matrix</param>
        /// <param name="columns">number of columns of the matrix</param>
        /// <param name="startRow">Initial row index</param>
        /// <param name="endRow">Final row index (inclusive)</param>
        /// <param name="startColumn">Initial column index</param>
        /// <param name="endColumn">Final column index (inclusive)</param>
        void start(int rows, int columns, int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>
        /// Visit one matrix entry.
        /// </summary>
        /// <param name="row">row index of the entry</param>
        /// <param name="column">column index of the entry</param>
        /// <param name="value">current value of the entry</param>
        /// <returns>the new value to be set for the entry</returns>
        void visit(int row, int column, T value);

        /// <summary>
        /// End visiting a matrix.
        /// <para>This method is called once after all entries of the matrix have been visited.
        /// </para>
        /// </summary>
        /// <returns>the value that the <c>walkInXxxOrder</c> must return</returns>
        T end();
    }
}
