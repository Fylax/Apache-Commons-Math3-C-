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

namespace Math3.linear
{
    /// <summary>
    /// Exception to be thrown when a symmetric matrix is expected.
    /// </summary>
    public class NonSymmetricMatrixException : MathIllegalArgumentException
    {
        /// <summary>
        /// Row.
        /// </summary>
        private readonly int row;
        
        /// <summary>
        /// Column.
        /// </summary>
        private readonly int column;
        
        /// <summary>
        /// Threshold.
        /// </summary>
        private readonly double threshold;

        /// <summary>
        /// Construct an exception.
        /// </summary>
        /// <param name="row">Row index.</param>
        /// <param name="column">Column index.</param>
        /// <param name="threshold">Relative symmetry threshold.</param>
        public NonSymmetricMatrixException(int row, int column, double threshold)
            : base(new LocalizedFormats("NON_SYMMETRIC_MATRIX"), row, column, threshold)
        {
            this.row = row;
            this.column = column;
            this.threshold = threshold;
        }

        /// <returns>the row index of the entry.</returns>
        public int getRow()
        {
            return row;
        }
        
        /// <returns>the column index of the entry.</returns>
        public int getColumn()
        {
            return column;
        }
        
        /// <returns>the relative symmetry threshold.</returns>
        public double getThreshold()
        {
            return threshold;
        }
    }
}
