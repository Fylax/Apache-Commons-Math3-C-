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
    /// Default implementation of the <see cref="FieldMatrixPreservingVisitor"/> interface.
    /// <para>
    /// This class is a convenience to create custom visitors without defining all
    /// methods. This class provides default implementations that do nothing.
    /// </para>
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public class DefaultFieldMatrixPreservingVisitor<T> : FieldMatrixPreservingVisitor<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Zero element of the field.
        /// </summary>
        private readonly T zero;

        /// <summary>
        /// Build a new instance.
        /// </summary>
        /// <param name="zero">additive identity of the field</param>
        public DefaultFieldMatrixPreservingVisitor(T zero)
        {
            this.zero = zero;
        }

        /// <inheritdoc/>
        public void start(int rows, int columns,
                          int startRow, int endRow, int startColumn, int endColumn)
        {
        }

        /// <inheritdoc/>
        public void visit(int row, int column, T value) { }

        /// <inheritdoc/>
        public T end()
        {
            return zero;
        }
    }
}
