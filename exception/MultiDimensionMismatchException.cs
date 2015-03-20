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
using Math3.exception.util;
using System;

namespace Math3.exception
{
    /// <summary>
    /// Exception to be thrown when two sets of dimensions differ.
    /// </summary>
    public class MultiDimensionMismatchException : MathIllegalArgumentException
    {
        /// <summary>
        /// Wrong dimensions.
        /// </summary>
        private readonly Int32[] wrong;
        
        /// <summary>
        /// Correct dimensions.
        /// </summary>
        private readonly Int32[] expected;

        /// <summary>
        /// Construct an exception from the mismatched dimensions.
        /// </summary>
        /// <param name="wrong">Wrong dimensions.</param>
        /// <param name="expected">Expected dimensions.</param>
        public MultiDimensionMismatchException(Int32[] wrong, Int32[] expected) : this(new LocalizedFormats("DIMENSIONS_MISMATCH"), wrong, expected) { }

        /// <summary>
        /// Construct an exception from the mismatched dimensions.
        /// </summary>
        /// <param name="specific">specific Message pattern providing the specific context of
        /// the error.</param>
        /// <param name="wrong">Wrong dimensions.</param>
        /// <param name="expected">Expected dimensions.</param>
        public MultiDimensionMismatchException(Localizable specific, Int32[] wrong, Int32[] expected) : base(specific, wrong, expected)
        {
            this.wrong = (Int32[])wrong.Clone();
            this.expected = (Int32[])expected.Clone();
        }

        /// <summary></summary>
        /// <returns>an array containing the wrong dimensions.</returns>
        public Int32[] getWrongDimensions()
        {
            return (Int32[])wrong.Clone();
        }
        
        /// <summary></summary>
        /// <returns>an array containing the expected dimensions.</returns>
        public Int32[] getExpectedDimensions()
        {
            return (Int32[])expected.Clone();
        }

        /// <summary></summary>
        /// <param name="index">Dimension index.</param>
        /// <returns>the wrong dimension stored at <c>index</c>.</returns>
        public int getWrongDimension(int index)
        {
            return wrong[index];
        }
        
        /// <summary></summary>
        /// <param name="index">Dimension index.</param>
        /// <returns>the expected dimension stored at <c>index</c>.</returns>
        public int getExpectedDimension(int index)
        {
            return expected[index];
        }
    }
}
