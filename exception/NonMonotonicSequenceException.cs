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
using Math3.util;
using System;

namespace Math3.exception
{
    /// <summary>
    /// Exception to be thrown when the a sequence of values is not monotonically
    /// increasing or decreasing.
    /// </summary>
    /// <typeparam name="T">A number</typeparam>
    /// <typeparam name="U">A number</typeparam>
    public class NonMonotonicSequenceException<T, U> : MathIllegalNumberException<U>
        where T : IComparable
        where U : IComparable
    {
        /// <summary>
        /// Direction (positive for increasing, negative for decreasing).
        /// </summary>
        private readonly MathArrays.OrderDirection direction;
        
        /// <summary>
        /// Whether the sequence must be strictly increasing or decreasing.
        /// </summary>
        private readonly Boolean strict;
        
        /// <summary>
        /// Index of the wrong value.
        /// </summary>
        private readonly int index;
       
        /// <summary>
        /// Previous value.
        /// </summary>
        private readonly T previous;

        /// <summary>
        /// Construct the exception.
        /// This constructor uses default values assuming that the sequence should
        /// have been strictly increasing.
        /// </summary>
        /// <param name="wrong">Value that did not match the requirements.</param>
        /// <param name="previous">Previous value in the sequence.</param>
        /// <param name="index">Index of the value that did not match the requirements.</param>
        public NonMonotonicSequenceException(U wrong, T previous, int index) :
            this(wrong, previous, index, MathArrays.OrderDirection.INCREASING, true) { }

        /// <summary>
        /// Construct the exception. 
        /// </summary>
        /// <param name="wrong">Value that did not match the requirements.</param>
        /// <param name="previous">Previous value in the sequence.</param>
        /// <param name="index">Index of the value that did not match the requirements.</param>
        /// <param name="direction">Strictly positive for a sequence required to be
        /// increasing, negative (or zero) for a decreasing sequence.</param>
        /// <param name="strict">Whether the sequence must be strictly increasing or
        /// decreasing.</param>
        public NonMonotonicSequenceException(U wrong, T previous, int index, MathArrays.OrderDirection direction, Boolean strict) :
            base((direction == MathArrays.OrderDirection.INCREASING) ?
                  (strict ? new LocalizedFormats("NOT_STRICTLY_INCREASING_SEQUENCE") :
                   new LocalizedFormats("NOT_INCREASING_SEQUENCE")) :
                  (strict ? new LocalizedFormats("NOT_STRICTLY_DECREASING_SEQUENCE") :
                   new LocalizedFormats("NOT_DECREASING_SEQUENCE")),
                  wrong, previous, index, index - 1)
        {
            this.direction = direction;
            this.strict = strict;
            this.index = index;
            this.previous = previous;
        }

        /// <summary></summary>
        /// <returns>the order direction.</returns>
        public MathArrays.OrderDirection getDirection()
        {
            return direction;
        }
        
        /// <summary></summary>
        /// <returns><c>true</c> is the sequence should be strictly monotonic.</returns>
        public Boolean getStrict()
        {
            return strict;
        }
        
        /// <summary>
        /// Get the index of the wrong value.
        /// </summary>
        /// <returns>the current index.</returns>
        public int getIndex()
        {
            return index;
        }
        
        /// <summary></summary>
        /// <returns>the previous value.</returns>
        public T getPrevious()
        {
            return previous;
        }
    }
}
