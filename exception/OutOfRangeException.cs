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
    /// Exception to be thrown when some argument is out of range.
    /// </summary>
    /// <typeparam name="T">A Number</typeparam>
    public class OutOfRangeException<T> : MathIllegalNumberException<T> where T : IComparable
    {

        /// <summary>
        /// Lower bound.
        /// </summary>
        private readonly T lo;
        
        /// <summary>
        /// Higher bound.
        /// </summary>
        private readonly T hi;

        /// <summary>
        /// Construct an exception from the mismatched dimensions.
        /// </summary>
        /// <param name="wrong">Requested value.</param>
        /// <param name="lo">Lower bound.</param>
        /// <param name="hi">Higher bound.</param>
        public OutOfRangeException(T wrong, T lo, T hi) : this(new LocalizedFormats("OUT_OF_RANGE_SIMPLE"), wrong, lo, hi) { }

        /// <summary>
        /// Construct an exception from the mismatched dimensions with a
        /// specific context information.
        /// </summary>
        /// <param name="specific">Context information.</param>
        /// <param name="wrong">Requested value.</param>
        /// <param name="lo">Lower bound.</param>
        /// <param name="hi">Higher bound.</param>
        public OutOfRangeException(Localizable specific, T wrong, T lo, T hi) : base(specific, wrong, lo, hi)
        {
            this.lo = lo;
            this.hi = hi;
        }

        /// <summary></summary>
        /// <returns>the lower bound.</returns>
        public T getLo()
        {
            return lo;
        }
        
        /// <summary></summary>
        /// <returns>the higher bound.</returns>
        public T getHi()
        {
            return hi;
        }
    }
}
