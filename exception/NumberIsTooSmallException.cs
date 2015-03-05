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
    /// Exception to be thrown when a number is too small.
    /// </summary>
    /// <typeparam name="T">A Number</typeparam>
    /// <typeparam name="U">A Number</typeparam>
    public class NumberIsTooSmallException<T, U> : MathIllegalNumberException<T>
        where T : IComparable
        where U : IComparable
    {

        /// <summary>
        /// Higher bound.
        /// </summary>
        private readonly U min;
        
        /// <summary>
        /// Whether the maximum is included in the allowed range.
        /// </summary>
        private readonly Boolean boundIsAllowed;

        /// <summary>
        /// Construct the exception.
        /// </summary>
        /// <param name="wrong">Value that is smaller than the minimum.</param>
        /// <param name="min"></param>
        /// <param name="boundMinimum.IsAllowed">boundIsAllowed Whether <c>min</c> is included in the allowed range.</param>
        public NumberIsTooSmallException(T wrong, U min, Boolean boundIsAllowed) : this(boundIsAllowed ? new LocalizedFormats("NUMBER_TOO_SMALL") : new LocalizedFormats("NUMBER_TOO_SMALL_BOUND_EXCLUDED"), wrong, min, boundIsAllowed) { }

        /// <summary>
        /// Construct the exception with a specific context. 
        /// </summary>
        /// <param name="specific">Specific context pattern.</param>
        /// <param name="wrong">Value that is smaller than the minimum.</param>
        /// <param name="min">Minimum.</param>
        /// <param name="boundIsAllowed">boundIsAllowed Whether <c>min</c> is included in the allowed range.</param>
        public NumberIsTooSmallException(Localizable specific, T wrong, U min, Boolean boundIsAllowed) : base(specific, wrong, min)
        {
            this.min = min;
            this.boundIsAllowed = boundIsAllowed;
        }

        /// <summary></summary>
        /// <returns><c>true</c> if the minimum is included in the allowed range.</returns>
        public Boolean getBoundIsAllowed()
        {
            return boundIsAllowed;
        }

        /// <summary></summary>
        /// <returns>the minimum.</returns>
        public U getMin()
        {
            return min;
        }
    }

}
