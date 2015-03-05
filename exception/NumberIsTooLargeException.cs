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
    /// Exception to be thrown when a number is too large.
    /// </summary>
    /// <typeparam name="T">A number</typeparam>
    /// <typeparam name="U">A number</typeparam>
    public class NumberIsTooLargeException<T, U> : MathIllegalNumberException<T>
        where T : IComparable
        where U : IComparable
    {
        /// <summary>
        /// Higher bound.
        /// </summary>
        private readonly U max;
        
        /// <summary>
        /// Whether the maximum is included in the allowed range.
        /// </summary>
        private readonly Boolean boundIsAllowed;

        /// <summary>
        /// Construct the exception.
        /// </summary>
        /// <param name="wrong">Value that is larger than the maximum.</param>
        /// <param name="max">Maximum.</param>
        /// <param name="boundIsAllowed">if true the maximum is included in the allowed range.</param>
        public NumberIsTooLargeException(T wrong, U max, Boolean boundIsAllowed) : this(boundIsAllowed ? new LocalizedFormats("NUMBER_TOO_LARGE") : new LocalizedFormats("NUMBER_TOO_LARGE_BOUND_EXCLUDED"), wrong, max, boundIsAllowed) { }
        
        /// <summary>
        /// Construct the exception with a specific context.
        /// </summary>
        /// <param name="specific">Specific context pattern.</param>
        /// <param name="wrong">Value that is larger than the maximum.</param>
        /// <param name="max">Maximum.</param>
        /// <param name="boundIsAllowed">if true the maximum is included in the allowed range.</param>
        public NumberIsTooLargeException(Localizable specific, T wrong, U max, Boolean boundIsAllowed) : base(specific, wrong, max)
        {
            this.max = max;
            this.boundIsAllowed = boundIsAllowed;
        }

        /// <summary>
        /// <c>true</c> if the maximum is included in the allowed range.
        /// </summary>
        /// <returns></returns>
        public Boolean getBoundIsAllowed()
        {
            return boundIsAllowed;
        }

        /// <summary></summary>
        /// <returns>the maximum.</returns>
        public U getMax()
        {
            return max;
        }
    }
}
