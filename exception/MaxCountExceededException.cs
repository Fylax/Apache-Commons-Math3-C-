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
    /// Exception to be thrown when some counter maximum value is exceeded.
    /// </summary>
    /// <typeparam name="T">A number</typeparam>
    public class MaxCountExceededException<T> : MathIllegalStateException where T : IComparable
    {

        /// <summary>
        /// Maximum number of evaluations.
        /// </summary>
        private readonly T max;

        /// <summary>
        /// Construct the exception.
        /// </summary>
        /// <param name="max">Maximum.</param>
        public MaxCountExceededException(T max) : this(new LocalizedFormats("MAX_COUNT_EXCEEDED"), max) { }

        /// <summary>
        /// Construct the exception with a specific context.
        /// </summary>
        /// <param name="specific">Specific context pattern.</param>
        /// <param name="max">Maximum.</param>
        /// <param name="args">Additional arguments.</param>
        public MaxCountExceededException(Localizable specific, T max, params Object[] args)
        {
            getContext().addMessage(specific, max, args);
            this.max = max;
        }

        /// <summary></summary>
        /// <returns>the maximum number of evaluations.</returns>
        public T getMax()
        {
            return max;
        }
    }
}
