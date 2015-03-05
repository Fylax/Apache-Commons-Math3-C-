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
    /// Base class for exceptions raised by a wrong number.
    /// This class is not intended to be instantiated directly: it should serve
    /// as a base class to create all the exceptions that are raised because some
    /// precondition is violated by a number argument.
    /// </summary>
    /// <typeparam name="T">A number.</typeparam>
    public class MathIllegalNumberException<T> : MathIllegalArgumentException where T : IComparable
    {
        /// <summary>
        /// Helper to avoid boxing warnings
        /// </summary>
        protected const Int32 INTEGER_ZERO = 0;

        /// <summary>
        /// Requested.
        /// </summary>
        protected readonly T argument;

        /// <summary>
        /// Construct an exception. 
        /// </summary>
        /// <param name="pattern">Localizable pattern.</param>
        /// <param name="wrong">Wrong number.</param>
        /// <param name="arguments">Arguments.</param>
        protected MathIllegalNumberException(Localizable pattern, T wrong, params Object[] arguments) : base(pattern, wrong, arguments)
        {
            argument = wrong;
        }

        /// <summary></summary>
        /// <returns>the requested value.</returns>
        public T getArgument()
        {
            return argument;
        }
    }
}
