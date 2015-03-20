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
using System;

namespace Math3.fraction
{
    /// <summary>
    /// Representation of the fractional numbers field.
    /// <para>
    /// This class is a singleton.
    /// </para>
    /// </summary>
    /// <remarks>
    /// See <see cref="Fraction"/>
    /// </remarks>
    public class FractionField : Field<Fraction>
    {
        /// <summary>
        /// Private constructor for the singleton.
        /// </summary>
        private FractionField() { }

        /// <summary>
        /// Get the unique instance.
        /// </summary>
        /// <returns>the unique instance</returns>
        public static FractionField getInstance()
        {
            return LazyHolder.INSTANCE;
        }

        /// <inheritdoc/>
        public Fraction getOne()
        {
            return Fraction.ONE;
        }

        /// <inheritdoc/>
        public Fraction getZero()
        {
            return Fraction.ZERO;
        }

        /// <inheritdoc/>
        public Type getRuntimeClass()
        {
            return this.GetType();
        }
        // CHECKSTYLE: stop HideUtilityClassConstructor
        
        /// <summary>
        /// Holder for the instance.
        /// <para>We use here the Initialization On Demand Holder Idiom.</para>
        /// </summary>
        private static class LazyHolder
        {
            /// <summary>
            /// Cached field instance.
            /// </summary>
            internal static readonly FractionField INSTANCE = new FractionField();
        }
        // CHECKSTYLE: resume HideUtilityClassConstructor

        /// <summary>
        /// Handle deserialization of the singleton.
        /// </summary>
        /// <returns>the singleton instance</returns>
        private Object readResolve()
        {
            // return the singleton instance
            return LazyHolder.INSTANCE;
        }
    }
}