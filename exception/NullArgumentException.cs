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
// limitations under the License
using Math3.exception.util;
using System;

namespace Math3.exception
{
    /// <summary>
    /// All conditions checks that fail due to a <c>null</c> argument must throw
    /// this exception.
    /// This class is meant to signal a precondition violation ("null is an illegal
    /// argument") and so does not extend the standard <c>NullReferenceException</c>.
    /// Propagation of <c>NullReferenceException</c> from within Commons-Math is
    /// construed to be a bug.
    /// </summary>
    public class NullArgumentException : MathIllegalArgumentException
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public NullArgumentException() : this(new LocalizedFormats("NULL_NOT_ALLOWED")) { }
        
        /// <summary></summary>
        /// <param name="pattern">Message pattern providing the specific context of
        /// the error.</param>
        /// <param name="arguments">Values for replacing the placeholders in <c>pattern</c>.</param>
        public NullArgumentException(Localizable pattern, params Object[] arguments) : base(pattern, arguments) { }
    }
}
