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
    /// Error thrown when a numerical computation can not be performed because the
    /// numerical result failed to converge to a finite value.
    /// </summary>
    public class ConvergenceException : MathIllegalStateException
    {
        /// <summary>
        /// Construct the exception.
        /// </summary>
        public ConvergenceException() : this(new LocalizedFormats("CONVERGENCE_FAILED")) { }

        /// <summary>
        /// Construct the exception with a specific context and arguments.
        /// </summary>
        /// <param name="pattern">Message pattern providing the specific context of
        /// the error.</param>
        /// <param name="args">Arguments.</param>
        public ConvergenceException(Localizable pattern, params Object[] args)
        {
            getContext().addMessage(pattern, args);
        }
    }
}
