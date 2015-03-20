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
using Math3.exception;
using Math3.exception.util;

namespace Math3.fraction
{
    /// <summary>
    /// Error thrown when a double value cannot be converted to a fraction
    /// in the allowed number of iterations.
    /// </summary>
    public class FractionConversionException : ConvergenceException
    {
        /// <summary>
        /// Constructs an exception with specified formatted detail message. 
        /// </summary>
        /// <param name="value">double value to convert</param>
        /// <param name="maxIterations">maximal number of iterations allowed</param>
        public FractionConversionException(double value, int maxIterations) : base(new LocalizedFormats("FAILED_FRACTION_CONVERSION"), value, maxIterations) { }

        /// <summary>
        /// Constructs an exception with specified formatted detail message.
        /// </summary>
        /// <param name="value">double value to convert</param>
        /// <param name="p">current numerator</param>
        /// <param name="q">current denominator</param>
        public FractionConversionException(double value, long p, long q) : base(new LocalizedFormats("FRACTION_CONVERSION_OVERFLOW"), value, p, q) { }
    }
}
