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
using System;
using System.Runtime.InteropServices;

namespace Math3.linear
{
    /// <summary>
    /// Exception to be thrown when a positive definite matrix is expected.
    /// </summary>
    [ComVisible(false)]
    public class NonPositiveDefiniteMatrixException : NumberIsTooSmallException<Double, Double>
    {
        /// <summary>
        /// Index (diagonal element).
        /// </summary>
        private readonly int index;
        
        /// <summary>
        /// Threshold.
        /// </summary>
        private readonly double threshold;

        /// <summary>
        /// Construct an exception.
        /// </summary>
        /// <param name="wrong">Value that fails the positivity check.</param>
        /// <param name="index">Row (and column) index.</param>
        /// <param name="threshold">Absolute positivity threshold.</param>
        public NonPositiveDefiniteMatrixException(double wrong, int index, double threshold)
            : base(wrong, threshold, false)
        {
            this.index = index;
            this.threshold = threshold;

            ExceptionContext context = getContext();
            context.addMessage(new LocalizedFormats("NOT_POSITIVE_DEFINITE_MATRIX"));
            context.addMessage(new LocalizedFormats("ARRAY_ELEMENT"), wrong, index);
        }

        /// <summary></summary>
        /// <returns>the row index.</returns>
        public int getRow()
        {
            return index;
        }
        
        /// <summary></summary>
        /// <returns>the column index.</returns>
        public int getColumn()
        {
            return index;
        }
        
        /// <summary></summary>
        /// <returns>the absolute positivity threshold.</returns>
        public double getThreshold()
        {
            return threshold;
        }
    }
}