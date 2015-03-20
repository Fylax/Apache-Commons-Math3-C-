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
using System.Runtime.InteropServices;

namespace Math3.linear
{
    /// <summary>
    /// Exception to be thrown when a square matrix is expected.
    /// </summary>
    [ComVisible(false)]
    public class NonSquareMatrixException : DimensionMismatchException
    {
        /// <summary>
        /// Construct an exception from the mismatched dimensions.
        /// </summary>
        /// <param name="wrong">Row dimension.</param>
        /// <param name="expected">Column dimension.</param>
        public NonSquareMatrixException(int wrong, int expected) : base(new LocalizedFormats("NON_SQUARE_MATRIX"), wrong, expected) { }
    }
}
