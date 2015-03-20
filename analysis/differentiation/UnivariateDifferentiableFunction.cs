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
namespace Math3.analysis.differentiation
{
    /// <summary>
    /// Interface for univariate functions derivatives.
    /// <para>This interface represents a simple function which computes
    /// both the value and the first derivative of a mathematical function.
    /// The derivative is computed with respect to the input variable.</para>
    /// </summary>
    /// <remarks>
    /// See <seealso cref="UnivariateDifferentiableFunction"/><para/>
    /// See <seealso cref="UnivariateFunctionDifferentiator"/>
    /// </remarks>
    public interface UnivariateDifferentiableFunction : UnivariateFunction
    {
        /// <summary>
        /// Simple mathematical function.
        /// <para><see cref="UnivariateDifferentiableFunction"/> classes compute both the
        /// value and the first derivative of the function.</para>
        /// </summary>
        /// <param name="t">function input value</param>
        /// <returns>function result</returns>
        /// <exception cref="DimensionMismatchException"> if t is inconsistent with the
        /// function's free parameters or order</exception>
        DerivativeStructure value(DerivativeStructure t);
    }
}
