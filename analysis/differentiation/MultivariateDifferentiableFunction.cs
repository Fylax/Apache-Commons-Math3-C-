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
    /// Extension of <see cref="MultivariateFunction"/> representing a
    /// multivariate differentiable real function.
    /// </summary>
    public interface MultivariateDifferentiableFunction : MultivariateFunction
    {
        /// <summary>
        /// Compute the value for the function at the given point.
        /// </summary>
        /// <param name="point">Point at which the function must be evaluated.</param>
        /// <returns>the function value for the given point.</returns>
        /// <exception cref="MathIllegalArgumentException"> if <c>point</c> does not
        /// satisfy the function's constraints (wrong dimension, argument out of bound,
        /// or unsupported derivative order for example)</exception>
        DerivativeStructure value(DerivativeStructure[] point);
    }
}