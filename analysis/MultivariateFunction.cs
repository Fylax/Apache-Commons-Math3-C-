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
namespace Math3.analysis
{
    /// <summary>
    /// An interface representing a multivariate real function.
    /// </summary>
    public interface MultivariateFunction
    {
        /// <summary>
        /// Compute the value for the function at the given point.
        /// </summary>
        /// <param name="point">Point at which the function must be evaluated.</param>
        /// <returns>function value for the given point.</returns>
        /// <exception cref="DimensionMismatchException"> if the parameter's dimension is
        /// wrong for the function being evaluated.</exception>
        /// <exception cref="MathIllegalArgumentException"> when the activated method itself 
        /// can ascertain that preconditions, specified in the API expressed at the level
        /// of the activated method, have been violated.  In the vast majority of cases 
        /// where Commons Math throws this exception, it is the result of argument checking 
        /// of actual parameters immediately passed to a method.</exception>
        double value(double[] point);
    }
}