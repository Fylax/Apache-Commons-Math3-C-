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
namespace Math3.analysis.solvers
{
    /// <summary>
    /// Interface for <see cref="UnivariateSolver (univariate real)"> root-finding
    /// algorithms</see> that maintain a bracketed solution. There are several advantages
    /// to having such root-finding algorithms:
    /// <list type="bullet">
    /// <item>The bracketed solution guarantees that the root is kept within the
    /// interval. As such, these algorithms generally also guarantee
    /// convergence.</item>
    /// <item>The bracketed solution means that we have the opportunity to only
    /// return roots that are greater than or equal to the actual root, or
    /// are less than or equal to the actual root. That is, we can control
    /// whether under-approximations and over-approximations are
    /// <see cref="AllowedSolution allowed solutions"/>. Other root-finding
    /// algorithms can usually only guarantee that the solution (the root that
    /// was found) is around the actual root.</item>
    /// </list>
    /// <para>For backwards compatibility, all root-finding algorithms must have
    /// <see cref="AllowedSolution#ANY_SIDE ANY_SIDE"/> as default for the allowed
    /// solutions.</para>
    /// </summary>
    /// <typeparam name="FUNC">Type of function to solve.</typeparam>
    /// <remarks>
    /// See <see cref="AllowedSolution"/>
    /// </remarks>
    public interface BracketedUnivariateSolver<FUNC> : BaseUnivariateSolver<FUNC> where FUNC : UnivariateFunction
    {
        /// <summary>
        /// Solve for a zero in the given interval.
        /// A solver may require that the interval brackets a single zero root.
        /// Solvers that do require bracketing should be able to handle the case
        /// where one of the endpoints is itself a root
        /// </summary>
        /// <param name="maxEval">Maximum number of evaluations.</param>
        /// <param name="f">Function to solve.</param>
        /// <param name="min">Lower bound for the interval.</param>
        /// <param name="max">Upper bound for the interval.</param>
        /// <param name="allowedSolution">The kind of solutions that the root-finding
        /// algorithm may accept as solutions.</param>
        /// <returns>A value where the function is zero.</returns>
        /// <exception cref="MathIllegalArgumentException"> if the arguments
        /// do not satisfy the requirements specified by the solver.</exception>
        /// <exception cref="TooManyEvaluationsException"> if the allowed number 
        /// of evaluations is exceeded.</exception>
        double solve(int maxEval, FUNC f, double min, double max, AllowedSolution allowedSolution);

        /// <summary>
        /// Solve for a zero in the given interval, start at <c>startValue</c>.
        /// A solver may require that the interval brackets a single zero root.
        /// Solvers that do require bracketing should be able to handle the case
        /// where one of the endpoints is itself a root.
        /// </summary>
        /// <param name="maxEval">Maximum number of evaluations.</param>
        /// <param name="f">Function to solve.</param>
        /// <param name="min">Lower bound for the interval.</param>
        /// <param name="max">Upper bound for the interval.</param>
        /// <param name="startValue">Start value to use.</param>
        /// <param name="allowedSolution">The kind of solutions that the root-finding 
        /// algorithm may accept as solutions.</param>
        /// <returns>A value where the function is zero.</returns>
        /// <exception cref="MathIllegalArgumentException"> if the arguments
        /// do not satisfy the requirements specified by the solver.</exception>
        /// <exception cref="TooManyEvaluationsException"> if the allowed number 
        /// of evaluations is exceeded.</exception>
        double solve(int maxEval, FUNC f, double min, double max, double startValue, AllowedSolution allowedSolution);
    }
}
