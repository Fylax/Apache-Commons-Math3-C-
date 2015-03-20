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
using Math3.util;
using System;

namespace Math3.analysis.solvers
{
    /// <summary>
    /// Provide a default implementation for several functions useful to generic
    /// solvers.
    /// </summary>
    /// <typeparam name="FUNC">Type of function to solve.</typeparam>
    public abstract class BaseAbstractUnivariateSolver<FUNC> : BaseUnivariateSolver<FUNC> where FUNC : UnivariateFunction
    {
        /// <summary>
        /// Default relative accuracy.
        /// </summary>
        private const double DEFAULT_RELATIVE_ACCURACY = 1e-14;

        /// <summary>
        /// Default function value accuracy.
        /// </summary>
        private const double DEFAULT_FUNCTION_VALUE_ACCURACY = 1e-15;

        /// <summary>
        /// Function value accuracy.
        /// </summary>
        private readonly double functionValueAccuracy;

        /// <summary>
        /// Absolute accuracy.
        /// </summary>
        private readonly double absoluteAccuracy;

        /// <summary>
        /// Relative accuracy.
        /// </summary>
        private readonly double relativeAccuracy;

        /// <summary>
        /// Evaluations counter.
        /// </summary>
        private readonly Incrementor evaluations = new Incrementor();

        /// <summary>
        /// Lower end of search interval.
        /// </summary>
        private double searchMin;

        /// <summary>
        /// Higher end of search interval.
        /// </summary>
        private double searchMax;

        /// <summary>
        /// Initial guess.
        /// </summary>
        private double searchStart;

        /// <summary>
        /// Function to solve.
        /// </summary>
        private FUNC function;

        /// <summary>
        /// Construct a solver with given absolute accuracy.
        /// </summary>
        /// <param name="absoluteAccuracy">Maximum absolute error.</param>
        protected BaseAbstractUnivariateSolver(double absoluteAccuracy) : this(DEFAULT_RELATIVE_ACCURACY, absoluteAccuracy, DEFAULT_FUNCTION_VALUE_ACCURACY) { }

        /// <summary>
        /// Construct a solver with given accuracies.
        /// </summary>
        /// <param name="relativeAccuracy">Maximum relative error.</param>
        /// <param name="absoluteAccuracy">Maximum absolute error.</param>
        protected BaseAbstractUnivariateSolver(double relativeAccuracy, double absoluteAccuracy) : this(relativeAccuracy, absoluteAccuracy, DEFAULT_FUNCTION_VALUE_ACCURACY) { }

        /// <summary>
        /// Construct a solver with given accuracies.
        /// </summary>
        /// <param name="relativeAccuracy">Maximum relative error.</param>
        /// <param name="absoluteAccuracy">Maximum absolute error.</param>
        /// <param name="functionValueAccuracy">Maximum function value error.</param>
        protected BaseAbstractUnivariateSolver(double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy)
        {
            this.absoluteAccuracy = absoluteAccuracy;
            this.relativeAccuracy = relativeAccuracy;
            this.functionValueAccuracy = functionValueAccuracy;
        }

        /// <inheritdoc/>
        public int getMaxEvaluations()
        {
            return evaluations.getMaximalCount();
        }

        /// <inheritdoc/>
        public int getEvaluations()
        {
            return evaluations.getCount();
        }

        /// <summary></summary>
        /// <returns>the lower end of the search interval.</returns>
        public double getMin()
        {
            return searchMin;
        }

        /// <summary></summary>
        /// <returns>the higher end of the search interval.</returns>
        public double getMax()
        {
            return searchMax;
        }

        /// <summary></summary>
        /// <returns>the initial guess.</returns>
        public double getStartValue()
        {
            return searchStart;
        }

        /// <inheritdoc/>
        public double getAbsoluteAccuracy()
        {
            return absoluteAccuracy;
        }

        /// <inheritdoc/>
        public double getRelativeAccuracy()
        {
            return relativeAccuracy;
        }

        /// <inheritdoc/>
        public double getFunctionValueAccuracy()
        {
            return functionValueAccuracy;
        }

        /// <summary>
        /// Compute the objective function value.
        /// </summary>
        /// <param name="point">Point at which the objective function must be evaluated.</param>
        /// <returns>the objective function value at specified point.</returns>
        /// <exception cref="TooManyEvaluationsException"> if the maximal number of evaluations
        /// is exceeded.</exception>
        protected double computeObjectiveValue(double point)
        {
            incrementEvaluationCount();
            return function.value(point);
        }

        /// <summary>
        /// Prepare for computation.
        /// Subclasses must call this method if they override any of the
        /// <c>solve</c> methods.
        /// </summary>
        /// <param name="maxEval"></param>
        /// <param name="f">Function to solve.</param>
        /// <param name="min">Lower bound for the interval</param>
        /// <param name="max">Upper bound for the interval.</param>
        /// <param name="startValue">Start value to use.</param>
        /// <exception cref="NullArgumentException"> if f is null</exception>
        protected void setup(int maxEval, FUNC f, double min, double max, double startValue)
        {
            // Checks.
            MathUtils.checkNotNull(f);

            // Reset.
            searchMin = min;
            searchMax = max;
            searchStart = startValue;
            function = f;
            evaluations.setMaximalCount(maxEval);
            evaluations.resetCount();
        }

        /// <inheritdoc/>
        public double solve(int maxEval, FUNC f, double min, double max, double startValue)
        {
            // Initialization.
            setup(maxEval, f, min, max, startValue);

            // Perform computation.
            return doSolve();
        }

        /// <inheritdoc/>
        public double solve(int maxEval, FUNC f, double min, double max)
        {
            return solve(maxEval, f, min, max, min + 0.5 * (max - min));
        }

        /// <inheritdoc/>
        public double solve(int maxEval, FUNC f, double startValue)
        {
            return solve(maxEval, f, Double.NaN, Double.NaN, startValue);
        }

        /// <summary>
        /// Method for implementing actual optimization algorithms in derived
        /// classes.
        /// </summary>
        /// <returns>the root.</returns>
        /// <exception cref="TooManyEvaluationsException"> if the maximal number of evaluations
        /// is exceeded.</exception>
        /// <exception cref="NoBracketingException"> if the initial search interval does
        /// not bracket a root and the solver requires it.</exception>
        protected abstract double doSolve();

        /// <summary>
        /// Check whether the function takes opposite signs at the endpoints.
        /// </summary>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <returns><c>true</c> if the function values have opposite signs at the
        /// given points.</returns>
        protected Boolean isBracketing(double lower, double upper)
        {
            return UnivariateSolverUtils.isBracketing(function, lower, upper);
        }

        /// <summary>
        /// Check whether the arguments form a (strictly) increasing sequence.
        /// </summary>
        /// <param name="start">First number.</param>
        /// <param name="mid">Second number.</param>
        /// <param name="end">Third number.</param>
        /// <returns><c>true</c> if the arguments form an increasing sequence.</returns>
        protected Boolean isSequence(double start, double mid, double end)
        {
            return UnivariateSolverUtils.isSequence(start, mid, end);
        }

        /// <summary>
        /// Check that the endpoints specify an interval. 
        /// </summary>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c></exception>
        protected void verifyInterval(double lower, double upper)
        {
            UnivariateSolverUtils.verifyInterval(lower, upper);
        }

        /// <summary>
        /// Check that <c>lower < initial < upper</c>.
        /// </summary>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="initial">Initial value.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= initial</c> or
        /// <c>initial >= upper</c>.</exception>
        protected void verifySequence(double lower, double initial, double upper)
        {
            UnivariateSolverUtils.verifySequence(lower, initial, upper);
        }

        /// <summary>
        /// Check that the endpoints specify an interval and the function takes
        /// opposite signs at the endpoints.
        /// </summary>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <exception cref="NullArgumentException"> if the function has not been set.</exception>
        /// <exception cref="NoBracketingException"> if the function has the same sign at
        /// the endpoints.</exception>
        protected void verifyBracketing(double lower, double upper)
        {
            UnivariateSolverUtils.verifyBracketing(function, lower, upper);
        }

        /// <summary>
        /// Increment the evaluation count by one.
        /// Method <see cref="computeObjectiveValue(double)"/> calls this method internally.
        /// It is provided for subclasses that do not exclusively use
        /// <c>computeObjectiveValue</c> to solve the function.
        /// See e.g. <see cref="AbstractUnivariateDifferentiableSolver"/>.
        /// </summary>
        /// <exception cref="TooManyEvaluationsException"> when the allowed number of function
        /// evaluations has been exhausted.</exception>
        protected void incrementEvaluationCount()
        {
            try
            {
                evaluations.incrementCount();
            }
            catch (MaxCountExceededException<Int32> e)
            {
                throw new TooManyEvaluationsException<Int32>(e.getMax());
            }
        }
    }
}
