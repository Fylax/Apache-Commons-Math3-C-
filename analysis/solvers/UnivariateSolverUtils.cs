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
using Math3.util;
using System;

namespace Math3.analysis.solvers
{
    /// <summary>
    /// Utility routines for <see cref="UnivariateSolver"/> objects.
    /// </summary>
    public class UnivariateSolverUtils
    {
        /// <summary>
        /// Class contains only static methods.
        /// </summary>
        private UnivariateSolverUtils() { }

        /// <summary>
        /// Convenience method to find a zero of a univariate real function.  A default
        /// solver is used.
        /// </summary>
        /// <param name="function">Function.</param>
        /// <param name="x0">Lower bound for the interval.</param>
        /// <param name="x1">Upper bound for the interval.</param>
        /// <returns>a value where the function is zero.</returns>
        /// <exception cref="NoBracketingException"> if the function has the same sign at the
        /// endpoints.</exception>
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static double solve(UnivariateFunction function, double x0, double x1)
        {
            if (function == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FUNCTION"));
            }
            UnivariateSolver solver = new BrentSolver();
            return solver.solve(Int32.MaxValue, function, x0, x1);
        }

        /// <summary>
        /// Convenience method to find a zero of a univariate real function.  A default
        /// solver is used.
        /// </summary>
        /// <param name="function">Function.</param>
        /// <param name="x0">Lower bound for the interval.</param>
        /// <param name="x1">Upper bound for the interval.</param>
        /// <param name="absoluteAccuracy">Accuracy to be used by the solver.</param>
        /// <returns>a value where the function is zero.</returns>
        /// <exception cref="NoBracketingException"> if the function has the same sign at the
        /// endpoints.</exception>
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static double solve(UnivariateFunction function, double x0, double x1, double absoluteAccuracy)
        {
            if (function == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FUNCTION"));
            }
            UnivariateSolver solver = new BrentSolver(absoluteAccuracy);
            return solver.solve(Int32.MaxValue, function, x0, x1);
        }

        /// <summary>
        /// Force a root found by a non-bracketing solver to lie on a specified side,
        /// as if the solver was a bracketing one.
        /// </summary>
        /// <param name="maxEval">maximal number of new evaluations of the function
        /// (evaluations already done for finding the root should have already been 
        /// subtracted from this number)</param>
        /// <param name="f">function to solve</param>
        /// <param name="bracketing">bracketing solver to use for shifting the root</param>
        /// <param name="baseRoot">original root found by a previous non-bracketing solver</param>
        /// <param name="min">minimal bound of the search interval</param>
        /// <param name="max">maximal bound of the search interval</param>
        /// <param name="allowedSolution">the kind of solutions that the root-finding
        /// algorithm may accept as solutions.</param>
        /// <returns>a root approximation, on the specified side of the exact root</returns>
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static double forceSide(int maxEval, UnivariateFunction f, BracketedUnivariateSolver<UnivariateFunction> bracketing, double baseRoot, double min, double max, AllowedSolution allowedSolution)
        {

            if (allowedSolution == AllowedSolution.ANY_SIDE)
            {
                // no further bracketing required
                return baseRoot;
            }

            // find a very small interval bracketing the root
            double step = FastMath.max(bracketing.getAbsoluteAccuracy(),
                                             FastMath.abs(baseRoot * bracketing.getRelativeAccuracy()));
            double xLo = FastMath.max(min, baseRoot - step);
            double fLo = f.value(xLo);
            double xHi = FastMath.min(max, baseRoot + step);
            double fHi = f.value(xHi);
            int remainingEval = maxEval - 2;
            while (remainingEval > 0)
            {

                if ((fLo >= 0 && fHi <= 0) || (fLo <= 0 && fHi >= 0))
                {
                    // compute the root on the selected side
                    return bracketing.solve(remainingEval, f, xLo, xHi, baseRoot, allowedSolution);
                }

                // try increasing the interval
                Boolean changeLo = false;
                Boolean changeHi = false;
                if (fLo < fHi)
                {
                    // increasing function
                    if (fLo >= 0)
                    {
                        changeLo = true;
                    }
                    else
                    {
                        changeHi = true;
                    }
                }
                else if (fLo > fHi)
                {
                    // decreasing function
                    if (fLo <= 0)
                    {
                        changeLo = true;
                    }
                    else
                    {
                        changeHi = true;
                    }
                }
                else
                {
                    // unknown variation
                    changeLo = true;
                    changeHi = true;
                }

                // update the lower bound
                if (changeLo)
                {
                    xLo = FastMath.max(min, xLo - step);
                    fLo = f.value(xLo);
                    remainingEval--;
                }

                // update the higher bound
                if (changeHi)
                {
                    xHi = FastMath.min(max, xHi + step);
                    fHi = f.value(xHi);
                    remainingEval--;
                }

            }

            throw new NoBracketingException(new LocalizedFormats("FAILED_BRACKETING"),
                                            xLo, xHi, fLo, fHi,
                                            maxEval - remainingEval, maxEval, baseRoot,
                                            min, max);

        }

        /// <summary>
        /// This method simply calls 
        /// <see cref="bracket(UnivariateFunction, double, double, double, double, double, int)"/> 
        /// with <c>q</c> and <c>r</c> set to 1.0 and <c>maximumIterations</c> set to
        /// <c>Int32.MaxValue</c>.
        /// <para>
        /// Note: this method can take
        /// <c>Int32.MaxValue</c> iterations to throw a
        /// <c>ConvergenceException.</c>  Unless you are confident that there
        /// is a root between <c>lowerBound</c> and <c>upperBound</c>
        /// near <c>initial,</c> it is better to use
        /// <see cref="bracket(UnivariateFunction, double, double, double, double, double, int)"/>
        /// , explicitly specifying the maximum number of iterations.</para>
        /// </summary>
        /// <param name="function">Function.</param>
        /// <param name="initial">Initial midpoint of interval being expanded to
        /// bracket a root.</param>
        /// <param name="lowerBound">Lower bound (a is never lower than this value)</param>
        /// <param name="upperBound">Upper bound (b never is greater than this value).</param>
        /// <returns>a two-element array holding a and b.</returns>
        /// <exception cref="NoBracketingException"> if a root cannot be bracketted.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if <c>maximumIterations <= 0</c></exception>.
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static double[] bracket(UnivariateFunction function, double initial, double lowerBound, double upperBound)
        {
            return bracket(function, initial, lowerBound, upperBound, 1.0, 1.0, Int32.MaxValue);
        }

        /// <summary>
        /// This method simply calls 
        /// <see cref="bracket(UnivariateFunction, double, double, double, double, double, int)"/>
        /// with <c>q</c> and <c>r</c> set to 1.0.
        /// </summary>
        /// <param name="function">Function.</param>
        /// <param name="initial">Initial midpoint of interval being expanded to
        /// bracket a root.</param>
        /// <param name="lowerBound">Lower bound (a is never lower than this value).</param>
        /// <param name="upperBound">Upper bound (b never is greater than this value).</param>
        /// <param name="maximumIterations">Maximum number of iterations to perform</param>
        /// <returns>a two element array holding a and b.</returns>
        /// <exception cref="NoBracketingException"> if the algorithm fails to find a and b
        /// satisfying the desired conditions.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if <c>maximumIterations <= 0</c></exception>.
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static double[] bracket(UnivariateFunction function, double initial, double lowerBound, double upperBound, int maximumIterations)
        {
            return bracket(function, initial, lowerBound, upperBound, 1.0, 1.0, maximumIterations);
        }

        /// <summary>
        /// This method attempts to find two values a and b satisfying 
        /// <list type="bullet">
        /// <item> <c>lowerBound <= a < initial < b <= upperBound</c> </item>
        /// <item> <c>f(a) * f(b) <= 0</c> </item>
        /// </list>
        /// If <c>f</c> is continuous on <c>[a,b]</c>, this means that <c>a</c>
        /// and <c>b</c> bracket a root of <c>f</c>.
        /// <para>
        /// The algorithm checks the sign of \( f(l_k) \) and \( f(u_k) \) for increasing
        /// values of k, where \( l_k = max(lower, initial - \delta_k) \),
        /// \( u_k = min(upper, initial + \delta_k) \), using recurrence
        /// \( \delta_{k+1} = r \delta_k + q, \delta_0 = 0\) and starting search with \( k=1 \).
        /// The algorithm stops when one of the following happens: 
        /// <list type="bullet">
        /// <item> at least one positive and one negative value have been found --  success!</item>
        /// <item> both endpoints have reached their respective limites -- NoBracketingException </item>
        /// <item> <c>maximumIterations</c> iterations elapse -- NoBracketingException </item>
        /// </list></para>
        /// <para>
        /// If different signs are found at first iteration (<c>k=1</c>), then the returned
        /// interval will be \( [a, b] = [l_1, u_1] \). If different signs are found at a later
        /// iteration (<c>k>1</c>, then the returned interval will be either
        /// \( [a, b] = [l_{k+1}, l_{k}] \) or \( [a, b] = [u_{k}, u_{k+1}] \). A root 
        /// solver called with these parameters will therefore start with the smallest 
        /// bracketing interval known at this step.
        /// </para>
        /// <para>
        /// Interval expansion rate is tuned by changing the recurrence parameters <c>r</c> and
        /// <c>q</c>. When the multiplicative factor <c>r</c> is set to 1, the sequence is a
        /// simple arithmetic sequence with linear increase. When the multiplicative factor 
        /// <c>r</c> is larger than 1, the sequence has an asymtotically exponential rate.
        /// Note than the additive parameter <c>q</c> should never be set to zero, otherwise
        /// the interval would degenerate to the single initial point for all values of <c>k</c>.
        /// </para>
        /// <para> 
        /// As a rule of thumb, when the location of the root is expected to be approximately 
        /// known within some error margin, <c>r</c> should be set to 1 and <c>q</c> should be 
        /// set to the order of magnitude of the error margin. When the location of the root 
        /// is really a wild guess, then <c>r</c> should be set to a value larger than 1 
        /// (typically 2 to double the interval length at each iteration) and <c>q</c> should
        /// be set according to half the initial search interval length.
        /// </para>
        /// <para>
        /// As an example, if we consider the trivial function <c>f(x) = 1 - x</c> and use
        /// <c>initial = 4</c>, <c>r = 1</c>, <c>q = 2</c>, the algorithm will compute
        /// <c>f(4-2) = f(2) = -1</c> and <c>f(4+2) = f(6) = -5</c> for <c>k = 1</c>, then
        /// <c>f(4-4) = f(0) = +1</c> and <c>f(4+4) = f(8) = -7</c> for <c>k = 2</c>. Then it will
        /// return the interval <c>[0, 2]</c> as the smallest one known to be bracketing the root.
        /// As shown by this example, the initial value (here <c>4</c>) may lie outside of the
        /// returned bracketing interval.
        /// </para>
        /// </summary>
        /// <param name="function">function to check</param>
        /// <param name="initial">Initial midpoint of interval being expanded to bracket a 
        /// root.</param>
        /// <param name="lowerBound">Lower bound (a is never lower than this value).</param>
        /// <param name="upperBound">Upper bound (b never is greater than this value).</param>
        /// <param name="q">additive offset used to compute bounds sequence (must be strictly 
        /// positive)</param>
        /// <param name="r">multiplicative factor used to compute bounds sequence</param>
        /// <param name="maximumIterations">Maximum number of iterations to perform</param>
        /// <returns>a two element array holding the bracketing values.</returns>
        /// <exception cref="NoBracketingException"> if function cannot be bracketed in the
        /// search interval</exception>
        public static double[] bracket(UnivariateFunction function, double initial, double lowerBound, double upperBound, double q, double r, int maximumIterations)
        {

            if (function == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FUNCTION"));
            }
            if (q <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(q);
            }
            if (maximumIterations <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("INVALID_MAX_ITERATIONS"), maximumIterations);
            }
            verifySequence(lowerBound, initial, upperBound);

            // initialize the recurrence
            double a = initial;
            double b = initial;
            double fa = Double.NaN;
            double fb = Double.NaN;
            double delta = 0;

            for (int numIterations = 0;
                 (numIterations < maximumIterations) && (a > lowerBound || b > upperBound);
                 ++numIterations)
            {

                double previousA = a;
                double previousFa = fa;
                double previousB = b;
                double previousFb = fb;

                delta = r * delta + q;
                a = FastMath.max(initial - delta, lowerBound);
                b = FastMath.min(initial + delta, upperBound);
                fa = function.value(a);
                fb = function.value(b);

                if (numIterations == 0)
                {
                    // at first iteration, we don't have a previous interval
                    // we simply compare both sides of the initial interval
                    if (fa * fb <= 0)
                    {
                        // the first interval already brackets a root
                        return new double[] { a, b };
                    }
                }
                else
                {
                    // we have a previous interval with constant sign and expand it,
                    // we expect sign changes to occur at boundaries
                    if (fa * previousFa <= 0)
                    {
                        // sign change detected at near lower bound
                        return new double[] { a, previousA };
                    }
                    else if (fb * previousFb <= 0)
                    {
                        // sign change detected at near upper bound
                        return new double[] { previousB, b };
                    }
                }

            }

            // no bracketing found
            throw new NoBracketingException(a, b, fa, fb);

        }

        /// <summary>
        /// Compute the midpoint of two values.
        /// </summary>
        /// <param name="a">first value.</param>
        /// <param name="b">second value.</param>
        /// <returns>the midpoint.</returns>
        public static double midpoint(double a, double b)
        {
            return (a + b) * 0.5;
        }

        /// <summary>
        /// Check whether the interval bounds bracket a root. That is, if the
        /// values at the endpoints are not equal to zero, then the function takes
        /// opposite signs at the endpoints.
        /// </summary>
        /// <param name="function">Function.</param>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <returns><c>true</c> if the function values have opposite signs at the
        /// given points.</returns>
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static Boolean isBracketing(UnivariateFunction function, double lower, double upper)
        {
            if (function == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FUNCTION"));
            }
            double fLo = function.value(lower);
            double fHi = function.value(upper);
            return (fLo >= 0 && fHi <= 0) || (fLo <= 0 && fHi >= 0);
        }

        /// <summary>
        /// Check whether the arguments form a (strictly) increasing sequence.
        /// </summary>
        /// <param name="start">First number.</param>
        /// <param name="mid">Second number.</param>
        /// <param name="end">Third number.</param>
        /// <returns><c>true</c> if the arguments form an increasing sequence.</returns>
        public static Boolean isSequence(double start, double mid, double end)
        {
            return (start < mid) && (mid < end);
        }

        /// <summary>
        /// Check that the endpoints specify an interval.
        /// </summary>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c>.</exception>
        public static void verifyInterval(double lower, double upper)
        {
            if (lower >= upper)
            {
                throw new NumberIsTooLargeException<Double, Double>(new LocalizedFormats("ENDPOINTS_NOT_AN_INTERVAL"), lower, upper, false);
            }
        }

        /// <summary>
        /// Check that <c>lower < initial < upper</c>.
        /// </summary>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="initial">Initial value.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= initial</c> or
        /// <c>initial >= upper</c>.</exception>
        public static void verifySequence(double lower, double initial, double upper)
        {
            verifyInterval(lower, initial);
            verifyInterval(initial, upper);
        }

        /// <summary>
        /// Check that the endpoints specify an interval and the end points
        /// bracket a root.
        /// </summary>
        /// <param name="function">Function.</param>
        /// <param name="lower">Lower endpoint.</param>
        /// <param name="upper">Upper endpoint.</param>
        /// <exception cref="NoBracketingException"> if the function has the same sign at the
        /// endpoints.</exception>
        /// <exception cref="NullArgumentException"> if <c>function</c> is <c>null</c>.</exception>
        public static void verifyBracketing(UnivariateFunction function, double lower, double upper)
        {
            if (function == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FUNCTION"));
            }
            verifyInterval(lower, upper);
            if (!isBracketing(function, lower, upper))
            {
                throw new NoBracketingException(lower, upper, function.value(lower), function.value(upper));
            }
        }
    }
}
