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

namespace Math3.util
{
    /// <summary>
    /// Provides a generic means to evaluate continued fractions.  Subclasses simply
    /// provided the a and b coefficients to evaluate the continued fraction.
    /// <para>
    /// References:
    /// <list type="bullet">
    /// <item><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
    /// Continued Fraction</a></item>
    /// </list>
    /// </para>
    /// </summary>
    public abstract class ContinuedFraction
    {
        /// <summary>
        /// Maximum allowed numerical error.
        /// </summary>
        private const double DEFAULT_EPSILON = 10e-9;

        /// <summary>
        /// Default constructor.
        /// </summary>
        protected ContinuedFraction() : base() { }

        /// <summary>
        /// Access the n-th a coefficient of the continued fraction.  Since a can be
        /// a function of the evaluation point, x, that is passed in as well.
        /// </summary>
        /// <param name="n">the coefficient index to retrieve.</param>
        /// <param name="x">the evaluation point.</param>
        /// <returns>the n-th a coefficient.</returns>
        protected abstract double getA(int n, double x);

        /// <summary>
        /// Access the n-th b coefficient of the continued fraction.  Since b can be
        /// a function of the evaluation point, x, that is passed in as well.
        /// </summary>
        /// <param name="n">the coefficient index to retrieve.</param>
        /// <param name="x">the evaluation point.</param>
        /// <returns>the n-th b coefficient.</returns>
        protected abstract double getB(int n, double x);

        /// <summary>
        /// Evaluates the continued fraction at the value x.
        /// </summary>
        /// <param name="x">the evaluation point.</param>
        /// <returns>the value of the continued fraction evaluated at x.</returns>
        /// <exception cref="ConvergenceException"> if the algorithm fails to converge.</exception>
        public double evaluate(double x)
        {
            return evaluate(x, DEFAULT_EPSILON, Int32.MaxValue);
        }

        /// <summary>
        /// Evaluates the continued fraction at the value x.
        /// </summary>
        /// <param name="x">the evaluation point.</param>
        /// <param name="epsilon">maximum error allowed.</param>
        /// <returns>the value of the continued fraction evaluated at x.</returns>
        /// <exception cref="ConvergenceException"> if the algorithm fails to converge.</exception>
        public double evaluate(double x, double epsilon)
        {
            return evaluate(x, epsilon, Int32.MaxValue);
        }

        /// <summary>
        /// Evaluates the continued fraction at the value x.
        /// </summary>
        /// <param name="x">the evaluation point.</param>
        /// <param name="maxIterations">maximum number of convergents</param>
        /// <returns>the value of the continued fraction evaluated at x.</returns>
        /// <exception cref="ConvergenceException"> if the algorithm fails to converge.</exception>
        /// <exception cref="MaxCountExceededException"> if maximal number of iterations
        /// is reached</exception>
        public double evaluate(double x, int maxIterations)
        {
            return evaluate(x, DEFAULT_EPSILON, maxIterations);
        }

        /// <summary>
        /// Evaluates the continued fraction at the value x.
        /// <para>
        /// The implementation of this method is based on the modified Lentz algorithm as described
        /// on page 18 ff. in:
        /// <list type="bullet">
        /// <item>
        /// I. J. Thompson,  A. R. Barnett. "Coulomb and Bessel Functions of Complex Arguments and Order."
        /// <a target="_blank" href="http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf">
        /// http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf </a>
        /// </item>
        /// </list>
        /// Note: the implementation uses the terms a_i and b_i as defined in
        /// <a href="http://mathworld.wolfram.com/ContinuedFraction.html">Continued Fraction @ MathWorld</a>.
        /// </para>
        /// </summary>
        /// <param name="x">the evaluation point.</param>
        /// <param name="epsilon">maximum error allowed.</param>
        /// <param name="maxIterations">maximum number of convergents</param>
        /// <returns>the value of the continued fraction evaluated at x.</returns>
        /// <exception cref="ConvergenceException"> if the algorithm fails to converge.</exception>
        /// <exception cref="MaxCountExceededException"> if maximal number of iterations
        /// is reached</exception>
        public double evaluate(double x, double epsilon, int maxIterations)
        {
            double small = 1e-50;
            double hPrev = getA(0, x);

            // use the value of small as epsilon criteria for zero checks
            if (Precision.equals(hPrev, 0.0, small))
            {
                hPrev = small;
            }

            int n = 1;
            double dPrev = 0.0;
            double cPrev = hPrev;
            double hN = hPrev;

            while (n < maxIterations)
            {
                double a = getA(n, x);
                double b = getB(n, x);

                double dN = a + b * dPrev;
                if (Precision.equals(dN, 0.0, small))
                {
                    dN = small;
                }
                double cN = a + b / cPrev;
                if (Precision.equals(cN, 0.0, small))
                {
                    cN = small;
                }

                dN = 1 / dN;
                double deltaN = cN * dN;
                hN = hPrev * deltaN;

                if (Double.IsInfinity(hN))
                {
                    throw new ConvergenceException(new LocalizedFormats("CONTINUED_FRACTION_INFINITY_DIVERGENCE"), x);
                }
                if (Double.IsNaN(hN))
                {
                    throw new ConvergenceException(new LocalizedFormats("CONTINUED_FRACTION_NAN_DIVERGENCE"), x);
                }

                if (FastMath.abs(deltaN - 1.0) < epsilon)
                {
                    break;
                }

                dPrev = dN;
                cPrev = cN;
                hPrev = hN;
                n++;
            }

            if (n >= maxIterations)
            {
                throw new MaxCountExceededException<Int32>(new LocalizedFormats("NON_CONVERGENT_CONTINUED_FRACTION"), maxIterations, x);
            }
            return hN;
        }
    }
}
