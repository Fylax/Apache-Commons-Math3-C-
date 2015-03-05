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
using Math3.random;
using Math3.util;
using System;

namespace Math3.distribution
{
    /// <summary>
    /// Base class for integer-valued discrete distributions.  Default
    /// implementations are provided for some of the methods that do not vary
    /// from distribution to distribution.
    /// </summary>
    public abstract class AbstractIntegerDistribution : IntegerDistribution
    {
        /// <summary>
        /// RandomData instance used to generate samples from the distribution.
        /// </summary>
        #pragma warning disable 0612
        [Obsolete("Please use the random instance variable instead")]
        protected readonly RandomDataImpl randomData = new RandomDataImpl();
        #pragma warning restore 0612

        /// <summary>
        /// RNG instance used to generate samples from the distribution.
        /// </summary>
        protected readonly RandomGenerator random;


        [Obsolete("Please use AbtractIntegerDistribution(RandomGeerator) intead")]
        protected AbstractIntegerDistribution()
        {
            // Legacy users are only allowed to access the deprecated "randomData".
            // New users are forbidden to use this constructor.
            random = null;
        }

        /// <param name="rng">Random number generator.</param>
        protected AbstractIntegerDistribution(RandomGenerator rng)
        {
            random = rng;
        }

        /// <inheritdoc/>
        /// <remarks>The default implementation uses the identity
        /// <para><c>P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)</c></para></remarks>
        public double cumulativeProbability(int x0, int x1)
        {
            if (x1 < x0)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT"), x0, x1, true);
            }
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        /// <inheritdoc/>
        /// <remarks>The default implementation returns
        /// <list type="bullet">
        /// <item><see cref="getSupportLowerBound()"/> for <c>p = 0</c>,</item>
        /// <item><see cref="getSupportUpperBound()"/> for <c>p = 1</c>, and</item>
        /// <item><see cref="solveInverseCumulativeProbability(double, int, int)"/> for
        /// <c>0 < p < 1</c>.</item>
        /// </list>
        public int inverseCumulativeProbability(double p)
        {
            if (p < 0.0 || p > 1.0)
            {
                throw new OutOfRangeException<Double>(p, 0, 1);
            }

            int lower = getSupportLowerBound();
            if (p == 0.0)
            {
                return lower;
            }
            if (lower == Int32.MinValue)
            {
                if (checkedCumulativeProbability(lower) >= p)
                {
                    return lower;
                }
            }
            else
            {
                lower -= 1; // this ensures cumulativeProbability(lower) < p, which
                // is important for the solving step
            }

            int upper = getSupportUpperBound();
            if (p == 1.0)
            {
                return upper;
            }

            // use the one-sided Chebyshev inequality to narrow the bracket
            // cf. AbstractRealDistribution.inverseCumulativeProbability(double)
            double mu = getNumericalMean();
            double sigma = FastMath.sqrt(getNumericalVariance());
            Boolean chebyshevApplies = !(Double.IsInfinity(mu) || Double.IsNaN(mu) ||
                    Double.IsInfinity(sigma) || Double.IsNaN(sigma) || sigma == 0.0);
            if (chebyshevApplies)
            {
                double k = FastMath.sqrt((1.0 - p) / p);
                double tmp = mu - k * sigma;
                if (tmp > lower)
                {
                    lower = ((int)FastMath.ceil(tmp)) - 1;
                }
                k = 1.0 / k;
                tmp = mu + k * sigma;
                if (tmp < upper)
                {
                    upper = ((int)FastMath.ceil(tmp)) - 1;
                }
            }

            return solveInverseCumulativeProbability(p, lower, upper);
        }

        /// <summary>
        /// This is a utility function used by
        /// <see cref="inverseCumulativeProbability(double)}"/>. It assumes <c>0 < p < 1</c> and
        /// that the inverse cumulative probability lies in the bracket <c>(lower, upper]</c>.
        /// The implementation does simple bisection to find the smallest <c>p</c>-quantile 
        /// <c>inf{x in Z | P(X<=x) >= p}</c>.
        /// </summary>
        /// <param name="p">the cumulative probability</param>
        /// <param name="lower">a value satisfying <c>cumulativeProbability(lower) < p</c></param>
        /// <param name="upper">a value satisfying <c>p <= cumulativeProbability(upper)</c></param>
        /// <returns>the smallest <c>p</c>-quantile of this distribution</returns>
        protected int solveInverseCumulativeProbability(double p, int lower, int upper)
        {
            while (lower + 1 < upper)
            {
                int xm = (lower + upper) / 2;
                if (xm < lower || xm > upper)
                {
                    /*
                     * Overflow.
                     * There will never be an overflow in both calculation methods
                     * for xm at the same time
                     */
                    xm = lower + (upper - lower) / 2;
                }

                double pm = checkedCumulativeProbability(xm);
                if (pm >= p)
                {
                    upper = xm;
                }
                else
                {
                    lower = xm;
                }
            }
            return upper;
        }

        #pragma warning disable 0618
        /// <inheritdoc/>
        public void reseedRandomGenerator(long seed)
        {
            random.setSeed(seed);
            randomData.reSeed(seed);
        }
        #pragma warning restore 0618

        /// <inheritdoc/>
        /// <remarks>The default implementation uses the
        /// <a href="http://en.wikipedia.org/wiki/Inverse_transform_sampling">
        /// inversion method</a>.</remarks>
        public int sample()
        {
            return inverseCumulativeProbability(random.nextDouble());
        }

        /// <inheritdoc/>
        /// <remarks>The default implementation generates the sample by calling
        /// <see cref="sample()"/> in a loop.</remarks>
        public int[] sample(int sampleSize)
        {
            if (sampleSize <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_SAMPLES"), sampleSize);
            }
            int[] outp = new int[sampleSize];
            for (int i = 0; i < sampleSize; i++)
            {
                outp[i] = sample();
            }
            return outp;
        }

        /// <summary>
        /// Computes the cumulative probability function and checks for <c>NaN</c>
        /// values returned. Throws <c>MathInternalError</c> if the value is
        /// <c>NaN</c>. Rethrows any exception encountered evaluating the cumulative
        /// probability function. Throws <c>MathInternalError</c> if the cumulative
        /// probability function returns <c>NaN</c>.
        /// </summary>
        /// <param name="argument">input value</param>
        /// <returns>the cumulative probability</returns>
        /// <exception cref="MathInternalError"> if the cumulative probability is <c>NaN</c></exception>
        private double checkedCumulativeProbability(int argument)
        {
            double result = Double.NaN;
            result = cumulativeProbability(argument);
            if (Double.IsNaN(result))
            {
                throw new MathInternalError(new LocalizedFormats("DISCRETE_CUMULATIVE_PROBABILITY_RETURNED_NAN"), argument);
            }
            return result;
        }

        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according to
        /// this distribution, this method returns <c>log(P(X = x))</c>, where
        /// <c>log</c> is the natural logarithm. In other words, this method
        /// represents the logarithm of the probability mass function (PMF) for the
        /// distribution. Note that due to the floating point precision and
        /// under/overflow issues, this method will for some distributions be more
        /// precise and faster than computing the logarithm of <see cref="probability(int)"/>.
        /// <para>
        /// The default implementation simply computes the logarithm of <see cref="probability(x)"/>.</para> 
        /// </summary>
        /// <param name="x">the point at which the PMF is evaluated</param>
        /// <returns>the logarithm of the value of the probability mass function at </c>x</c></returns>
        public double logProbability(int x)
        {
            return FastMath.log(probability(x));
        }

        /// <inheritdoc/>
        public abstract double probability(int x);

        /// <inheritdoc/>
        public abstract double cumulativeProbability(int x);

        /// <inheritdoc/>
        public abstract double getNumericalMean();

        /// <inheritdoc/>
        public abstract double getNumericalVariance();

        /// <inheritdoc/>
        public abstract int getSupportLowerBound();

        /// <inheritdoc/>
        public abstract int getSupportUpperBound();

        /// <inheritdoc/>
        public abstract Boolean isSupportConnected();
    }
}
