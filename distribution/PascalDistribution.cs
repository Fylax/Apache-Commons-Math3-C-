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
using Math3.special;
using Math3.util;
using System;

namespace Math3.distribution
{
    /// <summary>
    /// <para>
    /// Implementation of the Pascal distribution. The Pascal distribution is a
    /// special case of the Negative Binomial distribution where the number of
    /// successes parameter is an integer.
    /// </para>
    /// <para>
    /// There are various ways to express the probability mass and distribution
    /// functions for the Pascal distribution. The present implementation represents
    /// the distribution of the number of failures before <c>r</r> successes occur.
    /// This is the convention adopted in e.g.
    /// <a href="http://mathworld.wolfram.com/NegativeBinomialDistribution.html">MathWorld</a>,
    /// but not in
    /// <a href="http://en.wikipedia.org/wiki/Negative_binomial_distribution">Wikipedia</a>.
    /// </para><para>
    /// For a random variable <c>X</c> whose values are distributed according to this
    /// distribution, the probability mass function is given by<para/>
    /// <c>P(X = k) = C(k + r - 1, r - 1) * p^r * (1 - p)^k</c>,<para/>
    /// where <c>r</c> is the number of successes, <c>p</c> is the probability of
    /// success, and <c>X</c> is the total number of failures. <c>C(n, k)</c> is
    /// the binomial coefficient (<c>n</c> choose <c>k</c>). The mean and variance
    /// of <c>X</c> are<para/>
    /// <c>E(X) = (1 - p) * r / p, var(X) = (1 - p) * r / p^2.</c><para/>
    /// Finally, the cumulative distribution function is given by<para/>
    /// <c>P(X <= k) = I(p, r, k + 1)</c>,
    /// where I is the regularized incomplete Beta function.
    /// </para>
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Negative_binomial_distribution">
    /// Negative binomial distribution (Wikipedia)</a></para>
    /// See <a href="http://mathworld.wolfram.com/NegativeBinomialDistribution.html">
    /// Negative binomial distribution (MathWorld)</a></remarks>
    public class PascalDistribution : AbstractIntegerDistribution
    {
        /// <summary>
        /// The number of successes.
        /// </summary>
        private readonly int numberOfSuccesses;

        /// <summary>
        /// The probability of success.
        /// </summary>
        private readonly double probabilityOfSuccess;

        /// <summary>
        /// The value of <c>log(p)</c>, where <c>p</c> is the probability of success,
        /// stored for faster computation.
        /// </summary>
        private readonly double logProbabilityOfSuccess;

        /// <summary>
        /// The value of <c>log(1-p)</c>, where <c>p</c> is the probability of success,
        /// stored for faster computation.
        /// </summary>
        private readonly double log1mProbabilityOfSuccess;

        /// <summary>
        /// Create a Pascal distribution with the given number of successes and
        /// probability of success.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para> 
        /// </summary>
        /// <param name="r">Number of successes.</param>
        /// <param name="p">Probability of success.</param>
        /// <exception cref="NotStrictlyPositiveException"> if the number of successes is not positive</exception>
        /// <exception cref="OutOfRangeException">if the probability of success is not in the
        /// range <c>[0, 1]</c>.</exception>
        public PascalDistribution(int r, double p) : this(new Well19937c(), r, p) { }

        /// <summary>
        /// Create a Pascal distribution with the given number of successes and
        /// probability of success.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="r">Number of successes.</param>
        /// <param name="p">Probability of success.</param>
        /// <exception cref="NotStrictlyPositiveException"> if the number of successes is not positive</exception>
        /// <exception cref="OutOfRangeException">if the probability of success is not in the
        /// range <c>[0, 1]</c>.</exception>
        public PascalDistribution(RandomGenerator rng, int r, double p)
            : base(rng)
        {
            if (r <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_SUCCESSES"), r);
            }
            if (p < 0 || p > 1)
            {
                throw new OutOfRangeException<Double>(p, 0, 1);
            }

            numberOfSuccesses = r;
            probabilityOfSuccess = p;
            logProbabilityOfSuccess = FastMath.log(p);
            log1mProbabilityOfSuccess = FastMath.log1p(-p);
        }

        /// <summary>
        /// Access the number of successes for this distribution.
        /// </summary>
        /// <returns>the number of successes.</returns>
        public int getNumberOfSuccesses()
        {
            return numberOfSuccesses;
        }

        /// <summary>
        /// Access the probability of success for this distribution.
        /// </summary>
        /// <returns>the probability of success.</returns>
        public double getProbabilityOfSuccess()
        {
            return probabilityOfSuccess;
        }

        /// <inheritdoc/>
        public override double probability(int x)
        {
            double ret;
            if (x < 0)
            {
                ret = 0.0;
            }
            else
            {
                ret = CombinatoricsUtils.binomialCoefficientDouble(x +
                      numberOfSuccesses - 1, numberOfSuccesses - 1) *
                      FastMath.pow(probabilityOfSuccess, numberOfSuccesses) *
                      FastMath.pow(1.0 - probabilityOfSuccess, x);
            }
            return ret;
        }

        /// <inheritdoc/>
        public new double logProbability(int x)
        {
            double ret;
            if (x < 0)
            {
                ret = Double.NegativeInfinity;
            }
            else
            {
                ret = CombinatoricsUtils.binomialCoefficientLog(x +
                      numberOfSuccesses - 1, numberOfSuccesses - 1) +
                      logProbabilityOfSuccess * numberOfSuccesses +
                      log1mProbabilityOfSuccess * x;
            }
            return ret;
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(int x)
        {
            double ret;
            if (x < 0)
            {
                ret = 0.0;
            }
            else
            {
                ret = Beta.regularizedBeta(probabilityOfSuccess,
                        numberOfSuccesses, x + 1.0);
            }
            return ret;
        }

        /// <inheritdoc/>
        /// <remarks>For number of successes <c>r</c> and probability of success <c>p</c>,
        /// the mean is <c>r * (1 - p) / p</c>.</remarks>
        public override double getNumericalMean()
        {
            double p = getProbabilityOfSuccess();
            double r = getNumberOfSuccesses();
            return (r * (1 - p)) / p;
        }

        /// <inheritdoc/>
        /// <remarks>For number of successes <c>r</c> and probability of success <c>p</c>,
        /// the variance is <c>r * (1 - p) / p^"</c>.</remarks>
        public override double getNumericalVariance()
        {
            double p = getProbabilityOfSuccess();
            double r = getNumberOfSuccesses();
            return r * (1 - p) / (p * p);
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>The lower bound of the support is always 0 no matter the parameters.</remarks>
        public override int getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always <c>Integer.MaxValue</c>
        /// for positive infinity)</returns>
        /// <remarks>The upper bound of the support is always positive infinity no matter the
        /// parameters. Positive infinity is symbolized by <c>Integer.MaxValue</c></remarks>
        public override int getSupportUpperBound()
        {
            return Int32.MaxValue;
        }

        /// <inheritdoc/>
        /// <returns><c>true</c></returns>
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
