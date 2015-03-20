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
    /// Implementation of the binomial distribution.
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Binomial_distribution">Binomial distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/BinomialDistribution.html">Binomial Distribution (MathWorld)</a>
    /// </remarks>
    public class BinomialDistribution : AbstractIntegerDistribution
    {
        /// <summary>
        /// The number of trials.
        /// </summary>
        private readonly int numberOfTrials;

        /// <summary>
        /// The probability of success.
        /// </summary>
        private readonly double probabilityOfSuccess;

        /// <summary>
        /// Create a binomial distribution with the given number of trials and
        /// probability of success.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="trials">Number of trials.</param>
        /// <param name="p">Probability of success.</param>
        /// <exception cref="NotPositiveException"> if <c>trials < 0</c></exception>.
        /// <exception cref="OutOfRangeException"> if <c>p < 0</c> or <c>p > 1</c>.</exception>
        public BinomialDistribution(int trials, double p) : this(new Well19937c(), trials, p) { }

        /// <summary>
        /// Creates a binomial distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="trials">Number of trials.</param>
        /// <param name="p">Probability of success.</param>
        /// <exception cref="NotPositiveException"> if <c>trials < 0</c>.</exception>
        /// <exception cref="OutOfRangeException"> if <c>p < 0</c> or <c>p > 1</c>.</exception>
        public BinomialDistribution(RandomGenerator rng, int trials, double p)
            : base(rng)
        {
            if (trials < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_TRIALS"), trials);
            }
            if (p < 0 || p > 1)
            {
                throw new OutOfRangeException<Double>(p, 0, 1);
            }

            probabilityOfSuccess = p;
            numberOfTrials = trials;
        }

        /// <summary>
        /// Access the number of trials for this distribution. 
        /// </summary>
        /// <returns>the number of trials.</returns>
        public int getNumberOfTrials()
        {
            return numberOfTrials;
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
            double LogProbability = logProbability(x);
            return LogProbability == Double.NegativeInfinity ? 0 : FastMath.exp(LogProbability);
        }

        /// <inheritdoc/>
        public new double logProbability(int x)
        {
            if (numberOfTrials == 0)
            {
                return (x == 0) ? 0d : Double.NegativeInfinity;
            }
            double ret;
            if (x < 0 || x > numberOfTrials)
            {
                ret = Double.NegativeInfinity;
            }
            else
            {
                ret = SaddlePointExpansion.logBinomialProbability(x,
                        numberOfTrials, probabilityOfSuccess,
                        1.0 - probabilityOfSuccess);
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
            else if (x >= numberOfTrials)
            {
                ret = 1.0;
            }
            else
            {
                ret = 1.0 - Beta.regularizedBeta(probabilityOfSuccess, x + 1.0, numberOfTrials - x);
            }
            return ret;
        }

        /// <inheritdoc/>
        /// <remarks>For <c>n</c> trials and probability parameter <c>p</c>, the mean is
        /// <c>n * p</c>.</remarks>
        public override double getNumericalMean()
        {
            return numberOfTrials * probabilityOfSuccess;
        }

        /// <inheritdoc/>
        /// <remarks>For <c>n</c> trials and probability parameter <c>p</c>, the variance is
        /// <c>n * p * (1 - p)</c>.</remarks>
        public override double getNumericalVariance()
        {
            double p = probabilityOfSuccess;
            return numberOfTrials * p * (1 - p);
        }

        /// <inheritdoc/>
        /// <remarks>The lower bound of the support is always 0 except for the probability
        /// parameter <c>p = 1</c>.</remarks>
        /// <returns>lower bound of the support (0 or the number of trials)</returns>
        public override int getSupportLowerBound()
        {
            return probabilityOfSuccess < 1.0 ? 0 : numberOfTrials;
        }

        /// <inheritdoc/>
        /// <remarks>The upper bound of the support is the number of trials except for the
        /// probability parameter <c>p = 0</c>.</remarks>
        /// <returns>upper bound of the support (number of trials or 0)</returns>
        public override int getSupportUpperBound()
        {
            return probabilityOfSuccess > 0.0 ? numberOfTrials : 0;
        }

        /// <inheritdoc/>
        /// <remarks>he support of this distribution is connected.</remarks>
        /// <returns><c>true</c></returns>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
