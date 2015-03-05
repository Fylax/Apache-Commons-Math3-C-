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
    /// Implementation of the hypergeometric distribution.
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Hypergeometric_distribution">Hypergeometric distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/HypergeometricDistribution.html">Hypergeometric distribution (MathWorld)</a>
    /// </remarks>
    public class HypergeometricDistribution : AbstractIntegerDistribution
    {
        /// <summary>
        /// The number of successes in the population.
        /// </summary>
        private readonly int numberOfSuccesses;

        /// <summary>
        /// The population size.
        /// </summary>
        private readonly int populationSize;

        /// <summary>
        /// The sample size.
        /// </summary>
        private readonly int sampleSize;

        /// <summary>
        /// Cached numerical variance
        /// </summary>
        private double numericalVariance = Double.NaN;

        /// <summary>
        /// Whether or not the numerical variance has been calculated
        /// </summary>
        private Boolean numericalVarianceIsCalculated = false;

        /// <summary>
        /// Construct a new hypergeometric distribution with the specified population
        /// size, number of successes in the population, and sample size.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="populationSize">Population size.</param>
        /// <param name="numberOfSuccesses">Number of successes in the population.</param>
        /// <param name="sampleSize">Sample size.</param>
        /// <exception cref="NotPositiveException"> if <c>numberOfSuccesses < 0</c>.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>populationSize <= 0</c> or <c>scale <= 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException">if <c>numberOfSuccesses > populationSize</c>,
        /// or <c>sampleSize > populationSize</c>.</exception>
        public HypergeometricDistribution(int populationSize, int numberOfSuccesses, int sampleSize) : this(new Well19937c(), populationSize, numberOfSuccesses, sampleSize) { }

        /// <summary>
        /// Creates a new hypergeometric distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="populationSize">Population size.</param>
        /// <param name="numberOfSuccesses">Number of successes in the population.</param>
        /// <param name="sampleSize">Sample size.</param>
        /// <exception cref="NotPositiveException"> if <c>numberOfSuccesses < 0</c>.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>populationSize <= 0</c> or <c>scale <= 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException">if <c>numberOfSuccesses > populationSize</c>,
        /// or <c>sampleSize > populationSize</c>.</exception>
        public HypergeometricDistribution(RandomGenerator rng, int populationSize, int numberOfSuccesses, int sampleSize)
            : base(rng)
        {
            if (populationSize <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("POPULATION_SIZE"), populationSize);
            }
            if (numberOfSuccesses < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_SUCCESSES"), numberOfSuccesses);
            }
            if (sampleSize < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_SAMPLES"), sampleSize);
            }

            if (numberOfSuccesses > populationSize)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("NUMBER_OF_SUCCESS_LARGER_THAN_POPULATION_SIZE"), numberOfSuccesses, populationSize, true);
            }
            if (sampleSize > populationSize)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("SAMPLE_SIZE_LARGER_THAN_POPULATION_SIZE"), sampleSize, populationSize, true);
            }

            this.numberOfSuccesses = numberOfSuccesses;
            this.populationSize = populationSize;
            this.sampleSize = sampleSize;
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(int x)
        {
            double ret;

            int[] domain = getDomain(populationSize, numberOfSuccesses, sampleSize);
            if (x < domain[0])
            {
                ret = 0.0;
            }
            else if (x >= domain[1])
            {
                ret = 1.0;
            }
            else
            {
                ret = innerCumulativeProbability(domain[0], x, 1);
            }

            return ret;
        }

        /// <summary>
        /// Return the domain for the given hypergeometric distribution parameters.
        /// </summary>
        /// <param name="n">Population size.</param>
        /// <param name="m">Number of successes in the population.</param>
        /// <param name="k">Sample size.</param>
        /// <returns>a two element array containing the lower and upper bounds of the
        /// hypergeometric distribution.</returns>
        private int[] getDomain(int n, int m, int k)
        {
            return new int[] { getLowerDomain(n, m, k), getUpperDomain(m, k) };
        }

        /// <summary>
        /// Return the lowest domain value for the given hypergeometric distribution
        /// parameters.
        /// </summary>
        /// <param name="n">Population size.</param>
        /// <param name="m">Number of successes in the population.</param>
        /// <param name="k">Sample size.</param>
        /// <returns>the lowest domain value of the hypergeometric distribution.</returns>
        private int getLowerDomain(int n, int m, int k)
        {
            return FastMath.max(0, m - (n - k));
        }

        /// <summary>
        /// Access the number of successes.
        /// </summary>
        /// <returns>the number of successes.</returns>
        public int getNumberOfSuccesses()
        {
            return numberOfSuccesses;
        }

        /// <summary>
        /// Access the population size.
        /// </summary>
        /// <returns>the population size.</returns>
        public int getPopulationSize()
        {
            return populationSize;
        }

        /// <summary>
        /// Access the sample size.
        /// </summary>
        /// <returns>the sample size.</returns>
        public int getSampleSize()
        {
            return sampleSize;
        }

        /// <summary>
        /// Return the highest domain value for the given hypergeometric distribution
        /// parameters.
        /// </summary>
        /// <param name="m">Number of successes in the population.</param>
        /// <param name="k">Sample size.</param>
        /// <returns>the highest domain value of the hypergeometric distribution.</returns>
        private int getUpperDomain(int m, int k)
        {
            return FastMath.min(k, m);
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
            double ret;

            int[] domain = getDomain(populationSize, numberOfSuccesses, sampleSize);
            if (x < domain[0] || x > domain[1])
            {
                ret = Double.NegativeInfinity;
            }
            else
            {
                double p = (double)sampleSize / (double)populationSize;
                double q = (double)(populationSize - sampleSize) / (double)populationSize;
                double p1 = SaddlePointExpansion.logBinomialProbability(x,
                        numberOfSuccesses, p, q);
                double p2 =
                        SaddlePointExpansion.logBinomialProbability(sampleSize - x,
                                populationSize - numberOfSuccesses, p, q);
                double p3 =
                        SaddlePointExpansion.logBinomialProbability(sampleSize, populationSize, p, q);
                ret = p1 + p2 - p3;
            }

            return ret;
        }

        /// <summary>
        /// For this distribution, <c>X</c>, this method returns <c>P(X >= x)</c>.
        /// </summary>
        /// <param name="x">Value at which the CDF is evaluated.</param>
        /// <returns>the upper tail CDF for this distribution.</returns>
        public double upperCumulativeProbability(int x)
        {
            double ret;

            int[] domain = getDomain(populationSize, numberOfSuccesses, sampleSize);
            if (x <= domain[0])
            {
                ret = 1.0;
            }
            else if (x > domain[1])
            {
                ret = 0.0;
            }
            else
            {
                ret = innerCumulativeProbability(domain[1], x, -1);
            }

            return ret;
        }

        /// <summary>
        /// For this distribution, <c>X</c>, this method returns
        /// <c>P(x0 <= X <= x1)</c>.
        /// This probability is computed by summing the point probabilities for the
        /// values <c>x0, x0 + 1, x0 + 2, ..., x1</c>, in the order directed by
        /// <c>dx</c>.
        /// </summary>
        /// <param name="x0">Inclusive lower bound.</param>
        /// <param name="x1">Inclusive upper bound.</param>
        /// <param name="dx">Direction of summation (1 indicates summing from x0 to x1, and
        /// 0 indicates summing from x1 to x0).</param>
        /// <returns><c>P(x0 <= X <= x1)</c></returns>
        private double innerCumulativeProbability(int x0, int x1, int dx)
        {
            double ret = probability(x0);
            while (x0 != x1)
            {
                x0 += dx;
                ret += probability(x0);
            }
            return ret;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For population size <c>N</c>, number of successes <c>m</c>, and sample
        /// size <c>n</c>, the mean is <c>n * m / N</c>.
        /// </remarks>
        public override double getNumericalMean()
        {
            return getSampleSize() * (getNumberOfSuccesses() / (double)getPopulationSize());
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For population size <c>N</c>, number of successes <c>m</c>, and sample
        /// size <c>n</c>, the variance is
        /// <c>[n * m * (N - n) * (N - m)] / [N^2 * (N - 1)]</c>.
        /// </remarks>
        public override double getNumericalVariance()
        {
            if (!numericalVarianceIsCalculated)
            {
                numericalVariance = calculateNumericalVariance();
                numericalVarianceIsCalculated = true;
            }
            return numericalVariance;
        }

        /// <summary>
        /// Used by <see cref="getNumericalVariance()"/>.
        /// </summary>
        /// <returns>the variance of this distribution</returns>
        protected double calculateNumericalVariance()
        {
            double N = getPopulationSize();
            double m = getNumberOfSuccesses();
            double n = getSampleSize();
            return (n * m * (N - n) * (N - m)) / (N * N * (N - 1));
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support</returns>
        /// <remarks>
        /// For population size <c>N</c>, number of successes <c>m</c>, and sample
        /// <c>n</c>, the lower bound of the support is
        /// <c>max(0, n + m - N)</c>.
        /// </remarks>
        public override int getSupportLowerBound()
        {
            return FastMath.max(0, getSampleSize() + getNumberOfSuccesses() - getPopulationSize());
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support</returns>
        /// <remarks>
        /// For number of successes <c>m</c> and sample size <c>n</c>, the upper
        /// bound of the support is <c>min(m, n)</c>.
        /// </remarks>
        public override int getSupportUpperBound()
        {
            return FastMath.min(getNumberOfSuccesses(), getSampleSize());
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
