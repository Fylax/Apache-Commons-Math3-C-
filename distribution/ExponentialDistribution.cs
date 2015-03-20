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
    /// Implementation of the exponential distribution.
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Exponential_distribution">Exponential distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/ExponentialDistribution.html">Exponential distribution (MathWorld)</a>
    /// </remarks>
    public class ExponentialDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// Used when generating Exponential samples.
        /// Table containing the constants
        /// q_i = sum_{j=1}^i (ln 2)^j/j! = ln 2 + (ln 2)^2/2 + ... + (ln 2)^i/i!
        /// until the largest representable fraction below 1 is exceeded.
        /// Note that
        /// 1 = 2 - 1 = exp(ln 2) - 1 = sum_{n=1}^infty (ln 2)^n / n!
        /// thus q_i -> 1 as i -> +inf,
        /// so the higher i, the closer to one we get (the series is not alternating).
        /// By trying, n = 16 in Java is enough to reach 1.0.
        /// </summary>
        private static readonly double[] EXPONENTIAL_SA_QI;

        /// <summary>
        /// The mean of this distribution.
        /// </summary>
        private readonly double mean;

        /// <summary>
        /// The logarithm of the mean, stored to reduce computing time.
        /// </summary>
        private readonly double logMean;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Initialize tables.
        /// </summary>
        static ExponentialDistribution()
        {
            /*
             * Filling EXPONENTIAL_SA_QI table.
             * Note that we don't want qi = 0 in the table.
             */
            double LN2 = FastMath.log(2);
            double qi = 0;
            int i = 1;

            /*
             * ArithmeticUtils provides factorials up to 20, so let's use that
             * limit together with Precision.EPSILON to generate the following
             * code (a priori, we know that there will be 16 elements, but it is
             * better to not hardcode it).
             */
            ResizableDoubleArray ra = new ResizableDoubleArray(20);

            while (qi < 1)
            {
                qi += FastMath.pow(LN2, i) / CombinatoricsUtils.factorial(i);
                ra.addElement(qi);
                ++i;
            }

            EXPONENTIAL_SA_QI = ra.getElements();
        }

        /// <summary>
        /// Create an exponential distribution with the given mean.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="mean">mean of this distribution.</param>
        public ExponentialDistribution(double mean) : this(mean, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create an exponential distribution with the given mean.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="mean">Mean of this distribution.</param>
        /// <param name="inverseCumAccuracy">Maximum absolute error in inverse
        /// cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>mean <= 0</c>.</exception>
        public ExponentialDistribution(double mean, double inverseCumAccuracy) : this(new Well19937c(), mean, inverseCumAccuracy) { }

        /// <summary>
        /// Creates an exponential distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="mean">Mean of this distribution.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>mean <= 0</c>.</exception>
        public ExponentialDistribution(RandomGenerator rng, double mean) : this(rng, mean, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates an exponential distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="mean">Mean of this distribution.</param>
        /// <param name="inverseCumAccuracy">Maximum absolute error in inverse
        /// cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>mean <= 0</c>.</exception>
        public ExponentialDistribution(RandomGenerator rng, double mean, double inverseCumAccuracy)
            : base(rng)
        {
            if (mean <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("MEAN"), mean);
            }
            this.mean = mean;
            logMean = FastMath.log(mean);
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <summary>
        /// Access the mean. 
        /// </summary>
        /// <returns>the mean.</returns>
        public double getMean()
        {
            return mean;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            double LogDensity = logDensity(x);
            return LogDensity == Double.NegativeInfinity ? 0 : FastMath.exp(LogDensity);
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            if (x < 0)
            {
                return Double.NegativeInfinity;
            }
            return -x / mean - logMean;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The implementation of this method is based on:
        /// <list type="bullet">
        /// <item>
        /// <a href="http://mathworld.wolfram.com/ExponentialDistribution.html">
        /// Exponential Distribution</a>, equation (1).</item>
        /// </list>
        /// </remarks>
        public override double cumulativeProbability(double x)
        {
            double ret;
            if (x <= 0.0)
            {
                ret = 0.0;
            }
            else
            {
                ret = 1.0 - FastMath.exp(-x / mean);
            }
            return ret;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Returns <c>0</c> when <c>p= = 0</c> and
        /// <c>Double.PositiveInfinity</c> when <c>p == 1</c>.
        /// </remarks>
        public new double inverseCumulativeProbability(double p)
        {
            double ret;

            if (p < 0.0 || p > 1.0)
            {
                throw new OutOfRangeException<Double>(p, 0.0, 1.0);
            }
            else if (p == 1.0)
            {
                ret = Double.PositiveInfinity;
            }
            else
            {
                ret = -mean * FastMath.log(1.0 - p);
            }

            return ret;
        }

        /// <inheritdoc/>
        /// <returns>a random value.</returns>
        /// <remarks>
        /// Algorithm Description: this implementation uses the
        /// <a href="http://www.jesus.ox.ac.uk/~clifford/a5/chap1/node5.html">
        /// Inversion Method</a> to generate exponentially distributed random values
        /// from uniform deviates.
        /// </remarks>
        public new double sample()
        {
            // Step 1:
            double a = 0;
            double u = random.nextDouble();

            // Step 2 and 3:
            while (u < 0.5)
            {
                a += EXPONENTIAL_SA_QI[0];
                u *= 2;
            }

            // Step 4 (now u >= 0.5):
            u += u - 1;

            // Step 5:
            if (u <= EXPONENTIAL_SA_QI[0])
            {
                return mean * (a + u);
            }

            // Step 6:
            int i = 0; // Should be 1, be we iterate before it in while using 0
            double u2 = random.nextDouble();
            double umin = u2;

            // Step 7 and 8:
            do
            {
                ++i;
                u2 = random.nextDouble();

                if (u2 < umin)
                {
                    umin = u2;
                }

                // Step 8:
            } while (u > EXPONENTIAL_SA_QI[i]); // Ensured to exit since EXPONENTIAL_SA_QI[MAX] = 1

            return mean * (a + umin * EXPONENTIAL_SA_QI[0]);
        }

        /// <inheritdoc/>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For mean parameter <c>k</c>, the mean is <c>k</c>.
        /// </remarks>
        public override double getNumericalMean()
        {
            return getMean();
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For mean parameter <c>k</c>, the variance is <c>k^2</c>.
        /// </remarks>
        public override double getNumericalVariance()
        {
            double m = getMean();
            return m * m;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>
        /// The lower bound of the support is always 0 no matter the mean parameter.
        /// </remarks>
        public override double getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always Double.PositiveInfinity)</returns>
        /// <remarks>
        /// The upper bound of the support is always positive infinity
        /// no matter the mean parameter.
        /// </remarks>
        public override double getSupportUpperBound()
        {
            return Double.PositiveInfinity;
        }

        /// <inheritdoc/>
        public override Boolean isSupportLowerBoundInclusive()
        {
            return true;
        }

        /// <inheritdoc/>
        public override Boolean isSupportUpperBoundInclusive()
        {
            return false;
        }

        /// <inheritdoc/>
        /// <returns><c>true</c></returns>
        /// <remarks>
        /// The support of this distribution is connected.
        /// </remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
