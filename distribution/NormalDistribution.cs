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
    /// Implementation of the normal (gaussian) distribution. 
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Normal_distribution">Normal distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/NormalDistribution.html">Normal distribution (MathWorld)</a>
    /// </remarks>
    public class NormalDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// &radic;(2)
        /// </summary>
        private static readonly double SQRT2 = FastMath.sqrt(2.0);

        /// <summary>
        /// Mean of this distribution.
        /// </summary>
        private readonly double mean;

        /// <summary>
        /// Standard deviation of this distribution.
        /// </summary>
        private readonly double standardDeviation;

        /// <summary>
        /// The value of <c>log(sd) + 0.5*log(2*pi)</c> stored for faster computation.
        /// </summary>
        private readonly double logStandardDeviationPlusHalfLog2Pi;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Create a normal distribution with mean equal to zero and standard
        /// deviation equal to one.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        public NormalDistribution() : this(0, 1) { }

        /// <summary>
        /// Create a normal distribution with mean equal to zero and standard
        /// deviation equal to one.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="mean">Mean for this distribution.</param>
        /// <param name="sd">Standard deviation for this distribution.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sd <= 0</c>.</exception>
        public NormalDistribution(double mean, double sd) : this(mean, sd, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create a normal distribution with mean equal to zero and standard
        /// deviation equal to one.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="mean">Mean for this distribution.</param>
        /// <param name="sd">Standard deviation for this distribution.</param>
        /// <param name="inverseCumAccuracy">Inverse cumulative probability accuracy.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sd <= 0</c>.</exception>
        public NormalDistribution(double mean, double sd, double inverseCumAccuracy) : this(new Well19937c(), mean, sd, inverseCumAccuracy) { }

        /// <summary>Creates a normal distribution.</summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="mean">Mean for this distribution.</param>
        /// <param name="sd">Standard deviation for this distribution.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sd <= 0</c>.</exception>
        public NormalDistribution(RandomGenerator rng, double mean, double sd) : this(rng, mean, sd, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>Creates a normal distribution.</summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="mean">Mean for this distribution.</param>
        /// <param name="sd">Standard deviation for this distribution.</param>
        /// <param name="inverseCumAccuracy">Inverse cumulative probability accuracy.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sd <= 0</c>.</exception>
        public NormalDistribution(RandomGenerator rng, double mean, double sd, double inverseCumAccuracy)
            : base(rng)
        {
            if (sd <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("STANDARD_DEVIATION"), sd);
            }

            this.mean = mean;
            standardDeviation = sd;
            logStandardDeviationPlusHalfLog2Pi = FastMath.log(sd) + 0.5 * FastMath.log(2 * FastMath.PI);
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <summary>
        /// Access the mean.
        /// </summary>
        /// <returns>the mean for this distribution</returns>
        public double getMean()
        {
            return mean;
        }

        /// <summary>
        /// Access the standard deviation.
        /// </summary>
        /// <returns>the standard deviation for this distribution.</returns>
        public double getStandardDeviation()
        {
            return standardDeviation;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            return FastMath.exp(logDensity(x));
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            double x0 = x - mean;
            double x1 = x0 / standardDeviation;
            return -0.5 * x1 * x1 - logStandardDeviationPlusHalfLog2Pi;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// If <c>x</c> is more than 40 standard deviations from the mean, 0 or 1
        /// is returned, as in these cases the actual value is within
        /// <c>Double.MinValue</c> of 0 or 1.
        /// </remarks>
        public override double cumulativeProbability(double x)
        {
            double dev = x - mean;
            if (FastMath.abs(dev) > 40 * standardDeviation)
            {
                return dev < 0 ? 0.0d : 1.0d;
            }
            return 0.5 * (1 + Erf.erf(dev / (standardDeviation * SQRT2)));
        }

        /// <inheritdoc/>
        public new double inverseCumulativeProbability(double p)
        {
            if (p < 0.0 || p > 1.0)
            {
                throw new OutOfRangeException<Double>(p, 0, 1);
            }
            return mean + standardDeviation * SQRT2 * Erf.erfInv(2 * p - 1);
        }

        /// <inheritdoc/>
        [Obsolete("See RealDistribution#cumulativeProbability(double,double)")]
        public new double cumulativeProbability(double x0, double x1)
        {
            return probability(x0, x1);
        }

        /// <inheritdoc/>
        public new double probability(double x0, double x1)
        {
            if (x0 > x1)
            {
                throw new NumberIsTooLargeException<Double, Double>(new LocalizedFormats("LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT"), x0, x1, true);
            }
            double denom = standardDeviation * SQRT2;
            double v0 = (x0 - mean) / denom;
            double v1 = (x1 - mean) / denom;
            return 0.5 * Erf.erf(v0, v1);
        }

        /// <inheritdoc/>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For mean parameter <c>mu</c>, the mean is <c>mu</c>.
        /// </remarks>
        public override double getNumericalMean()
        {
            return getMean();
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For standard deviation parameter <c>s</c>, the variance is <c>s^2</c>.
        /// </remarks>
        public override double getNumericalVariance()
        {
            double s = getStandardDeviation();
            return s * s;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always
        /// <c>Double.NegativeInfinity</c></returns>
        /// <remarks>
        /// The lower bound of the support is always negative infinity
        /// no matter the parameters.
        /// </remarks>
        public override double getSupportLowerBound()
        {
            return Double.NegativeInfinity;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always
        /// <c>Double.PositiveInfinity</c></returns>
        /// <remarks>
        /// The upper bound of the support is always positive infinity
        /// no matter the parameters.
        /// </remarks>
        public override double getSupportUpperBound()
        {
            return Double.PositiveInfinity;
        }

        /// <inheritdoc/>
        public override Boolean isSupportLowerBoundInclusive()
        {
            return false;
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

        /// <inheritdoc/>
        public new double sample()
        {
            return standardDeviation * random.nextGaussian() + mean;
        }
    }
}
