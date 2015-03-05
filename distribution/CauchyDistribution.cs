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
    /// Implementation of the Cauchy distribution.
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Cauchy_distribution">Cauchy distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/CauchyDistribution.html">Cauchy Distribution (MathWorld)</a>
    /// </remarks>
    public class CauchyDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// The median of this distribution.
        /// </summary>
        private readonly double median;

        /// <summary>
        /// The scale of this distribution.
        /// </summary>
        private readonly double scale;

        /// <summary>
        /// Inverse cumulative probability accuracy
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Creates a Cauchy distribution with the median equal to zero and scale
        /// equal to one.
        /// </summary>
        public CauchyDistribution() : this(0, 1) { }

        /// <summary>
        /// Creates a Cauchy distribution using the given median and scale.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.
        /// </summary>
        /// <param name="median">Median for this distribution.</param>
        /// <param name="scale">Scale parameter for this distribution.</param>
        public CauchyDistribution(double median, double scale) : this(median, scale, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates a Cauchy distribution using the given median and scale.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="median">Median for this distribution.</param>
        /// <param name="scale">Scale parameter for this distribution.</param>
        /// <param name="inverseCumAccuracy">Maximum absolute error in inverse
        /// cumulative probability estimates
        /// (defaults to <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>scale <= 0</c>.</exception>
        public CauchyDistribution(double median, double scale,
                                  double inverseCumAccuracy)
            : this(new Well19937c(), median, scale, inverseCumAccuracy) { }

        /// <summary>
        /// Creates a Cauchy distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="median">Median for this distribution.</param>
        /// <param name="scale">Scale parameter for this distribution.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>scale <= 0</c>.</exception>
        public CauchyDistribution(RandomGenerator rng, double median, double scale) : this(rng, median, scale, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates a Cauchy distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="median">Median for this distribution.</param>
        /// <param name="scale">Scale parameter for this distribution.</param>
        /// <param name="inverseCumAccuracy">Maximum absolute error in inverse
        /// cumulative probability estimates
        /// (defaults to <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>scale <= 0</c>.</exception>s
        public CauchyDistribution(RandomGenerator rng, double median, double scale, double inverseCumAccuracy)
            : base(rng)
        {
            if (scale <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("SCALE"), scale);
            }
            this.scale = scale;
            this.median = median;
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(double x)
        {
            return 0.5 + (FastMath.atan((x - median) / scale) / FastMath.PI);
        }

        /// <summary>
        /// Access the median.
        /// </summary>
        /// <returns>the median for this distribution.</returns>
        public double getMedian()
        {
            return median;
        }

        /// <summary>
        /// Access the scale parameter.
        /// </summary>
        /// <returns>the scale parameter for this distribution.</returns>
        public double getScale()
        {
            return scale;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            double dev = x - median;
            return (1 / FastMath.PI) * (scale / (dev * dev + scale * scale));
        }

        /// <inheritdoc/>
        /// <remarks>Returns <c>Double.NegativeInfinity</c> when <c>p == 0</c>
        /// and <c>Double.PositiveInfinity</c> when <c>p == 1</c>.</remarks>
        public new double inverseCumulativeProbability(double p)
        {
            double ret;
            if (p < 0 || p > 1)
            {
                throw new OutOfRangeException<Double>(p, 0, 1);
            }
            else if (p == 0)
            {
                ret = Double.NegativeInfinity;
            }
            else if (p == 1)
            {
                ret = Double.PositiveInfinity;
            }
            else
            {
                ret = median + scale * FastMath.tan(FastMath.PI * (p - .5));
            }
            return ret;
        }

        /// <inheritdoc/>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>The mean is always undefined no matter the parameters.</remarks>
        /// <returns>mean (always Double.NaN)</returns>
        public override double getNumericalMean()
        {
            return Double.NaN;
        }

        /// <inheritdoc/>
        /// <remarks>The variance is always undefined no matter the parameters.</remarks>
        /// <returns>variance (always Double.NaN)</returns>
        public override double getNumericalVariance()
        {
            return Double.NaN;
        }

        /// <inheritdoc/>
        /// <remarks>The lower bound of the support is always negative infinity no matter
        /// the parameters.</remarks>
        /// <returns>lower bound of the support (always Double.NegativeInfinity)</returns>
        public override double getSupportLowerBound()
        {
            return Double.NegativeInfinity;
        }

        /// <inheritdoc/>
        /// <remarks>The upper bound of the support is always positive infinity no matter
        /// the parameters.</remarks>
        /// <returns>upper bound of the support (always Double.PositiveInfinity)</returns>
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
        /// <remarks>The support of this distribution is connected.</remarks>
        /// <returns><c>true</c></returns>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
