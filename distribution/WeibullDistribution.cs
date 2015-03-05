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
    /// Implementation of the Weibull distribution. This implementation uses the
    /// two parameter form of the distribution defined by
    /// <a href="http://mathworld.wolfram.com/WeibullDistribution.html">
    /// Weibull Distribution</a>, equations (1) and (2).
    /// </summary>
    public class WeibullDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// The shape parameter.
        /// </summary>
        private readonly double shape;

        /// <summary>
        /// The scale parameter.
        /// </summary>
        private readonly double scale;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Cached numerical mean
        /// </summary>
        private double numericalMean = Double.NaN;

        /// <summary>
        /// Whether or not the numerical mean has been calculated
        /// </summary>
        private Boolean numericalMeanIsCalculated = false;

        /// <summary>
        /// Cached numerical variance
        /// </summary>
        private double numericalVariance = Double.NaN;

        /// <summary>
        /// Whether or not the numerical variance has been calculated
        /// </summary>
        private Boolean numericalVarianceIsCalculated = false;

        /// <summary>
        /// Create a Weibull distribution with the given shape and scale and a
        /// location equal to zero.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>   
        /// </summary>
        /// <param name="alpha">Shape parameter.</param>
        /// <param name="beta">Scale parameter.</param>
        /// <exception cref="NotStrictlyPositiveException">if <c>alpha <= 0</c> or
        /// <c>beta <= 0</c>.<exception>
        public WeibullDistribution(double alpha, double beta) : this(alpha, beta, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create a Weibull distribution with the given shape, scale and inverse
        /// cumulative probability accuracy and a location equal to zero.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>   
        /// </summary>
        /// <param name="alpha">Shape parameter.</param>
        /// <param name="beta">Scale parameter.</param>
        /// <param name="inverseCumAccuracy">Maximum absolute error in inverse
        /// cumulative probability estimates
        /// (defaults to <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException">if <c>alpha <= 0</c> or
        /// <c>beta <= 0</c>.<exception>
        public WeibullDistribution(double alpha, double beta, double inverseCumAccuracy) : this(new Well19937c(), alpha, beta, inverseCumAccuracy) { }

        /// <summary>
        /// Creates a Weibull distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="alpha">Shape parameter.</param>
        /// <param name="beta">Scale parameter.</param>
        /// <exception cref="NotStrictlyPositiveException">if <c>alpha <= 0</c> or
        /// <c>beta <= 0</c>.<exception>
        public WeibullDistribution(RandomGenerator rng, double alpha, double beta) : this(rng, alpha, beta, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates a Weibull distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="alpha">Shape parameter.</param>
        /// <param name="beta">Scale parameter.</param>
        /// <param name="inverseCumAccuracy">Maximum absolute error in inverse
        /// cumulative probability estimates
        /// (defaults to <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException">if <c>alpha <= 0</c> or
        /// <c>beta <= 0</c>.<exception>
        public WeibullDistribution(RandomGenerator rng, double alpha, double beta, double inverseCumAccuracy)
            : base(rng)
        {
            if (alpha <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("SHAPE"), alpha);
            }
            if (beta <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("SCALE"), beta);
            }
            scale = beta;
            shape = alpha;
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <summary>
        /// Access the shape parameter, <c>alpha</c>. 
        /// </summary>
        /// <returns>the shape parameter, <c>alpha</c>.</returns>
        public double getShape()
        {
            return shape;
        }

        /// <summary>
        /// Access the scale parameter, <c>beta</c>. 
        /// </summary>
        /// <returns>the scale parameter, <c>beta</c>.</returns>
        public double getScale()
        {
            return scale;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            if (x < 0)
            {
                return 0;
            }

            double xscale = x / scale;
            double xscalepow = FastMath.pow(xscale, shape - 1);

            /*
             * FastMath.pow(x / scale, shape) =
             * FastMath.pow(xscale, shape) =
             * FastMath.pow(xscale, shape - 1) * xscale
             */
            double xscalepowshape = xscalepow * xscale;

            return (shape / scale) * xscalepow * FastMath.exp(-xscalepowshape);
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            if (x < 0)
            {
                return Double.NegativeInfinity;
            }

            double xscale = x / scale;
            double logxscalepow = FastMath.log(xscale) * (shape - 1);

            /*
             * FastMath.pow(x / scale, shape) =
             * FastMath.pow(xscale, shape) =
             * FastMath.pow(xscale, shape - 1) * xscale
             */
            double xscalepowshape = FastMath.exp(logxscalepow) * xscale;

            return FastMath.log(shape / scale) + logxscalepow - xscalepowshape;
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(double x)
        {
            double ret;
            if (x <= 0.0)
            {
                ret = 0.0;
            }
            else
            {
                ret = 1.0 - FastMath.exp(-FastMath.pow(x / scale, shape));
            }
            return ret;
        }

        /// <inheritdoc/>
        /// <remarks>Returns <c>0</c> when <c>p == 0</c> and
        /// <c>Double.PositiveInfinity</c> when <c>p == 1</c>.</remarks>
        public new double inverseCumulativeProbability(double p)
        {
            double ret;
            if (p < 0.0 || p > 1.0)
            {
                throw new OutOfRangeException<Double>(p, 0.0, 1.0);
            }
            else if (p == 0)
            {
                ret = 0.0;
            }
            else if (p == 1)
            {
                ret = Double.PositiveInfinity;
            }
            else
            {
                ret = scale * FastMath.pow(-FastMath.log1p(-p), 1.0 / shape);
            }
            return ret;
        }

        /// <summary>
        /// Return the absolute accuracy setting of the solver used to estimate
        /// inverse cumulative probabilities.
        /// </summary>
        /// <returns>the solver absolute accuracy.</returns>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The mean is <c>scale * Gamma(1 + (1 / shape))</c>, where <c>Gamma()</c>
        /// is the Gamma-function.
        /// </remarks>
        public override double getNumericalMean()
        {
            if (!numericalMeanIsCalculated)
            {
                numericalMean = calculateNumericalMean();
                numericalMeanIsCalculated = true;
            }
            return numericalMean;
        }

        /// <summary>
        /// used by <see cref="getNumericalMean()"/>
        /// </summary>
        /// <returns>the mean of this distribution</returns>
        protected double calculateNumericalMean()
        {
            double sh = getShape();
            double sc = getScale();

            return sc * FastMath.exp(Gamma.logGamma(1 + (1 / sh)));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The variance is <c>scale^2 * Gamma(1 + (2 / shape)) - mean^2</c>
        /// where <c>Gamma()</c> is the Gamma-function.
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
        /// used by <see cref="getNumericalVariance()"/>
        /// </summary>
        /// <returns>the variance of this distribution</returns>
        protected double calculateNumericalVariance()
        {
            double sh = getShape();
            double sc = getScale();
            double mn = getNumericalMean();

            return (sc * sc) * FastMath.exp(Gamma.logGamma(1 + (2 / sh))) -
                   (mn * mn);
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>The lower bound of the support is always 0 no matter the parameters.</remarks>
        public override double getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always
        /// <c>Double.PositiveInfinity</c>)</returns>
        /// <remarks>The upper bound of the support is always positive infinity
        /// no matter the parameters.</remarks>
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
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
