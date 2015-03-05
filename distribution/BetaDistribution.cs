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
    /// Implements the Beta distribution.
    /// </summary>
    /// <remarks>See <a href="http://en.wikipedia.org/wiki/Beta_distribution">Beta distribution</a></remarks>
    public class BetaDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// First shape parameter.
        /// </summary>
        private readonly double alpha;

        /// <summary>
        /// Second shape parameter.
        /// </summary>
        private readonly double beta;

        /// <summary>
        /// Normalizing factor used in density computations.
        /// updated whenever alpha or beta are changed.
        /// </summary>
        private double z;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Build a new instance.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="alpha">First shape parameter (must be positive).</param>
        /// <param name="beta">Second shape parameter (must be positive).</param>
        public BetaDistribution(double alpha, double beta) : this(alpha, beta, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Build a new instance.
        /// <para>
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="alpha">First shape parameter (must be positive).</param>
        /// <param name="beta">Second shape parameter (must be positive).</param>
        /// <param name="inverseCumAccuracy">inverseCumAccuracy Maximum absolute error in
        /// inverse cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        public BetaDistribution(double alpha, double beta, double inverseCumAccuracy) : this(new Well19937c(), alpha, beta, inverseCumAccuracy) { }

        /// <summary>
        /// Creates a &beta; distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="alpha">First shape parameter (must be positive).</param>
        /// <param name="beta">Second shape parameter (must be positive).</param>
        public BetaDistribution(RandomGenerator rng, double alpha, double beta) : this(rng, alpha, beta, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates a &beta; distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="alpha">First shape parameter (must be positive).</param>
        /// <param name="beta">Second shape parameter (must be positive).</param>
        /// <param name="inverseCumAccuracy">inverseCumAccuracy Maximum absolute error in
        /// inverse cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        public BetaDistribution(RandomGenerator rng, double alpha, double beta, double inverseCumAccuracy)
            : base(rng)
        {
            this.alpha = alpha;
            this.beta = beta;
            z = Double.NaN;
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <summary>
        /// Access the first shape parameter, <c>alpha</c>. 
        /// </summary>
        /// <returns>the first shape parameter.</returns>
        public double getAlpha()
        {
            return alpha;
        }

        /// <summary>
        /// Access the second shape parameter, <c>beta</c>.
        /// </summary>
        /// <returns>the second shape parameter.</returns>
        public double getBeta()
        {
            return beta;
        }

        /// <summary>
        /// Recompute the normalization factor.
        /// </summary>
        private void recomputeZ()
        {
            if (Double.IsNaN(z))
            {
                z = Gamma.logGamma(alpha) + Gamma.logGamma(beta) - Gamma.logGamma(alpha + beta);
            }
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
            recomputeZ();
            if (x < 0 || x > 1)
            {
                return Double.NegativeInfinity;
            }
            else if (x == 0)
            {
                if (alpha < 1)
                {
                    throw new NumberIsTooSmallException<Double, Int32>(new LocalizedFormats("CANNOT_COMPUTE_BETA_DENSITY_AT_0_FOR_SOME_ALPHA"), alpha, 1, false);
                }
                return Double.NegativeInfinity;
            }
            else if (x == 1)
            {
                if (beta < 1)
                {
                    throw new NumberIsTooSmallException<Double, Int32>(new LocalizedFormats("CANNOT_COMPUTE_BETA_DENSITY_AT_1_FOR_SOME_BETA"), beta, 1, false);
                }
                return Double.NegativeInfinity;
            }
            else
            {
                double logX = FastMath.log(x);
                double log1mX = FastMath.log1p(-x);
                return (alpha - 1) * logX + (beta - 1) * log1mX - z;
            }
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            else if (x >= 1)
            {
                return 1;
            }
            else
            {
                return Beta.regularizedBeta(x, alpha, beta);
            }
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
        /// <remarks>For first shape parameter <c>alpha</c> and second shape parameter
        /// <c>beta</c>, the mean is <c>alpha / (alpha + beta)</c>.</remarks>
        public override double getNumericalMean()
        {
            double a = getAlpha();
            return a / (a + getBeta());
        }

        /// <inheritdoc/>
        /// <remarks>For first shape parameter <c>alpha</c> and second shape parameter
        /// <c>beta</c>, the variance is <c>(alpha * beta) / [(alpha + beta)^2 * (alpha + beta + 1)]</c></remarks>
        public override double getNumericalVariance()
        {
            double a = getAlpha();
            double b = getBeta();
            double alphabetasum = a + b;
            return (a * b) / ((alphabetasum * alphabetasum) * (alphabetasum + 1));
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>The lower bound of the support is always 0 no matter the parameters.</remarks>
        public override double getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always 1)</returns>
        /// <remarks>The upper bound of the support is always 1 no matter the parameters.</remarks>
        public override double getSupportUpperBound()
        {
            return 1;
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
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
