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
using Math3.random;
using System;

namespace Math3.distribution
{
    /// <summary>
    /// Implementation of the chi-squared distribution.
    /// </summary>
    /// <remarks>See <a href="http://en.wikipedia.org/wiki/Chi-squared_distribution">Chi-squared distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">Chi-squared Distribution (MathWorld)</a></remarks>
    public class ChiSquaredDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// Internal Gamma distribution.
        /// </summary>
        private readonly GammaDistribution gamma;

        /// <summary>
        /// Inverse cumulative probability accuracy
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Create a Chi-Squared distribution with the given degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        public ChiSquaredDistribution(double degreesOfFreedom) : this(degreesOfFreedom, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create a Chi-Squared distribution with the given degrees of freedom and
        /// inverse cumulative probability accuracy.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        public ChiSquaredDistribution(double degreesOfFreedom, double inverseCumAccuracy) : this(new Well19937c(), degreesOfFreedom, inverseCumAccuracy) { }

        /// <summary>
        /// Create a Chi-Squared distribution with the given degrees of freedom.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        public ChiSquaredDistribution(RandomGenerator rng, double degreesOfFreedom) : this(rng, degreesOfFreedom, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create a Chi-Squared distribution with the given degrees of freedom and
        /// inverse cumulative probability accuracy.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <param name="inverseCumAccuracy">inverseCumAccuracy the maximum absolute error
        /// in inverse cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        public ChiSquaredDistribution(RandomGenerator rng, double degreesOfFreedom, double inverseCumAccuracy)
            : base(rng)
        {
            gamma = new GammaDistribution(degreesOfFreedom / 2, 2);
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <summary>
        /// Access the number of degrees of freedom.
        /// </summary>
        /// <returns>the degrees of freedom.</returns>
        public double getDegreesOfFreedom()
        {
            return gamma.getShape() * 2.0;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            return gamma.density(x);
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            return gamma.logDensity(x);
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(double x)
        {
            return gamma.cumulativeProbability(x);
        }

        /// <inheritdoc/>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>For <c>k</c> degrees of freedom, the mean is <c>k</c>.</remarks>
        public override double getNumericalMean()
        {
            return getDegreesOfFreedom();
        }

        /// <inheritdoc/>
        /// <returns><c>2 * k</c>, where <c>k</c> is the number of degrees of freedom.</returns>
        public override double getNumericalVariance()
        {
            return 2 * getDegreesOfFreedom();
        }

        /// <inheritdoc/>
        /// <remarks>The lower bound of the support is always 0 no matter the
        /// degrees of freedom.</remarks>
        /// <returns>zero</returns>
        public override double getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <remarks>The upper bound of the support is always positive infinity no matter the
        /// degrees of freedom.</remarks>
        /// <returns><c>Double.PositiveInfinity</c></returns>
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
        /// <remarks>The support of this distribution is connected..</remarks>
        /// <returns><c>true</c></returns>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
