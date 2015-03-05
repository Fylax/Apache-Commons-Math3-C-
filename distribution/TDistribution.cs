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
    /// Implementation of Student's t-distribution.
    /// </summary>
    /// <remarks>
    /// See <a href='http://en.wikipedia.org/wiki/Student&apos;s_t-distribution'>Student's t-distribution (Wikipedia)</a><para/>
    /// See <a href='http://mathworld.wolfram.com/Studentst-Distribution.html'>Student's t-distribution (MathWorld)</a>
    /// </remarks>
    public class TDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// The degrees of freedom.
        /// </summary>
        private readonly double degreesOfFreedom;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Static computation factor based on degreesOfFreedom.
        /// </summary>
        private readonly double factor;

        /// <summary>
        /// Create a t distribution using the given degrees of freedom.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>degreesOfFreedom <= 0</c></exception>
        public TDistribution(double degreesOfFreedom) : this(degreesOfFreedom, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create a t distribution using the given degrees of freedom and the
        /// specified inverse cumulative probability absolute accuracy.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates
        /// (defaults to <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>degreesOfFreedom <= 0</c></exception>
        public TDistribution(double degreesOfFreedom, double inverseCumAccuracy) : this(new Well19937c(), degreesOfFreedom, inverseCumAccuracy) { }

        /// <summary>
        /// Creates a t distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>degreesOfFreedom <= 0</c></exception>
        public TDistribution(RandomGenerator rng, double degreesOfFreedom) : this(rng, degreesOfFreedom, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates a t distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates
        /// (defaults to <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>degreesOfFreedom <= 0</c></exception>
        public TDistribution(RandomGenerator rng, double degreesOfFreedom, double inverseCumAccuracy)
            : base(rng)
        {
            if (degreesOfFreedom <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("DEGREES_OF_FREEDOM"), degreesOfFreedom);
            }
            this.degreesOfFreedom = degreesOfFreedom;
            solverAbsoluteAccuracy = inverseCumAccuracy;

            double n = degreesOfFreedom;
            double nPlus1Over2 = (n + 1) / 2;
            factor = Gamma.logGamma(nPlus1Over2) -
                     0.5 * (FastMath.log(FastMath.PI) + FastMath.log(n)) -
                     Gamma.logGamma(n / 2);
        }

        /// <summary>
        /// Access the degrees of freedom.
        /// </summary>
        /// <returns>the degrees of freedom.</returns>
        public double getDegreesOfFreedom()
        {
            return degreesOfFreedom;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            return FastMath.exp(logDensity(x));
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            double n = degreesOfFreedom;
            double nPlus1Over2 = (n + 1) / 2;
            return factor - nPlus1Over2 * FastMath.log(1 + x * x / n);
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(double x)
        {
            double ret;
            if (x == 0)
            {
                ret = 0.5;
            }
            else
            {
                double t =
                    Beta.regularizedBeta(
                        degreesOfFreedom / (degreesOfFreedom + (x * x)),
                        0.5 * degreesOfFreedom,
                        0.5);
                if (x < 0.0)
                {
                    ret = 0.5 * t;
                }
                else
                {
                    ret = 1.0 - 0.5 * t;
                }
            }

            return ret;
        }

        /// <inheritdoc/>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For degrees of freedom parameter <c>df</c>, the mean is
        /// <list type="bullet">
        /// <item>if <c>df > 1</c> then <c>0</c>,</item>
        /// <item>else undefined (<c>Double.NaN</c>).</item>
        /// </list>
        /// </remarks>
        public override double getNumericalMean()
        {
            double df = getDegreesOfFreedom();

            if (df > 1)
            {
                return 0;
            }

            return Double.NaN;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For degrees of freedom parameter <c>df</c>, the variance is
        /// <list type="bullet">
        /// <item>if <c>df > 2</c> then <c>df / (df - 2)</c>,</item>
        /// <item>if <c>1 < df <= 2</c> then positive infinity
        /// (<c>Double.PositiveInfinity</c>),</item>
        /// <item>else undefined (<c>Double.NaN</c>).</item>
        /// </list>
        /// </remarks>
        public override double getNumericalVariance()
        {
            double df = getDegreesOfFreedom();

            if (df > 2)
            {
                return df / (df - 2);
            }

            if (df > 1 && df <= 2)
            {
                return Double.PositiveInfinity;
            }

            return Double.NaN;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always
        /// <c>Double.NegativeInfinity</c>)</returns>
        /// <remarks>
        /// The lower bound of the support is always negative infinity no matter the
        /// parameters.
        /// </remarks>
        public override double getSupportLowerBound()
        {
            return Double.NegativeInfinity;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always
        /// <c>Double.PositiveInfinity</c>)</returns>
        /// <remarks>
        /// The upper bound of the support is always positive infinity no matter the
        /// parameters.
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
        /// <returns><c>true</c>)</returns>
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
