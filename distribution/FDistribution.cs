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
    /// Implementation of the F-distribution..
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/F-distribution">F-distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/F-Distribution.html">F-distribution (MathWorld)</a>
    /// </remarks>
    public class FDistribution : AbstractRealDistribution
    {
        /// <summary>
        /// Default inverse cumulative probability accuracy.
        /// </summary>
        public const double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;

        /// <summary>
        /// The numerator degrees of freedom.
        /// </summary>
        private readonly double numeratorDegreesOfFreedom;

        /// <summary>
        /// The numerator degrees of freedom.
        /// </summary>
        private readonly double denominatorDegreesOfFreedom;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Cached numerical variance
        /// </summary>
        private double numericalVariance = Double.NaN;

        /// <summary>
        /// Whether or not the numerical variance has been calculated
        /// </summary>
        private Boolean numericalVarianceIsCalculated = false;

        /// <summary>
        /// Create an F distribution using the given degrees of freedom.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="numeratorDegreesOfFreedom">Numerator degrees of freedom.</param>
        /// <param name="denominatorDegreesOfFreedom">Denominator degrees of freedom.</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>numeratorDegreesOfFreedom <= 0</c> or
        /// <c>denominatorDegreesOfFreedom <= 0</c>.</exception>
        public FDistribution(double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom) : this(numeratorDegreesOfFreedom, denominatorDegreesOfFreedom, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Create an F distribution using the given degrees of freedom.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="numeratorDegreesOfFreedom">Numerator degrees of freedom.</param>
        /// <param name="denominatorDegreesOfFreedom">Denominator degrees of freedom.</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates.</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>numeratorDegreesOfFreedom <= 0</c> or
        /// <c>denominatorDegreesOfFreedom <= 0</c>.</exception>
        public FDistribution(double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom, double inverseCumAccuracy) : this(new Well19937c(), numeratorDegreesOfFreedom, denominatorDegreesOfFreedom, inverseCumAccuracy) { }

        /// <summary>
        /// Creates an F distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="numeratorDegreesOfFreedom">Numerator degrees of freedom.</param>
        /// <param name="denominatorDegreesOfFreedom">Denominator degrees of freedom.</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>numeratorDegreesOfFreedom <= 0</c> or
        /// <c>denominatorDegreesOfFreedom <= 0</c>.</exception>
        public FDistribution(RandomGenerator rng, double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom) : this(rng, numeratorDegreesOfFreedom, denominatorDegreesOfFreedom, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates an F distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="numeratorDegreesOfFreedom">Numerator degrees of freedom.</param>
        /// <param name="denominatorDegreesOfFreedom">Denominator degrees of freedom.</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates.</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>numeratorDegreesOfFreedom <= 0</c> or
        /// <c>denominatorDegreesOfFreedom <= 0</c>.</exception>
        public FDistribution(RandomGenerator rng, double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom, double inverseCumAccuracy)
            : base(rng)
        {
            if (numeratorDegreesOfFreedom <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("DEGREES_OF_FREEDOM"), numeratorDegreesOfFreedom);
            }
            if (denominatorDegreesOfFreedom <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("DEGREES_OF_FREEDOM"), denominatorDegreesOfFreedom);
            }
            this.numeratorDegreesOfFreedom = numeratorDegreesOfFreedom;
            this.denominatorDegreesOfFreedom = denominatorDegreesOfFreedom;
            solverAbsoluteAccuracy = inverseCumAccuracy;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            return FastMath.exp(logDensity(x));
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            double nhalf = numeratorDegreesOfFreedom / 2;
            double mhalf = denominatorDegreesOfFreedom / 2;
            double logx = FastMath.log(x);
            double logn = FastMath.log(numeratorDegreesOfFreedom);
            double logm = FastMath.log(denominatorDegreesOfFreedom);
            double lognxm = FastMath.log(numeratorDegreesOfFreedom * x + denominatorDegreesOfFreedom);
            return nhalf * logn + nhalf * logx - logx +
                   mhalf * logm - nhalf * lognxm - mhalf * lognxm -
                   Beta.logBeta(nhalf, mhalf);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The implementation of this method is based on
        /// <list type="bullet">
        /// <item>
        /// <a href="http://mathworld.wolfram.com/F-Distribution.html">
        /// F-Distribution</a>, equation (4).
        /// </item>
        /// </list>
        /// </remarks>
        public override double cumulativeProbability(double x)
        {
            double ret;
            if (x <= 0)
            {
                ret = 0;
            }
            else
            {
                double n = numeratorDegreesOfFreedom;
                double m = denominatorDegreesOfFreedom;

                ret = Beta.regularizedBeta((n * x) / (m + n * x),
                    0.5 * n,
                    0.5 * m);
            }
            return ret;
        }

        /// <summary>
        /// Access the numerator degrees of freedom. 
        /// </summary>
        /// <returns>the numerator degrees of freedom.</returns>
        public double getNumeratorDegreesOfFreedom()
        {
            return numeratorDegreesOfFreedom;
        }

        /// <summary>
        /// Access the denominator degrees of freedom.
        /// </summary>
        /// <returns>the denominator degrees of freedom.</returns>
        public double getDenominatorDegreesOfFreedom()
        {
            return denominatorDegreesOfFreedom;
        }

        /// <inheritdoc/>
        protected new double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For denominator degrees of freedom parameter <c>b</c>, the mean is
        /// <list type="bullet">
        /// <item>if <c>b > 2</c> then <c>b / (b - 2)</c>,</item>
        /// <item>else undefined (<c>@code Double.NaN</>).</item>
        /// </list>
        /// </remarks>
        public override double getNumericalMean()
        {
            double denominatorDF = getDenominatorDegreesOfFreedom();

            if (denominatorDF > 2)
            {
                return denominatorDF / (denominatorDF - 2);
            }

            return Double.NaN;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For numerator degrees of freedom parameter <c>a</c> and denominator
        /// degrees of freedom parameter <c>b</c>, the variance is
        /// <list type="bullet">
        /// <item>if <c>b > 4</c> then
        /// <c>[2 * b^2 * (a + b - 2)] / [a * (b - 2)^2 * (b - 4)]</c>,</item>
        /// <item>else undefined (<c>@code Double.NaN</>).</item>
        /// </list>
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
            double denominatorDF = getDenominatorDegreesOfFreedom();

            if (denominatorDF > 4)
            {
                double numeratorDF = getNumeratorDegreesOfFreedom();
                double denomDFMinusTwo = denominatorDF - 2;

                return (2 * (denominatorDF * denominatorDF) * (numeratorDF + denominatorDF - 2)) / ((numeratorDF * (denomDFMinusTwo * denomDFMinusTwo) * (denominatorDF - 4)));
            }

            return Double.NaN;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>The lower bound of the support is always 0 no matter the parameters.</remarks>
        public override double getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always Double.PositiveInfinity)</returns>
        /// <remarks>The upper bound of the support is always positive infinity
        /// no matter the parameters.</remarks>
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
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }
    }
}
