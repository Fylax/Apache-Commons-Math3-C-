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
    /// Implementation of the Zipf distribution.
    /// </summary>
    /// <remarks>See <a href="http://mathworld.wolfram.com/ZipfDistribution.html">Zipf distribution (MathWorld)</a></remarks>
    public class ZipfDistribution : AbstractIntegerDistribution
    {
        /// <summary>
        /// Number of elements.
        /// </summary>
        private readonly int numberOfElements;

        /// <summary>
        /// Exponent parameter of the distribution.
        /// </summary>
        private readonly double exponent;

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
        /// Create a new Zipf distribution with the given number of elements and
        /// exponent.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="numberOfElements">Number of elements.</param>
        /// <param name="exponent">Exponent.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>numberOfElements <= 0</c>
        /// or <c>exponent <= 0</c>.<exception>
        public ZipfDistribution(int numberOfElements, double exponent) : this(new Well19937c(), numberOfElements, exponent) { }

        /// <summary>
        /// Creates a Zipf distribution.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="numberOfElements">Number of elements.</param>
        /// <param name="exponent">Exponent.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>numberOfElements <= 0</c>
        /// or <c>exponent <= 0</c>.<exception>
        public ZipfDistribution(RandomGenerator rng, int numberOfElements, double exponent)
            : base(rng)
        {
            if (numberOfElements <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("DIMENSION"),
                                                       numberOfElements);
            }
            if (exponent <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("EXPONENT"), exponent);
            }

            this.numberOfElements = numberOfElements;
            this.exponent = exponent;
        }

        /// <summary>
        /// Get the number of elements (e.g. corpus size) for the distribution.
        /// </summary>
        /// <returns>the number of elements</returns>
        public int getNumberOfElements()
        {
            return numberOfElements;
        }

        /// <summary>
        /// Get the exponent characterizing the distribution. 
        /// </summary>
        /// <returns>the exponent</returns>
        public double getExponent()
        {
            return exponent;
        }

        /// <inheritdoc/>
        public override double probability(int x)
        {
            if (x <= 0 || x > numberOfElements)
            {
                return 0.0;
            }

            return (1.0 / FastMath.pow(x, exponent)) / generalizedHarmonic(numberOfElements, exponent);
        }

        /// <inheritdoc/>
        public new double logProbability(int x)
        {
            if (x <= 0 || x > numberOfElements)
            {
                return Double.NegativeInfinity;
            }

            return -FastMath.log(x) * exponent - FastMath.log(generalizedHarmonic(numberOfElements, exponent));
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(int x)
        {
            if (x <= 0)
            {
                return 0.0;
            }
            else if (x >= numberOfElements)
            {
                return 1.0;
            }

            return generalizedHarmonic(x, exponent) / generalizedHarmonic(numberOfElements, exponent);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For number of elements <c>N</c> and exponent <c>s</c>, the mean is
        /// <c>Hs1 / Hs</c>, where
        /// <list type="bullet">
        /// <item><c>Hs1 = generalizedHarmonic(N, s - 1)</c>,</item>
        /// <item><c>Hs = generalizedHarmonic(N, s)</c>.</item>
        /// </list>
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
        /// Used by <see cref="getNumericalMean()"/>.
        /// </summary>
        /// <returns>the mean of this distribution</returns>
        protected double calculateNumericalMean()
        {
            int N = getNumberOfElements();
            double s = getExponent();

            double Hs1 = generalizedHarmonic(N, s - 1);
            double Hs = generalizedHarmonic(N, s);

            return Hs1 / Hs;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For number of elements <c>N</c> and exponent <c>s</c>, the mean is
        /// <c>(Hs2 / Hs)- (Hs1^2 / Hs^2)</c>, where
        /// <list type="bullet">
        /// <item><c>Hs2 = generalizedHarmonic(N, s - 2)</c>,</item>
        /// <item><c>Hs1 = generalizedHarmonic(N, s - 1)</c>,</item>
        /// <item><c>Hs = generalizedHarmonic(N, s)</c>.</item>
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
        /// Used by <see cref="getNumericalVariance()"/>.
        /// </summary>
        /// <returns>the variance of this distribution</returns>
        protected double calculateNumericalVariance()
        {
            int N = getNumberOfElements();
            double s = getExponent();

            double Hs2 = generalizedHarmonic(N, s - 2);
            double Hs1 = generalizedHarmonic(N, s - 1);
            double Hs = generalizedHarmonic(N, s);

            return (Hs2 / Hs) - ((Hs1 * Hs1) / (Hs * Hs));
        }

        /// <summary>
        /// Calculates the Nth generalized harmonic number. See
        /// <a href="http://mathworld.wolfram.com/HarmonicSeries.html">Harmonic
        /// Series</a>.
        /// </summary>
        /// <param name="n">Term in the series to calculate (must be larger than 1)</param>
        /// <param name="m">Exponent (special case <c>m = 1</c> is the harmonic series).</param>
        /// <returns>the n^{th} generalized harmonic number.</returns>
        private double generalizedHarmonic(int n, double m)
        {
            double value = 0;
            for (int k = n; k > 0; --k)
            {
                value += 1.0 / FastMath.pow(k, m);
            }
            return value;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 1)</returns>
        /// <remarks>The lower bound of the support is always 1 no matter the parameters.</remarks>
        public override int getSupportLowerBound()
        {
            return 1;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support</returns>
        /// <remarks>The upper bound of the support is the number of elements.</remarks>
        public override int getSupportUpperBound()
        {
            return getNumberOfElements();
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
