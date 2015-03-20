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
    /// Implementation of the Gamma distribution. 
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Gamma_distribution">Gamma distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/GammaDistribution.html">Gamma distribution (MathWorld)</a>
    /// </remarks>
    public class GammaDistribution : AbstractRealDistribution
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
        /// The constant value of <c>shape + g + 0.5</c>, where <c>g</c> is the
        /// Lanczos constant <see cref="Gamma.LANCZOS_G"/>.
        /// </summary>
        private readonly double shiftedShape;

        /// <summary>
        /// The constant value of
        /// <c>shape / scale * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape)</c>,
        /// where <c>L(shape)</c> is the Lanczos approximation returned by
        /// <see cref="Gamma.lanczos(double)"/>. This prefactor is used in
        /// <see cref="density(double)"/>, when no overflow occurs with the natural
        /// calculation.
        /// </summary>
        private readonly double densityPrefactor1;

        /// <summary>
        /// The constant value of
        /// <c>log(shape / scale * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape))</c>,
        /// where <c>L(shape)</c> is the Lanczos approximation returned by
        /// <see cref="Gamma.lanczos(double)"/>. This prefactor is used in
        /// <see cref="logDensity(double)"/>, when no overflow occurs with the natural
        /// calculation.
        /// </summary>
        private readonly double logDensityPrefactor1;
        
        /// <summary>
        /// The constant value of
        /// <c>shape * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape)</c>,
        /// where <c>L(shape)</c> is the Lanczos approximation returned by
        /// <see cref="Gamma.lanczos(double)"/>. This prefactor is used in
        /// <see cref="density(double)"/>, when overflow occurs with the natural
        /// calculation.
        /// </summary>
        private readonly double densityPrefactor2;

        /// <summary>
        /// The constant value of
        /// <c>log(shape * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape))</c>,
        /// where <c>L(shape)</c> is the Lanczos approximation returned by
        /// <see cref="Gamma.lanczos(double)"/>. This prefactor is used in
        /// <see cref="logDensity(double)"/>, when overflow occurs with the natural
        /// calculation.
        /// </summary>
        private readonly double logDensityPrefactor2;

        /// <summary>
        /// Lower bound on <c>y = x / scale</c> for the selection of the computation
        /// method in <see cref="density(double)"/>. For <c>y <= minY</c>, the natural
        /// calculation overflows.
        /// </summary>
        private readonly double minY;

        /// <summary>
        /// Upper bound on <c>log(y)</c> (<c>y = x / scale</c>) for the selection
        /// of the computation method in <see cref="density(double)"/>. For
        /// <see cref="log(y) >= maxLogY"/>, the natural calculation overflows.
        /// </summary>
        private readonly double maxLogY;

        /// <summary>
        /// Inverse cumulative probability accuracy.
        /// </summary>
        private readonly double solverAbsoluteAccuracy;

        /// <summary>
        /// Creates new gamma distribution with specified values of the shape and
        /// scale parameters.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="shape">the shape parameter</param>
        /// <param name="scale">the scale parameter</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>shape <= 0</c> or <c>scale <= 0</c>.</exception>
        public GammaDistribution(double shape, double scale) : this(shape, scale, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates new gamma distribution with specified values of the shape and
        /// scale parameters.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="shape">the shape parameter</param>
        /// <param name="scale">the scale parameter</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>shape <= 0</c> or <c>scale <= 0</c>.</exception>
        public GammaDistribution(double shape, double scale, double inverseCumAccuracy) : this(new Well19937c(), shape, scale, inverseCumAccuracy) { }

        /// <summary>
        /// Creates a Gamma distribution .
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="shape">the shape parameter</param>
        /// <param name="scale">the scale parameter</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>shape <= 0</c> or <c>scale <= 0</c>.</exception>
        public GammaDistribution(RandomGenerator rng, double shape, double scale) : this(rng, shape, scale, DEFAULT_INVERSE_ABSOLUTE_ACCURACY) { }

        /// <summary>
        /// Creates a Gamma distribution .
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="shape">the shape parameter</param>
        /// <param name="scale">the scale parameter</param>
        /// <param name="inverseCumAccuracy">the maximum absolute error in inverse
        /// cumulative probability estimates (defaults to
        /// <see cref="DEFAULT_INVERSE_ABSOLUTE_ACCURACY"/>).</param>
        /// <exception cref="NotStrictlyPositiveException"> if
        /// <c>shape <= 0</c> or <c>scale <= 0</c>.</exception>
        public GammaDistribution(RandomGenerator rng, double shape, double scale, double inverseCumAccuracy)
            : base(rng)
        {
            if (shape <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("SHAPE"), shape);
            }
            if (scale <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("SCALE"), scale);
            }

            this.shape = shape;
            this.scale = scale;
            this.solverAbsoluteAccuracy = inverseCumAccuracy;
            this.shiftedShape = shape + Gamma.LANCZOS_G + 0.5;
            double aux = FastMath.E / (2.0 * FastMath.PI * shiftedShape);
            this.densityPrefactor2 = shape * FastMath.sqrt(aux) / Gamma.lanczos(shape);
            this.logDensityPrefactor2 = FastMath.log(shape) + 0.5 * FastMath.log(aux) -
                                        FastMath.log(Gamma.lanczos(shape));
            this.densityPrefactor1 = this.densityPrefactor2 / scale *
                    FastMath.pow(shiftedShape, -shape) *
                    FastMath.exp(shape + Gamma.LANCZOS_G);
            this.logDensityPrefactor1 = this.logDensityPrefactor2 - FastMath.log(scale) -
                    FastMath.log(shiftedShape) * shape +
                    shape + Gamma.LANCZOS_G;
            this.minY = shape + Gamma.LANCZOS_G - FastMath.log(Double.MaxValue);
            this.maxLogY = FastMath.log(Double.MaxValue) / (shape - 1.0);
        }

        /// <summary>
        /// Returns the shape parameter of <c>this</c> distribution.
        /// </summary>
        /// <returns>the shape parameter</returns>
        [Obsolete("getShape() should be preferred.")]
        public double getAlpha()
        {
            return shape;
        }

        /// <summary>
        /// Returns the shape parameter of <c>this</c> distribution.
        /// </summary>
        /// <returns>the shape parameter</returns>
        public double getShape()
        {
            return shape;
        }

        /// <summary>
        /// Returns the scale parameter of <c>this</c> distribution.
        /// </summary>
        /// <returns>the scale parameter</returns>
        [Obsolete("getScale() should be preferred.")]
        public double getBeta()
        {
            return scale;
        }

        /// <summary>
        /// Returns the scale parameter of <c>this</c> distribution.
        /// </summary>
        /// <returns>the scale parameter</returns>
        public double getScale()
        {
            return scale;
        }

        /// <inheritdoc/>
        public override double density(double x)
        {
            /* The present method must return the value of
             *
             *     1       x a     - x
             * ---------- (-)  exp(---)
             * x Gamma(a)  b        b
             *
             * where a is the shape parameter, and b the scale parameter.
             * Substituting the Lanczos approximation of Gamma(a) leads to the
             * following expression of the density
             *
             * a              e            1         y      a
             * - sqrt(------------------) ---- (-----------)  exp(a - y + g),
             * x      2 pi (a + g + 0.5)  L(a)  a + g + 0.5
             *
             * where y = x / b. The above formula is the "natural" computation, which
             * is implemented when no overflow is likely to occur. If overflow occurs
             * with the natural computation, the following identity is used. It is
             * based on the BOOST library
             * http://www.boost.org/doc/libs/1_35_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html
             * Formula (15) needs adaptations, which are detailed below.
             *
             *       y      a
             * (-----------)  exp(a - y + g)
             *  a + g + 0.5
             *                              y - a - g - 0.5    y (g + 0.5)
             *               = exp(a log1pm(---------------) - ----------- + g),
             *                                a + g + 0.5      a + g + 0.5
             *
             *  where log1pm(z) = log(1 + z) - z. Therefore, the value to be
             *  returned is
             *
             * a              e            1
             * - sqrt(------------------) ----
             * x      2 pi (a + g + 0.5)  L(a)
             *                              y - a - g - 0.5    y (g + 0.5)
             *               * exp(a log1pm(---------------) - ----------- + g).
             *                                a + g + 0.5      a + g + 0.5
             */
            if (x < 0)
            {
                return 0;
            }
            double y = x / scale;
            if ((y <= minY) || (FastMath.log(y) >= maxLogY))
            {
                /*
                 * Overflow.
                 */
                double aux1 = (y - shiftedShape) / shiftedShape;
                double aux2 = shape * (FastMath.log1p(aux1) - aux1);
                double aux3 = -y * (Gamma.LANCZOS_G + 0.5) / shiftedShape +
                        Gamma.LANCZOS_G + aux2;
                return densityPrefactor2 / x * FastMath.exp(aux3);
            }
            /*
             * Natural calculation.
             */
            return densityPrefactor1 * FastMath.exp(-y) * FastMath.pow(y, shape - 1);
        }

        /// <inheritdoc/>
        public new double logDensity(double x)
        {
            /*
             * see the comment in <see cref="#density(double)"/> for computation details
             */
            if (x < 0)
            {
                return Double.NegativeInfinity;
            }
            double y = x / scale;
            if ((y <= minY) || (FastMath.log(y) >= maxLogY))
            {
                /*
                 * Overflow.
                 */
                double aux1 = (y - shiftedShape) / shiftedShape;
                double aux2 = shape * (FastMath.log1p(aux1) - aux1);
                double aux3 = -y * (Gamma.LANCZOS_G + 0.5) / shiftedShape +
                        Gamma.LANCZOS_G + aux2;
                return logDensityPrefactor2 - FastMath.log(x) + aux3;
            }
            /*
             * Natural calculation.
             */
            return logDensityPrefactor1 - y + FastMath.log(y) * (shape - 1);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The implementation of this method is based on:
        /// <list type="bullet">
        /// <item>
        /// <a href="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
        /// Chi-Squared Distribution</a>, equation (9).
        /// </item>
        /// <item>Casella, G., & Berger, R. (1990). <i>Statistical Inference</i>.
        /// Belmont, CA: Duxbury Press.
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
                ret = Gamma.regularizedGammaP(shape, x / scale);
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
        /// For shape parameter <c>alpha</c> and scale parameter <c>beta</c>, the
        /// mean is <c>alpha * beta</c>.
        /// </remarks>
        public override double getNumericalMean()
        {
            return shape * scale;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// For shape parameter <c>alpha</c> and scale parameter <c>beta</c>, the
        /// variance is <c>alpha * beta^2</c>.
        /// </remarks>
        public override double getNumericalVariance()
        {
            return shape * scale * scale;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>
        /// The lower bound of the support is always 0 no matter the parameters.
        /// </remarks>
        public override double getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always Double.PositiveInfinity)</returns>
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

        /// <summary>
        /// <para>This implementation uses the following algorithms: </para>
        /// <para>For 0 < shape < 1: <para/>
        /// Ahrens, J. H. and Dieter, U., Computer methods for
        /// sampling from gamma, beta, Poisson and binomial distributions.
        /// Computing, 12, 223-246, 1974.</para>
        /// <para>For shape >= 1: <ara/>
        /// Marsaglia and Tsang, A Simple Method for Generating
        /// Gamma Variables. ACM Transactions on Mathematical Software,
        /// Volume 26 Issue 3, September, 2000.</para>
        /// </summary>
        /// <returns>random value sampled from the Gamma(shape, scale) distribution</returns>
        public new double sample()
        {
            if (shape < 1)
            {
                // [1]: p. 228, Algorithm GS

                while (true)
                {
                    // Step 1:
                    double u = random.nextDouble();
                    double bGS = 1 + shape / FastMath.E;
                    double p = bGS * u;

                    if (p <= 1)
                    {
                        // Step 2:

                        double x = FastMath.pow(p, 1 / shape);
                        double u2 = random.nextDouble();

                        if (u2 > FastMath.exp(-x))
                        {
                            // Reject
                            continue;
                        }
                        else
                        {
                            return scale * x;
                        }
                    }
                    else
                    {
                        // Step 3:

                        double x = -1 * FastMath.log((bGS - p) / shape);
                        double u2 = random.nextDouble();

                        if (u2 > FastMath.pow(x, shape - 1))
                        {
                            // Reject
                            continue;
                        }
                        else
                        {
                            return scale * x;
                        }
                    }
                }
            }

            // Now shape >= 1

            double d = shape - 0.333333333333333333;
            double c = 1 / (3 * FastMath.sqrt(d));

            while (true)
            {
                double x = random.nextGaussian();
                double v = (1 + c * x) * (1 + c * x) * (1 + c * x);

                if (v <= 0)
                {
                    continue;
                }

                double x2 = x * x;
                double u = random.nextDouble();

                // Squeeze
                if (u < 1 - 0.0331 * x2 * x2)
                {
                    return scale * d * v;
                }

                if (FastMath.log(u) < 0.5 * x2 + d * (1 - v + FastMath.log(v)))
                {
                    return scale * d * v;
                }
            }
        }
    }
}
