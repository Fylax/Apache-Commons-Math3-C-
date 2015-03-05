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
    /// Implementation of the Poisson distribution. 
    /// </summary>
    /// <remarks>
    /// See <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution (Wikipedia)</a><para/>
    /// See <a href="http://mathworld.wolfram.com/PoissonDistribution.html">Poisson distribution (MathWorld)</a>
    /// </remarks>
    public class PoissonDistribution : AbstractIntegerDistribution
    {
        /// <summary>
        /// Default maximum number of iterations for cumulative probability calculations.
        /// </summary>
        public const int DEFAULT_MAX_ITERATIONS = 10000000;

        /// <summary>
        /// Default convergence criterion.
        /// </summary>
        public const double DEFAULT_EPSILON = 1e-12;

        /// <summary>
        /// Distribution used to compute normal approximation.
        /// </summary>
        private readonly NormalDistribution normal;

        /// <summary>
        /// Distribution needed for the <see cref="sample()"/> method.
        /// </summary>
        private readonly ExponentialDistribution exponential;

        /// <summary>
        /// Mean of the distribution.
        /// </summary>
        private readonly double mean;

        /// <summary>
        /// Maximum number of iterations for cumulative probability. Cumulative
        /// probabilities are estimated using either Lanczos series approximation
        /// of <see cref="Gamma.regularizedGammaP(double, double, double, int)"/>
        /// or continued fraction approximation of
        /// <see cref="Gamma.regularizedGammaQ(double, double, double, int)"/>.
        /// </summary>
        private readonly int maxIterations;

        /// <summary>
        /// Convergence criterion for cumulative probability.
        /// </summary>
        private readonly double epsilon;

        /// <summary>
        /// Creates a new Poisson distribution with specified mean.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="p">the Poisson mean</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>p <= 0</c>.</exception>
        public PoissonDistribution(double p) : this(p, DEFAULT_EPSILON, DEFAULT_MAX_ITERATIONS) { }

        /// <summary>
        /// Creates a new Poisson distribution with specified mean, convergence
        /// criterion and maximum number of iterations.
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <see cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para>
        /// </summary>
        /// <param name="p">the Poisson mean</param>
        /// <param name="epsilon">Convergence criterion for cumulative probabilities.</param>
        /// <param name="maxIterations">the maximum number of iterations for cumulative
        /// probabilities.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>p <= 0</c>.</exception>
        public PoissonDistribution(double p, double epsilon, int maxIterations) : this(new Well19937c(), p, epsilon, maxIterations) { }

        /// <summary>
        /// Creates a new Poisson distribution with specified mean, convergence
        /// criterion and maximum number of iterations..
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="p">the Poisson mean</param>
        /// <param name="epsilon">Convergence criterion for cumulative probabilities.</param>
        /// <param name="maxIterations">the maximum number of iterations for cumulative
        /// probabilities.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>p <= 0</c>.</exception>
        public PoissonDistribution(RandomGenerator rng, double p, double epsilon, int maxIterations)
            : base(rng)
        {
            if (p <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("MEAN"), p);
            }
            mean = p;
            this.epsilon = epsilon;
            this.maxIterations = maxIterations;

            // Use the same RNG instance as the parent class.
            normal = new NormalDistribution(rng, p, FastMath.sqrt(p),
                                            NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
            exponential = new ExponentialDistribution(rng, 1, ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
        }

        /// <summary>
        /// Creates a new Poisson distribution with the specified mean and
        /// convergence criterion.
        /// </summary>
        /// <param name="p">the Poisson mean</param>
        /// <param name="epsilon">Convergence criterion for cumulative probabilities.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>p <= 0</c>.</exception>
        public PoissonDistribution(double p, double epsilon) : this(p, epsilon, DEFAULT_MAX_ITERATIONS) { }

        /// <summary>
        /// Creates a new Poisson distribution with the specified mean and maximum
        /// number of iterations.
        /// </summary>
        /// <param name="p">the Poisson mean</param>
        /// <param name="maxIterations">the maximum number of iterations for cumulative
        /// probabilities.</param>
        /// <exception cref="NotStrictlyPositiveException"> if <c>p <= 0</c>.</exception>
        public PoissonDistribution(double p, int maxIterations) : this(p, DEFAULT_EPSILON, maxIterations) { }

        /// <summary>
        /// Get the mean for the distribution.
        /// </summary>
        /// <returns>the mean for the distribution.</returns>
        public double getMean()
        {
            return mean;
        }

        /// <inheritdoc/>
        public override double probability(int x)
        {
            double LogProbability = logProbability(x);
            return LogProbability == Double.NegativeInfinity ? 0 : FastMath.exp(LogProbability);
        }

        /// <inheritdoc/>
        public new double logProbability(int x)
        {
            double ret;
            if (x < 0 || x == Int32.MaxValue)
            {
                ret = Double.NegativeInfinity;
            }
            else if (x == 0)
            {
                ret = -mean;
            }
            else
            {
                ret = -SaddlePointExpansion.getStirlingError(x) -
                      SaddlePointExpansion.getDeviancePart(x, mean) -
                      0.5 * FastMath.log(MathUtils.TWO_PI) - 0.5 * FastMath.log(x);
            }
            return ret;
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(int x)
        {
            if (x < 0)
            {
                return 0;
            }
            if (x == Int32.MaxValue)
            {
                return 1;
            }
            return Gamma.regularizedGammaQ((double)x + 1, mean, epsilon, maxIterations);
        }

        /// <summary>
        /// Calculates the Poisson distribution function using a normal
        /// approximation. The <c>N(mean, sqrt(mean))</c> distribution is used
        /// to approximate the Poisson distribution. The computation uses
        /// "half-correction" (evaluating the normal distribution function at
        /// <c>x + 0.5</c>).
        /// </summary>
        /// <param name="x">Upper bound, inclusive.</param>
        /// <returns>the distribution function value calculated using a normal
        /// approximation.</returns>
        public double normalApproximateProbability(int x)
        {
            // calculate the probability using half-correction
            return normal.cumulativeProbability(x + 0.5);
        }

        /// <inheritdoc/>
        /// <remarks>For mean parameter <c>p</c>, the mean is <c>p</c>.</remarks>
        public override double getNumericalMean()
        {
            return getMean();
        }

        /// <inheritdoc/>
        /// <remarks>For mean parameter <c>p</c>, the variance is <c>p</c>.</remarks>
        public override double getNumericalVariance()
        {
            return getMean();
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support (always 0)</returns>
        /// <remarks>The lower bound of the support is always 0 no matter the mean parameter.</remarks>
        public override int getSupportLowerBound()
        {
            return 0;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support (always <c>Int32.MaxValue</c> for
        /// positive infinity)</returns>
        /// <remarks>The upper bound of the support is positive infinity,
        /// regardless of the parameter values. There is no integer infinity,
        /// so this method returns <c>Int32.MaxValue</c>.</remarks>
        public override int getSupportUpperBound()
        {
            return Int32.MaxValue;
        }

        /// <inheritdoc/>
        /// <returns><c>true</c></returns>
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }

        /// <inheritdoc/>
        /// <returns>a random value.</returns>
        /// <remarks>Algorithm Description:
        /// <list type="bullet">
        /// <item>For small means, uses simulation of a Poisson process
        /// using Uniform deviates, as described
        /// <a href="http://irmi.epfl.ch/cmos/Pmmi/interactive/rng7.htm"> here</a>.
        /// The Poisson process (and hence value returned) is bounded by 1000 * mean.
        /// </item>
        /// <item>For large means, uses the rejection algorithm described in
        /// <para/>
        /// Devroye, Luc. (1981). The Computer Generation of Poisson Random Variables
        /// Computing vol. 26 pp. 197-207.
        /// </item>
        /// </list></remarks>
        public new int sample()
        {
            return (int)FastMath.min(nextPoisson(mean), Int32.MaxValue);
        }

        /// <param name="meanPoisson">Mean of the Poisson distribution.</param>
        /// <returns>the next sample.</returns>
        private long nextPoisson(double meanPoisson)
        {
            double pivot = 40.0d;
            if (meanPoisson < pivot)
            {
                double p = FastMath.exp(-meanPoisson);
                long n = 0;
                double r = 1.0d;
                double rnd = 1.0d;

                while (n < 1000 * meanPoisson)
                {
                    rnd = random.nextDouble();
                    r *= rnd;
                    if (r >= p)
                    {
                        n++;
                    }
                    else
                    {
                        return n;
                    }
                }
                return n;
            }
            else
            {
                double lambda = FastMath.floor(meanPoisson);
                double lambdaFractional = meanPoisson - lambda;
                double logLambda = FastMath.log(lambda);
                double logLambdaFactorial = CombinatoricsUtils.factorialLog((int)lambda);
                long y2 = lambdaFractional < Double.MinValue ? 0 : nextPoisson(lambdaFractional);
                double delta = FastMath.sqrt(lambda * FastMath.log(32 * lambda / FastMath.PI + 1));
                double halfDelta = delta / 2;
                double twolpd = 2 * lambda + delta;
                double a1 = FastMath.sqrt(FastMath.PI * twolpd) * FastMath.exp(1 / (8 * lambda));
                double a2 = (twolpd / delta) * FastMath.exp(-delta * (1 + delta) / twolpd);
                double aSum = a1 + a2 + 1;
                double p1 = a1 / aSum;
                double p2 = a2 / aSum;
                double c1 = 1 / (8 * lambda);

                double x = 0;
                double y = 0;
                double v = 0;
                int a = 0;
                double t = 0;
                double qr = 0;
                double qa = 0;
                for (; ; )
                {
                    double u = random.nextDouble();
                    if (u <= p1)
                    {
                        double n = random.nextGaussian();
                        x = n * FastMath.sqrt(lambda + halfDelta) - 0.5d;
                        if (x > delta || x < -lambda)
                        {
                            continue;
                        }
                        y = x < 0 ? FastMath.floor(x) : FastMath.ceil(x);
                        double e = exponential.sample();
                        v = -e - (n * n / 2) + c1;
                    }
                    else
                    {
                        if (u > p1 + p2)
                        {
                            y = lambda;
                            break;
                        }
                        else
                        {
                            x = delta + (twolpd / delta) * exponential.sample();
                            y = FastMath.ceil(x);
                            v = -exponential.sample() - delta * (x + 1) / twolpd;
                        }
                    }
                    a = x < 0 ? 1 : 0;
                    t = y * (y + 1) / (2 * lambda);
                    if (v < -t && a == 0)
                    {
                        y = lambda + y;
                        break;
                    }
                    qr = t * ((2 * y + 1) / (6 * lambda) - 1);
                    qa = qr - (t * t) / (3 * (lambda + a * (y + 1)));
                    if (v < qa)
                    {
                        y = lambda + y;
                        break;
                    }
                    if (v > qr)
                    {
                        continue;
                    }
                    if (v < y * logLambda - CombinatoricsUtils.factorialLog((int)(y + lambda)) + logLambdaFactorial)
                    {
                        y = lambda + y;
                        break;
                    }
                }
                return y2 + (long)y;
            }
        }
    }
}
