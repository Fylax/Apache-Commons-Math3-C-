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
using Math3.analysis;
using Math3.analysis.solvers;
using Math3.exception;
using Math3.exception.util;
using Math3.random;
using Math3.util;
using System;

namespace Math3.distribution
{
    /// <summary>
    /// Base class for probability distributions on the reals.
    /// Default implementations are provided for some of the methods
    /// that do not vary from distribution to distribution.
    /// </summary>
    public abstract class AbstractRealDistribution : RealDistribution
    {
        /// <summary>
        /// Default accuracy.
        /// </summary>
        public const double SOLVER_DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

        /// <summary>
        /// RandomData instance used to generate samples from the distribution.
        /// </summary>
#pragma warning disable 0612
        [Obsolete("Please use random instance variable instead")]
        protected RandomDataImpl randomData = new RandomDataImpl();
#pragma warning restore 0612

        /// <summary>
        /// RNG instance used to generate samples from the distribution.
        /// </summary>
        protected readonly RandomGenerator random;

        /// <summary>
        /// Solver absolute accuracy for inverse cumulative computation
        /// </summary>
        private double solverAbsoluteAccuracy = SOLVER_DEFAULT_ABSOLUTE_ACCURACY;

        [Obsolete("Please use AbstractRealDistribution(RandomGenerator) instead")]
        protected AbstractRealDistribution()
        {
            // Legacy users are only allowed to access the deprecated "randomData".
            // New users are forbidden to use this constructor.
            random = null;
        }

        /// <param name="rng">Random number generator.</param>
        protected AbstractRealDistribution(RandomGenerator rng)
        {
            random = rng;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The default implementation uses the identity
        /// <para><c>P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)</c></para>
        /// </remarks>
        [Obsolete("Please use probability(double,double) instead")]
        public double cumulativeProbability(double x0, double x1)
        {
            return probability(x0, x1);
        }

        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(x0 < X <= x1)</c>.
        /// </summary>
        /// <param name="x0">Lower bound (excluded).</param>
        /// <param name="x1">Upper bound (included).</param>
        /// <returns>the probability that a random variable with this distribution
        /// takes a value between <c>x0</c> and <c>x1</c>, excluding the lower
        /// and including the upper endpoint.</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>x0 > x1</c>.
        /// The default implementation uses the identity
        /// <c>P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)</c>
        public double probability(double x0,
                                  double x1)
        {
            if (x0 > x1)
            {
                throw new NumberIsTooLargeException<Double, Double>(new LocalizedFormats("LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT"), x0, x1, true);
            }
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The default implementation returns
        /// <list type="bullet">
        /// <item><see cref="getSupportLowerBound()"/> for <c>p = 0</c>,</item>
        /// <item><see cref="getSupportUpperBound()"/> for <c>p = 1</c>.</item>
        /// </list>
        /// </remarks>
        public double inverseCumulativeProbability(double p)
        {
            /*
             * IMPLEMENTATION NOTES
             * --------------------
             * Where applicable, use is made of the one-sided Chebyshev inequality
             * to bracket the root. This inequality states that
             * P(X - mu >= k * sig) <= 1 / (1 + k^2),
             * mu: mean, sig: standard deviation. Equivalently
             * 1 - P(X < mu + k * sig) <= 1 / (1 + k^2),
             * F(mu + k * sig) >= k^2 / (1 + k^2).
             *
             * For k = sqrt(p / (1 - p)), we find
             * F(mu + k * sig) >= p,
             * and (mu + k * sig) is an upper-bound for the root.
             *
             * Then, introducing Y = -X, mean(Y) = -mu, sd(Y) = sig, and
             * P(Y >= -mu + k * sig) <= 1 / (1 + k^2),
             * P(-X >= -mu + k * sig) <= 1 / (1 + k^2),
             * P(X <= mu - k * sig) <= 1 / (1 + k^2),
             * F(mu - k * sig) <= 1 / (1 + k^2).
             *
             * For k = sqrt((1 - p) / p), we find
             * F(mu - k * sig) <= p,
             * and (mu - k * sig) is a lower-bound for the root.
             *
             * In cases where the Chebyshev inequality does not apply, geometric
             * progressions 1, 2, 4, ... and -1, -2, -4, ... are used to bracket
             * the root.
             */
            if (p < 0.0 || p > 1.0)
            {
                throw new OutOfRangeException<Double>(p, 0, 1);
            }

            double lowerBound = getSupportLowerBound();
            if (p == 0.0)
            {
                return lowerBound;
            }

            double upperBound = getSupportUpperBound();
            if (p == 1.0)
            {
                return upperBound;
            }

            double mu = getNumericalMean();
            double sig = FastMath.sqrt(getNumericalVariance());
            Boolean chebyshevApplies;
            chebyshevApplies = !(Double.IsInfinity(mu) || Double.IsNaN(mu) ||
                                 Double.IsInfinity(sig) || Double.IsNaN(sig));

            if (lowerBound == Double.NegativeInfinity)
            {
                if (chebyshevApplies)
                {
                    lowerBound = mu - sig * FastMath.sqrt((1d - p) / p);
                }
                else
                {
                    lowerBound = -1.0;
                    while (cumulativeProbability(lowerBound) >= p)
                    {
                        lowerBound *= 2.0;
                    }
                }
            }

            if (upperBound == Double.PositiveInfinity)
            {
                if (chebyshevApplies)
                {
                    upperBound = mu + sig * FastMath.sqrt(p / (1d - p));
                }
                else
                {
                    upperBound = 1.0;
                    while (cumulativeProbability(upperBound) < p)
                    {
                        upperBound *= 2.0;
                    }
                }
            }

            UnivariateFunction toSolve = new UnivariateFunctionHelper(p);

            double x = UnivariateSolverUtils.solve(toSolve,
                                                       lowerBound,
                                                       upperBound,
                                                       getSolverAbsoluteAccuracy());

            if (!isSupportConnected())
            {
                /* Test for plateau. */
                double dx = getSolverAbsoluteAccuracy();
                if (x - dx >= getSupportLowerBound())
                {
                    double px = cumulativeProbability(x);
                    if (cumulativeProbability(x - dx) == px)
                    {
                        upperBound = x;
                        while (upperBound - lowerBound > dx)
                        {
                            double midPoint = 0.5 * (lowerBound + upperBound);
                            if (cumulativeProbability(midPoint) < px)
                            {
                                lowerBound = midPoint;
                            }
                            else
                            {
                                upperBound = midPoint;
                            }
                        }
                        return upperBound;
                    }
                }
            }
            return x;
        }

#pragma warning disable 0618
        private class UnivariateFunctionHelper : AbstractRealDistribution, UnivariateFunction
        {
            private double p;
            internal UnivariateFunctionHelper(double p)
            {
                this.p = p;
            }
            public double value(double x)
            {
                return cumulativeProbability(x) - p;
            }

            public override double density(double x)
            {
                throw new NotImplementedException();
            }

            public override double cumulativeProbability(double x)
            {
                throw new NotImplementedException();
            }

            public override double getNumericalMean()
            {
                throw new NotImplementedException();
            }

            public override double getNumericalVariance()
            {
                throw new NotImplementedException();
            }

            public override double getSupportLowerBound()
            {
                throw new NotImplementedException();
            }

            public override double getSupportUpperBound()
            {
                throw new NotImplementedException();
            }

            public override bool isSupportLowerBoundInclusive()
            {
                throw new NotImplementedException();
            }

            public override bool isSupportUpperBoundInclusive()
            {
                throw new NotImplementedException();
            }

            public override bool isSupportConnected()
            {
                throw new NotImplementedException();
            }
        }
#pragma warning restore 0618

        /// <summary>
        /// Returns the solver absolute accuracy for inverse cumulative computation.
        /// You can override this method in order to use a Brent solver with an
        /// absolute accuracy different from the default.
        /// </summary>
        /// <returns>the maximum absolute error in inverse cumulative probability estimates</returns>
        protected double getSolverAbsoluteAccuracy()
        {
            return solverAbsoluteAccuracy;
        }

#pragma warning disable 0618

        /// <inheritdoc/>
        public void reseedRandomGenerator(long seed)
        {
            random.setSeed(seed);
            randomData.reSeed(seed);
        }
#pragma warning restore 0618

        /// <inheritdoc/>
        /// The default implementation uses the
        /// <a href="http://en.wikipedia.org/wiki/Inverse_transform_sampling">
        /// inversion method.
        /// </remarks>
        public double sample()
        {
            return inverseCumulativeProbability(random.nextDouble());
        }

        /// <inheritdoc/>
        /// The default implementation generates the sample by calling
        /// <see cref="sample()"/> in a loop.
        /// </remarks>
        public double[] sample(int sampleSize)
        {
            if (sampleSize <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_SAMPLES"), sampleSize);
            }
            double[] outp = new double[sampleSize];
            for (int i = 0; i < sampleSize; i++)
            {
                outp[i] = sample();
            }
            return outp;
        }

        /// <inheritdoc/>
        /// <returns>zero.</returns>
        public double probability(double x)
        {
            return 0d;
        }

        /// <summary>
        /// Returns the natural logarithm of the probability density function (PDF) of
        /// this distribution evaluated at the specified point <c>x</c>. In general, the PDF
        /// is the derivative of the <see cref="cumulativeProbability(double)">CDF</see>. 
        /// If the derivative does not exist at <c>x<c>, then an appropriate replacement should
        /// be returned, e.g. <c>Double.PositiveInfinity</c>, <c>Double.NaN</c>, or the limit 
        /// inferior or limit superior of the difference quotient. Note that due to the floating
        /// point precision and under/overflow issues, this method will for some distributions 
        /// be more precise and faster than computing the logarithm of 
        /// <see cref="density(double)"/>. The default implementation simply computes the 
        /// logarithm of <c>density(x)</c>.
        /// </summary>
        /// <param name="x">the point at which the PDF is evaluated</param>
        /// <returns>the logarithm of the value of the probability density function at point <c>x</c></returns>
        public double logDensity(double x)
        {
            return FastMath.log(density(x));
        }

        /// <inheritdoc/>
        public abstract double density(double x);

        /// <inheritdoc/>
        public abstract double cumulativeProbability(double x);

        /// <inheritdoc/>
        public abstract double getNumericalMean();

        /// <inheritdoc/>
        public abstract double getNumericalVariance();

        /// <inheritdoc/>
        public abstract double getSupportLowerBound();

        /// <inheritdoc/>
        public abstract double getSupportUpperBound();

        /// <inheritdoc/>
        public abstract bool isSupportLowerBoundInclusive();

        /// <inheritdoc/>
        public abstract bool isSupportUpperBoundInclusive();

        /// <inheritdoc/>
        public abstract bool isSupportConnected();
    }
}
