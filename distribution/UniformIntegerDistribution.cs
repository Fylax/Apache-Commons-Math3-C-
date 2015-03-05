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
using System;

namespace Math3.distribution
{
    /// <summary>
    /// Implementation of the uniform integer distribution.
    /// </summary>
    /// See <a href="http://en.wikipedia.org/wiki/Uniform_distribution_(discrete)"
    /// >Uniform distribution (discrete), at Wikipedia</a>
    public class UniformIntegerDistribution : AbstractIntegerDistribution
    {
        /// <summary>
        /// Lower bound (inclusive) of this distribution.
        /// </summary>
        private readonly int lower;

        /// <summary>
        /// Upper bound (inclusive) of this distribution.
        /// </summary>
        private readonly int upper;

        /// <summary>
        /// Creates a new uniform integer distribution using the given lower and
        /// upper bounds (both inclusive).
        /// <para>
        /// Note: this constructor will implicitly create an instance of
        /// <see cref="Well19937c"/> as random generator to be used for sampling only (see
        /// <see cref="sample()"/> and <se cref="sample(int)"/>). In case no sampling is
        /// needed for the created distribution, it is advised to pass <c>null</c>
        /// as random generator via the appropriate constructors to avoid the
        /// additional initialisation overhead.</para> 
        /// </summary>
        /// <param name="lower">Lower bound (inclusive) of this distribution.</param>
        /// <param name="upper">Upper bound (inclusive) of this distribution.</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c>.<exception>
        public UniformIntegerDistribution(int lower, int upper) : this(new Well19937c(), lower, upper) { }

        /// <summary>
        /// Creates a new uniform integer distribution using the given lower and
        /// upper bounds (both inclusive).
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="lower">Lower bound (inclusive) of this distribution.</param>
        /// <param name="upper">Upper bound (inclusive) of this distribution.</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c>.<exception>
        public UniformIntegerDistribution(RandomGenerator rng, int lower, int upper)
            : base(rng)
        {
            if (lower > upper)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("LOWER_BOUND_NOT_BELOW_UPPER_BOUND"), lower, upper, true);
            }
            this.lower = lower;
            this.upper = upper;
        }

        /// <inheritdoc/>
        public override double probability(int x)
        {
            if (x < lower || x > upper)
            {
                return 0;
            }
            return 1.0 / (upper - lower + 1);
        }

        /// <inheritdoc/>
        public override double cumulativeProbability(int x)
        {
            if (x < lower)
            {
                return 0;
            }
            if (x > upper)
            {
                return 1;
            }
            return (x - lower + 1.0) / (upper - lower + 1.0);
        }

        /// <inheritdoc/>
        /// <remarks>For lower bound <c>lower</c> and upper bound <c>upper</c>, the mean is
        /// <c>0.5 * (lower + upper)</c>.</remarks>
        public override double getNumericalMean()
        {
            return 0.5 * (lower + upper);
        }

        /// <inheritdoc/>
        /// <remarks>For lower bound <c>lower</c> and upper bound <c>upper</c>, and
        /// <c>n = upper - lower + 1</c>,the variance is <c>(n^2 -1) / 12</c>.</remarks>
        public override double getNumericalVariance()
        {
            double n = upper - lower + 1;
            return (n * n - 1) / 12.0;
        }

        /// <inheritdoc/>
        /// <returns>lower bound of the support</returns>
        /// <remarks>The lower bound of the support is equal to the lower bound parameter
        /// of the distribution.</remarks>
        public override int getSupportLowerBound()
        {
            return lower;
        }

        /// <inheritdoc/>
        /// <returns>upper bound of the support</returns>
        /// <remarks>The upper bound of the support is equal to the upper bound parameter
        /// of the distribution.</remarks>
        public override int getSupportUpperBound()
        {
            return upper;
        }

        /// <inheritdoc/>
        /// <returns><c>true</c></returns>
        /// <remarks>The support of this distribution is connected.</remarks>
        public override Boolean isSupportConnected()
        {
            return true;
        }

        /// <inheritdoc/>
        public new int sample()
        {
            int max = (upper - lower) + 1;
            if (max <= 0)
            {
                // The range is too wide to fit in a positive int (larger
                // than 2^31); as it covers more than half the integer range,
                // we use a simple rejection method.
                while (true)
                {
                    int r = random.nextInt();
                    if (r >= lower &&
                        r <= upper)
                    {
                        return r;
                    }
                }
            }
            else
            {
                // We can shift the range and directly generate a positive int.
                return lower + random.nextInt(max);
            }
        }
    }
}
