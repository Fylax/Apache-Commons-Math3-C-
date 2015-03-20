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
using System;

namespace Math3.random
{
    /// <summary>
    /// Utilities for creating <see cref="RandomGenerator"/> instances.
    /// </summary>
    public class RandomGeneratorFactory
    {
        /// <summary>
        /// Class contains only static methods.
        /// </summary>
        private RandomGeneratorFactory() { }

        /// <summary>
        /// Creates a <see cref="RandomDataGenerator"/> instance that wraps a
        /// Random instance.
        /// </summary>
        /// <param name="rng">Random instance that will generate the
        /// the random data.</param>
        /// <returns>the given RNG, wrapped in a <see cref="RandomGenerator"/>.</returns>
        public static RandomGenerator createRandomGenerator(Random rng)
        {
            return new RandomGeneratorInternal(rng);
        }

        private class RandomGeneratorInternal : RandomGenerator
        {
            private Random rng;

            internal RandomGeneratorInternal(Random rng)
            {
                this.rng = rng;
            }

            /// <inheritdoc/>
            public void setSeed(int seed)
            {
                this.rng = new Random(seed);
            }

            /// <inheritdoc/>
            public void setSeed(int[] seed)
            {
                this.rng = new Random((Int32)convertToLong(seed));
            }

            /// <inheritdoc/>
            public void setSeed(long seed)
            {
                this.rng = new Random((Int32)seed);
            }

            /// <inheritdoc/>
            public void nextBytes(byte[] bytes)
            {
                rng.NextBytes(bytes);
            }

            /// <inheritdoc/>
            public int nextInt()
            {
                return rng.Next();
            }

            /// <inheritdoc/>
            public int nextInt(int n)
            {
                if (n <= 0)
                {
                    throw new NotStrictlyPositiveException<Int32>(n);
                }
                return rng.Next(n);
            }

            /// <inheritdoc/>
            public long nextLong()
            {
                return (Int32)rng.Next();
            }

            /// <inheritdoc/>
            public Boolean nextBoolean()
            {
                return Convert.ToBoolean(rng.Next());
            }

            /// <inheritdoc/>
            public float nextFloat()
            {
                return (Single)rng.NextDouble();
            }

            /// <inheritdoc/>
            public double nextDouble()
            {
                return rng.NextDouble();
            }

            /// <inheritdoc/>
            public double nextGaussian()
            {
                return rng.NextDouble();
            }
        }

        /// <summary>
        /// Converts seed from one representation to another.
        /// </summary>
        /// <param name="seed">Original seed.</param>
        /// <returns>the converted seed.</returns>
        public static long convertToLong(int[] seed)
        {
            // The following number is the largest prime that fits
            // in 32 bits (i.e. 2^32 - 5).
            long prime = 4294967291L;

            long combined = 0L;
            foreach (int s in seed)
            {
                combined = combined * prime + s;
            }

            return combined;
        }
    }
}
