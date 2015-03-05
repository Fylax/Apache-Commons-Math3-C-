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
using Math3.util;
using System;

namespace Math3.random
{
    /// <summary>
    /// Base class for random number generators that generates bits streams.
    /// </summary>
    public abstract class BitsStreamGenerator : RandomGenerator
    {
        /// <summary>
        /// Next gaussian.
        /// </summary>
        private double NextGaussian;

        /// <summary>
        /// Creates a new random number generator.
        /// </summary>
        public BitsStreamGenerator()
        {
            NextGaussian = Double.NaN;
        }

        /// <inheritdoc/>
        public abstract void setSeed(int seed);

        /// <inheritdoc/>
        public abstract void setSeed(int[] seed);

        /// <inheritdoc/>
        public abstract void setSeed(long seed);

        /// <summary>
        /// Generate next pseudorandom number.
        /// <para>This method is the core generation algorithm. It is used by all the
        /// public generation methods for the various primitive types 
        /// <see cref="nextBoolean()"/>, <see cref="nextBytes(byte[])"/>, <see cref="nextDouble()"/>,
        /// <see cref="nextFloat()"/>, <see cref="nextGaussian()"/>, <see cref="nextInt()"/>,
        /// <see cref="next(int)"/> and <see cref="nextLong()"/>.</para>
        /// </summary>
        /// <param name="bits">number of random bits to produce</param>
        /// <returns>random bits generated</returns>
        protected abstract int next(int bits);

        /// <inheritdoc/>
        public Boolean nextBoolean()
        {
            return next(1) != 0;
        }

        /// <inheritdoc/>
        public void nextBytes(byte[] bytes)
        {
            int i = 0;
            int iEnd = bytes.Length - 3;
            int random;
            while (i < iEnd)
            {
                random = next(32);
                bytes[i] = (byte)(random & 0xff);
                bytes[i + 1] = (byte)((random >> 8) & 0xff);
                bytes[i + 2] = (byte)((random >> 16) & 0xff);
                bytes[i + 3] = (byte)((random >> 24) & 0xff);
                i += 4;
            }
            random = next(32);
            while (i < bytes.Length)
            {
                bytes[i++] = (byte)(random & 0xff);
                random >>= 8;
            }
        }

        /// <inheritdoc/>
        public double nextDouble()
        {
            long high = ((long)next(26)) << 26;
            long low = next(26);
            return (high | low) * (0x1 * Math.Pow(2, -52d));
        }

        /// <inheritdoc/>
        public float nextFloat()
        {
            return (Single)next(23) * (0x1f * (Single)Math.Pow(2, -23f));
        }

        /// <inheritdoc/>
        public double nextGaussian()
        {
            double random;
            if (Double.IsNaN(NextGaussian))
            {
                // generate a new pair of gaussian numbers
                double x = nextDouble();
                double y = nextDouble();
                double alpha = 2 * FastMath.PI * x;
                double r = FastMath.sqrt(-2 * FastMath.log(y));
                random = r * FastMath.cos(alpha);
                NextGaussian = r * FastMath.sin(alpha);
            }
            else
            {
                // use the second element of the pair already generated
                random = NextGaussian;
                NextGaussian = Double.NaN;
            }

            return random;

        }

        /// <inheritdoc/>
        public int nextInt()
        {
            return next(32);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>This default implementation is copied from Apache Harmony</para>
        /// <para>Implementation notes: 
        /// <list type="bullet">
        /// <item>If n is a power of 2, this method returns
        /// <c>(int) ((n * (long) next(31)) >> 31)</c>.</item>
        /// <item>If n is not a power of 2, what is returned is <c>next(31) % n</c>
        /// with <c>next(31)</c> values rejected (i.e. regenerated) until a
        /// value that is larger than the remainder of <c>Int32.MaxValue / n</c>
        /// is generated. Rejection of this initial segment is necessary to ensure
        /// a uniform distribution.</item>
        /// </list></para>
        /// </remarks>
        public int nextInt(int n)
        {
            if (n > 0)
            {
                if ((n & -n) == n)
                {
                    return (int)((n * (long)next(31)) >> 31);
                }
                int bits;
                int val;
                do
                {
                    bits = next(31);
                    val = bits % n;
                } while (bits - val + (n - 1) < 0);
                return val;
            }
            throw new NotStrictlyPositiveException<Int32>(n);
        }

        /// <inheritdoc/>
        public long nextLong()
        {
            long high = ((long)next(32)) << 32;
            long low = ((long)next(32)) & 0xffffffffL;
            return high | low;
        }

        /// <summary>
        /// Returns a pseudorandom, uniformly distributed <c>long</c> value
        /// between 0 (inclusive) and the specified value (exclusive), drawn from
        /// this random number generator's sequence. 
        /// </summary>
        /// <param name="n">the bound on the random number to be returned.  Must be
        /// positive.</param>
        /// <returns>a pseudorandom, uniformly distributed <c>long</c>
        /// value between 0 (inclusive) and n (exclusive).</returns>
        /// <exception cref="IllegalArgumentException"> if n is not positive.<exception>
        public long nextLong(long n)
        {
            if (n > 0)
            {
                long bits;
                long val;
                do
                {
                    bits = ((long)next(31)) << 32;
                    bits |= ((long)next(32)) & 0xffffffffL;
                    val = bits % n;
                } while (bits - val + (n - 1) < 0);
                return val;
            }
            throw new NotStrictlyPositiveException<Int64>(n);
        }

        /// <summary>
        /// Clears the cache used by the default implementation of
        /// <see cref="nextGaussian"/>.
        /// </summary>
        public void clear()
        {
            NextGaussian = Double.NaN;
        }
    }
}
