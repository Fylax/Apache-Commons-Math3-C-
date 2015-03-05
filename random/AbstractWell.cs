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
using Math3.util;
using System;

namespace Math3.random
{
    /// <summary>
    /// This abstract class implements the WELL class of pseudo-random number generator
    /// from Fran&ccedil;ois Panneton, Pierre L'Ecuyer and Makoto Matsumoto.
    /// <para>This generator is described in a paper by Fran&ccedil;ois Panneton,
    /// Pierre L'Ecuyer and Makoto Matsumoto <a
    /// href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng.pdf">Improved
    /// Long-Period Generators Based on Linear Recurrences Modulo 2</a> ACM
    /// Transactions on Mathematical Software, 32, 1 (2006). The errata for the paper
    /// are in <a href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng-errata.txt">wellrng-errata.txt</a>.</para>
    /// </summary>
    /// <remarks>See <a href="http://www.iro.umontreal.ca/~panneton/WELLRNG.html">WELL Random number generator</a></remarks>
    public abstract class AbstractWell : BitsStreamGenerator
    {
        /// <summary>
        /// Current index in the bytes pool.
        /// </summary>
        protected int index;

        /// <summary>
        /// Bytes pool.
        /// </summary>
        protected readonly int[] v;

        /// <summary>
        /// Index indirection table giving for each index its predecessor taking table size into account.
        /// </summary>
        protected readonly int[] iRm1;

        /// <summary>
        /// Index indirection table giving for each index its second predecessor taking table size into account.
        /// </summary>
        protected readonly int[] iRm2;

        /// <summary>
        /// Index indirection table giving for each index the value index + m1 taking table size into account.
        /// </summary>
        protected readonly int[] i1;

        /// <summary>
        /// Index indirection table giving for each index the value index + m2 taking table size into account.
        /// </summary>
        protected readonly int[] i2;

        /// <summary>
        /// Index indirection table giving for each index the value index + m3 taking table size into account.
        /// </summary>
        protected readonly int[] i3;

        /// <summary>
        /// Creates a new random number generator.
        /// <para>The instance is initialized using the current time plus the
        /// system identity hash code of this instance as the seed.</para>
        /// </summary>
        /// <param name="k">number of bits in the pool (not necessarily a multiple of 32)</param>
        /// <param name="m1">first parameter of the algorithm</param>
        /// <param name="m2">second parameter of the algorithm</param>
        /// <param name="m3">third parameter of the algorithm</param>
        protected AbstractWell(int k, int m1, int m2, int m3) : this(k, m1, m2, m3, null) { }

        /// <summary>
        /// Creates a new random number generator using a single int seed. 
        /// </summary>
        /// <param name="k">number of bits in the pool (not necessarily a multiple of 32)</param>
        /// <param name="m1">first parameter of the algorithm</param>
        /// <param name="m2">second parameter of the algorithm</param>
        /// <param name="m3">third parameter of the algorithm</param>
        /// <param name="seed">the initial seed (32 bits integer)</param>
        protected AbstractWell(int k, int m1, int m2, int m3, int seed) : this(k, m1, m2, m3, new int[] { seed }) { }

        /// <summary>
        /// Creates a new random number generator using an int array seed.
        /// </summary>
        /// <param name="k">number of bits in the pool (not necessarily a multiple of 32)</param>
        /// <param name="m1">first parameter of the algorithm</param>
        /// <param name="m2">second parameter of the algorithm</param>
        /// <param name="m3">third parameter of the algorithm</param>
        /// <param name="seed">the initial seed (32 bits integers array), if null
        /// the seed of the generator will be related to the current time</param>
        protected AbstractWell(int k, int m1, int m2, int m3, int[] seed)
        {
            // the bits pool contains k bits, k = r w - p where r is the number
            // of w bits blocks, w is the block size (always 32 in the original paper)
            // and p is the number of unused bits in the last block
            int w = 32;
            int r = (k + w - 1) / w;
            this.v = new int[r];
            this.index = 0;

            // precompute indirection index tables. These tables are used for optimizing access
            // they allow saving computations like "(j + r - 2) % r" with costly modulo operations
            iRm1 = new int[r];
            iRm2 = new int[r];
            i1 = new int[r];
            i2 = new int[r];
            i3 = new int[r];
            for (int j = 0; j < r; ++j)
            {
                iRm1[j] = (j + r - 1) % r;
                iRm2[j] = (j + r - 2) % r;
                i1[j] = (j + m1) % r;
                i2[j] = (j + m2) % r;
                i3[j] = (j + m3) % r;
            }

            // initialize the pool content: causes a CA2214, not my fault
            this.setSeed(seed);

        }

        /// <summary>
        /// Creates a new random number generator using a single long seed.
        /// </summary>
        /// <param name="k">number of bits in the pool (not necessarily a multiple of 32)</param>
        /// <param name="m1">first parameter of the algorithm</param>
        /// <param name="m2">second parameter of the algorithm</param>
        /// <param name="m3">third parameter of the algorithm</param>
        /// <param name="seed">the initial seed (64 bits integer)</param>
        protected AbstractWell(int k, int m1, int m2, int m3, long seed) : this(k, m1, m2, m3, new int[] { unchecked((int)((UInt64)seed >> 32)), (int)(seed & 0xffffffffL) }) { }

        /// <summary>
        /// Reinitialize the generator as if just built with the given int seed.
        /// <para>The state of the generator is exactly the same as a new
        /// generator built with the same seed.</para>
        /// </summary>
        /// <param name="seed">the initial seed (32 bits integer)</param>
        public override void setSeed(int seed)
        {
            setSeed(new int[] { seed });
        }

        /// <summary>
        /// Reinitialize the generator as if just built with the given int array seed.
        /// <para>The state of the generator is exactly the same as a new
        /// generator built with the same seed.</para>
        /// </summary>
        /// <param name="seed">seed the initial seed (32 bits integers array). If null
        /// the seed of the generator will be the system time plus the system identity
        /// hash code of the instance.</param>
        public override void setSeed(int[] seed)
        {
            if (seed == null)
            {
                //Find unix timestamp (seconds since 01/01/1970)
                long ticks = DateTime.UtcNow.Ticks - DateTime.Parse("01/01/1970 00:00:00").Ticks;
                ticks /= 10000; //current timestamp in millis
                setSeed(ticks + this.GetHashCode());
                return;
            }

            Array.Copy(seed, 0, v, 0, FastMath.min(seed.Length, v.Length));

            if (seed.Length < v.Length)
            {
                for (int i = seed.Length; i < v.Length; ++i)
                {
                    long l = v[i - seed.Length];
                    v[i] = (int)((1812433253L * (l ^ (l >> 30)) + i) & 0xffffffffL);
                }
            }

            index = 0;
            clear();  // Clear normal deviate cache
        }

        /// <summary>
        /// Reinitialize the generator as if just built with the given long seed.
        /// <para>The state of the generator is exactly the same as a new
        /// generator built with the same seed.</para> 
        /// </summary>
        /// <param name="seed">the initial seed (64 bits integer)</param>
        public override void setSeed(long seed)
        {
            setSeed(new int[] { unchecked((int)((UInt64)seed >> 32)), (int)(seed & 0xffffffffL) });
        }

        /// <inheritdoc/>
        protected override abstract int next(int bits);
    }
}
