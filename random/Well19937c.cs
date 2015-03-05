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
using System;

namespace Math3.random
{
    /// <summary>
    /// This class implements the WELL19937c pseudo-random number generator
    /// from Fran&ccedil;ois Panneton, Pierre L'Ecuyer and Makoto Matsumoto.
    /// <para>This generator is described in a paper by Fran&ccedil;ois Panneton,
    /// Pierre L'Ecuyer and Makoto Matsumoto <a
    /// href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng.pdf">Improved
    /// Long-Period Generators Based on Linear Recurrences Modulo 2</a> ACM
    /// Transactions on Mathematical Software, 32, 1 (2006). The errata for the paper
    /// are in <a href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng-errata.txt">wellrng-errata.txt</a>.</para>
    /// </summary>
    /// <remarks>See <a href="http://www.iro.umontreal.ca/~panneton/WELLRNG.html">WELL Random 
    /// number generator</a></remarks>
    public class Well19937c : AbstractWell
    {
        /// <summary>
        /// Number of bits in the pool.
        /// </summary>
        private const int K = 19937;

        /// <summary>
        /// First parameter of the algorithm.
        /// </summary>
        private const int M1 = 70;

        /// <summary>
        /// Second parameter of the algorithm.
        /// </summary>
        private const int M2 = 179;

        /// <summary>
        /// Third parameter of the algorithm.
        /// </summary>
        private const int M3 = 449;

        /// <summary>
        /// Creates a new random number generator.
        /// <para>The instance is initialized using the current time as the
        /// seed.</para>
        /// </summary>
        public Well19937c() : base(K, M1, M2, M3) { }

        /// <summary>
        /// Creates a new random number generator using a single int seed.
        /// </summary>
        /// <param name="seed">the initial seed (32 bits integer)</param>
        public Well19937c(int seed) : base(K, M1, M2, M3, seed) { }

        /// <summary>
        /// Creates a new random number generator using an int array seed.
        /// </summary>
        /// <param name="seed">the initial seed (32 bits integers array), if null
        /// the seed of the generator will be related to the current time</param>
        public Well19937c(int[] seed) : base(K, M1, M2, M3, seed) { }

        /// <summary>
        /// Creates a new random number generator using a single long seed.
        /// </summary>
        /// <param name="seed">the initial seed (64 bits integer)</param>
        public Well19937c(long seed) : base(K, M1, M2, M3, seed) { }

        /// <inheritdoc/>
        protected override int next(int bits)
        {

            int indexRm1 = iRm1[index];
            int indexRm2 = iRm2[index];

            int v0 = v[index];
            int vM1 = v[i1[index]];
            int vM2 = v[i2[index]];
            int vM3 = v[i3[index]];

            int z0 = (Int32.MinValue & v[indexRm1]) ^ (0x7FFFFFFF & v[indexRm2]);
            int z1 = (v0 ^ (v0 << 25)) ^ (vM1 ^ (unchecked((Int32)((UInt32)vM1 >> 27))));
            int z2 = unchecked((Int32)((UInt32)vM2 >> 9)) ^ (vM3 ^ unchecked((Int32)((UInt32)vM3 >> 1)));
            int z3 = z1 ^ z2;
            int z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ unchecked((Int32)((UInt32)z3 >> 21)));

            v[index] = z3;
            v[indexRm1] = z4;
            v[indexRm2] &= Int32.MinValue;
            index = indexRm1;

            // add Matsumoto-Kurita tempering
            // to get a maximally-equidistributed generator
            z4 ^= (z4 << 7) & -0x646e1700;
            z4 ^= (z4 << 15) & -0x1b868000;

            return unchecked((Int32)((UInt32)z4 >> (32 - bits)));
        }
    }
}
