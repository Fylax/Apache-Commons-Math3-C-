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
    /// Interface extracted from <c>java.util.Random</c>.  This interface is
    /// implemented by <see cref="AbstractRandomGenerator"/>.
    /// </summary>
    public interface RandomGenerator
    {
        /// <summary>
        /// Sets the seed of the underlying random number generator using an
        /// <c>int</c> seed.
        /// <para>Sequences of values generated starting with the same seeds
        /// should be identical.
        /// </para>
        /// </summary>
        /// <param name="seed">the seed value</param>
        void setSeed(int seed);

        /// <summary>
        /// Sets the seed of the underlying random number generator using an
        /// <c>int</c> array seed.
        /// <para>Sequences of values generated starting with the same seeds
        /// should be identical.
        /// </para>
        /// </summary>
        /// <param name="seed">the seed value</param>
        void setSeed(int[] seed);

        /// <summary>
        /// Sets the seed of the underlying random number generator using a
        /// <c>long</c> seed.
        /// <para>Sequences of values generated starting with the same seeds
        /// should be identical.
        /// </para>
        /// </summary>
        /// <param name="seed">the seed value</param>
        void setSeed(long seed);

        /// <summary>
        /// Generates random bytes and places them into a user-supplied
        /// byte array.  The number of random bytes produced is equal to
        /// the length of the byte array.
        /// </summary>
        /// <param name="bytes">bytes the non-null byte array in which to put the
        /// random bytes</param>
        void nextBytes(byte[] bytes);

        /// <summary>
        /// Returns the next pseudorandom, uniformly distributed <c>int</c>
        /// value from this random number generator's sequence.
        /// All 2<^32 possible <c>int</c> values
        /// should be produced with  (approximately) equal probability.
        /// </summary>
        /// <returns>the next pseudorandom, uniformly distributed <c>int</c>
        /// value from this random number generator's sequence</returns>
        int nextInt();

        /// <summary>
        /// Returns a pseudorandom, uniformly distributed <c>int</c> value
        /// between 0 (inclusive) and the specified value (exclusive), drawn from
        /// this random number generator's sequence.
        /// </summary>
        /// <param name="n">the bound on the random number to be returned.  Must be
        /// positive.</param>
        /// <returns>a pseudorandom, uniformly distributed <c>int</c>
        /// value between 0 (inclusive) and n (exclusive).</returns>
        /// <exception cref="IllegalArgumentException"> if n is not positive.</exception>
        int nextInt(int n);

        /// <summary>
        /// Returns the next pseudorandom, uniformly distributed <c>long</c>
        /// value from this random number generator's sequence.  All
        /// 2^64 possible <c>long</c> values
        /// should be produced with (approximately) equal probability.
        /// </summary>
        /// <returns>the next pseudorandom, uniformly distributed <c>long</c>
        /// value from this random number generator's sequence</returns>
        long nextLong();

        /// <summary>
        /// Returns the next pseudorandom, uniformly distributed
        /// <c>boolean</c> value from this random number generator's
        /// sequence.
        /// </summary>
        /// <returns>the next pseudorandom, uniformly distributed
        /// <c>boolean</c> value from this random number generator's
        /// sequence</returns>
        Boolean nextBoolean();

        /// <summary>
        /// Returns the next pseudorandom, uniformly distributed <c>float</c>
        /// value between <c>0.0</c> and <c>1.0</c> from this random
        /// number generator's sequence. 
        /// </summary>
        /// <returns>the next pseudorandom, uniformly distributed <c>float</c>
        /// value between <c>0.0</c> and <c>1.0</c> from this
        /// random number generator's sequence</returns>
        float nextFloat();

        /// <summary>
        /// Returns the next pseudorandom, uniformly distributed
        /// <c>double</c> value between <c>0.0</c> and
        /// <c>1.0</c> from this random number generator's sequence.
        /// </summary>
        /// <returns>the next pseudorandom, uniformly distributed
        /// <c>double</c> value between <c>0.0</c> and
        /// <c>1.0</c> from this random number generator's sequence</returns>
        double nextDouble();

        /// <summary>
        /// Returns the next pseudorandom, Gaussian ("normally") distributed
        /// <c>double</c> value with mean <c>0.0</c> and standard
        /// deviation <c>1.0</c> from this random number generator's sequence. 
        /// </summary>
        /// <returns>the next pseudorandom, Gaussian ("normally") distributed
        /// <c>double</c> value with mean <c>0.0</c> and
        /// standard deviation <c>1.0</c> from this random number
        /// generator's sequence</returns>
        double nextGaussian();
    }
}
