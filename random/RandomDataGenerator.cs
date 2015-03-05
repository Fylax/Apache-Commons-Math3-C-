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
using Math3.distribution;
using Math3.exception;
using Math3.exception.util;
using Math3.util;
using Org.BouncyCastle.Crypto.Digests;
using Org.BouncyCastle.Security;
using System;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;

namespace Math3.random
{
    /// <summary>
    /// Implements the <see cref="RandomData"/> interface using a <see cref="RandomGenerator"/>
    /// instance to generate non-secure data and a SecureRandom
    /// instance to provide data for the <c>nextSecureXxx</c> methods. If no
    /// <c>RandomGenerator</c> is provided in the constructor, the default is
    /// to use a <see cref="Well19937c"/> generator. To plug in a different
    /// implementation, either implement <c>RandomGenerator</c> directly or
    /// extend <see cref="AbstractRandomGenerator"/>.
    /// <para>
    /// Supports reseeding the underlying pseudo-random number generator (PRNG). The
    /// <c>SecurityProvider</c> and <c>Algorithm</c> used by the
    /// <c>SecureRandom</c> instance can also be reset.
    /// </para>
    /// <para>
    /// For details on the default PRNGs, see <see cref="System.Math.Random"/> and SecureRandom".
    /// </para>
    /// <para>
    /// <strong>Usage Notes</strong>:
    /// <list type="bullet">
    /// <item>
    /// Instance variables are used to maintain <c>RandomGenerator</c> and
    /// <c>SecureRandom</c> instances used in data generation. Therefore, to
    /// generate a random sequence of values or strings, you should use just
    /// one <c>RandomDataImpl</c> instance repeatedly.</item>
    /// <item>
    /// The "secure" methods are *much* slower. These should be used only when a
    /// cryptographically secure random sequence is required. A secure random
    /// sequence is a sequence of pseudo-random values which, in addition to being
    /// well-dispersed (so no subsequence of values is an any more likely than other
    /// subsequence of the the same length), also has the additional property that
    /// knowledge of values generated up to any point in the sequence does not make
    /// it any easier to predict subsequent values.</item>
    /// <item>
    /// When a new <c>RandomDataImpl</c> is created, the underlying random
    /// number generators are not initialized. If you do not
    /// explicitly seed the default non-secure generator, it is seeded with the
    /// current time in milliseconds plus the system identity hash code on first use.
    /// The same holds for the secure generator. If you provide a <code>RandomGenerator</code>
    /// to the constructor, however, this generator is not reseeded by the constructor
    /// nor is it reseeded on first use.</item>
    /// <item>
    /// The <c>reSeed</c> and <c>reSeedSecure</c> methods delegate to the
    /// corresponding methods on the underlying <c>RandomGenerator</c> and
    /// <c>SecureRandom</c> instances. Therefore, <c>reSeed(long)</c>
    /// fully resets the initial state of the non-secure random number generator (so
    /// that reseeding with a specific value always results in the same subsequent
    /// random sequence); whereas reSeedSecure(long) does not
    /// reinitialize the secure random number generator (so secure sequences started
    /// with calls to reseedSecure(long) won't be identical).</item>
    /// <item>
    /// This implementation is not synchronized. The underlying <code>RandomGenerator</code>
    /// or <c>SecureRandom</c> instances are not protected by synchronization and
    /// are not guaranteed to be thread-safe.  Therefore, if an instance of this class
    /// is concurrently utilized by multiple threads, it is the responsibility of
    /// client code to synchronize access to seeding and data generation methods.
    /// </item>
    /// </list>
    /// </para>
    /// </summary>
#pragma warning disable 0612
    public class RandomDataGenerator : RandomData
    {
        /// <summary>
        /// underlying random number generator
        /// </summary>
        private RandomGenerator rand = null;

        /// <summary>
        /// underlying secure random number generator
        /// </summary>
        private RandomGenerator secRand = null;

        /// <summary>
        /// Construct a RandomDataGenerator, using a default random generator as the source
        /// of randomness.
        /// <para>The default generator is a <see cref="Well19937c"/> seeded
        /// with <c>System.currentTimeMillis() + System.identityHashCode(this))</c>.
        /// The generator is initialized and seeded on first use.</para>
        /// </summary>
        public RandomDataGenerator() { }

        /// <summary>
        /// Construct a RandomDataGenerator using the supplied <see cref="RandomGenerator"/> as
        /// the source of (non-secure) random data.
        /// </summary>
        /// <param name="rand">the source of (non-secure) random data
        /// (may be null, resulting in the default generator)</param>
        public RandomDataGenerator(RandomGenerator rand)
        {
            this.rand = rand;
        }

        /// <inheritdoc/>
        /// <param name="len">the desired string length.</param>
        /// <returns>the random string.</returns>
        /// <remarks>
        /// <para>
        /// Algorithm Description: hex strings are generated using a
        /// 2-step process.
        /// <list type="numbers">
        /// <item><c>len / 2 + 1</c> binary bytes are generated using the underlying
        /// Random</item>
        /// <item>Each binary byte is translated into 2 hex digits</item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <exception cref="NotStrictlyPositiveException"> if <c>len <= 0</c>.</exception>
        public String nextHexString(int len)
        {
            if (len <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("LENGTH"), len);
            }

            // Get a random number generator
            RandomGenerator ran = getRandomGenerator();

            // Initialize output buffer
            StringBuilder outBuffer = new StringBuilder();

            // Get int(len/2)+1 random bytes
            byte[] randomBytes = new byte[(len / 2) + 1];
            ran.nextBytes(randomBytes);

            // Convert each byte to 2 hex digits
            for (int i = 0; i < randomBytes.Length; i++)
            {
                Int32 c = randomBytes[i];

                /*
               /// Add 128 to byte value to make interval 0-255 before doing hex
               /// conversion. This guarantees <= 2 hex digits from toHexString()
               /// toHexString would otherwise add 2^32 to negative arguments.
                 */
                String hex = (c + 128).ToString("X");

                // Make sure we add 2 hex digits for each byte
                if (hex.Length == 1)
                {
                    hex = "0" + hex;
                }
                outBuffer.Append(hex);
            }
            return outBuffer.ToString().Substring(0, len);
        }

        ///<inheritdoc/>
        public int nextInt(int lower, int upper)
        {
            return new UniformIntegerDistribution(getRandomGenerator(), lower, upper).sample();
        }

        /// <inheritdoc/>
        public long nextLong(long lower, long upper)
        {
            if (lower >= upper)
            {
                throw new NumberIsTooLargeException<Int64, Int64>(new LocalizedFormats("LOWER_BOUND_NOT_BELOW_UPPER_BOUND"), lower, upper, false);
            }
            long max = (upper - lower) + 1;
            if (max <= 0)
            {
                // the range is too wide to fit in a positive long (larger than 2^63); as it covers
                // more than half the long range, we use directly a simple rejection method
                RandomGenerator rng = getRandomGenerator();
                while (true)
                {
                    long r = rng.nextLong();
                    if (r >= lower && r <= upper)
                    {
                        return r;
                    }
                }
            }
            else if (max < Int32.MaxValue)
            {
                // we can shift the range and generate directly a positive int
                return lower + getRandomGenerator().nextInt((int)max);
            }
            else
            {
                // we can shift the range and generate directly a positive long
                return lower + nextLong(getRandomGenerator(), max);
            }
        }

        /// <param name="rng">random generator to use</param>
        /// <param name="n">the bound on the random number to be returned.  Must be
        /// positive.</param>
        /// <returns>a pseudorandom, uniformly distributed <c>long</c>
        /// value between 0 (inclusive) and n (exclusive).</returns>
        /// <remarks>
        /// Returns a pseudorandom, uniformly distributed <c>long</c> value
        /// between 0 (inclusive) and the specified value (exclusive), drawn from
        /// this random number generator's sequence.
        /// </remarks>
        /// <exception cref="NotStrictlyPositiveException">if n is not positive.</exception>
        private static long nextLong(RandomGenerator rng, long n)
        {
            if (n > 0)
            {
                byte[] byteArray = new byte[8];
                long bits;
                long val;
                do
                {
                    rng.nextBytes(byteArray);
                    bits = 0;
                    foreach (byte b in byteArray)
                    {
                        bits = (bits << 8) | (((long)b) & 0xffL);
                    }
                    bits &= 0x7fffffffffffffffL;
                    val = bits % n;
                } while (bits - val + (n - 1) < 0);
                return val;
            }
            throw new NotStrictlyPositiveException<Int64>(n);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// <strong>Algorithm Description:</strong> hex strings are generated in
        /// 40-byte segments using a 3-step process.
        /// <list type="numbers">
        /// <item>
        /// 20 random bytes are generated using the underlying
        /// <code>SecureRandom</code>.</item>
        /// <item>
        /// SHA-1 hash is applied to yield a 20-byte binary digest.</item>
        /// <item>
        /// Each byte of the binary digest is converted to 2 hex digits.</item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <exception cref="NotStrictlyPositiveException"> if <c>len <= 0</c></exception>
        public String nextSecureHexString(int len)
        {
            if (len <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("LENGTH"), len);
            }

            // Get SecureRandom and setup Digest provider
            RandomGenerator secRan = getSecRan();
            Sha1Digest alg = new Sha1Digest(); //The buildin HashAlgorithm is not supported as portable, so I used the BouncyCastle, from NuGet
            alg.Reset();

            // Compute number of iterations required (40 bytes each)
            int numIter = (len / 40) + 1;

            StringBuilder outBuffer = new StringBuilder();
            for (int iter = 1; iter < numIter + 1; iter++)
            {
                byte[] randomBytes = new byte[40];
                secRan.nextBytes(randomBytes);
                alg.BlockUpdate(randomBytes, 0, randomBytes.Length - 1);

                // Compute hash -- will create 20-byte binary hash
                byte[] hash = new byte[20];
                alg.DoFinal(hash, 0);

                // Loop over the hash, converting each byte to 2 hex digits
                for (int i = 0; i < hash.Length; i++)
                {
                    Int32 c = hash[i];

                    /*
                   /// Add 128 to byte value to make interval 0-255 This guarantees
                   /// <= 2 hex digits from toHexString() toHexString would
                   /// otherwise add 2^32 to negative arguments
                     */
                    String hex = (c + 128).ToString("X");

                    // Keep strings uniform length -- guarantees 40 bytes
                    if (hex.Length == 1)
                    {
                        hex = "0" + hex;
                    }
                    outBuffer.Append(hex);
                }
            }
            return outBuffer.ToString().Substring(0, len);
        }

        /// <inheritdoc/>
        public int nextSecureInt(int lower, int upper)
        {
            return new UniformIntegerDistribution(getSecRan(), lower, upper).sample();
        }

        /// <inheritdoc/>
        public long nextSecureLong(long lower, long upper)
        {
            if (lower >= upper)
            {
                throw new NumberIsTooLargeException<Int64, Int64>(new LocalizedFormats("LOWER_BOUND_NOT_BELOW_UPPER_BOUND"), lower, upper, false);
            }
            RandomGenerator rng = getSecRan();
            long max = (upper - lower) + 1;
            if (max <= 0)
            {
                // the range is too wide to fit in a positive long (larger than 2^63); as it covers
                // more than half the long range, we use directly a simple rejection method
                while (true)
                {
                    long r = rng.nextLong();
                    if (r >= lower && r <= upper)
                    {
                        return r;
                    }
                }
            }
            else if (max < Int32.MaxValue)
            {
                // we can shift the range and generate directly a positive int
                return lower + rng.nextInt((int)max);
            }
            else
            {
                // we can shift the range and generate directly a positive long
                return lower + nextLong(rng, max);
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Algorithm Description:
        /// <list type="bullet"><item> For small means, uses simulation of a Poisson process
        /// using Uniform deviates, as described
        /// <a href="http://irmi.epfl.ch/cmos/Pmmi/interactive/rng7.htm"> here.</a>
        /// The Poisson process (and hence value returned) is bounded by 1000 * mean.</item>
        /// <item> For large means, uses the rejection algorithm described in <para/>
        /// Devroye, Luc. (1981).The Computer Generation of Poisson Random Variables
        /// <strong>Computing</strong> vol. 26 pp. 197-207.</item></list></para>
        /// </remarks>
        /// <exception cref="NotStrictlyPositiveException"> if <c>len <= 0</c>.</exception>
        public long nextPoisson(double mean)
        {
            return new PoissonDistribution(getRandomGenerator(), mean,
                    PoissonDistribution.DEFAULT_EPSILON,
                    PoissonDistribution.DEFAULT_MAX_ITERATIONS).sample();
        }

        /// <inheritdoc/>
        public double nextGaussian(double mu, double sigma)
        {
            if (sigma <= 0)
            {
                throw new NotStrictlyPositiveException<Double>(new LocalizedFormats("STANDARD_DEVIATION"), sigma);
            }
            return sigma * getRandomGenerator().nextGaussian() + mu;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Algorithm Description: Uses the Algorithm SA (Ahrens)
        /// from p. 876 in:
        /// [1]: Ahrens, J. H. and Dieter, U. (1972). Computer methods for
        /// sampling from the exponential and normal distributions.
        /// Communications of the ACM, 15, 873-882.
        /// </para>
        /// </remarks>
        public double nextExponential(double mean)
        {
            return new ExponentialDistribution(getRandomGenerator(), mean,
                    ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// <para>Generates a random value from the
        /// <see cref="org.apache.commons.math3.distribution.GammaDistribution Gamma Distribution"/>.</para>
        /// <para>This implementation uses the following algorithms: </para>
        /// <para>For 0 < shape < 1: <br/>
        /// Ahrens, J. H. and Dieter, U., <i>Computer methods for
        /// sampling from gamma, beta, Poisson and binomial distributions.</i>
        /// Computing, 12, 223-246, 1974.</para>
        /// <para>For shape >= 1: <br/>
        /// Marsaglia and Tsang, <i>A Simple Method for Generating
        /// Gamma Variables.</i> ACM Transactions on Mathematical Software,
        /// Volume 26 Issue 3, September, 2000.</para> 
        /// </summary>
        /// <param name="shape">the median of the Gamma distribution</param>
        /// <param name="scale">the scale parameter of the Gamma distribution</param>
        /// <returns>random value sampled from the Gamma(shape, scale) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException">if <c>shape <= 0</c> or
        /// <c>scale <= 0</c>.</exception>
        public double nextGamma(double shape, double scale)
        {
            return new GammaDistribution(getRandomGenerator(), shape, scale,
                    GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="HypergeometricDistribution Hypergeometric Distribution"/>.  
        /// </summary>
        /// <param name="populationSize">the population size of the Hypergeometric 
        /// distribution</param>
        /// <param name="numberOfSuccesses">number of successes in the population of 
        /// the Hypergeometric distribution</param>
        /// <param name="sampleSize">he sample size of the Hypergeometric distribution</param>
        /// <returns>random value sampled from the Hypergeometric(numberOfSuccesses, sampleSize) 
        /// distribution</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>numberOfSuccesses > populationSize</c>,
        /// or <c>sampleSize > populationSize</c>.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if <c>populationSize <= 0</c>.</exception>
        /// <exception cref="NotPositiveException"> if <c>numberOfSuccesses < 0</c></exception>
        public int nextHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize)
        {
            return new HypergeometricDistribution(getRandomGenerator(), populationSize,
                    numberOfSuccesses, sampleSize).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="PascalDistribution Pascal Distribution"/>.
        /// </summary>
        /// <param name="r">the number of successes of the Pascal distribution</param>
        /// <param name="p">the probability of success of the Pascal distribution</param>
        /// <returns>random value sampled from the Pascal(r, p) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if the number of successes is not positive</exception>
        /// <exception cref="OutOfRangeException"> if the probability of success is not in the
        /// range <c>[0, 1]</c>.</exception>
        public int nextPascal(int r, double p)
        {
            return new PascalDistribution(getRandomGenerator(), r, p).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="TDistribution T Distribution"/>. 
        /// </summary>
        /// <param name="df">the degrees of freedom of the T distribution</param>
        /// <returns>random value from the T(df) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>df <= 0</c></exception>
        public double nextT(double df)
        {
            return new TDistribution(getRandomGenerator(), df,
                    TDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="WeibullDistribution Weibull Distribution"/>. 
        /// </summary>
        /// <param name="shape">the shape parameter of the Weibull distribution</param>
        /// <param name="scale">the scale parameter of the Weibull distribution</param>
        /// <returns>value sampled from the Weibull(shape, size) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException">if <c>shape <= 0</c> or
        /// <c>scale <= 0</c>.</exception>
        public double nextWeibull(double shape, double scale)
        {
            return new WeibullDistribution(getRandomGenerator(), shape, scale,
                    WeibullDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="ZipfDistribution Zipf Distribution"/>. 
        /// </summary>
        /// <param name="numberOfElements">the number of elements of the ZipfDistribution</param>
        /// <param name="exponent">the exponent of the ZipfDistribution</param>
        /// <returns>random value sampled from the Zipf(numberOfElements, exponent) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>numberOfElements <= 0</c>
        /// or <c>exponent <= 0</c>.</exception>
        public int nextZipf(int numberOfElements, double exponent)
        {
            return new ZipfDistribution(getRandomGenerator(), numberOfElements, exponent).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="BetaDistribution Beta Distribution"/>.
        /// </summary>
        /// <param name="alpha">first distribution shape parameter</param>
        /// <param name="beta">second distribution shape parameter</param>
        /// <returns>random value sampled from the beta(alpha, beta) distribution</returns>
        public double nextBeta(double alpha, double beta)
        {
            return new BetaDistribution(getRandomGenerator(), alpha, beta,
                    BetaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="BinomialDistribution Binomial Distribution"/>.
        /// </summary>
        /// <param name="numberOfTrials">number of trials of the Binomial distribution</param>
        /// <param name="probabilityOfSuccess">probability of success of the Binomial distribution</param>
        /// <returns>random value sampled from the Binomial(numberOfTrials, probabilityOfSuccess) distribution</returns>
        public int nextBinomial(int numberOfTrials, double probabilityOfSuccess)
        {
            return new BinomialDistribution(getRandomGenerator(), numberOfTrials, probabilityOfSuccess).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="CauchyDistribution Cauchy Distribution"/>.
        /// </summary>
        /// <param name="median">the median of the Cauchy distribution</param>
        /// <param name="scale">the scale parameter of the Cauchy distribution</param>
        /// <returns>random value sampled from the Cauchy(median, scale) distribution</returns>
        public double nextCauchy(double median, double scale)
        {
            return new CauchyDistribution(getRandomGenerator(), median, scale,
                    CauchyDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="ChiSquaredDistribution ChiSquare Distribution"/>.
        /// </summary>
        /// <param name="df">the degrees of freedom of the ChiSquare distribution</param>
        /// <returns>random value sampled from the ChiSquare(df) distribution</returns>
        public double nextChiSquare(double df)
        {
            return new ChiSquaredDistribution(getRandomGenerator(), df,
                    ChiSquaredDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <summary>
        /// Generates a random value from the <see cref="FDistribution F Distribution"/>.  
        /// </summary>
        /// <param name="numeratorDf">the numerator degrees of freedom of the F 
        /// distribution</param>
        /// <param name="denominatorDf">the denominator degrees of freedom of the F
        /// distribution</param>
        /// <returns>random value sampled from the F(numeratorDf, denominatorDf) 
        /// distribution</returns>
        /// <exception cref="NotStrictlyPositiveException">if
        /// <c>numeratorDf <= 0</c> or <c>denominatorDf <= 0</c>.</exception>
        public double nextF(double numeratorDf, double denominatorDf)
        {
            return new FDistribution(getRandomGenerator(), numeratorDf, denominatorDf,
                    FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY).sample();
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Algorithm Description: scales the output of
        /// Random.NextDouble(), but rejects 0 values (i.e., will generate another
        /// random double if Random.NextDouble() returns 0). This is necessary to
        /// provide a symmetric output interval (both endpoints excluded).
        /// </para>
        /// </remarks>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c></exception>
        /// <exception cref="NotFiniteNumberException"> if one of the bounds is infinite
        /// </exception>
        /// <exception cref="NotANumberException"> if one of the bounds is NaN</exception>
        public double nextUniform(double lower, double upper)
        {
            return nextUniform(lower, upper, false);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Algorithm Description: if the lower bound is excluded,
        /// scales the output of Random.nextDouble(), but rejects 0 values (i.e.,
        /// will generate another random double if Random.nextDouble() returns 0).
        /// This is necessary to provide a symmetric output interval (both
        /// endpoints excluded).
        /// </para>
        /// </remarks>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c></exception>
        /// <exception cref="NotFiniteNumberException"> if one of the bounds is infinite
        /// </exception>
        /// <exception cref="NotANumberException"> if one of the bounds is NaN</exception>
        public double nextUniform(double lower, double upper, Boolean lowerInclusive)
        {
            if (lower >= upper)
            {
                throw new NumberIsTooLargeException<Double, Double>(new LocalizedFormats("LOWER_BOUND_NOT_BELOW_UPPER_BOUND"), lower, upper, false);
            }

            if (Double.IsInfinity(lower))
            {
                throw new NotFiniteNumberException<Double>(new LocalizedFormats("INFINITE_BOUND"), lower);
            }
            if (Double.IsInfinity(upper))
            {
                throw new NotFiniteNumberException<Double>(new LocalizedFormats("INFINITE_BOUND"), upper);
            }

            if (Double.IsNaN(lower) || Double.IsNaN(upper))
            {
                throw new NotANumberException();
            }

            RandomGenerator generator = getRandomGenerator();

            // ensure nextDouble() isn't 0.0
            double u = generator.nextDouble();
            while (!lowerInclusive && u <= 0.0)
            {
                u = generator.nextDouble();
            }

            return u * upper + (1.0 - u) * lower;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This method calls <see cref="MathArrays.shuffle(int[],RandomGenerator)"/> 
        /// in order to create a random shuffle of the set
        /// of natural numbers <c>{ 0, 1, ..., n - 1 </c>}.
        /// </remarks>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if <c>k <= 0</c>.</exception>
        public int[] nextPermutation(int n, int k)
        {
            if (k > n)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("PERMUTATION_EXCEEDS_N"), k, n, true);
            }
            if (k <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("PERMUTATION_SIZE"), k);
            }

            int[] index = MathArrays.natural(n);
            MathArrays.shuffle(index, getRandomGenerator());

            // Return a new array containing the first "k" entries of "index".
            return MathArrays.copyOf(index, k);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This method calls <see cref="nextPermutation(int,int)">nextPermutation(c.size(), k)
        /// </see>
        /// </remarks>
        /// in order to sample the collection.
        public Object[] nextSample<T>(Collection<T> c, int k)
        {

            int len = c.Count;
            if (k > len)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE"), k, len, true);
            }
            if (k <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_SAMPLES"), k);
            }

            T[] objects = c.ToArray<T>();
            int[] index = nextPermutation(len, k);
            Object[] result = new Object[k];
            for (int i = 0; i < k; i++)
            {
                result[i] = objects[index[i]];
            }
            return result;
        }



        /// <summary>
        /// Reseeds the random number generator with the supplied seed.
        /// <para>
        /// Will create and initialize if null.
        /// </para>
        /// </summary>
        /// <param name="seed">seed the seed value to use</param>
        public void reSeed(long seed)
        {
            getRandomGenerator().setSeed(seed);
        }

        /// <summary>
        /// Reseeds the secure random number generator with the current time in
        /// milliseconds.
        /// <para>
        /// Will create and initialize if null.
        /// </para>
        /// </summary>
        public void reSeedSecure()
        {
            //Find unix timestamp (seconds since 01/01/1970)
            long ticks = DateTime.UtcNow.Ticks - DateTime.Parse("01/01/1970 00:00:00").Ticks;
            ticks /= 10000; //current timestamp in millis
            getSecRan().setSeed(ticks);
        }

        /// <summary>
        /// Reseeds the secure random number generator with the supplied seed.
        /// <para>
        /// Will create and initialize if null.
        /// </para>
        /// </summary>
        /// <param name="seed">the seed value to use</param>
        public void reSeedSecure(long seed)
        {
            getSecRan().setSeed(seed);
        }

        /// <summary>
        /// Reseeds the random number generator with
        /// <c>System.currentTimeMillis() + System.identityHashCode(this))</c>.
        /// </summary>
        public void reSeed()
        {
            //Find unix timestamp (seconds since 01/01/1970)
            long ticks = DateTime.UtcNow.Ticks - DateTime.Parse("01/01/1970 00:00:00").Ticks;
            ticks /= 10000; //current timestamp in millis
            getRandomGenerator().setSeed(ticks + this.GetHashCode());
        }

        /// <summary>
        /// Sets the PRNG algorithm for the underlying SecureRandom instance using
        /// the Security Provider API. The Security Provider API is defined in <a
        /// href =
        /// "http://java.sun.com/j2se/1.3/docs/guide/security/CryptoSpec.html#AppA">
        /// Java Cryptography Architecture API Specification & Reference.</a>
        /// <para>
        /// <strong>USAGE NOTE:</strong> This method carries <i>significant</i>
        /// overhead and may take several seconds to execute.
        /// </para>
        /// </summary>
        /// <param name="algorithm">the name of the PRNG algorithm</param>
        public void setSecureAlgorithm(String algorithm)
        {
            secRand = RandomGeneratorFactory.createRandomGenerator(SecureRandom.GetInstance(algorithm));
        }

        /// <summary>
        /// Returns the RandomGenerator used to generate non-secure random data.
        /// <para>
        /// Creates and initializes a default generator if null. Uses a <see cref="Well19937c"/>
        /// generator with <c>System.currentTimeMillis() + System.identityHashCode(this))</c>
        /// as the default seed.
        /// </para>
        /// </summary>
        /// <returns>the Random used to generate random data</returns>
        public RandomGenerator getRandomGenerator()
        {
            if (rand == null)
            {
                initRan();
            }
            return rand;
        }

        /// <summary>
        /// Sets the default generator to a <see cref="Well19937c"/> generator seeded with
        /// <c>System.currentTimeMillis() + System.identityHashCode(this))</c>.
        /// </summary>
        private void initRan()
        {
            //Find unix timestamp (seconds since 01/01/1970)
            long ticks = DateTime.UtcNow.Ticks - DateTime.Parse("01/01/1970 00:00:00").Ticks;
            ticks /= 10000; //current timestamp in millis
            rand = new Well19937c(ticks + this.GetHashCode());
        }

        /// <summary>
        /// Returns the SecureRandom used to generate secure random data.
        /// <para>
        /// Creates and initializes if null.  Uses
        /// <c>System.currentTimeMillis() + System.identityHashCode(this)</c> as the default seed.
        /// </para>
        /// </summary>
        /// <returns>the SecureRandom used to generate secure random data, wrapped in a
        /// <see cref="RandomGenerator"/>.</returns>
        private RandomGenerator getSecRan()
        {
            if (secRand == null)
            {
                secRand = RandomGeneratorFactory.createRandomGenerator(new SecureRandom());
                //Find unix timestamp (seconds since 01/01/1970)
                long ticks = DateTime.UtcNow.Ticks - DateTime.Parse("01/01/1970 00:00:00").Ticks;
                ticks /= 10000; //current timestamp in millis
                secRand.setSeed(ticks + this.GetHashCode());
            }
            return secRand;
        }
    }
}
