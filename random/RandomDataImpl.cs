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
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Math3.random
{
    /// <summary>
    /// Generates random deviates and other random data using a <see cref="RandomGenerator"/>
    /// instance to generate non-secure data and a <see cref="java.security.SecureRandom"/>
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
    /// For details on the default PRNGs, see <see cref="java.util.Random"/> and
    /// <see cref="java.security.SecureRandom"/>.
    /// </para>
    /// <para>
    /// <strong>Usage Notes</strong>:
    /// <list type="bullet">
    /// <item>
    /// Instance variables are used to maintain <c>RandomGenerator</c> and
    /// <c>SecureRandom</c> instances used in data generation. Therefore, to
    /// generate a random sequence of values or strings, you should use just
    /// one <c>RandomDataGenerator</c> instance repeatedly.</item>
    /// <item>
    /// The "secure" methods are *much* slower. These should be used only when a
    /// cryptographically secure random sequence is required. A secure random
    /// sequence is a sequence of pseudo-random values which, in addition to being
    /// well-dispersed (so no subsequence of values is an any more likely than other
    /// subsequence of the the same length), also has the additional property that
    /// knowledge of values generated up to any point in the sequence does not make
    /// it any easier to predict subsequent values.</item>
    /// <item>
    /// When a new <c>RandomDataGenerator</c> is created, the underlying random
    /// number generators are <strong>not</strong> initialized. If you do not
    /// explicitly seed the default non-secure generator, it is seeded with the
    /// current time in milliseconds plus the system identity hash code on first use.
    /// The same holds for the secure generator. If you provide a <c>RandomGenerator</c>
    /// to the constructor, however, this generator is not reseeded by the constructor
    /// nor is it reseeded on first use.</item>
    /// <item>
    /// The <c>reSeed</c> and <c>reSeedSecure</c> methods delegate to the
    /// corresponding methods on the underlying <c>RandomGenerator</c> and
    /// <c>SecureRandom</c> instances. Therefore, <c>reSeed(long)</c>
    /// fully resets the initial state of the non-secure random number generator (so
    /// that reseeding with a specific value always results in the same subsequent
    /// random sequence); whereas reSeedSecure(long) does <strong>not</strong>
    /// reinitialize the secure random number generator (so secure sequences started
    /// with calls to reseedSecure(long) won't be identical).</item>
    /// <item>
    /// This implementation is not synchronized. The underlying <c>RandomGenerator</c>
    /// or <c>SecureRandom</c> instances are not protected by synchronization and
    /// are not guaranteed to be thread-safe.  Therefore, if an instance of this class
    /// is concurrently utilized by multiple threads, it is the responsibility of
    /// client code to synchronize access to seeding and data generation methods.
    /// </item>
    /// </list>
    /// </para>
    /// </summary>
    [Obsolete]
    public class RandomDataImpl : RandomData
    {
        /// <summary>
        /// RandomDataGenerator delegate
        /// </summary>
        private readonly RandomDataGenerator delegates;

        /// <summary>
        /// Construct a RandomDataImpl, using a default random generator as the source
        /// of randomness.
        /// <para>The default generator is a <see cref="Well19937c"/> seeded
        /// with <c>System.currentTimeMillis() + System.identityHashCode(this))</c>.
        /// The generator is initialized and seeded on first use.</para>
        /// </summary>
        public RandomDataImpl()
        {
            delegates = new RandomDataGenerator();
        }

        /// <summary>
        /// Construct a RandomDataImpl using the supplied <see cref="RandomGenerator"/> as
        /// the source of (non-secure) random data.
        /// </summary>
        /// <param name="rand">rand the source of (non-secure) random data
        /// (may be null, resulting in the default generator)</param>
        public RandomDataImpl(RandomGenerator rand)
        {
            delegates = new RandomDataGenerator(rand);
        }

        /// <returns>the delegate object.</returns>
        [Obsolete]
        RandomDataGenerator getDelegate()
        {
            return delegates;
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
            return delegates.nextHexString(len);
        }

        /// <inheritdoc/>
        public int nextInt(int lower, int upper)
        {
            return delegates.nextInt(lower, upper);
        }

        /// <inheritdoc/>
        public long nextLong(long lower, long upper)
        {
            return delegates.nextLong(lower, upper);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Algorithm Description: hex strings are generated in
        /// 40-byte segments using a 3-step process.
        /// <list type="numbers">
        /// <item>
        /// 20 random bytes are generated using the underlying
        /// <c>SecureRandom</c>.</item>
        /// <item>
        /// SHA-1 hash is applied to yield a 20-byte binary digest.</item>
        /// <item>
        /// Each byte of the binary digest is converted to 2 hex digits.</item>
        /// </list>
        /// </para>
        /// </remarks>
        public String nextSecureHexString(int len)
        {
            return delegates.nextSecureHexString(len);
        }

        /// <inheritdoc/>
        public int nextSecureInt(int lower, int upper)
        {
            return delegates.nextSecureInt(lower, upper);
        }

        /// <inheritdoc/>
        public long nextSecureLong(long lower, long upper)
        {
            return delegates.nextSecureLong(lower, upper);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// <strong>Algorithm Description</strong>:
        /// <list type="bullet"><item> For small means, uses simulation of a Poisson process
        /// using Uniform deviates, as described
        /// <a href="http://irmi.epfl.ch/cmos/Pmmi/interactive/rng7.htm"> here.</a>
        /// The Poisson process (and hence value returned) is bounded by 1000 * mean.</item>
        /// <item> For large means, uses the rejection algorithm described in <para/>
        /// Devroye, Luc. (1981).<i>The Computer Generation of Poisson Random Variables
        /// Computing vol. 26 pp. 197-207.</item></list></para>
        /// </remarks>
        public long nextPoisson(double mean)
        {
            return delegates.nextPoisson(mean);
        }

        /// <inheritdoc/>
        public double nextGaussian(double mu, double sigma)
        {
            return delegates.nextGaussian(mu, sigma);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// <strong>Algorithm Description</strong>: Uses the Algorithm SA (Ahrens)
        /// from p. 876 in:
        /// [1]: Ahrens, J. H. and Dieter, U. (1972). Computer methods for
        /// sampling from the exponential and normal distributions.
        /// Communications of the ACM, 15, 873-882.
        /// </para>
        /// </remarks>
        public double nextExponential(double mean)
        {
            return delegates.nextExponential(mean);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Algorithm Description: scales the output of
        /// Random.nextDouble(), but rejects 0 values (i.e., will generate another
        /// random double if Random.nextDouble() returns 0). This is necessary to
        /// provide a symmetric output interval (both endpoints excluded).
        /// </para>
        /// </remarks>
        public double nextUniform(double lower, double upper)
        {
            return delegates.nextUniform(lower, upper);
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
        public double nextUniform(double lower, double upper, Boolean lowerInclusive)
        {
            return delegates.nextUniform(lower, upper, lowerInclusive);
        }

        /// <summary>
        /// Generates a random value from the <see cref="BetaDistribution">Beta Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(RealDistribution)">inversion</see>
        /// to generate random values.
        /// </summary>
        /// <param name="alpha">first distribution shape parameter</param>
        /// <param name="beta">second distribution shape parameter</param>
        /// <returns>random value sampled from the beta(alpha, beta) distribution</returns>
        public double nextBeta(double alpha, double beta)
        {
            return delegates.nextBeta(alpha, beta);
        }

        /// <summary>
        /// Generates a random value from the <see cref="BinomialDistribution">Binomial 
        /// Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(RealDistribution)">inversion
        /// </see> to generate random values.
        /// </summary>
        /// <param name="numberOfTrials">number of trials of the Binomial distribution</param>
        /// <param name="probabilityOfSuccess">probability of success of the Binomial 
        /// distribution</param>
        /// <returns>random value sampled from the Binomial(numberOfTrials, probabilityOfSuccess)
        /// distribution</returns>
        public int nextBinomial(int numberOfTrials, double probabilityOfSuccess)
        {
            return delegates.nextBinomial(numberOfTrials, probabilityOfSuccess);
        }

        /// <summary>
        /// Generates a random value from the <see cref="CauchyDistribution">Cauchy 
        /// Distribution<see/>.
        /// This implementation uses <see cref="nextInversionDeviate(RealDistribution)">inversion
        /// </see> to generate random values.
        /// </summary>
        /// <param name="median">the median of the Cauchy distribution</param>
        /// <param name="scale">the scale parameter of the Cauchy distribution</param>
        /// <returns>random value sampled from the Cauchy(median, scale) distribution</returns>
        public double nextCauchy(double median, double scale)
        {
            return delegates.nextCauchy(median, scale);
        }

        /// <summary>
        /// Generates a random value from the <see cref="ChiSquaredDistribution">ChiSquare 
        /// Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(RealDistribution)"> inversion"</see> to generate random values.
        /// </summary>
        /// <param name="df">the degrees of freedom of the ChiSquare distribution</param>
        /// <returns>random value sampled from the ChiSquare(df) distribution</returns>
        public double nextChiSquare(double df)
        {
            return delegates.nextChiSquare(df);
        }

        /// <summary>
        /// Generates a random value from the <see cref="FDistribution">F-Distribution</see>.
        ///This implementation uses <see cref="nextInversionDeviate(RealDistribution)">inversion
        ///</see> to generate random values.
        /// </summary>
        /// <param name="numeratorDf">the numerator degrees of freedom of the F distribution</param>
        /// <param name="denominatorDf">the denominator degrees of freedom of the F distribution</param>
        /// <returns>value sampled from the F(numeratorDf, denominatorDf) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if
        ///<c>numeratorDf <= 0</c> or <c>denominatorDf <= 0</c>.</exception>
        public double nextF(double numeratorDf, double denominatorDf)
        {
            return delegates.nextF(numeratorDf, denominatorDf);
        }

        /// <summary>
        /// <para>Generates a random value from the
        /// <see cref="GammaDistribution">Gamma Distribution</see>.</para>
        ///<para>This implementation uses the following algorithms: </para>
        ///<para>For 0 < shape < 1: <para/>
        ///Ahrens, J. H. and Dieter, U., Computer methods for
        ///sampling from gamma, beta, Poisson and binomial distributions.
        ///Computing, 12, 223-246, 1974.</para>
        ///<para>For shape >= 1: <para/>
        ///Marsaglia and Tsang, A Simple Method for Generating
        ///Gamma Variables. ACM Transactions on Mathematical Software,
        ///Volume 26 Issue 3, September, 2000.</para>
        /// </summary>
        /// <param name="shape">the median of the Gamma distribution</param>
        /// <param name="scale">the scale parameter of the Gamma distribution</param>
        /// <returns>random value sampled from the Gamma(shape, scale) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>shape <= 0</c> or
        /// <c>scale <= 0</c>.</exception>
        public double nextGamma(double shape, double scale)
        {
            return delegates.nextGamma(shape, scale);
        }

        /// <summary>
        /// Generates a random value from the <see cref="HypergeometricDistribution">
        /// Hypergeometric Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(IntegerDistribution)">
        /// inversion</see> to generate random values.
        ///@throws  
        /// </summary>
        /// <param name="populationSize">the population size of the Hypergeometric distribution</param>
        /// <param name="numberOfSuccesses">number of successes in the population of the
        /// Hypergeometric distribution</param>
        /// <param name="sampleSize">the sample size of the Hypergeometric distribution</param>
        /// <returns>random value sampled from the Hypergeometric(numberOfSuccesses, sampleSize)
        /// distribution</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>numberOfSuccesses > populationSize
        /// </c>, or <c>sampleSize > populationSize</c>.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if <c>populationSize <= 0</c>.</exception>
        /// <exception cref="NotPositiveException"> if <c>numberOfSuccesses < 0</c>.</exception>
        public int nextHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize)
        {
            return delegates.nextHypergeometric(populationSize, numberOfSuccesses, sampleSize);
        }

        /// <summary>
        /// Generates a random value from the <see cref="PascalDistribution">
        /// Pascal Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(IntegerDistribution)">
        /// inversion</see> to generate random values.
        /// </summary>
        /// <param name="r">the number of successes of the Pascal distribution</param>
        /// <param name="p">the probability of success of the Pascal distribution</param>
        /// <returns>random value sampled from the Pascal(r, p) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if the number of successes is not 
        /// positive</exception>
        /// <exception cref="OutOfRangeException"> if the probability of success is not in the
        /// range <c>[0, 1]</c>.</exception>
        public int nextPascal(int r, double p)
        {
            return delegates.nextPascal(r, p);
        }

        /// <summary>
        /// Generates a random value from the <see cref="TDistribution">T-Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(RealDistribution)">inversion
        /// </see> to generate random values.
        /// </summary>
        /// <param name="df">the degrees of freedom of the T distribution</param>
        /// <returns>random value from the T(df) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>df <= 0</c></exception>
        public double nextT(double df)
        {
            return delegates.nextT(df);
        }

        /// <summary>
        /// Generates a random value from the <see cref="WeibullDistribution">Weibull 
        /// Distribution</see>.
        /// This implementation uses <see cref="nextInversionDeviate(RealDistribution)">inversion"
        /// </see> to generate random values.
        /// </summary>
        /// <param name="shape">the shape parameter of the Weibull distribution</param>
        /// <param name="scale">the scale parameter of the Weibull distribution</param>
        /// <returns>random value sampled from the Weibull(shape, size) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>shape <= 0</c> or
        /// <c>scale <= 0</c>.</exception>
        public double nextWeibull(double shape, double scale)
        {
            return delegates.nextWeibull(shape, scale);
        }

        /// <summary>
        /// Generates a random value from the <see cref="ZipfDistribution">Zipf Distribution
        /// </see>.
        /// This implementation uses <see cref="nextInversionDeviate(IntegerDistribution)">
        /// inversion</see> to generate random values.
        /// </summary>
        /// <param name="numberOfElements">the number of elements of the ZipfDistribution</param>
        /// <param name="exponent">the exponent of the ZipfDistribution</param>
        /// <returns>random value sampled from the Zipf(numberOfElements, exponent) distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>numberOfElements <= 0</c>
        /// or <c>exponent <= 0</c>.</exception>
        public int nextZipf(int numberOfElements, double exponent)
        {
            return delegates.nextZipf(numberOfElements, exponent);
        }


        /// <summary>
        /// Reseeds the random number generator with the supplied seed.
        /// <para>
        /// Will create and initialize if null.
        /// </para>       
        /// </summary>
        /// <param name="seed">the seed value to use</param>
        public void reSeed(long seed)
        {
            delegates.reSeed(seed);
        }

        /// <summary>
        /// Reseeds the secure random number generator with the current time in
        /// milliseconds.
        /// <para>
        /// Will create and initialize if null.
        /// </summary>
        public void reSeedSecure()
        {
            delegates.reSeedSecure();
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
            delegates.reSeedSecure(seed);
        }

        /// <summary>
        /// Reseeds the random number generator with
        /// <c>System.currentTimeMillis() + System.identityHashCode(this))</c>.
        /// </summary>
        public void reSeed()
        {
            delegates.reSeed();
        }

        /// <summary>
        /// Sets the PRNG algorithm for the underlying SecureRandom instance using
        /// the Security Provider API. The Security Provider API is defined in <a
        /// href =
        /// "http://java.sun.com/j2se/1.3/docs/guide/security/CryptoSpec.html#AppA">
        /// Java Cryptography Architecture API Specification & Reference.</a>
        /// <para>
        /// USAGE NOTE: This method carries significant
        /// overhead and may take several seconds to execute.
        /// </para>    
        /// </summary>
        /// <param name="algorithm">the name of the PRNG algorithm</param>
        public void setSecureAlgorithm(String algorithm)
        {
            delegates.setSecureAlgorithm(algorithm);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Uses a 2-cycle permutation shuffle. The shuffling process is described <a
        /// href="http://www.maths.abdn.ac.uk/~igc/tch/mx4002/notes/node83.html">
        /// here</a>.
        /// </para>
        /// </remarks>
        public int[] nextPermutation(int n, int k)
        {
            return delegates.nextPermutation(n, k);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// <strong>Algorithm Description</strong>: Uses a 2-cycle permutation
        /// shuffle to generate a random permutation of <c>c.size()</c> and
        /// then returns the elements whose indexes correspond to the elements of the
        /// generated permutation. This technique is described, and proven to
        /// generate random samples <a
        /// href="http://www.maths.abdn.ac.uk/~igc/tch/mx4002/notes/node83.html">
        /// here</a>
        /// </para>
        /// </remarks>
        public Object[] nextSample<T>(Collection<T> c, int k)
        {
            return delegates.nextSample(c, k);
        }

        /// <summary>
        /// Generate a random deviate from the given distribution using the
        ///<a href="http://en.wikipedia.org/wiki/Inverse_transform_sampling"> inversion method.</a>
        /// </summary>
        /// <param name="distribution">Continuous distribution to generate a random value from</param>
        /// <returns>a random value sampled from the given distribution</returns>
        /// <exception cref="MathIllegalArgumentException"> if the underlynig distribution throws one</exception>
        [Obsolete("use the distribution's sample() method")]
        public double nextInversionDeviate(RealDistribution distribution)
        {
            return distribution.inverseCumulativeProbability(nextUniform(0, 1));

        }

        /// <summary>
        /// Generate a random deviate from the given distribution using the
        /// <a href="http://en.wikipedia.org/wiki/Inverse_transform_sampling"> inversion method.</a>
        /// </summary>
        /// <param name="distribution">Integer distribution to generate a random value from</param>
        /// <returns>a random value sampled from the given distribution</returns>
        /// <exception cref="MathIllegalArgumentException"> if the underlynig distribution throws 
        /// one</exception>
        [Obsolete("use the distribution's sample() method")]
        public int nextInversionDeviate(IntegerDistribution distribution)
        {
            return distribution.inverseCumulativeProbability(nextUniform(0, 1));
        }
    }
}
