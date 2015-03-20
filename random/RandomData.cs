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
using System.Collections.ObjectModel;

namespace Math3.random
{
    /// <summary>
    /// Random data generation utilities.
    /// </summary>
    [Obsolete]
    public interface RandomData
    {
        /// <summary>
        /// Generates a random string of hex characters of length <c>len</c>.
        /// <para>
        /// The generated string will be random, but not cryptographically
        /// secure. To generate cryptographically secure strings, use
        /// <see cref="nextSecureHexString(int)"/>.
        /// </para>
        /// </summary>
        /// <param name="len">the length of the string to be generated</param>
        /// <returns>a random string of hex characters of length <c>len</c></returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>len <= 0</c><exception>
        String nextHexString(int len);

        /// <summary>
        /// Generates a uniformly distributed random integer between <c>lower</c>
        /// and <c>upper</c> (endpoints included).
        /// <para>
        /// The generated integer will be random, but not cryptographically secure.
        /// To generate cryptographically secure integer sequences, use
        /// <see cref="nextSecureInt(int, int)"/>.
        /// </para>
        /// </summary>
        /// <param name="lower">lower bound for generated integer</param>
        /// <param name="upper">upper bound for generated integer</param>
        /// <returns>a random integer greater than or equal to <c>lower</c>
        /// and less than or equal to <c>upper</c></returns>
        /// <exception cref="NumberIsTooLargeException">if <c>lower >= upper</c></exception>
        int nextInt(int lower, int upper);

        /// <summary>
        /// Generates a uniformly distributed random long integer between
        /// <c>lower</c> and <c>upper</c> (endpoints included).
        /// <para>
        /// The generated long integer values will be random, but not
        /// cryptographically secure. To generate cryptographically secure sequences
        /// of longs, use <see cref="nextSecureLong(long, long)"/>.
        /// </para>
        /// </summary>
        /// <param name="lower">lower bound for generated long integer</param>
        /// <param name="upper">upper bound for generated long integer</param>
        /// <returns>a random long integer greater than or equal to <c>lower</c> and
        /// less than or equal to <c>upper</c></returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c></exception>
        long nextLong(long lower, long upper);

        /// <summary>
        /// Generates a random string of hex characters from a secure random
        /// sequence.
        /// <para>
        /// If cryptographic security is not required, use
        /// <see cref="nextHexString(int)"/>.
        /// </para>
        /// </summary>
        /// <param name="len">the length of the string to be generated</param>
        /// <returns>a random string of hex characters of length <c>len</c></returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>len <= 0</c></exception>
        String nextSecureHexString(int len);

        /// <summary>
        /// Generates a uniformly distributed random integer between <c>lower</c>
        /// and <c>upper</c> (endpoints included) from a secure random sequence.
        /// <para>
        /// Sequences of integers generated using this method will be
        /// cryptographically secure. If cryptographic security is not required,
        /// <see cref="nextInt(int, int)"/> should be used instead of this method.</para>
        /// <para>
        /// Definition:
        /// a href="http://en.wikipedia.org/wiki/Cryptographically_secure_pseudo-random_number_generator">
        /// Secure Random Sequence</a></para> 
        /// </summary>
        /// <param name="lower">lower bound for generated integer</param>
        /// <param name="upper">upper bound for generated integer</param>
        /// <returns>a random integer greater than or equal to <c>lower</c> and less
        /// than or equal to <c>upper</c>.</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c>.</exception>
        int nextSecureInt(int lower, int upper);

        /// <summary>
        /// Generates a uniformly distributed random long integer between
        /// <c>lower</c> and <c>upper</c> (endpoints included) from a secure random
        /// sequence.
        /// <para>
        /// Sequences of long values generated using this method will be
        /// cryptographically secure. If cryptographic security is not required,
        /// <see cref="nextLong(long, long)"/> should be used instead of this method.</para>
        /// <para>
        /// Definition:
        /// <a href="http://en.wikipedia.org/wiki/Cryptographically_secure_pseudo-random_number_generator">
        /// Secure Random Sequence</a></para>
        /// </summary>
        /// <param name="lower">lower bound for generated integer</param>
        /// <param name="upper">upper bound for generated integer</param>
        /// <returns>a random long integer greater than or equal to <c>lower</c> and
        /// less than or equal to <c>upper</c>.</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c>.</exception>
        long nextSecureLong(long lower, long upper);

        /// <summary>
        /// Generates a random value from the Poisson distribution with the given
        /// mean.
        /// <para>
        /// Definition:
        /// <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda366j.htm">
        /// Poisson Distribution</a></para> 
        /// </summary>
        /// <param name="mean">the mean of the Poisson distribution</param>
        /// <returns>a random value following the specified Poisson distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>mean <= 0</c>.</exception>
        long nextPoisson(double mean);

        /// <summary>
        /// Generates a random value from the Normal (or Gaussian) distribution with
        /// specified mean and standard deviation.
        /// <para>
        /// Definition:
        /// <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3661.htm">
        /// Normal Distribution</a></para>
        /// </summary>
        /// <param name="mu">the mean of the distribution</param>
        /// <param name="sigma">the standard deviation of the distribution</param>
        /// <returns>a random value following the specified Gaussian distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sigma <= 0</c>.</exception>
        double nextGaussian(double mu, double sigma);

        /// <summary>
        /// Generates a random value from the exponential distribution
        /// with specified mean.
        /// <para>
        /// Definition:
        /// <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3667.htm">
        /// Exponential Distribution</a></para> 
        /// </summary>
        /// <param name="mean">the mean of the distribution</param>
        /// <returns>a random value following the specified exponential distribution</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>mean <= 0</c>.</ecxeption>
        double nextExponential(double mean);

        /// <summary>
        /// Generates a uniformly distributed random value from the open interval
        /// <c>(lower, upper)</c> (i.e., endpoints excluded).
        /// <para>
        /// Definition:
        /// <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm">
        /// Uniform Distribution</a> <c>lower</c> and <c>upper - lower</c> are the
        /// <a href = "http://www.itl.nist.gov/div898/handbook/eda/section3/eda364.htm">
        /// location and scale parameters</a>, respectively.</para>
        /// </summary>
        /// <param name="lower">the exclusive lower bound of the support</param>
        /// <param name="upper">the exclusive upper bound of the support</param>
        /// <returns>a uniformly distributed random value between lower and upper
        /// (exclusive)</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c></exception>
        /// <exception cref="NotFiniteNumberException"> if one of the bounds is infinite</exception>
        /// <exception cref="NotANumberException"> if one of the bounds is NaN</exception>
        double nextUniform(double lower, double upper);

        /// <summary>
        /// Generates a uniformly distributed random value from the interval
        /// <c>(lower, upper)</c> or the interval <c>[lower, upper)</c>. The lower
        /// bound is thus optionally included, while the upper bound is always
        /// excluded.
        /// <para>
        /// Definition:
        /// <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm">
        /// Uniform Distribution</a> <c>lower</c> and <c>upper - lower</c> are the
        /// <a href = "http://www.itl.nist.gov/div898/handbook/eda/section3/eda364.htm">
        /// location and scale parameters</a>, respectively.</para>
        /// </summary>
        /// <param name="lower">the lower bound of the support</param>
        /// <param name="upper">the exclusive upper bound of the support</param>
        /// <param name="lowerInclusive">lowerInclusive <c>true</c> if the lower bound is inclusive</param>
        /// <returns>uniformly distributed random value in the <c>(lower, upper)</c>
        /// interval, if <c>lowerInclusive</c> is <c>false</c>, or in the
        /// <c>[lower, upper)</c> interval, if <c>lowerInclusive</c> is
        /// <c>true</c></returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>lower >= upper</c></exception>
        /// <exception cref="NotFiniteNumberException"> if one of the bounds is infinite</exception>
        /// <exception cref="NotANumberException"> if one of the bounds is NaN</exception>
        double nextUniform(double lower, double upper, Boolean lowerInclusive);

        /// <summary>
        /// Generates an integer array of length <c>k</c> whose entries are selected
        /// randomly, without repetition, from the integers <c>0, ..., n - 1</c>
        /// (inclusive).
        /// <para>
        /// Generated arrays represent permutations of <c>n</c> taken <c>k</c> at a
        /// time.</para>
        /// </summary>
        /// <param name="n">the domain of the permutation</param>
        /// <param name="k">the size of the permutation</param>
        /// <returns>a random <c>-permutation of <c>n</c>, as an array of
        /// integers</returns>
        /// <exception cref="NumberIsTooLargeException">if <c>k > n</c>.</cref>
        /// <exception cref="NotStrictlyPositiveException"> if <c>k <= 0</c>.</cref>
        int[] nextPermutation(int n, int k);

        /// <summary>
        /// Returns an array of <c>k</c> objects selected randomly from the
        /// Collection <c>c</c>.
        /// <para>
        /// Sampling from <c>c</c> is without replacement; but if <c>c</c> contains
        /// identical objects, the sample may include repeats.  If all elements of
        /// <c>c</c> are distinct, the resulting object array represents a
        /// <a href="http://rkb.home.cern.ch/rkb/AN16pp/node250.html#SECTION0002500000000000000000">
        /// Simple Random Sample</a> of size <c>k</c> from the elements of
        /// <c>c</c>.</para>
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="c">the collection to be sampled</param>
        /// <param name="k">the size of the sample</param>
        /// <returns>a random sample of <c>k</c> elements from <c>c</c></returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > c.size()</c>.<exception>
        /// <exception cref="NotStrictlyPositiveException"> if <c>k <= 0</c>.<exception>
        Object[] nextSample<T>(Collection<T> c, int k);
    }
}
