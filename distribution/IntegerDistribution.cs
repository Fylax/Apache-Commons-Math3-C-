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

namespace Math3.distribution
{
    /// <summary>
    /// Interface for distributions on the integers.
    /// </summary>
    public interface IntegerDistribution
    {
        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(X = x)</c>. In other
        /// words, this method represents the probability mass function (PMF)
        /// for the distribution.
        /// </summary>
        /// <param name="x">the point at which the PMF is evaluated</param>
        /// <returns>the value of the probability mass function at <c>x</c></returns>
        double probability(int x);

        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(X <= x)</c>.  In other
        /// words, this method represents the (cumulative) distribution function
        /// (CDF) for this distribution.
        /// </summary>
        /// <param name="x">the point at which the CDF is evaluated</param>
        /// <returns>the probability that a random variable with this
        /// distribution takes a value less than or equal to <c>x</c></returns>
        double cumulativeProbability(int x);

        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(x0 < X <= x1)</c>.
        /// </summary>
        /// <param name="x0">the exclusive lower bound</param>
        /// <param name="x1">the inclusive upper bound</param>
        /// <returns>the probability that a random variable with this distribution
        /// will take a value between <c>x0</c> and <c>x1</c>,
        /// excluding the lower and including the upper endpoint</returns>
        /// <exception cref="NumberIsTooLargeException"> if <c>x0 > x1</c></exception>
        double cumulativeProbability(int x0, int x1);

        /// <summary>
        /// Computes the quantile function of this distribution.
        /// For a random variable <c>X</c> distributed according to this distribution,
        /// the returned value is
        /// <list type="bullet">
        /// <item><c>inf{x in Z | P(X<=x) >= p}</c> for <c>0 < p <= 1</c>,</item>
        /// <item><c>inf{x in Z | P(X<=x) > 0}</code> for <c>p = 0</c>.</item>
        /// </list>
        /// If the result exceeds the range of the data type <c>int</c>,
        /// then <c>Int32.MimValue</c> or <c>Int32.MaxValue</c> is returned.
        /// </summary>
        /// <param name="p">the cumulative probability</param>
        /// <returns>the smallest <c>p</c>-quantile of this distribution
        /// (largest 0-quantile for <c>p = 0</c>)</returns>
        /// <exception cref="OutOfRangeException">if <c>p < 0</c> or <c>p > 1</c></exception>
        int inverseCumulativeProbability(double p);

        /// <summary>
        /// Use this method to get the numerical value of the mean of this
        /// distribution.
        /// </summary>
        /// <returns>the mean or <c>Double.NaN</c> if it is not defined</returns>
        double getNumericalMean();

        /// <summary>
        /// Use this method to get the numerical value of the variance of this
        /// distribution.
        /// </summary>
        /// <returns>the variance (possibly <c>Double.PositiveInfinity</c> or
        /// <c>Double.NaN</c> if it is not defined)</returns>
        double getNumericalVariance();

        /// <summary>
        /// Access the lower bound of the support. This method must return the same
        /// value as <c>inverseCumulativeProbability(0)</c>. In other words, this
        /// method must return
        /// <para><c>inf {x in Z | P(X <= x) > 0}</c>.</para>
        /// </summary>
        /// <returns>lower bound of the support (<c>Int32.MinValue</c>
        /// for negative infinity)</returns>
        int getSupportLowerBound();

        /// <summary>
        /// Access the upper bound of the support. This method must return the same
        /// value as <c>inverseCumulativeProbability(1)</c>. In other words, this
        /// method must return
        /// <para><c>inf {x in R | P(X <= x) = 1}</c>.</para>
        /// </summary>
        /// <returns>upper bound of the support (<c>Int32.MaxValue</c>
        /// for positive infinity)</returns>
        int getSupportUpperBound();

        /// <summary>
        /// Use this method to get information about whether the support is
        /// connected, i.e. whether all integers between the lower and upper bound of
        /// the support are included in the support.
        /// </summary>
        /// <returns>whether the support is connected or not</returns>
        Boolean isSupportConnected();

        /// <summary>
        /// Reseed the random generator used to generate samples.
        /// </summary>
        /// <param name="seed">seed the new seed</param>
        void reseedRandomGenerator(long seed);

        /// <summary>
        /// Generate a random value sampled from this distribution. 
        /// </summary>
        /// <returns>a random value</returns>
        int sample();

        /// <summary>
        /// Generate a random sample from the distribution.
        /// </summary>
        /// <param name="sampleSize">the number of random values to generate</param>
        /// <returns>an array representing the random sample</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sampleSize</c> is not positive</exception>
        int[] sample(int sampleSize);
    }
}
