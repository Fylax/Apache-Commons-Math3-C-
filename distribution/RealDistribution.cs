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
    /// Base interface for distributions on the reals.
    /// </summary>
    public interface RealDistribution
    {
        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(X = x)</c>. In other
        /// words, this method represents the probability mass function (PMF)
        /// for the distribution.
        /// </summary>
        /// <param name="x">the point at which the PMF is evaluated</param>
        /// <returns>the value of the probability mass function at point <c>x</c></returns>
        double probability(double x);

        /// <summary>
        /// Returns the probability density function (PDF) of this distribution
        /// evaluated at the specified point <c>x</c>. In general, the PDF is
        /// the derivative of the <see cref="cumulativeProbability(double)">CDF</see>.
        /// If the derivative does not exist at <c>x</c>, then an appropriate
        /// replacement should be returned, e.g. <c>Double.PositiveInfinity</c>,
        /// <c>Double.NaN</c>, or  the limit inferior or limit superior of the
        /// difference quotient.
        /// </summary>
        /// <param name="x">the point at which the PDF is evaluated</param>
        /// <returns>the value of the probability density function at point <c>x</c></returns>
        double density(double x);

        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(X <= x)</c>. In other
        /// words, this method represents the (cumulative) distribution function
        /// (CDF) for this distribution.
        /// </summary>
        /// <param name="x">the point at which the CDF is evaluated</param>
        /// <returns>the probability that a random variable with this
        /// distribution takes a value less than or equal to <c>x</c></returns>
        double cumulativeProbability(double x);

        /// <summary>
        /// For a random variable <c>X</c> whose values are distributed according
        /// to this distribution, this method returns <c>P(x0 < X <= x1)</c>.
        /// </summary>
        /// <param name="x0">the exclusive lower bound</param>
        /// <param name="x1">the inclusive upper bound</param>
        /// <returns>the probability that a random variable with this distribution
        /// takes a value between <c>x0</c> and <c>x1</c>,
        /// excluding the lower and including the upper endpoint</returns>
        /// <exception cref="NumberIsTooLargeException">if <c>x0 > x1</c></exception>
        [Obsolete]
        double cumulativeProbability(double x0, double x1);

        /// <summary>
        /// Computes the quantile function of this distribution. For a random
        /// variable <c>X</c> distributed according to this distribution, the
        /// returned value is
        /// <list type="bullet">
        /// <item><c>inf{x in R | P(X<=x) >= p}</c> for <c>0 < p <= 1</c>,</item>
        /// <item><c>inf{x in R | P(X<=x) > 0}</c> for <c>p = 0</c>.</item>
        /// </list>
        /// </summary>
        /// <param name="p">the cumulative probability</param>
        /// <returns>the smallest <c>p</c>-quantile of this distribution
        /// (largest 0-quantile for <c>p = 0</c>)</returns>
        /// <exception cref="OutOfRangeException"> if <c>p < 0</c> or <c>p > 1</c></exception>
        double inverseCumulativeProbability(double p);

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
        /// <returns>the variance (possibly <c>Double.PositiveInfinity</c> as
        /// for certain cases in <see cref="TDistribution"/>) or <c>Double.NaN</c> if it
        /// is not defined</returns>
        double getNumericalVariance();

        /// <summary>
        /// Access the lower bound of the support. This method must return the same
        /// value as <c>inverseCumulativeProbability(0)</c>. In other words, this
        /// method must return
        /// <para><c>inf {x in R | P(X <= x) > 0}</c>.</para>
        /// </summary>
        /// <returns>lower bound of the support (might be
        /// <c>Double.NegativeInfinity</c>)</returns>
        double getSupportLowerBound();

        /// <summary>
        /// Access the upper bound of the support. This method must return the same
        /// value as <c>inverseCumulativeProbability(1)</c>. In other words, this
        /// method must return
        /// <para><c>inf {x in R | P(X <= x) = 1}</c>.</para>
        /// </summary>
        /// <returns>upper bound of the support (might be
        /// <c>Double.PositiveInfinity</c>)</returns>
        double getSupportUpperBound();

        /// <summary>
        /// Whether or not the lower bound of support is in the domain of the density
        /// function.  Returns true iff <c>getSupporLowerBound()</c> is finite and
        /// <c>density(getSupportLowerBound())</c> returns a non-NaN, non-infinite
        /// value.
        /// </summary>
        /// <returns>true if the lower bound of support is finite and the density
        /// function returns a non-NaN, non-infinite value there</returns>
        [Obsolete]
        Boolean isSupportLowerBoundInclusive();

        /// <summary>
        /// Whether or not the upper bound of support is in the domain of the density
        /// function.  Returns true iff <c>getSupportUpperBound()</c> is finite and
        /// <c>density(getSupportUpperBound())</c> returns a non-NaN, non-infinite
        /// value.
        /// </summary>
        /// <returns>true if the upper bound of support is finite and the density
        /// function returns a non-NaN, non-infinite value there</returns>
        [Obsolete]
        Boolean isSupportUpperBoundInclusive();

        /// <summary>
        /// Use this method to get information about whether the support is connected,
        /// i.e. whether all values between the lower and upper bound of the support
        /// are included in the support.
        /// </summary>
        /// <returns>whether the support is connected or not</returns>
        Boolean isSupportConnected();

        /// <summary>
        /// Reseed the random generator used to generate samples.
        /// </summary>
        /// <param name="seed">the new seed</param>
        void reseedRandomGenerator(long seed);

        /// <summary>
        /// Generate a random value sampled from this distribution.
        /// </summary>
        /// <returns>a random value.</returns>
        double sample();

        /// <summary>
        /// Generate a random sample from the distribution.
        /// </summary>
        /// <param name="sampleSize">the number of random values to generate</param>
        /// <returns>an array representing the random sample</returns>
        /// <exception cref="NotStrictlyPositiveException"> if <c>sampleSize</c> is not positive</exception>
        double[] sample(int sampleSize);
    }
}
