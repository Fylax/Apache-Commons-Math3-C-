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
using Math3.exception.util;
using Math3.util;
using System;

namespace Math3.complex
{
    /// <summary>
    /// A helper class for the computation and caching of the <c>n-th</c> roots of
    /// unity.
    /// </summary>
    /// <remarks>Altough in Java version methods where syncronized, I decided
    /// to make them just standard method as syncronized methods are not supported
    /// by mobile.</remarks>
    public class RootsOfUnity
    {

        /// <summary>
        /// Number of roots of unity.
        /// </summary>
        private int omegaCount;

        /// <summary>
        /// Real part of the roots.
        /// </summary>
        private double[] omegaReal;

        /// <summary>
        /// Imaginary part of the <c>n-th</c> roots of unity, for positive values
        /// of <c>n</c>. In this array, the roots are stored in counter-clockwise
        /// order.
        /// </summary>
        private double[] omegaImaginaryCounterClockwise;

        /// <summary>
        /// Imaginary part of the <c>n-th</c> roots of unity, for negative values
        /// of <c>n</c>. In this array, the roots are stored in clockwise order.
        /// </summary>
        private double[] omegaImaginaryClockwise;

        /// <summary>
        /// <c>true</c> if <see cref="computeRoots(int)"/> was called with a positive
        /// value of its argument <c>n</c>. In this case, counter-clockwise ordering
        /// of the roots of unity should be used.
        /// </summary>
        private Boolean isCounterClockWise;

        /// <summary>
        /// Build an engine for computing the <c>n-th</c> roots of unity.
        /// </summary>
        public RootsOfUnity()
        {
            omegaCount = 0;
            omegaReal = null;
            omegaImaginaryCounterClockwise = null;
            omegaImaginaryClockwise = null;
            isCounterClockWise = true;
        }

        /// <summary>
        /// Returns <c>true</c> if <see cref="computeRoots(int)"/> was called with a
        /// positive value of its argument <c>n</c>. If <c>true</c>, then
        /// counter-clockwise ordering of the roots of unity should be used.
        /// </summary>
        /// <returns><c>true</c> if the roots of unity are stored in
        /// counter-clockwise order</returns>
        /// <exception cref="MathIllegalStateException"> if no roots of unity have been computed
        /// yet</exception>
        public Boolean IsCounterClockWise()
        {
            lock (this)
            {
                if (omegaCount == 0)
                {
                    throw new MathIllegalStateException(new LocalizedFormats("ROOTS_OF_UNITY_NOT_COMPUTED_YET"));
                }
                return isCounterClockWise;
            }
        }

        /// <summary>
        /// Computes the <c>n-th</c> roots of unity. The roots are stored in
        /// <c>omega[]</c>, such that <c>omega[k] = w ^ k</c>, where
        /// <c>k = 0, ..., n - 1</c>, <c>w = exp(2 * pi * i / n)</c> and
        /// <c>i = sqrt(-1)</c>.
        /// <para>
        /// Note that <c>n</c> can be positive of negative
        /// </para>
        /// <list type="bullet">
        /// <item><c>abs(n)</c> is always the number of roots of unity.</item>
        /// <item>If <c>n > 0</c>, then the roots are stored in counter-clockwise order.</item>
        /// <item>If <c>n < 0</c>, then the roots are stored in clockwise order.</item>
        /// </para>
        /// </list>
        /// </summary>
        /// <param name="n">the (signed) number of roots of unity to be computed</param>
        /// <exception cref="ZeroException"> if <c>n = 0</c></exception>
        public void computeRoots(int n)
        {
            lock (this)
            {
                if (n == 0)
                {
                    throw new ZeroException(new LocalizedFormats("CANNOT_COMPUTE_0TH_ROOT_OF_UNITY"));
                }
                isCounterClockWise = n > 0;
                // avoid repetitive calculations
                int absN = FastMath.abs(n);
                if (absN == omegaCount)
                {
                    return;
                }
                // calculate everything from scratch
                double t = 2.0 * FastMath.PI / absN;
                double cosT = FastMath.cos(t);
                double sinT = FastMath.sin(t);
                omegaReal = new double[absN];
                omegaImaginaryCounterClockwise = new double[absN];
                omegaImaginaryClockwise = new double[absN];
                omegaReal[0] = 1.0;
                omegaImaginaryCounterClockwise[0] = 0.0;
                omegaImaginaryClockwise[0] = 0.0;
                for (int i = 1; i < absN; i++)
                {
                    omegaReal[i] = omegaReal[i - 1] * cosT - omegaImaginaryCounterClockwise[i - 1] * sinT;
                    omegaImaginaryCounterClockwise[i] = omegaReal[i - 1] * sinT + omegaImaginaryCounterClockwise[i - 1] * cosT;
                    omegaImaginaryClockwise[i] = -omegaImaginaryCounterClockwise[i];
                }
                omegaCount = absN;
            }
        }

        /// <summary>
        /// Get the real part of the <c>k-th</c> <c>n-th</c> root of unity.
        /// </summary>
        /// <param name="k">index of the <c>n-th</c> root of unity</param>
        /// <returns>real part of the <c>k-th</c> <c>n-th</c> root of unity</returns>
        /// <exception cref="MathIllegalStateException"> if no roots of unity have been
        /// computed yet</exception>
        /// <exception cref="MathIllegalArgumentException"> if <c>k</c> is out of
        /// range</exception>
        public double getReal(int k)
        {
            lock (this)
            {
                if (omegaCount == 0)
                {
                    throw new MathIllegalStateException(new LocalizedFormats("ROOTS_OF_UNITY_NOT_COMPUTED_YET"));
                }
                if ((k < 0) || (k >= omegaCount))
                {
                    throw new OutOfRangeException<Int32>(new LocalizedFormats("OUT_OF_RANGE_ROOT_OF_UNITY_INDEX"), k, 0, omegaCount - 1);
                }

                return omegaReal[k];
            }
        }

        /// <summary>
        /// Get the imaginary part of the <c>k-th</c> <c>n-th</c> root of unity.
        /// </summary>
        /// <param name="k">index of the <c>n-th</c> root of unity</param>
        /// <returns>imaginary part of the <c>k-th</c> <c>n-th</c> root of unity</returns>
        /// <exception cref="MathIllegalStateException"> if no roots of unity have been
        /// computed yet</exception>
        /// <exception cref="MathIllegalArgumentException"> if <c>k</c> is out of
        /// range</exception>
        public double getImaginary(int k)
        {
            lock (this)
            {
                if (omegaCount == 0)
                {
                    throw new MathIllegalStateException(new LocalizedFormats("ROOTS_OF_UNITY_NOT_COMPUTED_YET"));
                }
                if ((k < 0) || (k >= omegaCount))
                {
                    throw new OutOfRangeException<Int32>(new LocalizedFormats("OUT_OF_RANGE_ROOT_OF_UNITY_INDEX"), k, 0, omegaCount - 1);
                }

                return isCounterClockWise ? omegaImaginaryCounterClockwise[k] :
                    omegaImaginaryClockwise[k];
            }
        }

        /// <summary>
        /// Returns the number of roots of unity currently stored. If
        /// <see cref="computeRoots(int)"/> was called with <c>n</c>, then this method
        /// returns <c>abs(n)</c>. If no roots of unity have been computed yet, this
        /// method returns 0.
        /// </summary>
        /// <returns>the number of roots of unity currently stored</returns>
        public int getNumberOfRoots()
        {
            lock (this)
            {
                return omegaCount;
            }
        }
    }
}
