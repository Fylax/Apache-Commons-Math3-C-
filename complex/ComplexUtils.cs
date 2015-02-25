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

namespace Math3.complex
{
    /// <summary>
    /// Static implementations of common
    /// <see cref="Math3.complex.Complex"/> utilities functions.
    /// </summary>
    public class ComplexUtils
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        private ComplexUtils() { }

        /// <summary>
        /// Creates a complex number from the given polar representation.
        /// <para>
        /// The value returned is <c>r&middot;e^{i&middot;theta}</c>,
        /// computed as <c>r&middot;cos(theta) + r&middot;sin(theta)i</c></para>
        /// <para>
        /// If either <c>r</c> or <c>theta</c> is NaN, or
        /// <c>theta</c> is infinite, <see cref="Complex.NaN"/> is returned.</para>
        /// <para>
        /// If <c>r</c> is infinite and <c>theta</c> is finite,
        /// infinite or NaN values may be returned in parts of the result, following
        /// the rules for double arithmetic.
        /// </para>
        /// </summary>
        /// <param name="r">the modulus of the complex number to create</param>
        /// <param name="theta">the argument of the complex number to create</param>
        /// <returns><c>r&middot;e^{i&middot;theta}</x></returns>
        /// <example>
        /// <code>
        /// polar2Complex(INFINITY, &pi;/4) = INFINITY + INFINITY i
        /// polar2Complex(INFINITY, 0) = INFINITY + NaN i
        /// polar2Complex(INFINITY, -&pi;/4) = INFINITY - INFINITY i
        /// polar2Complex(INFINITY, 5&pi;/4) = -INFINITY - INFINITY i
        /// </code></example>
        /// <exception cref="MathIllegalArgumentException">if <c>r</c> is negative.</exception>
        public static Complex polar2Complex(double r, double theta)
        {
            if (r < 0)
            {
                throw new MathIllegalArgumentException(
                      new LocalizedFormats("NEGATIVE_COMPLEX_MODULE"), r);
            }
            return new Complex(r * FastMath.cos(theta), r * FastMath.sin(theta));
        }

        /// <summary>
        /// Convert an array of primitive doubles to an array of <c>Complex</c> objects.
        /// </summary>
        /// <param name="real">Array of numbers to be converted to their
        /// <c>Complex</c> equivalent</param>
        /// <returns>an array of <c>Complex</c>} objects.</returns>
        public static Complex[] convertToComplex(double[] real)
        {
            Complex[] c = new Complex[real.Length];
            for (int i = 0; i < real.Length; i++)
            {
                c[i] = new Complex(real[i], 0);
            }

            return c;
        }
    }
}
