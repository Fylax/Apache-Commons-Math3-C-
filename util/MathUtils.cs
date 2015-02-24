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
using System;

namespace Math3.util
{
    /// <summary>
    /// Miscellaneous utility functions.
    /// </summary>
    /// <remarks>
    /// See <seealso cref="ArithmeticUtils"/>, <seealso cref="Precision"/>,
    /// <seealso cref="MathArrays"/>
    /// </remarks>
    public class MathUtils
    {
        /// <summary>
        /// \(2\pi\)
        /// </summary>
        public const double TWO_PI = 2 * FastMath.PI;

        /// <summary>
        /// \(\pi^2\)
        /// </summary>
        public const double PI_SQUARED = FastMath.PI * FastMath.PI;


        /// <summary>
        /// Class contains only static methods.
        /// </summary>
        private MathUtils() { }


        /// <summary>
        /// Returns an integer hash code representing the given double value.
        /// </summary>
        /// <param name="value">the value to be hashed</param>
        /// <returns>the hash code</returns>
        public static int hash(double value)
        {
            return value.GetHashCode();
        }

        /// <summary>
        /// Returns <c>true</c> if the values are equal according to semantics of
        /// <see cref="System.Double.Equals(Object, Object)"/>.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="y">y Value</param>
        /// <returns><c>Double.Equals(x, y);</c></returns>
        public static Boolean equals(double x, double y)
        {
            return Double.Equals(x, y);
        }

        /// <summary>
        /// Returns an integer hash code representing the given double array.
        /// </summary>
        /// <param name="value">the value to be hashed (may be null)</param>
        /// <returns>the hash code</returns>
        public static int hash(double[] value)
        {
            return value.GetHashCode();
        }

        /// <summary>
        /// Normalize an angle in a 2&pi; wide interval around a center value.
        /// <para>This method has three main uses:</para>
        /// <list type="bullet">
        /// <item>normalize an angle between 0 and 2&pi;:<para/>
        /// <c>a = MathUtils.normalizeAngle(a, FastMath.PI);</c></item>
        /// <item>normalize an angle between -&pi; and +&pi;<para/>
        /// <c>a = MathUtils.normalizeAngle(a, 0.0);</c></item>
        /// <item>compute the angle between two defining angular positions:<para/>
        /// <c>angle = MathUtils.normalizeAngle(end, start) - start;</c></item>
        /// </list>
        /// <para>Note that due to numerical accuracy and since &pi; cannot be represented
        /// exactly, the result interval is closed, it cannot be half-closed
        /// as would be more satisfactory in a purely mathematical view.</para>
        /// </summary>
        /// <param name="a">angle to normalize</param>
        /// <param name="center">center of the desired 2&pi; interval for the result</param>
        /// <returns>a-2k&pi; with integer k and center-&pi; &lt;= a-2k&pi; &lt;= center+&pi;</returns>
        public static double normalizeAngle(double a, double center)
        {
            return a - TWO_PI * FastMath.floor((a + FastMath.PI - center) / TWO_PI);
        }

        /// <summary>
        /// <para>Reduce <c>|a - offset|</c> to the primary interval
        /// <c>[0, |period|)</c>.</para>
        /// <para>Specifically, the value returned is <aara/>
        /// <c>a - |period| * floor((a - offset) / |period|) - offset</c>.</para>
        /// <para>If any of the parameters are <c>NaN</c> or infinite, the result is
        /// <c>NaN</c>.</para>
        /// </summary>
        /// <param name="a">Value to reduce.</param>
        /// <param name="period">Value that will be mapped to <c>0</c>.Period.</param>
        /// <param name="offset"></param>
        /// <returns>the value, within the interval <c>[0 |period|)</c>,
        /// that corresponds to <c>a</c>.</returns>
        public static double reduce(double a,
                                    double period,
                                    double offset)
        {
            double p = FastMath.abs(period);
            return a - p * FastMath.floor((a - offset) / p) - offset;
        }

        /// <summary>
        /// Returns the first argument with the sign of the second argument.
        /// </summary>
        /// <param name="magnitude">Magnitude of the returned value.</param>
        /// <param name="sign">ign of the returned value.</param>
        /// <returns>a value with magnitude equal to <c>magnitude</c> and with the
        /// same sign as the <c>sign</c> argument.</returns>
        /// <exception cref="MathArithmeticException"> if <c>magnitude == Byte.MinValue</c>
        /// and <c>sign >= 0</c>.</exception>
        public static byte copySign(byte magnitude, byte sign)
        {
            if ((magnitude >= 0 && sign >= 0) ||
                (magnitude < 0 && sign < 0))
            { // Sign is OK.
                return magnitude;
            }
            else if (sign >= 0 &&
                     magnitude == Byte.MinValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("OVERFLOW"));
            }
            else
            {
                return (byte)-magnitude; // Flip sign.
            }
        }

        /// <summary>
        /// Returns the first argument with the sign of the second argument.
        /// </summary>
        /// <param name="magnitude">Magnitude of the returned value.</param>
        /// <param name="sign">ign of the returned value.</param>
        /// <returns>a value with magnitude equal to <c>magnitude</c> and with the
        /// same sign as the <c>sign</c> argument.</returns>
        /// <exception cref="MathArithmeticException"> if <c>magnitude == Int16.MinValue</c>
        /// and <c>sign >= 0</c>.</exception>
        public static short copySign(short magnitude, short sign)
        {
            if ((magnitude >= 0 && sign >= 0) ||
                (magnitude < 0 && sign < 0))
            { // Sign is OK.
                return magnitude;
            }
            else if (sign >= 0 &&
                     magnitude == Int16.MinValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("OVERFLOW"));
            }
            else
            {
                return (short)-magnitude; // Flip sign.
            }
        }

        /// <summary>
        /// Returns the first argument with the sign of the second argument.
        /// </summary>
        /// <param name="magnitude">Magnitude of the returned value.</param>
        /// <param name="sign">ign of the returned value.</param>
        /// <returns>a value with magnitude equal to <c>magnitude</c> and with the
        /// same sign as the <c>sign</c> argument.</returns>
        /// <exception cref="MathArithmeticException"> if <c>magnitude == Int32.MinValue</c>
        /// and <c>sign >= 0</c>.</exception>
        public static int copySign(int magnitude, int sign)
        {
            if ((magnitude >= 0 && sign >= 0) ||
                (magnitude < 0 && sign < 0))
            { // Sign is OK.
                return magnitude;
            }
            else if (sign >= 0 &&
                     magnitude == Int32.MinValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("OVERFLOW"));
            }
            else
            {
                return -magnitude; // Flip sign.
            }
        }

        /// <summary>
        /// Returns the first argument with the sign of the second argument.
        /// </summary>
        /// <param name="magnitude">Magnitude of the returned value.</param>
        /// <param name="sign">ign of the returned value.</param>
        /// <returns>a value with magnitude equal to <c>magnitude</c> and with the
        /// same sign as the <c>sign</c> argument.</returns>
        /// <exception cref="MathArithmeticException"> if <c>magnitude == Int64.MinValue</c>
        /// and <c>sign >= 0</c>.</exception>
        public static long copySign(long magnitude, long sign)
        {
            if ((magnitude >= 0 && sign >= 0) ||
                (magnitude < 0 && sign < 0))
            { // Sign is OK.
                return magnitude;
            }
            else if (sign >= 0 &&
                     magnitude == Int16.MinValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("OVERFLOW"));
            }
            else
            {
                return -magnitude; // Flip sign.
            }
        }

        /// <summary>
        /// Check that the argument is a real number.
        /// </summary>
        /// <param name="x">Argument.</param>
        /// <exception cref="NotFiniteNumberException">if <c>x</c> is not a
        /// finite real number.</exception>
        public static void checkFinite(double x)
        {
            if (Double.IsInfinity(x) || Double.IsNaN(x))
            {
                throw new NotFiniteNumberException<Double>(x);
            }
        }

        /// <summary>
        /// Check that all the elements are real numbers.
        /// </summary>
        /// <param name="val">Arguments.</param>
        /// <exception cref="NotFiniteNumberException">if any values of the array is not a
        /// finite real number.</exception>
        public static void checkFinite(double[] val)
        {
            for (int i = 0; i < val.Length; i++)
            {
                double x = val[i];
                if (Double.IsInfinity(x) || Double.IsNaN(x))
                {
                    throw new NotFiniteNumberException<Double>(new LocalizedFormats("ARRAY_ELEMENT"), x, i);
                }
            }
        }

        /// <summary>
        /// Checks that an object is not null.
        /// </summary>
        /// <param name="o">Object to be checked.</param>
        /// <param name="pattern">Message pattern.</param>
        /// <param name="args">Arguments to replace the placeholders in <c>pattern</c>.</param>
        /// <exception cref="NullArgumentException">if <c>o</c> is <c>null</c>.</exception>
        public static void checkNotNull(Object o, Localizable pattern, params Object[] args)
        {
            if (o == null)
            {
                throw new NullArgumentException(pattern, args);
            }
        }

        /// <summary>
        /// Checks that an object is not null.
        /// </summary>
        /// <param name="o">Object to be checked.</param>
        /// <exception cref="NullArgumentException">if <c>o</c> is <c>null</c>.</exception>
        public static void checkNotNull(Object o)
        {
            if (o == null)
            {
                throw new NullArgumentException();
            }
        }
    }
}
