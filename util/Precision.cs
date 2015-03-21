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
    /// Utilities for comparing numbers.
    /// </summary>
    public class Precision
    {
        /// <summary>
        /// <para>
        /// Largest double-precision floating-point number such that
        /// <c>1 + EPSILON</c> is numerically equal to 1. This value is an upper
        /// bound on the relative error due to rounding real numbers to double
        /// precision floating-point numbers.
        /// </para>
        /// <para>
        /// In IEEE 754 arithmetic, this is 2<sup>-53</sup>.
        /// </para>
        /// </summary>
        /// <remarks>
        /// See <a href="http://en.wikipedia.org/wiki/Machine_epsilon">Machine epsilon</a>
        /// </remarks>
        public static readonly double EPSILON;

        /// <summary>
        /// Safe minimum, such that <c>1 / SAFE_MIN</c> does not overflow.
        /// <para>
        /// In IEEE 754 arithmetic, this is also the smallest normalized
        /// number 2^-1022.
        /// </para>
        /// </summary>
        public static readonly Double SAFE_MIN = BitConverter.Int64BitsToDouble((EXPONENT_OFFSET - 1022L) << 52);

        /// <summary>
        /// Exponent offset in IEEE754 representation.
        /// </summary>
        private const Int64 EXPONENT_OFFSET = 1023L;

        /// <summary>
        /// Offset to order signed double numbers lexicographically.
        /// </summary>
        private const Int64 SGN_MASK = Int64.MinValue;

        /// <summary>
        /// Offset to order signed double numbers lexicographically.
        /// </summary>
        private const Int32 SGN_MASK_FLOAT = Int32.MinValue;

        /// <summary>
        /// Positive zero.
        /// </summary>
        private const Double POSITIVE_ZERO = 0d;

        /// <summary>
        /// Positive zero bits.
        /// </summary>
        private static readonly Int64 POSITIVE_ZERO_DOUBLE_BITS = BitConverter.DoubleToInt64Bits(+0.0);

        /// <summary>
        /// Negative zero bits.
        /// </summary>
        private static readonly Int64 NEGATIVE_ZERO_DOUBLE_BITS = BitConverter.DoubleToInt64Bits(-0.0);

        /// <summary>
        /// Positive zero bits.
        /// </summary>
        private static readonly Int32 POSITIVE_ZERO_FLOAT_BITS = (Int32)BitConverter.DoubleToInt64Bits(+0.0f);

        /// <summary>
        /// Negative zero bits.
        /// </summary>
        private static readonly Int32 NEGATIVE_ZERO_FLOAT_BITS = (Int32)BitConverter.DoubleToInt64Bits(-0.0f);


        /// <summary>
        /// Private constructor.
        /// </summary>
        private Precision() { }

        /// <summary>
        /// Compares two numbers given some amount of allowed error.
        /// </summary>
        /// <param name="x">the first number</param>
        /// <param name="y">the second number</param>
        /// <param name="eps"></param>
        /// <returns><list type="bullet">
        /// <item>0 if <see cref="equals(double, double, double)"/></item>
        /// <item>&lt; 0 if <see cref="equals(double, double, double)"/>
        /// &amp;&amp; x &lt; y</item>
        /// <item>&gt; 0 if <see cref="equals(double, double, double)"/>
        /// &amp;&amp; x &gt; y</item>
        /// </list></returns>
        public static Int32 compareTo(double x, double y, double eps)
        {
            if (equals(x, y, eps))
            {
                return 0;
            }
            else if (x < y)
            {
                return -1;
            }
            return 1;
        }

        /// <summary>
        /// Compares two numbers given some amount of allowed error.
        /// Two float numbers are considered equal if there are <c>(maxUlps - 1)</c>
        /// (or fewer) floating point numbers between them, i.e. two adjacent floating
        /// point numbers are considered equal.
        /// Adapted from <a
        /// href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
        /// Bruce Dawson</a>
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="maxUlps"><c>(maxUlps - 1)</c> is the number of floating point
        /// values between <c>x</c> and <c>y</c>.</param>
        /// <returns><list>
        /// <item>0 if  <see cref="equals(double, double, int)"/></item>
        /// <li>&lt; 0 if <see cref="equals(double, double, int)"/> &amp;&amp; x &lt; y</item>
        /// <li>&gt; 0 if <see cref="equals(double, double, int)"/> &amp;&amp; x &gt y</item>
        /// </list></returns>
        public static Int32 compareTo(double x, double y, int maxUlps)
        {
            if (equals(x, y, maxUlps))
            {
                return 0;
            }
            else if (x < y)
            {
                return -1;
            }
            return 1;
        }

        /// <summary>
        /// Returns true iff they are equal as defined by
        /// <see cref="equals(float,float,int)"/>.
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <returns><c>true</c> if the values are equal.</returns>
        public static Boolean equals(float x, float y)
        {
            return equals(x, y, 1);
        }

        /// <summary>
        /// Returns true if both arguments are NaN or neither is NaN and they are
        /// equal as defined by <see cref="#equals(float,float) equals(x, y, 1)"/>.
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <returns><c>true</c> if the values are equal or both are NaN.</returns>
        public static Boolean equalsIncludingNaN(float x, float y)
        {
            return (Single.IsNaN(x) || Single.IsNaN(y)) ? !(Single.IsNaN(x) ^ Single.IsNaN(y)) : equals(x, y, 1);
        }

        /// <summary>
        /// Returns true if both arguments are equal or within the range of allowed
        /// error (inclusive).
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="eps">the amount of absolute error to allow.</param>
        /// <returns><c>true</c> if the values are equal or within range of each other.
        /// </returns>
        public static Boolean equals(float x, float y, float eps)
        {
            return equals(x, y, 1) || FastMath.abs(y - x) <= eps;
        }

        /// <summary>
        /// Returns true if both arguments are NaN or are equal or within the range
        /// of allowed error (inclusive).
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="eps">the amount of absolute error to allow.</param>
        /// <returns><c>true</c> if the values are equal or within range of each other,
        /// or both are NaN.</returns>
        public static Boolean equalsIncludingNaN(float x, float y, float eps)
        {
            return equalsIncludingNaN(x, y) || (FastMath.abs(y - x) <= eps);
        }

        /// <summary>
        /// Returns true if both arguments are equal or within the range of allowed
        /// error (inclusive).
        /// Two float numbers are considered equal if there are <c>(maxUlps - 1)</c>
        /// (or fewer) floating point numbers between them, i.e. two adjacent floating
        /// point numbers are considered equal.
        /// Adapted from <a
        /// href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
        /// Bruce Dawson</a> 
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// values between <c>x</c> and <c>y</c>.</param>
        /// <returns><c>true</c> if there are fewer than <c>maxUlps</c> floating
        /// point values between <c>x</c> and <c>y</c>.</returns>
        public static Boolean equals(Single x, Single y, Int32 maxUlps)
        {

            int xInt = (Int32)BitConverter.DoubleToInt64Bits(x);
            int yInt = (Int32)BitConverter.DoubleToInt64Bits(y);

            Boolean isEqual;
            if (((xInt ^ yInt) & SGN_MASK_FLOAT) == 0)
            {
                // number have same sign, there is no risk of overflow
                isEqual = FastMath.abs(xInt - yInt) <= maxUlps;
            }
            else
            {
                // number have opposite signs, take care of overflow
                int deltaPlus;
                int deltaMinus;
                if (xInt < yInt)
                {
                    deltaPlus = yInt - POSITIVE_ZERO_FLOAT_BITS;
                    deltaMinus = xInt - NEGATIVE_ZERO_FLOAT_BITS;
                }
                else
                {
                    deltaPlus = xInt - POSITIVE_ZERO_FLOAT_BITS;
                    deltaMinus = yInt - NEGATIVE_ZERO_FLOAT_BITS;
                }

                if (deltaPlus > maxUlps)
                {
                    isEqual = false;
                }
                else
                {
                    isEqual = deltaMinus <= (maxUlps - deltaPlus);
                }

            }

            return isEqual && !Single.IsNaN(x) && !Single.IsNaN(y);

        }

        /// <summary>
        /// Returns true if both arguments are NaN or if they are equal as defined
        /// by <see cref="equals(float,float,int)"/>. 
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="maxUlps"><c>(maxUlps - 1)</c> is the number of floating point
        /// values between <c>x</c> and <c>y</c>.</param>
        /// <returns><c>true</c> if both arguments are NaN or if there are less than
        /// <c>maxUlps</c> floating point values between <c>x</c> and <c>y</c>.</returns>
        public static Boolean equalsIncludingNaN(float x, float y, int maxUlps)
        {
            return (Single.IsNaN(x) || Single.IsNaN(y)) ? !(Single.IsNaN(x) ^ Single.IsNaN(y)) : equals(x, y, maxUlps);
        }

        /// <summary>
        /// Returns true iff they are equal as defined by
        /// <see cref="equals(double,double,int)"/>. 
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <returns><c>true</c> if the values are equal.</returns>
        public static Boolean equals(double x, double y)
        {
            return equals(x, y, 1);
        }

        /// <summary>
        /// Returns true if both arguments are NaN or neither is NaN and they are
        /// equal as defined by <see cref="equals(double,double)"/>.
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <returns><c>true</c> if the values are equal or both are NaN.</returns>
        public static Boolean equalsIncludingNaN(double x, double y)
        {
            return (Double.IsNaN(x) || Double.IsNaN(y)) ? !(Double.IsNaN(x) ^ Double.IsNaN(y)) : equals(x, y, 1);
        }

        /// <summary>
        /// Returns <c>true</c> if there is no double value strictly between the
        /// arguments or the difference between them is within the range of allowed
        /// error (inclusive).
        /// </summary>
        /// <param name="x">First value.</param>
        /// <param name="y">Second value.</param>
        /// <param name="eps">Amount of allowed absolute error.</param>
        /// <returns><c>true</c> if the values are two adjacent floating point
        /// numbers or they are within range of each other.</returns>
        public static Boolean equals(double x, double y, double eps)
        {
            return equals(x, y, 1) || FastMath.abs(y - x) <= eps;
        }

        /// <summary>
        /// Returns <c>true</c> if there is no double value strictly between the
        /// arguments or the relative difference between them is smaller or equal
        /// to the given tolerance.
        /// </summary>
        /// <param name="x">First value.</param>
        /// <param name="y">Second value.</param>
        /// <param name="eps">Amount of allowed relative error.</param>
        /// <returns><c>true</c> if the values are two adjacent floating point
        /// numbers or they are within range of each other.</returns>
        public static Boolean equalsWithRelativeTolerance(double x, double y, double eps)
        {
            if (equals(x, y, 1))
            {
                return true;
            }

            double absoluteMax = FastMath.max(FastMath.abs(x), FastMath.abs(y));
            double relativeDifference = FastMath.abs((x - y) / absoluteMax);

            return relativeDifference <= eps;
        }

        /// <summary>
        /// Returns true if both arguments are NaN or are equal or within the range
        /// of allowed error (inclusive).
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="eps">the amount of absolute error to allow.</param>
        /// <returns><c>true</c> if the values are equal or within range of each other,
        /// or both are NaN.</returns>
        public static Boolean equalsIncludingNaN(double x, double y, double eps)
        {
            return equalsIncludingNaN(x, y) || (FastMath.abs(y - x) <= eps);
        }

        /// <summary>
        /// Returns true if both arguments are equal or within the range of allowed
        /// error (inclusive).
        /// <para>
        /// Two float numbers are considered equal if there are <c>(maxUlps - 1)</c>
        /// (or fewer) floating point numbers between them, i.e. two adjacent
        /// floating point numbers are considered equal.
        /// </para>
        /// <para>
        /// Adapted from <a
        /// href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
        /// Bruce Dawson</a>
        /// </para> 
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="maxUlps"><c>(maxUlps - 1)</c> is the number of floating point
        /// values between <c>x</> and <c>y</c>.</param>
        /// <returns><c>true</c> if there are fewer than <c>maxUlps</c> floating
        /// point values between <c>x</c> and <c>y</c>.</returns>
        public static Boolean equals(double x, double y, int maxUlps)
        {

            long xInt = BitConverter.DoubleToInt64Bits(x);
            long yInt = BitConverter.DoubleToInt64Bits(y);

            Boolean isEqual;
            if (((xInt ^ yInt) & SGN_MASK) == 0L)
            {
                // number have same sign, there is no risk of overflow
                isEqual = FastMath.abs(xInt - yInt) <= maxUlps;
            }
            else
            {
                // number have opposite signs, take care of overflow
                long deltaPlus;
                long deltaMinus;
                if (xInt < yInt)
                {
                    deltaPlus = yInt - POSITIVE_ZERO_DOUBLE_BITS;
                    deltaMinus = xInt - NEGATIVE_ZERO_DOUBLE_BITS;
                }
                else
                {
                    deltaPlus = xInt - POSITIVE_ZERO_DOUBLE_BITS;
                    deltaMinus = yInt - NEGATIVE_ZERO_DOUBLE_BITS;
                }

                if (deltaPlus > maxUlps)
                {
                    isEqual = false;
                }
                else
                {
                    isEqual = deltaMinus <= (maxUlps - deltaPlus);
                }

            }

            return isEqual && !Double.IsNaN(x) && !Double.IsNaN(y);

        }

        /// <summary>
        /// Returns true if both arguments are NaN or if they are equal as defined
        /// by <see cref="equals(double,double,int)"/>.
        /// </summary>
        /// <param name="x">first value</param>
        /// <param name="y">second value</param>
        /// <param name="maxUlps"><c>(maxUlps - 1)</c> is the number of floating point
        /// values between <c>x</c> and <c>y</c>.</param>
        /// <returns><c>true</c> if both arguments are NaN or if there are less than
        /// <c>maxUlps</c> floating point values between <c>x</c> and <c>y</c>.</returns>
        public static Boolean equalsIncludingNaN(double x, double y, int maxUlps)
        {
            return (Double.IsNaN(x) || Double.IsNaN(y)) ? !(Double.IsNaN(x) ^ Double.IsNaN(y)) : equals(x, y, maxUlps);
        }

        /// <summary>
        /// Rounds the given value to the specified number of decimal places.
        /// The value is rounded using the <see cref="BigDecimal.ROUND_HALF_UP"/> method.
        /// </summary>
        /// <param name="x">Value to round.</param>
        /// <param name="scale"> Number of digits to the right of the decimal point.</param>
        /// <returns>the rounded value.</returns>
        public static double round(double x, int scale)
        {
            return round(x, scale, BigDecimal.ROUND_HALF_UP);
        }

        /// <summary>
        /// Rounds the given value to the specified number of decimal places.
        /// The value is rounded using the given method which is any method defined
        /// in <see cref="BigDecimal"/>.
        /// If <c>x</c> is infinite or <c>NaN</c>, then the value of <c>x</c> is
        /// returned unchanged, regardless of the other parameters. 
        /// </summary>
        /// <param name="x">Value to round</param>
        /// <param name="scale">Number of digits to the right of the decimal point.</param>
        /// <param name="roundingMethod">Rounding method as defined in
        /// <see cref="BigDecimal"/>.</param>
        /// <returns>the rounded value.</returns>
        /// <exception cref="ArithmeticException">if <c>roundingMethod == ROUND_UNNECESSARY</c>
        /// and the specified scaling operation would require rounding.</exception>
        /// <exception cref="ArgumentOutOfRangeException">if <c>roundingMethod</c> does not
        /// represent a valid rounding mode.</exception>
        public static double round(double x, int scale, int roundingMethod)
        {
            try
            {

                double rounded = new BigDecimal(x)
                       .SetScale(scale, (Byte)roundingMethod)
                       .DoubleValue();
                // MATH-1089: negative values rounded to zero should result in negative zero
                return rounded == POSITIVE_ZERO ? POSITIVE_ZERO * x : rounded;
            }
            catch (FormatException)
            {
                if (Double.IsInfinity(x))
                {
                    return x;
                }
                else
                {
                    return Double.NaN;
                }
            }
        }

        /// <summary>
        /// Rounds the given value to the specified number of decimal places.
        /// The value is rounded using the <see cref="BigDecimal.ROUND_HALF_UP"/> method. 
        /// </summary>
        /// <param name="x">Value to round.</param>
        /// <param name="scale">Number of digits to the right of the decimal point.</param>
        /// <returns>the rounded value.</returns>
        public static float round(float x, int scale)
        {
            return round(x, scale, BigDecimal.ROUND_HALF_UP);
        }

        /// <summary>
        /// Rounds the given value to the specified number of decimal places.
        /// The value is rounded using the given method which is any method defined
        /// in <see cref="BigDecimal"/>.
        /// </summary>
        /// <param name="x">Value to round.</param>
        /// <param name="scale">Number of digits to the right of the decimal point</param>
        /// <param name="roundingMethod">Rounding method as defined in <see cref="BigDecimal"/>.
        /// </param>
        /// <returns>the rounded value.</returns>
        /// <exception cref="MathArithmeticException"> if an exact operation is required but 
        /// result is not exact</exception>
        /// <exception cref="MathIllegalArgumentException"> if <c>roundingMethod</c> is not a
        /// valid rounding method.</exception>
        public static float round(float x, int scale, int roundingMethod)
        {
            float sign = FastMath.copySign(1f, x);
            float factor = (float)FastMath.pow(10.0f, scale) * sign;
            return (float)roundUnscaled(x * factor, sign, roundingMethod) / factor;
        }

        /// <summary>
        /// Rounds the given non-negative value to the "nearest" integer. Nearest is
        /// determined by the rounding method specified. Rounding methods are defined
        /// in <see cref="BigDecimal"/>.
        /// </summary>
        /// <param name="unscaled">Value to round.</param>
        /// <param name="sign">Sign of the original, scaled value.</param>
        /// <param name="roundingMethod">Rounding method, as defined in <see cref="BigDecimal"/>.
        /// </param>
        /// <returns>the rounded value.</returns>
        /// <exception cref="MathArithmeticException"> if an exact operation is required but 
        /// result is not exact</exception>
        /// <exception cref="MathIllegalArgumentException"> if <c>roundingMethod</c> is not a
        /// valid rounding method.</exception>
        private static double roundUnscaled(double unscaled, double sign, int roundingMethod)
        {
            switch (roundingMethod)
            {
                case BigDecimal.ROUND_CEILING:
                    if (sign == -1)
                    {
                        unscaled = FastMath.floor(FastMath.nextAfter(unscaled, Double.NegativeInfinity));
                    }
                    else
                    {
                        unscaled = FastMath.ceil(FastMath.nextAfter(unscaled, Double.PositiveInfinity));
                    }
                    break;
                case BigDecimal.ROUND_DOWN:
                    unscaled = FastMath.floor(FastMath.nextAfter(unscaled, Double.NegativeInfinity));
                    break;
                case BigDecimal.ROUND_FLOOR:
                    if (sign == -1)
                    {
                        unscaled = FastMath.ceil(FastMath.nextAfter(unscaled, Double.PositiveInfinity));
                    }
                    else
                    {
                        unscaled = FastMath.floor(FastMath.nextAfter(unscaled, Double.NegativeInfinity));
                    }
                    break;
                case BigDecimal.ROUND_HALF_DOWN:
                    {
                        unscaled = FastMath.nextAfter(unscaled, Double.NegativeInfinity);
                        double fraction = unscaled - FastMath.floor(unscaled);
                        if (fraction > 0.5)
                        {
                            unscaled = FastMath.ceil(unscaled);
                        }
                        else
                        {
                            unscaled = FastMath.floor(unscaled);
                        }
                        break;
                    }
                case BigDecimal.ROUND_HALF_EVEN:
                    {
                        double fraction = unscaled - FastMath.floor(unscaled);
                        if (fraction > 0.5)
                        {
                            unscaled = FastMath.ceil(unscaled);
                        }
                        else if (fraction < 0.5)
                        {
                            unscaled = FastMath.floor(unscaled);
                        }
                        else
                        {
                            // The following equality test is intentional and needed for rounding purposes
                            if (FastMath.floor(unscaled) / 2.0 == FastMath.floor(FastMath.floor(unscaled) / 2.0))
                            { // even
                                unscaled = FastMath.floor(unscaled);
                            }
                            else
                            { // odd
                                unscaled = FastMath.ceil(unscaled);
                            }
                        }
                        break;
                    }
                case BigDecimal.ROUND_HALF_UP:
                    {
                        unscaled = FastMath.nextAfter(unscaled, Double.PositiveInfinity);
                        double fraction = unscaled - FastMath.floor(unscaled);
                        if (fraction >= 0.5)
                        {
                            unscaled = FastMath.ceil(unscaled);
                        }
                        else
                        {
                            unscaled = FastMath.floor(unscaled);
                        }
                        break;
                    }
                case BigDecimal.ROUND_UNNECESSARY:
                    if (unscaled != FastMath.floor(unscaled))
                    {
                        throw new MathArithmeticException();
                    }
                    break;
                case BigDecimal.ROUND_UP:
                    // do not round if the discarded fraction is equal to zero
                    if (unscaled != FastMath.floor(unscaled))
                    {
                        unscaled = FastMath.ceil(FastMath.nextAfter(unscaled, Double.PositiveInfinity));
                    }
                    break;
                default:
                    throw new MathIllegalArgumentException(new LocalizedFormats("INVALID_ROUNDING_METHOD"),
                                                           roundingMethod,
                                                           "ROUND_CEILING", BigDecimal.ROUND_CEILING,
                                                           "ROUND_DOWN", BigDecimal.ROUND_DOWN,
                                                           "ROUND_FLOOR", BigDecimal.ROUND_FLOOR,
                                                           "ROUND_HALF_DOWN", BigDecimal.ROUND_HALF_DOWN,
                                                           "ROUND_HALF_EVEN", BigDecimal.ROUND_HALF_EVEN,
                                                           "ROUND_HALF_UP", BigDecimal.ROUND_HALF_UP,
                                                           "ROUND_UNNECESSARY", BigDecimal.ROUND_UNNECESSARY,
                                                           "ROUND_UP", BigDecimal.ROUND_UP);
            }
            return unscaled;
        }

        /// <summary>
        /// Computes a number <c>delta</c> close to <c>originalDelta</c> with
        /// the property that <code>
        ///   x + delta - x
        /// </code>
        /// is exactly machine-representable.
        /// This is useful when computing numerical derivatives, in order to reduce
        /// roundoff errors.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="originalDelta">Offset value.</param>
        /// <returns>a number <c>delta</c> so that <c>x + delta</c> and <c>x</c>
        /// differ by a representable floating number.</returns>
        public static double representableDelta(double x,
                                                double originalDelta)
        {
            return x + originalDelta - x;
        }
    }
}
