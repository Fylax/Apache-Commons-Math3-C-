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
using System.Numerics;

namespace Math3.util
{
    /// <summary>
    /// Some useful, arithmetics related, additions to the built-in functions in
    /// <see cref="System.Math"/>.
    /// </summary>
    public class ArithmeticUtils
    {
        /// <summary>
        /// Private constructor.
        /// </summary>
        private ArithmeticUtils() { }

        /// <summary>
        /// Add two integers, checking for overflow.
        /// </summary>
        /// <param name="x">an addend</param>
        /// <param name="y">an addend</param>
        /// <returns>the sum <c>x+y</c></returns>
        /// <exception cref=">MathArithmeticException" if the result can not be represented
        /// as an <c>int</c>.</exception>
        public static int addAndCheck(int x, int y)
        {
            long s = (long)x + (long)y;
            if (s < Int32.MinValue || s > Int32.MaxValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("OVERFLOW_IN_ADDITION"), x, y);
            }
            return (int)s;
        }

        /// <summary>
        /// Add two long integers, checking for overflow.
        /// </summary>
        /// <param name="a">an addend</param>
        /// <param name="b">an addend</param>
        /// <returns>the sum <c>a+b</c></returns>
        /// <exception cref="MathArithmeticException"> if the result can not be represented
        /// as an long</exception>
        public static long addAndCheck(long a, long b)
        {
            return addAndCheck(a, b, new LocalizedFormats("OVERFLOW_IN_ADDITION"));
        }

        /// <summary>
        /// Returns an exact representation of the <a
        /// href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
        /// Coefficient</a>, "<c>n choose k</c>", the number of
        /// <c>k</c>-element subsets that can be selected from an
        /// <c>n</c>-element set.
        /// <para>
        /// Preconditions:
        /// <list type="bullet">
        /// <item> <c>0 <= k <= n</c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// <item> The result is small enough to fit into a <c>long</c>. The
        /// largest value of <c>n</c> for which all coefficients are
        /// <c> < Long.MAX_VALUE</c> is 66. If the computed value exceeds
        /// <c>Int64.MaxValue</c> an <c>ArithMeticException</c> is
        /// thrown.</item>
        /// </list></para>
        /// </summary>
        /// <param name="n">the size of the set</param>
        /// <param name="k">the size of the subsets to be counted</param>
        /// <returns><c>n choose k</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result is too large to be
        /// represented by a long integer.</exception>
        [Obsolete("use CombinatoricsUtils.binomialCoefficient(int, int)")]
        public static long binomialCoefficient(int n, int k)
        {
            return CombinatoricsUtils.binomialCoefficient(n, k);
        }

        /// <summary>
        /// Returns a <c>double</c> representation of the <a
        /// href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
        /// Coefficient</a>, "<c>n choose k</c>", the number of
        /// <c>k</c>-element subsets that can be selected from an
        /// <c>n</c>-element set.
        /// <para>
        /// Preconditions:
        /// <list type="bullet">
        /// <item> <c>0 <= k <= n </c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// <item> The result is small enough to fit into a <c>double</c>. The
        /// largest value of <c>n</c> for which all coefficients are <
        /// Double.MaxValue is 1029. If the computed value exceeds Double.MaxValue,
        /// Double.PositiveInfinity is returned</item>
        /// </list></para>
        /// </summary>
        /// <param name="n">the size of the set</param>
        /// <param name="k">the size of the subsets to be counted</param>
        /// <returns><c>n choose k</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result is too large to be
        /// represented by a long integer.</exception>
        [Obsolete("use CombinatoricsUtils.binomialCoefficientDouble(int, int)")]
        public static double binomialCoefficientDouble(int n, int k)
        {
            return CombinatoricsUtils.binomialCoefficientDouble(n, k);
        }

        /// <summary>
        /// Returns the natural <c>log</c> of the <a
        /// href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
        /// Coefficient</a>, "<c>n choose k</c>", the number of
        /// <c>k</c>-element subsets that can be selected from an
        /// <c>n</c>-element set.
        /// <para>
        /// Preconditions:
        /// <list type="bullet">
        /// <item> <c> 0 <= k <= n </c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// </list></para>
        /// </summary>
        /// <param name="n">the size of the set</param>
        /// <param name="k">the size of the subsets to be counted</param>
        /// <returns><<c>n choose k</c>/returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result is too large to be
        /// represented by a long integer.</exception>
        [Obsolete("use CombinatoricsUtils.binomialCoefficientLong(int, int)")]
        public static double binomialCoefficientLog(int n, int k)
        {
            return CombinatoricsUtils.binomialCoefficientLog(n, k);
        }

        /// <summary>
        /// Returns n!. Shorthand for <c>n</c> <a
        /// href="http://mathworld.wolfram.com/Factorial.html"> Factorial</a>, the
        /// product of the numbers <c>1,...,n</c>.
        /// <para>
        /// Preconditions:
        /// <list type="bullet">
        /// <item> <c> n >= 0</c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// <item> The result is small enough to fit into a <c>long</c>. The
        /// largest value of <c>n</c> for which <c>n! <
        /// Int64.MaxValue</c> is 20. If the computed value exceeds <c>Int64.MaxValue</c>
        /// an <c>ArithMeticException</c> is thrown.</item>
        /// </list>
        /// </para>
        /// </summary>
        /// <param name="n">argument</param>
        /// <returns><c>n!</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result is too large to be
        /// represented by a long integer.</exception>
        [Obsolete("use CombinatoricsUtils.factorial(int)")]
        public static long factorial(int n)
        {
            return CombinatoricsUtils.factorial(n);
        }

        /// <summary>
        /// Compute n!, the<a href="http://mathworld.wolfram.com/Factorial.html">
        /// factorial</a> of <c>n</c> (the product of the numbers 1 to n), as a
        /// <c>double</c>.
        /// The result should be small enough to fit into a <c>double</c>: The
        /// largest <c>n</c> for which <c>n! < Double.MaxValue</c> is 170.
        /// If the computed value exceeds <c>Double.MaxValue</c>,
        /// <c>Double.PositiveInfinity</c> is returned.
        /// </summary>
        /// <param name="n">Argument.</param>
        /// <returns><c>n!</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c></exception>
        [Obsolete("use CombinatoricsUtils.factorialDouble(int)")]
        public static double factorialDouble(int n)
        {
            return CombinatoricsUtils.factorialDouble(n);
        }

        /// <summary>
        /// Compute the natural logarithm of the factorial of <c>n</c>.
        /// </summary>
        /// <param name="n">Argument.</param>
        /// <returns><c>n!</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c></exception>
        [Obsolete("use CombinatoricsUtils.factorialLog(int)")]
        public static double factorialLog(int n)
        {
            return CombinatoricsUtils.factorialLog(n);
        }

        /// <summary>
        /// Computes the greatest common divisor of the absolute value of two
        /// numbers, using a modified version of the "binary gcd" method.
        /// See Knuth 4.5.2 algorithm B.
        /// The algorithm is due to Josef Stein (1961).
        /// <para/>
        /// Special cases:
        /// <list type="bullet">
        /// <item>The invocations
        /// <c>gcd(Int32.MinValue, Int32.MinValue)</c>,
        /// <c>gcd(Int32.MinValue, 0)</c> and
        /// <c>gcd(0, Int32.MinValue)</c> throw an
        /// <c>ArithmeticException</c>, because the result would be 2^31, which
        /// is too large for an int value.</item>
        /// <item>The result of <c>gcd(x, x)</c>, <c>gcd(0, x)</c> and
        /// <c>gcd(x, 0)</c> is the absolute value of <c>x</c>, except
        /// for the special cases above.</item>
        /// <item>The invocation <c>gcd(0, 0)</c> is the only one which returns
        /// <c>0</c>.</item>
        /// </list>
        /// </summary>
        /// <param name="p">Number.</param>
        /// <param name="q">Number.</param>
        /// <returns>the greatest common divisor (never negative).</returns>
        /// <exception cref="MathArithmeticException"> if the result cannot be represented as
        /// a non-negative <c>int</c> value.</exception>
        public static int gcd(int p, int q)
        {
            int a = p;
            int b = q;
            if (a == 0 ||
                b == 0)
            {
                if (a == Int32.MinValue ||
                    b == Int32.MinValue)
                {
                    throw new MathArithmeticException(new LocalizedFormats("GCD_OVERFLOW_32_BITS"), p, q);
                }
                return FastMath.abs(a + b);
            }

            long al = a;
            long bl = b;
            Boolean useLong = false;
            if (a < 0)
            {
                if (Int32.MinValue == a)
                {
                    useLong = true;
                }
                else
                {
                    a = -a;
                }
                al = -al;
            }
            if (b < 0)
            {
                if (Int32.MinValue == b)
                {
                    useLong = true;
                }
                else
                {
                    b = -b;
                }
                bl = -bl;
            }
            if (useLong)
            {
                if (al == bl)
                {
                    throw new MathArithmeticException(new LocalizedFormats("GCD_OVERFLOW_32_BITS"), p, q);
                }
                long blbu = bl;
                bl = al;
                al = blbu % al;
                if (al == 0)
                {
                    if (bl > Int32.MaxValue)
                    {
                        throw new MathArithmeticException(new LocalizedFormats("GCD_OVERFLOW_32_BITS"), p, q);
                    }
                    return (int)bl;
                }
                blbu = bl;

                // Now "al" and "bl" fit in an "int".
                b = (int)al;
                a = (int)(blbu % al);
            }

            return gcdPositive(a, b);
        }

        /// <summary>
        /// Computes the greatest common divisor of two positive numbers
        /// (this precondition is not checked and the result is undefined
        /// if not fulfilled) using the "binary gcd" method which avoids division
        /// and modulo operations.
        /// See Knuth 4.5.2 algorithm B.
        /// The algorithm is due to Josef Stein (1961).
        /// <para/>
        /// Special cases:
        /// <list type="bullet">
        /// <item>The result of <c>gcd(x, x)</c>, <c>gcd(0, x)</c> and
        /// <c>gcd(x, 0)</c> is the value of <c>x</c>.</item>
        /// <item>The invocation <c>gcd(0, 0)</c> is the only one which returns
        /// <c>0.</item>
        /// </list>
        /// </summary>
        /// <param name="a">Positive number.</param>
        /// <param name="b">Positive number.</param>
        /// <returns>the greatest common divisor.</returns>
        private static int gcdPositive(int a, int b)
        {
            if (a == 0)
            {
                return b;
            }
            else if (b == 0)
            {
                return a;
            }

            // Make "a" and "b" odd, keeping track of common power of 2.
            int mask = 1;
            int aTwos = 32;
            for (int i = 0; i < 32; ++i, mask <<= 1) //number of trailing zeros
            {
                if ((a & mask) != 0)
                {
                    aTwos = i;
                }
            }
            a >>= aTwos;
            mask = 1;
            int bTwos = 32;
            for (int i = 0; i < 32; ++i, mask <<= 1) //number of trailing zeros
            {
                if ((b & mask) != 0)
                {
                    bTwos = i;
                }
            }
            b >>= bTwos;
            int shift = FastMath.min(aTwos, bTwos);

            // "a" and "b" are positive.
            // If a > b then "gdc(a, b)" is equal to "gcd(a - b, b)".
            // If a < b then "gcd(a, b)" is equal to "gcd(b - a, a)".
            // Hence, in the successive iterations:
            //  "a" becomes the absolute difference of the current values,
            //  "b" becomes the minimum of the current values.
            while (a != b)
            {
                int delta = a - b;
                b = Math.Min(a, b);
                a = Math.Abs(delta);

                // Remove any power of 2 in "a" ("b" is guaranteed to be odd).
                mask = 1;
                int aTZ = 32;
                for (int i = 0; i < 32; ++i, mask <<= 1) //number of trailing zeros
                {
                    if ((a & mask) != 0)
                    {
                        aTZ = i;
                    }
                }
                a >>= aTZ;
            }

            // Recover the common power of 2.
            return a << shift;
        }

        /// <summary>
        /// <para>
        /// Gets the greatest common divisor of the absolute value of two numbers,
        /// using the "binary gcd" method which avoids division and modulo
        /// operations. See Knuth 4.5.2 algorithm B. This algorithm is due to Josef
        /// Stein (1961).
        /// </para>
        /// Special cases:
        /// <list type="bullet">
        /// <item>The invocations
        /// <c>gcd(Int64.MinValue, Int64.MinValue)</c>,
        /// <c>gcd(Int64.MinValue, 0L)</c> and
        /// <c>gcd(0L, Int64.MinValue)</c> throw an
        /// <c>ArithmeticException</c>, because the result would be 2^63, which
        /// is too large for a long value.</item>
        /// <item>The result of <c>gcd(x, x)</c>, <c>gcd(0L, x)</c> and
        /// <c>gcd(x, 0L)</c> is the absolute value of <c>x</c>, except
        /// for the special cases above.</item>
        /// <item>The invocation <c>gcd(0L, 0L)</c> is the only one which returns
        /// <c>0L</c>.</item>
        /// </list>
        /// </summary>
        /// <param name="p">Number.</param>
        /// <param name="q">Number.</param>
        /// <returns>the greatest common divisor, never negative.</returns>
        /// <exception cref="MathArithmeticException"> if the result cannot be represented as
        /// a non-negative <c>long</c> value.</exception>
        public static long gcd(long p, long q)
        {
            long u = p;
            long v = q;
            if ((u == 0) || (v == 0))
            {
                if ((u == Int64.MinValue) || (v == Int64.MinValue))
                {
                    throw new MathArithmeticException(new LocalizedFormats("GCD_OVERFLOW_64_BITS"), p, q);
                }
                return FastMath.abs(u) + FastMath.abs(v);
            }
            // keep u and v negative, as negative integers range down to
            // -2^63, while positive numbers can only be as large as 2^63-1
            // (i.e. we can't necessarily negate a negative number without
            // overflow)
            /* assert u!=0 && v!=0; */
            if (u > 0)
            {
                u = -u;
            } // make u negative
            if (v > 0)
            {
                v = -v;
            } // make v negative
            // B1. [Find power of 2]
            int k = 0;
            while ((u & 1) == 0 && (v & 1) == 0 && k < 63)
            { // while u and v are
                // both even...
                u /= 2;
                v /= 2;
                k++; // cast out twos.
            }
            if (k == 63)
            {
                throw new MathArithmeticException(new LocalizedFormats("GCD_OVERFLOW_64_BITS"), p, q);
            }
            // B2. Initialize: u and v have been divided by 2^k and at least
            // one is odd.
            long t = ((u & 1) == 1) ? v : -(u / 2)/* B3 */;
            // t negative: u was odd, v may be even (t replaces v)
            // t positive: u was even, v is odd (t replaces u)
            do
            {
                /* assert u<0 && v<0; */
                // B4/B3: cast out twos from t.
                while ((t & 1) == 0)
                { // while t is even..
                    t /= 2; // cast out twos
                }
                // B5 [reset max(u,v)]
                if (t > 0)
                {
                    u = -t;
                }
                else
                {
                    v = t;
                }
                // B6/B3. at this point both u and v should be odd.
                t = (v - u) / 2;
                // |u| larger: t positive (replace u)
                // |v| larger: t negative (replace v)
            } while (t != 0);
            return -u * (1L << k); // gcd is u*2^k
        }

        /// <summary>
        /// <para>
        /// Returns the least common multiple of the absolute value of two numbers,
        /// using the formula <c>lcm(a,b) = (a / gcd(a,b)) * b</c>.
        /// </para>
        /// Special cases:
        /// <list type="bullet">
        /// <item>The invocations <c>lcm(Int32.MinValue, n)</c> and
        /// <c>lcm(n, Int32.MinValue)</c>, where <c>abs(n)</c> is a
        /// power of 2, throw an <c>ArithmeticException</c>, because the result
        /// would be 2^31, which is too large for an int value.</item>
        /// <item>The result of <c>lcm(0, x)</c> and <c>lcm(x, 0)</c> is
        /// <c>0</c> for any <c>x</c>.</item>
        /// </list>
        /// </summary>
        /// <param name="a">Number.</param>
        /// <param name="b">Number.</param>
        /// <returns>the least common multiple, never negative.</returns>
        /// <exception cref="MathArithmeticException"> if the result cannot be represented as
        /// a non-negative <c>int</c> value.</exception>
        public static int lcm(int a, int b)
        {
            if (a == 0 || b == 0)
            {
                return 0;
            }
            int lcm = FastMath.abs(ArithmeticUtils.mulAndCheck(a / gcd(a, b), b));
            if (lcm == Int32.MinValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("LCM_OVERFLOW_32_BITS"), a, b);
            }
            return lcm;
        }

        /// <summary>
        /// <para>
        /// Returns the least common multiple of the absolute value of two numbers,
        /// using the formula <c>lcm(a,b) = (a / gcd(a,b)) * b</c>.
        /// </para>
        /// Special cases:
        /// <list type="bullet">
        /// <item>The invocations <c>lcm(Int64.MinValue, n)</c> and
        /// <c>lcm(n, Int64.MinValue)</c>, where <c>abs(n)</c> is a
        /// power of 2, throw an <c>ArithmeticException</c>, because the result
        /// would be 2^63, which is too large for an int value.</item>
        /// <item>The result of <c>lcm(0L, x)</c> and <c>lcm(x, 0L)</c> is
        /// <c>0L</c> for any <c>x</c>.</item>
        /// </list>
        /// </summary>
        /// <param name="a">Number.</param>
        /// <param name="b">Number.</param>
        /// <returns>the least common multiple, never negative.</returns>
        /// <exception cref="MathArithmeticException"> if the result cannot be represented
        /// as a non-negative <c>long</c> value.</exception>
        public static long lcm(long a, long b)
        {
            if (a == 0 || b == 0)
            {
                return 0;
            }
            long lcm = FastMath.abs(ArithmeticUtils.mulAndCheck(a / gcd(a, b), b));
            if (lcm == Int64.MinValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("LCM_OVERFLOW_64_BITS"), a, b);
            }
            return lcm;
        }

        /// <summary>
        /// Multiply two integers, checking for overflow.
        /// </summary>
        /// <param name="x">Factor.</param>
        /// <param name="y">Factor.</param>
        /// <returns>the product <c>x * y</c>.</returns>
        /// <exception cref="MathArithmeticException"> if the result can not be
        /// represented as an <c>int</c>.</exception>
        public static int mulAndCheck(int x, int y)
        {
            long m = ((long)x) * ((long)y);
            if (m < Int32.MinValue || m > Int32.MaxValue)
            {
                throw new MathArithmeticException();
            }
            return (int)m;
        }

        /// <summary>
        /// Multiply two long integers, checking for overflow.
        /// </summary>
        /// <param name="a">Factor.</param>
        /// <param name="b">Factor.</param>
        /// <returns>the product <c>a * b</c>.</returns>
        /// <exception cref="MathArithmeticException"> if the result can not be represented
        /// as a <c>long</c>.</exception>
        public static long mulAndCheck(long a, long b)
        {
            long ret;
            if (a > b)
            {
                // use symmetry to reduce boundary cases
                ret = mulAndCheck(b, a);
            }
            else
            {
                if (a < 0)
                {
                    if (b < 0)
                    {
                        // check for positive overflow with negative a, negative b
                        if (a >= Int64.MaxValue / b)
                        {
                            ret = a * b;
                        }
                        else
                        {
                            throw new MathArithmeticException();
                        }
                    }
                    else if (b > 0)
                    {
                        // check for negative overflow with negative a, positive b
                        if (Int64.MinValue / b <= a)
                        {
                            ret = a * b;
                        }
                        else
                        {
                            throw new MathArithmeticException();
                        }
                    }
                    else
                    {
                        // assert b == 0
                        ret = 0;
                    }
                }
                else if (a > 0)
                {
                    // assert a > 0
                    // assert b > 0

                    // check for positive overflow with positive a, positive b
                    if (a <= Int64.MaxValue / b)
                    {
                        ret = a * b;
                    }
                    else
                    {
                        throw new MathArithmeticException();
                    }
                }
                else
                {
                    // assert a == 0
                    ret = 0;
                }
            }
            return ret;
        }

        /// <summary>
        /// Subtract two integers, checking for overflow.
        /// </summary>
        /// <param name="x">Minuend.</param>
        /// <param name="y">Subtrahend.</param>
        /// <returns>the difference <c>x - y</c>.</returns>
        /// <exception cref="MathArithmeticException"> if the result can not be represented
        /// as an <c>int</c>.</exception>
        public static int subAndCheck(int x, int y)
        {
            long s = (long)x - (long)y;
            if (s < Int32.MinValue || s > Int32.MaxValue)
            {
                throw new MathArithmeticException(new LocalizedFormats("OVERFLOW_IN_SUBTRACTION"), x, y);
            }
            return (int)s;
        }

        /// <summary>
        /// Subtract two long integers, checking for overflow.
        /// </summary>
        /// <param name="x">Minuend.</param>
        /// <param name="y">Subtrahend.</param>
        /// <returns>the difference <c>x - y</c>.</returns>
        /// <exception cref="MathArithmeticException"> if the result can not be represented
        /// as a <c>long</c>.</exception>
        public static long subAndCheck(long a, long b)
        {
            long ret;
            if (b == Int64.MinValue)
            {
                if (a < 0)
                {
                    ret = a - b;
                }
                else
                {
                    throw new MathArithmeticException(new LocalizedFormats("OVERFLOW_IN_ADDITION"), a, -b);
                }
            }
            else
            {
                // use additive inverse
                ret = addAndCheck(a, -b, new LocalizedFormats("OVERFLOW_IN_ADDITION"));
            }
            return ret;
        }

        /// <summary>
        /// Raise an int to an int power.
        /// </summary>
        /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result would overflow.</exception>
        public static int pow(int k, int e)
        {
            if (e < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("EXPONENT"), e);
            }

            try
            {
                int exp = e;
                int result = 1;
                int k2p = k;
                while (true)
                {
                    if ((exp & 0x1) != 0)
                    {
                        result = mulAndCheck(result, k2p);
                    }

                    exp >>= 1;
                    if (exp == 0)
                    {
                        break;
                    }

                    k2p = mulAndCheck(k2p, k2p);
                }

                return result;
            }
            catch (MathArithmeticException mae)
            {
                // Add context information.
                mae.getContext().addMessage(new LocalizedFormats("OVERFLOW"));
                mae.getContext().addMessage(new LocalizedFormats("BASE"), k);
                mae.getContext().addMessage(new LocalizedFormats("EXPONENT"), e);

                // Rethrow.
                throw;
            }
        }

        /// <summary>
        /// Raise an int to a long power.
        /// </summary>
        /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        [Obsolete("Please use pow(int,int) instead.")]
        public static int pow(int k, long e)
        {
            if (e < 0)
            {
                throw new NotPositiveException<Int64>(new LocalizedFormats("EXPONENT"), e);
            }

            int result = 1;
            int k2p = k;
            while (e != 0)
            {
                if ((e & 0x1) != 0)
                {
                    result *= k2p;
                }
                k2p *= k2p;
                e >>= 1;
            }

            return result;
        }

        /// <summary>
        /// Raise a long to an int power.
        /// </summary>
        /// /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result would overflow.</exception>
        public static long pow(long k, int e)
        {
            if (e < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("EXPONENT"), e);
            }

            try
            {
                int exp = e;
                long result = 1;
                long k2p = k;
                while (true)
                {
                    if ((exp & 0x1) != 0)
                    {
                        result = mulAndCheck(result, k2p);
                    }

                    exp >>= 1;
                    if (exp == 0)
                    {
                        break;
                    }

                    k2p = mulAndCheck(k2p, k2p);
                }

                return result;
            }
            catch (MathArithmeticException mae)
            {
                // Add context information.
                mae.getContext().addMessage(new LocalizedFormats("OVERFLOW"));
                mae.getContext().addMessage(new LocalizedFormats("BASE"), k);
                mae.getContext().addMessage(new LocalizedFormats("EXPONENT"), e);

                // Rethrow.
                throw;
            }
        }

        /// <summary>
        /// Raise a long to a long power.
        /// </summary>
        /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        [Obsolete("Please use pow(long,int) instead.")]
        public static long pow(long k, long e)
        {
            if (e < 0)
            {
                throw new NotPositiveException<Int64>(new LocalizedFormats("EXPONENT"), e);
            }
            long result = 1L;
            long k2p = k;
            while (e != 0)
            {
                if ((e & 0x1) != 0)
                {
                    result *= k2p;
                }
                k2p *= k2p;
                e >>= 1;
            }
            return result;
        }

        /// <summary>
        /// Raise an BigInteger to an int power.
        /// </summary>
        /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        public static BigInteger pow(BigInteger k, int e)
        {
            if (e < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("EXPONENT"), e);
            }
            return BigInteger.Pow(k, e);
        }

        /// <summary>
        /// Raise an BigInteger to an long power.
        /// </summary>
        /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        public static BigInteger pow(BigInteger k, long e)
        {
            if (e < 0)
            {
                throw new NotPositiveException<Int64>(new LocalizedFormats("EXPONENT"), e);
            }

            BigInteger result = BigInteger.One;
            BigInteger k2p = k;
            while (e != 0)
            {
                if ((e & 0x1) != 0)
                {
                    result = BigInteger.Multiply(result, k2p);
                }
                k2p = BigInteger.Multiply(k2p, k2p);
                e >>= 1;
            }
            return result;
        }

        /// <summary>
        /// Raise an BigInteger to a BigInteger power.
        /// </summary>
        /// <param name="k">Number to raise.</param>
        /// <param name="e">Exponent (must be positive or zero).</param>
        /// <returns>\( k^e \)</returns>
        /// <exception cref="NotPositiveException"> if <c>e < 0</c>.</exception>
        public static BigInteger pow(BigInteger k, BigInteger e)
        {
            if (e.CompareTo(BigInteger.Zero) < 0)
            {
                throw new NotPositiveException<BigInteger>(new LocalizedFormats("EXPONENT"), e);
            }

            BigInteger result = BigInteger.One;
            BigInteger k2p = k;
            while (!BigInteger.Zero.Equals(e))
            {
                if ((e & (1 << 0)) != 0) //equivalent of BigInteger.testBit(0) 
                {
                    result = BigInteger.Multiply(result, k2p);
                }
                k2p = BigInteger.Multiply(k2p, k2p);
                e >>= 1;
            }
            return result;
        }

        /// <summary>
        /// Returns the <a
        /// href="http://mathworld.wolfram.com/StirlingNumberoftheSecondKind.html">
        /// Stirling number of the second kind</a>, "<c>S(n,k)</c>", the number of
        /// ways of partitioning an <c>n</c>-element set into <c>k</c> non-empty
        /// subsets.
        /// <para>
        /// The preconditions are <c>0 <= k <= n</c> (otherwise
        /// <c>NotPositiveException</c> is thrown)
        /// </para>
        /// </summary>
        /// <param name="n">the size of the set</param>
        /// <param name="k">the number of non-empty subsets</param>
        /// <returns><c>S(n,k)</c></returns>
        /// <exception cref="NotPositiveException"> if <c>k < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if some overflow happens, typically for
        /// n exceeding 25 and k between 20 and n-2 (S(n,n-1) is handled specifically and does 
        /// not overflow)</exception>
        [Obsolete("CombinatoricsUtils.stirlingS2(int, int)")]
        public static long stirlingS2(int n, int k)
        {
            return CombinatoricsUtils.stirlingS2(n, k);
        }

        /// <summary>
        /// Add two long integers, checking for overflow.
        /// </summary>
        /// <param name="a">Addend.</param>
        /// <param name="b">Addend.</param>
        /// <param name="pattern">Pattern to use for any thrown exception.</param>
        /// <returns>the sum <c>a + b</c>.</returns>
        /// <exception cref="MathArithmeticException"> if the result cannot be represented
        /// as a <c>long</c>.</exception>
        private static long addAndCheck(long a, long b, Localizable pattern)
        {
            long result = a + b;
            if (!((a ^ b) < 0 | (a ^ result) >= 0))
            {
                throw new MathArithmeticException(pattern, a, b);
            }
            return result;
        }

        /// <summary>
        /// Returns true if the argument is a power of two.
        /// </summary>
        /// <param name="n">the number to test</param>
        /// <returns>true if the argument is a power of two</returns>
        public static Boolean isPowerOfTwo(long n)
        {
            return (n > 0) && ((n & (n - 1)) == 0);
        }
    }
}
