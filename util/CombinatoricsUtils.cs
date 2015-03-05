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
using System.Collections.Generic;

namespace Math3.util
{
    /// <summary>
    /// Combinatorial utilities.
    /// </summary>
    public class CombinatoricsUtils
    {
        /// <summary>
        /// All long-representable factorials
        /// </summary>
        internal static readonly long[] FACTORIALS = new long[]
        {
                           1L,                  1L,                   2L,
                           6L,                 24L,                 120L,
                         720L,               5040L,               40320L,
                      362880L,            3628800L,            39916800L,
                   479001600L,         6227020800L,         87178291200L,
               1307674368000L,     20922789888000L,     355687428096000L,
            6402373705728000L, 121645100408832000L, 2432902008176640000L 
        };

        /// <summary>
        /// Stirling numbers of the second kind.
        /// </summary>
        internal static long[][] STIRLING_S2 = null;

        /// <summary>
        /// Private constructor (class contains only static methods).
        /// </summary>
        private CombinatoricsUtils() { }

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
        /// <c>MathIllegalArgumentException</c> is thrown)</item>
        /// <item> The result is small enough to fit into a <c>long</c>. The
        /// largest value of <c>n</c> for which all coefficients are
        /// <c> < Int64.MaxValue</c> is 66. If the computed value exceeds
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
        public static long binomialCoefficient(int n, int k)
        {
            CombinatoricsUtils.checkBinomial(n, k);
            if ((n == k) || (k == 0))
            {
                return 1;
            }
            if ((k == 1) || (k == n - 1))
            {
                return n;
            }
            // Use symmetry for large k
            if (k > n / 2)
            {
                return binomialCoefficient(n, n - k);
            }

            // We use the formula
            // (n choose k) = n! / (n-k)! / k!
            // (n choose k) == ((n-k+1)*...*n) / (1*...*k)
            // which could be written
            // (n choose k) == (n-1 choose k-1) * n / k
            long result = 1;
            if (n <= 61)
            {
                // For n <= 61, the naive implementation cannot overflow.
                int i = n - k + 1;
                for (int j = 1; j <= k; j++)
                {
                    result = result * i / j;
                    i++;
                }
            }
            else if (n <= 66)
            {
                // For n > 61 but n <= 66, the result cannot overflow,
                // but we must take care not to overflow intermediate values.
                int i = n - k + 1;
                for (int j = 1; j <= k; j++)
                {
                    // We know that (result * i) is divisible by j,
                    // but (result * i) may overflow, so we split j:
                    // Filter out the gcd, d, so j/d and i/d are integer.
                    // result is divisible by (j/d) because (j/d)
                    // is relative prime to (i/d) and is a divisor of
                    // result * (i/d).
                    long d = ArithmeticUtils.gcd(i, j);
                    result = (result / (j / d)) * (i / d);
                    i++;
                }
            }
            else
            {
                // For n > 66, a result overflow might occur, so we check
                // the multiplication, taking care to not overflow
                // unnecessary.
                int i = n - k + 1;
                for (int j = 1; j <= k; j++)
                {
                    long d = ArithmeticUtils.gcd(i, j);
                    result = ArithmeticUtils.mulAndCheck(result / (j / d), i / d);
                    i++;
                }
            }
            return result;
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
        /// <item> <c>0 <= k <= n</c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// <item> The result is small enough to fit into a <c>double</c>. The
        /// largest value of <c>n</c> for which all coefficients are <
        /// Double.MaxValue is 1029. If the computed value exceeds Double.MaxValue,
        /// Double.PositiveInvinifty is returned</item>
        /// </list></para>
        /// </summary>
        /// <param name="n">the size of the set</param>
        /// <param name="k">the size of the subsets to be counted</param>
        /// <returns><c>n choose k</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result is too large to be
        /// represented by a double.</exception>
        public static double binomialCoefficientDouble(int n, int k)
        {
            CombinatoricsUtils.checkBinomial(n, k);
            if ((n == k) || (k == 0))
            {
                return 1d;
            }
            if ((k == 1) || (k == n - 1))
            {
                return n;
            }
            if (k > n / 2)
            {
                return binomialCoefficientDouble(n, n - k);
            }
            if (n < 67)
            {
                return binomialCoefficient(n, k);
            }

            double result = 1d;
            for (int i = 1; i <= k; i++)
            {
                result *= (double)(n - k + i) / (double)i;
            }

            return FastMath.floor(result + 0.5);
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
        /// <item> <c>0 <= k <= n</c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// </list></para>
        /// </summary>
        /// <param name="n">the size of the set</param>
        /// <param name="k">the size of the subsets to be counted</param>
        /// <returns><c>n choose k</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        /// <exception cref="MathArithmeticException"> if the result is too large to be
        /// represented by a double.</exception>
        public static double binomialCoefficientLog(int n, int k)
        {
            CombinatoricsUtils.checkBinomial(n, k);
            if ((n == k) || (k == 0))
            {
                return 0;
            }
            if ((k == 1) || (k == n - 1))
            {
                return FastMath.log(n);
            }

            /*
             * For values small enough to do exact integer computation,
             * return the log of the exact value
             */
            if (n < 67)
            {
                return FastMath.log(binomialCoefficient(n, k));
            }

            /*
             * Return the log of binomialCoefficientDouble for values that will not
             * overflow binomialCoefficientDouble
             */
            if (n < 1030)
            {
                return FastMath.log(binomialCoefficientDouble(n, k));
            }

            if (k > n / 2)
            {
                return binomialCoefficientLog(n, n - k);
            }

            /*
             * Sum logs for values that could overflow
             */
            double logSum = 0;

            // n!/(n-k)!
            for (int i = n - k + 1; i <= n; i++)
            {
                logSum += FastMath.log(i);
            }

            // divide by k!
            for (int i = 2; i <= k; i++)
            {
                logSum -= FastMath.log(i);
            }

            return logSum;
        }

        /// <summary>
        /// Returns n!. Shorthand for <c>n</c> <a
        /// href="http://mathworld.wolfram.com/Factorial.html">Factorial</a>, the
        /// product of the numbers <c>1,...,n</c>.
        /// <para>
        /// Preconditions:
        /// <list type="bullet">
        /// <item> <c>n >= 0</c> (otherwise
        /// <c>IllegalArgumentException</c> is thrown)</item>
        /// <item> The result is small enough to fit into a <c>long</c>. The
        /// largest value of <c>n</c> for which <c>n! < Int64.MaxValue</c>
        /// is 20. If the computed value exceeds <c>Int64.MaxValue</c>
        /// an <c>ArithMeticException</c> is thrown.</item>
        /// </list>
        /// </para>
        /// </summary>
        /// <param name="n">argument</param>
        /// <returns><c>n!</c></returns>
        /// <exception cref="MathArithmeticException"> if the result is too large to be represented
        /// by a <c>long</c>.</exception>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="MathArithmeticException"> if <c>n > 20</c>: The factorial value is too
        /// large to fit in a <c>long</c>.</exception>
        public static long factorial(int n)
        {
            if (n < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("FACTORIAL_NEGATIVE_PARAMETER"), n);
            }
            if (n > 20)
            {
                throw new MathArithmeticException();
            }
            return FACTORIALS[n];
        }

        /// <summary>
        /// Compute n!, the <a href="http://mathworld.wolfram.com/Factorial.html">
        /// factorial</a> of <c>n</c> (the product of the numbers 1 to n), as a
        /// <c>double</c>.
        /// The result should be small enough to fit into a <c>double</c>: The
        /// largest <c>n</c> for which <c>n! < Double.MaxValue</c> is 170.
        /// If the computed value exceeds <c>Double.MaxValue</c>,
        /// <c>Double.PositiveInfinity</c> is returned.
        /// </summary>
        /// <param name="n">Argument.</param>
        /// <returns><c>n!</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        public static double factorialDouble(int n)
        {
            if (n < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("FACTORIAL_NEGATIVE_PARAMETER"), n);
            }
            if (n < 21)
            {
                return FACTORIALS[n];
            }
            return FastMath.floor(FastMath.exp(CombinatoricsUtils.factorialLog(n)) + 0.5);
        }

        /// <summary>
        /// Compute the natural logarithm of the factorial of <c>n</c>.
        /// </summary>
        /// <param name="n">Argument.</param>
        /// <returns><c>n!</c></returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        public static double factorialLog(int n)
        {
            if (n < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("FACTORIAL_NEGATIVE_PARAMETER"), n);
            }
            if (n < 21)
            {
                return FastMath.log(FACTORIALS[n]);
            }
            double logSum = 0;
            for (int i = 2; i <= n; i++)
            {
                logSum += FastMath.log(i);
            }
            return logSum;
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
        /// <exception cref="MathArithmeticException"> if some overflow happens, typically
        /// for n exceeding 25 and k between 20 and n-2 (S(n,n-1) is handled specifically
        /// and does not overflow)</exception>
        public static long stirlingS2(int n, int k)
        {
            if (k < 0)
            {
                throw new NotPositiveException<Int32>(k);
            }
            if (k > n)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(k, n, true);
            }

            long[][] stirlingS2;
            lock (STIRLING_S2)
            {
                stirlingS2 = STIRLING_S2;
                if (stirlingS2 == null)
                {
                    // the cache has never been initialized, compute the first numbers
                    // by direct recurrence relation

                    // as S(26,9) = 11201516780955125625 is larger than Long.MAX_VALUE
                    // we must stop computation at row 26
                    int maxIndex = 26;
                    stirlingS2 = new long[maxIndex][];
                    stirlingS2[0] = new long[] { 1L };
                    for (int i = 1; i < stirlingS2.Length; ++i)
                    {
                        stirlingS2[i] = new long[i + 1];
                        stirlingS2[i][0] = 0;
                        stirlingS2[i][1] = 1;
                        stirlingS2[i][i] = 1;
                        for (int j = 2; j < i; ++j)
                        {
                            STIRLING_S2[i][j] = j * stirlingS2[i - 1][j] + stirlingS2[i - 1][j - 1];
                        }
                    }

                    // atomically save the cache, thread-safe
                    STIRLING_S2 = (long[][])stirlingS2.Clone();
                }
            }

            if (n < stirlingS2.Length)
            {
                // the number is in the small cache
                return stirlingS2[n][k];
            }
            else
            {
                // use explicit formula to compute the number without caching it
                if (k == 0)
                {
                    return 0;
                }
                else if (k == 1 || k == n)
                {
                    return 1;
                }
                else if (k == 2)
                {
                    return (1L << (n - 1)) - 1L;
                }
                else if (k == n - 1)
                {
                    return binomialCoefficient(n, 2);
                }
                else
                {
                    // definition formula: note that this may trigger some overflow
                    long sum = 0;
                    long sign = ((k & 0x1) == 0) ? 1 : -1;
                    for (int j = 1; j <= k; ++j)
                    {
                        sign = -sign;
                        sum += sign * binomialCoefficient(k, j) * ArithmeticUtils.pow(j, n);
                        if (sum < 0)
                        {
                            // there was an overflow somewhere
                            throw new MathArithmeticException(new LocalizedFormats("ARGUMENT_OUTSIDE_DOMAIN"), n, 0, stirlingS2.Length - 1);
                        }
                    }
                    return sum / factorial(k);
                }
            }

        }

        /// <summary>
        /// Returns an iterator whose range is the k-element subsets of {0, ..., n - 1}
        /// represented as <c>int[]</c> arrays.
        /// <para>
        /// The arrays returned by the iterator are sorted in descending order and
        /// they are visited in lexicographic order with significance from right to
        /// left. For example, combinationsIterator(4, 2) returns an Iterator that
        /// will generate the following sequence of arrays on successive calls to
        /// <c>next()</c>:<para/>
        /// <c>[0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]</c>
        /// </para>
        /// If <c>k == 0</c> an Iterator containing an empty array is returned and
        /// if <c>k == n</c> an Iterator containing [0, ..., n -1] is returned.
        /// </summary>
        /// <param name="n">Size of the set from which subsets are selected.</param>
        /// <param name="k">Size of the subsets to be enumerated.</param>
        /// <returns>an <c>Iterator iterator</c> over the k-sets in n.</returns>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        public static IEnumerator<int[]> combinationsIterator(int n, int k)
        {
            return new Combinations(n, k).iterator();
        }

        /// <summary>
        /// Check binomial preconditions.
        /// </summary>
        /// <param name="n">Size of the set.</param>
        /// <param name="k">Size of the subsets to be counted.</param>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        public static void checkBinomial(int n, int k)
        {
            if (n < k)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("BINOMIAL_INVALID_PARAMETERS_ORDER"),
                                                    k, n, true);
            }
            if (n < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("BINOMIAL_NEGATIVE_PARAMETER"), n);
            }
        }
    }
}
