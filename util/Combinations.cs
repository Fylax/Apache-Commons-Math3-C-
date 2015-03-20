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
using System;
using System.Collections;
using System.Collections.Generic;

namespace Math3.util
{
    /// <summary>
    /// Utility to create <a href="http://en.wikipedia.org/wiki/Combination">
    /// combinations</a> <c>(n, k)</c> of <c>k</c> elements in a set of
    /// <c>n</c> elements.
    /// </summary>
    public class Combinations : IEnumerable<int[]>
    {
        /// <summary>
        /// Size of the set from which combinations are drawn.
        /// </summary>
        private readonly int n;

        /// <summary>
        /// Number of elements in each combination.
        /// </summary>
        private readonly int k;

        /// <summary>
        /// Iteration order.
        /// </summary>
        private readonly IterationOrder iterationOrder;

        /// <summary>
        /// Describes the type of iteration performed by the
        /// <see cref="iterator()"/>.
        /// </summary>
        private enum IterationOrder
        {
            /// <summary>
            /// Lexicographic order.
            /// </summary>
            LEXICOGRAPHIC
        }

        /// <summary>
        /// Creates an instance whose range is the k-element subsets of
        /// {0, ..., n - 1} represented as <c>int[]</c> arrays.
        /// <para>
        /// The iteration order is lexicographic: the arrays returned by the
        /// <c>iterator()</c> are sorted in descending order and
        /// they are visited in lexicographic order with significance from
        /// right to left.
        /// For example, <c>new Combinations(4, 2).iterator()</c> returns
        /// an iterator that will generate the following sequence of arrays
        /// on successive calls to
        /// <c>next()</c>:<para/>
        /// <c>[0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]</c>
        /// </para>
        /// If <c>k == 0</c> an iterator containing an empty array is returned;
        /// if <c>k == n</c> an iterator containing [0, ..., n - 1] is returned.
        /// 
        /// </summary>
        /// <param name="n">Size of the set from which subsets are selected.</param>
        /// <param name="k"Size of the subsets to be enumerated.></param>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        public Combinations(int n, int k) : this(n, k, IterationOrder.LEXICOGRAPHIC) { }

        /// <summary>
        /// Creates an instance whose range is the k-element subsets of
        /// {0, ..., n - 1} represented as <c>int[]</c> arrays.
        /// <para>
        /// If the <c>iterationOrder</c> argument is set to
        /// <see cref="IterationOrder.LEXICOGRAPHIC"/>, the arrays returned by the
        /// <see cref="iterator()"/> are sorted in descending order and
        /// they are visited in lexicographic order with significance from
        /// right to left.
        /// For example, <c>new Combinations(4, 2).iterator()</c> returns
        /// an iterator that will generate the following sequence of arrays
        /// on successive calls to
        /// <c>next()</c>:<para/>
        /// <c>[0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]</c>
        /// </para>
        /// If <c>k == 0</c> an iterator containing an empty array is returned;
        /// if <c>k == n</c> an iterator containing [0, ..., n - 1] is returned.
        /// 
        /// </summary>
        /// <param name="n">Size of the set from which subsets are selected.</param>
        /// <param name="k">Size of the subsets to be enumerated.></param>
        /// <param name="iterationOrder">Specifies the <see cref="iterator()"/>.</param>
        /// <exception cref="NotPositiveException"> if <c>n < 0</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if <c>k > n</c>.</exception>
        private Combinations(int n, int k, IterationOrder iterationOrder)
        {
            CombinatoricsUtils.checkBinomial(n, k);
            this.n = n;
            this.k = k;
            this.iterationOrder = iterationOrder;
        }

        /// <summary>
        /// Gets the size of the set from which combinations are drawn.
        /// </summary>
        /// <returns>the size of the universe.</returns>
        public int getN()
        {
            return n;
        }

        /// <summary>
        /// Gets the number of elements in each combination. 
        /// </summary>
        /// <returns>the size of the subsets to be enumerated.</returns>
        public int getK()
        {
            return k;
        }

        /// <inheritdoc cref="GetEnumerator"/>
        public IEnumerator<int[]> iterator()
        {
            if (k == 0 ||
                k == n)
            {
                return new SingletonIterator(MathArrays.natural(k));
            }

            switch (iterationOrder)
            {
                case IterationOrder.LEXICOGRAPHIC:
                    return new LexicographicIterator(n, k);
                default:
                    throw new MathInternalError(); // Should never happen.
            }
        }

        /// <inheritdoc/>
        IEnumerator IEnumerable.GetEnumerator()
        {
            return (IEnumerator)iterator();
        }

        /// <inheritdoc cref="GetEnumerator"/>
        public IEnumerator<int[]> GetEnumerator()
        {
            return iterator();
        }

        /// <summary>
        /// Defines a lexicographic ordering of combinations.
        /// The returned comparator allows to compare any two combinations
        /// that can be produced by this instance's <c>iterator()</c>.
        /// Its <c>ompare(int[],int[])</c> method will throw exceptions if
        /// passed combinations that are inconsistent with this instance:
        /// <list type="bullet">
        /// <item><c>DimensionMismatchException</c> if the array lengths are not
        /// equal to <c>k</c></item>
        /// <item><c>OutOfRangeException</c> if an element of the array is not
        /// within the interval [0, <c>n</c>).</item>
        /// </list>
        /// </summary>
        /// <returns>a lexicographic comparator.</returns>
        public IComparer<int[]> comparator()
        {
            return new LexicographicComparator(n, k);
        }

        /// <summary>
        /// Lexicographic combinations iterator.
        /// <para>
        /// Implementation follows Algorithm T in The Art of Computer Programming
        /// Internet Draft (PRE-FASCICLE 3A), "A Draft of Section 7.2.1.3 Generating All
        /// Combinations, D. Knuth, 2004.</para>
        /// <para>
        /// The degenerate cases <c>k == 0</c> and <c>k == n</c> are NOT handled by this
        /// implementation.  If constructor arguments satisfy <c>k == 0</c>
        /// or <c>k >= n</c>, no exception is generated, but the iterator is empty.
        /// </para>
        /// </summary>
        private class LexicographicIterator : IEnumerator<int[]>
        {
            /// <summary>
            /// Size of subsets returned by the iterator
            /// </summary>
            private readonly int k;

            /// <summary>
            /// c[1], ..., c[k] stores the next combination; c[k + 1], c[k + 2] are
            /// sentinels.
            /// <para>
            /// Note that c[0] is "wasted" but this makes it a little easier to
            /// follow the code.
            /// </para>
            /// </summary>
            private readonly int[] c;

            /// <summary>
            /// Return value for <see cref="hasNext()"/>
            /// </summary>
            private Boolean more = true;

            /// <summary>
            /// Marker: smallest index such that c[j + 1] > j
            /// </summary>
            private int j;

            /// <summary>
            /// Construct a CombinationIterator to enumerate k-sets from n.
            /// <para>
            /// NOTE: If <c>k === 0</c> or <c>k >= n</c>, the Iterator will be empty
            /// (that is, <see cref="hasNext()"/> will return <c>false</c> immediately.
            /// </para>
            /// </summary>
            /// <param name="n">size of the set from which subsets are enumerated</param>
            /// <param name="k">size of the subsets to enumerate</param>
            public LexicographicIterator(int n, int k)
            {
                this.k = k;
                c = new int[k + 3];
                if (k == 0 || k >= n)
                {
                    more = false;
                    return;
                }
                // Initialize c to start with lexicographically first k-set
                for (int i = 1; i <= k; i++)
                {
                    c[i] = i - 1;
                }
                // Initialize sentinels
                c[k + 1] = n;
                c[k + 2] = 0;
                j = k; // Set up invariant: j is smallest index such that c[j + 1] > j
            }

            /// <inheritdoc cref="SingletonIterator.hasNext()"/>
            public Boolean hasNext()
            {
                return more;
            }

            /// <inheritdoc cref="MoveNext()"/>
            public int[] next()
            {
                // Copy return value (prepared by last activation)
                int[] ret = new int[k];
                Array.Copy(c, 1, ret, 0, k);

                // Prepare next iteration
                // T2 and T6 loop
                int x = 0;
                if (j > 0)
                {
                    x = j;
                    c[j] = x;
                    j--;
                    return ret;
                }
                // T3
                if (c[1] + 1 < c[2])
                {
                    c[1]++;
                    return ret;
                }
                else
                {
                    j = 2;
                }
                // T4
                Boolean stepDone = false;
                while (!stepDone)
                {
                    c[j - 1] = j - 2;
                    x = c[j] + 1;
                    if (x == c[j + 1])
                    {
                        j++;
                    }
                    else
                    {
                        stepDone = true;
                    }
                }
                // T5
                if (j > k)
                {
                    more = false;
                    return ret;
                }
                // T6
                c[j] = x;
                j--;
                return ret;
            }

            /// <inheritdoc/>
            public bool MoveNext()
            {
                this.next();
                return this.more;
            }

            /// <summary>
            /// Not supported.
            /// </summary>
            public void remove()
            {
                throw new NotSupportedException();
            }

            /// <summary>
            /// Not supported.
            /// </summary>
            public void Reset()
            {
                this.remove();
            }

            public void Dispose()
            {
                return;
            }

            object IEnumerator.Current
            {
                get
                {
                    return Current;
                }
            }

            public int[] Current
            {
                get
                {
                    try
                    {
                        return c;
                    }
                    catch (IndexOutOfRangeException)
                    {
                        throw new InvalidOperationException();
                    }
                }
            }
        }

        /// <summary>
        /// Iterator with just one element to handle degenerate cases (full array,
        /// empty array) for combination iterator.
        /// </summary>
        private class SingletonIterator : IEnumerator<int[]>
        {
            /// <summary>
            /// Singleton array
            /// </summary>
            private readonly int[] singleton;
            
            /// <summary>
            /// True on initialization, false after first call to next
            /// </summary>
            private Boolean more = true;
            
            /// <summary>
            /// Create a singleton iterator providing the given array.
            /// </summary>
            /// <param name="singleton">array returned by the iterator</param>
            public SingletonIterator(int[] singleton)
            {
                this.singleton = singleton;
            }
            
            /// <returns>True until next is called the first time, then false</returns>
            public Boolean hasNext()
            {
                return more;
            }
            
            /// <summary></summary>
            /// <returnsthe singleton in first activation; throws NSEE thereafterreturns>
            public int[] next()
            {
                more = false;
                return singleton;
            }

            /// <inheritdoc cref="next()"/>
            public bool MoveNext()
            {
                this.next();
                return this.more;
            }
            
            /// <summary>
            /// Not supported
            /// </summary>
            public void remove()
            {
                throw new NotSupportedException();
            }

            /// <summary>
            /// Not supported
            /// </summary>
            public void Reset()
            {
                this.remove();
            }

            public void Dispose()
            {
                return;
            }

            object IEnumerator.Current
            {
                get
                {
                    return Current;
                }
            }

            public int[] Current
            {
                get
                {
                    try
                    {
                        return singleton;
                    }
                    catch (IndexOutOfRangeException)
                    {
                        throw new InvalidOperationException();
                    }
                }
            }
        }

        /// <summary>
        /// Defines the lexicographic ordering of combinations, using
        /// the <see cref="exNorm(int[])"/> method.
        /// </summary>
        private class LexicographicComparator : IComparer<int[]>
        {
            /// <summary>
            /// Size of the set from which combinations are drawn.
            /// </summary>
            private int n;

            /// <summary>
            /// Number of elements in each combination.
            /// </summary>
            private int k;

            /// <summary></summary>
            /// <param name="n">Size of the set from which subsets are selected.</param>
            /// <param name="k">Size of the subsets to be enumerated.</param>
            public LexicographicComparator(int n, int k)
            {
                this.n = n;
                this.k = k;
            }

            /// <inheritdoc cref="Compare"/>
            public int compare(int[] c1, int[] c2)
            {
                if (c1.Length != k)
                {
                    throw new DimensionMismatchException(c1.Length, k);
                }
                if (c2.Length != k)
                {
                    throw new DimensionMismatchException(c2.Length, k);
                }

                // Method "lexNorm" works with ordered arrays.
                int[] c1s = MathArrays.copyOf(c1);
                Array.Sort(c1s);
                int[] c2s = MathArrays.copyOf(c2);
                Array.Sort(c2s);

                long v1 = lexNorm(c1s);
                long v2 = lexNorm(c2s);

                if (v1 < v2)
                {
                    return -1;
                }
                else if (v1 > v2)
                {
                    return 1;
                }
                else
                {
                    return 0;
                }
            }

            /// <inheritdoc/>
            /// <exception cref="DimensionMismatchException"> if the array lengths are not
            /// equal to <c>k</c>.</exception>
            /// <exception cref="OutOfRangeException"> if an element of the array is not
            /// within the interval [0, <c>n</c>).</exception>
            public int Compare(int[] c1, int[] c2)
            {
                return this.compare(c1, c2);
            }

            /// <summary>
            /// Computes the value (in base 10) represented by the digit
            /// (interpreted in base <c>n</c>) in the input array in reverse
            /// order.
            /// For example if <c>c</c> is <c>{3, 2, 1</c>}, and <c>n</c>
            /// is 3, the method will return 18.
            /// </summary>
            /// <param name="c">Input array.</param>
            /// <returns>the lexicographic norm.</returns>
            /// <exception cref="OutOfRangeException"> if an element of the array is not
            /// within the interval [0, <c>n</c>).</exception>
            private long lexNorm(int[] c)
            {
                long ret = 0;
                for (int i = 0; i < c.Length; i++)
                {
                    int digit = c[i];
                    if (digit < 0 ||
                        digit >= n)
                    {
                        throw new OutOfRangeException<Int32>(digit, 0, n - 1);
                    }
                    ret += c[i] * ArithmeticUtils.pow(n, i);
                }
                return ret;
            }
        }
    }
}
