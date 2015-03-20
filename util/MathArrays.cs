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
using Math3.distribution;
using Math3.exception;
using Math3.exception.util;
using Math3.random;
using System;
using System.Collections.Generic;

namespace Math3.util
{
    /// <summary>
    /// Arrays utilities.
    /// </summary>
    public class MathArrays
    {
        /// <summary>
        /// Factor used for splitting double numbers: n = 2^27 + 1 (i.e. <value/>).
        /// </summary>
        private const int SPLIT_FACTOR = 0x8000001;

        /// <summary>
        /// Private constructor.
        /// </summary>
        private MathArrays() { }

        /// <summary>
        /// Real-valued function that operate on an array or a part of it.
        /// </summary>
        public interface Function
        {
            /// <summary>
            /// Operates on an entire array.
            /// </summary>
            /// <param name="array">Array to operate on.</param>
            /// <returns>the result of the operation.</returns>
            double evaluate(double[] array);

            /// <param name="array">Array to operate on.</param>
            /// <param name="startIndex">Index of the first element to take into account.</param>
            /// <param name="numElements">Number of elements to take into account.</param>
            /// <returns>the result of the operation.</returns>
            double evaluate(double[] array, int startIndex, int numElements);
        }

        /// <summary>
        /// Create a copy of an array scaled by a value.
        /// </summary>
        /// <param name="val">Array to scale.</param>
        /// <param name="arr">Scalar.</param>
        /// <returns>scaled copy of array with each entry multiplied by val.</returns>
        public static double[] scale(double val, double[] arr)
        {
            double[] newArr = new double[arr.Length];
            for (int i = 0; i < arr.Length; i++)
            {
                newArr[i] = arr[i] * val;
            }
            return newArr;
        }

        /// <summary>
        /// <para>Multiply each element of an array by a value.</para>
        /// <para>The array is modified in place (no copy is created).</para>
        /// </summary>
        /// <param name="val">Array to scale</param>
        /// <param name="arr">Scalar</param>
        public static void scaleInPlace(double val, double[] arr)
        {
            for (int i = 0; i < arr.Length; i++)
            {
                arr[i] *= val;
            }
        }

        /// <summary>
        /// Creates an array whose contents will be the element-by-element
        /// addition of the arguments.
        /// </summary>
        /// <param name="a">First term of the addition.</param>
        /// <param name="b">Second term of the addition.</param>
        /// <returns>a new array <c>r</c> where <c>r[i] = a[i] + b[i]</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the array lengths differ.</exception>
        public static double[] ebeAdd(double[] a, double[] b)
        {
            if (a.Length != b.Length)
            {
                throw new DimensionMismatchException(a.Length, b.Length);
            }

            double[] result = (Double[])a.Clone();
            for (int i = 0; i < a.Length; i++)
            {
                result[i] += b[i];
            }
            return result;
        }

        /// <summary>
        /// Creates an array whose contents will be the element-by-element
        /// subtraction of the second argument from the first.
        /// </summary>
        /// <param name="a">First term.</param>
        /// <param name="b">Element to be subtracted.</param>
        /// <returns>a new array <c>r</c> where <c>r[i] = a[i] - b[i]</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the array lengths differ.</exception>
        public static double[] ebeSubtract(double[] a, double[] b)
        {
            if (a.Length != b.Length)
            {
                throw new DimensionMismatchException(a.Length, b.Length);
            }

            double[] result = (Double[])a.Clone();
            for (int i = 0; i < a.Length; i++)
            {
                result[i] -= b[i];
            }
            return result;
        }

        /// <summary>
        /// Creates an array whose contents will be the element-by-element
        /// multiplication of the arguments.
        /// </summary>
        /// <param name="a">First factor of the multiplication.</param>
        /// <param name="b">Second factor of the multiplication.</param>
        /// <returns>a new array <c>r</c> where <c>r[i] = a[i] * b[i]</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the array lengths differ.</exception>
        public static double[] ebeMultiply(double[] a, double[] b)
        {
            if (a.Length != b.Length)
            {
                throw new DimensionMismatchException(a.Length, b.Length);
            }
            double[] result = (Double[])a.Clone();
            for (int i = 0; i < a.Length; i++)
            {
                result[i] *= b[i];
            }
            return result;
        }

        /// <summary>
        /// Creates an array whose contents will be the element-by-element
        /// division of the first argument by the second.
        /// </summary>
        /// <param name="a">Numerator of the division.</param>
        /// <param name="b">Denominator of the division.</param>
        /// <returns>a new array <c>r</c> where <c>r[i] = a[i] / b[i]</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if the array lengths differ.</exception>
        public static double[] ebeDivide(double[] a, double[] b)
        {
            if (a.Length != b.Length)
            {
                throw new DimensionMismatchException(a.Length, b.Length);
            }
            double[] result = (Double[])a.Clone();
            for (int i = 0; i < a.Length; i++)
            {
                result[i] /= b[i];
            }
            return result;
        }

        /// <summary>
        /// Calculates the L_1 (sum of abs) distance between two points.
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns>the L_1 distance between the two points</returns>
        public static double distance1(double[] p1, double[] p2)
        {
            double sum = 0;
            for (int i = 0; i < p1.Length; i++)
            {
                sum += FastMath.abs(p1[i] - p2[i]);
            }
            return sum;
        }

        /// <summary>
        /// Calculates the L_1 (sum of abs) distance between two points.
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns>the L_1 distance between the two points</returns>
        public static int distance1(int[] p1, int[] p2)
        {
            int sum = 0;
            for (int i = 0; i < p1.Length; i++)
            {
                sum += FastMath.abs(p1[i] - p2[i]);
            }
            return sum;
        }

        /// <summary>
        /// Calculates the L_2 (Euclidean) distance between two points.
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns>the L_2 distance between the two points</returns>
        public static double distance(double[] p1, double[] p2)
        {
            double sum = 0;
            for (int i = 0; i < p1.Length; i++)
            {
                double dp = p1[i] - p2[i];
                sum += dp * dp;
            }
            return FastMath.sqrt(sum);
        }

        /// <summary>
        /// Calculates the L_2 (Euclidean) distance between two points.
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns>the L_2 distance between the two points</returns>
        public static double distance(int[] p1, int[] p2)
        {
            double sum = 0;
            for (int i = 0; i < p1.Length; i++)
            {
                double dp = p1[i] - p2[i];
                sum += dp * dp;
            }
            return FastMath.sqrt(sum);
        }

        /// <summary>
        /// Calculates the L_&infin; (max of abs) distance between two points.
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns>the L_&infin; distance between the two points</returns>
        public static double distanceInf(double[] p1, double[] p2)
        {
            double max = 0;
            for (int i = 0; i < p1.Length; i++)
            {
                max = FastMath.max(max, FastMath.abs(p1[i] - p2[i]));
            }
            return max;
        }

        /// <summary>
        /// Calculates the L_&infin; (max of abs) distance between two points.
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns>the L_&infin; distance between the two points</returns>
        public static int distanceInf(int[] p1, int[] p2)
        {
            int max = 0;
            for (int i = 0; i < p1.Length; i++)
            {
                max = FastMath.max(max, FastMath.abs(p1[i] - p2[i]));
            }
            return max;
        }

        /// <summary>
        /// Specification of ordering direction.
        /// </summary>
        public enum OrderDirection
        {
            /// <summary>
            /// Constant for increasing direction.
            /// </summary>
            INCREASING,

            /// <summary>
            /// Constant for decreasing direction.
            /// </summary>
            DECREASING
        }

        /// <summary>
        /// Check that an array is monotonically increasing or decreasing.
        /// </summary>
        /// <typeparam name="T">the type of the elements in the specified array</typeparam>
        /// <param name="val">Values.</param>
        /// <param name="dir">Ordering direction.</param>
        /// <param name="strict">Whether the order should be strict.</param>
        /// <returns><c>true</c> if sorted, <c>false</c> otherwise.</returns>
        public static Boolean isMonotonic<T>(T[] val, OrderDirection dir, Boolean strict) where T : IComparable<T>
        {
            T previous = val[0];
            int max = val.Length;
            for (int i = 1; i < max; i++)
            {
                int comp;
                switch (dir)
                {
                    case OrderDirection.INCREASING:
                        comp = previous.CompareTo(val[i]);
                        if (strict)
                        {
                            if (comp >= 0)
                            {
                                return false;
                            }
                        }
                        else
                        {
                            if (comp > 0)
                            {
                                return false;
                            }
                        }
                        break;
                    case OrderDirection.DECREASING:
                        comp = val[i].CompareTo(previous);
                        if (strict)
                        {
                            if (comp >= 0)
                            {
                                return false;
                            }
                        }
                        else
                        {
                            if (comp > 0)
                            {
                                return false;
                            }
                        }
                        break;
                    default:
                        // Should never happen.
                        throw new MathInternalError();
                }
                previous = val[i];
            }
            return true;
        }

        /// <summary>
        /// Check that an array is monotonically increasing or decreasing.
        /// </summary>
        /// <param name="val">Values.</param>
        /// <param name="dir">Ordering direction.</param>
        /// <param name="strict">Whether the order should be strict.</param>
        /// <returns><c>true</c> if sorted, <c>false</c> otherwise.</returns>
        public static Boolean isMonotonic(double[] val, OrderDirection dir, Boolean strict)
        {
            return checkOrder(val, dir, strict, false);
        }

        /// <summary>
        /// Check that the given array is sorted.
        /// </summary>
        /// <param name="val">Values.</param>
        /// <param name="dir">Ordering direction.</param>
        /// <param name="strict">Whether the order should be strict.</param>
        /// <param name="abort">Whether to throw an exception if the check fails.</param>
        /// <returns><c>true</c> if the array is sorted.</returns>
        /// and <c>abort</c> is <c>true</c>.</exception>
        public static Boolean checkOrder(double[] val, OrderDirection dir, Boolean strict, Boolean abort)
        {
            double previous = val[0];
            int max = val.Length;
            int index;
            for (index = 1; index < max; ++index)
            {
                switch (dir)
                {
                    case OrderDirection.INCREASING:
                        if (strict)
                        {
                            if (val[index] <= previous)
                            {
                                goto ITEM;
                            }
                        }
                        else
                        {
                            if (val[index] < previous)
                            {
                                goto ITEM;
                            }
                        }
                        break;
                    case OrderDirection.DECREASING:
                        if (strict)
                        {
                            if (val[index] >= previous)
                            {
                                goto ITEM;
                            }
                        }
                        else
                        {
                            if (val[index] > previous)
                            {
                                goto ITEM;
                            }
                        }
                        break;
                    default:
                        // Should never happen.
                        throw new MathInternalError();
                }

                previous = val[index];
            }
        ITEM:
            if (index == max)
            {
                // Loop completed.
                return true;
            }

            // Loop early exit means wrong ordering.
            if (abort)
            {
                throw new NonMonotonicSequenceException<Double, Double>(val[index], previous, index, dir, strict);
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// Check that the given array is sorted.
        /// </summary>
        /// <param name="val">Values.</param>
        /// <param name="dir">Ordering direction.</param>
        /// <param name="strict">Whether the order should be strict.</param>
        /// <exception cref="NonMonotonicSequenceException"> if the array is not 
        /// sorted.</exception>
        public static void checkOrder(double[] val, OrderDirection dir, Boolean strict)
        {
            checkOrder(val, dir, strict, true);
        }

        /// <summary>
        /// Check that the given array is sorted in strictly increasing order.
        /// </summary>
        /// <param name="val">Values.</param>
        /// <exception cref="NonMonotonicSequenceException"> if the array is not
        /// sorted.</exception>
        public static void checkOrder(double[] val)
        {
            checkOrder(val, OrderDirection.INCREASING, true);
        }

        /// <summary>
        /// Throws DimensionMismatchException if the input array is not rectangular.
        /// </summary>
        /// <param name="inp">array to be tested</param>
        /// <exception cref="NullArgumentException"> if input array is null</exception>
        /// <exception cref="DimensionMismatchException"> if input array is not
        /// rectangular</exception>
        public static void checkRectangular(long[][] inp)
        {
            MathUtils.checkNotNull(inp);
            for (int i = 1; i < inp.Length; i++)
            {
                if (inp[i].Length != inp[0].Length)
                {
                    throw new DimensionMismatchException(new LocalizedFormats("DIFFERENT_ROWS_LENGTHS"), inp[i].Length, inp[0].Length);
                }
            }
        }

        /// <summary>
        /// Check that all entries of the input array are strictly positive.
        /// </summary>
        /// <param name="inp">Array to be tested</param>
        /// <exception cref="NotStrictlyPositiveException"> if any entries of the array
        /// are not strictly positive.</exception>
        public static void checkPositive(double[] inp)
        {
            for (int i = 0; i < inp.Length; i++)
            {
                if (inp[i] <= 0)
                {
                    throw new NotStrictlyPositiveException<Double>(inp[i]);
                }
            }
        }

        /// <summary>
        /// Check that no entry of the input array is <c>NaN</c>.
        /// </summary>
        /// <param name="inp">Array to be tested.</param>
        /// <exception cref="NotANumberException"> if an entry is <c>NaN</c>.</exception>
        public static void checkNotNaN(double[] inp)
        {
            for (int i = 0; i < inp.Length; i++)
            {
                if (Double.IsNaN(inp[i]))
                {
                    throw new NotANumberException();
                }
            }
        }

        /// <summary>
        /// Check that all entries of the input array are >= 0.
        /// </summary>
        /// <param name="inp">Array to be tested</param>
        /// <exception cref="NotPositiveException"> if any array entries are
        /// less than 0.</exception>
        public static void checkNonNegative(long[] inp)
        {
            for (int i = 0; i < inp.Length; i++)
            {
                if (inp[i] < 0)
                {
                    throw new NotPositiveException<Int64>(inp[i]);
                }
            }
        }

        /// <summary>
        /// Check that all entries of the input array are >= 0.
        /// </summary>
        /// <param name="inp">Array to be tested</param>
        /// <exception cref="NotPositiveException"> if any array entries are
        /// less than 0.</exception>
        public static void checkNonNegative(long[][] inp)
        {
            for (int i = 0; i < inp.Length; i++)
            {
                for (int j = 0; j < inp[i].Length; j++)
                {
                    if (inp[i][j] < 0)
                    {
                        throw new NotPositiveException<Int64>(inp[i][j]);
                    }
                }
            }
        }

        /// <summary>
        /// Returns the Cartesian norm (2-norm), handling both overflow and underflow.
        /// Translation of the minpack enorm subroutine.
        /// The redistribution policy for MINPACK is available
        /// <a href="http://www.netlib.org/minpack/disclaimer">here</a>, for
        /// convenience, it is reproduced below.<para/>
        /// <list type="table">
        /// <item>
        /// Minpack Copyright Notice (1999) University of Chicago.
        /// All rights reserved
        /// </item>
        /// <item>
        /// Redistribution and use in source and binary forms, with or without
        /// modification, are permitted provided that the following conditions
        /// are met:
        /// <list type="number">
        /// <item>Redistributions of source code must retain the above copyright
        /// notice, this list of conditions and the following disclaimer.</item>
        /// <item>Redistributions in binary form must reproduce the above
        /// copyright notice, this list of conditions and the following
        /// disclaimer in the documentation and/or other materials provided
        /// with the distribution.</item>
        /// <item>The end-user documentation included with the redistribution, if any,
        /// must include the following acknowledgment:
        /// <c>This product includes software developed by the University of
        /// Chicago, as Operator of Argonne National Laboratory.</c>
        /// Alternately, this acknowledgment may appear in the software itself,
        /// if and wherever such third-party acknowledgments normally appear.</item>
        /// <item>WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
        /// WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
        /// UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
        /// THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
        /// IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
        /// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
        /// OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
        /// OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
        /// USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
        /// THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
        /// DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
        /// UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
        /// BE CORRECTED.</item>
        /// <item>LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
        /// HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
        /// ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
        /// INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
        /// ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
        /// PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
        /// SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
        /// (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
        /// EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
        /// POSSIBILITY OF SUCH LOSS OR DAMAGES.</item>
        /// </list></item>
        /// </list>
        /// </summary>
        /// <param name="v">Vector of doubles.</param>
        /// <returns>the 2-norm of the vector.</returns>
        public static double safeNorm(double[] v)
        {
            double rdwarf = 3.834e-20;
            double rgiant = 1.304e+19;
            double s1 = 0;
            double s2 = 0;
            double s3 = 0;
            double x1max = 0;
            double x3max = 0;
            double floatn = v.Length;
            double agiant = rgiant / floatn;
            for (int i = 0; i < v.Length; i++)
            {
                double xabs = FastMath.abs(v[i]);
                if (xabs < rdwarf || xabs > agiant)
                {
                    if (xabs > rdwarf)
                    {
                        if (xabs > x1max)
                        {
                            double r = x1max / xabs;
                            s1 = 1 + s1 * r * r;
                            x1max = xabs;
                        }
                        else
                        {
                            double r = xabs / x1max;
                            s1 += r * r;
                        }
                    }
                    else
                    {
                        if (xabs > x3max)
                        {
                            double r = x3max / xabs;
                            s3 = 1 + s3 * r * r;
                            x3max = xabs;
                        }
                        else
                        {
                            if (xabs != 0)
                            {
                                double r = xabs / x3max;
                                s3 += r * r;
                            }
                        }
                    }
                }
                else
                {
                    s2 += xabs * xabs;
                }
            }
            double norm;
            if (s1 != 0)
            {
                norm = x1max * Math.Sqrt(s1 + (s2 / x1max) / x1max);
            }
            else
            {
                if (s2 == 0)
                {
                    norm = x3max * Math.Sqrt(s3);
                }
                else
                {
                    if (s2 >= x3max)
                    {
                        norm = Math.Sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
                    }
                    else
                    {
                        norm = Math.Sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
                    }
                }
            }
            return norm;
        }

        /// <summary>
        /// Sort an array in ascending order in place and perform the same reordering
        /// of entries on other arrays. For example, if
        /// <c>x = [3, 1, 2], y = [1, 2, 3]</c> and <c>z = [0, 5, 7]</c>, then
        /// <c>sortInPlace(x, y, z)</c> will update <c>x</c> to <c>[1, 2, 3]</c>,
        /// <c>y</c> to <c>[2, 3, 1]</c> and <c>z</c> to <c>[5, 7, 0]</c>.
        /// </summary>
        /// <param name="x">Array to be sorted and used as a pattern for permutation
        /// of the other arrays.</param>
        /// <param name="yList">Set of arrays whose permutations of entries will follow
        /// those performed on <c>x</c>.</param>
        /// <exception cref="">DimensionMismatchException if any <c>y</c> is not the same
        /// size as <c>x</c>.</exception>
        /// <exception cref="NullArgumentException"> if <c>x</c> or any <c>y</c> is null.
        /// </exception>
        public static void sortInPlace(double[] x, params double[][] yList)
        {
            sortInPlace(x, OrderDirection.INCREASING, yList);
        }

        /// <summary>
        /// Sort an array in place and perform the same reordering of entries on
        /// other arrays.  This method works the same as the other
        /// <see cref="sortInPlace(double[], double[][])">sortInPlace</see> method, but
        /// allows the order of the sort to be provided in the <c>dir</c>
        /// parameter.
        /// </summary>
        /// <param name="x">Array to be sorted and used as a pattern for permutation
        /// of the other arrays.</param>
        /// <param name="dir">Order direction.</param>
        /// <param name="yList">Set of arrays whose permutations of entries will follow
        /// those performed on <c>x</c>.</param>
        /// <exception cref="DimensionMismatchException"> if any <c>y</c> is not the same
        /// size as <c>x</c>.</exception>
        /// <exception cref="NullArgumentException"> if <c>x</c> or any <c>y</c> is null
        /// </exception>
        public static void sortInPlace(double[] x, OrderDirection dir, params double[][] yList)
        {
            // Consistency checks.
            if (x == null)
            {
                throw new NullArgumentException();
            }

            int yListLen = yList.Length;
            int len = x.Length;

            for (int j = 0; j < yListLen; j++)
            {
                double[] y = yList[j];
                if (y == null)
                {
                    throw new NullArgumentException();
                }
                if (y.Length != len)
                {
                    throw new DimensionMismatchException(y.Length, len);
                }
            }

            // Associate each abscissa "x[i]" with its index "i".
            List<KeyValuePair<Double, Int32>> list = new List<KeyValuePair<Double, Int32>>(len);
            for (int i = 0; i < len; i++)
            {
                list.Add(new KeyValuePair<Double, Int32>(x[i], i));
            }

            // Sort.
            if (dir == MathArrays.OrderDirection.INCREASING)
            {
                list.Sort((o1, o2) => o1.Key.CompareTo(o2.Key));
            }
            else
            {
                list.Sort((o1, o2) => o2.Key.CompareTo(o1.Key));
            }

            // Modify the original array so that its elements are in
            // the prescribed order.
            // Retrieve indices of original locations.
            int[] indices = new int[len];
            for (int i = 0; i < len; i++)
            {
                KeyValuePair<Double, Int32> e = list[i];
                x[i] = e.Key;
                indices[i] = e.Value;
            }

            // In each of the associated arrays, move the
            // elements to their new location.
            for (int j = 0; j < yListLen; j++)
            {
                // Input array will be modified in place.
                double[] yInPlace = yList[j];
                double[] yOrig = (Double[])yInPlace.Clone();

                for (int i = 0; i < len; i++)
                {
                    yInPlace[i] = yOrig[indices[i]];
                }
            }
        }

        /// <summary>
        /// Creates a copy of the <c>source</c> array.
        /// </summary>
        /// <param name="source">Array to be copied.</param>
        /// <returns>the copied array.</returns>
        public static int[] copyOf(int[] source)
        {
            return copyOf(source, source.Length);
        }

        /// <summary>
        /// Creates a copy of the <c>source<c> array.
        /// </summary>
        /// <param name="source">Array to be copied.</param>
        /// <returns>the copied array.</returns>
        public static double[] copyOf(double[] source)
        {
            return copyOf(source, source.Length);
        }

        /// <summary>
        /// Creates a copy of the <c>source</c> array.
        /// </summary>
        /// <param name="source">Array to be copied.</param>
        /// <param name="len">Number of entries to copy. If smaller then the source
        /// length, the copy will be truncated, if larger it will padded with
        /// zeroes.</param>
        /// <returns>the copied array.</returns>
        public static int[] copyOf(int[] source, int len)
        {
            int[] output = new int[len];
            Array.Copy(source, 0, output, 0, FastMath.min(len, source.Length));
            return output;
        }

        /// <summary>
        /// Creates a copy of the <c>source</c> array.
        /// </summary>
        /// <param name="source">Array to be copied.</param>
        /// <param name="len">Number of entries to copy. If smaller then the source
        /// length, the copy will be truncated, if larger it will padded with
        /// zeroes.</param>
        /// <returns>the copied array.</returns>
        public static double[] copyOf(double[] source, int len)
        {
            double[] output = new double[len];
            Array.Copy(source, 0, output, 0, FastMath.min(len, source.Length));
            return output;
        }

        /// <summary>
        /// Creates a copy of the <c>source</c> array.
        /// </summary>
        /// <param name="source">Array to be copied.</param>
        /// <param name="from">Initial index of the range to be copied, inclusive.</param>
        /// <param name="to">Final index of the range to be copied, exclusive. 
        /// (This index may lie outside the array.)</param>
        /// <returns>the copied array.</returns>
        public static double[] copyOfRange(double[] source, int from, int to)
        {
            int len = to - from;
            double[] output = new double[len];
            Array.Copy(source, from, output, 0, FastMath.min(len, source.Length - from));
            return output;
        }

        /// <summary>
        /// Compute a linear combination accurately.
        /// This method computes the sum of the products
        /// <c>a_i b_i</c> to high accuracy.
        /// It does so by using specific multiplication and addition algorithms to
        /// preserve accuracy and reduce cancellation effects.
        /// <para/>
        /// It is based on the 2005 paper
        /// <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
        /// Accurate Sum and Dot Product</a> by Takeshi Ogita, Siegfried M. Rump,
        /// and Shin'ichi Oishi published in SIAM J. Sci. Comput.
        /// </summary>
        /// <param name="a">Factors.</param>
        /// <param name="b">Factors.</param>
        /// <returns><c>&Sigma;_i a_i b_i</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if arrays dimensions don't 
        /// match</exception>
        public static double linearCombination(double[] a, double[] b)
        {
            int len = a.Length;
            if (len != b.Length)
            {
                throw new DimensionMismatchException(len, b.Length);
            }

            if (len == 1)
            {
                // Revert to scalar multiplication.
                return a[0] * b[0];
            }

            double[] prodHigh = new double[len];
            double prodLowSum = 0;

            for (int i = 0; i < len; i++)
            {
                double ai = a[i];
                double ca = SPLIT_FACTOR * ai;
                double aHigh = ca - (ca - ai);
                double aLow = ai - aHigh;

                double bi = b[i];
                double cb = SPLIT_FACTOR * bi;
                double bHigh = cb - (cb - bi);
                double bLow = bi - bHigh;
                prodHigh[i] = ai * bi;
                double prodLow = aLow * bLow - (((prodHigh[i] -
                                                        aHigh * bHigh) -
                                                       aLow * bHigh) -
                                                      aHigh * bLow);
                prodLowSum += prodLow;
            }


            double prodHighCur = prodHigh[0];
            double prodHighNext = prodHigh[1];
            double sHighPrev = prodHighCur + prodHighNext;
            double sPrime = sHighPrev - prodHighNext;
            double sLowSum = (prodHighNext - (sHighPrev - sPrime)) + (prodHighCur - sPrime);

            int lenMinusOne = len - 1;
            for (int i = 1; i < lenMinusOne; i++)
            {
                prodHighNext = prodHigh[i + 1];
                double sHighCur = sHighPrev + prodHighNext;
                sPrime = sHighCur - prodHighNext;
                sLowSum += (prodHighNext - (sHighCur - sPrime)) + (sHighPrev - sPrime);
                sHighPrev = sHighCur;
            }

            double result = sHighPrev + (prodLowSum + sLowSum);

            if (Double.IsNaN(result))
            {
                // either we have split infinite numbers or some coefficients were NaNs,
                // just rely on the naive implementation and let IEEE754 handle this
                result = 0;
                for (int i = 0; i < len; ++i)
                {
                    result += a[i] * b[i];
                }
            }

            return result;
        }

        /// <summary>
        /// Compute a linear combination accurately.
        /// <para>
        /// This method computes a_1 &times; b_1 +
        /// a_2 &times; b_2 to high accuracy. It does
        /// so by using specific multiplication and addition algorithms to
        /// preserve accuracy and reduce cancellation effects. It is based
        /// on the 2005 paper <a
        /// href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
        /// Accurate Sum and Dot Product</a> by Takeshi Ogita,
        /// Siegfried M. Rump, and Shin'ichi Oishi published in SIAM J. Sci. Comput.
        /// </para>
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2"second factor of the second term></param>
        /// <returns>a_1 &times; b_1 +
        /// a_2 &times; b_2</returns>
        /// <remarks>
        /// See <see cref="linearCombination(double, double, double, double, double, double)"/>
        /// <para/>
        /// <see cref="linearCombination(double, double, double, double, double, double, double, double)"/>
        /// </remarks>
        public static double linearCombination(double a1, double b1, double a2, double b2)
        {

            // the code below is split in many additions/subtractions that may
            // appear redundant. However, they should NOT be simplified, as they
            // use IEEE754 floating point arithmetic rounding properties.
            // as an example, the expression "ca1 - (ca1 - a1)" is NOT the same as "a1"
            // The variable naming conventions are that xyzHigh contains the most significant
            // bits of xyz and xyzLow contains its least significant bits. So theoretically
            // xyz is the sum xyzHigh + xyzLow, but in many cases below, this sum cannot
            // be represented in only one double precision number so we preserve two numbers
            // to hold it as long as we can, combining the high and low order bits together
            // only at the end, after cancellation may have occurred on high order bits

            // split a1 and b1 as two 26 bits numbers
            double ca1 = SPLIT_FACTOR * a1;
            double a1High = ca1 - (ca1 - a1);
            double a1Low = a1 - a1High;
            double cb1 = SPLIT_FACTOR * b1;
            double b1High = cb1 - (cb1 - b1);
            double b1Low = b1 - b1High;

            // accurate multiplication a1 * b1
            double prod1High = a1 * b1;
            double prod1Low = a1Low * b1Low - (((prod1High - a1High * b1High) - a1Low * b1High) - a1High * b1Low);

            // split a2 and b2 as two 26 bits numbers
            double ca2 = SPLIT_FACTOR * a2;
            double a2High = ca2 - (ca2 - a2);
            double a2Low = a2 - a2High;
            double cb2 = SPLIT_FACTOR * b2;
            double b2High = cb2 - (cb2 - b2);
            double b2Low = b2 - b2High;

            // accurate multiplication a2 * b2
            double prod2High = a2 * b2;
            double prod2Low = a2Low * b2Low - (((prod2High - a2High * b2High) - a2Low * b2High) - a2High * b2Low);

            // accurate addition a1 * b1 + a2 * b2
            double s12High = prod1High + prod2High;
            double s12Prime = s12High - prod2High;
            double s12Low = (prod2High - (s12High - s12Prime)) + (prod1High - s12Prime);

            // final rounding, s12 may have suffered many cancellations, we try
            // to recover some bits from the extra words we have saved up to now
            double result = s12High + (prod1Low + prod2Low + s12Low);

            if (Double.IsNaN(result))
            {
                // either we have split infinite numbers or some coefficients were NaNs,
                // just rely on the naive implementation and let IEEE754 handle this
                result = a1 * b1 + a2 * b2;
            }

            return result;
        }

        /// <summary>
        /// Compute a linear combination accurately.
        /// <para>
        /// This method computes a_1 &times; b_1 +
        /// a_2 &times; b_2 + a_3 &times; b_3
        /// to high accuracy. It does so by using specific multiplication and
        /// addition algorithms to preserve accuracy and reduce cancellation effects.
        /// It is based on the 2005 paper <a
        /// href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
        /// Accurate Sum and Dot Product</a> by Takeshi Ogita,
        /// Siegfried M. Rump, and Shin'ichi Oishi published in SIAM J. Sci. Comput.
        /// </para>
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <param name="a3">first factor of the third term</param>
        /// <param name="b3">second factor of the third term</param>
        /// <returns>a_1 &times; b_1 +
        /// a_2 &times; b_2 + a_3 &times; b_3</returns>
        /// <remarks>See <see cref="linearCombination(double, double, double, double)"/><para/>
        /// <see cref="linearCombination(double, double, double, double, double, double, double, double)"/>
        /// </remarks>
        public static double linearCombination(double a1, double b1, double a2, double b2, double a3, double b3)
        {

            // the code below is split in many additions/subtractions that may
            // appear redundant. However, they should NOT be simplified, as they
            // do use IEEE754 floating point arithmetic rounding properties.
            // as an example, the expression "ca1 - (ca1 - a1)" is NOT the same as "a1"
            // The variables naming conventions are that xyzHigh contains the most significant
            // bits of xyz and xyzLow contains its least significant bits. So theoretically
            // xyz is the sum xyzHigh + xyzLow, but in many cases below, this sum cannot
            // be represented in only one double precision number so we preserve two numbers
            // to hold it as long as we can, combining the high and low order bits together
            // only at the end, after cancellation may have occurred on high order bits

            // split a1 and b1 as two 26 bits numbers
            double ca1 = SPLIT_FACTOR * a1;
            double a1High = ca1 - (ca1 - a1);
            double a1Low = a1 - a1High;
            double cb1 = SPLIT_FACTOR * b1;
            double b1High = cb1 - (cb1 - b1);
            double b1Low = b1 - b1High;

            // accurate multiplication a1 * b1
            double prod1High = a1 * b1;
            double prod1Low = a1Low * b1Low - (((prod1High - a1High * b1High) - a1Low * b1High) - a1High * b1Low);

            // split a2 and b2 as two 26 bits numbers
            double ca2 = SPLIT_FACTOR * a2;
            double a2High = ca2 - (ca2 - a2);
            double a2Low = a2 - a2High;
            double cb2 = SPLIT_FACTOR * b2;
            double b2High = cb2 - (cb2 - b2);
            double b2Low = b2 - b2High;

            // accurate multiplication a2 * b2
            double prod2High = a2 * b2;
            double prod2Low = a2Low * b2Low - (((prod2High - a2High * b2High) - a2Low * b2High) - a2High * b2Low);

            // split a3 and b3 as two 26 bits numbers
            double ca3 = SPLIT_FACTOR * a3;
            double a3High = ca3 - (ca3 - a3);
            double a3Low = a3 - a3High;
            double cb3 = SPLIT_FACTOR * b3;
            double b3High = cb3 - (cb3 - b3);
            double b3Low = b3 - b3High;

            // accurate multiplication a3 * b3
            double prod3High = a3 * b3;
            double prod3Low = a3Low * b3Low - (((prod3High - a3High * b3High) - a3Low * b3High) - a3High * b3Low);

            // accurate addition a1 * b1 + a2 * b2
            double s12High = prod1High + prod2High;
            double s12Prime = s12High - prod2High;
            double s12Low = (prod2High - (s12High - s12Prime)) + (prod1High - s12Prime);

            // accurate addition a1 * b1 + a2 * b2 + a3 * b3
            double s123High = s12High + prod3High;
            double s123Prime = s123High - prod3High;
            double s123Low = (prod3High - (s123High - s123Prime)) + (s12High - s123Prime);

            // final rounding, s123 may have suffered many cancellations, we try
            // to recover some bits from the extra words we have saved up to now
            double result = s123High + (prod1Low + prod2Low + prod3Low + s12Low + s123Low);

            if (Double.IsNaN(result))
            {
                // either we have split infinite numbers or some coefficients were NaNs,
                // just rely on the naive implementation and let IEEE754 handle this
                result = a1 * b1 + a2 * b2 + a3 * b3;
            }

            return result;
        }

        /// <summary>
        /// Compute a linear combination accurately.
        /// <para>
        /// This method computes a_1 &times; b_1 +
        /// a_2 &times; b_2 + a_3 &times; b_3 +
        /// a_4 &times; b_4
        /// to high accuracy. It does so by using specific multiplication and
        /// addition algorithms to preserve accuracy and reduce cancellation effects.
        /// It is based on the 2005 paper <a
        /// href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
        /// Accurate Sum and Dot Product</a> by Takeshi Ogita,
        /// Siegfried M. Rump, and Shin'ichi Oishi published in SIAM J. Sci. Comput.
        /// </para>
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <param name="a3">first factor of the third term</param>
        /// <param name="b3">second factor of the third term</param>
        /// <param name="a4">first factor of the third term</param>
        /// <param name="b4">second factor of the third term</param>
        /// <returns>a_1 &times; b_1 +
        /// a_2 &times; b_2 + a_3 &times; b_3 +
        /// a_4 &times; b_4</returns>
        /// <remarks>
        /// See <see cref="linearCombination(double, double, double, double)"/><para>
        /// See <see cref="linearCombination(double, double, double, double, double, double)"/>
        /// </remarks>
        public static double linearCombination(double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4)
        {

            // the code below is split in many additions/subtractions that may
            // appear redundant. However, they should NOT be simplified, as they
            // do use IEEE754 floating point arithmetic rounding properties.
            // as an example, the expression "ca1 - (ca1 - a1)" is NOT the same as "a1"
            // The variables naming conventions are that xyzHigh contains the most significant
            // bits of xyz and xyzLow contains its least significant bits. So theoretically
            // xyz is the sum xyzHigh + xyzLow, but in many cases below, this sum cannot
            // be represented in only one double precision number so we preserve two numbers
            // to hold it as long as we can, combining the high and low order bits together
            // only at the end, after cancellation may have occurred on high order bits

            // split a1 and b1 as two 26 bits numbers
            double ca1 = SPLIT_FACTOR * a1;
            double a1High = ca1 - (ca1 - a1);
            double a1Low = a1 - a1High;
            double cb1 = SPLIT_FACTOR * b1;
            double b1High = cb1 - (cb1 - b1);
            double b1Low = b1 - b1High;

            // accurate multiplication a1 * b1
            double prod1High = a1 * b1;
            double prod1Low = a1Low * b1Low - (((prod1High - a1High * b1High) - a1Low * b1High) - a1High * b1Low);

            // split a2 and b2 as two 26 bits numbers
            double ca2 = SPLIT_FACTOR * a2;
            double a2High = ca2 - (ca2 - a2);
            double a2Low = a2 - a2High;
            double cb2 = SPLIT_FACTOR * b2;
            double b2High = cb2 - (cb2 - b2);
            double b2Low = b2 - b2High;

            // accurate multiplication a2 * b2
            double prod2High = a2 * b2;
            double prod2Low = a2Low * b2Low - (((prod2High - a2High * b2High) - a2Low * b2High) - a2High * b2Low);

            // split a3 and b3 as two 26 bits numbers
            double ca3 = SPLIT_FACTOR * a3;
            double a3High = ca3 - (ca3 - a3);
            double a3Low = a3 - a3High;
            double cb3 = SPLIT_FACTOR * b3;
            double b3High = cb3 - (cb3 - b3);
            double b3Low = b3 - b3High;

            // accurate multiplication a3 * b3
            double prod3High = a3 * b3;
            double prod3Low = a3Low * b3Low - (((prod3High - a3High * b3High) - a3Low * b3High) - a3High * b3Low);

            // split a4 and b4 as two 26 bits numbers
            double ca4 = SPLIT_FACTOR * a4;
            double a4High = ca4 - (ca4 - a4);
            double a4Low = a4 - a4High;
            double cb4 = SPLIT_FACTOR * b4;
            double b4High = cb4 - (cb4 - b4);
            double b4Low = b4 - b4High;

            // accurate multiplication a4 * b4
            double prod4High = a4 * b4;
            double prod4Low = a4Low * b4Low - (((prod4High - a4High * b4High) - a4Low * b4High) - a4High * b4Low);

            // accurate addition a1 * b1 + a2 * b2
            double s12High = prod1High + prod2High;
            double s12Prime = s12High - prod2High;
            double s12Low = (prod2High - (s12High - s12Prime)) + (prod1High - s12Prime);

            // accurate addition a1 * b1 + a2 * b2 + a3 * b3
            double s123High = s12High + prod3High;
            double s123Prime = s123High - prod3High;
            double s123Low = (prod3High - (s123High - s123Prime)) + (s12High - s123Prime);

            // accurate addition a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
            double s1234High = s123High + prod4High;
            double s1234Prime = s1234High - prod4High;
            double s1234Low = (prod4High - (s1234High - s1234Prime)) + (s123High - s1234Prime);

            // final rounding, s1234 may have suffered many cancellations, we try
            // to recover some bits from the extra words we have saved up to now
            double result = s1234High + (prod1Low + prod2Low + prod3Low + prod4Low + s12Low + s123Low + s1234Low);

            if (Double.IsNaN(result))
            {
                // either we have split infinite numbers or some coefficients were NaNs,
                // just rely on the naive implementation and let IEEE754 handle this
                result = a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
            }

            return result;
        }

        /// <summary>
        /// Returns true iff both arguments are null or have same dimensions and all
        /// their elements are equal as defined by
        /// <see cref="Precision.equals(float,float)"/>.
        /// </summary>
        /// <param name="x">first array</param>
        /// <param name="y">second array</param>
        /// <returns>true if the values are both null or have same dimension
        /// and equal elements.</returns>
        public static Boolean equals(float[] x, float[] y)
        {
            if ((x == null) || (y == null))
            {
                return !((x == null) ^ (y == null));
            }
            if (x.Length != y.Length)
            {
                return false;
            }
            for (int i = 0; i < x.Length; ++i)
            {
                if (!Precision.equals(x[i], y[i]))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Returns true iff both arguments are null or have same dimensions and all
        /// their elements are equal as defined by
        /// <see cref="Precision.equalsIncludingNaN(double,double)"/>.
        /// </summary>
        /// <param name="x">first array</param>
        /// <param name="y">second array</param>
        /// <returns>true if the values are both null or have same dimension and
        /// equal elements</returns>
        public static Boolean equalsIncludingNaN(float[] x, float[] y)
        {
            if ((x == null) || (y == null))
            {
                return !((x == null) ^ (y == null));
            }
            if (x.Length != y.Length)
            {
                return false;
            }
            for (int i = 0; i < x.Length; ++i)
            {
                if (!Precision.equalsIncludingNaN(x[i], y[i]))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Returns <c>true</c> iff both arguments are <c>null</c> or have same
        /// dimensions and all their elements are equal as defined by
        /// <see cref="Precision.equals(double,double)"/>.
        /// </summary>
        /// <param name="x">First array.</param>
        /// <param name="y">Second array.</param>
        /// <returns><c>true</c> if the values are both <c>null</c> or have same
        /// dimension and equal elements.</returns>
        public static Boolean equals(double[] x, double[] y)
        {
            if ((x == null) || (y == null))
            {
                return !((x == null) ^ (y == null));
            }
            if (x.Length != y.Length)
            {
                return false;
            }
            for (int i = 0; i < x.Length; ++i)
            {
                if (!Precision.equals(x[i], y[i]))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Returns <c>true</c> iff both arguments are <c>null</c> or have same
        /// dimensions and all their elements are equal as defined by
        /// <see cref="Precision.equalsIncludingNaN(double,double)"/>.
        /// </summary>
        /// <param name="x">First array.</param>
        /// <param name="y">Second array.</param>
        /// <returns><c>true</c> if the values are both <c>null</c> or have same
        /// dimension and equal elements.</returns>
        public static Boolean equalsIncludingNaN(double[] x, double[] y)
        {
            if ((x == null) || (y == null))
            {
                return !((x == null) ^ (y == null));
            }
            if (x.Length != y.Length)
            {
                return false;
            }
            for (int i = 0; i < x.Length; ++i)
            {
                if (!Precision.equalsIncludingNaN(x[i], y[i]))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Normalizes an array to make it sum to a specified value.
        /// Returns the result of the transformation
        /// <code>
        ///    x |-> x * normalizedSum / sum
        /// </code>
        /// applied to each non-NaN element x of the input array, where sum is the
        /// sum of the non-NaN entries in the input array.
        /// <para>
        /// Throws IllegalArgumentException if <c>normalizedSum</c> is infinite
        /// or NaN and ArithmeticException if the input array contains any infinite elements
        /// or sums to 0.
        /// </para>
        /// <para>
        /// Ignores (i.e., copies unchanged to the output array) NaNs in the input array.
        /// </para>
        /// </summary>
        /// <param name="values">Input array to be normalized</param>
        /// <param name="normalizedSum">Target sum for the normalized array</param>
        /// <returns>the normalized array.</returns>
        /// <exception cref="MathArithmeticException"> if the input array contains infinite
        /// elements or sums to zero.</exception>
        /// <exception cref="MathIllegalArgumentException"> if the target sum is infinite 
        /// or <c>NaN</c>.</exception>
        public static double[] normalizeArray(double[] values, double normalizedSum)
        {
            if (Double.IsInfinity(normalizedSum))
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("NORMALIZE_INFINITE"));
            }
            if (Double.IsNaN(normalizedSum))
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("NORMALIZE_NAN"));
            }
            double sum = 0d;
            int len = values.Length;
            double[] outp = new double[len];
            for (int i = 0; i < len; i++)
            {
                if (Double.IsInfinity(values[i]))
                {
                    throw new MathIllegalArgumentException(new LocalizedFormats("INFINITE_ARRAY_ELEMENT"), values[i], i);
                }
                if (!Double.IsNaN(values[i]))
                {
                    sum += values[i];
                }
            }
            if (sum == 0)
            {
                throw new MathArithmeticException(new LocalizedFormats("ARRAY_SUMS_TO_ZERO"));
            }
            for (int i = 0; i < len; i++)
            {
                if (Double.IsNaN(values[i]))
                {
                    outp[i] = Double.NaN;
                }
                else
                {
                    outp[i] = values[i] * normalizedSum / sum;
                }
            }
            return outp;
        }

        /// <summary>
        /// Build an array of elements.
        /// <para>
        /// Arrays are filled with field.getZero()
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="field">field to which array elements belong</param>
        /// <param name="length">length of the array</param>
        /// <returns>a new array</returns>
        public static T[] buildArray<T>(Field<T> field, int length)
        {
            //@SuppressWarnings("unchecked") // OK because field must be correct class
            T[] array = (T[])Array.CreateInstance(field.getRuntimeClass(), length);
            for (int i = 0; i < array.Length; ++i)
            {
                array[i] = field.getZero();
            }
            return array;
        }

        /// <summary>
        /// Build a double dimension  array of elements.
        /// <para>
        /// Arrays are filled with field.getZero()
        /// </summary>
        /// <typeparam name="T">the type of the field elements</typeparam>
        /// <param name="field">field to which array elements belong</param>
        /// <param name="rows">number of rows in the array</param>
        /// <param name="columns">number of columns (may be negative to build partial
        /// arrays in the same way <c>new Field[rows][]</c> works)</param>
        /// <returns>a new array</returns>
        public static T[][] buildArray<T>(Field<T> field, int rows, int columns)
        {
            T[][] array;
            if (columns < 0)
            {
                T[] dummyRow = buildArray(field, 0);
                array = (T[][])Array.CreateInstance(dummyRow.GetType(), rows);
            }
            else
            {
                array = (T[][])Array.CreateInstance(field.getRuntimeClass(), new int[] { rows, columns });
                for (int i = 0; i < rows; ++i)
                {
                    for (int j = 0; j < array.Length; ++j)
                    {
                        array[i][j] = field.getZero();
                    }
                }
            }
            return array;
        }

        /// <summary>
        /// Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
        /// convolution</a> between two sequences.
        /// <para>
        /// The solution is obtained via straightforward computation of the
        /// convolution sum (and not via FFT). Whenever the computation needs
        /// an element that would be located at an index outside the input arrays,
        /// the value is assumed to be zero.
        /// </para>
        /// </summary>
        /// <param name="x">First sequence.
        /// Typically, this sequence will represent an input signal to a system.</param>
        /// <param name="h">Second sequence.
        /// Typically, this sequence will represent the impulse response of the system.</param>
        /// <returns>the convolution of <c>x</c> and <c>h</c>.
        /// This array's length will be <c>x.length + h.length - 1</c>.</returns>
        /// <exception cref="NullArgumentException"> if either <c>x</c> or <c>h</c> is 
        /// <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if either <c>x</c> or <c>h</c> is empty.</exception>
        public static double[] convolve(double[] x, double[] h)
        {
            MathUtils.checkNotNull(x);
            MathUtils.checkNotNull(h);

            int xLen = x.Length;
            int hLen = h.Length;

            if (xLen == 0 || hLen == 0)
            {
                throw new NoDataException();
            }

            // initialize the output array
            int totalLength = xLen + hLen - 1;
            double[] y = new double[totalLength];

            // straightforward implementation of the convolution sum
            for (int n = 0; n < totalLength; n++)
            {
                double yn = 0;
                int k = FastMath.max(0, n + 1 - xLen);
                int j = n - k;
                while (k < hLen && j >= 0)
                {
                    yn += x[j--] * h[k++];
                }
                y[n] = yn;
            }

            return y;
        }

        /// <summary>
        /// Specification for indicating that some operation applies
        /// before or after a given index.
        /// </summary>
        public enum Position
        {
            /// <summary>
            /// Designates the beginning of the array (near index 0).
            /// </summary>
            HEAD,

            /// <summary>
            /// Designates the end of the array.
            /// </summary>
            TAIL
        }

        /// <summary>
        /// Shuffle the entries of the given array.
        /// The <c>start</c> and <c>pos</c> parameters select which portion
        /// of the array is randomized and which is left untouched.
        /// </summary>
        /// <param name="list">Array whose entries will be shuffled (in-place).</param>
        /// <param name="start">Index at which shuffling begins.</param>
        /// <param name="pos">Shuffling is performed for index positions between
        /// <c>start</c> and either the end (if <see cref="Position.TAIL"/>)
        /// or the beginning (if <see cref="Position.HEAD"/>) of the array.</param>
        /// <remarks>See <see cref="shuffle(int[],int,Position,RandomGenerator)"/></remarks>
        public static void shuffle(int[] list, int start, Position pos)
        {
            shuffle(list, start, pos, new Well19937c());
        }

        /// <summary>
        /// Shuffle the entries of the given array, using the
        /// <a href="http://en.wikipedia.org/wiki/Fisher–Yates_shuffle#The_modern_algorithm">
        /// Fisher–Yates</a> algorithm.
        /// The <c>start</c> and <c>pos</c> parameters select which portion
        /// of the array is randomized and which is left untouched.
        /// </summary>
        /// <param name="list">Array whose entries will be shuffled (in-place).</param>
        /// <param name="start">Index at which shuffling begins.</param>
        /// <param name="pos">Shuffling is performed for index positions between
        /// <c>start</c> and either the end (if <see cref="Position.TAIL"/>)
        /// or the beginning (if <see cref="Position.HEAD"/>) of the array.</param>
        /// <param name="rng">Random number generator.</param>
        public static void shuffle(int[] list, int start, Position pos, RandomGenerator rng)
        {
            switch (pos)
            {
                case Position.TAIL:
                    {
                        for (int i = list.Length - 1; i >= start; i--)
                        {
                            int target;
                            if (i == start)
                            {
                                target = start;
                            }
                            else
                            {
                                // NumberIsTooLargeException cannot occur.
                                target = new UniformIntegerDistribution(rng, start, i).sample();
                            }
                            int temp = list[target];
                            list[target] = list[i];
                            list[i] = temp;
                        }
                    }
                    break;
                case Position.HEAD:
                    {
                        for (int i = 0; i <= start; i++)
                        {
                            int target;
                            if (i == start)
                            {
                                target = start;
                            }
                            else
                            {
                                // NumberIsTooLargeException cannot occur.
                                target = new UniformIntegerDistribution(rng, i, start).sample();
                            }
                            int temp = list[target];
                            list[target] = list[i];
                            list[i] = temp;
                        }
                    }
                    break;
                default:
                    throw new MathInternalError(); // Should never happen.
            }
        }

        /// <summary>
        /// Shuffle the entries of the given array.
        /// </summary>
        /// <param name="list">Array whose entries will be shuffled (in-place).</param>
        /// <param name="rng">Random number generator.</param>
        /// <remarks>See <see cref="shuffle(int[],int,Position,RandomGenerator)"/></remarks>
        public static void shuffle(int[] list, RandomGenerator rng)
        {
            shuffle(list, 0, Position.TAIL, rng);
        }

        /// <summary>
        /// Shuffle the entries of the given array.
        /// </summary>
        /// <param name="list">Array whose entries will be shuffled (in-place).</param>
        /// <remarks>See <see cref="shuffle(int[],int,Position,RandomGenerator)"/></remarks>
        public static void shuffle(int[] list)
        {
            shuffle(list, new Well19937c());
        }

        /// <summary>
        /// Returns an array representing the natural number <c>n</c>.
        /// </summary>
        /// <param name="n">Natural number.</param>
        /// <returns>an array whose entries are the numbers 0, 1, ..., <c>n</c>-1.
        /// If <c>n == 0</c>, the returned array is empty.</returns>
        public static int[] natural(int n)
        {
            return sequence(n, 0, 1);
        }

        /// <summary>
        /// Returns an array of <c>size</c> integers starting at <c>start</c>,
        /// skipping <c>stride</c> numbers.
        /// </summary>
        /// <param name="size">Natural number.</param>
        /// <param name="start">Natural number.</param>
        /// <param name="stride">Natural number.</param>
        /// <returns>an array whose entries are the numbers
        /// <c>start, start + stride, ..., start + (size - 1) * stride</c>.
        /// If <c>size == 0</c>, the returned array is empty.</returns>
        public static int[] sequence(int size, int start, int stride)
        {
            int[] a = new int[size];
            for (int i = 0; i < size; i++)
            {
                a[i] = start + i * stride;
            }
            return a;
        }

        /// <summary>
        /// This method is used
        /// to verify that the input parameters designate a subarray of positive length.
        /// <para>
        /// <list type="bullet">
        /// <item>returns <c>true</c> iff the parameters designate a subarray of
        /// positive length</item>
        /// <item>throws <c>MathIllegalArgumentException</c> if the array is null or
        /// or the indices are invalid</item>
        /// <item>returns <c>false</c> if the array is non-null, but
        /// <c>length</c> is 0.</item>
        /// </list></para>
        /// </summary>
        /// <param name="values">the input array</param>
        /// <param name="begin">index of the first array element to include</param>
        /// <param name="length">the number of elements to include</param>
        /// <returns>true if the parameters are valid and designate a subarray of positive length
        /// </returns>
        /// <exception cref="MathIllegalArgumentException"> if the indices are invalid or
        /// the array is null</exception>
        public static Boolean verifyValues(double[] values, int begin, int length)
        {
            return verifyValues(values, begin, length, false);
        }

        /// <summary>
        /// This method is used
        /// to verify that the input parameters designate a subarray of positive length.
        /// <para>
        /// <list type="bullet">
        /// <item>returns <c>true</c> iff the parameters designate a subarray of
        /// non-negative length</item>
        /// <item>throws <c>IllegalArgumentException</c> if the array is null or
        /// or the indices are invalid</item>
        /// <item>returns <c>false</c> if the array is non-null, but
        /// <c>length</c> is 0 unless <c>allowEmpty</c> is <c>true</c></item>
        /// </list></para>
        /// </summary>
        /// <param name="values">the input array</param>
        /// <param name="begin">index of the first array element to include</param>
        /// <param name="length">the number of elements to include</param>
        /// <param name="allowEmpty">if <c>true</c> then zero length arrays are allowed</param>
        /// <returns>true if the parameters are valid</returns>
        /// <exception cref="MathIllegalArgumentException"> if the indices are invalid or
        /// the array is null</exception>
        public static Boolean verifyValues(double[] values, int begin, int length, Boolean allowEmpty)
        {
            if (values == null)
            {
                throw new NullArgumentException(new LocalizedFormats("INPUT_ARRAY"));
            }

            if (begin < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("START_POSITION"), begin);
            }

            if (length < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("LENGTH"), length);
            }

            if (begin + length > values.Length)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(new LocalizedFormats("SUBARRAY_ENDS_AFTER_ARRAY_END"), begin + length, values.Length, true);
            }

            if (length == 0 && !allowEmpty)
            {
                return false;
            }

            return true;

        }

        /// <summary>
        /// This method is used
        /// to verify that the begin and length parameters designate a subarray of positive length
        /// and the weights are all non-negative, non-NaN, finite, and not all zero.
        /// <para>
        /// <list type="bullet">
        /// <item>returns <c>true</c> iff the parameters designate a subarray of
        /// positive length and the weights array contains legitimate values.</item>
        /// <item>throws <c>IllegalArgumentException</c> if any of the following are true:
        /// <list type="bullet"><item>the values array is null</item>
        /// <item>the weights array is null</item>
        /// <item>the weights array does not have the same length as the values array</item>
        /// <item>the weights array contains one or more infinite values</item>
        /// <item>the weights array contains one or more NaN values</item>
        /// <item>the weights array contains negative values</item>
        /// <item>the start and length arguments do not determine a valid array</item></list>
        /// </item>
        /// <item>returns <c>false</c> if the array is non-null, but
        /// <c>length</c> is 0.</item>
        /// </list></para>
        /// </summary>
        /// <param name="values">the input array</param>
        /// <param name="weights">the weights array</param>
        /// <param name="begin">index of the first array element to include</param>
        /// <param name="length">the number of elements to include</param>
        /// <returns>true if the parameters are valid and designate a subarray of positive 
        /// length</returns>
        /// <exception cref="MathIllegalArgumentException"> if the indices are invalid
        /// or the array is null</exception>
        public static Boolean verifyValues(double[] values, double[] weights, int begin, int length)
        {
            return verifyValues(values, weights, begin, length, false);
        }

        /// <summary>
        /// This method is used
        /// to verify that the begin and length parameters designate a subarray of positive length
        /// and the weights are all non-negative, non-NaN, finite, and not all zero.
        /// <para>
        /// <list type="bullet">
        /// <item>returns <c>true</c> iff the parameters designate a subarray of
        /// non-negative length and the weights array contains legitimate values.</item>
        /// <item>throws <c>MathIllegalArgumentException</c> if any of the following are true:
        /// <list type="bullet"><item>the values array is null</item>
        /// <item>the weights array is null</item>
        /// <item>the weights array does not have the same length as the values array</item>
        /// <item>the weights array contains one or more infinite values</item>
        /// <item>the weights array contains one or more NaN values</item>
        /// <item>the weights array contains negative values</item>
        /// <item>the start and length arguments do not determine a valid array</item>
        /// </list>
        /// </item>
        /// <item>returns <c>false</c> if the array is non-null, but
        /// <c>length</c> is 0 unless <c>allowEmpty</c> is <c>true</c>.</item>
        /// </list></para>
        /// </summary>
        /// <param name="values">the input array.</param>
        /// <param name="weights">the weights array.</param>
        /// <param name="begin">index of the first array element to include.</param>
        /// <param name="length">the number of elements to include.</param>
        /// <param name="allowEmpty">if <c>true</c> than allow zero length arrays to pass.</param>
        /// <returns><c>true</c> if the parameters are valid.</returns>
        /// <exception cref="NullArgumentException"> if either of the arrays are null</exception>
        /// <exception cref="MathIllegalArgumentException"> if the array indices are not valid,
        /// the weights array contains NaN, infinite or negative elements, or there
        /// are no positive weights.</exception>
        public static Boolean verifyValues(double[] values, double[] weights, int begin, int length, Boolean allowEmpty)
        {

            if (weights == null || values == null)
            {
                throw new NullArgumentException(new LocalizedFormats("INPUT_ARRAY"));
            }

            if (weights.Length != values.Length)
            {
                throw new DimensionMismatchException(weights.Length, values.Length);
            }

            Boolean containsPositiveWeight = false;
            for (int i = begin; i < begin + length; i++)
            {
                double weight = weights[i];
                if (Double.IsNaN(weight))
                {
                    throw new MathIllegalArgumentException(new LocalizedFormats("NAN_ELEMENT_AT_INDEX"), i);
                }
                if (Double.IsInfinity(weight))
                {
                    throw new MathIllegalArgumentException(new LocalizedFormats("INFINITE_ARRAY_ELEMENT"), weight, i);
                }
                if (weight < 0)
                {
                    throw new MathIllegalArgumentException(new LocalizedFormats("NEGATIVE_ELEMENT_AT_INDEX"), i, weight);
                }
                if (!containsPositiveWeight && weight > 0.0)
                {
                    containsPositiveWeight = true;
                }
            }

            if (!containsPositiveWeight)
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("WEIGHT_AT_LEAST_ONE_NON_ZERO"));
            }
            return verifyValues(values, begin, length, allowEmpty);
        }
    }
}