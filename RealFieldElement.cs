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
namespace Math3
{
    /// <summary>
    /// Interface representing a <a href="http://mathworld.wolfram.com/RealNumber.html">real</a>
    /// <a href="http://mathworld.wolfram.com/Field.html">field</a>.
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    /// <remarks>
    /// See <see cref="FieldElement"/>
    /// </remarks>
    public interface RealFieldElement<T> : FieldElement<T>
    {
        /// <summary>
        /// Get the real value of the number.
        /// </summary>
        /// <returns>real value</returns>
        double getReal();

        /// <summary>
        /// '+' operator.
        /// </summary>
        /// <param name="a">right hand side parameter of the operator</param>
        /// <returns>this+a</returns>
        T add(double a);

        /// <summary>
        /// '-' operator.
        /// </summary>
        /// <param name="a">right hand side parameter of the operator</param>
        /// <returns>this-a</returns>
        T subtract(double a);

        /// <summary>
        /// '&times;' operator.
        /// </summary>
        /// <param name="a">right hand side parameter of the operator</param>
        /// <returns>this&times;a</returns>
        T multiply(double a);

        /// <summary>
        /// '&divide;' operator.
        /// </summary>
        /// <param name="a">right hand side parameter of the operator</param>
        /// <returns>this&divides;a</returns>
        T divide(double a);

        /// <summary>
        /// IEEE remainder operator.
        /// </summary>
        /// <param name="a">right hand side parameter of the operator</param>
        /// <returns>this - n &times; a where n is the closest integer to this/a
        /// (the even integer is chosen for n if this/a is halfway between two integers)</returns>
        T remainder(double a);

        /// <summary>
        /// IEEE remainder operator.
        /// </summary>
        /// <param name="a">right hand side parameter of the operator</param>
        /// <returns>this - n &times; a where n is the closest integer to this/a
        /// (the even integer is chosen for n if this/a is halfway between two integers)</returns>
        /// <exception cref="DimensionMismatchException"> if number of free parameters or orders are inconsistent</exception>
        T remainder(T a);

        /// <summary>
        /// absolute value.
        /// </summary>
        /// <returns>abs(this)</returns>
        T abs();

        /// <summary>
        /// Get the smallest whole number larger than instance.
        /// </summary>
        /// <returns>ceil(this)</returns>
        T ceil();

        /// <summary>
        /// Get the largest whole number smaller than instance. 
        /// </summary>
        /// <returns>floor(this)</returns>
        T floor();

        /// <summary>
        /// Get the whole number that is the nearest to the instance, or the even one if x
        /// is exactly half way between two integers.
        /// </summary>
        /// <returns>a double number r such that r is an integer r - 0.5 <= this <= r + 0.5</returns>
        T rint();

        /// <summary>
        /// Get the closest long to instance value.
        /// </summary>
        /// <returns>closest long to <see cref="getReal()"/></returns>
        long round();

        /// <summary>
        /// Compute the signum of the instance.
        /// The signum is -1 for negative numbers, +1 for positive numbers and 0 otherwise
        /// </summary>
        /// <returns>-1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a</returns>
        T signum();

        /// <summary>
        /// Returns the instance with the sign of the argument.
        /// A NaN <c>sign</c> argument is treated as positive.
        /// </summary>
        /// <param name="sign">the sign for the returned value</param>
        /// <returns>the instance with the same sign as the <c>sign</c> argument</returns>
        T copySign(T sign);

        /// <summary>
        /// Returns the instance with the sign of the argument.
        /// A NaN <c>sign</c> argument is treated as positive.
        /// </summary>
        /// <param name="sign">the sign for the returned value</param>
        /// <returns>the instance with the same sign as the <c>sign</c> argument</returns>
        T copySign(double sign);

        /// <summary>
        /// Multiply the instance by a power of 2.
        /// </summary>
        /// <param name="n">power of 2</param>
        /// <returns>this &times; 2^n</returns>
        T scalb(int n);

        /// <summary>
        /// Returns the hypotenuse of a triangle with sides <c>this</c> and <c>y</c>
        /// - sqrt(this^2&nbsp;+y^2)<para/>
        /// avoiding intermediate overflow or underflow.
        /// <list type="bullet">
        /// <item> If either argument is infinite, then the result is positive infinity.</item>
        /// <item> else, if either argument is NaN then the result is NaN.</item>
        /// </list>
        /// </summary>
        /// <param name="y">a value</param>
        /// <returns>sqrt(this^2&nbsp;+y^2)</returns>
        /// <exception cref="DimensionMismatchException"> if number of free parameters or orders are inconsistent</exception>
        T hypot(T y);

        /// <summary>
        /// Square root.
        /// </summary>
        /// <returns>square root of the instance</returns>
        T sqrt();

        /// <summary>
        /// Cubic root.
        /// </summary>
        /// <returns>cubic root of the instance</returns>
        T cbrt();

        /// <summary>
        /// N^th root.
        /// </summary>
        /// <param name="n">order of the root</param>
        /// <returns>n^th root of the instance</returns>
        T rootN(int n);

        /// <summary>
        /// Power operation.
        /// </summary>
        /// <param name="p">power to apply</param>
        /// <returns>this^p</returns>
        T pow(double p);

        /// <summary>
        /// Integer power operation. 
        /// </summary>
        /// <param name="n">power to apply</param>
        /// <returns>this^n</returns>
        T pow(int n);

        /// <summary>
        /// Power operation.
        /// </summary>
        /// <param name="e">exponent</param>
        /// <returns>this<sup>e</sup></returns>
        /// <exception cref="DimensionMismatchException"> if number of free parameters or orders are inconsistent</exception>
        T pow(T e);

        /// <summary>
        /// Exponential.
        /// </summary>
        /// <returns>exponential of the instance</returns>
        T exp();

        /// <summary>
        /// Exponential minus 1.
        /// </summary>
        /// <returns>exponential minus one of the instance</returns>
        T expm1();

        /// <summary>
        /// Natural logarithm.
        /// </summary>
        /// <returns>logarithm of the instance</returns>
        T log();

        /// <summary>
        /// Shifted natural logarithm.
        /// </summary>
        /// <returns>logarithm of one plus the instance</returns>
        T log1p();

        //    TODO: add this method in 4.0, as it is not possible to do it in 3.2
        //          due to incompatibility of the return type in the Dfp class
        //    /* Base 10 logarithm.
        //     * @return base 10 logarithm of the instance
        //     */
        //    T log10();

        /// <summary>
        /// Cosine operation. 
        /// </summary>
        /// <returns>cos(this)</returns>
        T cos();

        /// <summary>
        /// Sine operation.
        /// </summary>
        /// <returns>sin(this)</returns>
        T sin();

        /// <summary>
        /// Tangent operation.
        /// </summary>
        /// <returns>tan(this)</returns>
        T tan();

        /// <summary>
        /// Arc cosine operation.
        /// </summary>
        /// <returns>acos(this)</returns>
        T acos();

        /// <summary>
        /// Arc sine operation.
        /// </summary>
        /// <returns>asin(this)</returns>
        T asin();

        /// <summary>
        /// Arc tangent operation.
        /// </summary>
        /// <returns>atan(this)</returns>
        T atan();

        /// <summary>
        /// Two arguments arc tangent operation.
        /// </summary>
        /// <param name="x">second argument of the arc tangent</param>
        /// <returns>atan2(this, x)</returns>
        /// <exception cref="DimensionMismatchException"> if number of free parameters or orders are inconsistent</exception>
        T atan2(T x);

        /// <summary>
        /// Hyperbolic cosine operation.
        /// </summary>
        /// <returns>cosh(this)</returns>
        T cosh();

        /// <summary>
        /// Hyperbolic sine operation.
        /// </summary>
        /// <returns>sinh(this)</returns>
        T sinh();

        /// <summary>
        /// Hyperbolic tangent operation.
        /// </summary>
        /// <returns>tanh(this)</returns>
        T tanh();

        /// <summary>
        /// Inverse hyperbolic cosine operation.
        /// </summary>
        /// <returns>acosh(this)</returns>
        T acosh();

        /// <summary>
        /// Inverse hyperbolic sine operation.
        /// </summary>
        /// <returns>asin(this)</returns>
        T asinh();

        /// <summary>
        /// Inverse hyperbolic  tangent operation.
        /// </summary>
        /// <returns>atanh(this)</returns>
        T atanh();

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a">Factors.</param>
        /// <param name="b">Factors.</param>
        /// <returns><c>&Sigma;_i a_i b_i</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if arrays dimensions don't match
        /// </exception>
        T linearCombination(T[] a, T[] b);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a">Factors.</param>
        /// <param name="b">Factors.</param>
        /// <returns><c>&Sigma;_i a_i b_i</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if arrays dimensions don't match
        /// </exception>
        T linearCombination(double[] a, T[] b);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <returns>a_1&times;b_1 +
        /// a_2&times;b_2</returns>
        /// <remarks>
        /// See <see cref="linearCombination(Object, Object, Object, Object, Object, Object)"/>
        /// <para/>
        /// See <see cref="linearCombination(Object, Object, Object, Object, Object, Object, Object, Object)"/>
        /// </remarks>
        T linearCombination(T a1, T b1, T a2, T b2);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <returns>a_1&times;b_1 +
        /// a_2_&times;b_2</returns>
        /// <remarks>
        /// See <see cref="linearCombination(double, Object, double, Object, double, Object)"/>
        /// <para/>
        /// See <see cref="linearCombination(double, Object, double, Object, double, Object, double, Object)"/>
        /// </remarks>
        T linearCombination(double a1, T b1, double a2, T b2);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <param name="a3">first factor of the third term</param>
        /// <param name="b3">second factor of the third term</param>
        /// <returns>a_1&times;b_1 +
        /// a_2&times;b_2 + a_3&times;b_3</returns>
        /// <remarks>
        /// See <see cref="linearCombination(Object, Object, Object, Object)"/><para/>
        /// See <see cref="linearCombination(Object, Object, Object, Object, Object, Object, Object, Object)"/>
        /// </remarks>
        T linearCombination(T a1, T b1, T a2, T b2, T a3, T b3);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <param name="a3">first factor of the third term</param>
        /// <param name="b3">second factor of the third term</param>
        /// <returns>a_1&times;b_1 +
        /// a_2&times;b_2 + a_3>&times;b_3</returns>
        /// <remarks>
        /// See <see cref="linearCombination(double, Object, double, Object)"/><para/>
        /// See <see cref="linearCombination(double, Object, double, Object, double, Object, double, Object)"/>
        /// </remarks>
        T linearCombination(double a1, T b1, double a2, T b2, double a3, T b3);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <param name="a3">first factor of the third term</param>
        /// <param name="b3">second factor of the third term</param>
        /// <param name="a4">first factor of the third term</param>
        /// <param name="b4">second factor of the third term</param>
        /// <returns>a_1&times;b_1 +
        /// a_2&times;b_2 + a_3&times;b_3 +
        /// a<sub>4</sub>&times;b<sub>4</sub></returns>
        /// <remarks>
        /// See <see cref="linearCombination(Object, Object, Object, Object)"/><para/>
        /// See <see cref="linearCombination(Object, Object, Object, Object, Object, Object)"/>
        /// </remarks>
        T linearCombination(T a1, T b1, T a2, T b2, T a3, T b3, T a4, T b4);

        /// <summary>
        /// Compute a linear combination.
        /// </summary>
        /// <param name="a1">first factor of the first term</param>
        /// <param name="b1">second factor of the first term</param>
        /// <param name="a2">first factor of the second term</param>
        /// <param name="b2">second factor of the second term</param>
        /// <param name="a3">first factor of the third term</param>
        /// <param name="b3">second factor of the third term</param>
        /// <param name="a4">first factor of the third term</param>
        /// <param name="b4">second factor of the third term</param>
        /// <returns>a_1&times;b_1 + 
        /// a_2&times;b_2 + a_3&times;b_3 +
        /// a_4&times;b_4</returns>
        /// <remarks>
        /// See <see cref="linearCombination(double, Object, double, Object)"/><para/>
        /// See <see cref="linearCombination(double, Object, double, Object, double, Object)"/>
        /// </remarks>
        T linearCombination(double a1, T b1, double a2, T b2, double a3, T b3, double a4, T b4);
    }
}