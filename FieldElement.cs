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
    /// Interface representing <a href="http://mathworld.wolfram.com/Field.html">field</a> elements.
    /// </summary>
    /// <remarks>See <see cref="Field"/></remarks>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public interface FieldElement<T>
    {
        /// <summary>
        /// Compute this + a. 
        /// </summary>
        /// <param name="a">element to add</param>
        /// <returns>new element representing this + a</returns>
        /// <exception cref="NullArgumentException">if <c>a</c> is <c>null</c>.</exception>
        T add(T a);

        /// <summary>
        /// Compute this - a. 
        /// </summary>
        /// <param name="a">element to subctract</param>
        /// <returns>new element representing this - a</returns>
        /// <exception cref="NullArgumentException">if <c>a</c> is <c>null</c>.</exception>
        T subtract(T a);

        /// <summary>
        /// Returns the additive inverse of <c>this</c> element.
        /// </summary>
        /// <returns>the opposite of <c>this</c>.</returns>
        T negate();

        /// <summary>
        /// Compute n &times; this. Multiplication by an integer number is defined
        /// as the following sum
        /// n &times; this = &sum;_{i=1}^n this.
        /// </summary>
        /// <param name="n">Number of times <c>this</c> must be added to itself.</param>
        /// <returns>A new element representing n &times; this.</returns>
        T multiply(int n);

        /// <summary>
        /// Compute this &times; a. 
        /// </summary>
        /// <param name="a">element to multiply</param>
        /// <returns>new element representing this &times; a</returns>
        /// <exception cref="NullArgumentException">if <c>a</c> is <c>null</c>.</exception>
        T multiply(T a);

        /// <summary>
        /// Compute this &divide; a. 
        /// </summary>
        /// <param name="a">element to add</param>
        /// <returns>a new element representing this &divide; a</returns>
        /// <exception cref="NullArgumentException">if <c>addend</c> is <c>null</c>.</exception>
        /// <exception cref="MathArithmeticException">if <c>a</c> is zero</exception>
        T divide(T a);

        /// <summary>
        /// Returns the multiplicative inverse of <c>this</c> element.
        /// </summary>
        /// <returns>the inverse of <c>this</c>.</returns>
        /// <exception cref="MathArithmeticException">if <c>this</c> is zero</exception>
        T reciprocal();

        /// <summary>
        /// Get the <see cref="Field"/> to which the instance belongs.
        /// </summary>
        /// <returns><see cref="Field"/> to which the instance belongs</returns>
        Field<T> getField();
    }
}
