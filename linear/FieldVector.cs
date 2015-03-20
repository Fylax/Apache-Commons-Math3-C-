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
using System;

namespace Math3.linear
{
    /// <summary>
    /// Interface defining a field-valued vector with basic algebraic operations.
    /// <para>
    /// vector element indexing is 0-based -- e.g., <c>getEntry(0)</code>
    /// returns the first element of the vector.
    /// </para>
    /// <para>
    /// The various <c>mapXxx</c> and <c>mapXxxToSelf</c> methods operate
    /// on vectors element-wise, i.e. they perform the same operation (adding a scalar,
    /// applying a function ...) on each element in turn. The <c>mapXxx</c>
    /// versions create a new vector to hold the result and do not change the instance.
    /// The <c>mapXxxToSelf</c> versions use the instance itself to store the
    /// results, so the instance is changed by these methods. In both cases, the result
    /// vector is returned by the methods, this allows to use the <i>fluent API</i>
    /// style, like this:
    /// </para>
    /// <code>
    ///   RealVector result = v.mapAddToSelf(3.0).mapTanToSelf().mapSquareToSelf();
    /// </code>
    /// <para>
    /// Note that as almost all operations on <see cref="FieldElement"/> throw 
    /// <see cref="NullArgumentException"/> when operating on a null element, it is 
    /// the responsibility of <c>FieldVector</c> implementations to make sure no null 
    /// elements are inserted into the vector. This must be done in all constructors and
    /// all setters.
    /// <para>
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
    public interface FieldVector<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Get the type of field elements of the vector.
        /// </summary>
        /// <returns>type of field elements of the vector</returns>
        Field<T> getField();

        /// <summary>
        /// Returns a (deep) copy of this.
        /// </summary>
        /// <returns>vector copy</returns>
        FieldVector<T> copy();

        /// <summary>
        /// Compute the sum of <c>this</c> and <c>v</c>.
        /// </summary>
        /// <param name="v">vector to be added</param>
        /// <returns><c>this + v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size 
        /// as <c>this</c></exception>
        FieldVector<T> add(FieldVector<T> v);

        /// <summary>
        /// Compute <c>this</c> minus <c>v</c>.
        /// </summary>
        /// <param name="v">vector to be subtracted</param>
        /// <returns><c>this - v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size
        /// as <c>this</c></exception>
        FieldVector<T> subtract(FieldVector<T> v);

        /// <summary>
        /// Map an addition operation to each entry.
        /// </summary>
        /// <param name="d">value to be added to each entry</param>
        /// <returns><c>this + d</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        FieldVector<T> mapAdd(T d);

        /// <summary>
        /// Map an addition operation to each entry.
        /// <para>The instance is changed by this method.</para>
        /// </summary>
        /// <param name="d">value to be added to each entry</param>
        /// <returns>for convenience, return <c>this</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        FieldVector<T> mapAddToSelf(T d);

        /// <summary>
        /// Map a subtraction operation to each entry.
        /// </summary>
        /// <param name="d">value to be subtracted to each entry</param>
        /// <returns><c>this - d</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c></exception>
        FieldVector<T> mapSubtract(T d);

        /// <summary>
        /// Map a subtraction operation to each entry.
        /// <para>The instance <strong>is</strong> changed by this method.</para>
        /// </summary>
        /// <param name="d">value to be subtracted to each entry</param>
        /// <returns>for convenience, return <c>this</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c></exception>
        FieldVector<T> mapSubtractToSelf(T d);

        /// <summary>
        /// Map a multiplication operation to each entry.
        /// </summary>
        /// <param name="d">value to multiply all entries by</param>
        /// <returns><c>this * d</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        FieldVector<T> mapMultiply(T d);

        /// <summary>
        /// Map a multiplication operation to each entry.
        /// <para>The instance <strong>is</strong> changed by this method.</para>
        /// </summary>
        /// <param name="d">value to multiply all entries by</param>
        /// <returns>for convenience, return <c>this</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        FieldVector<T> mapMultiplyToSelf(T d);

        /// <summary>
        /// Map a division operation to each entry.
        /// </summary>
        /// <param name="d">value to divide all entries by</param>
        /// <returns><c>this / d</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        FieldVector<T> mapDivide(T d);

        /// <summary>
        /// Map a division operation to each entry.
        /// <para>The instance <strong>is</strong> changed by this method.</para>
        /// </summary>
        /// <param name="d">value to divide all entries by</param>
        /// <returns>for convenience, return <c>this</c></returns>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="MathArithmeticException"> if <c>d</c> is zero.</exception>
        FieldVector<T> mapDivideToSelf(T d);

        /// <summary>
        /// Map the 1/x function to each entry.
        /// </summary>
        /// <returns>a vector containing the result of applying the function to each entry.
        /// </returns>
        /// <exception cref="MathArithmeticException"> if one of the entries is zero.</exception>
        FieldVector<T> mapInv();

        /// <summary>
        /// Map the 1/x function to each entry.
        /// <para>The instance <strong>is</strong> changed by this method.</para>
        /// </summary>
        /// <returns>for convenience, return <c>this</c></returns>
        /// <exception cref="MathArithmeticException"> if one of the entries is zero.</exception>
        FieldVector<T> mapInvToSelf();

        /// <summary>
        /// Element-by-element multiplication.
        /// </summary>
        /// <param name="v">vector by which instance elements must be multiplied</param>
        /// <returns>a vector containing <c>this[i] * v[i]</c> for all <c>i</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size
        /// as <c>this</c></exception>
        FieldVector<T> ebeMultiply(FieldVector<T> v);

        /// <summary>
        /// Element-by-element division.
        /// </summary>
        /// <param name="v">vector by which instance elements must be divided</param>
        /// <returns>a vector containing <c>this[i] / v[i]</c> for all <c>i</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size 
        /// as <c>this</c></exception>
        /// <exception cref="MathArithmeticException"> if one entry of <c>v</c> is zero.
        /// </exception>
        FieldVector<T> ebeDivide(FieldVector<T> v);

        /// <summary>
        /// Returns vector entries as a T array.
        /// </summary>
        /// <returns>T array of entries</returns>
        [Obsolete("Please use the toArray() method instead.")]
        T[] getData();

        /// <summary>
        /// Compute the dot product.
        /// </summary>
        /// <param name="v">vector with which dot product should be computed</param>
        /// <returns>the scalar dot product of <c>this</c> and <c>v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size
        /// as <c>this</c></exception>
        T dotProduct(FieldVector<T> v);

        /// <summary>
        /// Find the orthogonal projection of this vector onto another vector.
        /// </summary>
        /// <param name="v">vector onto which <c>this</c> must be projected</param>
        /// <returns>projection of <c>this</c> onto <c>v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size
        /// as <c>this</c></exception>
        /// <exception cref="MathArithmeticException"> if <c>v</c> is the null vector.</exception>
        FieldVector<T> projection(FieldVector<T> v);

        /// <summary>
        /// Compute the outer product.
        /// </summary>
        /// <param name="v">vector with which outer product should be computed</param>
        /// <returns>the matrix outer product between instance and v</returns>
        FieldMatrix<T> outerProduct(FieldVector<T> v);

        /// <summary>
        /// Returns the entry in the specified index.
        /// </summary>
        /// <param name="index">Index location of entry to be fetched.</param>
        /// <returns>the vector entry at <c>index</c>.</returns>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        /// <remarks>
        /// See <see cref="setEntry(int, FieldElement)"/>
        /// </remarks>
        T getEntry(int index);

        /// <summary>
        /// Set a single element.
        /// </summary>
        /// <param name="index">element index.</param>
        /// <param name="value">new value for the element.</param>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        /// <remarks>
        /// See <see cref="getEntry(int)"/>
        /// </remarks>
        void setEntry(int index, T value);

        /// <summary>
        /// Returns the size of the vector.
        /// </summary>
        /// <returns>size</returns>
        int getDimension();

        /// <summary>
        /// Construct a vector by appending a vector to this vector.
        /// </summary>
        /// <param name="v">vector to append to this one.</param>
        /// <returns>a new vector</returns>
        FieldVector<T> append(FieldVector<T> v);

        /// <summary>
        /// Construct a vector by appending a T to this vector.
        /// </summary>
        /// <param name="d">T to append.</param>
        /// <returns>a new vector</returns>
        FieldVector<T> append(T d);

        /// <summary>
        /// Get a subvector from consecutive elements.
        /// </summary>
        /// <param name="index">index of first element.</param>
        /// <param name="n">number of elements to be retrieved.</param>
        /// <returns>a vector containing n elements.</returns>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        /// <exception cref="NotPositiveException"> if the number of elements if not positive.
        /// </exception>
        FieldVector<T> getSubVector(int index, int n);

        /// <summary>
        /// Set a set of consecutive elements.
        /// </summary>
        /// <param name="index">index of first element to be set.</param>
        /// <param name="v">vector containing the values to set.</param>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        void setSubVector(int index, FieldVector<T> v);

        /// <summary>
        /// Set all elements to a single value.
        /// </summary>
        /// <param name="value">single value to set for all elements</param>
        void set(T value);

        /// <summary>
        /// Convert the vector to a T array.
        /// <para>The array is independent from vector data, it's elements
        /// are copied.</para>
        /// </summary>
        /// <returns>array containing a copy of vector elements</returns>
        T[] toArray();
    }
}