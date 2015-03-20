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
using System;

namespace Math3.linear
{
    /// <summary>
    /// This class implements the <see cref="FieldVector"/> interface with a 
    /// <see cref="FieldElement"/> array.
    /// </summary>
    /// <typeparam name="T">the type of the field elements</typeparam>
#pragma warning disable 0618
    public class ArrayFieldVector<T> : FieldVector<T> where T : FieldElement<T>
    {
        /// <summary>
        /// Entries of the vector.
        /// </summary>
        private T[] data;

        /// <summary>
        /// Field to which the elements belong.
        /// </summary>
        private readonly Field<T> field;

        /// <summary>
        /// Build a 0-length vector.
        /// Zero-length vectors may be used to initialize construction of vectors
        /// by data gathering. We start with zero-length and use either the 
        /// <see cref="ArrayFieldVector(ArrayFieldVector, ArrayFieldVector)">constructor</see>
        /// or one of the <c>append</c> methods (<see cref="add(FieldVector)"/> or
        /// <see cref="append(ArrayFieldVector)"/>) to gather data into this vector.
        /// </summary>
        /// <param name="field">field to which the elements belong</param>
        public ArrayFieldVector(Field<T> field) : this(field, 0) { }

        /// <summary>
        /// Construct a vector of zeroes.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="size">Size of the vector.</param>
        public ArrayFieldVector(Field<T> field, int size)
        {
            this.field = field;
            this.data = MathArrays.buildArray(field, size);
        }

        /// <summary>
        /// Construct a vector with preset values.
        /// </summary>
        /// <param name="size">Size of the vector.</param>
        /// <param name="preset">All entries will be set with this value.</param>
        public ArrayFieldVector(int size, T preset)
            : this(preset.getField(), size)
        {
            for (int i = 0; i < data.Length; ++i)
            {
                data[i] = preset;
            }
        }

        /// <summary>
        /// Construct a vector from an array, copying the input array.
        /// This constructor needs a non-empty <c>d</c> array to retrieve
        /// the field from its first element. This implies it cannot build
        /// 0 length vectors. To build vectors from any size, one should
        /// use the <see cref="ArrayFieldVector(Field, FieldElement[])"/> constructor.
        /// </summary>
        /// <param name="d">Array.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="ZeroException"> if <c>d</c> is empty.</exception>
        /// <remarks>
        /// See <see cref="ArrayFieldVector(Field, FieldElement[])"/>
        /// </remarks>
        public ArrayFieldVector(T[] d)
        {
            MathUtils.checkNotNull(d);
            try
            {
                field = d[0].getField();
                data = (T[])d.Clone();
            }
            catch (IndexOutOfRangeException)
            {
                throw new ZeroException(new LocalizedFormats("VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT"));
            }
        }

        /// <summary>
        /// Construct a vector from an array, copying the input array.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="d">Array.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <remarks>
        /// See <see cref="ArrayFieldVector(FieldElement[])"/>
        /// </remarks>
        public ArrayFieldVector(Field<T> field, T[] d)
        {
            MathUtils.checkNotNull(d);
            this.field = field;
            data = (T[])d.Clone();
        }

        /// <summary>
        /// Create a new ArrayFieldVector using the input array as the underlying
        /// data array.
        /// If an array is built specially in order to be embedded in a
        /// ArrayFieldVector and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>. This will prevent the copying and improve
        /// performance as no new array will be built and no data will be copied.
        /// This constructor needs a non-empty <c>d</c> array to retrieve
        /// the field from its first element. This implies it cannot build
        /// 0 length vectors. To build vectors from any size, one should
        /// use the <see cref="ArrayFieldVector(Field, FieldElement[], boolean)"/>
        /// constructor.
        /// </summary>
        /// <param name="d">Data for the new vector.</param>
        /// <param name="copyArray">If <c>true</c>, the input array will be copied,
        /// otherwise it will be referenced.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="ZeroException"> if <c>d</c> is empty.</exception>
        /// <remarks>
        /// See <see cref="ArrayFieldVector(FieldElement[])"/><para/>
        /// See <see cref="ArrayFieldVector(Field, FieldElement[], boolean)"/>
        /// </remarks>
        public ArrayFieldVector(T[] d, Boolean copyArray)
        {
            MathUtils.checkNotNull(d);
            if (d.Length == 0)
            {
                throw new ZeroException(new LocalizedFormats("VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT"));
            }
            field = d[0].getField();
            data = copyArray ? (T[])d.Clone() : d;
        }

        /// <summary>
        /// Create a new ArrayFieldVector using the input array as the underlying
        /// data array.
        /// If an array is built specially in order to be embedded in a
        /// ArrayFieldVector and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>. This will prevent the copying and improve
        /// performance as no new array will be built and no data will be copied.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="d">Data for the new vector.</param>
        /// <param name="copyArray">If <c>true</c>, the input array will be copied,
        /// otherwise it will be referenced.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <remarks>
        /// See <see cref="ArrayFieldVector(FieldElement[], boolean)"/>
        /// </remarks>
        public ArrayFieldVector(Field<T> field, T[] d, Boolean copyArray)
        {
            MathUtils.checkNotNull(d);
            this.field = field;
            data = copyArray ? (T[])d.Clone() : d;
        }

        /// <summary>
        /// Construct a vector from part of a array.
        /// </summary>
        /// <param name="d">Array.</param>
        /// <param name="pos">Position of the first entry.</param>
        /// <param name="size">Number of entries to copy.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if the size of <c>d</c> is less
        /// than <c>pos + size</c>.</exception>
        public ArrayFieldVector(T[] d, int pos, int size)
        {
            MathUtils.checkNotNull(d);
            if (d.Length < pos + size)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(pos + size, d.Length, true);
            }
            field = d[0].getField();
            data = MathArrays.buildArray(field, size);
            Array.Copy(d, pos, data, 0, size);
        }

        /// <summary>
        /// Construct a vector from part of a array.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="d">Array.</param>
        /// <param name="pos">Position of the first entry.</param>
        /// <param name="size">Number of entries to copy.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if the size of <c>d</c> is less
        /// than <c>pos + size</c>.</exception>
        public ArrayFieldVector(Field<T> field, T[] d, int pos, int size)
        {
            MathUtils.checkNotNull(d);
            if (d.Length < pos + size)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(pos + size, d.Length, true);
            }
            this.field = field;
            data = MathArrays.buildArray(field, size);
            Array.Copy(d, pos, data, 0, size);
        }

        /// <summary>
        /// Construct a vector from another vector, using a deep copy.
        /// </summary>
        /// <param name="v">Vector to copy.</param>
        /// <exception cref="NullArgumentException"> if <c>v</c> is <c>null</c>.</exception>
        public ArrayFieldVector(FieldVector<T> v)
        {
            MathUtils.checkNotNull(v);
            field = v.getField();
            data = MathArrays.buildArray(field, v.getDimension());
            for (int i = 0; i < data.Length; ++i)
            {
                data[i] = v.getEntry(i);
            }
        }

        /// <summary>
        /// Construct a vector from another vector, using a deep copy.
        /// </summary>
        /// <param name="v">Vector to copy.</param>
        /// <exception cref="NullArgumentException"> if <c>v</c> is <c>null</c>.</exception>
        public ArrayFieldVector(ArrayFieldVector<T> v)
        {
            MathUtils.checkNotNull(v);
            field = v.getField();
            data = (T[])v.data.Clone();
        }

        /// <summary>
        /// Construct a vector from another vector.
        /// </summary>
        /// <param name="v">Vector to copy.</param>
        /// <param name="deep">If <c>true</c> perform a deep copy, otherwise perform
        /// a shallow copy</param>
        /// <exception cref="NullArgumentException"> if <c>v</c> is <c>null</c>.</exception>
        public ArrayFieldVector(ArrayFieldVector<T> v, Boolean deep)
        {
            MathUtils.checkNotNull(v);
            field = v.getField();
            data = deep ? (T[])v.data.Clone() : v.data;
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        [Obsolete("replaced by ArrayFieldVector(FieldVector, FieldVector)")]
        public ArrayFieldVector(ArrayFieldVector<T> v1, ArrayFieldVector<T> v2) : this((FieldVector<T>)v1, (FieldVector<T>)v2) { }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        public ArrayFieldVector(FieldVector<T> v1, FieldVector<T> v2)
        {
            MathUtils.checkNotNull(v1);
            MathUtils.checkNotNull(v2);
            field = v1.getField();
            T[] v1Data =
                    (v1 is ArrayFieldVector<T>) ? ((ArrayFieldVector<T>)v1).data : v1.toArray();
            T[] v2Data =
                    (v2 is ArrayFieldVector<T>) ? ((ArrayFieldVector<T>)v2).data : v2.toArray();
            data = MathArrays.buildArray(field, v1Data.Length + v2Data.Length);
            Array.Copy(v1Data, 0, data, 0, v1Data.Length);
            Array.Copy(v2Data, 0, data, v1Data.Length, v2Data.Length);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        [Obsolete("replaced by ArrayFieldVector(FieldVector, FieldElement[])")]
        public ArrayFieldVector(ArrayFieldVector<T> v1, T[] v2) : this((FieldVector<T>)v1, v2) { }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        public ArrayFieldVector(FieldVector<T> v1, T[] v2)
        {
            MathUtils.checkNotNull(v1);
            MathUtils.checkNotNull(v2);
            field = v1.getField();
            T[] v1Data = (v1 is ArrayFieldVector<T>) ? ((ArrayFieldVector<T>)v1).data : v1.toArray();
            data = MathArrays.buildArray(field, v1Data.Length + v2.Length);
            Array.Copy(v1Data, 0, data, 0, v1Data.Length);
            Array.Copy(v2, 0, data, v1Data.Length, v2.Length);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        [Obsolete("replaced by ArrayFieldVector(FieldElement[], FieldVector)")]
        public ArrayFieldVector(T[] v1, ArrayFieldVector<T> v2) : this(v1, (FieldVector<T>)v2) { }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        public ArrayFieldVector(T[] v1, FieldVector<T> v2)
        {
            MathUtils.checkNotNull(v1);
            MathUtils.checkNotNull(v2);
            field = v2.getField();
            T[] v2Data = (v2 is ArrayFieldVector<T>) ? ((ArrayFieldVector<T>)v2).data : v2.toArray();
            data = MathArrays.buildArray(field, v1.Length + v2Data.Length);
            Array.Copy(v1, 0, data, 0, v1.Length);
            Array.Copy(v2Data, 0, data, v1.Length, v2Data.Length);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// This constructor needs at least one non-empty array to retrieve
        /// the field from its first element. This implies it cannot build
        /// 0 length vectors. To build vectors from any size, one should
        /// use the <see cref="ArrayFieldVector(Field, FieldElement[], FieldElement[])"/>
        /// constructor.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        /// <exception cref="ZeroException"> if both arrays are empty.</exception>
        /// <remarks>
        /// See <see cref="ArrayFieldVector(Field, FieldElement[], FieldElement[])"/>
        /// </remarks>
        public ArrayFieldVector(T[] v1, T[] v2)
        {
            MathUtils.checkNotNull(v1);
            MathUtils.checkNotNull(v2);
            if (v1.Length + v2.Length == 0)
            {
                throw new ZeroException(new LocalizedFormats("VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT"));
            }
            data = MathArrays.buildArray(v1[0].getField(), v1.Length + v2.Length);
            Array.Copy(v1, 0, data, 0, v1.Length);
            Array.Copy(v2, 0, data, v1.Length, v2.Length);
            field = data[0].getField();
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="field">Field to which the elements belong.</param>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        /// <exception cref="NullArgumentException"> if <c>v1</c> or <c>v2</c> is
        /// <c>null</c>.</exception>
        /// <exception cref="ZeroException"> if both arrays are empty.</exception>
        /// <remarks>
        /// See <see cref="ArrayFieldVector(FieldElement[], FieldElement[])"/>
        /// </remarks>
        public ArrayFieldVector(Field<T> field, T[] v1, T[] v2)
        {
            MathUtils.checkNotNull(v1);
            MathUtils.checkNotNull(v2);
            if (v1.Length + v2.Length == 0)
            {
                throw new ZeroException(new LocalizedFormats("VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT"));
            }
            data = MathArrays.buildArray(field, v1.Length + v2.Length);
            Array.Copy(v1, 0, data, 0, v1.Length);
            Array.Copy(v2, 0, data, v1.Length, v2.Length);
            this.field = field;
        }

        /// <inheritdoc/>
        public Field<T> getField()
        {
            return field;
        }

        /// <inheritdoc/>
        public FieldVector<T> copy()
        {
            return new ArrayFieldVector<T>(this, true);
        }

        /// <inheritdoc/>
        public FieldVector<T> add(FieldVector<T> v)
        {
            try
            {
                return add((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                checkVectorDimensions(v);
                T[] outp = MathArrays.buildArray(field, data.Length);
                for (int i = 0; i < data.Length; i++)
                {
                    outp[i] = data[i].add(v.getEntry(i));
                }
                return new ArrayFieldVector<T>(field, outp, false);
            }
        }

        /// <summary>
        /// Compute the sum of <c>this</c> and <c>v</c>.
        /// </summary>
        /// <param name="v">vector to be added</param>
        /// <returns><c>this + v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c></exception>
        public ArrayFieldVector<T> add(ArrayFieldVector<T> v)
        {
            checkVectorDimensions(v.data.Length);
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].add(v.data[i]);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> subtract(FieldVector<T> v)
        {
            try
            {
                return subtract((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                checkVectorDimensions(v);
                T[] outp = MathArrays.buildArray(field, data.Length);
                for (int i = 0; i < data.Length; i++)
                {
                    outp[i] = data[i].subtract(v.getEntry(i));
                }
                return new ArrayFieldVector<T>(field, outp, false);
            }
        }

        /// <summary>
        /// Compute <c>this</c> minus <c>v</c>.
        /// </summary>
        /// <param name="v">vector to be subtracted</param>
        /// <returns><c>this - v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c></exception>
        public ArrayFieldVector<T> subtract(ArrayFieldVector<T> v)
        {
            checkVectorDimensions(v.data.Length);
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].subtract(v.data[i]);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> mapAdd(T d)
        {
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].add(d);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> mapAddToSelf(T d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = data[i].add(d);
            }
            return this;
        }

        /// <inheritdoc/>
        public FieldVector<T> mapSubtract(T d)
        {
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].subtract(d);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> mapSubtractToSelf(T d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = data[i].subtract(d);
            }
            return this;
        }

        /// <inheritdoc/>
        public FieldVector<T> mapMultiply(T d)
        {
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].multiply(d);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> mapMultiplyToSelf(T d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = data[i].multiply(d);
            }
            return this;
        }

        /// <inheritdoc/>
        public FieldVector<T> mapDivide(T d)
        {
            MathUtils.checkNotNull(d);
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].divide(d);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> mapDivideToSelf(T d)
        {
            MathUtils.checkNotNull(d);
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = data[i].divide(d);
            }
            return this;
        }

        /// <inheritdoc/>
        public FieldVector<T> mapInv()
        {
            T[] outp = MathArrays.buildArray(field, data.Length);
            T one = field.getOne();
            for (int i = 0; i < data.Length; i++)
            {
                try
                {
                    outp[i] = one.divide(data[i]);
                }
                catch (MathArithmeticException)
                {
                    throw new MathArithmeticException(new LocalizedFormats("INDEX"), i);
                }
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> mapInvToSelf()
        {
            T one = field.getOne();
            for (int i = 0; i < data.Length; i++)
            {
                try
                {
                    data[i] = one.divide(data[i]);
                }
                catch (MathArithmeticException)
                {
                    throw new MathArithmeticException(new LocalizedFormats("INDEX"), i);
                }
            }
            return this;
        }

        /// <inheritdoc/>
        public FieldVector<T> ebeMultiply(FieldVector<T> v)
        {
            try
            {
                return ebeMultiply((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                checkVectorDimensions(v);
                T[] outp = MathArrays.buildArray(field, data.Length);
                for (int i = 0; i < data.Length; i++)
                {
                    outp[i] = data[i].multiply(v.getEntry(i));
                }
                return new ArrayFieldVector<T>(field, outp, false);
            }
        }

        /// <summary>
        /// Element-by-element multiplication.
        /// </summary>
        /// <param name="v">vector by which instance elements must be multiplied</param>
        /// <returns>vector containing <c>this[i] * v[i]</c> for all <c>i</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c></exception>
        public ArrayFieldVector<T> ebeMultiply(ArrayFieldVector<T> v)
        {
            checkVectorDimensions(v.data.Length);
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                outp[i] = data[i].multiply(v.data[i]);
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> ebeDivide(FieldVector<T> v)
        {
            try
            {
                return ebeDivide((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                checkVectorDimensions(v);
                T[] outp = MathArrays.buildArray(field, data.Length);
                for (int i = 0; i < data.Length; i++)
                {
                    try
                    {
                        outp[i] = data[i].divide(v.getEntry(i));
                    }
                    catch (MathArithmeticException)
                    {
                        throw new MathArithmeticException(new LocalizedFormats("INDEX"), i);
                    }
                }
                return new ArrayFieldVector<T>(field, outp, false);
            }
        }

        /// <summary>
        /// Element-by-element division.
        /// </summary>
        /// <param name="v">vector by which instance elements must be divided</param>
        /// <returns>a vector containing <c>this[i] / v[i]</c> for all <c>i</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c></exception>
        /// <exception cref="MathArithmeticException"> if one entry of <c>v</c> is zero.
        /// </exception>
        public ArrayFieldVector<T> ebeDivide(ArrayFieldVector<T> v)
        {
            checkVectorDimensions(v.data.Length);
            T[] outp = MathArrays.buildArray(field, data.Length);
            for (int i = 0; i < data.Length; i++)
            {
                try
                {
                    outp[i] = data[i].divide(v.data[i]);
                }
                catch (MathArithmeticException)
                {
                    throw new MathArithmeticException(new LocalizedFormats("INDEX"), i);
                }
            }
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public T[] getData()
        {
            return (T[])data.Clone();
        }

        /// <summary>
        /// Returns a reference to the underlying data array.
        /// <para>Does not make a fresh copy of the underlying data.</para>
        /// </summary>
        /// <returns>array of entries</returns>
        public T[] getDataRef()
        {
            return data;
        }

        /// <inheritdoc/>
        public T dotProduct(FieldVector<T> v)
        {
            try
            {
                return dotProduct((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                checkVectorDimensions(v);
                T dot = field.getZero();
                for (int i = 0; i < data.Length; i++)
                {
                    dot = dot.add(data[i].multiply(v.getEntry(i)));
                }
                return dot;
            }
        }

        /// <summary>
        /// Compute the dot product.
        /// </summary>
        /// <param name="v">vector with which dot product should be computed</param>
        /// <returns>the scalar dot product of <c>this</c> and <c>v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// mo<c>this</c></exception>
        public T dotProduct(ArrayFieldVector<T> v)
        {
            checkVectorDimensions(v.data.Length);
            T dot = field.getZero();
            for (int i = 0; i < data.Length; i++)
            {
                dot = dot.add(data[i].multiply(v.data[i]));
            }
            return dot;
        }

        /// <inheritdoc/>
        public FieldVector<T> projection(FieldVector<T> v)
        {
            return v.mapMultiply(dotProduct(v).divide(v.dotProduct(v)));
        }

        /// <summary>
        /// Find the orthogonal projection of this vector onto another vector.
        /// </summary>
        /// <param name="v">vector onto which <c>this</c> must be projected</param>
        /// <returns>projection of <c>this</c> onto <c>v</c></returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c></exception>
        /// <exception cref="MathArithmeticException"> if <c>v</c> is the null vector.</exception>
        public ArrayFieldVector<T> projection(ArrayFieldVector<T> v)
        {
            return (ArrayFieldVector<T>)v.mapMultiply(dotProduct(v).divide(v.dotProduct(v)));
        }

        /// <inheritdoc/>
        public FieldMatrix<T> outerProduct(FieldVector<T> v)
        {
            try
            {
                return outerProduct((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                int m = data.Length;
                int n = v.getDimension();
                FieldMatrix<T> outp = new Array2DRowFieldMatrix<T>(field, m, n);
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        outp.setEntry(i, j, data[i].multiply(v.getEntry(j)));
                    }
                }
                return outp;
            }
        }

        /// <summary>
        /// Compute the outer product.
        /// </summary>
        /// <param name="v">vector with which outer product should be computed</param>
        /// <returns>the matrix outer product between instance and v</returns>
        public FieldMatrix<T> outerProduct(ArrayFieldVector<T> v)
        {
            int m = data.Length;
            int n = v.data.Length;
            FieldMatrix<T> outp = new Array2DRowFieldMatrix<T>(field, m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    outp.setEntry(i, j, data[i].multiply(v.data[j]));
                }
            }
            return outp;
        }

        /// <inheritdoc/>
        public T getEntry(int index)
        {
            return data[index];
        }

        /// <inheritdoc/>
        public int getDimension()
        {
            return data.Length;
        }

        /// <inheritdoc/>
        public FieldVector<T> append(FieldVector<T> v)
        {
            try
            {
                return append((ArrayFieldVector<T>)v);
            }
            catch (InvalidCastException)
            {
                return new ArrayFieldVector<T>(this, new ArrayFieldVector<T>(v));
            }
        }

        /// <summary>
        /// Construct a vector by appending a vector to this vector.
        /// </summary>
        /// <param name="v">vector to append to this one.</param>
        /// <returns>a new vector</returns>
        public ArrayFieldVector<T> append(ArrayFieldVector<T> v)
        {
            return new ArrayFieldVector<T>(this, v);
        }

        /// <inheritdoc/>
        public FieldVector<T> append(T inp)
        {
            T[] outp = MathArrays.buildArray(field, data.Length + 1);
            Array.Copy(data, 0, outp, 0, data.Length);
            outp[data.Length] = inp;
            return new ArrayFieldVector<T>(field, outp, false);
        }

        /// <inheritdoc/>
        public FieldVector<T> getSubVector(int index, int n)
        {
            if (n < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE"), n);
            }
            ArrayFieldVector<T> outp = new ArrayFieldVector<T>(field, n);
            try
            {
                Array.Copy(data, index, outp.data, 0, n);
            }
            catch (IndexOutOfRangeException)
            {
                checkIndex(index);
                checkIndex(index + n - 1);
            }
            return outp;
        }

        /// <inheritdoc/>
        public void setEntry(int index, T value)
        {
            try
            {
                data[index] = value;
            }
            catch (IndexOutOfRangeException)
            {
                checkIndex(index);
            }
        }

        /// <inheritdoc/>
        public void setSubVector(int index, FieldVector<T> v)
        {
            try
            {
                try
                {
                    set(index, (ArrayFieldVector<T>)v);
                }
                catch (InvalidCastException)
                {
                    for (int i = index; i < index + v.getDimension(); ++i)
                    {
                        data[i] = v.getEntry(i - index);
                    }
                }
            }
            catch (IndexOutOfRangeException)
            {
                checkIndex(index);
                checkIndex(index + v.getDimension() - 1);
            }
        }

        /// <summary>
        /// Set a set of consecutive elements.
        /// </summary>
        /// <param name="index">index of first element to be set.</param>
        /// <param name="v">vector containing the values to set.</param>
        /// <exception cref="OutOfRangeException"> if the index is invalid.</exception>
        public void set(int index, ArrayFieldVector<T> v)
        {
            try
            {
                Array.Copy(v.data, 0, data, index, v.data.Length);
            }
            catch (IndexOutOfRangeException)
            {
                checkIndex(index);
                checkIndex(index + v.data.Length - 1);
            }
        }

        /// <inheritdoc/>
        public void set(T value)
        {
            for (int i = 0; i < data.Length; ++i)
            {
                data[i] = value;
            }
        }

        /// <inheritdoc/>
        public T[] toArray()
        {
            return (T[])data.Clone();
        }

        /// <summary>
        /// Check if instance and specified vectors have the same dimension.
        /// </summary>
        /// <param name="v">vector to compare instance with</param>
        /// <exception cref="DimensionMismatchException"> if the vectors do not
        /// have the same dimensions</exception>
        protected void checkVectorDimensions(FieldVector<T> v)
        {
            checkVectorDimensions(v.getDimension());
        }

        /// <summary>
        /// Check if instance dimension is equal to some expected value.
        /// </summary>
        /// <param name="n">Expected dimension.</param>
        /// <exception cref="DimensionMismatchException"> if the dimension is not equal to the
        /// size of <c>this</c> vector.</exception>
        protected void checkVectorDimensions(int n)
        {
            if (data.Length != n)
            {
                throw new DimensionMismatchException(data.Length, n);
            }
        }

        /// <summary>
        /// Visits (but does not alter) all entries of this vector in default order
        /// (increasing index).
        /// </summary>
        /// <param name="visitor">the visitor to be used to process the entries of this
        /// vector</param>
        /// <returns>the value returned by <see cref="FieldVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        public T walkInDefaultOrder(FieldVectorPreservingVisitor<T> visitor)
        {
            int dim = getDimension();
            visitor.start(dim, 0, dim - 1);
            for (int i = 0; i < dim; i++)
            {
                visitor.visit(i, getEntry(i));
            }
            return visitor.end();
        }

        /// <summary>
        /// Visits (but does not alter) some entries of this vector in default order
        /// (increasing index).
        /// </summary>
        /// <param name="visitor">visitor to be used to process the entries of this vector</param>
        /// <param name="start">the index of the first entry to be visited</param>
        /// <param name="end">the index of the last entry to be visited (inclusive)</param>
        /// <returns>the value returned by <see cref="FieldVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public T walkInDefaultOrder(FieldVectorPreservingVisitor<T> visitor, int start, int end)
        {
            checkIndices(start, end);
            visitor.start(getDimension(), start, end);
            for (int i = start; i <= end; i++)
            {
                visitor.visit(i, getEntry(i));
            }
            return visitor.end();
        }

        /// <summary>
        /// Visits (but does not alter) all entries of this vector in optimized
        /// order. The order in which the entries are visited is selected so as to
        /// lead to the most efficient implementation; it might depend on the
        /// concrete implementation of this abstract class.
        /// </summary>
        /// <param name="visitor">the visitor to be used to process the entries of this
        /// vector</param>
        /// <returns>the value returned by <see cref="FieldVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        public T walkInOptimizedOrder(FieldVectorPreservingVisitor<T> visitor)
        {
            return walkInDefaultOrder(visitor);
        }

        /// <summary>
        /// Visits (but does not alter) some entries of this vector in optimized
        /// order. The order in which the entries are visited is selected so as to
        /// lead to the most efficient implementation; it might depend on the
        /// concrete implementation of this abstract class.
        /// </summary>
        /// <param name="visitor">visitor to be used to process the entries of this vector</param>
        /// <param name="start">the index of the first entry to be visited</param>
        /// <param name="end">the index of the last entry to be visited (inclusive)</param>
        /// <returns>the value returned by <see cref="FieldVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public T walkInOptimizedOrder(FieldVectorPreservingVisitor<T> visitor, int start, int end)
        {
            return walkInDefaultOrder(visitor, start, end);
        }

        /// <summary>
        /// Visits (and possibly alters) all entries of this vector in default order
        /// (increasing index).
        /// </summary>
        /// <param name="visitor">the visitor to be used to process and modify the entries
        /// of this vector</param>
        /// <returns>the value returned by <see cref="FieldVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        public T walkInDefaultOrder(FieldVectorChangingVisitor<T> visitor)
        {
            int dim = getDimension();
            visitor.start(dim, 0, dim - 1);
            for (int i = 0; i < dim; i++)
            {
                setEntry(i, visitor.visit(i, getEntry(i)));
            }
            return visitor.end();
        }

        /// <summary>
        /// Visits (and possibly alters) some entries of this vector in default order
        /// (increasing index).
        /// </summary>
        /// <param name="visitor">visitor to be used to process the entries of this vector</param>
        /// <param name="start">the index of the first entry to be visited</param>
        /// <param name="end">the index of the last entry to be visited (inclusive)</param>
        /// <returns>the value returned by <see cref="FieldVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public T walkInDefaultOrder(FieldVectorChangingVisitor<T> visitor, int start, int end)
        {
            checkIndices(start, end);
            visitor.start(getDimension(), start, end);
            for (int i = start; i <= end; i++)
            {
                setEntry(i, visitor.visit(i, getEntry(i)));
            }
            return visitor.end();
        }

        /// <summary>
        /// Visits (and possibly alters) all entries of this vector in optimized
        /// order. The order in which the entries are visited is selected so as to
        /// lead to the most efficient implementation; it might depend on the
        /// concrete implementation of this abstract class.
        /// </summary>
        /// <param name="visitor">the visitor to be used to process the entries of this
        /// vector</param>
        /// <returns>the value returned by <see cref="FieldVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        public T walkInOptimizedOrder(FieldVectorChangingVisitor<T> visitor)
        {
            return walkInDefaultOrder(visitor);
        }

        /// <summary>
        /// Visits (and possibly change) some entries of this vector in optimized
        /// order. The order in which the entries are visited is selected so as to
        /// lead to the most efficient implementation; it might depend on the
        /// concrete implementation of this abstract class.
        /// </summary>
        /// <param name="visitor">visitor to be used to process the entries of this vector</param>
        /// <param name="start">the index of the first entry to be visited</param>
        /// <param name="end">the index of the last entry to be visited (inclusive)</param>
        /// <returns>the value returned by <see cref="FieldVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public T walkInOptimizedOrder(FieldVectorChangingVisitor<T> visitor, int start, int end)
        {
            return walkInDefaultOrder(visitor, start, end);
        }

        /// <summary>
        /// Test for the equality of two vectors.
        /// </summary>
        /// <param name="other">Object to test for equality.</param>
        /// <returns><c>true</c> if two vector objects are equal, <c>false</c>
        /// otherwise.</returns>
        public override Boolean Equals(Object other)
        {
            if (this == other)
            {
                return true;
            }
            if (other == null)
            {
                return false;
            }

            try
            {
                FieldVector<T> rhs = (FieldVector<T>)other;
                if (data.Length != rhs.getDimension())
                {
                    return false;
                }

                for (int i = 0; i < data.Length; ++i)
                {
                    if (!data[i].Equals(rhs.getEntry(i)))
                    {
                        return false;
                    }
                }
                return true;
            }
            catch (InvalidCastException)
            {
                // ignore exception
                return false;
            }
        }

        /// <summary>
        /// Get a hashCode for the real vector.
        /// <para>All NaN values have the same hash code.</para>
        /// </summary>
        /// <returns>a hash code value for this object</returns>
        public override int GetHashCode()
        {
            int h = 3542;
            foreach (T a in data)
            {
                h ^= a.GetHashCode();
            }
            return h;
        }

        /// <summary>
        /// Check if an index is valid.
        /// </summary>
        /// <param name="index">Index to check.</param>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        private void checkIndex(int index)
        {
            if (index < 0 || index >= getDimension())
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("INDEX"), index, 0, getDimension() - 1);
            }
        }

        /// <summary>
        /// Checks that the indices of a subvector are valid.
        /// </summary>
        /// <param name="start">the index of the first entry of the subvector</param>
        /// <param name="end">the index of the last entry of the subvector (inclusive)</param>
        /// <exception cref="OutOfRangeException"> if <c>start</c> of <c>end</c> are not valid
        /// </exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c></exception>
        private void checkIndices(int start, int end)
        {
            int dim = getDimension();
            if ((start < 0) || (start >= dim))
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("INDEX"), start, 0, dim - 1);
            }
            if ((end < 0) || (end >= dim))
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("INDEX"), end, 0, dim - 1);
            }
            if (end < start)
            {
                throw new NumberIsTooSmallException<Int32, Int32>(new LocalizedFormats("INITIAL_ROW_AFTER_FINAL_ROW"), end, start, false);
            }
        }
    }
}
