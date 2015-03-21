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
using Math3.analysis;
using Math3.analysis.function;
using Math3.exception;
using Math3.exception.util;
using Math3.util;
using System;
using System.Collections.Generic;

namespace Math3.linear
{
    /// <summary>
    /// Class defining a real-valued vector with basic algebraic operations.
    /// <para>
    /// vector element indexing is 0-based -- e.g., <c>getEntry(0)</c>
    /// returns the first element of the vector.
    /// </para>
    /// <para>
    /// The <c>code map</c> and <c>mapToSelf</c> methods operate
    /// on vectors element-wise, i.e. they perform the same operation (adding a scalar,
    /// applying a function ...) on each element in turn. The <c>map</c>
    /// versions create a new vector to hold the result and do not change the instance.
    /// The <c>mapToSelf</c> version uses the instance itself to store the
    /// results, so the instance is changed by this method. In all cases, the result
    /// vector is returned by the methods, allowing the fluent API
    /// style, like this:
    /// </para>
    /// <code>
    ///   RealVector result = v.mapAddToSelf(3.4).mapToSelf(new Tan()).mapToSelf(new Power(2.3));
    /// </code>
    /// </summary>
    public abstract class RealVector
    {
        /// <summary>
        /// Returns the size of the vector. 
        /// </summary>
        /// <returns>the size of this vector.</returns>
        public abstract int getDimension();

        /// <summary>
        /// Return the entry at the specified index.
        /// </summary>
        /// <param name="index">Index location of entry to be fetched.</param>
        /// <returns>the vector entry at <c>index</c>.</returns>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        /// <remarks>
        /// See <see cref="setEntry(int, double)"/>
        /// </remarks>
        public abstract double getEntry(int index);

        /// <summary>
        /// Set a single element.
        /// </summary>
        /// <param name="index">element index.</param>
        /// <param name="value">new value for the element.</param>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        /// <remarks>
        /// See <see cref="setEntry(int, double)"/>
        /// </remarks>
        public abstract void setEntry(int index, double value);

        /// <summary>
        /// Change an entry at the specified index.
        /// </summary>
        /// <param name="index">Index location of entry to be set.</param>
        /// <param name="increment">Value to add to the vector entry.</param>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        public void addToEntry(int index, double increment)
        {
            setEntry(index, getEntry(index) + increment);
        }

        /// <summary>
        /// Construct a new vector by appending a vector to this vector.
        /// </summary>
        /// <param name="v">vector to append to this one.</param>
        /// <returns>a new vector.</returns>
        public abstract RealVector append(RealVector v);

        /// <summary>
        /// Construct a new vector by appending a double to this vector.
        /// </summary>
        /// <param name="d">double to append.</param>
        /// <returns>a new vector.</returns>
        public abstract RealVector append(double d);

        /// <summary>
        /// Get a subvector from consecutive elements.
        /// </summary>
        /// <param name="index">index of first element.</param>
        /// <param name="n">number of elements to be retrieved.</param>
        /// <returns>a vector containing n elements.</returns>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        /// <exception cref="NotPositiveException"> if the number of elements is not positive
        /// </exception>
        public abstract RealVector getSubVector(int index, int n);

        /// <summary>
        /// Set a sequence of consecutive elements.
        /// </summary>
        /// <param name="index">index of first element to be set.</param>
        /// <param name="v">vector containing the values to set.</param>
        /// <exception cref="OutOfRangeException"> if the index is not valid.</exception>
        public abstract void setSubVector(int index, RealVector v);

        /// <summary>
        /// Check whether any coordinate of this vector is <c>NaN</c>.
        /// </summary>
        /// <returns><c>true</c> if any coordinate of this vector is <c>NaN</c>,
        /// <c>false</c> otherwise.</returns>
        public abstract Boolean isNaN();

        /// <summary>
        /// Check whether any coordinate of this vector is infinite and none are <c>NaN</c>.
        /// </summary>
        /// <returns><c>true</c> if any coordinate of this vector is infinite and
        /// none are <c>NaN</c>, <c>false</c> otherwise.</returns>
        public abstract Boolean isInfinite();

        /// <summary>
        /// Check if instance and specified vectors have the same dimension.
        /// </summary>
        /// <param name="v">Vector to compare instance with.</param>
        /// <exception cref="DimensionMismatchException"> if the vectors do not
        /// have the same dimension.</exception>
        protected void checkVectorDimensions(RealVector v)
        {
            checkVectorDimensions(v.getDimension());
        }

        /// <summary>
        /// Check if instance dimension is equal to some expected value.
        /// </summary>
        /// <param name="n">Expected dimension.</param>
        /// <exception cref="DimensionMismatchException"> if the dimension is
        /// inconsistent with the vector size.</exception>
        protected void checkVectorDimensions(int n)
        {
            int d = getDimension();
            if (d != n)
            {
                throw new DimensionMismatchException(d, n);
            }
        }

        /// <summary>
        /// Check if an index is valid.
        /// </summary>
        /// <param name="index">Index to check.</param>
        /// <exception cref="OutOfRangeException"> if <c>index</c> is not valid.</exception>
        protected void checkIndex(int index)
        {
            if (index < 0 ||
                index >= getDimension())
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
        protected void checkIndices(int start, int end)
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
                // TODO Use more specific error message
                throw new NumberIsTooSmallException<Int32, Int32>(new LocalizedFormats("INITIAL_ROW_AFTER_FINAL_ROW"), end, start, false);
            }
        }

        /// <summary>
        /// Compute the sum of this vector and <c>v</c>.
        /// Returns a new vector. Does not change instance data.
        /// </summary>
        /// <param name="v">Vector to be added.</param>
        /// <returns><c>this</c> + <c>v</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public RealVector add(RealVector v)
        {
            checkVectorDimensions(v);
            RealVector result = v.copy();
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                int index = e.getIndex();
                result.setEntry(index, e.getValue() + result.getEntry(index));
            }
            return result;
        }

        /// <summary>
        /// Subtract <c>v</c> from this vector.
        /// Returns a new vector. Does not change instance data.
        /// </summary>
        /// <param name="v">Vector to be subtracted.</param>
        /// <returns><c>this</c> - <c>v</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public RealVector subtract(RealVector v)
        {
            checkVectorDimensions(v);
            RealVector result = v.mapMultiply(-1d);
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                int index = e.getIndex();
                result.setEntry(index, e.getValue() + result.getEntry(index));
            }
            return result;
        }

        /// <summary>
        /// Add a value to each entry.
        /// Returns a new vector. Does not change instance data.
        /// </summary>
        /// <param name="d">Value to be added to each entry.</param>
        /// <returns><c>this</c> + <c>d</c>.</returns>
        public RealVector mapAdd(double d)
        {
            return copy().mapAddToSelf(d);
        }

        /// <summary>
        /// Add a value to each entry.
        /// The instance is changed in-place.
        /// </summary>
        /// <param name="d">Value to be added to each entry.</param>
        /// <returns><c>this</c>.</returns>
        public RealVector mapAddToSelf(double d)
        {
            if (d != 0)
            {
                return mapToSelf(FunctionUtils.fix2ndArgument(new Add(), d));
            }
            return this;
        }

        /// <summary>
        /// Returns a (deep) copy of this vector.
        /// </summary>
        /// <returns>a vector copy.</returns>
        public abstract RealVector copy();

        /// <summary>
        /// Compute the dot product of this vector with <c>v</c>.
        /// </summary>
        /// <param name="v">Vector with which dot product should be computed</param>
        /// <returns>the scalar dot product between this instance and <c>v</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public double dotProduct(RealVector v)
        {
            checkVectorDimensions(v);
            double d = 0;
            int n = getDimension();
            for (int i = 0; i < n; i++)
            {
                d += getEntry(i) * v.getEntry(i);
            }
            return d;
        }

        /// <summary>
        /// Computes the cosine of the angle between this vector and the
        /// argument.
        /// </summary>
        /// <param name="v">Vector.</param>
        /// <returns>the cosine of the angle between this vector and <c>v</c>.</returns>
        /// <exception cref="MathArithmeticException"> if <c>this</c> or <c>v</c> is the null
        /// vector</exception>
        /// <exception cref="DimensionMismatchException"> if the dimensions of <c>this</c> and
        /// <c>v</c> do not match</exception>
        public double cosine(RealVector v)
        {
            double norm = getNorm();
            double vNorm = v.getNorm();

            if (norm == 0 ||
                vNorm == 0)
            {
                throw new MathArithmeticException(new LocalizedFormats("ZERO_NORM"));
            }
            return dotProduct(v) / (norm * vNorm);
        }

        /// <summary>
        /// Element-by-element division.
        /// </summary>
        /// <param name="v">Vector by which instance elements must be divided.</param>
        /// <returns>a vector containing this[i] / v[i] for all i.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public abstract RealVector ebeDivide(RealVector v);

        /// <summary>
        /// Element-by-element multiplication.
        /// </summary>
        /// <param name="v">Vector by which instance elements must be multiplied</param>
        /// <returns>a vector containing this[i] * v[i] for all i.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public abstract RealVector ebeMultiply(RealVector v);

        /// <summary>
        /// Distance between two vectors.
        /// <para>This method computes the distance consistent with the
        /// L^2 norm, i.e. the square root of the sum of
        /// element differences, or Euclidean distance.</para>
        /// </summary>
        /// <param name="v">Vector to which distance is requested.</param>
        /// <returns>the distance between two vectors.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        /// <remarks>
        /// See <see cref="getL1Distance(RealVector)"/><para/>
        /// See <see cref="getLInfDistance(RealVector)"/><para/>
        /// See <see cref="getNorm()"/>
        /// </remarks>
        public double getDistance(RealVector v)
        {
            checkVectorDimensions(v);
            double d = 0;
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                double diff = e.getValue() - v.getEntry(e.getIndex());
                d += diff * diff;
            }
            return FastMath.sqrt(d);
        }

        /// <summary>
        /// Returns the L^2 norm of the vector.
        /// <para>The L^2 norm is the root of the sum of
        /// the squared elements.</para>
        /// </summary>
        /// <returns>the norm.</returns>
        /// <remarks>
        /// See <see cref="getL1Norm()"/><para/>
        /// See <see cref="getLInfNorm()"/><para/>
        /// See <see cref="getDistance(RealVector)"/>
        /// </remarks>
        public double getNorm()
        {
            double sum = 0;
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                double value = e.getValue();
                sum += value * value;
            }
            return FastMath.sqrt(sum);
        }

        /// <summary>
        /// Returns the L<sub>1</sub> norm of the vector.
        /// <para>The L<sub>1</sub> norm is the sum of the absolute
        /// values of the elements.</para>
        /// </summary>
        /// <returns>the norm.</returns>
        /// <remarks>
        /// See<see cref="getNorm()"/><para/>
        /// See <see cref="getLInfNorm()"/>para/>
        /// See <see cref="getL1Distance(RealVector)"/>
        /// </remarks>
        public double getL1Norm()
        {
            double norm = 0;
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                norm += FastMath.abs(e.getValue());
            }
            return norm;
        }

        /// <summary>
        /// Returns the L<sub>&infin;</sub> norm of the vector.
        /// <para>The L^&infin; norm is the max of the absolute
        /// values of the elements.</para>
        /// </summary>
        /// <returns>the norm.</returns>
        /// <remarks>
        /// See <see cref="getNorm()"/><para/>
        /// See <see cref="getL1Norm()"/><para/>
        /// See <see cref="getLInfDistance(RealVector)"/>
        /// </remarks>
        public double getLInfNorm()
        {
            double norm = 0;
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                norm = FastMath.max(norm, FastMath.abs(e.getValue()));
            }
            return norm;
        }

        /// <summary>
        /// Distance between two vectors.
        /// <para>This method computes the distance consistent with
        /// L^1 norm, i.e. the sum of the absolute values of
        /// the elements differences.</para>
        /// </summary>
        /// <param name="v">Vector to which distance is requested.</param>
        /// <returns>the distance between two vectors.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public double getL1Distance(RealVector v)
        {
            checkVectorDimensions(v);
            double d = 0;
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                d += FastMath.abs(e.getValue() - v.getEntry(e.getIndex()));
            }
            return d;
        }

        /// <summary>
        /// Distance between two vectors.
        /// <para>This method computes the distance consistent with
        /// L^&infin; norm, i.e. the max of the absolute values of
        /// element differences.</para>
        /// </summary>
        /// <param name="v">Vector to which distance is requested.</param>
        /// <returns>the distance between two vectors.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        /// <remarks>
        /// See <see cref="getDistance(RealVector)"/><para/>
        /// See <see cref="getL1Distance(RealVector)"/><para/>
        /// See <see cref="getLInfNorm()"/> 
        /// </remarks>
        public double getLInfDistance(RealVector v)
        {
            checkVectorDimensions(v);
            double d = 0;
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                d = FastMath.max(FastMath.abs(e.getValue() - v.getEntry(e.getIndex())), d);
            }
            return d;
        }

        /// <summary>
        /// Get the index of the minimum entry.
        /// </summary>
        /// <returns>the index of the minimum entry or -1 if vector length is 0
        /// or all entries are <c>NaN</c>.</returns>
        public int getMinIndex()
        {
            int minIndex = -1;
            double minValue = Double.PositiveInfinity;
            IEnumerator<Entry> Iterator = iterator();
            while (Iterator.MoveNext())
            {
                Entry entry = Iterator.Current;
                if (entry.getValue() <= minValue)
                {
                    minIndex = entry.getIndex();
                    minValue = entry.getValue();
                }
            }
            return minIndex;
        }

        /// <summary>
        /// Get the value of the minimum entry.
        /// </summary>
        /// <returns>the value of the minimum entry or <c>NaN</c> if all
        /// entries are <c>NaN</c>.</returns>
        public double getMinValue()
        {
            int minIndex = getMinIndex();
            return minIndex < 0 ? Double.NaN : getEntry(minIndex);
        }

        /// <summary>
        /// Get the index of the maximum entry.
        /// </summary>
        /// <returns>the index of the maximum entry or -1 if vector length is 0
        /// or all entries are <c>NaN</c></returns>
        public int getMaxIndex()
        {
            int maxIndex = -1;
            double maxValue = Double.NegativeInfinity;
            IEnumerator<Entry> Iterator = iterator();
            while (Iterator.MoveNext())
            {
                Entry entry = Iterator.Current;
                if (entry.getValue() >= maxValue)
                {
                    maxIndex = entry.getIndex();
                    maxValue = entry.getValue();
                }
            }
            return maxIndex;
        }

        /// <summary>
        /// Get the value of the maximum entry.
        /// </summary>
        /// <returns>the value of the maximum entry or <c>NaN</c> if all
        /// entries are <c>NaN</c>.</returns>
        public double getMaxValue()
        {
            int maxIndex = getMaxIndex();
            return maxIndex < 0 ? Double.NaN : getEntry(maxIndex);
        }


        /// <summary>
        /// Multiply each entry by the argument. Returns a new vector.
        /// Does not change instance data.
        /// </summary>
        /// <param name="d">Multiplication factor.</param>
        /// <returns><c>this</c> * <c>d</c>.</returns>
        public RealVector mapMultiply(double d)
        {
            return copy().mapMultiplyToSelf(d);
        }

        /// <summary>
        /// Multiply each entry.
        /// The instance is changed in-place.
        /// </summary>
        /// <param name="d">Multiplication factor.</param>
        /// <returns><c>this</c>.</returns>
        public RealVector mapMultiplyToSelf(double d)
        {
            return mapToSelf(FunctionUtils.fix2ndArgument(new Multiply(), d));
        }

        /// <summary>
        /// Subtract a value from each entry. Returns a new vector.
        /// Does not change instance data.
        /// </summary>
        /// <param name="d">Value to be subtracted.</param>
        /// <returns><c>this</c> - <c>d</c>.</returns>
        public RealVector mapSubtract(double d)
        {
            return copy().mapSubtractToSelf(d);
        }

        /// <summary>
        /// Subtract a value from each entry.
        /// The instance is changed in-place.
        /// </summary>
        /// <param name="d">Value to be subtracted.</param>
        /// <returns><c>this</c>.</returns>
        public RealVector mapSubtractToSelf(double d)
        {
            return mapAddToSelf(-d);
        }

        /// <summary>
        /// Divide each entry by the argument. Returns a new vector.
        /// Does not change instance data.
        /// </summary>
        /// <param name="d">Value to divide by.</param>
        /// <returns><c>this</c> / <c>d</c>.</returns>
        public RealVector mapDivide(double d)
        {
            return copy().mapDivideToSelf(d);
        }

        /// <summary>
        /// Divide each entry by the argument.
        /// The instance is changed in-place.
        /// </summary>
        /// <param name="d">Value to divide by.</param>
        /// <returns><c>this</c>.</returns>
        public RealVector mapDivideToSelf(double d)
        {
            return mapToSelf(FunctionUtils.fix2ndArgument(new Divide(), d));
        }

        /// <summary>
        /// Compute the outer product.
        /// </summary>
        /// <param name="v">Vector with which outer product should be computed.</param>
        /// <returns>the matrix outer product between this instance and <c>v</c>.</returns>
        public RealMatrix outerProduct(RealVector v)
        {
            int m = this.getDimension();
            int n = v.getDimension();
            RealMatrix product;
            if (v is SparseRealVector || this is SparseRealVector)
            {
                product = new OpenMapRealMatrix(m, n);
            }
            else
            {
                product = new Array2DRowRealMatrix(m, n);
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    product.setEntry(i, j, this.getEntry(i) * v.getEntry(j));
                }
            }
            return product;
        }

        /// <summary>
        /// Find the orthogonal projection of this vector onto another vector.
        /// </summary>
        /// <param name="v">vector onto which instance must be projected.</param>
        /// <returns>projection of the instance onto <c>v</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>v</c> is not the same size as
        /// <c>this</c> vector.</exception>
        /// <exception cref="MathArithmeticException"> if <c>this</c> or <c>v</c> is the null
        /// vector</exception>
        public RealVector projection(RealVector v)
        {
            double norm2 = v.dotProduct(v);
            if (norm2 == 0.0)
            {
                throw new MathArithmeticException(new LocalizedFormats("ZERO_NORM"));
            }
            return v.mapMultiply(dotProduct(v) / v.dotProduct(v));
        }

        /// <summary>
        /// Set all elements to a single value.
        /// </summary>
        /// <param name="value">Single value to set for all elements.</param>
        public void set(double value)
        {
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                e.setValue(value);
            }
        }

        /// <summary>
        /// Convert the vector to an array of <c>double</c>s.
        /// The array is independent from this vector data: the elements
        /// are copied.
        /// </summary>
        /// <returns>an array containing a copy of the vector elements.</returns>
        public double[] toArray()
        {
            int dim = getDimension();
            double[] values = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                values[i] = getEntry(i);
            }
            return values;
        }

        /// <summary>
        /// Creates a unit vector pointing in the direction of this vector.
        /// The instance is not changed by this method.
        /// </summary>
        /// <returns>a unit vector pointing in direction of this vector.</returns>
        /// <exception cref="MathArithmeticException"> if the norm is zero.</exception>
        public RealVector unitVector()
        {
            double norm = getNorm();
            if (norm == 0)
            {
                throw new MathArithmeticException(new LocalizedFormats("ZERO_NORM"));
            }
            return mapDivide(norm);
        }

        /// <summary>
        /// Converts this vector into a unit vector.
        /// The instance itself is changed by this method.
        /// </summary>
        /// <exception cref="MathArithmeticException"> if the norm is zero.</exception>
        public void unitize()
        {
            double norm = getNorm();
            if (norm == 0)
            {
                throw new MathArithmeticException(new LocalizedFormats("ZERO_NORM"));
            }
            mapDivideToSelf(getNorm());
        }

        /// <summary>
        /// Create a sparse iterator over the vector, which may omit some entries.
        /// The ommitted entries are either exact zeroes (for dense implementations)
        /// or are the entries which are not stored (for real sparse vectors).
        /// No guarantees are made about order of iteration.
        /// <para>Note: derived classes are required to return an <see cref="IEnumerator"/> that
        /// returns non-null <see cref="Entry"/> objects as long as <see cref="IEnumerator.Next()"/>
        /// returns <c>true</c>.</para>
        /// </summary>
        /// <returns>a sparse iterator.</returns>
        public IEnumerator<Entry> sparseIterator()
        {
            return new SparseEntryIterator(this);
        }

        /// <summary>
        /// Generic dense iterator. Iteration is in increasing order
        /// of the vector index.
        /// <para>Note: derived classes are required to return an <see cref="IEnumerator"/> that
        /// returns non-null <see cref="Entry"/> objects as long as <see cref="IEnumerator.Next()"/>
        /// returns <c>true</c>.</para>
        /// </summary>
        /// <returns>a dense iterator.</returns>
        public IEnumerator<Entry> iterator()
        {
            int dim = getDimension();
            return new EntryIterator(dim, this);
        }

        private class EntryIterator : IEnumerator<Entry>
        {
            private Int32 dim;

            /// <summary>
            /// Current index.
            /// </summary>
            private int i = 0;

            /// <summary>
            /// Current entry.
            /// </summary>
            private Entry e;
            public EntryIterator(Int32 dim, RealVector v)
            {
                this.dim = dim;
                this.e = new Entry(v);
            }

            /// <inheritdoc/>
            public Boolean hasNext()
            {
                return i < dim;
            }

            /// <inheritdoc/>
            public Entry next()
            {
                e.setIndex(i++);
                return e;
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all circumstances.
            /// </exception>
            public void remove()
            {
                throw new MathUnsupportedOperationException();
            }

            public Entry Current
            {
                get { throw new NotImplementedException(); }
            }

            object System.Collections.IEnumerator.Current
            {
                get { return this; }
            }

            public bool MoveNext()
            {
                if (this.hasNext())
                {
                    this.next();
                    return true;
                }
                return false;
            }

            public void Reset()
            {
                return;
            }

            public void Dispose()
            {
                return;
            }
        }

        /// <summary>
        /// Acts as if implemented as:
        /// <code>
        ///  return copy().mapToSelf(function);
        /// </code>
        /// Returns a new vector. Does not change instance data.
        /// </summary>
        /// <param name="function">Function to apply to each entry.</param>
        /// <returns>a new vector.</returns>
        public RealVector map(UnivariateFunction function)
        {
            return copy().mapToSelf(function);
        }

        /// <summary>
        /// Acts as if it is implemented as:
        /// <code>
        ///  Entry e = null;
        ///  for(IEnumerator<Entry> it = iterator(); it.Next(); e = it.Current) {
        ///      e.setValue(function.value(e.getValue()));
        ///  }
        /// </code>
        /// Entries of this vector are modified in-place by this method.
        /// </summary>
        /// <param name="function">Function to apply to each entry.</param>
        /// <returns>a reference to this vector.</returns>
        public RealVector mapToSelf(UnivariateFunction function)
        {
            IEnumerator<Entry> it = iterator();
            while (it.MoveNext())
            {
                Entry e = it.Current;
                e.setValue(function.value(e.getValue()));
            }
            return this;
        }

        /// <summary>
        /// Returns a new vector representing <c>a * this + b * y</c>, the linear
        /// combination of <c>this</c> and <c>y</c>.
        /// Returns a new vector. Does not change instance data.
        /// </summary>
        /// <param name="a">Coefficient of <c>this</c>.</param>
        /// <param name="b">Coefficient of <c>y</c>.</param>
        /// <param name="y">Vector with which <c>this</c> is linearly combined</param>
        /// <returns>a vector containing <c>a * this[i] + b * y[i]</c> for all
        /// <c>i</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>y</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public RealVector combine(double a, double b, RealVector y)
        {
            return copy().combineToSelf(a, b, y);
        }

        /// <summary>
        /// Updates <c>this</c> with the linear combination of <c>this</c> and
        /// <c>y</c>.
        /// </summary>
        /// <param name="a">Weight of <c>this</c>.</param>
        /// <param name="b">Weight of <c>y</c>.</param>
        /// <param name="y">Vector with which <c>this</c> is linearly combined.</param>
        /// <returns><c>this</c>, with components equal to
        /// <c>a * this[i] + b * y[i]</c> for all <c>i</c>.</returns>
        /// <exception cref="DimensionMismatchException"> if <c>y</c> is not the same size as
        /// <c>this</c> vector.</exception>
        public RealVector combineToSelf(double a, double b, RealVector y)
        {
            checkVectorDimensions(y);
            for (int i = 0; i < getDimension(); i++)
            {
                double xi = getEntry(i);
                double yi = y.getEntry(i);
                setEntry(i, a * xi + b * yi);
            }
            return this;
        }

        /// <summary>
        /// Visits (but does not alter) all entries of this vector in default order
        /// (increasing index).
        /// </summary>
        /// <param name="visitor">the visitor to be used to process the entries of this
        /// vector</param>
        /// <returns>the value returned by <see cref="RealVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        public double walkInDefaultOrder(RealVectorPreservingVisitor visitor)
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
        /// <returns>the value returned by <see cref="RealVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public double walkInDefaultOrder(RealVectorPreservingVisitor visitor, int start, int end)
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
        /// <returns>the value returned by <see cref="RealVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        public double walkInOptimizedOrder(RealVectorPreservingVisitor visitor)
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
        /// <returns>the value returned by <see cref="RealVectorPreservingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public double walkInOptimizedOrder(RealVectorPreservingVisitor visitor, int start, int end)
        {
            return walkInDefaultOrder(visitor, start, end);
        }

        /// <summary>
        /// Visits (and possibly alters) all entries of this vector in default order
        /// (increasing index).
        /// </summary>
        /// <param name="visitor">the visitor to be used to process and modify the entries
        /// of this vector</param>
        /// <returns>the value returned by <see cref="RealVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        public double walkInDefaultOrder(RealVectorChangingVisitor visitor)
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
        /// <returns>value returned by <see cref="RealVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public double walkInDefaultOrder(RealVectorChangingVisitor visitor, int start, int end)
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
        /// <returns>the value returned by <see cref="RealVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        public double walkInOptimizedOrder(RealVectorChangingVisitor visitor)
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
        /// <returns>the value returned by <see cref="RealVectorChangingVisitor.end()"/>
        /// at the end of the walk</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>end < start</c>.</exception>
        /// <exception cref="OutOfRangeException"> if the indices are not valid.</exception>
        public double walkInOptimizedOrder(RealVectorChangingVisitor visitor, int start, int end)
        {
            return walkInDefaultOrder(visitor, start, end);
        }

        /// <summary>
        /// An entry in the vector.
        /// </summary>
        public class Entry
        {
            /// <summary>
            /// Index of this entry.
            /// </summary>
            private int index;
            protected RealVector v;

            /// <summary>
            /// Simple constructor.
            /// </summary>
            /// <param name="v"></param>
            public Entry(RealVector v)
            {
                this.v = v;
                setIndex(0);
            }

            /// <summary>
            /// Get the value of the entry.
            /// </summary>
            /// <returns>the value of the entry.</returns>
            public double getValue()
            {
                return v.getEntry(getIndex());
            }

            /// <summary>
            /// Set the value of the entry.
            /// </summary>
            /// <param name="value">New value for the entry.</param>
            public void setValue(double value)
            {
                v.setEntry(getIndex(), value);
            }

            /// <summary>
            /// Get the index of the entry.
            /// </summary>
            /// <returns>the index of the entry.</returns>
            public int getIndex()
            {
                return index;
            }

            /// <summary>
            /// Set the index of the entry.
            /// </summary>
            /// <param name="index">New index for the entry.</param>
            public void setIndex(int index)
            {
                this.index = index;
            }
        }

        /// <summary>
        /// <para>
        /// Test for the equality of two real vectors. If all coordinates of two real
        /// vectors are exactly the same, and none are <c>NaN</c>, the two real
        /// vectors are considered to be equal. <c>NaN</c> coordinates are
        /// considered to affect globally the vector and be equals to each other -
        /// i.e, if either (or all) coordinates of the real vector are equal to
        /// <c>NaN</c>, the real vector is equal to a vector with all <c>NaN</c>
        /// coordinates.
        /// </para>
        /// <para>
        /// This method <em>must</em> be overriden by concrete subclasses of
        /// <see cref="RealVector"/> (the current implementation throws an exception).
        /// </para>
        /// </summary>
        /// <param name="other">Object to test for equality.</param>
        /// <returns><c>true</c> if two vector objects are equal, <c>false</c> if
        /// <c>other</c> is null, not an instance of <c>RealVector</c>, or
        /// not equal to this <c>RealVector</c> instance.</returns>
        /// <exception cref="MathUnsupportedOperationException"> if this method is not
        /// overridden.</exception>
        public override Boolean Equals(Object other)
        {
            throw new MathUnsupportedOperationException();
        }

        /// <inheritdoc/>. 
        /// <remarks>
        /// This method <em>must</em> be overriden by concrete
        /// subclasses of <see cref="RealVector"/> (current implementation throws an
        /// exception).
        /// </remarks>
        /// <exception cref="MathUnsupportedOperationException"> if this method is not
        /// overridden.</exception>
        public override int GetHashCode()
        {
            throw new MathUnsupportedOperationException();
        }

        /// <summary>
        /// This class should rarely be used, but is here to provide
        /// a default implementation of sparseIterator(), which is implemented
        /// by walking over the entries, skipping those that are zero.
        /// Concrete subclasses which are SparseVector implementations should
        /// make their own sparse iterator, rather than using this one.
        /// This implementation might be useful for ArrayRealVector, when expensive
        /// operations which preserve the default value are to be done on the entries,
        /// and the fraction of non-default values is small (i.e. someone took a
        /// SparseVector, and passed it into the copy-constructor of ArrayRealVector)
        /// </summary>
        public class SparseEntryIterator : IEnumerator<Entry>
        {
            /// <summary>
            /// Dimension of the vector.
            /// </summary>
            private int dim;

            /// <summary>
            /// Last entry returned by <see cref="Next()"/>.
            /// </summary>
            private Entry current;

            /// <summary>
            /// Next entry for <see cref="Next()"/> to return.
            /// </summary>
            private Entry Next;
            private Boolean disposed;

            /// <summary>
            /// Simple constructor.
            /// </summary>
            /// <param name="v"></param>
            public SparseEntryIterator(RealVector v)
            {
                dim = v.getDimension();
                current = new Entry(v);
                Next = new Entry(v);
                if (Next.getValue() == 0)
                {
                    advance(Next);
                }
            }

            /// <summary>
            /// Advance an entry up to the next nonzero one.
            /// </summary>
            /// <param name="e">entry to advance.</param>
            public void advance(Entry e)
            {
                if (e == null)
                {
                    return;
                }
                do
                {
                    e.setIndex(e.getIndex() + 1);
                } while (e.getIndex() < dim && e.getValue() == 0);
                if (e.getIndex() >= dim)
                {
                    e.setIndex(-1);
                }
            }

            /// <inheritdoc/>
            public Boolean hasNext()
            {
                return Next.getIndex() >= 0;
            }

            /// <inheritdoc/>
            public Entry next()
            {
                int index = Next.getIndex();
                current.setIndex(index);
                advance(Next);
                return current;
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all circumstances.</exception>
            public void remove()
            {
                throw new MathUnsupportedOperationException();
            }

            public Entry Current
            {
                get { return this.Current; }
            }

            object System.Collections.IEnumerator.Current
            {
                get { return this.Current; }
            }

            public bool MoveNext()
            {
                if (this.hasNext())
                {
                    int index = Next.getIndex();
                    current.setIndex(index);
                    advance(Next);
                    return true;
                }
                return false;
            }

            public void Reset()
            {
                throw new NotImplementedException();
            }

            protected virtual void Dispose(Boolean disposing)
            {
                if (!this.disposed)
                {
                    if (disposing)
                    {
                        current = null;
                        Next = null;
                        disposed = true;
                    }
                }
            }

            public void Dispose()
            {
                Dispose(true);
                // This object will be cleaned up by the Dispose method.
                // Therefore, you should call GC.SupressFinalize to
                // take this object off the finalization queue
                // and prevent finalization code for this object
                // from executing a second time.
                GC.SuppressFinalize(this);
            }
        }

        /// <summary>
        /// Returns an unmodifiable view of the specified vector.
        /// The returned vector has read-only access. An attempt to modify it will
        /// result in a <see cref="MathUnsupportedOperationException"/>. However, the
        /// returned vector is <em>not</em> immutable, since any modification of
        /// <c>v</c> will also change the returned view.
        /// For example, in the following piece of code
        /// <code>
        ///     RealVector v = new ArrayRealVector(2);
        ///     RealVector w = RealVector.unmodifiableRealVector(v);
        ///     v.setEntry(0, 1.2);
        ///     v.setEntry(1, -3.4);
        /// </code>
        /// the changes will be seen in the <c>w</c> view of <c>v</c>.
        /// </summary>
        /// <param name="v">Vector for which an unmodifiable view is to be returned.</param>
        /// <returns>an unmodifiable view of <c>v</c>.</returns>
        public static RealVector unmodifiableRealVector(RealVector v)
        {
            return new RealVectorHelper(v);
        }

        /// <summary>
        /// This anonymous class is an implementation of <see cref="RealVector"/>
        /// with read-only access.
        /// It wraps any <see cref="RealVector"/>, and exposes all methods which
        /// do not modify it. Invoking methods which should normally result
        /// in the modification of the calling <see cref="RealVector"/> results in
        /// a <see cref="MathUnsupportedOperationException"/>. It should be noted
        /// that <see cref="UnmodifiableVector"/> is <em>not</em> immutable.
        /// </summary>
        private class RealVectorHelper : RealVector
        {
            private RealVector v;
            public RealVectorHelper(RealVector v)
            {
                this.v = v;
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all circumstances.
            /// </exception>
            public new RealVector mapToSelf(UnivariateFunction function)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public new RealVector map(UnivariateFunction function)
            {
                return v.map(function);
            }

            /// <inheritdoc/>
            public new IEnumerator<Entry> iterator()
            {
                IEnumerator<Entry> i = v.iterator();
                return new AnotherEntryIterator(i, v);
            }

            private class AnotherEntryIterator : IEnumerator<Entry>
            {
                /// <summary>
                /// The current entry.
                /// </summary>
                private readonly UnmodifiableEntry e;
                private IEnumerator<Entry> i;
                public AnotherEntryIterator(IEnumerator<Entry> i, RealVector v)
                {
                    this.i = i;
                    this.e = new UnmodifiableEntry(v);
                }

                /// <inheritdoc/>
                /// <exception cref="MathUnsupportedOperationException"> in all
                /// circumstances.</exception>
                public void remove()
                {
                    throw new MathUnsupportedOperationException();
                }

                public Entry Current
                {
                    get { return e; }
                }

                object System.Collections.IEnumerator.Current
                {
                    get { return e; }
                }

                public bool MoveNext()
                {
                    if (i.MoveNext())
                    {
                        e.setIndex(i.Current.getIndex());
                        return true;
                    }
                    return false;
                }

                public void Reset()
                {
                    throw new NotImplementedException();
                }

                public void Dispose()
                {
                    return;
                }
            }

            /// <inheritdoc/>
            public new IEnumerator<Entry> sparseIterator()
            {
                IEnumerator<Entry> i = v.sparseIterator();
                return new AnotherEntryIterator(i, v);
            }

            /// <inheritdoc/>
            public override RealVector copy()
            {
                return v.copy();
            }

            /// <inheritdoc/>
            public new RealVector add(RealVector w)
            {
                return v.add(w);
            }

            /// <inheritdoc/>
            public new RealVector subtract(RealVector w)
            {
                return v.subtract(w);
            }

            /// <inheritdoc/>
            public new RealVector mapAdd(double d)
            {
                return v.mapAdd(d);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new RealVector mapAddToSelf(double d)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public new RealVector mapSubtract(double d)
            {
                return v.mapSubtract(d);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new RealVector mapSubtractToSelf(double d)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public new RealVector mapMultiply(double d)
            {
                return v.mapMultiply(d);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new RealVector mapMultiplyToSelf(double d)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public new RealVector mapDivide(double d)
            {
                return v.mapDivide(d);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new RealVector mapDivideToSelf(double d)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public override RealVector ebeMultiply(RealVector w)
            {
                return v.ebeMultiply(w);
            }

            /// <inheritdoc/>
            public override RealVector ebeDivide(RealVector w)
            {
                return v.ebeDivide(w);
            }

            /// <inheritdoc/>
            public new double dotProduct(RealVector w)
            {
                return v.dotProduct(w);
            }

            /// <inheritdoc/>
            public new double cosine(RealVector w)
            {
                return v.cosine(w);
            }

            /// <inheritdoc/>
            public new double getNorm()
            {
                return v.getNorm();
            }

            /// <inheritdoc/>
            public new double getL1Norm()
            {
                return v.getL1Norm();
            }

            /// <inheritdoc/>
            public new double getLInfNorm()
            {
                return v.getLInfNorm();
            }

            /// <inheritdoc/>
            public new double getDistance(RealVector w)
            {
                return v.getDistance(w);
            }

            /// <inheritdoc/>
            public new double getL1Distance(RealVector w)
            {
                return v.getL1Distance(w);
            }

            /// <inheritdoc/>
            public new double getLInfDistance(RealVector w)
            {
                return v.getLInfDistance(w);
            }

            /// <inheritdoc/>
            public new RealVector unitVector()
            {
                return v.unitVector();
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new void unitize()
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public new RealMatrix outerProduct(RealVector w)
            {
                return v.outerProduct(w);
            }

            /// <inheritdoc/>
            public override double getEntry(int index)
            {
                return v.getEntry(index);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public override void setEntry(int index, double value)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new void addToEntry(int index, double value)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public override int getDimension()
            {
                return v.getDimension();
            }

            /// <inheritdoc/>
            public override RealVector append(RealVector w)
            {
                return v.append(w);
            }

            /// <inheritdoc/>
            public override RealVector append(double d)
            {
                return v.append(d);
            }

            /// <inheritdoc/>
            public override RealVector getSubVector(int index, int n)
            {
                return v.getSubVector(index, n);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public override void setSubVector(int index, RealVector w)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new void set(double value)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <inheritdoc/>
            public new double[] toArray()
            {
                return v.toArray();
            }

            /// <inheritdoc/>
            public override Boolean isNaN()
            {
                return v.isNaN();
            }

            /// <inheritdoc/>
            public override Boolean isInfinite()
            {
                return v.isInfinite();
            }

            /// <inheritdoc/>
            public new RealVector combine(double a, double b, RealVector y)
            {
                return v.combine(a, b, y);
            }

            /// <inheritdoc/>
            /// <exception cref="MathUnsupportedOperationException"> in all
            /// circumstances.</exception>
            public new RealVector combineToSelf(double a, double b, RealVector y)
            {
                throw new MathUnsupportedOperationException();
            }

            /// <summary>
            /// An entry in the vector.
            /// </summary>
            internal class UnmodifiableEntry : Entry
            {
                internal UnmodifiableEntry(RealVector v) : base(v) { }
                /// <inheritdoc/>
                public new double getValue()
                {
                    return v.getEntry(getIndex());
                }

                /// <inheritdoc/>
                /// <exception cref="MathUnsupportedOperationException"> in all
                /// circumstances.</exception>
                public new void setValue(double value)
                {
                    throw new MathUnsupportedOperationException();
                }
            }
        }
    }
}