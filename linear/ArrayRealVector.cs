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
using Math3.exception;
using Math3.exception.util;
using Math3.util;
using System;
using System.Collections.Generic;

namespace Math3.linear
{
    /// <summary>
    /// This class implements the <see cref="RealVector"/> interface with a double array.
    /// </summary>
    public class ArrayRealVector : RealVector
    {
        /// <summary>
        /// Default format.
        /// </summary>
        private static readonly RealVectorFormat DEFAULT_FORMAT = RealVectorFormat.getInstance();

        /// <summary>
        /// Entries of the vector.
        /// </summary>
        private double[] data;

        /// <summary>
        /// Build a 0-length vector.
        /// Zero-length vectors may be used to initialized construction of vectors
        /// by data gathering. We start with zero-length and use either the
        /// <see cref="ArrayRealVector(ArrayRealVector, ArrayRealVector)">constructor</see>
        /// or one of the <c>append</c> method (<see cref="append(double)"/>,
        /// <see cref="append(ArrayRealVector)"/>) to gather data into this vector.
        /// </summary>
        public ArrayRealVector()
        {
            data = new double[0];
        }

        /// <summary>
        /// Construct a vector of zeroes. 
        /// </summary>
        /// <param name="size">Size of the vector.</param>
        public ArrayRealVector(int size)
        {
            data = new double[size];
        }

        /// <summary>
        /// Construct a vector with preset values.
        /// </summary>
        /// <param name="size">Size of the vector</param>
        /// <param name="preset">All entries will be set with this value.</param>
        public ArrayRealVector(int size, double preset)
        {
            data = new double[size];
            for (int i = 0; i < data.Length; ++i)
            {
                data[i] = preset;
            }
        }

        /// <summary>
        /// Construct a vector from an array, copying the input array.
        /// </summary>
        /// <param name="d">Array.</param>
        public ArrayRealVector(double[] d)
        {
            data = (Double[])d.Clone();
        }

        /// <summary>
        /// Create a new ArrayRealVector using the input array as the underlying
        /// data array.
        /// If an array is built specially in order to be embedded in a
        /// ArrayRealVector and not used directly, the <c>copyArray</c> may be
        /// set to <c>false</c>. This will prevent the copying and improve
        /// performance as no new array will be built and no data will be copied.
        /// </summary>
        /// <param name="d">Data for the new vector.</param>
        /// <param name="copyArray">if <c>true</c>, the input array will be copied,
        /// otherwise it will be referenced.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <remarks>
        /// See <see cref="ArrayRealVector(double[])"/>
        /// </remarks>
        public ArrayRealVector(double[] d, Boolean copyArray)
        {
            if (d == null)
            {
                throw new NullArgumentException();
            }
            data = copyArray ? (Double[])d.Clone() : d;
        }

        /// <summary>
        /// Construct a vector from part of a array.
        /// </summary>
        /// <param name="d">Array.</param>
        /// <param name="pos">Position of first entry.</param>
        /// <param name="size">Number of entries to copy.</param>
        /// <exception cref="NullArgumentException"> if <c>d</c> is <c>null</c>.</exception>
        /// <exception cref="NumberIsTooLargeException"> if the size of <c>d</c> is less
        /// than <c>pos + size</c>.</exception>
        public ArrayRealVector(double[] d, int pos, int size)
        {
            if (d == null)
            {
                throw new NullArgumentException();
            }
            if (d.Length < pos + size)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(pos + size, d.Length, true);
            }
            data = new double[size];
            Array.Copy(d, pos, data, 0, size);
        }

        /// <summary>
        /// Construct a vector from another vector, using a deep copy.
        /// </summary>
        /// <param name="v">vector to copy.</param>
        /// <exception cref="NullArgumentException"> if <c>v</c> is <c>null</c>.</exception>
        public ArrayRealVector(RealVector v)
        {
            if (v == null)
            {
                throw new NullArgumentException();
            }
            data = new double[v.getDimension()];
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
        public ArrayRealVector(ArrayRealVector v) : this(v, true) { }

        /// <summary>
        /// Construct a vector from another vector.
        /// </summary>
        /// <param name="v">Vector to copy.</param>
        /// <param name="deep">If <c>true</c> perform a deep copy, otherwise perform a
        /// shallow copy.</param>
        public ArrayRealVector(ArrayRealVector v, Boolean deep)
        {
            data = deep ? (Double[])v.data.Clone() : v.data;
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        public ArrayRealVector(ArrayRealVector v1, ArrayRealVector v2)
        {
            data = new double[v1.data.Length + v2.data.Length];
            Array.Copy(v1.data, 0, data, 0, v1.data.Length);
            Array.Copy(v2.data, 0, data, v1.data.Length, v2.data.Length);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        public ArrayRealVector(ArrayRealVector v1, RealVector v2)
        {
            int l1 = v1.data.Length;
            int l2 = v2.getDimension();
            data = new double[l1 + l2];
            Array.Copy(v1.data, 0, data, 0, l1);
            for (int i = 0; i < l2; ++i)
            {
                data[l1 + i] = v2.getEntry(i);
            }
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        public ArrayRealVector(RealVector v1, ArrayRealVector v2)
        {
            int l1 = v1.getDimension();
            int l2 = v2.data.Length;
            data = new double[l1 + l2];
            for (int i = 0; i < l1; ++i)
            {
                data[i] = v1.getEntry(i);
            }
            Array.Copy(v2.data, 0, data, l1, l2);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        public ArrayRealVector(ArrayRealVector v1, double[] v2)
        {
            int l1 = v1.getDimension();
            int l2 = v2.Length;
            data = new double[l1 + l2];
            Array.Copy(v1.data, 0, data, 0, l1);
            Array.Copy(v2, 0, data, l1, l2);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        public ArrayRealVector(double[] v1, ArrayRealVector v2)
        {
            int l1 = v1.Length;
            int l2 = v2.getDimension();
            data = new double[l1 + l2];
            Array.Copy(v1, 0, data, 0, l1);
            Array.Copy(v2.data, 0, data, l1, l2);
        }

        /// <summary>
        /// Construct a vector by appending one vector to another vector.
        /// </summary>
        /// <param name="v1">First vector (will be put in front of the new vector).</param>
        /// <param name="v2">Second vector (will be put at back of the new vector).</param>
        public ArrayRealVector(double[] v1, double[] v2)
        {
            int l1 = v1.Length;
            int l2 = v2.Length;
            data = new double[l1 + l2];
            Array.Copy(v1, 0, data, 0, l1);
            Array.Copy(v2, 0, data, l1, l2);
        }

        /// <inheritdoc/>
        public override RealVector copy()
        {
            return new ArrayRealVector(this, true);
        }

        /// <inheritdoc/>
        public new ArrayRealVector add(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                int dim = vData.Length;
                checkVectorDimensions(dim);
                ArrayRealVector result = new ArrayRealVector(dim);
                double[] resultData = result.data;
                for (int i = 0; i < dim; i++)
                {
                    resultData[i] = data[i] + vData[i];
                }
                return result;
            }
            else
            {
                checkVectorDimensions(v);
                double[] outp = (Double[])data.Clone();
                IEnumerator<Entry> it = v.iterator();
                while (it.MoveNext())
                {
                    Entry e = it.Current;
                    outp[e.getIndex()] += e.getValue();
                }
                return new ArrayRealVector(outp, false);
            }
        }

        /// <inheritdoc/>
        public new ArrayRealVector subtract(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                int dim = vData.Length;
                checkVectorDimensions(dim);
                ArrayRealVector result = new ArrayRealVector(dim);
                double[] resultData = result.data;
                for (int i = 0; i < dim; i++)
                {
                    resultData[i] = data[i] - vData[i];
                }
                return result;
            }
            else
            {
                checkVectorDimensions(v);
                double[] outp = (Double[])data.Clone();
                IEnumerator<Entry> it = v.iterator();
                while (it.MoveNext())
                {
                    Entry e = it.Current;
                    outp[e.getIndex()] -= e.getValue();
                }
                return new ArrayRealVector(outp, false);
            }
        }

        /// <inheritdoc/>
        public new ArrayRealVector map(UnivariateFunction function)
        {
            return (ArrayRealVector)copy().mapToSelf(function);
        }

        /// <inheritdoc/>
        public new ArrayRealVector mapToSelf(UnivariateFunction function)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = function.value(data[i]);
            }
            return this;
        }

        /// <inheritdoc/>
        public new RealVector mapAddToSelf(double d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] += d;
            }
            return this;
        }

        /// <inheritdoc/>
        public new RealVector mapSubtractToSelf(double d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] -= d;
            }
            return this;
        }

        /// <inheritdoc/>
        public new RealVector mapMultiplyToSelf(double d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] *= d;
            }
            return this;
        }

        /// <inheritdoc/>
        public new RealVector mapDivideToSelf(double d)
        {
            for (int i = 0; i < data.Length; i++)
            {
                data[i] /= d;
            }
            return this;
        }

        /// <inheritdoc/>
        public override RealVector ebeMultiply(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                int dim = vData.Length;
                checkVectorDimensions(dim);
                ArrayRealVector result = new ArrayRealVector(dim);
                double[] resultData = result.data;
                for (int i = 0; i < dim; i++)
                {
                    resultData[i] = data[i] * vData[i];
                }
                return result;
            }
            else
            {
                checkVectorDimensions(v);
                double[] outp = (Double[])data.Clone();
                for (int i = 0; i < data.Length; i++)
                {
                    outp[i] *= v.getEntry(i);
                }
                return new ArrayRealVector(outp, false);
            }
        }

        /// <inheritdoc/>
        public override RealVector ebeDivide(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                int dim = vData.Length;
                checkVectorDimensions(dim);
                ArrayRealVector result = new ArrayRealVector(dim);
                double[] resultData = result.data;
                for (int i = 0; i < dim; i++)
                {
                    resultData[i] = data[i] / vData[i];
                }
                return result;
            }
            else
            {
                checkVectorDimensions(v);
                double[] outp = (Double[])data.Clone();
                for (int i = 0; i < data.Length; i++)
                {
                    outp[i] /= v.getEntry(i);
                }
                return new ArrayRealVector(outp, false);
            }
        }

        /// <summary>
        /// Get a reference to the underlying data array.
        /// This method does not make a fresh copy of the underlying data. 
        /// </summary>
        /// <returns>the array of entries.</returns>
        public double[] getDataRef()
        {
            return data;
        }

        /// <inheritdoc/>
        public new double dotProduct(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                checkVectorDimensions(vData.Length);
                double dot = 0;
                for (int i = 0; i < data.Length; i++)
                {
                    dot += data[i] * vData[i];
                }
                return dot;
            }
            return base.dotProduct(v);
        }

        /// <inheritdoc/>
        public new double getNorm()
        {
            double sum = 0;
            foreach (double a in data)
            {
                sum += a * a;
            }
            return FastMath.sqrt(sum);
        }

        /// <inheritdoc/>
        public new double getL1Norm()
        {
            double sum = 0;
            foreach (double a in data)
            {
                sum += FastMath.abs(a);
            }
            return sum;
        }

        /// <inheritdoc/>
        public new double getLInfNorm()
        {
            double max = 0;
            foreach (double a in data)
            {
                max = FastMath.max(max, FastMath.abs(a));
            }
            return max;
        }

        /// <inheritdoc/>
        public new double getDistance(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                checkVectorDimensions(vData.Length);
                double sum = 0;
                for (int i = 0; i < data.Length; ++i)
                {
                    double delta = data[i] - vData[i];
                    sum += delta * delta;
                }
                return FastMath.sqrt(sum);
            }
            else
            {
                checkVectorDimensions(v);
                double sum = 0;
                for (int i = 0; i < data.Length; ++i)
                {
                    double delta = data[i] - v.getEntry(i);
                    sum += delta * delta;
                }
                return FastMath.sqrt(sum);
            }
        }

        /// <inheritdoc/>
        public new double getL1Distance(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                checkVectorDimensions(vData.Length);
                double sum = 0;
                for (int i = 0; i < data.Length; ++i)
                {
                    double delta = data[i] - vData[i];
                    sum += FastMath.abs(delta);
                }
                return sum;
            }
            else
            {
                checkVectorDimensions(v);
                double sum = 0;
                for (int i = 0; i < data.Length; ++i)
                {
                    double delta = data[i] - v.getEntry(i);
                    sum += FastMath.abs(delta);
                }
                return sum;
            }
        }

        /// <inheritdoc/>
        public new double getLInfDistance(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                checkVectorDimensions(vData.Length);
                double max = 0;
                for (int i = 0; i < data.Length; ++i)
                {
                    double delta = data[i] - vData[i];
                    max = FastMath.max(max, FastMath.abs(delta));
                }
                return max;
            }
            else
            {
                checkVectorDimensions(v);
                double max = 0;
                for (int i = 0; i < data.Length; ++i)
                {
                    double delta = data[i] - v.getEntry(i);
                    max = FastMath.max(max, FastMath.abs(delta));
                }
                return max;
            }
        }

        /// <inheritdoc/>
        public new RealMatrix outerProduct(RealVector v)
        {
            if (v is ArrayRealVector)
            {
                double[] vData = ((ArrayRealVector)v).data;
                int m = data.Length;
                int n = vData.Length;
                RealMatrix outp = MatrixUtils.createRealMatrix(m, n);
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        outp.setEntry(i, j, data[i] * vData[j]);
                    }
                }
                return outp;
            }
            else
            {
                int m = data.Length;
                int n = v.getDimension();
                RealMatrix outp = MatrixUtils.createRealMatrix(m, n);
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        outp.setEntry(i, j, data[i] * v.getEntry(j));
                    }
                }
                return outp;
            }
        }

        /// <inheritdoc/>
        public override double getEntry(int index)
        {
            try
            {
                return data[index];
            }
            catch (IndexOutOfRangeException)
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("INDEX"), index, 0, getDimension() - 1);
            }
        }

        /// <inheritdoc/>
        public override int getDimension()
        {
            return data.Length;
        }

        /// <inheritdoc/>
        public override RealVector append(RealVector v)
        {
            try
            {
                return new ArrayRealVector(this, (ArrayRealVector)v);
            }
            catch (InvalidCastException)
            {
                return new ArrayRealVector(this, v);
            }
        }

        /// <summary>
        /// Construct a vector by appending a vector to this vector.
        /// </summary>
        /// <param name="v">Vector to append to this one.</param>
        /// <returns>a new vector.</returns>
        public ArrayRealVector append(ArrayRealVector v)
        {
            return new ArrayRealVector(this, v);
        }

        /// <inheritdoc/>
        public override RealVector append(double inp)
        {
            double[] outp = new double[data.Length + 1];
            Array.Copy(data, 0, outp, 0, data.Length);
            outp[data.Length] = inp;
            return new ArrayRealVector(outp, false);
        }

        /// <inheritdoc/>
        public override RealVector getSubVector(int index, int n)
        {
            if (n < 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE"), n);
            }
            ArrayRealVector outp = new ArrayRealVector(n);
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
        public override void setEntry(int index, double value)
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
        public new void addToEntry(int index, double increment)
        {
            try
            {
                data[index] += increment;
            }
            catch (IndexOutOfRangeException)
            {
                throw new OutOfRangeException<Int32>(new LocalizedFormats("INDEX"), index, 0, data.Length - 1);
            }
        }

        /// <inheritdoc/>
        public override void setSubVector(int index, RealVector v)
        {
            if (v is ArrayRealVector)
            {
                setSubVector(index, ((ArrayRealVector)v).data);
            }
            else
            {
                try
                {
                    for (int i = index; i < index + v.getDimension(); ++i)
                    {
                        data[i] = v.getEntry(i - index);
                    }
                }
                catch (IndexOutOfRangeException)
                {
                    checkIndex(index);
                    checkIndex(index + v.getDimension() - 1);
                }
            }
        }

        /// <summary>
        /// Set a set of consecutive elements.
        /// </summary>
        /// <param name="index">Index of first element to be set.</param>
        /// <param name="v">Vector containing the values to set.</param>
        /// <exception cref="OutOfRangeException"> if the index is inconsistent with the vector
        /// size.</exception>
        public void setSubVector(int index, double[] v)
        {
            try
            {
                Array.Copy(v, 0, data, index, v.Length);
            }
            catch (IndexOutOfRangeException)
            {
                checkIndex(index);
                checkIndex(index + v.Length - 1);
            }
        }

        /// <inheritdoc/>
        public new void set(double value)
        {
            for (int i = 0; i < data.Length; ++i)
            {
                data[i] = value;
            }
        }

        /// <inheritdoc/>
        public new double[] toArray()
        {
            return (Double[])data.Clone();
        }

        /// <inheritdoc/>
        public override String ToString()
        {
            return DEFAULT_FORMAT.format(this);
        }

        /// <summary>
        /// Check if instance and specified vectors have the same dimension.
        /// </summary>
        /// <param name="v">Vector to compare instance with.</param>
        /// <exception cref="DimensionMismatchException"> if the vectors do not
        /// have the same dimension.</exception>
        protected new void checkVectorDimensions(RealVector v)
        {
            checkVectorDimensions(v.getDimension());
        }

        /// <summary>
        /// Check if instance dimension is equal to some expected value.
        /// </summary>
        /// <param name="n">Expected dimension.</param>
        /// <exception cref="DimensionMismatchException"> if the dimension is
        /// inconsistent with vector size.</exception>
        protected new void checkVectorDimensions(int n)
        {
            if (data.Length != n)
            {
                throw new DimensionMismatchException(data.Length, n);
            }
        }

        /// <summary>
        /// Check if any coordinate of this vector is <c>NaN</c>.
        /// </summary>
        /// <returns><c>true</c> if any coordinate of this vector is <c>NaN</c>,
        /// <c>false</c> otherwise.</returns>
        public override Boolean isNaN()
        {
            foreach (double v in data)
            {
                if (Double.IsNaN(v))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// Check whether any coordinate of this vector is infinite and none
        /// are <c>NaN</c>.
        /// </summary>
        /// <returns><c>true</c> if any coordinate of this vector is infinite and
        /// none are <c>NaN</c>, <c>false</c> otherwise.</returns>
        public override Boolean isInfinite()
        {
            if (isNaN())
            {
                return false;
            }

            foreach (double v in data)
            {
                if (Double.IsInfinity(v))
                {
                    return true;
                }
            }

            return false;
        }

        /// <inheritdoc/>
        public override Boolean Equals(Object other)
        {
            if (this == other)
            {
                return true;
            }

            if (!(other is RealVector))
            {
                return false;
            }

            RealVector rhs = (RealVector)other;
            if (data.Length != rhs.getDimension())
            {
                return false;
            }

            if (rhs.isNaN())
            {
                return this.isNaN();
            }

            for (int i = 0; i < data.Length; ++i)
            {
                if (data[i] != rhs.getEntry(i))
                {
                    return false;
                }
            }
            return true;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// All <c>NaN</c> values have the same hash code.
        /// </remarks>
        public override int GetHashCode()
        {
            if (isNaN())
            {
                return 9;
            }
            return MathUtils.hash(data);
        }

        /// <inheritdoc/>
        public new ArrayRealVector combine(double a, double b, RealVector y)
        {
            return (ArrayRealVector)copy().combineToSelf(a, b, y);
        }

        /// <inheritdoc/>
        public new ArrayRealVector combineToSelf(double a, double b, RealVector y)
        {
            if (y is ArrayRealVector)
            {
                double[] yData = ((ArrayRealVector)y).data;
                checkVectorDimensions(yData.Length);
                for (int i = 0; i < this.data.Length; i++)
                {
                    data[i] = a * data[i] + b * yData[i];
                }
            }
            else
            {
                checkVectorDimensions(y);
                for (int i = 0; i < this.data.Length; i++)
                {
                    data[i] = a * data[i] + b * y.getEntry(i);
                }
            }
            return this;
        }

        /// <inheritdoc/>
        public new double walkInDefaultOrder(RealVectorPreservingVisitor visitor)
        {
            visitor.start(data.Length, 0, data.Length - 1);
            for (int i = 0; i < data.Length; i++)
            {
                visitor.visit(i, data[i]);
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInDefaultOrder(RealVectorPreservingVisitor visitor, int start, int end)
        {
            checkIndices(start, end);
            visitor.start(data.Length, start, end);
            for (int i = start; i <= end; i++)
            {
                visitor.visit(i, data[i]);
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        /// <remarks>
        /// In this implementation, the optimized order is the default order.
        /// </remarks>
        public new double walkInOptimizedOrder(RealVectorPreservingVisitor visitor)
        {
            return walkInDefaultOrder(visitor);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// In this implementation, the optimized order is the default order.
        /// </remarks>
        public new double walkInOptimizedOrder(RealVectorPreservingVisitor visitor, int start, int end)
        {
            return walkInDefaultOrder(visitor, start, end);
        }

        /// <inheritdoc/>
        public new double walkInDefaultOrder(RealVectorChangingVisitor visitor)
        {
            visitor.start(data.Length, 0, data.Length - 1);
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = visitor.visit(i, data[i]);
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        public new double walkInDefaultOrder(RealVectorChangingVisitor visitor, int start, int end)
        {
            checkIndices(start, end);
            visitor.start(data.Length, start, end);
            for (int i = start; i <= end; i++)
            {
                data[i] = visitor.visit(i, data[i]);
            }
            return visitor.end();
        }

        /// <inheritdoc/>
        /// <remarks>
        /// In this implementation, the optimized order is the default order.
        /// </remarks>
        public new double walkInOptimizedOrder(RealVectorChangingVisitor visitor)
        {
            return walkInDefaultOrder(visitor);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// In this implementation, the optimized order is the default order.
        /// </remarks>
        public new double walkInOptimizedOrder(RealVectorChangingVisitor visitor, int start, int end)
        {
            return walkInDefaultOrder(visitor, start, end);
        }
    }
}
