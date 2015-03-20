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
using Math3.util;
using System;

namespace Math3.analysis.differentiation
{
    /// <summary>
    /// Class representing both the value and the differentials of a function.
    /// <para>This class is the workhorse of the differentiation package.</para>
    /// <para>This class is an implementation of the extension to Rall's
    /// numbers described in Dan Kalman's paper <a
    /// href="http://www1.american.edu/cas/mathstat/People/kalman/pdffiles/mmgautodiff.pdf">Doubly
    /// Recursive Multivariate Automatic Differentiation</a>, Mathematics Magazine, vol. 75,
    /// no. 3, June 2002.</para>. Rall's numbers are an extension to the real numbers used
    /// throughout mathematical expressions; they hold the derivative together with the
    /// value of a function. Dan Kalman's derivative structures hold all partial derivatives
    /// up to any specified order, with respect to any number of free parameters. Rall's
    /// numbers therefore can be seen as derivative structures for order one derivative and
    /// one free parameter, and real numbers can be seen as derivative structures with zero
    /// order derivative and no free parameters.</para>
    /// <para><see cref="DerivativeStructure"/> instances can be used directly thanks to
    /// the arithmetic operators to the mathematical functions provided as
    /// methods by this class (+, -, *, /, %, sin, cos ...).</para>
    /// <para>Implementing complex expressions by hand using these classes is
    /// a tedious and error-prone task but has the advantage of having no limitation
    /// on the derivation order despite no requiring users to compute the derivatives by
    /// themselves. Implementing complex expression can also be done by developing computation
    /// code using standard primitive double values and to use 
    /// <see cref="UnivariateFunctionDifferentiator">differentiators</see> to create the
    /// <see cref="DerivativeStructure"/>-based instances. This method is simpler but may
    /// be limited in the accuracy and derivation orders and may be computationally intensive
    /// (this is typically the case for <see cref="FiniteDifferencesDifferentiator">
    /// finite differences differentiator</see>.</para>
    /// <para>Instances of this class are guaranteed to be immutable.</para>
    /// </summary>
    /// <remarks>
    /// See <see cref="DSCompiler"/>
    /// </remarks>
    public class DerivativeStructure : RealFieldElement<DerivativeStructure>
    {
        /// <summary>
        /// Compiler for the current dimensions.
        /// </summary>
        private DSCompiler compiler;

        /// <summary>
        /// Combined array holding all values.
        /// </summary>
        private readonly double[] data;

        /// <summary>
        /// Build an instance with all values and derivatives set to 0.
        /// </summary>
        /// <param name="compiler">compiler to use for computation</param>
        private DerivativeStructure(DSCompiler compiler)
        {
            this.compiler = compiler;
            this.data = new double[compiler.getSize()];
        }

        /// <summary>
        /// Build an instance with all values and derivatives set to 0.
        /// </summary>
        /// <param name="parameters">number of free parameters</param>
        /// <param name="order">derivation order</param>
        /// <exception cref="NumberIsTooLargeException"> if order is too large</exception>
        public DerivativeStructure(int parameters, int order) : this(DSCompiler.getCompiler(parameters, order)) { }

        /// <summary>
        /// Build an instance representing a constant value.
        /// </summary>
        /// <param name="parameters">number of free parameters</param>
        /// <param name="order">derivation order</param>
        /// <param name="value">value of the constant</param>
        /// <exception cref="NumberIsTooLargeException"> if order is too large</exception>
        /// <remarks>
        /// See <see cref="DerivativeStructure(int, int, int, double)"/>
        /// </remarks>
        public DerivativeStructure(int parameters, int order, double value)
            : this(parameters, order)
        {
            this.data[0] = value;
        }

        /// <summary>
        /// Build an instance representing a variable.
        /// <para>Instances built using this constructor are considered
        /// to be the free variables with respect to which differentials
        /// are computed. As such, their differential with respect to
        /// themselves is +1.</para>
        /// </summary>
        /// <param name="parameters">number of free parameters</param>
        /// <param name="order">derivation order</param>
        /// <param name="index">index of the variable (from 0 to <c>parameters - 1</c>)</param>
        /// <param name="value">value of the variable</param>
        /// <exception cref="NumberIsTooLargeException"> if <c>index >= parameters</c>.</exception>
        /// <remarks>
        /// See <see cref="DerivativeStructure(int, int, double)"/>
        /// </remarks>
        public DerivativeStructure(int parameters, int order, int index, double value)
            : this(parameters, order, value)
        {
            if (index >= parameters)
            {
                throw new NumberIsTooLargeException<Int32, Int32>(index, parameters, false);
            }

            if (order > 0)
            {
                // the derivative of the variable with respect to itself is 1.
                data[DSCompiler.getCompiler(index, order).getSize()] = 1.0;
            }

        }

        /// <summary>
        /// Linear combination constructor.
        /// The derivative structure built will be a1 * ds1 + a2 * ds2
        /// </summary>
        /// <param name="a1">first scale factor</param>
        /// <param name="ds1">first base (unscaled) derivative structure</param>
        /// <param name="a2">second scale factor</param>
        /// <param name="ds2">second base (unscaled) derivative structure</param>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders are inconsistent</exception>
        public DerivativeStructure(double a1, DerivativeStructure ds1, double a2, DerivativeStructure ds2)
            : this(ds1.compiler)
        {
            compiler.checkCompatibility(ds2.compiler);
            compiler.linearCombination(a1, ds1.data, 0, a2, ds2.data, 0, data, 0);
        }

        /// <summary>
        /// Linear combination constructor.
        /// The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3
        /// </summary>
        /// <param name="a1">first scale factor</param>
        /// <param name="ds1">first base (unscaled) derivative structure</param>
        /// <param name="a2">second scale factor</param>
        /// <param name="ds2">second base (unscaled) derivative structure</param>
        /// <param name="a3">third scale factor</param>
        /// <param name="ds3">third base (unscaled) derivative structure</param>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders are inconsistent</exception>
        public DerivativeStructure(double a1, DerivativeStructure ds1, double a2, DerivativeStructure ds2, double a3, DerivativeStructure ds3)
            : this(ds1.compiler)
        {
            compiler.checkCompatibility(ds2.compiler);
            compiler.checkCompatibility(ds3.compiler);
            compiler.linearCombination(a1, ds1.data, 0, a2, ds2.data, 0, a3, ds3.data, 0, data, 0);
        }

        /// <summary>
        /// Linear combination constructor.
        /// The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
        /// </summary>
        /// <param name="a1">first scale factor</param>
        /// <param name="ds1">first base (unscaled) derivative structure</param>
        /// <param name="a2">second scale factor</param>
        /// <param name="ds2">second base (unscaled) derivative structure</param>
        /// <param name="a3">third scale factor</param>
        /// <param name="ds3">third base (unscaled) derivative structure</param>
        /// <param name="a4">fourth scale factor</param>
        /// <param name="ds4">fourth base (unscaled) derivative structure</param>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders are inconsistent</exception>
        public DerivativeStructure(double a1, DerivativeStructure ds1, double a2, DerivativeStructure ds2, double a3, DerivativeStructure ds3, double a4, DerivativeStructure ds4)
            : this(ds1.compiler)
        {
            compiler.checkCompatibility(ds2.compiler);
            compiler.checkCompatibility(ds3.compiler);
            compiler.checkCompatibility(ds4.compiler);
            compiler.linearCombination(a1, ds1.data, 0, a2, ds2.data, 0,
                                       a3, ds3.data, 0, a4, ds4.data, 0,
                                       data, 0);
        }

        /// <summary>
        /// Build an instance from all its derivatives.
        /// </summary>
        /// <param name="parameters">number of free parameters</param>
        /// <param name="order">derivation order</param>
        /// <param name="derivatives">derivatives sorted according to
        /// <see cref="DSCompiler.getPartialDerivativeIndex(params int)"/></param>
        /// <exception cref="DimensionMismatchException"> if derivatives array does not match the
        /// <see cref="DSCompiler.getSize()">size</see> expected by the compiler</exception>
        /// <exception cref="NumberIsTooLargeException"> if order is too large</exception>
        /// <remarks>
        /// See <see cref="getAllDerivatives()"/>
        /// </remarks>
        public DerivativeStructure(int parameters, int order, params double[] derivatives)
            : this(parameters, order)
        {
            if (derivatives.Length != data.Length)
            {
                throw new DimensionMismatchException(derivatives.Length, data.Length);
            }
            Array.Copy(derivatives, 0, data, 0, data.Length);
        }

        /// <summary>
        /// Copy constructor.
        /// </summary>
        /// <param name="ds">instance to copy</param>
        private DerivativeStructure(DerivativeStructure ds)
        {
            this.compiler = ds.compiler;
            this.data = (Double[])ds.data.Clone();
        }

        /// <summary>
        /// Get the number of free parameters.
        /// </summary>
        /// <returns>number of free parameters</returns>
        public int getFreeParameters()
        {
            return compiler.getFreeParameters();
        }

        /// <summary>
        /// Get the derivation order.
        /// </summary>
        /// <returns>derivation order</returns>
        public int getOrder()
        {
            return compiler.getOrder();
        }

        /// <summary>
        /// Create a constant compatible with instance order and number of parameters.
        /// <para>
        /// This method is a convenience factory method, it simply calls
        /// <c>new DerivativeStructure(getFreeParameters(), getOrder(), c)</c>
        /// </para>
        /// </summary>
        /// <param name="c">value of the constant</param>
        /// <returns>a constant compatible with instance order and number of parameters</returns>
        /// <remarks>See <see cref="DerivativeStructure(int, int, double)"/></remarks>
        public DerivativeStructure createConstant(double c)
        {
            return new DerivativeStructure(getFreeParameters(), getOrder(), c);
        }

        /// <inheritdoc/>
        public double getReal()
        {
            return data[0];
        }

        /// <summary>
        /// Get the value part of the derivative structure.
        /// </summary>
        /// <returns>value part of the derivative structure</returns>
        /// <remarks>
        /// See <see cref="getPartialDerivative(params int)"/>
        /// </remarks>
        public double getValue()
        {
            return data[0];
        }

        /// <summary>
        /// Get a partial derivative.
        /// </summary>
        /// <param name="orders">derivation orders with respect to each variable (if all 
        /// orders are 0, the value is returned)</param>
        /// <returns>partial derivative</returns>
        /// <exception cref="DimensionMismatchException">if the numbers of variables does not
        /// match the instance</exception>
        /// <exception cref="NumberIsTooLargeException"> if sum of derivation orders is larger
        /// than the instance limits</exception>
        /// <remarks>
        /// See <see cref="getValue()"/>
        /// </remarks>
        public double getPartialDerivative(params int[] orders)
        {
            return data[compiler.getPartialDerivativeIndex(orders)];
        }

        /// <summary>
        /// Get all partial derivatives.
        /// </summary>
        /// <returns>a fresh copy of partial derivatives, in an array sorted according to
        /// <see cref="DSCompiler#getPartialDerivativeIndex(params int)"/></returns>
        public double[] getAllDerivatives()
        {
            return (Double[])data.Clone();
        }

        /// <inheritdoc/>
        public DerivativeStructure add(double a)
        {
            DerivativeStructure ds = new DerivativeStructure(this);
            ds.data[0] += a;
            return ds;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure add(DerivativeStructure a)
        {
            compiler.checkCompatibility(a.compiler);
            DerivativeStructure ds = new DerivativeStructure(this);
            compiler.add(data, 0, a.data, 0, ds.data, 0);
            return ds;
        }

        /// <inheritdoc/>
        public DerivativeStructure subtract(double a)
        {
            return add(-a);
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure subtract(DerivativeStructure a)
        {
            compiler.checkCompatibility(a.compiler);
            DerivativeStructure ds = new DerivativeStructure(this);
            compiler.subtract(data, 0, a.data, 0, ds.data, 0);
            return ds;
        }

        /// <inheritdoc/>
        public DerivativeStructure multiply(int n)
        {
            return multiply((double)n);
        }

        /// <inheritdoc/>
        public DerivativeStructure multiply(double a)
        {
            DerivativeStructure ds = new DerivativeStructure(this);
            for (int i = 0; i < ds.data.Length; ++i)
            {
                ds.data[i] *= a;
            }
            return ds;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure multiply(DerivativeStructure a)
        {
            compiler.checkCompatibility(a.compiler);
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.multiply(data, 0, a.data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure divide(double a)
        {
            DerivativeStructure ds = new DerivativeStructure(this);
            for (int i = 0; i < ds.data.Length; ++i)
            {
                ds.data[i] /= a;
            }
            return ds;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure divide(DerivativeStructure a)
        {
            compiler.checkCompatibility(a.compiler);
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.divide(data, 0, a.data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure remainder(double a)
        {
            DerivativeStructure ds = new DerivativeStructure(this);
            ds.data[0] = FastMath.IEEEremainder(ds.data[0], a);
            return ds;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure remainder(DerivativeStructure a)
        {
            compiler.checkCompatibility(a.compiler);
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.remainder(data, 0, a.data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure negate()
        {
            DerivativeStructure ds = new DerivativeStructure(compiler);
            for (int i = 0; i < ds.data.Length; ++i)
            {
                ds.data[i] = -data[i];
            }
            return ds;
        }

        /// <inheritdoc/>
        public DerivativeStructure abs()
        {
            if (BitConverter.DoubleToInt64Bits(data[0]) < 0)
            {
                // we use the bits representation to also handle -0.0
                return negate();
            }
            else
            {
                return this;
            }
        }

        /// <inheritdoc/>
        public DerivativeStructure ceil()
        {
            return new DerivativeStructure(compiler.getFreeParameters(),
                                           compiler.getOrder(),
                                           FastMath.ceil(data[0]));
        }

        /// <inheritdoc/>
        public DerivativeStructure floor()
        {
            return new DerivativeStructure(compiler.getFreeParameters(),
                                           compiler.getOrder(),
                                           FastMath.floor(data[0]));
        }

        /// <inheritdoc/>
        public DerivativeStructure rint()
        {
            return new DerivativeStructure(compiler.getFreeParameters(),
                                           compiler.getOrder(),
                                           FastMath.rint(data[0]));
        }

        /// <inheritdoc/>
        public long round()
        {
            return FastMath.round(data[0]);
        }

        /// <inheritdoc/>
        public DerivativeStructure signum()
        {
            return new DerivativeStructure(compiler.getFreeParameters(),
                                           compiler.getOrder(),
                                           FastMath.signum(data[0]));
        }

        /// <inheritdoc/>
        public DerivativeStructure copySign(DerivativeStructure sign)
        {
            long m = BitConverter.DoubleToInt64Bits(data[0]);
            long s = BitConverter.DoubleToInt64Bits(sign.data[0]);
            if ((m >= 0 && s >= 0) || (m < 0 && s < 0))
            { // Sign is currently OK
                return this;
            }
            return negate(); // flip sign
        }

        /// <inheritdoc/>
        public DerivativeStructure copySign(double sign)
        {
            long m = BitConverter.DoubleToInt64Bits(data[0]);
            long s = BitConverter.DoubleToInt64Bits(sign);
            if ((m >= 0 && s >= 0) || (m < 0 && s < 0))
            { // Sign is currently OK
                return this;
            }
            return negate(); // flip sign
        }

        /// <summary>
        /// Return the exponent of the instance value, removing the bias.
        /// <para>
        /// For double numbers of the form 2^x, the unbiased
        /// exponent is exactly x.
        /// </para>
        /// </summary>
        /// <returns>exponent for instance in IEEE754 representation, without bias</returns>
        public int getExponent()
        {
            return FastMath.getExponent(data[0]);
        }

        /// <inheritdoc/>
        public DerivativeStructure scalb(int n)
        {
            DerivativeStructure ds = new DerivativeStructure(compiler);
            for (int i = 0; i < ds.data.Length; ++i)
            {
                ds.data[i] = FastMath.scalb(data[i], n);
            }
            return ds;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException">if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure hypot(DerivativeStructure y)
        {
            compiler.checkCompatibility(y.compiler);

            if (Double.IsInfinity(data[0]) || Double.IsInfinity(y.data[0]))
            {
                return new DerivativeStructure(compiler.getFreeParameters(),
                                               compiler.getFreeParameters(),
                                               Double.PositiveInfinity);
            }
            else if (Double.IsNaN(data[0]) || Double.IsNaN(y.data[0]))
            {
                return new DerivativeStructure(compiler.getFreeParameters(),
                                               compiler.getFreeParameters(),
                                               Double.NaN);
            }
            else
            {
                int expX = getExponent();
                int expY = y.getExponent();
                if (expX > expY + 27)
                {
                    // y is neglectible with respect to x
                    return abs();
                }
                else if (expY > expX + 27)
                {
                    // x is neglectible with respect to y
                    return y.abs();
                }
                else
                {

                    // find an intermediate scale to avoid both overflow and underflow
                    int middleExp = (expX + expY) / 2;

                    // scale parameters without losing precision
                    DerivativeStructure scaledX = scalb(-middleExp);
                    DerivativeStructure scaledY = y.scalb(-middleExp);

                    // compute scaled hypotenuse
                    DerivativeStructure scaledH =
                            scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                    // remove scaling
                    return scaledH.scalb(middleExp);

                }

            }
        }

        /// <summary>
        /// Returns the hypotenuse of a triangle with sides <c>x</c> and <c>y</c>
        /// - sqrt(x^2&nbsp;+y^2)<br/>
        /// avoiding intermediate overflow or underflow.
        /// <list type="bullet">
        /// <item> If either argument is infinite, then the result is positive infinity.</item>
        /// <item> else, if either argument is NaN then the result is NaN.</item>
        /// </list>
        /// </summary>
        /// <param name="x">a value</param>
        /// <param name="y">a value</param>
        /// <returns>sqrt(x^2&nbsp;+y^2)</returns>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public static DerivativeStructure hypot(DerivativeStructure x, DerivativeStructure y)
        {
            return x.hypot(y);
        }

        /// <summary>
        /// Compute composition of the instance by a univariate function.
        /// </summary>
        /// <param name="f">array of value and derivatives of the function at
        /// the current point (i.e. [f(<see cref="getValue()"/>),
        /// f'(<see cref="getValue()"/>), f''(<see cref="getValue()"/>)...]).</param>
        /// <returns>f(this)</returns>
        /// <exception cref="DimensionMismatchException"> if the number of derivatives
        /// in the array is not equal to <see cref="getOrder() order"/> + 1</exception>
        public DerivativeStructure compose(params double[] f)
        {
            if (f.Length != getOrder() + 1)
            {
                throw new DimensionMismatchException(f.Length, getOrder() + 1);
            }
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.compose(data, 0, f, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure reciprocal()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.pow(data, 0, -1, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure sqrt()
        {
            return rootN(2);
        }

        /// <inheritdoc/>
        public DerivativeStructure cbrt()
        {
            return rootN(3);
        }

        /// <inheritdoc/>
        public DerivativeStructure rootN(int n)
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.rootN(data, 0, n, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public Field<DerivativeStructure> getField()
        {
            return new DerivateStructureField(compiler);
        }

        private class DerivateStructureField : Field<DerivativeStructure>
        {
            private DSCompiler compiler;

            public DerivateStructureField(DSCompiler compiler)
            {
                this.compiler = compiler;
            }
            /// <inheritdoc/>
            public DerivativeStructure getZero()
            {
                return new DerivativeStructure(compiler.getFreeParameters(), compiler.getOrder(), 0.0);
            }

            /// <inheritdoc/>
            public DerivativeStructure getOne()
            {
                return new DerivativeStructure(compiler.getFreeParameters(), compiler.getOrder(), 1.0);
            }

            /// <inheritdoc/>
            public Type getRuntimeClass()
            {
                return this.GetType();
            }
        }

        /// <summary>
        /// Compute a^x where a is a double and x a <see cref="DerivativeStructure"/>
        /// </summary>
        /// <param name="a">number to exponentiate</param>
        /// <param name="x">power to apply</param>
        /// <returns>a^x<returns>
        public static DerivativeStructure pow(double a, DerivativeStructure x)
        {
            DerivativeStructure result = new DerivativeStructure(x.compiler);
            x.compiler.pow(a, x.data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure pow(double p)
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.pow(data, 0, p, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure pow(int n)
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.pow(data, 0, n, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException">if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure pow(DerivativeStructure e)
        {
            compiler.checkCompatibility(e.compiler);
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.pow(data, 0, e.data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure exp()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.exp(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure expm1()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.expm1(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure log()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.log(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure log1p()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.log1p(data, 0, result.data, 0);
            return result;
        }

        /// <summary>
        /// Base 10 logarithm.
        /// </summary>
        /// <returns>base 10 logarithm of the instance</returns>
        public DerivativeStructure log10()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.log10(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure cos()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.cos(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure sin()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.sin(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure tan()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.tan(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure acos()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.acos(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure asin()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.asin(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure atan()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.atan(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure atan2(DerivativeStructure x)
        {
            compiler.checkCompatibility(x.compiler);
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.atan2(data, 0, x.data, 0, result.data, 0);
            return result;
        }

        /// <summary>
        /// Two arguments arc tangent operation.
        /// </summary>
        /// <param name="y">first argument of the arc tangent</param>
        /// <param name="x">second argument of the arc tangent</param>
        /// <returns>atan2(y, x)</returns>
        /// <exception cref="DimensionMismatchException">if number of free parameters
        /// orders do not match</exception>
        public static DerivativeStructure atan2(DerivativeStructure y, DerivativeStructure x)
        {
            return y.atan2(x);
        }

        /// <inheritdoc/>
        public DerivativeStructure cosh()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.cosh(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure sinh()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.sinh(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure tanh()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.tanh(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure acosh()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.acosh(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure asinh()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.asinh(data, 0, result.data, 0);
            return result;
        }

        /// <inheritdoc/>
        public DerivativeStructure atanh()
        {
            DerivativeStructure result = new DerivativeStructure(compiler);
            compiler.atanh(data, 0, result.data, 0);
            return result;
        }

        /// <summary>
        /// Convert radians to degrees, with error of less than 0.5 ULP
        /// </summary>
        /// <returns>instance converted into degrees</returns>
        public DerivativeStructure toDegrees()
        {
            DerivativeStructure ds = new DerivativeStructure(compiler);
            for (int i = 0; i < ds.data.Length; ++i)
            {
                ds.data[i] = FastMath.toDegrees(data[i]);
            }
            return ds;
        }

        /// <summary>
        /// Convert degrees to radians, with error of less than 0.5 ULP
        /// </summary>
        /// <returns>instance converted into radians</returns>
        public DerivativeStructure toRadians()
        {
            DerivativeStructure ds = new DerivativeStructure(compiler);
            for (int i = 0; i < ds.data.Length; ++i)
            {
                ds.data[i] = FastMath.toRadians(data[i]);
            }
            return ds;
        }

        /// <summary>
        /// Evaluate Taylor expansion a derivative structure.
        /// </summary>
        /// <param name="delta">parameters offsets (&Delta;x, &Delta;y, ...)</param>
        /// <returns>value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...</returns>
        /// <exception cref="MathArithmeticException"> if factorials becomes too large</exception>
        public double taylor(params double[] delta)
        {
            return compiler.taylor(data, 0, delta);
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(DerivativeStructure[] a, DerivativeStructure[] b)
        {

            // compute an accurate value, taking care of cancellations
            double[] aDouble = new double[a.Length];
            for (int i = 0; i < a.Length; ++i)
            {
                aDouble[i] = a[i].getValue();
            }
            double[] bDouble = new double[b.Length];
            for (int i = 0; i < b.Length; ++i)
            {
                bDouble[i] = b[i].getValue();
            }
            double accurateValue = MathArrays.linearCombination(aDouble, bDouble);

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = a[0].getField().getZero();
            for (int i = 0; i < a.Length; ++i)
            {
                simpleValue = simpleValue.add(a[i].multiply(b[i]));
            }

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(simpleValue.getFreeParameters(), simpleValue.getOrder(), all);

        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(double[] a, DerivativeStructure[] b)
        {

            // compute an accurate value, taking care of cancellations
            double[] bDouble = new double[b.Length];
            for (int i = 0; i < b.Length; ++i)
            {
                bDouble[i] = b[i].getValue();
            }
            double accurateValue = MathArrays.linearCombination(a, bDouble);

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = b[0].getField().getZero();
            for (int i = 0; i < a.Length; ++i)
            {
                simpleValue = simpleValue.add(b[i].multiply(a[i]));
            }

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(simpleValue.getFreeParameters(), simpleValue.getOrder(), all);

        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(DerivativeStructure a1, DerivativeStructure b1, DerivativeStructure a2, DerivativeStructure b2)
        {

            // compute an accurate value, taking care of cancellations
            double accurateValue = MathArrays.linearCombination(a1.getValue(), b1.getValue(),
                                                                      a2.getValue(), b2.getValue());

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = a1.multiply(b1).add(a2.multiply(b2));

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(getFreeParameters(), getOrder(), all);

        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(double a1, DerivativeStructure b1, double a2, DerivativeStructure b2)
        {
            // compute an accurate value, taking care of cancellations
            double accurateValue = MathArrays.linearCombination(a1, b1.getValue(),
                                                                      a2, b2.getValue());

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = b1.multiply(a1).add(b2.multiply(a2));

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(getFreeParameters(), getOrder(), all);
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(DerivativeStructure a1, DerivativeStructure b1, DerivativeStructure a2, DerivativeStructure b2, DerivativeStructure a3, DerivativeStructure b3)
        {
            // compute an accurate value, taking care of cancellations
            double accurateValue = MathArrays.linearCombination(a1.getValue(), b1.getValue(),
                                                                      a2.getValue(), b2.getValue(),
                                                                      a3.getValue(), b3.getValue());

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3));

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(getFreeParameters(), getOrder(), all);
        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(double a1, DerivativeStructure b1, double a2, DerivativeStructure b2, double a3, DerivativeStructure b3)
        {

            // compute an accurate value, taking care of cancellations
            double accurateValue = MathArrays.linearCombination(a1, b1.getValue(),
                                                                      a2, b2.getValue(),
                                                                      a3, b3.getValue());

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = b1.multiply(a1).add(b2.multiply(a2)).add(b3.multiply(a3));

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(getFreeParameters(), getOrder(), all);

        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(DerivativeStructure a1, DerivativeStructure b1, DerivativeStructure a2, DerivativeStructure b2, DerivativeStructure a3, DerivativeStructure b3, DerivativeStructure a4, DerivativeStructure b4)
        {
            // compute an accurate value, taking care of cancellations
            double accurateValue = MathArrays.linearCombination(a1.getValue(), b1.getValue(),
                                                                      a2.getValue(), b2.getValue(),
                                                                      a3.getValue(), b3.getValue(),
                                                                      a4.getValue(), b4.getValue());

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3)).add(a4.multiply(b4));

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(getFreeParameters(), getOrder(), all);

        }

        /// <inheritdoc/>
        /// <exception cref="DimensionMismatchException"> if number of free parameters
        /// or orders do not match</exception>
        public DerivativeStructure linearCombination(double a1, DerivativeStructure b1, double a2, DerivativeStructure b2, double a3, DerivativeStructure b3, double a4, DerivativeStructure b4)
        {

            // compute an accurate value, taking care of cancellations
            double accurateValue = MathArrays.linearCombination(a1, b1.getValue(),
                                                                      a2, b2.getValue(),
                                                                      a3, b3.getValue(),
                                                                      a4, b4.getValue());

            // compute a simple value, with all partial derivatives
            DerivativeStructure simpleValue = b1.multiply(a1).add(b2.multiply(a2)).add(b3.multiply(a3)).add(b4.multiply(a4));

            // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
            double[] all = simpleValue.getAllDerivatives();
            all[0] = accurateValue;
            return new DerivativeStructure(getFreeParameters(), getOrder(), all);
        }

        /// <summary>
        /// Test for the equality of two derivative structures.
        /// <para>
        /// Derivative structures are considered equal if they have the same number
        /// of free parameters, the same derivation order, and the same derivatives.
        /// </para>
        /// </summary>
        /// <param name="other">Object to test for equality to this</param>
        /// <returns>true if two derivative structures are equal</returns>
        public override Boolean Equals(Object other)
        {
            if (this == other)
            {
                return true;
            }

            if (other is DerivativeStructure)
            {
                DerivativeStructure rhs = (DerivativeStructure)other;
                return (getFreeParameters() == rhs.getFreeParameters()) &&
                       (getOrder() == rhs.getOrder()) &&
                       MathArrays.equals(data, rhs.data);
            }
            return false;
        }

        /// <summary>
        /// Get a hashCode for the derivative structure. 
        /// </summary>
        /// <returns>a hash code value for this object</returns>
        public override int GetHashCode()
        {
            return 227 + 229 * getFreeParameters() + 233 * getOrder() + 239 * MathUtils.hash(data);
        }
    }
}