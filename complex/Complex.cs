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
using System.Collections.Generic;

namespace Math3.complex
{
    /// <summary>
    /// Representation of a Complex number, i.e. a number which has both a
    /// real and imaginary part.
    /// <para/>
    /// Implementations of arithmetic operations handle <c>NaN</c> and
    /// infinite values according to the rules for <c>System.Double</c>, i.e.
    /// equals is an equivalence relation for all instances that have
    /// a <c>NaN</c> in either real or imaginary part, e.g. the following are
    /// considered equal:
    /// <list type="bullet">
    /// <item><c>1 + NaNi</c></item>
    /// <item><c>NaN + i</c></item>
    /// <item><c>NaN + NaNi</c></item>
    /// </list>
    /// Note that this is in contradiction with the IEEE-754 standard for floating
    /// point numbers (according to which the test <c>x == x</c> must fail if
    /// <c>x</c> is <c>NaN</c>). The method
    /// <see cref="Math3.util.Precision.equals(double,double,int)"/>
    /// equals for primitive double in <see cref="Math3.util.Precision"/>
    /// conforms with IEEE-754 while this class conforms with the standard behavior
    /// for C# object types.
    /// </summary>
    /// <remarks>I arbitrarily decided to add the operators +, -, *, / 
    /// which use the corrispondent method. This maybe was wanted
    /// by the original developers too, but it's impossible in Java,
    /// as far as I know.</remarks>
    public class Complex : FieldElement<Complex>
    {
        /// <summary>
        /// The square root of -1. A number representing "0.0 + 1.0i"
        /// </summary>
        public static readonly Complex I = new Complex(0.0, 1.0);

        // CHECKSTYLE: stop ConstantName
        /// <summary>
        /// A complex number representing "NaN + NaNi"
        /// </summary>
        public static readonly Complex NaN = new Complex(Double.NaN, Double.NaN);

        // CHECKSTYLE: resume ConstantName
        /// <summary>
        /// A complex number representing "+INF + INFi"
        /// </summary>
        public static readonly Complex INF = new Complex(Double.PositiveInfinity, Double.PositiveInfinity);

        /// <summary>
        /// A complex number representing "1.0 + 0.0i"
        /// </summary>
        public static readonly Complex ONE = new Complex(1.0, 0.0);

        /// <summary>
        /// A complex number representing "0.0 + 0.0i"
        /// </summary>
        public static readonly Complex ZERO = new Complex(0.0, 0.0);

        /// <summary>
        /// The imaginary part.
        /// </summary>
        private readonly double imaginary;

        /// <summary>
        /// The real part.
        /// </summary>
        private readonly double real;

        /// <summary>
        /// Record whether this complex number is equal to NaN.
        /// </summary>
        private readonly Boolean isNaN;

        /// <summary>
        /// Record whether this complex number is infinite.
        /// </summary>
        private readonly Boolean isInfinite;

        /// <summary>
        /// Create a complex number given only the real part.
        /// </summary>
        /// <param name="real">Real part.</param>
        public Complex(double real) : this(real, 0.0) { }

        /// <summary>
        /// Create a complex number given only the real part.
        /// </summary>
        /// <param name="real">Real part.</param>
        /// <param name="imaginary">Imaginary part.</param>
        public Complex(double real, double imaginary)
        {
            this.real = real;
            this.imaginary = imaginary;

            isNaN = Double.IsNaN(real) || Double.IsNaN(imaginary);
            isInfinite = !isNaN &&
                (Double.IsInfinity(real) || Double.IsInfinity(imaginary));
        }

        /// <summary>
        /// Return the absolute value of this complex number.
        /// Returns <c>NaN</c> if either real or imaginary part is <c>NaN</c>
        /// and <c>Double.PositiveInfinity</c> if neither part is <c>NaN</c>,
        /// but at least one part is infinite.
        /// </summary>
        /// <returns>the absolute value.</returns>
        public double abs()
        {
            if (isNaN)
            {
                return Double.NaN;
            }
            if (IsInfinity())
            {
                return Double.PositiveInfinity;
            }
            if (FastMath.abs(real) < FastMath.abs(imaginary))
            {
                if (imaginary == 0.0)
                {
                    return FastMath.abs(real);
                }
                double q = real / imaginary;
                return FastMath.abs(imaginary) * FastMath.sqrt(1 + q * q);
            }
            else
            {
                if (real == 0.0)
                {
                    return FastMath.abs(imaginary);
                }
                double q = imaginary / real;
                return FastMath.abs(real) * FastMath.sqrt(1 + q * q);
            }
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is
        /// <c>(this + addend)</c>.
        /// Uses the definitional formula
        /// <code>
        /// (a + bi) + (c + di) = (a+c) + (b+d)i
        /// </code>
        /// <para/>
        /// If either <c>this</c> or <c>addend</c> has a <c>NaN</c> value in
        /// either part, <see cref="NaN"/> is returned; otherwise <c>Infinite</c>
        /// and <c>NaN</c> values are returned in the parts of the result
        /// according to the rules for <see cref="System.Double"/> arithmetic.
        /// </summary>
        /// <param name="addend">Value to be added to this <c>Complex</c>.</param>
        /// <returns><c>this + addend</c>.</returns>
        /// <exception cref="NullArgumentException">if <c>addend</c> is <c>null</c>.</exception>
        public Complex add(Complex addend)
        {
            MathUtils.checkNotNull(addend);
            if (isNaN || addend.isNaN)
            {
                return NaN;
            }

            return createComplex(real + addend.getReal(),
                                 imaginary + addend.getImaginary());
        }

        public static Complex operator +(Complex Complex1, Complex Complex2)
        {
            return Complex1.add(Complex2);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is <c>(this + addend)</c>,
        /// with <c>addend</c> interpreted as a real number.
        /// </summary>
        /// <param name="addend">Value to be added to this <c>Complex</c>.</param>
        /// <returns><c>this + addend</c>.</returns>
        /// <remarks>See <see cref="add(Complex)"/></remarks>
        public Complex add(double addend)
        {
            if (isNaN || Double.IsNaN(addend))
            {
                return NaN;
            }

            return createComplex(real + addend, imaginary);
        }

        public static Complex operator +(Complex Complex, Double Double)
        {
            return Complex.add(Double);
        }

        public static Complex operator +(Double Double, Complex Complex)
        {
            return Complex.add(Double);
        }

        /// <summary>
        /// Return the conjugate of this complex number.
        /// The conjugate of <c>a + bi</c> is <c>a - bi</c>.
        /// <para/>
        /// <see cref="NaN"/> is returned if either the real or imaginary
        /// part of this Complex number equals <c>Double.NaN</c>.
        /// <para/>
        /// If the imaginary part is infinite, and the real part is not
        /// <c>NaN</c>, the returned value has infinite imaginary part
        /// of the opposite sign, e.g. the conjugate of
        /// <c>1 + POSITIVE_INFINITY i</c> is <c>1 - NEGATIVE_INFINITY i</c>.
        /// </summary>
        /// <returns>the conjugate of this Complex object.</returns>
        public Complex conjugate()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(real, -imaginary);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is
        /// <c>(this / divisor)</c>.
        /// Implements the definitional formula
        /// <code>
        ///    a + bi          ac + bd + (bc - ad)i
        ///    ----------- = -------------------------
        ///    c + di              c^2 + d^2
        /// </code>
        /// but uses
        /// <a href="http://doi.acm.org/10.1145/1039813.1039814">
        /// prescaling of operands</a> to limit the effects of overflows and
        /// underflows in the computation.
        /// <para/>
        /// <c>Infinite</c> and <c>NaN</c> values are handled according to the
        /// following rules, applied in the order presented:
        /// <list type="bullet">
        /// <item>If either <c>this</c> or <c>divisor</c> has a <c>NaN</c> value
        /// in either part, <see cref="NaN"/> is returned.
        /// </item>
        /// <item>If <c>divisor</c> equals <see cref="ZERO"/>, <see cref="NaN"/> is returned.
        /// </item>
        /// <item>If <c>this</c> and <c>divisor</c> are both infinite,
        /// <see cref="NaN"/> is returned.
        /// </item>
        /// <item>If <c>this</c> is finite (i.e., has no <c>Infinite</c> or
        /// <c>NaN</c> parts) and <c>divisor</c> is infinite (one or both parts
        /// infinite), <see cref="ZERO"/> is returned.
        /// </item>
        /// <item>If <c>this</c> is infinite and <c>divisor</c> is finite,
        /// <c>NaN</c> values are returned in the parts of the result if the
        /// <see cref="Sytem.Double"/> rules applied to the definitional formula
        /// force <c>NaN</c> results.
        /// </item>
        /// </list>
        /// </summary>
        /// <param name="divisor">Value by which this <c>Complex</c> is to be divided.</param>
        /// <returns><c>this / divisor</c>.</returns>
        /// <exception cref="NullArgumentException">if <c>divisor</c> is <c>null</c>.</exception>
        public Complex divide(Complex divisor)
        {
            MathUtils.checkNotNull(divisor);
            if (isNaN || divisor.isNaN)
            {
                return NaN;
            }

            double c = divisor.getReal();
            double d = divisor.getImaginary();
            if (c == 0.0 && d == 0.0)
            {
                return NaN;
            }

            if (divisor.IsInfinity() && !IsInfinity())
            {
                return ZERO;
            }

            if (FastMath.abs(c) < FastMath.abs(d))
            {
                double q = c / d;
                double denominator = c * q + d;
                return createComplex((real * q + imaginary) / denominator,
                    (imaginary * q - real) / denominator);
            }
            else
            {
                double q = d / c;
                double denominator = d * q + c;
                return createComplex((imaginary * q + real) / denominator,
                    (imaginary - real * q) / denominator);
            }
        }

        public static Complex operator /(Complex Complex1, Complex Complex2)
        {
            return Complex1.divide(Complex2);
        }
        

        /// <summary>
        /// Returns a <c>Complex</c> whose value is <c>(this / divisor)</c>,
        /// with <c>divisor</c> interpreted as a real number.
        /// </summary>
        /// <param name="divisor">Value by which this <c>Complex</c> is to be divided.</param>
        /// <returns><c>this / divisor</c>.</returns>
        /// <remarks><see cref="divide(Complex)"/></remarks>
        public Complex divide(double divisor)
        {
            if (isNaN || Double.IsNaN(divisor))
            {
                return NaN;
            }
            if (divisor == 0d)
            {
                return NaN;
            }
            if (Double.IsInfinity(divisor))
            {
                return !IsInfinity() ? ZERO : NaN;
            }
            return createComplex(real / divisor,
                                 imaginary / divisor);
        }

        public static Complex operator /(Complex Complex, Double Double)
        {
            return Complex.divide(Double);
        }

        public static Complex operator /(Double Double, Complex Complex)
        {
            return Complex.divide(Double);
        }

        /// <inheritdoc/>
        public Complex reciprocal()
        {
            if (isNaN)
            {
                return NaN;
            }

            if (real == 0.0 && imaginary == 0.0)
            {
                return INF;
            }

            if (isInfinite)
            {
                return ZERO;
            }

            if (FastMath.abs(real) < FastMath.abs(imaginary))
            {
                double q = real / imaginary;
                double scale = 1d / (real * q + imaginary);
                return createComplex(scale * q, -scale);
            }
            else
            {
                double q = imaginary / real;
                double scale = 1d / (imaginary * q + real);
                return createComplex(scale, -scale * q);
            }
        }

        /// <summary>
        /// Test for the floating-point equality between Complex objects.
        /// It returns <c>true</c> if both arguments are equal or within the
        /// range of allowed error (inclusive).
        /// </summary>
        /// <param name="x">First value (cannot be <c>null</c>).</param>
        /// <param name="y">Second value (cannot be <c>null</c>).</param>
        /// <param name="maxUlps">maxUlps <c>(maxUlps - 1)</c> is the number of floating point
        /// values between the real (resp. imaginary) parts of <c>x</c> and
        /// <c>y</c>.</param>
        /// <returns><c>true</c> if there are fewer than <c>maxUlps</c> floating
        /// point values between the real (resp. imaginary) parts of <c>x</c>
        /// and <c>y</c>.</returns>
        /// <remarks><see cref="Math3.util.Precision.equals(double,double,int)"/></remarks>
        public static Boolean equals(Complex x, Complex y, int maxUlps)
        {
            return Precision.equals(x.real, y.real, maxUlps) &&
                Precision.equals(x.imaginary, y.imaginary, maxUlps);
        }

        /// <summary>
        /// Returns <c>true</c> iff the values are equal as defined by
        /// <see cref="equals(Complex,Complex,int)"/>.
        /// </summary>
        /// <param name="x">First value (cannot be <c>null</c>).</param>
        /// <param name="y">Second value (cannot be <c>null</c>).</param>
        /// <returns> <c>true</c> if the values are equal.</returns>
        public override Boolean Equals(Object ToBeCompared)
        {
            if (ToBeCompared is Complex)
            {
                return Complex.equals(this, (Complex)ToBeCompared, 1);
            }
            return false;
        }

        public static Boolean operator ==(Complex Complex1, Complex Complex2)
        {
            return Complex.equals(Complex1, Complex2, 1);
        }

        public static Boolean operator !=(Complex Complex1, Complex Complex2)
        {
            return !Complex.equals(Complex1, Complex2, 1);
        }

        /// <summary>
        /// Returns <c>true</c> if, both for the real part and for the imaginary
        /// part, there is no double value strictly between the arguments or the
        /// difference between them is within the range of allowed error
        /// (inclusive).
        /// </summary>
        /// <param name="x">First value (cannot be <c>null</c>).</param>
        /// <param name="y">Second value (cannot be <c>null</c>).</param>
        /// <param name="eps">Amount of allowed absolute error.</param>
        /// <returns><c>true</c> if the values are two adjacent floating point
        /// numbers or they are within range of each other.</returns>
        /// <remarks><see cref="Math3.util.Precision.equals(double,double,double)"/></remarks>
        public static Boolean equals(Complex x, Complex y, double eps)
        {
            return Precision.equals(x.real, y.real, eps) &&
                Precision.equals(x.imaginary, y.imaginary, eps);
        }

        /// <summary>
        /// Returns <c>true</c> if, both for the real part and for the imaginary
        /// part, there is no double value strictly between the arguments or the
        /// relative difference between them is smaller or equal to the given
        /// tolerance. 
        /// </summary>
        /// <param name="x">First value (cannot be <c>null</c>).</param>
        /// <param name="y">Second value (cannot be <c>null</c>).</param>
        /// <param name="eps">Amount of allowed relative error.</param>
        /// <returns<c>true</c> if the values are two adjacent floating point
        /// numbers or they are within range of each other.></returns>
        /// <remarks><see cref="Precision.equalsWithRelativeTolerance(double,double,double)"/></remarks>
        public static Boolean equalsWithRelativeTolerance(Complex x, Complex y,
                                                          double eps)
        {
            return Precision.equalsWithRelativeTolerance(x.real, y.real, eps) &&
                Precision.equalsWithRelativeTolerance(x.imaginary, y.imaginary, eps);
        }

        /// <summary>
        /// Get a hashCode for the complex number.
        /// Any <c>Double.NaN</c> value in real or imaginary part produces
        /// the same hash code <c>7</c>.
        /// </summary>
        /// <returns>a hash code value for this object.</returns>
        public override int GetHashCode()
        {
            if (isNaN)
            {
                return 7;
            }
            return 37 * (17 * MathUtils.hash(imaginary) +
                MathUtils.hash(real));
        }

        /// <summary>
        /// Access the imaginary part.
        /// </summary>
        /// <returns>the imaginary part.</returns>
        public double getImaginary()
        {
            return imaginary;
        }

        /// <summary>
        /// Access the real part. 
        /// </summary>
        /// <returns>the real part.</returns>
        public double getReal()
        {
            return real;
        }

        /// <summary>
        /// Checks whether either or both parts of this complex number is
        /// <c>NaN</c>.
        /// </summary>
        /// <returns>true if either or both parts of this complex number is
        /// <c>NaN</c>; false otherwise.</returns>
        public Boolean IsNaN()
        {
            return isNaN;
        }

        /// <summary>
        /// Checks whether either the real or imaginary part of this complex number
        /// takes an infinite value (either <c>Double.PositiveInfinity</c> or
        /// <c>Double.NegativeInfinity</c>) and neither part
        /// is <c>NaN</c>.
        /// </summary>
        /// <returns>true if one or both parts of this complex number are infinite
        /// and neither part is <c>NaN</c>.</returns>
        public Boolean IsInfinity()
        {
            return isInfinite;
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is <c>this * factor</c>.
        /// Implements preliminary checks for <c>NaN</c> and infinity followed by
        /// the definitional formula:
        /// <code>
        ///   (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        /// </code>
        /// Returns <see cref="#NaN"/> if either <c>this</c> or <c>factor</c> has one or
        /// more <c>NaN</c> parts.
        /// <para/>
        /// Returns <see cref="INF"/> if neither <c>this</c> nor <c>factor</c> has one
        /// or more <c>NaN</c> parts and if either <c>this</c> or <c>factor</c>
        /// has one or more infinite parts (same result is returned regardless of
        /// the sign of the components).
        /// <para/>
        /// Returns finite values in components of the result per the definitional
        /// formula in all remaining cases.
        /// </summary>
        /// <param name="factor">value to be multiplied by this <c>Complex</c>.</param>
        /// <returns><c>this * factor</c>.</returns>
        /// <exception cref="NullArgumentException">if <c>factor</c> is <c>null</c>.</exception>
        public Complex multiply(Complex factor)
        {
            MathUtils.checkNotNull(factor);
            if (isNaN || factor.isNaN)
            {
                return NaN;
            }
            if (Double.IsInfinity(real) ||
                Double.IsInfinity(imaginary) ||
                Double.IsInfinity(factor.real) ||
                Double.IsInfinity(factor.imaginary))
            {
                // we don't use isInfinite() to avoid testing for NaN again
                return INF;
            }
            return createComplex(real * factor.real - imaginary * factor.imaginary,
                                 real * factor.imaginary + imaginary * factor.real);
        }

        public static Complex operator *(Complex Complex1, Complex Complex2)
        {
            return Complex1.multiply(Complex2);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is <c>this * factor</c>, with <c>factor</c>
        /// interpreted as a integer number.
        /// </summary>
        /// <param name="factor">value to be multiplied by this <c>Complex</c>.</param>
        /// <returns><c>this * factor</c>.</returns>
        /// <see cref="multiply(Complex)"/>
        public Complex multiply(int factor)
        {
            if (isNaN)
            {
                return NaN;
            }
            if (Double.IsInfinity(real) ||
                Double.IsInfinity(imaginary))
            {
                return INF;
            }
            return createComplex(real * factor, imaginary * factor);
        }

        public static Complex operator *(Complex Complex, Int32 Integer)
        {
            return Complex.multiply(Integer);
        }

        public static Complex operator *(Int32 Integer, Complex Complex)
        {
            return Complex.multiply(Integer);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is <c>this * factor</c>, with <c>factor</c>
        /// interpreted as a real number.
        /// </summary>
        /// <param name="factor">value to be multiplied by this <c>Complex</c>.</param>
        /// <returns><c>this * factor</c>.</returns>
        /// <see cref="multiply(Complex)"/>
        public Complex multiply(double factor)
        {
            if (isNaN || Double.IsNaN(factor))
            {
                return NaN;
            }
            if (Double.IsInfinity(real) ||
                Double.IsInfinity(imaginary) ||
                Double.IsInfinity(factor))
            {
                // we don't use isInfinite() to avoid testing for NaN again
                return INF;
            }
            return createComplex(real * factor, imaginary * factor);
        }

        public static Complex operator *(Complex Complex, Double Double)
        {
            return Complex.multiply(Double);
        }

        public static Complex operator *(Double Double, Complex Complex)
        {
            return Complex.multiply(Double);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is <c>(-this)</c>.
        /// Returns <c>NaN</c> if either real or imaginary
        /// part of this Complex number equals <c>Double.NaN</c>.
        /// </summary>
        /// <returns><c>-this</c>.</returns>
        public Complex negate()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(-real, -imaginary);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is
        /// <c>(this - subtrahend)</c>.
        /// Uses the definitional formula
        /// <code>
        ///     (a + bi) - (c + di) = (a-c) + (b-d)i
        /// </code>
        /// If either <c>this</c> or <c>subtrahend</c> has a <c>NaN]</c> value in either part,
        /// <see cref="NaN"/> is returned; otherwise infinite and <c>NaN</c> values are
        /// returned in the parts of the result according to the rules for
        /// <see cref="System.Double"/> arithmetic.
        /// </summary>
        /// <param name="subtrahend">value to be subtracted from this <c>Complex</c>.</param>
        /// <returns><c>this - subtrahend</c>.</returns>
        /// <exception cref="NullArgumentException">if <c>subtrahend</c> is <c>null</c>.</exception>
        public Complex subtract(Complex subtrahend)
        {
            MathUtils.checkNotNull(subtrahend);
            if (isNaN || subtrahend.isNaN)
            {
                return NaN;
            }

            return createComplex(real - subtrahend.getReal(),
                                 imaginary - subtrahend.getImaginary());
        }

        public static Complex operator -(Complex Complex1, Complex Complex2)
        {
            return Complex1.subtract(Complex2);
        }

        /// <summary>
        /// Returns a <c>Complex</c> whose value is
        /// <c>(this - subtrahend)</c>.
        /// </summary>
        /// <param name="subtrahend">value to be subtracted from this <c>Complex</c>.</param>
        /// <returns><c>this - subtrahend</c>.</returns>
        /// <remarks><see cref="subtract(Complex)"/></remarks>
        public Complex subtract(double subtrahend)
        {
            if (isNaN || Double.IsNaN(subtrahend))
            {
                return NaN;
            }
            return createComplex(real - subtrahend, imaginary);
        }

        public static Complex operator -(Complex Complex, Double Double)
        {
            return Complex.multiply(Double);
        }

        public static Complex operator -(Double Double, Complex Complex)
        {
            return Complex.multiply(Double);
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
        /// inverse cosine</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///   acos(z) = -i (log(z + i (sqrt(1 - z^2))))
        /// </code>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c> or infinite.
        /// </summary>
        /// <returns>the inverse cosine of this complex number.</returns>
        public Complex acos()
        {
            if (isNaN)
            {
                return NaN;
            }

            return this.add(this.sqrt1z().multiply(I)).log().multiply(I.negate());
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/InverseSine.html" TARGET="_top">
        /// inverse sine</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///  asin(z) = -i (log(sqrt(1 - z<sup>2</sup>) + iz))
        /// </code>
        /// Returns <see cref="Complex#NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c> or infinite.
        /// </summary>
        /// <returns>the inverse sine of this complex number.</returns>
        public Complex asin()
        {
            if (isNaN)
            {
                return NaN;
            }

            return sqrt1z().add(this.multiply(I)).log().multiply(I.negate());
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/InverseTangent.html" TARGET="_top">
        /// Implements the formula:
        /// <code>
        ///   atan(z) = (i/2) log((i + z)/(i - z))
        /// </code>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c> or infinite.
        /// </summary>
        /// <returns>the inverse tangent of this complex number</returns>
        public Complex atan()
        {
            if (isNaN)
            {
                return NaN;
            }

            return this.add(I).divide(I.subtract(this)).log()
                    .multiply(I.divide(createComplex(2.0, 0.0)));
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/Cosine.html" TARGET="_top">
        /// cosine</a>
        /// of this complex number.
        /// Implements the formula:
        /// <code>
        ///  cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <see cref="FastMath.sin"/>, <see cref="FastMath.cos"/>,
        /// <see cref="FastMath.cosh"/> and <see cref="FastMath.sinh"/>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result.
        /// </summary>
        /// <example>
        /// <code>
        ///   cos(1 &plusmn; INFINITY i) = 1 &#x2213; INFINITY i
        ///   cos(&plusmn;INFINITY + i) = NaN + NaN i
        ///   cos(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        /// </code>
        /// </example>
        /// <returns>the cosine of this complex number.</returns>
        public Complex cos()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(FastMath.cos(real) * FastMath.cosh(imaginary),
                                 -FastMath.sin(real) * FastMath.sinh(imaginary));
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/HyperbolicCosine.html" TARGET="_top">
        /// hyperbolic cosine</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///   cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i}
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <see cref="FastMath.sin"/>, <see cref="FastMath.cos"/>,
        /// <see cref="FastMath.cosh"/> and <see cref="FastMath.sinh"/>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result.
        /// </summary>
        /// <example>
        /// <code>
        ///   cosh(1 &plusmn; INFINITY i) = NaN + NaN i
        ///   cosh(&plusmn;INFINITY + i) = INFINITY &plusmn; INFINITY i
        ///   cosh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        /// </code>
        /// </example>
        /// <returns>the hyperbolic cosine of this complex number.</returns>
        public Complex cosh()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(FastMath.cosh(real) * FastMath.cos(imaginary),
                                 FastMath.sinh(real) * FastMath.sin(imaginary));
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/ExponentialFunction.html" TARGET="_top">
        /// exponential function</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///   exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <see cref="FastMath.exp"/>, <see cref="FastMath.cos"/>, and
        /// <see cref="FastMath.sin"/>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result.
        /// </summary>
        /// <example>
        /// <code>
        ///   exp(1 &plusmn; INFINITY i) = NaN + NaN i
        ///   exp(INFINITY + i) = INFINITY + INFINITY i
        ///   exp(-INFINITY + i) = 0 + 0i
        ///   exp(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        /// </code>
        /// </example>
        /// <returns><code><i>e</i><sup>this</sup></code>.</returns>
        public Complex exp()
        {
            if (isNaN)
            {
                return NaN;
            }

            double expReal = FastMath.exp(real);
            return createComplex(expReal * FastMath.cos(imaginary),
                                 expReal * FastMath.sin(imaginary));
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/NaturalLogarithm.html" TARGET="_top">
        /// natural logarithm</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///   log(a + bi) = ln(|a + bi|) + arg(a + bi)i
        /// </code>
        /// where ln on the right hand side is <see cref="FastMath.log"/>,
        /// <c>|a + bi|</c> is the modulus, <see cref="abs"/>,  and
        /// <c>arg(a + bi) = </c><see cref="FastMath.atan2(double, double)"/>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite (or critical) values in real or imaginary parts of the input may
        /// result in infinite or NaN values returned in parts of the result.
        /// <example>
        /// <code>
        ///   log(1 &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/2)i
        ///   log(INFINITY + i) = INFINITY + 0i
        ///   log(-INFINITY + i) = INFINITY + &pi;i
        ///   log(INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/4)i
        ///   log(-INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (3&pi;/4)i
        ///   log(0 + 0i) = -INFINITY + 0i
        /// </code>
        /// </example>
        /// </summary>
        /// <returns>the value <code>ln &nbsp; this</code>, the natural logarithm
        /// of <c>this</c>.</returns>
        public Complex log()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(FastMath.log(abs()),
                                 FastMath.atan2(imaginary, real));
        }

        /// <summary>
        /// Returns of value of this complex number raised to the power of <c>x</c>.
        /// Implements the formula:
        /// <code>
        ///   y<sup>x</sup> = exp(x&middot;log(y))
        /// </code>
        /// where <c>exp</c> and <c>log</c> are <see cref="exp"/> and
        /// <see cref="log"/>, respectively.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c> or infinite, or if <c>y</c>
        /// equals <see cref="ZERO"/>.
        /// </summary>
        /// <param name="x">exponent to which this <c>Complex</c> is to be raised.</param>
        /// <returns><code> <c>this^x</c></code>.</returns>
        /// <exception cref="NullArgumentException">if x is <c>null</c>.</exception>
        public Complex pow(Complex x)
        {
            MathUtils.checkNotNull(x);
            return this.log().multiply(x).exp();
        }

        /// <summary>
        /// Returns of value of this complex number raised to the power of <c>x</c>.
        /// </summary>
        /// <param name="x">exponent to which this <c>Complex</c> is to be raised.</param>
        /// <returns><code>this^x</code>.</returns>
        /// <remarks>
        /// See <see cref="pow(Complex)"/>
        /// </remarks>
        public Complex pow(double x)
        {
            return this.log().multiply(x).exp();
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/Sine.html" TARGET="_top">
        /// sine</a>
        /// of this complex number.
        /// Implements the formula:
        /// <code>
        ///   sin(a + bi) = sin(a)cosh(b) - cos(a)sinh(b)i
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <c>FastMath.sin</c>, <c>FastMath.cos</c>,
        /// <c>FastMath.cosh</c> and <c>FastMath.sinh</c>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or <c>NaN</c> values returned in parts of the result.
        /// </summary>
        /// <example>
        /// <code>
        ///   sin(1 &plusmn; INFINITY i) = 1 &plusmn; INFINITY i
        ///   sin(&plusmn;INFINITY + i) = NaN + NaN i
        ///   sin(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        /// </code>
        /// </example>
        /// <returns>the sine of this complex number.</returns>
        public Complex sin()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(FastMath.sin(real) * FastMath.cosh(imaginary),
                                 FastMath.cos(real) * FastMath.sinh(imaginary));
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/HyperbolicSine.html" TARGET="_top">
        /// Implements the formula:
        /// <code>
        ///   sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <see cref="FastMath.sin"/>, <see cref="FastMath.cos"/>,
        /// <see cref="FastMath.cosh"/> and <see cref="FastMath.sinh"/>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result.
        /// <example>
        /// <code>
        ///   sinh(1 &plusmn; INFINITY i) = NaN + NaN i
        ///   sinh(&plusmn;INFINITY + i) = &plusmn; INFINITY + INFINITY i
        ///   sinh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        /// </code>
        /// </example>
        /// </summary>
        /// <returns>the hyperbolic sine of <c>this</c>.</returns>
        public Complex sinh()
        {
            if (isNaN)
            {
                return NaN;
            }

            return createComplex(FastMath.sinh(real) * FastMath.cos(imaginary),
                FastMath.cosh(real) * FastMath.sin(imaginary));
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
        /// square root</a> of this complex number.
        /// Implements the following algorithm to compute <c>sqrt(a + bi)</c>:
        /// <list type="number">
        /// <item>Let <c>t = sqrt((|a| + |a + bi|) / 2)</c></item>
        /// <item>
        /// <pre>if <c> a &#8805; 0</c> return <c>t + (b/2t)i</c>
        ///  else return <c>|b|/2t + sign(b)t i </c></pre>
        ///  </item>
        /// </list>
        /// where <list type="bullet">
        /// <item><c>|a| = </c><see cref="FastMath.abs(double)"/></item>
        /// <item><c>|a + bi| = </c><see cref="abs"/>(a + bi)</item>
        /// <item><c>sign(b) =  </c><see cref="FastMath.copySign(double,double)"/>
        /// </list>
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result.
        /// <example>
        /// <code>
        ///   sqrt(1 &plusmn; INFINITY i) = INFINITY + NaN i
        ///   sqrt(INFINITY + i) = INFINITY + 0i
        ///   sqrt(-INFINITY + i) = 0 + INFINITY i
        ///   sqrt(INFINITY &plusmn; INFINITY i) = INFINITY + NaN i
        ///   sqrt(-INFINITY &plusmn; INFINITY i) = NaN &plusmn; INFINITY i
        /// </code>
        /// </example>
        /// </summary>
        /// <returns>the square root of <c>this</c>.</returns>
        public Complex sqrt()
        {
            if (isNaN)
            {
                return NaN;
            }

            if (real == 0.0 && imaginary == 0.0)
            {
                return createComplex(0.0, 0.0);
            }

            double t = FastMath.sqrt((FastMath.abs(real) + abs()) / 2.0);
            if (real >= 0.0)
            {
                return createComplex(t, imaginary / (2.0 * t));
            }
            else
            {
                return createComplex(FastMath.abs(imaginary) / (2.0 * t),
                                     FastMath.copySign(1d, imaginary) * t);
            }
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
        /// square root</a> of <c>1 - this^2</c> for this complex
        /// number.
        /// Computes the result directly as
        /// <c>sqrt(ONE.subtract(z.multiply(z)))</c>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result. 
        /// </summary>
        /// <returns>the square root of <c>1 - this^2</c>.</returns>
        public Complex sqrt1z()
        {
            return createComplex(1.0, 0.0).subtract(this.multiply(this)).sqrt();
        }

        /// <summary>
        /// Compute the
        /// a href="http://mathworld.wolfram.com/Tangent.html" TARGET="_top">
        /// tangent</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///   tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <see cref="FastMath.sin"/>, <see cref="FastMath.cos"/>,
        /// <see cref="FastMathcosh"/> and <see cref="FastMath.sinh"/>.
        /// <para/>
        /// Returns <see cref="Complex#NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite (or critical) values in real or imaginary parts of the input may
        /// result in infinite or NaN values returned in parts of the result.
        /// </summary>
        /// <example>
        /// <code>
        ///   tan(a &plusmn; INFINITY i) = 0 &plusmn; i
        ///   tan(&plusmn;INFINITY + bi) = NaN + NaN i
        ///   tan(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        ///   tan(&plusmn;&pi;/2 + 0 i) = &plusmn;INFINITY + NaN i
        /// </code>
        /// </example>
        /// <returns>the tangent of <c>this</c>.</returns>
        public Complex tan()
        {
            if (isNaN || Double.IsInfinity(real))
            {
                return NaN;
            }
            if (imaginary > 20.0)
            {
                return createComplex(0.0, 1.0);
            }
            if (imaginary < -20.0)
            {
                return createComplex(0.0, -1.0);
            }

            double real2 = 2.0 * real;
            double imaginary2 = 2.0 * imaginary;
            double d = FastMath.cos(real2) + FastMath.cosh(imaginary2);

            return createComplex(FastMath.sin(real2) / d,
                                 FastMath.sinh(imaginary2) / d);
        }

        /// <summary>
        /// Compute the
        /// <a href="http://mathworld.wolfram.com/HyperbolicTangent.html" TARGET="_top">
        /// hyperbolic tangent</a> of this complex number.
        /// Implements the formula:
        /// <code>
        ///   tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
        /// </code>
        /// where the (real) functions on the right-hand side are
        /// <see cref="FastMath.sin"/>, <see cref="FastMath.cos"/>, <see cref="FastMath.cosh"/>
        /// and <see cref="FastMath.sinh"/>.
        /// <para/>
        /// Returns <see cref="NaN"/> if either real or imaginary part of the
        /// input argument is <c>NaN</c>.
        /// <para/>
        /// Infinite values in real or imaginary parts of the input may result in
        /// infinite or NaN values returned in parts of the result.
        /// <example>
        /// <code>
        ///   tanh(a &plusmn; INFINITY i) = NaN + NaN i
        ///   tanh(&plusmn;INFINITY + bi) = &plusmn;1 + 0 i
        ///   tanh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
        ///   tanh(0 + (&pi;/2)i) = NaN + INFINITY i
        /// </code>
        /// </example>
        /// </summary>
        /// <returns>the hyperbolic tangent of <c>this</c>.</returns>
        public Complex tanh()
        {
            if (isNaN || Double.IsInfinity(imaginary))
            {
                return NaN;
            }
            if (real > 20.0)
            {
                return createComplex(1.0, 0.0);
            }
            if (real < -20.0)
            {
                return createComplex(-1.0, 0.0);
            }
            double real2 = 2.0 * real;
            double imaginary2 = 2.0 * imaginary;
            double d = FastMath.cosh(real2) + FastMath.cos(imaginary2);

            return createComplex(FastMath.sinh(real2) / d,
                                 FastMath.sin(imaginary2) / d);
        }



        /// <summary>
        /// Compute the argument of this complex number.
        /// The argument is the angle phi between the positive real axis and
        /// the point representing this number in the complex plane.
        /// The value returned is between -PI (not inclusive)
        /// and PI (inclusive), with negative values returned for numbers with
        /// negative imaginary parts.
        /// <para/>
        /// If either real or imaginary part (or both) is NaN, NaN is returned.
        /// Infinite parts are handled as <c>Math.atan2</c> handles them,
        /// essentially treating finite parts as zero in the presence of an
        /// infinite coordinate and returning a multiple of pi/4 depending on
        /// the signs of the infinite parts.
        /// See <see cref="System.Math.atan2"/> for full details.
        /// </summary>
        /// <returns>the argument of <c>this</c>.</returns>
        public double getArgument()
        {
            return FastMath.atan2(getImaginary(), getReal());
        }

        /// <summary>
        ///  Computes the n-th roots of this complex number.
        ///  The nth roots are defined by the formula:
        ///  <code>
        ///   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
        /// </code>
        /// for <i><c>k=0, 1, ..., n-1</c></i>, where <c>abs</c> and <c>phi</c>
        /// are respectively the <see cref="abs()"/> and
        /// <see cref="getArgument()"/> of this complex number.
        /// <para/>
        /// If one or both parts of this complex number is NaN, a list with just
        /// one element, <see cref="NaN"/> is returned.
        /// if neither part is NaN, but at least one part is infinite, the result
        /// is a one-element list containing <see cref="INF"/>.
        /// </summary>
        /// <param name="n">Degree of root.</param>
        /// <returns>a List<Complex> of all <c>n</c>-th roots of <c>this</c>.</returns>
        /// <exception cref="NotPositiveException">if <c>n <= 0</c>.</exception>
        public List<Complex> nthRoot(int n)
        {

            if (n <= 0)
            {
                throw new NotPositiveException<Int32>(new LocalizedFormats("CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N"),
                                               n);
            }

            List<Complex> result = new List<Complex>();

            if (isNaN)
            {
                result.Add(NaN);
                return result;
            }
            if (IsInfinity())
            {
                result.Add(INF);
                return result;
            }

            // nth root of abs -- faster / more accurate to use a solver here?
            double nthRootOfAbs = FastMath.pow(abs(), 1.0 / n);

            // Compute nth roots of complex number with k = 0, 1, ... n-1
            double nthPhi = getArgument() / n;
            double slice = 2 * FastMath.PI / n;
            double innerPart = nthPhi;
            for (int k = 0; k < n; k++)
            {
                // inner part
                double realPart = nthRootOfAbs * FastMath.cos(innerPart);
                double imaginaryPart = nthRootOfAbs * FastMath.sin(innerPart);
                result.Add(createComplex(realPart, imaginaryPart));
                innerPart += slice;
            }

            return result;
        }

        /// <summary>
        /// Create a complex number given the real and imaginary parts.
        /// </summary>
        /// <param name="realPart">Real part.</param>
        /// <param name="imaginaryPart">Imaginary part.</param>
        /// <returns>a new complex number instance.</returns>
        /// <remarks>See <see cref="valueOf(double, double)"/></remarks>
        protected Complex createComplex(double realPart,
                                        double imaginaryPart)
        {
            return new Complex(realPart, imaginaryPart);
        }

        /// <summary>
        /// Create a complex number given the real and imaginary parts.
        /// </summary>
        /// <param name="realPart">Real part.</param>
        /// <param name="imaginaryPart">Imaginary part.</param>
        /// <returns>a Complex instance.</returns>
        public static Complex valueOf(double realPart,
                                      double imaginaryPart)
        {
            if (Double.IsNaN(realPart) ||
                Double.IsNaN(imaginaryPart))
            {
                return NaN;
            }
            return new Complex(realPart, imaginaryPart);
        }

        /// <summary>
        /// Create a complex number given only the real part.
        /// </summary>
        /// <param name="realPart">Real part.</param>
        /// <returns>a Complex instance.</returns>
        public static Complex valueOf(double realPart)
        {
            if (Double.IsNaN(realPart))
            {
                return NaN;
            }
            return new Complex(realPart);
        }

        /// <summary>
        /// Resolve the transient fields in a deserialized Complex Object.
        /// Subclasses will need to override <see cref="createComplex"/> to
        /// deserialize properly.
        /// </summary>
        /// <returns>A Complex instance with all fields resolved.</returns>
        protected Object readResolve()
        {
            return createComplex(real, imaginary);
        }

        /// <inheritdoc/>
        Field<Complex> FieldElement<Complex>.getField()
        {
            return ComplexField.getInstance();
        }

        /// <inheritdoc/>
        public ComplexField getField()
        {
            return ComplexField.getInstance();
        }

        /// <inheritdoc/>
        public override String ToString()
        {
            return (String.Format("({0},{1})", this.real, this.imaginary));
        }
    }
}
