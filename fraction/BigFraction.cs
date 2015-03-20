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
using System.Numerics;

namespace Math3.fraction
{
    /// <summary>
    /// Representation of a rational number without any overflow. This class is
    /// immutable.
    /// </summary>
    public class BigFraction : FieldElement<BigFraction>, IComparable<BigFraction>
    {
        /// <summary>
        /// A fraction representing "2 / 1".
        /// </summary>
        public static readonly BigFraction TWO = new BigFraction(2);

        /// <summary>
        /// A fraction representing "1".
        /// </summary>
        public static readonly BigFraction ONE = new BigFraction(1);

        /// <summary>
        /// A fraction representing "0".
        /// </summary>
        public static readonly BigFraction ZERO = new BigFraction(0);

        /// <summary>
        /// A fraction representing "-1 / 1".
        /// </summary>
        public static readonly BigFraction MINUS_ONE = new BigFraction(-1);

        /// <summary>
        /// A fraction representing "4/5".
        /// </summary>
        public static readonly BigFraction FOUR_FIFTHS = new BigFraction(4, 5);

        /// <summary>
        /// A fraction representing "1/5".
        /// </summary>
        public static readonly BigFraction ONE_FIFTH = new BigFraction(1, 5);

        /// <summary>
        /// A fraction representing "1/2".
        /// </summary>
        public static readonly BigFraction ONE_HALF = new BigFraction(1, 2);

        /// <summary>
        /// A fraction representing "1/4".
        /// </summary>
        public static readonly BigFraction ONE_QUARTER = new BigFraction(1, 4);

        /// <summary>
        /// A fraction representing "1/3".
        /// </summary>
        public static readonly BigFraction ONE_THIRD = new BigFraction(1, 3);

        /// <summary>
        /// A fraction representing "3/5".
        /// </summary>
        public static readonly BigFraction THREE_FIFTHS = new BigFraction(3, 5);

        /// <summary>
        /// A fraction representing "3/4".
        /// </summary>
        public static readonly BigFraction THREE_QUARTERS = new BigFraction(3, 4);

        /// <summary>
        /// A fraction representing "2/5".
        /// </summary>
        public static readonly BigFraction TWO_FIFTHS = new BigFraction(2, 5);

        /// <summary>
        /// A fraction representing "2/4".
        /// </summary>
        public static readonly BigFraction TWO_QUARTERS = new BigFraction(2, 4);

        /// <summary>
        /// A fraction representing "2/3".
        /// </summary>
        public static readonly BigFraction TWO_THIRDS = new BigFraction(2, 3);

        /// <summary>
        /// <c>BigInteger</c> representation of 100.
        /// </summary>
        private static readonly BigInteger ONE_HUNDRED = new BigInteger(100);

        /// <summary>
        /// The numerator.
        /// </summary>
        private readonly BigInteger numerator;

        /// <summary>
        /// The denominator.
        /// </summary>
        private readonly BigInteger denominator;

        /// <summary>
        /// <para>
        /// Create a <see cref="BigFraction"/> equivalent to the passed <c>BigInteger</c>, ie
        /// "num / 1".
        /// </para>
        /// </summary>
        /// <param name="num"></param>
        public BigFraction(BigInteger num) : this(num, BigInteger.One) { }

        /// <summary>
        /// Create a <see cref="BigFraction"/> given the numerator and denominator as
        /// <c>BigInteger</c>. The <see cref="BigFraction"/> is reduced to lowest terms.
        /// </summary>
        /// <param name="num">the numerator, must not be <c>null</c>.</param>
        /// <param name="den">the denominator, must not be <c>null</c>.</param>
        /// <exception cref="ZeroException"> if the denominator is zero.</exception>
        /// <exception cref="NullArgumentException"> if either of the arguments is null</exception>
        public BigFraction(BigInteger num, BigInteger den)
        {
            MathUtils.checkNotNull(num, new LocalizedFormats("NUMERATOR"));
            MathUtils.checkNotNull(den, new LocalizedFormats("DENOMINATOR"));
            if (BigInteger.Zero.Equals(den))
            {
                throw new ZeroException(new LocalizedFormats("ZERO_DENOMINATOR"));
            }
            if (BigInteger.Zero.Equals(num))
            {
                numerator = BigInteger.Zero;
                denominator = BigInteger.One;
            }
            else
            {

                // reduce numerator and denominator by greatest common denominator
                BigInteger gcd = BigInteger.GreatestCommonDivisor(num, den);
                if (BigInteger.One.CompareTo(gcd) < 0)
                {
                    num = num / gcd;
                    den = den / gcd;
                }

                // move sign to numerator
                if (BigInteger.Zero.CompareTo(den) > 0)
                {
                    num = BigInteger.Negate(num);
                    den = BigInteger.Negate(den);
                }

                // store the values in the final fields
                numerator = num;
                denominator = den;

            }
        }

        /// <summary>
        /// Create a fraction given the double value.
        /// <para>
        /// This constructor behaves differently from
        /// <see cref="BigFraction(double, double, int)"/>. It converts the double value
        /// exactly, considering its internal bits representation. This works for all
        /// values except NaN and infinities and does not requires any loop or
        /// convergence threshold.
        /// </para>
        /// <para>
        /// Since this conversion is exact and since double numbers are sometimes
        /// approximated, the fraction created may seem strange in some cases. For example,
        /// calling <c>new BigFraction(1.0 / 3.0)</c> does not create
        /// the fraction 1/3, but the fraction 6004799503160661 / 18014398509481984
        /// because the double number passed to the constructor is not exactly 1/3
        /// (this number cannot be stored exactly in IEEE754).
        /// </para>
        /// </summary>
        /// <param name="value">the double value to convert to a fraction.</param>
        /// <exception cref="MathIllegalArgumentException"> if value is NaN or infinite</exception>
        public BigFraction(double value)
        {
            if (Double.IsNaN(value))
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("NAN_VALUE_CONVERSION"));
            }
            if (Double.IsInfinity(value))
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("INFINITE_VALUE_CONVERSION"));
            }

            // compute m and k such that value = m * 2^k
            long bits = BitConverter.DoubleToInt64Bits(value);
            long sign = bits & Int64.MaxValue;
            long exponent = bits & 0x7ff0000000000000L;
            long m = bits & 0x000fffffffffffffL;
            if (exponent != 0)
            {
                // this was a normalized number, add the implicit most significant bit
                m |= 0x0010000000000000L;
            }
            if (sign != 0)
            {
                m = -m;
            }
            int k = ((int)(exponent >> 52)) - 1075;
            while (((m & 0x001ffffffffffffeL) != 0) && ((m & 0x1) == 0))
            {
                m >>= 1;
                ++k;
            }

            if (k < 0)
            {
                numerator = new BigInteger(m);
                denominator = BigInteger.Pow(BigInteger.Zero, (1 << -k));
            }
            else
            {
                numerator = new BigInteger(m) * BigInteger.Pow(BigInteger.Zero, (1 << -k));
                denominator = BigInteger.One;
            }

        }

        /// <summary>
        /// Create a fraction given the double value and maximum error allowed.
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
        /// Continued Fraction</a> equations (11) and (22)-(26)</item>
        /// </list>
        /// </para>
        /// </summary>
        /// <param name="value">the double value to convert to a fraction.</param>
        /// <param name="epsilon">maximum error allowed. The resulting fraction is within
        /// <c>epsilon</c> of <c>value</c>, in absolute terms.</param>
        /// <param name="maxIterations">maximum number of convergents.</param>
        /// <exception cref="FractionConversionException">if the continued fraction failed to converge.</exception>
        /// <remarks>
        /// See <see cref="BigFraction(double)"/>
        /// </remarks>
        public BigFraction(double value, double epsilon, int maxIterations) : this(value, epsilon, Int32.MaxValue, maxIterations) { }

        /// <summary>
        /// Create a fraction given the double value and either the maximum error
        /// allowed or the maximum number of denominator digits.
        /// <para>
        /// NOTE: This constructor is called with EITHER - a valid epsilon value and
        /// the maxDenominator set to Int32.MaxValue (that way the maxDenominator
        /// has no effect). OR - a valid maxDenominator value and the epsilon value
        /// set to zero (that way epsilon only has effect if there is an exact match
        /// before the maxDenominator value is reached).
        /// </para>
        /// <para>
        /// It has been done this way so that the same code can be (re)used for both
        /// scenarios. However this could be confusing to users if it were part of
        /// the public API and this constructor should therefore remain PRIVATE.
        /// </para>
        /// See JIRA issue ticket MATH-181 for more details:
        /// https://issues.apache.org/jira/browse/MATH-181
        /// </summary>
        /// <param name="value">the double value to convert to a fraction.</param>
        /// <param name="epsilon">maximum error allowed. The resulting fraction is within
        /// <c>epsilon</c> of <code>value</code>, in absolute terms.</param>
        /// <param name="maxDenominator">maximum denominator value allowed.</param>
        /// <param name="maxIterations">maximum number of convergents.</param>
        /// <exception cref="FractionConversionException">
        /// if the continued fraction failed to converge.</exception>
        private BigFraction(double value, double epsilon, int maxDenominator, int maxIterations)
        {
            long overflow = Int32.MaxValue;
            double r0 = value;
            long a0 = (long)FastMath.floor(r0);

            if (FastMath.abs(a0) > overflow)
            {
                throw new FractionConversionException(value, a0, 1L);
            }

            // check for (almost) integer arguments, which should not go
            // to iterations.
            if (FastMath.abs(a0 - value) < epsilon)
            {
                numerator = new BigInteger(a0);
                denominator = BigInteger.One;
                return;
            }

            long p0 = 1;
            long q0 = 0;
            long p1 = a0;
            long q1 = 1;

            long p2 = 0;
            long q2 = 1;

            int n = 0;
            Boolean stop = false;
            do
            {
                ++n;
                double r1 = 1.0 / (r0 - a0);
                long a1 = (long)FastMath.floor(r1);
                p2 = (a1 * p1) + p0;
                q2 = (a1 * q1) + q0;
                if ((p2 > overflow) || (q2 > overflow))
                {
                    // in maxDenominator mode, if the last fraction was very close to the actual value
                    // q2 may overflow in the next iteration; in this case return the last one.
                    if (epsilon == 0.0 && FastMath.abs(q1) < maxDenominator)
                    {
                        break;
                    }
                    throw new FractionConversionException(value, p2, q2);
                }

                double convergent = (double)p2 / (double)q2;
                if ((n < maxIterations) &&
                    (FastMath.abs(convergent - value) > epsilon) &&
                    (q2 < maxDenominator))
                {
                    p0 = p1;
                    p1 = p2;
                    q0 = q1;
                    q1 = q2;
                    a0 = a1;
                    r0 = r1;
                }
                else
                {
                    stop = true;
                }
            } while (!stop);

            if (n >= maxIterations)
            {
                throw new FractionConversionException(value, maxIterations);
            }

            if (q2 < maxDenominator)
            {
                numerator = new BigInteger(p2);
                denominator = new BigInteger(q2);
            }
            else
            {
                numerator = new BigInteger(p1);
                denominator = new BigInteger(q1);
            }
        }

        /// <summary>
        /// Create a fraction given the double value and maximum denominator.
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
        /// Continued Fraction</a> equations (11) and (22)-(26)</item>
        /// </list>
        /// </para>
        /// </summary>
        /// <param name="value">the double value to convert to a fraction.</param>
        /// <param name="maxDenominator">The maximum allowed value for denominator.</param>
        /// <exception cref="FractionConversionException"> if the continued fraction failed
        /// to converge.</exception>
        public BigFraction(double value, int maxDenominator) : this(value, 0, maxDenominator, 100) { }

        /// <summary>
        /// <para>
        /// Create a <see cref="BigFraction"/> equivalent to the passed <c>int</c>, ie
        /// "num / 1".
        /// </para>
        /// </summary>
        /// <param name="num">the numerator.</param>
        public BigFraction(int num) : this(new BigInteger(num), BigInteger.One) { }

        /// <summary>
        /// <para>
        /// Create a <see cref="BigFraction"/> given the numerator and denominator as simple
        /// <c>int</c>. The <see cref="BigFraction"/> is reduced to lowest terms.
        /// </para>  
        /// </summary>
        /// <param name="num">the numerator.</param>
        /// <param name="den">the denominator.</param>
        public BigFraction(int num, int den) : this(new BigInteger(num), new BigInteger(den)) { }

        /// <summary>
        /// <para>
        /// Create a <see cref="BigFraction"/> equivalent to the passed long, ie "num / 1".
        /// </para>    
        /// </summary>
        /// <param name="num">the numerator.</param>
        public BigFraction(long num) : this(new BigInteger(num), BigInteger.One) { }

        /// <summary>
        /// <para>
        /// Create a <see cref="BigFraction"/> given the numerator and denominator as simple
        /// <c>long</c>. The <see cref="BigFraction"/> is reduced to lowest terms.
        /// </para> 
        /// </summary>
        /// <param name="num">the numerator.</param>
        /// <param name="den">the denominator.</param>
        public BigFraction(long num, long den) : this(new BigInteger(num), new BigInteger(den)) { }

        /// <summary>
        /// <para>
        /// Creates a <c>BigFraction</c> instance with the 2 parts of a fraction
        /// Y/Z.
        /// </para>
        /// <para>
        /// Any negative signs are resolved to be on the numerator.
        /// </para>
        /// </summary>
        /// <param name="numerator">the numerator, for example the three in 'three sevenths'.</param>
        /// <param name="denominator">the denominator, for example the seven in 'three sevenths'.</param>
        /// <returns>a new fraction instance, with the numerator and denominator
        /// reduced.</returns>
        /// <exception cref="ArithmeticException"> if the denominator is <c>zero</c>.</exception>
        public static BigFraction getReducedFraction(int numerator, int denominator)
        {
            if (numerator == 0)
            {
                return ZERO; // normalize zero.
            }

            return new BigFraction(numerator, denominator);
        }

        /// <summary>
        /// <para>
        /// Returns the absolute value of this <see cref="BigFraction"/>.
        /// </para>
        /// </summary>
        /// <returns>the absolute value as a <see cref="BigFraction"/>.</returns>
        public BigFraction abs()
        {
            return (BigInteger.Zero.CompareTo(numerator) <= 0) ? this : negate();
        }

        /// <summary>
        /// <para>
        /// Adds the value of this fraction to the passed <see cref="BigInteger"/>,
        /// returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="bg">the <see cref="BigInteger"/> to add, must'nt be <c>null</c>.</param>
        /// <returns>a <c>BigFraction</c> instance with the resulting values.</returns>
        /// <exception cref="NullArgumentException"> if the <see cref="BigInteger"/> is
        /// <c>null</c>.</exception>
        public BigFraction add(BigInteger bg)
        {
            MathUtils.checkNotNull(bg);
            return new BigFraction(numerator + (denominator * bg), denominator);
        }

        /// <summary>
        /// <para>
        /// Adds the value of this fraction to the passed <c>integer</c>, returning
        /// the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="i">the <c>integer</c> to add.</param>
        /// <returns>a <c>BigFraction</c> instance with the resulting values.</returns>
        public BigFraction add(int i)
        {
            return add(new BigInteger(i));
        }

        /// <summary>
        /// <para>
        /// Adds the value of this fraction to the passed <c>long</c>, returning
        /// the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="l">the <c>long</c> to add.</param>
        /// <returns>a <c>BigFraction</c> instance with the resulting values.</returns>
        public BigFraction add(long l)
        {
            return add(new BigInteger(l));
        }

        /// <summary>
        /// <para>
        /// Adds the value of this fraction to another, returning the result in
        /// reduced form.
        /// </para>
        /// </summary>
        /// <param name="fraction">the <see cref="BigFraction"/> to add, must not be <c>null</c>.</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values.</returns>
        /// <exception cref="NullArgumentException"> if the <see cref="BigFraction"/> is 
        /// <c>null</c>.</exception>
        public BigFraction add(BigFraction fraction) {
        if (fraction == null) {
            throw new NullArgumentException(new LocalizedFormats("FRACTION"));
        }
        if (ZERO.Equals(fraction)) {
            return this;
        }

        BigInteger num = default(BigInteger);
        BigInteger den = default(BigInteger);

        if (denominator.Equals(fraction.denominator)) {
            num = numerator + fraction.numerator;
            den = denominator;
        } else {
            num = (numerator * fraction.denominator) + (fraction.numerator * denominator);
            den = denominator *fraction.denominator;
        }
        return new BigFraction(num, den);

    }

        /// <summary>
        /// <para>
        /// Gets the fraction as a <c>BigDecimal</c>. This calculates the
        /// fraction as the numerator divided by denominator.
        /// </para>
        /// </summary>
        /// <returns>the fraction as a <c>BigDecimal</c>.</returns>
        public BigDecimal bigDecimalValue()
        {
            return new BigDecimal(numerator).divide(new BigDecimal(denominator));
        }

        /// <summary>
        /// <para>
        /// Gets the fraction as a <c>BigDecimal</c> following the passed
        /// rounding mode. This calculates the fraction as the numerator divided by
        /// denominator.
        /// </para>
        /// </summary>
        /// <param name="roundingMode">rounding mode to apply. see <see cref="BigDecimal"/> constants.</param>
        /// <returns>the fraction as a <c>BigDecimal</c>.</returns>
        /// <exception cref="ArgumentException">if <c>roundingMode</c> does not represent a
        /// valid rounding mode.</exception>
        /// <remarks>
        /// See <see cref="BigDecimal"/>
        /// </remarks>
        public BigDecimal bigDecimalValue(int roundingMode)
        {
            return new BigDecimal(numerator).divide(new BigDecimal(denominator), roundingMode);
        }

        /// <summary>
        /// <para>
        /// Gets the fraction as a <c>BigDecimal</c> following the passed scale
        /// and rounding mode. This calculates the fraction as the numerator divided
        /// by denominator.
        /// </para>
        /// </summary>
        /// <param name="scale">scale of the <c>BigDecimal</c> quotient to be returned.
        /// see <see cref="BigDecimal"/> for more information.</param>
        /// <param name="roundingMode">rounding mode to apply. see <see cref="BigDecimal"/> 
        /// constants.</param>
        /// <returns><the fraction as a <c>BigDecimal</c>./returns>
        /// <remarks>
        /// See <see cref="BigDecimal"/>
        /// </remarks>
        public BigDecimal bigDecimalValue(int scale, int roundingMode)
        {
            return new BigDecimal(numerator).divide(new BigDecimal(denominator), scale, roundingMode);
        }

        /// <summary>
        /// <para>
        /// Compares this object to another based on size.
        /// </para>
        /// </summary>
        /// <param name="obj">the object to compare to, must not be <c>null</c>.</param>
        /// <returns>-1 if this is less than <c>object</c>, +1 if this is greater
        /// than <c>object</c>, 0 if they are equal.</returns>
        public int CompareTo(BigFraction obj)
        {
            BigInteger nOd = numerator * (obj.denominator);
            BigInteger dOn = denominator * (obj.numerator);
            return nOd.CompareTo(dOn);
        }

        /// <summary>
        /// <para>
        /// Divide the value of this fraction by the passed <c>BigInteger</c>,
        /// ie <c>this * 1 / bg</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="bg">the <c>BigInteger</c> to divide by, must not be <c>null</c></param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values</returns>
        /// <exception cref="NullArgumentException"> if the <c>BigInteger</c> is <c>null</c>
        /// </exception>
        /// <exception cref="MathArithmeticException"> if the fraction to divide by is zero
        /// </exception>
        public BigFraction divide(BigInteger bg)
        {
            if (bg == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FRACTION"));
            }
            if (BigInteger.Zero.Equals(bg))
            {
                throw new MathArithmeticException(new LocalizedFormats("ZERO_DENOMINATOR"));
            }
            return new BigFraction(numerator, denominator * bg);
        }

        /// <summary>
        /// <para>
        /// Divide the value of this fraction by the passed <c>int</c>, ie
        /// <c>this * 1 / i</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="i">the <c>int</c> to divide by</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values</returns>
        /// <exception cref="MathArithmeticException"> if the fraction to divide by is zero
        /// </exception>
        public BigFraction divide(int i)
        {
            return divide(new BigInteger(i));
        }

        /// <summary>
        /// <para>
        /// Divide the value of this fraction by the passed <c>long</c>, ie
        /// <c>this * 1 / l</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="l">the <c>long</c> to divide by</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values</returns>
        /// <exception cref="MathArithmeticException"> if the fraction to divide by is zero
        /// </exception>
        public BigFraction divide(long l)
        {
            return divide(new BigInteger(l));
        }

        /// <summary>
        /// <para>
        /// Divide the value of this fraction by another, returning the result in
        /// reduced form.
        /// </para>
        /// </summary>
        /// <param name="fraction">Fraction to divide by, must not be <c>null</c>.</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values.</returns>
        /// <exception cref="NullArgumentException"> if the <c>fraction</c> is <c>null</c>.
        /// </exception>
        /// <exception cref="MathArithmeticException"> if the fraction to divide by is zero
        /// </exception>
        public BigFraction divide(BigFraction fraction)
        {
            if (fraction == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FRACTION"));
            }
            if (BigInteger.Zero.Equals(fraction.numerator))
            {
                throw new MathArithmeticException(new LocalizedFormats("ZERO_DENOMINATOR"));
            }

            return multiply(fraction.reciprocal());
        }

        /// <summary>
        /// <para>
        /// Gets the fraction as a <c>double</c>. This calculates the fraction as
        /// the numerator divided by denominator.
        /// </para>
        /// </summary>
        /// <returns>the fraction as a <c>double</c></returns>
        public double doubleValue()
        {
            double result = Convert.ToDouble(numerator) / Convert.ToDouble(denominator);
            if (Double.IsNaN(result))
            {
                // Numerator and/or denominator must be out of range:
                // Calculate how far to shift them to put them in range.
                int shift = FastMath.max(numerator.ToByteArray().Length / 8, denominator.ToByteArray().Length / 8) - FastMath.getExponent(Double.MaxValue);
                result = Convert.ToDouble(numerator >> shift) /
                    Convert.ToDouble(denominator >> shift);
            }
            return result;
        }

        /// <summary>
        /// <para>
        /// Test for the equality of two fractions. If the lowest term numerator and
        /// denominators are the same for both fractions, the two fractions are
        /// considered to be equal.
        /// </para>
        /// </summary>
        /// <param name="other">fraction to test for equality to this fraction, can be
        /// <c>null</c>.</param>
        /// <returns>true if two fractions are equal, false if object is
        /// <c>null</c>, not an instance of <see cref="BigFraction"/>, or not
        /// equal to this fraction instance.</returns>
        public override Boolean Equals(Object other)
        {
            Boolean ret = false;

            if (this == other)
            {
                ret = true;
            }
            else if (other is BigFraction)
            {
                BigFraction rhs = ((BigFraction)other).reduce();
                BigFraction thisOne = this.reduce();
                ret = thisOne.numerator.Equals(rhs.numerator) && thisOne.denominator.Equals(rhs.denominator);
            }
            return ret;
        }

        /// <summary>
        /// <para>
        /// Gets the fraction as a <c>float</c>. This calculates the fraction as
        /// the numerator divided by denominator.
        /// </para>
        /// </summary>
        /// <returns>the fraction as a <c>float</c>.</returns>
        public float floatValue()
        {
            float result = Convert.ToSingle(numerator) / Convert.ToSingle(denominator);
            if (Single.IsNaN(result))
            {
                // Numerator and/or denominator must be out of range:
                // Calculate how far to shift them to put them in range.
                int shift = FastMath.max(numerator.ToByteArray().Length / 8,
                                         denominator.ToByteArray().Length / 8) - FastMath.getExponent(Single.MaxValue);
                result = Convert.ToSingle(numerator >> shift) / Convert.ToSingle(denominator >> shift);
            }
            return result;
        }

        /// <summary>
        /// <para>
        /// Access the denominator as a <c>BigInteger</c>.
        /// </para>
        /// </summary>
        /// <returns>the denominator as a <c>BigInteger</c>.</returns>
        public BigInteger getDenominator()
        {
            return denominator;
        }

        /// <summary>
        /// <para>
        /// Access the denominator as a <c>int</c>.
        /// </para>
        /// </summary>
        /// <returns>the denominator as a <c>int</c>.</returns>
        public int getDenominatorAsInt()
        {
            return Convert.ToInt32(denominator);
        }

        /// <summary>
        /// <para>
        /// Access the denominator as a <c>long</c>.
        /// </para>
        /// </summary>
        /// <returns>the denominator as a <c>long</c>.</returns>
        public long getDenominatorAsLong()
        {
            return Convert.ToInt64(denominator);
        }

        /// <summary>
        /// <para>
        /// Access the numerator as a <c>BigInteger</c>.
        /// </para>
        /// </summary>
        /// <returns>the numerator as a <c>BigInteger</c>.</returns>
        public BigInteger getNumerator()
        {
            return numerator;
        }

        /// <summary>
        /// <para>
        /// Access the numerator as a <c>int</c>.
        /// </para>
        /// </summary>
        /// <returns>the numerator as a <c>int</c>.</returns>
        public int getNumeratorAsInt()
        {
            return Convert.ToInt32(numerator);
        }

        /// <summary>
        /// <para>
        /// Access the numerator as a <c>long</c>.
        /// </para>
        /// </summary>
        /// <returns>the numerator as a <c>long</c>.</returns>
        public long getNumeratorAsLong()
        {
            return Convert.ToInt64(numerator);
        }

        /// <summary>
        /// <para>
        /// Gets a HashCode for the fraction.
        /// </para>
        /// </summary>
        /// <returns>a hash code value for this object.</returns>
        public override int GetHashCode()
        {
            return 37 * (37 * 17 + numerator.GetHashCode()) + denominator.GetHashCode();
        }

        /// <summary>
        /// <para>
        /// Gets the fraction as an <c>int</c>. This returns the whole number part
        /// of the fraction.
        /// </para>
        /// </summary>
        /// <returns>the whole number fraction part.</returns>
        public int intValue()
        {
            return Convert.ToInt32(numerator / denominator);
        }

        /// <summary>
        /// <para>
        /// Gets the fraction as a <c>long</c>. This returns the whole number part
        /// of the fraction.
        /// </para>
        /// </summary>
        /// <returns>the whole number fraction part.</returns>
        public long longValue()
        {
            return Convert.ToInt64(numerator / denominator);
        }

        /// <summary>
        /// <para>
        /// Multiplies the value of this fraction by the passed
        /// <c>BigInteger</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="bg">the <c>BigInteger</c> to multiply by.</param>
        /// <returns>a <c>BigFraction</c> instance with the resulting values.</returns>
        /// <exception cref="NullArgumentException"> if <c>bg</c> is <c>null</c>.</exception>
        public BigFraction multiply(BigInteger bg)
        {
            if (bg == null)
            {
                throw new NullArgumentException();
            }
            return new BigFraction(bg * numerator, denominator);
        }

        /// <summary>
        /// <para>
        /// Multiply the value of this fraction by the passed <c>int</c>, returning
        /// the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="i">the <c>int</c> to multiply by.</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values.</returns>
        public BigFraction multiply(int i)
        {
            return multiply(new BigInteger(i));
        }

        /// <summary>
        /// <para>
        /// Multiply the value of this fraction by the passed <c>long</c>,
        /// returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="l">the <c>long</c> to multiply by.</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values.</returns>
        public BigFraction multiply(long l)
        {
            return multiply(new BigInteger(l));
        }

        /// <summary>
        /// <para>
        /// Multiplies the value of this fraction by another, returning the result in
        /// reduced form.
        /// </para>
        /// </summary>
        /// <param name="fraction">Fraction to multiply by, must not be <c>null</c>.</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values.</returns>
        /// <exception cref="NullArgumentException"> if <c>fraction</c> is <c>null</c>.</exception>
        public BigFraction multiply(BigFraction fraction)
        {
            if (fraction == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FRACTION"));
            }
            if (numerator.Equals(BigInteger.Zero) ||
                fraction.numerator.Equals(BigInteger.Zero))
            {
                return ZERO;
            }
            return new BigFraction(numerator * fraction.numerator, denominator * fraction.denominator);
        }

        /// <summary>
        /// <para>
        /// Return the additive inverse of this fraction, returning the result in
        /// reduced form.
        /// </para>
        /// </summary>
        /// <returns>the negation of this fraction.</returns>
        public BigFraction negate()
        {
            return new BigFraction(BigInteger.Negate(numerator), denominator);
        }

        /// <summary>
        /// <para>
        /// Gets the fraction percentage as a <c>double</c>. This calculates the
        /// fraction as the numerator divided by denominator multiplied by 100.
        /// </para>
        /// </summary>
        /// <returns>the fraction percentage as a <c>double</c>.</returns>
        public double percentageValue()
        {
            return multiply(ONE_HUNDRED).doubleValue();
        }

        /// <summary>
        /// <para>
        /// Returns a <c>BigFraction</c> whose value is
        /// <c>(this<sup>exponent</sup>)</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="exponent">exponent to which this <c>BigFraction</c> is to be raised.
        /// </param>
        /// <returns>this^exponent</returns>
        public BigFraction pow(int exponent)
        {
            if (exponent < 0)
            {
                return new BigFraction(BigInteger.Pow(denominator, -exponent), BigInteger.Pow(numerator, -exponent));
            }
            return new BigFraction(BigInteger.Pow(numerator, exponent), BigInteger.Pow(denominator, exponent));
        }

        /// <summary>
        /// <para>
        /// Returns a <c>BigFraction</c> whose value is
        /// (this^exponent), returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="exponent">exponent to which this <c>BigFraction</c> is to be raised.
        /// </param>
        /// <returns>this^exponent as a <c>BigFraction</c>.</returns>
        public BigFraction pow(long exponent)
        {
            if (exponent < 0)
            {
                return new BigFraction(ArithmeticUtils.pow(denominator, -exponent),
                                       ArithmeticUtils.pow(numerator, -exponent));
            }
            return new BigFraction(ArithmeticUtils.pow(numerator, exponent),
                                   ArithmeticUtils.pow(denominator, exponent));
        }

        /// <summary>
        /// <para>
        /// Returns a <c>BigFraction</c> whose value is
        /// (this^exponent), returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="exponent">exponent to which this <c>BigFraction</c> is to be raised.
        /// </param>
        /// <returns>this^exponent as a <c>BigFraction</c>.</returns>
        public BigFraction pow(BigInteger exponent)
        {
            if (exponent.CompareTo(BigInteger.Zero) < 0)
            {
                BigInteger eNeg = BigInteger.Negate(exponent);
                return new BigFraction(ArithmeticUtils.pow(denominator, eNeg),
                                       ArithmeticUtils.pow(numerator, eNeg));
            }
            return new BigFraction(ArithmeticUtils.pow(numerator, exponent),
                                   ArithmeticUtils.pow(denominator, exponent));
        }

        /// <summary>
        /// <para>
        /// Returns a <c>double</c> whose value is
        /// (this^exponent), returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="exponent">exponent to which this <c>BigFraction</c> is to be raised.
        /// </param>
        /// <returns>this^exponent.</returns>
        public double pow(double exponent)
        {
            return FastMath.pow(Convert.ToDouble(numerator), exponent) /
                   FastMath.pow(Convert.ToDouble(denominator), exponent);
        }

        /// <summary>
        /// <para>
        /// Return the multiplicative inverse of this fraction.
        /// </para>
        /// </summary>
        /// <returns>the reciprocal fraction.</returns>
        public BigFraction reciprocal()
        {
            return new BigFraction(denominator, numerator);
        }

        /// <summary>
        /// <para>
        /// Reduce this <c>BigFraction</c> to its lowest terms.
        /// </para>
        /// </summary>
        /// <returns>the reduced <c>BigFraction</c>. It doesn't change anything if
        /// the fraction can be reduced.</returns>
        public BigFraction reduce()
        {
            BigInteger gcd = BigInteger.GreatestCommonDivisor(numerator, denominator);
            return new BigFraction(numerator / gcd, denominator / gcd);
        }

        /// <summary>
        /// <para>
        /// Subtracts the value of an <see cref="BigInteger"/> from the value of this
        /// <c>BigFraction</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="bg">the <see cref="BigInteger"/> to subtract.</param>
        /// <returns><c>BigFraction</c> instance with the resulting values.</returns>
        public BigFraction subtract(BigInteger bg)
        {
            return new BigFraction(numerator - (denominator * bg), denominator);
        }

        /// <summary>
        /// <para>
        /// Subtracts the value of an <c>integer</c> from the value of this
        /// <c>BigFraction</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="i">the <c>integer</c> to subtract.</param>
        /// <returns>a <c>BigFraction</c> instance with the resulting values.</returns>
        public BigFraction subtract(int i)
        {
            return subtract(new BigInteger(i));
        }

        /// <summary>
        /// <para>
        /// Subtracts the value of a <c>long</c> from the value of this
        /// <c>BigFraction</c>, returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="l">the <c>long</c> to subtract.</param>
        /// <returns>a <c>BigFraction</c> instance with the resulting values.</returns>
        public BigFraction subtract(long l)
        {
            return subtract(new BigInteger(l));
        }

        /// <summary>
        /// <para>
        /// Subtracts the value of another fraction from the value of this one,
        /// returning the result in reduced form.
        /// </para>
        /// </summary>
        /// <param name="fraction"><see cref="BigFraction"/> to subtract, must not be <c>null
        /// </c>.</param>
        /// <returns>a <see cref="BigFraction"/> instance with the resulting values</returns>
        /// <exception cref="NullArgumentException"> if the <c>fraction</c> is <c>null</c>.
        /// </exception>
        public BigFraction subtract(BigFraction fraction)
        {
            if (fraction == null)
            {
                throw new NullArgumentException(new LocalizedFormats("FRACTION"));
            }
            if (ZERO.Equals(fraction))
            {
                return this;
            }

            BigInteger num = default(BigInteger);
            BigInteger den = default(BigInteger);
            if (denominator.Equals(fraction.denominator))
            {
                num = numerator - fraction.numerator;
                den = denominator;
            }
            else
            {
                num = (numerator * fraction.denominator) - (fraction.numerator * denominator);
                den = denominator * fraction.denominator;
            }
            return new BigFraction(num, den);

        }

        /// <summary>
        /// <para>
        /// Returns the <code>String</code> representing this fraction, ie
        /// "num / dem" or just "num" if the denominator is one.
        /// </para>
        /// </summary>
        /// <returns>a string representation of the fraction.</returns>
        public override String ToString()
        {
            String str;
            if (BigInteger.One.Equals(denominator))
            {
                str = numerator.ToString();
            }
            else if (BigInteger.Zero.Equals(numerator))
            {
                str = "0";
            }
            else
            {
                str = numerator + " / " + denominator;
            }
            return str;
        }

        /// <inheritdoc/>
        public Field<BigFraction> getField()
        {
            return BigFractionField.getInstance();
        }
    }
}
