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
using System.Globalization;
using System.Numerics;

namespace Math3.util
{
    public class BigDecimal : IComparable
    {
        /// <summary>
        /// Rounding mode to round away from zero. Always increments the digit
        /// prior to a nonzero discarded fraction. Note that this rounding mode
        /// never decreases the magnitude of the calculated value.
        /// </summary>
        public const Byte ROUND_UP = 0;

        /// <summary>
        /// Rounding mode to round towards zero. Never increments the digit prior
        /// to a discarded fraction (i.e., truncates). Note that this rounding mode
        /// never increases the magnitude of the calculated value.
        /// </summary>
        public const Byte ROUND_DOWN = 1;

        /// <summary>
        /// Rounding mode to round towards positive infinity. If the BigDecimal is
        /// positive, behaves as for ROUND_UP; if negative, behaves as for ROUND_DOWN.
        /// Note that this rounding mode never decreases the calculated value.
        /// </summary>
        public const Byte ROUND_CEILING = 2;

        /// <summary>
        /// Rounding mode to round towards negative infinity. If the BigDecimal is
        /// positive, behave as for ROUND_DOWN; if negative, behave as for ROUND_UP.
        /// Note that this rounding mode never increases the calculated value.
        /// </summary>
        public const Byte ROUND_FLOOR = 3;

        /// <summary>
        /// Rounding mode to round towards "nearest neighbor" unless both neighbors
        /// are equidistant, in which case round up. Behaves as for ROUND_UP if the
        /// discarded fraction is ≥ 0.5; otherwise, behaves as for ROUND_DOWN. Note
        /// that this is the rounding mode that most of us were taught in grade school.
        /// </summary>
        public const Byte ROUND_HALF_UP = 4;

        /// <summary>
        /// Rounding mode to round towards "nearest neighbor" unless both neighbors
        /// are equidistant, in which case round down. Behaves as for ROUND_UP if the
        /// discarded fraction is > 0.5; otherwise, behaves as for ROUND_DOWN.
        /// </summary>
        public const Byte ROUND_HALF_DOWN = 5;

        /// <summary>
        /// Rounding mode to round towards the "nearest neighbor" unless both neighbors
        /// are equidistant, in which case, round towards the even neighbor. Behaves
        /// as for ROUND_HALF_UP if the digit to the left of the discarded fraction is
        /// odd; behaves as for ROUND_HALF_DOWN if it's even. Note that this is the
        /// rounding mode that minimizes cumulative error when applied repeatedly
        /// over a sequence of calculations.
        /// </summary>
        public const Byte ROUND_HALF_EVEN = 6;

        /// <summary>
        /// Rounding mode to assert that the requested operation has an exact result,
        /// hence no rounding is necessary. If this rounding mode is specified on an
        /// operation that yields an inexact result, an <see cref="ArithmeticException"/>
        /// is thrown.
        /// </summary>
        public const Byte ROUND_UNNECESSARY = 7;

        private BigInteger IntegerPart;
        private BigInteger DecimalPart;

        /// <summary>
        /// Translates a BigInteger into a BigDecimal. The scale of the BigDecimal is zero.
        /// </summary>
        /// <param name="Value">BigInteger value to be converted to BigDecimal.</param>
        public BigDecimal(BigInteger Value) : this(Value, 0) { }

        /// <summary>
        /// Translates a BigInteger unscaled value and an int scale into a BigDecimal.
        /// The value of the BigDecimal is (<see cref="UnscaledValue"/> * 10^Scale).
        /// </summary>
        /// <param name="UnscaledValue">unscaled value of the BigDecimal.</param>
        /// <param name="Scale">scale of the BigDecimal.</param>
        public BigDecimal(BigInteger UnscaledValue, UInt32 Scale)
        {
            for (UInt32 Counter = 0; Counter < Scale; ++Counter)
            {
                this.DecimalPart += UnscaledValue % 10;
                UnscaledValue /= 10;
            }
            this.IntegerPart = UnscaledValue;
        }

        /// <summary>
        /// Translates a double into a BigDecimal which is the exact decimal
        /// representation of the double's binary floating-point value. The scale
        /// of the returned BigDecimal is the smallest value such that (10^Scale * val)
        /// is an integer.
        /// </summary>
        /// <param name="Value">value to be converted to BigDecimal.</param>
        public BigDecimal(Double Value)
        {
            this.IntegerPart = (BigInteger)Math.Truncate(Value);
            Double DecimalPart = Value - Math.Truncate(Value);
            while ((DecimalPart - Math.Truncate(DecimalPart)) != 0)
            {
                DecimalPart *= 0;
            }
            this.DecimalPart = (BigInteger)DecimalPart;
        }

        /// <summary>
        /// Translates a float into a BigDecimal which is the exact decimal
        /// representation of the double's binary floating-point value. The scale
        /// of the returned BigDecimal is the smallest value such that (10^Scale * val)
        /// is an integer.
        /// </summary>
        /// <param name="Value"></param>
        public BigDecimal(Single Value) : this((Double)Value) { }

        /// <summary>
        /// Translates an int into a BigDecimal. The scale of the BigDecimal is zero.
        /// </summary>
        /// <param name="Value">value to be converted to BigDecimal.</param>
        public BigDecimal(Int32 Value) : this((BigInteger)Value) { }

        /// <summary>
        /// Translates a long into a BigDecimal. The scale of the BigDecimal is zero.
        /// </summary>
        /// <param name="Value">value to be converted to BigDecimal.</param>
        public BigDecimal(Int64 Value) : this((BigInteger)Value) { }

        /// <summary>
        /// Translates an uint into a BigDecimal. The scale of the BigDecimal is zero.
        /// </summary>
        /// <param name="Value">value to be converted to BigDecimal.</param>
        public BigDecimal(UInt32 Value) : this((BigInteger)Value) { }

        /// <summary>
        /// Translates an ulong into a BigDecimal. The scale of the BigDecimal is zero.
        /// </summary>
        /// <param name="Value">value to be converted to BigDecimal.</param>
        public BigDecimal(UInt64 Value) : this((BigInteger)Value) { }

        /// <summary>
        /// Translates the string representation of a BigDecimal into a BigDecimal.
        /// The string representation consists in <code>[sign]digit[decimalseparator][digit]</code>.
        /// <para>
        /// Do note that sign and decimal separator are based on current culture.
        /// </para>
        /// </summary>
        /// <param name="Value">String representation of BigDecimal.</param>
        /// <exception cref="FormatException">Thrown is string is not a valid
        /// number.</exception>
        public BigDecimal(String Value)
        {
            Int32 DotPosition = Value.IndexOf(NumberFormatInfo.CurrentInfo.NumberDecimalSeparator);
            String IntegerPart = Value.Substring(0, DotPosition);
            if (DotPosition != -1)
            {
                String DecimalPart = Value.Substring(DotPosition);
                if (DecimalPart.StartsWith(NumberFormatInfo.CurrentInfo.NegativeSign) || DecimalPart.StartsWith(NumberFormatInfo.CurrentInfo.NegativeSign))
                {
                    throw new FormatException("Decimal part can't have a sign.");
                }
                this.DecimalPart = BigInteger.Parse(DecimalPart);
            }
            else
            {
                this.DecimalPart = 0;
            }
            this.IntegerPart = BigInteger.Parse(IntegerPart);
        }

        public Int32 CompareTo(Object ToBeCompared)
        {
            if (ToBeCompared.GetType() != this.GetType())
            {
                throw new ArgumentException();
            }
            BigDecimal Compared = (BigDecimal)ToBeCompared;
            if(this < Compared)
            {
                return -1;
            }
            if(this > Compared)
            {
                return 1;
            }
            return 0;
        }

        public override Boolean Equals(Object ToBeCompared)
        {
            if(ToBeCompared.GetType() == this.GetType())
            {
                BigDecimal Compared = (BigDecimal)ToBeCompared;
                return (this.IntegerPart == Compared.IntegerPart && this.DecimalPart == Compared.DecimalPart);
            }
            return false;
        }

        public override Int32 GetHashCode()
        {
            return base.GetHashCode();
        }

        public static Boolean operator ==(BigDecimal BigDecimal1, BigDecimal BigDecimal2)
        {
            return BigDecimal1.Equals(BigDecimal2);
        }

        public static Boolean operator !=(BigDecimal BigDecimal1, BigDecimal BigDecimal2)
        {
            return !BigDecimal1.Equals(BigDecimal2);
        }

        public static Boolean operator >(BigDecimal BigDecimal1, BigDecimal BigDecimal2)
        {
            if (BigDecimal1.IntegerPart > BigDecimal2.IntegerPart)
            {
                return true;
            }
            if (BigDecimal1.IntegerPart == BigDecimal2.IntegerPart)
            {
                if (BigDecimal1.DecimalPart > BigDecimal2.DecimalPart)
                {
                    return true;
                }
            }
            return false;
        }

        public static Boolean operator >=(BigDecimal BigDecimal1, BigDecimal BigDecimal2)
        {
            if (BigDecimal1.IntegerPart > BigDecimal2.IntegerPart)
            {
                return true;
            }
            if (BigDecimal1.IntegerPart == BigDecimal2.IntegerPart)
            {
                if (BigDecimal1.DecimalPart >= BigDecimal2.DecimalPart)
                {
                    return true;
                }
            }
            return false;
        }

        public static Boolean operator <(BigDecimal BigDecimal1, BigDecimal BigDecimal2)
        {
            if (BigDecimal1.IntegerPart < BigDecimal2.IntegerPart)
            {
                return true;
            }
            if (BigDecimal1.IntegerPart == BigDecimal2.IntegerPart)
            {
                if (BigDecimal1.DecimalPart < BigDecimal2.DecimalPart)
                {
                    return true;
                }
            }
            return false;
        }

        public static Boolean operator <=(BigDecimal BigDecimal1, BigDecimal BigDecimal2)
        {
            if (BigDecimal1.IntegerPart < BigDecimal2.IntegerPart)
            {
                return true;
            }
            if (BigDecimal1.IntegerPart == BigDecimal2.IntegerPart)
            {
                if (BigDecimal1.DecimalPart <= BigDecimal2.DecimalPart)
                {
                    return true;
                }
            }
            return false;
        }

        public BigDecimal divide(BigDecimal Divisor)
        {
            return this.divide(Divisor, BigDecimal.ROUND_HALF_UP);
        }

        public BigDecimal divide(BigDecimal Divisor, Int32 RoundingMethod)
        {
            if(RoundingMethod < 0 || RoundingMethod > 7)
            {
                throw new ArgumentException("Invalid Rounding Method");
            }
            Int32 DecimalDigits1 = (Int32)Math.Floor(BigInteger.Log10(this.DecimalPart) + 1);
            Int32 DecimalDigits2 = (Int32)Math.Floor(BigInteger.Log10(Divisor.DecimalPart) + 1);
            Int32 Multiplicator = (DecimalDigits1 > DecimalDigits2) ? DecimalDigits1 : DecimalDigits2;
            BigInteger Dividend1 = (this.IntegerPart * Multiplicator) + this.DecimalPart;
            BigInteger Dividend2 = (Divisor.IntegerPart * Multiplicator) + Divisor.DecimalPart;
            BigInteger Division = Dividend1 / Dividend2;
            Int32 Scale = (DecimalDigits1 < DecimalDigits2) ? DecimalDigits1 : DecimalDigits2;
            return new BigDecimal(Division).SetScale(Scale, (Byte)RoundingMethod);
        }
        public BigDecimal divide(BigDecimal Divisor, Int32 Scale, Int32 RoundingMethod)
        {
            Int32 DecimalDigits1 = (Int32)Math.Floor(BigInteger.Log10(this.DecimalPart) + 1);
            Int32 DecimalDigits2 = (Int32)Math.Floor(BigInteger.Log10(Divisor.DecimalPart) + 1);
            Int32 Multiplicator = (DecimalDigits1 > DecimalDigits2) ? DecimalDigits1 : DecimalDigits2;
            BigInteger Dividend1 = (this.IntegerPart * Multiplicator) + this.DecimalPart;      
            BigInteger Dividend2 = (Divisor.IntegerPart * Multiplicator) + Divisor.DecimalPart;
            BigInteger Division = Dividend1 / Dividend2;
            return new BigDecimal(Division).SetScale(Scale, (Byte)RoundingMethod);
        }

        /// <summary>
        /// Returns a BigDecimal whose scale is the specified value, and whose
        /// unscaled value is determined by multiplying or dividing this BigDecimal's
        /// unscaled value by the appropriate power of ten to maintain its overall value.
        /// If the scale is reduced by the operation, the unscaled value must be divided
        /// (rather than multiplied), and the value may be changed; in this case, the
        /// specified rounding mode is applied to the division.
        /// <para>
        /// Note that since BigDecimal objects are immutable, calls of this method do
        /// not result in the original object being modified, contrary to the usual
        /// convention of having methods named setX mutate field X. Instead, setScale
        /// returns an object with the proper scale; the returned object may or may not
        /// be newly allocated.
        /// </para>
        /// </summary>
        /// <param name="Scale">scale of the BigDecimal value to be returned.</param>
        /// <param name="RoundingMethod">The rounding mode to apply.</param>
        /// <returns>a BigDecimal whose scale is the specified value, and whose unscaled
        /// value is determined by multiplying or dividing this BigDecimal's unscaled
        /// value by the appropriate power of ten to maintain its overall value.</returns>
        /// <exception cref="ArithmeticException">if <code>RoundingMethod==ROUND_UNNECESSARY
        /// </code> and the specified scaling operation would require rounding.</exception>
        /// <exception cref="ArgumentOutOfRangeException">if <see cref="RoungindMethod"/>
        /// does not represent a valid rounding mode.</exception>
        public BigDecimal SetScale(Int32 Scale, Byte RoundingMethod)
        {
            if (this.DecimalPart != 0)
            {
                Int32 DecimalDigits = (Int32)Math.Floor(BigInteger.Log10(this.DecimalPart) + 1);
                if (DecimalDigits > Scale)
                {
                    this.DecimalPart /= BigInteger.Pow(10, DecimalDigits - Scale - 1);
                    Byte LastDigit = (Byte)(this.DecimalPart % 10);
                    Boolean Increase = false;
                    switch (RoundingMethod)
                    {
                        case 0: //ROUND_UP
                            Increase = true;
                            break;
                        case 1: //ROUND_DOWN
                            Increase = false;
                            break;
                        case 2: //ROUND_CEILING
                            Increase = (this.IntegerPart > 0);
                            break;
                        case 3: //ROUND_FLOOR
                            Increase = (this.IntegerPart < 0);
                            break;
                        case 4: //ROUND_HALF-UP
                            Increase = (LastDigit >= 5);
                            break;
                        case 5: //ROUND_HALF-DOWN
                            Increase = (LastDigit > 5);
                            break;
                        case 6: //ROUND_HALF-EVEN
                            if (LastDigit > 5)
                            {
                                Increase = true;
                            }
                            if (LastDigit < 5)
                            {
                                Increase = false;
                            }
                            if(LastDigit == 5)
                            {
                                Byte SecondLastDigit = (Byte)((this.DecimalPart / 10) % 10);
                                Increase = (SecondLastDigit % 2 == 1);
                            }
                            break;
                        case 7: //ROUND_UNNECESSARY
                            throw new ArithmeticException("Rounding necessary");
                        default: throw new ArgumentOutOfRangeException("RoundingMethod", "Invalid rounding mode");
                    }
                    this.DecimalPart /= 10;
                    if (Increase)
                    {
                        ++this.DecimalPart;
                        if(DecimalDigits != (Int32)Math.Floor(BigInteger.Log10(this.DecimalPart) + 1))
                        {
                            this.DecimalPart = 0;
                            ++this.IntegerPart;
                        }
                    }
                }
            }
            return this;
        }

        /// <summary>
        /// Converts this BigDecimal to a double. This conversion is similar to the
        /// narrowing primitive conversion from double to float: if this BigDecimal
        /// has too great a magnitude represent as a double, it will be converted
        /// to <code>Double.NegativeInfinity</code> or <code>Double.PositiveInfinity</code>
        /// as appropriate.
        /// Note that even when the return value is finite,this conversion can lose
        /// information about the precision of the BigDecimal value.
        /// </summary>
        /// <returns>this BigDecimal converted to a double.</returns>
        public Double DoubleValue()
        {
            if (this > new BigDecimal(Double.MaxValue))
            {
                return Double.PositiveInfinity;
            }
            if (this < new BigDecimal(Double.MinValue))
            {
                return Double.NegativeInfinity;
            }
            Int32 DecimalDigits = (Int32)Math.Floor(BigInteger.Log10(this.DecimalPart) + 1);
            Double DecimalPart = (Double)this.DecimalPart / Math.Pow(10, DecimalDigits);
            return (Double)this.IntegerPart + DecimalPart;
        }
    }
}
