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
using System.Globalization;
using System.Text;

namespace Math3.complex
{
    /// <summary>
    /// Formats a Complex number in cartesian format "Re(c) + Im(c)i".  'i' can
    /// be replaced with 'j' (or anything else), and the number format for both real
    /// and imaginary parts can be configured.
    /// </summary>
    public class ComplexFormat
    {
        /// <summary>
        /// The default imaginary character.
        /// </summary>
        private const String DEFAULT_IMAGINARY_CHARACTER = "i";
        
        /// <summary>
        /// The notation used to signify the imaginary part of the complex number.
        /// </summary>
        private readonly String imaginaryCharacter;
        
        /// <summary>
        /// The format used for the imaginary part.
        /// </summary>
        private readonly NumberFormatInfo imaginaryFormat;
        
        /// <summary>
        /// The format used for the real part.
        /// </summary>
        private readonly NumberFormatInfo realFormat;

        /// <summary>
        /// Create an instance with the default imaginary character, 'i', and the
        /// default number format for both real and imaginary parts.
        /// </summary>
        public ComplexFormat()
        {
            this.imaginaryCharacter = DEFAULT_IMAGINARY_CHARACTER;
            this.imaginaryFormat = CompositeFormat.getDefaultNumberFormat();
            this.realFormat = imaginaryFormat;
        }

        /// <summary>
        /// Create an instance with a custom number format for both real and
        /// imaginary parts.
        /// </summary>
        /// <param name="format">the custom format for both real and imaginary parts.</param>
        /// <exception cref="NullArgumentException"> if <c>realFormat</c> is <c>null</c>.</exception>
        public ComplexFormat(NumberFormatInfo format)
        {
            if (format == null)
            {
                throw new NullArgumentException(new LocalizedFormats("IMAGINARY_FORMAT"));
            }
            this.imaginaryCharacter = DEFAULT_IMAGINARY_CHARACTER;
            this.imaginaryFormat = format;
            this.realFormat = format;
        }

        /// <summary>
        /// Create an instance with a custom number format for the real part and a
        /// custom number format for the imaginary part.
        /// </summary>
        /// <param name="realFormat">the custom format for the real part.</param>
        /// <param name="imaginaryFormat">the custom format for the imaginary part.</param>
        /// <exception cref="NullArgumentException"> if <c>imaginaryFormat</c> is <c>null</c>.</exception>
        /// <exception cref="NullArgumentException"> if <c>realFormat</c> is <c>null</c>.</exception>
        public ComplexFormat(NumberFormatInfo realFormat, NumberFormatInfo imaginaryFormat)
        {
            if (imaginaryFormat == null)
            {
                throw new NullArgumentException(new LocalizedFormats("IMAGINARY_FORMAT"));
            }
            if (realFormat == null)
            {
                throw new NullArgumentException(new LocalizedFormats("REAL_FORMAT"));
            }

            this.imaginaryCharacter = DEFAULT_IMAGINARY_CHARACTER;
            this.imaginaryFormat = imaginaryFormat;
            this.realFormat = realFormat;
        }

        /// <summary>
        /// Create an instance with a custom imaginary character, and the default
        /// number format for both real and imaginary parts.
        /// </summary>
        /// <param name="imaginaryCharacter">The custom imaginary character.</param>
        /// <exception cref="NullArgumentException"> if <c>imaginaryCharacter</c> is
        /// <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if <c>imaginaryCharacter</c> is an
        /// empty string.</exception>
        public ComplexFormat(String imaginaryCharacter) : this(imaginaryCharacter, CompositeFormat.getDefaultNumberFormat()) { }

        /// <summary>
        /// Create an instance with a custom imaginary character, and a custom number
        /// format for both real and imaginary parts.
        /// </summary>
        /// <param name="imaginaryCharacter">imaginaryCharacter The custom imaginary character.</param>
        /// <param name="format">the custom format for both real and imaginary parts.</param>
        /// <exception cref="NullArgumentException"> if <c>imaginaryCharacter</c> is
        /// <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if <c>imaginaryCharacter</c> is an
        /// empty string.</exception>
        /// <exception cref="NullArgumentException"> if <c>format</c> is <c>null</c>.</exception>
        public ComplexFormat(String imaginaryCharacter, NumberFormatInfo format) : this(imaginaryCharacter, format, format) { }

        /// <summary>
        /// Create an instance with a custom imaginary character, a custom number
        /// format for the real part, and a custom number format for the imaginary
        /// part.
        /// </summary>
        /// <param name="imaginaryCharacter">The custom imaginary character.</param>
        /// <param name="realFormat">the custom format for the real part.</param>
        /// <param name="imaginaryFormat">the custom format for the imaginary part.</param>
        /// <exception cref="NullArgumentException"> if <c>imaginaryCharacter</c> is
        /// <c>null</c>.</exception>
        /// <exception cref="NullArgumentException"> if <c>realFormat</c> is <c>null</c>.</exception>
        /// <exception cref="NullArgumentException"> if <c>imaginaryFormat</c> is <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if <c>imaginaryCharacter</c> is an
        /// empty string.</exception>
        public ComplexFormat(String imaginaryCharacter, NumberFormatInfo realFormat, NumberFormatInfo imaginaryFormat)
        {
            if (imaginaryCharacter == null)
            {
                throw new NullArgumentException();
            }
            if (imaginaryCharacter == String.Empty)
            {
                throw new NoDataException();
            }
            if (imaginaryFormat == null)
            {
                throw new NullArgumentException(new LocalizedFormats("IMAGINARY_FORMAT"));
            }
            if (realFormat == null)
            {
                throw new NullArgumentException(new LocalizedFormats("REAL_FORMAT"));
            }

            this.imaginaryCharacter = imaginaryCharacter;
            this.imaginaryFormat = imaginaryFormat;
            this.realFormat = realFormat;
        }

        /// <summary>
        /// This method calls <see cref="format(Object,StringBuilder)"/>.
        /// </summary>
        /// <param name="c">Complex object to format.</param>
        /// <returns>A formatted number in the form "Re(c) + Im(c)i".</returns>
        public String format(Complex c)
        {
            return format(c, new StringBuilder()).ToString();
        }

        /// <summary>
        /// This method calls <see cref="format(Object,StringBuilder)"/>.
        /// </summary>
        /// <param name="c">Double object to format.</param>
        /// <returns>A formatted number</returns>
        public String format(Double c)
        {
            return format(new Complex(c, 0), new StringBuilder()).ToString();
        }

        /// <summary>
        /// Formats a <see cref="Complex"/> object to produce a string.
        /// </summary>
        /// <param name="complex">the object to format.</param>
        /// <param name="toAppendTo">where the text is to be appended</param>
        /// <returns>the value passed in as toAppendTo.</returns>
        public StringBuilder format(Complex complex, StringBuilder toAppendTo)
        {
            CultureInfo Real = CultureInfo.CurrentCulture;
            Real.NumberFormat = getRealFormat();
            // format real
            double re = complex.getReal();
            CompositeFormat.formatDouble(re, Real, toAppendTo);

            // format sign and imaginary
            double im = complex.getImaginary();
            StringBuilder imAppendTo;
            if (im < 0.0)
            {
                toAppendTo.Append(" - ");
                imAppendTo = formatImaginary(-im, new StringBuilder());
                toAppendTo.Append(imAppendTo);
                toAppendTo.Append(getImaginaryCharacter());
            }
            else if (im > 0.0 || Double.IsNaN(im))
            {
                toAppendTo.Append(" + ");
                imAppendTo = formatImaginary(im, new StringBuilder());
                toAppendTo.Append(imAppendTo);
                toAppendTo.Append(getImaginaryCharacter());
            }
            return toAppendTo;
        }

        /// <summary>
        /// Format the absolute value of the imaginary part.
        /// </summary>
        /// <param name="absIm">Absolute value of the imaginary part of a complex number.</param>
        /// <param name="toAppendTo">where the text is to be appended.</param>
        /// <returns>the value passed in as toAppendTo.</returns>
        private StringBuilder formatImaginary(double absIm, StringBuilder toAppendTo)
        {
            CultureInfo Imag = CultureInfo.CurrentCulture;
            Imag.NumberFormat = getRealFormat();
            CompositeFormat.formatDouble(absIm, Imag, toAppendTo);
            if (toAppendTo.ToString().Equals("1"))
            {
                // Remove the character "1" if it is the only one.
                toAppendTo.Length = 0;
            }
            return toAppendTo;
        }

        /// <summary>
        /// Formats a object to produce a string.  <c>obj</c> must be either a
        /// <see cref="Complex"/> object or a Number. Any other type of
        /// object will result in an <see cref="MathIllegalArgumentException"/> being thrown.
        /// </summary>
        /// <param name="obj">the object to format.</param>
        /// <param name="toAppendTo">where the text is to be appended</param>
        /// <returns>the value passed in as toAppendTo.</returns>
        /// <exception cref="MathIllegalArgumentException"> is <c>obj</c> is not a valid type</exception>
        public StringBuilder format(Object obj, StringBuilder toAppendTo)
        {
            StringBuilder ret = null;
            if (obj is Complex)
            {
                ret = format((Complex)obj, toAppendTo);
            }
            else if (obj is Byte || obj is SByte || obj is Int16 || obj is UInt16
                || obj is Int32 || obj is UInt32 || obj is Int64 || obj is UInt64
                || obj is Single || obj is Double || obj is Decimal)
            {
                ret = format(new Complex((Double)obj, 0.0), toAppendTo);
            }
            else
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("CANNOT_FORMAT_INSTANCE_AS_COMPLEX"), obj.GetType().Name);
            }

            return ret;
        }

        /// <summary>
        /// Access the imaginaryCharacter.
        /// </summary>
        /// <returns>the imaginaryCharacter.</returns>
        public String getImaginaryCharacter()
        {
            return imaginaryCharacter;
        }

        /// <summary>
        /// Access the imaginaryFormat.
        /// </summary>
        /// <returns>the imaginaryFormat.</returns>
        public NumberFormatInfo getImaginaryFormat()
        {
            return imaginaryFormat;
        }

        /// <summary>
        /// Returns the default complex format for the current locale.
        /// </summary>
        /// <returns>the default complex format.</returns>
        public static ComplexFormat getInstance()
        {
            return getInstance(CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Returns the default complex format for the given locale.
        /// </summary>
        /// <param name="locale">the specific locale used by the format.</param>
        /// <returns>the complex format specific to the given locale.</returns>
        public static ComplexFormat getInstance(CultureInfo locale)
        {
            NumberFormatInfo f = CompositeFormat.getDefaultNumberFormat(locale);
            return new ComplexFormat(f);
        }

        /// <summary>
        /// Returns the default complex format for the given locale.
        /// </summary>
        /// <param name="imaginaryCharacter">Imaginary character.</param>
        /// <param name="locale">the specific locale used by the format.</param>
        /// <returns>the complex format specific to the given locale.</returns>
        /// <exception cref="NullArgumentException"> if <c>imaginaryCharacter</c> is
        /// <c>null</c>.</exception>
        /// <exception cref="NoDataException"> if <c>imaginaryCharacter</c> is an
        /// empty string.</exception>
        public static ComplexFormat getInstance(String imaginaryCharacter, CultureInfo locale)
        {
            NumberFormatInfo f = CompositeFormat.getDefaultNumberFormat(locale);
            return new ComplexFormat(imaginaryCharacter, f);
        }

        /// <summary>
        /// Access the realFormat.
        /// </summary>
        /// <returns>the realFormat.</returns>
        public NumberFormatInfo getRealFormat()
        {
            return realFormat;
        }

        /// <summary>
        /// Parses a string to produce a <c>Complex</c> object.
        /// </summary>
        /// <param name="source">source the string to parse.</param>
        /// <returns>the parsed <c>Complex</c> object.</returns>
        /// <exception cref="MathParseException"> if the beginning of the specified string
        /// cannot be parsed.</exception>
        public Complex parse(String source)
        {
            Int32 parsePosition = 0;
            Complex result = parse(source, parsePosition);
            if (result == null)
            {
                throw new MathParseException(source, parsePosition, typeof(Complex));
            }
            return result;
        }

        /// <summary>
        /// Parses a string to produce a <see cref="Complex"/> object. 
        /// </summary>
        /// <param name="source">the string to parse</param>
        /// <param name="pos">input/ouput parsing parameter.</param>
        /// <returns>the parsed <c>Complex</c> object.</returns>
        public Complex parse(String source, Int32 pos)
        {
            int initialIndex = pos;

            // parse whitespace
            Int32 nonSpace = CompositeFormat.parseAndIgnoreWhitespace(source, pos);

            // parse real
            Double re = CompositeFormat.parseNumber(source, CultureInfo.CurrentCulture);
            if (Double.IsNaN(re))
            {
                // invalid real number
                // set index back to initial, error index should already be set
                return null;
            }

            // parse sign
            Int32 p = CompositeFormat.parseNextCharacter(source, nonSpace);
            char c = source[p];
            int sign = 0;
            switch (c)
            {
                case '0':
                    // no sign
                    // return real only complex number
                    return new Complex(re, 0.0);
                case '-':
                    sign = -1;
                    break;
                case '+':
                    sign = 1;
                    break;
                default:
                    // invalid sign
                    // set index back to initial, error index should be the last
                    // character examined.
                    return null;
            }

            // parse whitespace
            CompositeFormat.parseAndIgnoreWhitespace(source, pos);

            // parse imaginary
            Double im = CompositeFormat.parseNumber(source, CultureInfo.CurrentCulture);
            if (Double.IsNaN(im))
            {
                // invalid imaginary number
                // set index back to initial, error index should already be set
                return null;
            }

            // parse imaginary character
            if (!CompositeFormat.parseFixedstring(source, getImaginaryCharacter()))
            {
                return null;
            }
            return new Complex(re, im * sign);
        }
    }
}
