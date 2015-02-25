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
using System.Text;

namespace Math3.util
{
    /// <summary>
    /// Base class for formatters of composite objects (complex numbers, vectors ...).
    /// </summary>
    public class CompositeFormat
    {

        /// <summary>
        /// Class contains only static methods.
        /// </summary>
        private CompositeFormat() { }

        /// <summary>
        /// Create a default number format.  The default number format is based on
        /// <see cref="NumberFormatInfo"/> with the only customizing that the
        /// maximum number of decimal digits is set to 10.
        /// </summary>
        /// <returns>the default number format.</returns>
        public static NumberFormatInfo getDefaultNumberFormat()
        {
            return getDefaultNumberFormat(CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Create a default number format.  The default number format is based on
        /// <see cref="CultureInfo.NumberFormat"/> with the only
        /// customizing that the maximum number of fraction digits is set to 10. 
        /// </summary>
        /// <param name="locale">the specific locale used by the format.</param>
        /// <returns>the default number format specific to the given locale.</returns>
        public static NumberFormatInfo getDefaultNumberFormat(CultureInfo locale)
        {
            NumberFormatInfo nf = locale.NumberFormat;
            nf.NumberDecimalDigits = 10;
            return nf;
        }

        /// <summary>
        /// Parses <c>source</c> until a non-whitespace character is found.
        /// </summary>
        /// <param name="source">the string to parse</param>
        /// <param name="pos">parsing parameter.</param>
        /// <returns>holds the index of the next non-whitespace character.</returns>
        public static Int32 parseAndIgnoreWhitespace(String source, Int32 pos)
        {
            return parseNextCharacter(source, pos);
        }

        /// <summary>
        /// Parses <c>source</c> until a non-whitespace character is found.
        /// </summary>
        /// <param name="source">the string to parse</param>
        /// <param name="pos">parsing parameter.</param>
        /// <returns>the position of the first non-whitespace character.</returns>
        public static Int32 parseNextCharacter(String source, Int32 pos)
        {
            int index = pos;
            int n = source.Length;
            if (index < n)
            {
                char c;
                do
                {
                    c = source[index++];
                } while (Char.IsWhiteSpace(c) && index < n);
                pos = index;
                if (index < n)
                {
                    return index;
                }
            }
            return n;
        }

        /// <summary>
        /// Parses <c>source</c> for special double values.  These values
        /// include Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity.
        /// </summary>
        /// <param name="source">the string to parse</param>
        /// <param name="value">the special value to parse.</param>
        /// <param name="pos">parsing parameter.</param>
        /// <returns>the special number.</returns>
        private static Double parseNumber(String source, double value, Int32 pos)
        {
            Double ret = default(Double);

            StringBuilder sb = new StringBuilder();
            sb.Append('(');
            sb.Append(value);
            sb.Append(')');

            int n = sb.Length;
            int startIndex = pos;
            int endIndex = startIndex + n;
            if (endIndex < source.Length &&
                source.Substring(startIndex, endIndex).CompareTo(sb.ToString()) == 0)
            {
                ret = value;
                pos = endIndex;
            }
            return ret;
        }

        /// <summary>
        /// Parses <c>source</c> for a number.  This method can parse normal,
        /// numeric values as well as special values.  These special values include
        /// Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity.
        /// </summary>
        /// <param name="source">the string to parse</param>
        /// <param name="format"><c>CultureInfo</c> for getting the number format used
        /// to parse normal, numeric values.</param>
        /// <returns>the parsed number.</returns>
        public static Double parseNumber(String source, CultureInfo format)
        {
            try
            {
                return Double.Parse(source, format.NumberFormat);
            }
            catch (FormatException)
            {
                return Double.NaN; //null is not accettable
            }
        }

        /// <summary>
        /// Parse <code>source</code> for an expected fixed string. 
        /// </summary>
        /// <param name="source">the string to parse</param>
        /// <param name="expected">expected string</param>
        /// <returns>true if the expected string was there</returns>
        public static Boolean parseFixedstring(String source, String expected)
        {
            return source.Contains(expected);
        }

        /// <summary>
        /// Formats a double value to produce a string.  In general, the value is
        /// formatted using the formatting rules of <c>format</c>.  There are
        /// three exceptions to this:
        /// <list type="number">
        /// <item>NaN is formatted as '(NaN)'</item>
        /// <item>Positive infinity is formatted as '(Infinity)'</item>
        /// <item>Negative infinity is formatted as '(-Infinity)'</item>
        /// </list> 
        /// </summary>
        /// <param name="value">the double to format.</param>
        /// <param name="format"><c>CultureInfo</c> to get the format used.</param>
        /// <param name="toAppendTo">where the text is to be appended</param>
        /// <returns>the value passed in as toAppendTo.</returns>
        public static StringBuilder formatDouble(double value, CultureInfo format, StringBuilder toAppendTo)
        {
            if (Double.IsNaN(value) || Double.IsInfinity(value))
            {
                toAppendTo.Append('(');
                toAppendTo.Append(value);
                toAppendTo.Append(')');
            }
            else
            {
                toAppendTo.Append(value.ToString(format.NumberFormat));
            }
            return toAppendTo;
        }
    }

}
