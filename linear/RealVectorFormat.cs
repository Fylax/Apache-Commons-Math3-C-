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
using Math3.util;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace Math3.linear
{
    /// <summary>
    /// Formats a vector in components list format "{v0; v1; ...; vk-1}".
    /// <para>The prefix and suffix "{" and "}" and the separator "; " can be replaced by
    /// any user-defined strings. The number format for components can be configured.</para>
    /// <para>White space is ignored at parse time, even if it is in the prefix, suffix
    /// or separator specifications. So even if the default separator does include a space
    /// character that is used at format time, both input string "{1;1;1}" and
    /// " { 1 ; 1 ; 1 } " will be parsed without error and the same vector will be
    /// returned. In the second case, however, the parse position after parsing will be
    /// just after the closing curly brace, i.e. just before the trailing space.</par
    /// </summary>
    public class RealVectorFormat
    {
        /// <summary>
        /// The default prefix: "{".
        /// </summary>
        private const String DEFAULT_PREFIX = "{";

        /// <summary>
        /// The default suffix: "}".
        /// </summary>
        private const String DEFAULT_SUFFIX = "}";

        /// <summary>
        /// The default separator: ", ".
        /// </summary>
        private const String DEFAULT_SEPARATOR = "; ";

        /// <summary>
        /// Prefix.
        /// </summary>
        private readonly String prefix;

        /// <summary>
        /// Suffix.
        /// </summary>
        private readonly String suffix;

        /// <summary>
        /// Separator.
        /// </summary>
        private readonly String separator;

        /// <summary>
        /// Trimmed prefix.
        /// </summary>
        private readonly String trimmedPrefix;

        /// <summary>
        /// Trimmed suffix.
        /// </summary>
        private readonly String trimmedSuffix;

        /// <summary>
        /// Trimmed separator.
        /// </summary>
        private readonly String trimmedSeparator;

        /// <summary>
        /// The format used for components.
        /// </summary>
        private readonly NumberFormatInfo Format;

        /// <summary>
        /// Create an instance with default settings.
        /// <para>The instance uses the default prefix, suffix and separator:
        /// "{", "}", and "; " and the default number format for components.</para>
        /// </summary>
        public RealVectorFormat() : this(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, CompositeFormat.getDefaultNumberFormat()) { }

        /// <summary>
        /// Create an instance with a custom number format for components.
        /// </summary>
        /// <param name="format">the custom format for components.</param>
        public RealVectorFormat(NumberFormatInfo format) : this(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, format) { }

        /// <summary>
        /// Create an instance with custom prefix, suffix and separator.
        /// </summary>
        /// <param name="prefix">prefix to use instead of the default "{"</param>
        /// <param name="suffix">suffix to use instead of the default "}"</param>
        /// <param name="separator">separator to use instead of the default "; "</param>
        public RealVectorFormat(String prefix, String suffix, String separator) : this(prefix, suffix, separator, CompositeFormat.getDefaultNumberFormat()) { }

        /// <summary>
        /// Create an instance with custom prefix, suffix, separator and format
        /// for components.
        /// </summary>
        /// <param name="prefix">prefix to use instead of the default "{"</param>
        /// <param name="suffix">suffix to use instead of the default "}"</param>
        /// <param name="separator">separator to use instead of the default "; "</param>
        /// <param name="format">the custom format for components.</param>
        public RealVectorFormat(String prefix, String suffix, String separator, NumberFormatInfo format)
        {
            this.prefix = prefix;
            this.suffix = suffix;
            this.separator = separator;
            trimmedPrefix = prefix.Trim();
            trimmedSuffix = suffix.Trim();
            trimmedSeparator = separator.Trim();
            this.Format = format;
        }

        /// <summary>
        /// Get the set of locales for which real vectors formats are available.
        /// <para>This is the same set as the <see cref="NumberFormatInfo"/> set.</para>
        /// </summary>
        /// <returns>available real vector format locales.</returns>
        public static CultureInfo[] getAvailableLocales()
        {
            // cannot fetch all cultures for portability :(
            CultureInfo[] cul = new CultureInfo[1];
            cul[0] = CultureInfo.CurrentCulture;
            return cul;
        }

        /// <summary>
        /// Get the format prefix.
        /// </summary>
        /// <returns>format prefix.</returns>
        public String getPrefix()
        {
            return prefix;
        }

        /// <summary>
        /// Get the format suffix.
        /// </summary>
        /// <returns>format suffix.</returns>
        public String getSuffix()
        {
            return suffix;
        }

        /// <summary>
        /// Get the format separator between components.
        /// </summary>
        /// <returns>format separator.</returns>
        public String getSeparator()
        {
            return separator;
        }

        /// <summary>
        /// Get the components format.
        /// </summary>
        /// <returns>components format.</returns>
        public NumberFormatInfo getFormat()
        {
            return Format;
        }

        /// <summary>
        /// Returns the default real vector format for the current locale.
        /// </summary>
        /// <returns>the default real vector format.</returns>
        public static RealVectorFormat getInstance()
        {
            return getInstance(CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Returns the default real vector format for the given locale.
        /// </summary>
        /// <param name="locale">the specific locale used by the format.</param>
        /// <returns>the real vector format specific to the given locale.</returns>
        public static RealVectorFormat getInstance(CultureInfo locale)
        {
            return new RealVectorFormat(CompositeFormat.getDefaultNumberFormat(locale));
        }

        /// <summary>
        /// This method calls <see cref="#format(RealVector,StringBuffer,FieldPosition)"/>.
        /// </summary>
        /// <param name="v">RealVector object to format.</param>
        /// <returns>a formatted vector.</returns>
        public String format(RealVector v)
        {
            return format(v, new StringBuilder()).ToString();
        }

        /// <summary>
        /// Formats a <see cref="RealVector"/> object to produce a string.
        /// </summary>
        /// <param name="vector">the object to format.</param>
        /// <param name="toAppendTo">where the text is to be appended</param>
        /// <returns>the value passed in as toAppendTo.</returns>
        public StringBuilder format(RealVector vector, StringBuilder toAppendTo)
        {

            // format prefix
            toAppendTo.Append(prefix);

            // format components
            for (int i = 0; i < vector.getDimension(); ++i)
            {
                if (i > 0)
                {
                    toAppendTo.Append(separator);
                }
                CompositeFormat.formatDouble(vector.getEntry(i), CultureInfo.CurrentCulture, toAppendTo);
            }

            // format suffix
            toAppendTo.Append(suffix);

            return toAppendTo;
        }

        /// <summary>
        /// Parse a string to produce a <see cref="RealVector"/> object.
        /// </summary>
        /// <param name="source">String to parse.</param>
        /// <returns>the parsed <see cref="RealVector"/> object.</returns>
        /// <exception cref="MathParseException"> if the beginning of the specified string
        /// cannot be parsed.</exception>
        public ArrayRealVector parse(String source)
        {
            Int32 parsePosition = 0;
            ArrayRealVector result = parse(source, parsePosition);
            if (parsePosition == 0)
            {
                throw new MathParseException(source, parsePosition, typeof(ArrayRealVector));
            }
            return result;
        }

        /// <summary>
        /// Parse a string to produce a <see cref="RealVector"/> object.
        /// </summary>
        /// <param name="source">String to parse.</param>
        /// <param name="pos">input/ouput parsing parameter.</param>
        /// <returns>the parsed <see cref="RealVector"/> object.</returns>
        public ArrayRealVector parse(String source, Int32 pos)
        {
            int initialIndex = pos;

            // parse prefix
            CompositeFormat.parseAndIgnoreWhitespace(source, pos);
            if (!CompositeFormat.parseFixedstring(source, trimmedPrefix))
            {
                return null;
            }

            // parse components
            List<double> components = new List<double>();
            for (Boolean loop = true; loop; )
            {

                if (!(components.Count == 0))
                {
                    CompositeFormat.parseAndIgnoreWhitespace(source, pos);
                    if (!CompositeFormat.parseFixedstring(source, trimmedSeparator))
                    {
                        loop = false;
                    }
                }

                if (loop)
                {
                    CompositeFormat.parseAndIgnoreWhitespace(source, pos);
                    double component = CompositeFormat.parseNumber(source, CultureInfo.CurrentCulture);
                    if (Double.IsNaN(component))
                    {
                        components.Add(component);
                    }
                    else
                    {
                        // invalid component
                        // set index back to initial, error index should already be set
                        pos = initialIndex;
                        return null;
                    }
                }

            }

            // parse suffix
            CompositeFormat.parseAndIgnoreWhitespace(source, pos);
            if (!CompositeFormat.parseFixedstring(source, trimmedSuffix))
            {
                return null;
            }

            // build vector
            double[] data = new double[components.Count];
            for (int i = 0; i < data.Length; ++i)
            {
                data[i] = components[i];
            }
            return new ArrayRealVector(data, false);
        }
    }
}
