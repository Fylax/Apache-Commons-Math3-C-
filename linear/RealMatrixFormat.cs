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
    /// Formats a <c>nxm</c> matrix in components list format
    /// "{{a_0_0,a_0_1, ...,
    /// a_0_{m-1}},{a_1_0},
    /// a_1_1, ..., a_1_{m-1}},{...},{
    /// a_{n-1}_0, a_{n-1}_1, ...,
    /// a_{n-1}_{m-1}}}".
    /// <para>The prefix and suffix "{" and "}", the row prefix and suffix "{" and "}",
    /// the row separator "," and the column separator "," can be replaced by any
    /// user-defined strings. The number format for components can be configured.</para>
    /// <para>White space is ignored at parse time, even if it is in the prefix, suffix
    /// or separator specifications. So even if the default separator does include a space
    /// character that is used at format time, both input string "{{1,1,1}}" and
    /// " { { 1 , 1 , 1 } } " will be parsed without error and the same matrix will be
    /// returned. In the second case, however, the parse position after parsing will be
    /// just after the closing curly brace, i.e. just before the trailing space.</para>
    /// <para>Note: the grouping functionality of the used <see cref="NumberFormatInfo"/> is
    /// disabled to prevent problems when parsing (e.g. 1,345.34 would be a valid number
    /// but conflicts with the default column separator).</para>
    /// </summary>
    public class RealMatrixFormat
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
        /// The default row prefix: "{".
        /// </summary>
        private const String DEFAULT_ROW_PREFIX = "{";

        /// <summary>
        /// The default row suffix: "}".
        /// </summary>
        private const String DEFAULT_ROW_SUFFIX = "}";

        /// <summary>
        /// The default row separator: ",".
        /// </summary>
        private const String DEFAULT_ROW_SEPARATOR = ",";

        /// <summary>
        /// The default column separator: ",".
        /// </summary>
        private const String DEFAULT_COLUMN_SEPARATOR = ",";

        /// <summary>
        /// Prefix.
        /// </summary>
        private readonly String prefix;

        /// <summary>
        /// Suffix.
        /// </summary>
        private readonly String suffix;

        /// <summary>
        /// Row prefix.
        /// </summary>
        private readonly String rowPrefix;

        /// <summary>
        /// Row suffix.
        /// </summary>
        private readonly String rowSuffix;

        /// <summary>
        /// Row separator.
        /// </summary>
        private readonly String rowSeparator;

        /// <summary>
        /// Column separator.
        /// </summary>
        private readonly String columnSeparator;

        /// <summary>
        /// The format used for components.
        /// </summary>
        private readonly NumberFormatInfo Format;

        /// <summary>
        /// Create an instance with default settings.
        /// <para>The instance uses the default prefix, suffix and row/column separator:
        /// "[", "]", ";" and ", " and the default number format for components.</para>
        /// </summary>
        public RealMatrixFormat() : this(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_ROW_PREFIX, DEFAULT_ROW_SUFFIX, DEFAULT_ROW_SEPARATOR, DEFAULT_COLUMN_SEPARATOR, CompositeFormat.getDefaultNumberFormat()) { }

        /// <summary>
        /// Create an instance with a custom number format for components.
        /// </summary>
        /// <param name="format">the custom format for components.</param>
        public RealMatrixFormat(NumberFormatInfo format) : this(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_ROW_PREFIX, DEFAULT_ROW_SUFFIX, DEFAULT_ROW_SEPARATOR, DEFAULT_COLUMN_SEPARATOR, format) { }

        /// <summary>
        /// Create an instance with custom prefix, suffix and separator.
        /// </summary>
        /// <param name="prefix">prefix to use instead of the default "{"</param>
        /// <param name="suffix">suffix to use instead of the default "}"</param>
        /// <param name="rowPrefix">row prefix to use instead of the default "{"</param>
        /// <param name="rowSuffix">row suffix to use instead of the default "}"</param>
        /// <param name="rowSeparator">tow separator to use instead of the default ";"</param>
        /// <param name="columnSeparator">column separator to use instead of the default ", "</param>
        public RealMatrixFormat(String prefix, String suffix, String rowPrefix, String rowSuffix, String rowSeparator, String columnSeparator) : this(prefix, suffix, rowPrefix, rowSuffix, rowSeparator, columnSeparator, CompositeFormat.getDefaultNumberFormat()) { }

        /// <summary>
        /// Create an instance with custom prefix, suffix, separator and format
        /// for components.
        /// </summary>
        /// <param name="prefix">prefix to use instead of the default "{"</param>
        /// <param name="suffix">suffix to use instead of the default "}"</param>
        /// <param name="rowPrefix">row prefix to use instead of the default "{"</param>
        /// <param name="rowSuffix">row suffix to use instead of the default "}"</param>
        /// <param name="rowSeparator">tow separator to use instead of the default ";"</param>
        /// <param name="columnSeparator">column separator to use instead of the default ", "
        /// </param>
        /// <param name="format">the custom format for components.</param>
        public RealMatrixFormat(String prefix, String suffix, String rowPrefix, String rowSuffix, String rowSeparator, String columnSeparator, NumberFormatInfo format)
        {
            this.prefix = prefix;
            this.suffix = suffix;
            this.rowPrefix = rowPrefix;
            this.rowSuffix = rowSuffix;
            this.rowSeparator = rowSeparator;
            this.columnSeparator = columnSeparator;
            this.Format = format;
            // disable grouping to prevent parsing problems
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
        /// Get the format prefix.
        /// </summary>
        /// <returns>format prefix.</returns>
        public String getRowPrefix()
        {
            return rowPrefix;
        }

        /// <summary>
        /// Get the format suffix.
        /// </summary>
        /// <returns>format suffix.</returns>
        public String getRowSuffix()
        {
            return rowSuffix;
        }

        /// <summary>
        /// Get the format separator between rows of the matrix.
        /// </summary>
        /// <returns>format separator for rows.</returns>
        public String getRowSeparator()
        {
            return rowSeparator;
        }

        /// <summary>
        /// Get the format separator between components.
        /// </summary>
        /// <returns>format separator between components.</returns>
        public String getColumnSeparator()
        {
            return columnSeparator;
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
        public static RealMatrixFormat getInstance()
        {
            return getInstance(CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Returns the default real vector format for the given locale.
        /// </summary>
        /// <param name="locale">the specific locale used by the format.</param>
        /// <returns>the real vector format specific to the given locale.</returns>
        public static RealMatrixFormat getInstance(CultureInfo locale)
        {
            return new RealMatrixFormat(CompositeFormat.getDefaultNumberFormat(locale));
        }

        /// <summary>
        /// This method calls <see cref="#format(RealMatrix,StringBuffer,FieldPosition)"/>.
        /// </summary>
        /// <param name="m">RealMatrix object to format.</param>
        /// <returns>a formatted matrix.</returns>
        public String format(RealMatrix m)
        {
            return format(m, new StringBuilder()).ToString();
        }

        /// <summary>
        /// Formats a <see cref="RealMatrix"/> object to produce a string.
        /// </summary>
        /// <param name="matrix">the object to format.</param>
        /// <param name="toAppendTo">where the text is to be appended</param>
        /// <returns>the value passed in as toAppendTo.</returns>
        public StringBuilder format(RealMatrix matrix, StringBuilder toAppendTo)
        {

            // format prefix
            toAppendTo.Append(prefix);

            // format rows
            int rows = matrix.getRowDimension();
            for (int i = 0; i < rows; ++i)
            {
                toAppendTo.Append(rowPrefix);
                for (int j = 0; j < matrix.getColumnDimension(); ++j)
                {
                    if (j > 0)
                    {
                        toAppendTo.Append(columnSeparator);
                    }
                    CompositeFormat.formatDouble(matrix.getEntry(i, j), CultureInfo.CurrentCulture, toAppendTo);
                }
                toAppendTo.Append(rowSuffix);
                if (i < rows - 1)
                {
                    toAppendTo.Append(rowSeparator);
                }
            }

            // format suffix
            toAppendTo.Append(suffix);

            return toAppendTo;
        }

        /// <summary>
        /// Parse a string to produce a <see cref="RealMatrix"/> object.
        /// </summary>
        /// <param name="source">String to parse.</param>
        /// <returns>the parsed <see cref="RealMatrix"/> object.</returns>
        /// <exception cref="MathParseException"> if the beginning of the specified string
        /// cannot be parsed.</exception>
        public RealMatrix parse(String source)
        {
            Int32 parsePosition = 0;
            RealMatrix result = parse(source, ref parsePosition);
            if (parsePosition == 0)
            {
                throw new MathParseException(source, parsePosition, typeof(Array2DRowRealMatrix));
            }
            return result;
        }

        /// <summary>
        /// Parse a string to produce a <see cref="RealMatrix"/> object.
        /// </summary>
        /// <param name="source">String to parse.</param>
        /// <param name="pos">input/ouput parsing parameter.</param>
        /// <returns>the parsed <see cref="RealMatrix"/> object.</returns>
        public RealMatrix parse(String source, ref Int32 pos)
        {
            int initialIndex = pos;

            String trimmedPrefix = prefix.Trim();
            String trimmedSuffix = suffix.Trim();
            String trimmedRowPrefix = rowPrefix.Trim();
            String trimmedRowSuffix = rowSuffix.Trim();
            String trimmedColumnSeparator = columnSeparator.Trim();
            String trimmedRowSeparator = rowSeparator.Trim();

            // parse prefix
            CompositeFormat.parseAndIgnoreWhitespace(source, pos);
            if (!CompositeFormat.parseFixedstring(source, trimmedPrefix))
            {
                return null;
            }


            // parse components
            List<List<double>> matrix = new List<List<double>>();
            List<double> rowComponents = new List<double>();
            for (Boolean loop = true; loop; )
            {
                if (rowComponents.Count != 0)
                {
                    CompositeFormat.parseAndIgnoreWhitespace(source, pos);
                    if (!CompositeFormat.parseFixedstring(source, trimmedColumnSeparator))
                    {
                        if (trimmedRowSuffix.Length != 0 &&
                            !CompositeFormat.parseFixedstring(source, trimmedRowSuffix))
                        {
                            return null;
                        }
                        else
                        {
                            CompositeFormat.parseAndIgnoreWhitespace(source, pos);
                            if (CompositeFormat.parseFixedstring(source, trimmedRowSeparator))
                            {
                                matrix.Add(rowComponents);
                                rowComponents = new List<double>();
                                continue;
                            }
                            else
                            {
                                loop = false;
                            }
                        }
                    }
                }
                else
                {
                    CompositeFormat.parseAndIgnoreWhitespace(source, pos);
                    if (trimmedRowPrefix.Length != 0 &&
                        !CompositeFormat.parseFixedstring(source, trimmedRowPrefix))
                    {
                        return null;
                    }
                }

                if (loop)
                {
                    CompositeFormat.parseAndIgnoreWhitespace(source, pos);
                    double component = CompositeFormat.parseNumber(source, CultureInfo.CurrentCulture);
                    if (Double.IsNaN(component))
                    {
                        rowComponents.Add(component);
                    }
                    else
                    {
                        if (rowComponents.Count == 0)
                        {
                            loop = false;
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

            }

            if (rowComponents.Count != 0)
            {
                matrix.Add(rowComponents);
            }

            // parse suffix
            CompositeFormat.parseAndIgnoreWhitespace(source, pos);
            if (!CompositeFormat.parseFixedstring(source, trimmedSuffix))
            {
                return null;
            }

            // do not allow an empty matrix
            if (matrix.Count == 0)
            {
                pos = initialIndex;
                return null;
            }

            // build vector
            double[][] data = new double[matrix.Count][];
            int row = 0;
            foreach (List<double> rowList in matrix)
            {
                data[row] = new double[rowList.Count];
                for (int i = 0; i < rowList.Count; i++)
                {
                    data[row][i] = rowList[i];
                }
                row++;
            }
            return MatrixUtils.createRealMatrix(data);
        }
    }
}
