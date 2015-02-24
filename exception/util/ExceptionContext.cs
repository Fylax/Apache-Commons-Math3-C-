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
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace Math3.exception.util
{
    /// <summary>
    /// Class that contains the actual implementation of the functionality mandated
    /// by the <see cref="ExceptionContext"/> interface. All Commons Math exceptions
    /// delegate the interface's methods to this class.
    /// </summary>
    public class ExceptionContext 
    {
        /// <summary>
        /// The throwable to which this context refers to.
        /// </summary>
        private Exception throwable;

        /// <summary>
        /// Various informations that enrich the informative message.
        /// </summary>
        private List<Localizable> msgPatterns;

        /// <summary>
        /// Various informations that enrich the informative message.
        /// The arguments will replace the corresponding place-holders in
        /// <see cref="msgPatterns"/>.
        /// </summary>
        private List<Object[]> msgArguments;

        /// <summary>
        /// Arbitrary context information.
        /// </summary>
        private Dictionary<String, Object> context;

        /// <summary>
        /// Simple constructor.
        /// </summary>
        /// <param name="throwable">the exception this context refers too</param>
        public ExceptionContext(Exception throwable) {
            this.throwable = throwable;
            this.msgPatterns = new List<Localizable>();
            this.msgArguments = new List<Object[]>();
            this.context = new Dictionary<String, Object>();
        }

        /// <summary>
        /// Get a reference to the exception to which the context relates.
        /// </summary>
        /// <returns>a reference to the exception to which the context relates</returns>
        public Exception getThrowable()
        {
            return throwable;
        }

        /// <summary>
        /// Adds a message.
        /// </summary>
        /// <param name="pattern">Message pattern.</param>
        /// <param name="arguments">Values for replacing the placeholders in the message
        /// pattern.</param>
        public void addMessage(Localizable pattern, params Object[] arguments) {
            this.msgPatterns.Add(pattern);
            this.msgArguments.Add(ArgUtils.flatten(arguments));
        }

        /// <summary>
        /// Sets the context (key, value) pair.
        /// Keys are assumed to be unique within an instance. If the same key is
        /// assigned a new value, the previous one will be lost.
        /// </summary>
        /// <param name="key">Context key (not null).</param>
        /// <param name="value">Context value.</param>
        public void setValue(String key, Object value)
        {
            this.context.Add(key, value);
        }

        /// <summary>
        /// Gets the value associated to the given context key.
        /// </summary>
        /// <param name="key">Context key.</param>
        /// <returns>the context value or <code>null</code> if the key does not exist.</returns>
        public Object getValue(String key)
        {
            Object value;
            this.context.TryGetValue(key, out value);
            return value;
        }

        /// <summary>
        /// Gets all the keys stored in the exception
        /// </summary>
        /// <returns>the set of keys.</returns>
        public HashSet<String> getKeys()
        {
            //All this craziness just to maintain the Java Set, converted to HashSet.
            String[] KeysArray = new String[this.context.Keys.Count];
            this.context.Keys.CopyTo(KeysArray, 0);
            return new HashSet<String>(KeysArray);
        }

        /// <summary>
        /// Gets the default message.
        /// </summary>
        /// <returns>the message.</returns>
        public String getMessage()
        {
            return getMessage(new CultureInfo("en-US"));
        }

        /// <summary>
        /// Gets the message in the default locale.
        /// </summary>
        /// <returns>the localized message.</returns>
        public String getLocalizedMessage()
        {
            return getMessage(CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Gets the message in a specified locale.
        /// </summary>
        /// <param name="locale">Locale in which the message should be translated.</param>
        /// <returns>the localized message.</returns>
        public String getMessage(CultureInfo locale) {
            return buildMessage(locale, ": ");
        }

        /// <summary>
        /// Gets the message in a specified locale.
        /// </summary>
        /// <param name="locale">Locale in which the message should be translated.</param>
        /// <param name="separator">Separator inserted between the message parts.</param>
        /// <returns>the localized message.</returns>
        public String getMessage(CultureInfo locale, String separator) {
            return buildMessage(locale, separator);
        }

        /// <summary>
        /// Builds a message string.
        /// </summary>
        /// <param name="locale">Locale in which the message should be translated.</param>
        /// <param name="separator">Message separator.</param>
        /// <returns>a localized message string.</returns>
        private String buildMessage(CultureInfo locale, String separator)
        {
            StringBuilder sb = new StringBuilder();
            int count = 0;
            int len = msgPatterns.Count;
            for (int i = 0; i < len; i++)
            {
                Localizable pat = msgPatterns[i];
                Object[] args = msgArguments[i];
                //MessageFormat fmt = new MessageFormat(pat.getLocalizedString(locale), locale);
                //The line above shows what actually happens on Java,
                //I tried my best to implement it the similiar way I could in C#
                sb.Append(String.Format(pat.getLocalizedString(locale),args));
                if (++count < len)
                {
                    // Add a separator if there are other messages.
                    sb.Append(separator);
                }
            }
            return sb.ToString();
        }
    }
}
