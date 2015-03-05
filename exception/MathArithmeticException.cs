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
using Math3.exception.util;
using System;

namespace Math3.exception
{
    /// <summary>
    ///  Base class for arithmetic exceptions.
    ///  It is used for all the exceptions that have the semantics of the standard
    ///  <see cref="ArithmeticException"/>, but must also
    ///  provide a localized message.
    /// </summary>
    public class MathArithmeticException : ArithmeticException, ExceptionContextProvider
    {
        /// <summary>
        /// Context
        /// </summary>
        private ExceptionContext context;

        /// <summary>
        /// Default constructor.
        /// </summary>
        public MathArithmeticException() : this(new LocalizedFormats("ARITHMETIC_EXCEPTION"), null) { }

        /// <summary>
        /// Simple constructor.
        /// </summary>
        /// <param name="pattern">Message pattern explaining the cause of the error.</param>
        /// <param name="args">Arguments.</param>
        public MathArithmeticException(Localizable pattern, params Object[] args)
        {
            context = new ExceptionContext(this);
            context.addMessage(pattern, args);
        }

        /// <summary>
        /// Gets a reference to the "rich context" data structure that allows to
        /// customize error messages and store key, value pairs in exceptions.
        /// </summary>
        /// <returns>a reference to the exception context.</returns>
        public ExceptionContext getContext()
        {
            return context;
        }

        /// <summary>
        /// Gets the default message.
        /// </summary>
        /// <returns>the message.</returns>
        public String getMessage()
        {
            return context.getMessage();
        }

        /// <summary>
        /// Gets the message in the default locale.
        /// </summary>
        /// <returns>the localized message.</returns>
        public String getLocalizedMessage()
        {
            return context.getLocalizedMessage();
        }
    }
}
