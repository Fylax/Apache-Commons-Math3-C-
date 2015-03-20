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
    /// Base class for all unsupported features.
    /// It is used for all the exceptions that have the semantics of the standard
    /// <see cref="UnsupportedOperationException"/>, but must also provide a localized
    /// message.
    /// </summary>
    public class MathUnsupportedOperationException : NotSupportedException, ExceptionContextProvider
    {
        /// <summary>
        /// Context.
        /// </summary>
        private readonly ExceptionContext context;

        /// <summary>
        /// Default constructor.
        /// </summary>
        public MathUnsupportedOperationException() : this(new LocalizedFormats("UNSUPPORTED_OPERATION")) { }
        
        /// <summary></summary>
        /// <param name="pattern">Message pattern providing the specific context of
        /// the error.</param>
        /// <param name="args">Arguments.</param>
        public MathUnsupportedOperationException(Localizable pattern, params Object[] args)
        {
            context = new ExceptionContext(this);
            context.addMessage(pattern, args);
        }

        /// <inheritdoc/>
        public ExceptionContext getContext()
        {
            return context;
        }

        /// <inheritdoc/>
        public String getMessage()
        {
            return context.getMessage();
        }

        /// <inheritdoc/>
        public String getLocalizedMessage()
        {
            return context.getLocalizedMessage();
        }
    }
}