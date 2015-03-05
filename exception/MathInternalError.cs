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
    /// Exception triggered when something that shouldn't happen does happen.
    /// </summary>
    public class MathInternalError : MathIllegalStateException
    {
        /// <summary>
        /// URL for reporting problems.
        /// </summary>
        private const String REPORT_URL = "https://issues.apache.org/jira/browse/MATH";

        /// <summary>
        /// Simple constructor.
        /// </summary>
        public MathInternalError()
        {
            getContext().addMessage(new LocalizedFormats("INTERNAL_ERROR"), REPORT_URL);
        }

        /// <summary>
        /// Simple constructor.
        /// </summary>
        /// <param name="cause">root cause</param>
        public MathInternalError(Exception cause) : base(cause, new LocalizedFormats("INTERNAL_ERROR"), REPORT_URL) { }

        /// <summary>
        /// Constructor accepting a localized message. 
        /// </summary>
        /// <param name="pattern">Message pattern explaining the cause of the error.</param>
        /// <param name="args">Arguments.</param>
        public MathInternalError(Localizable pattern, params Object[] args) : base(pattern, args) { }
    }
}
