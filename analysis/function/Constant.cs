﻿// Licensed to the Apache Software Foundation (ASF) under one or more
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
using Math3.analysis.differentiation;
using System;

namespace Math3.analysis.function
{
    /// <summary>
    /// Constant function.
    /// </summary>
#pragma warning disable 0618
    public class Constant : UnivariateDifferentiableFunction, DifferentiableUnivariateFunction
    {
        /// <summary>
        /// Constant.
        /// </summary>
        private readonly double c;

        /// <summary></summary>
        /// <param name="c">Constant.</param>
        public Constant(double c)
        {
            this.c = c;
        }

        /// <inheritdoc/>
        public double value(double x)
        {
            return c;
        }

        /// <inheritdoc/>
        [Obsolete("replaced by value(DerivativeStructure)")]
        public UnivariateFunction derivative()
        {
            return new Constant(0);
        }

        /// <inheritdoc/>
        public DerivativeStructure value(DerivativeStructure t)
        {
            return new DerivativeStructure(t.getFreeParameters(), t.getOrder(), c);
        }

        /// <inheritdoc/>
        double UnivariateFunction.value(double x)
        {
            throw new NotImplementedException();
        }
    }
}