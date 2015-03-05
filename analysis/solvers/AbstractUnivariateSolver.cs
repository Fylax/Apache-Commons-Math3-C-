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
using System.Runtime.InteropServices;

namespace Math3.analysis.solvers
{
    /// <summary>
    /// Base class for solvers.
    /// </summary>
    [ComVisible(false)]
    public abstract class AbstractUnivariateSolver : BaseAbstractUnivariateSolver<UnivariateFunction>, UnivariateSolver
    {
        /// <summary>
        /// Construct a solver with given absolute accuracy.
        /// </summary>
        /// <param name="absoluteAccuracy">Maximum absolute error.</param>
        protected AbstractUnivariateSolver(double absoluteAccuracy) : base(absoluteAccuracy) { }
        
        /// <summary>
        /// Construct a solver with given accuracies.
        /// </summary>
        /// <param name="relativeAccuracy">Maximum relative error.</param>
        /// <param name="absoluteAccuracy">Maximum absolute error.</param>
        protected AbstractUnivariateSolver(double relativeAccuracy, double absoluteAccuracy) : base(relativeAccuracy, absoluteAccuracy) { }
        
        /// <summary>
        /// Construct a solver with given accuracies.
        /// </summary>
        /// <param name="relativeAccuracy">Maximum relative error.</param>
        /// <param name="absoluteAccuracy">Maximum absolute error.</param>
        /// <param name="functionValueAccuracy">Maximum function value error.</param>
        protected AbstractUnivariateSolver(double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy) : base(relativeAccuracy, absoluteAccuracy, functionValueAccuracy) { }
    }
}
