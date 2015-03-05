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
namespace Math3.analysis
{
    /// <summary>
    /// An interface representing a univariate real function.
    /// <para/>
    /// When a user-defined function encounters an error during
    /// evaluation, the <see cref="value(double)"/> method should throw a
    /// >user-defined unchecked exception.
    /// <para/>
    /// The following code excerpt shows the recommended way to do that using
    /// a root solver as an example, but the same construct is applicable to
    /// ODE integrators or optimizers.
    /// <code>
    /// private static class LocalException extends RuntimeException {
    ///     // The x value that caused the problem.
    ///     private final double x;
    ///
    ///     public LocalException(double x) {
    ///         this.x = x;
    ///     }
    ///
    ///     public double getX() {
    ///         return x;
    ///     }
    /// }
    ///
    /// private static class MyFunction implements UnivariateFunction {
    ///     public double value(double x) {
    ///         double y = hugeFormula(x);
    ///         if (somethingBadHappens) {
    ///           throw new LocalException(x);
    ///         }
    ///         return y;
    ///     }
    /// }
    ///
    /// public void compute() {
    ///     try {
    ///         solver.solve(maxEval, new MyFunction(a, b, c), min, max);
    ///     } catch (LocalException le) {
    ///         // Retrieve the x value.
    ///     }
    /// }
    /// </code>
    ///
    /// As shown, the exception is local to the user's code and it is guaranteed
    /// that Apache Commons Math will not catch it.
    /// </summary>
    public interface UnivariateFunction
    {
        /// <summary>
        /// Compute the value of the function.
        /// </summary>
        /// <param name="x">Point at which the function value should be computed.</param>
        /// <returns>the value of the function.</returns>
        /// <exception cref="IllegalArgumentException"> when the activated method itself can
        /// ascertain that a precondition, specified in the API expressed at the
        /// level of the activated method, has been violated.
        /// When Commons Math throws an <c>IllegalArgumentException</c>, it is
        /// usually the consequence of checking the actual parameters passed to
        /// the method.</exception>
        double value(double x);
    }
}
