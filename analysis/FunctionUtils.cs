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
using Math3.analysis.differentiation;
using Math3.analysis.function;
using Math3.exception;
using Math3.exception.util;
using System;

namespace Math3.analysis
{
    /// <summary>
    /// Utilities for manipulating function objects.
    /// </summary>
    #pragma warning disable 0618
    public class FunctionUtils
    {
        /// <summary>
        /// Class only contains static methods.
        /// </summary>
        private FunctionUtils() { }

        /// <summary>
        /// Composes functions.
        /// <para/>
        /// The functions in the argument list are composed sequentially, in the
        /// given order.  For example, compose(f1,f2,f3) acts like f1(f2(f3(x))).
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>the composite function.</returns>
        public static UnivariateFunction compose(params UnivariateFunction[] f)
        {
            return new UnivariateFunctionAnonymous1(f);
        }
        
        /// <summary>
        /// Composes functions.
        /// <para/>
        /// The functions in the argument list are composed sequentially, in the
        /// given order.  For example, compose(f1,f2,f3) acts like f1(f2(f3(x))).
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>the composite function.</returns>
        public static UnivariateDifferentiableFunction compose(params UnivariateDifferentiableFunction[] f)
        {
            return new UnivariateDifferentiableFunctionAnonymous1(f);
        }

        /// <summary>
        /// Composes functions.
        /// <para/>
        /// The functions in the argument list are composed sequentially, in the
        /// given order.  For example, compose(f1,f2,f3) acts like f1(f2(f3(x))).
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>the composite function.</returns>
        [Obsolete("replaced by compose(UnivariateDifferentiableFunction)")]
        public static DifferentiableUnivariateFunction compose(params DifferentiableUnivariateFunction[] f)
        {
            return new DifferentiableUnivariateFunctionAnonymous1(f);
        }

        /// <summary>
        /// Adds functions.
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>a function that computes the sum of the functions.</returns>
        public static UnivariateFunction add(params UnivariateFunction[] f)
        {
            return new UnivariateFunctionAnonymous3(f);
        }

        /// <summary>
        /// Adds functions.
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>a function that computes the sum of the functions.</returns>
        public static UnivariateDifferentiableFunction add(params UnivariateDifferentiableFunction[] f)
        {
            return new UnivariateDifferentiableFunctionAnonymous2(f);
        }

        /// <summary>
        /// Adds functions.
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>a function that computes the sum of the functions.</returns>
        [Obsolete("replaced by add(params UnivariateDifferentiableFunction)")]
        public static DifferentiableUnivariateFunction add(params DifferentiableUnivariateFunction[] f)
        {
            return new DifferentiableUnivariateFunctionAnonymous2(f);
        }

        /// <summary>
        /// Multiplies functions.
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>a function that computes the product of the functions.</returns>
        public static UnivariateFunction multiply(params UnivariateFunction[] f)
        {
            return new UnivariateFunctionAnonymous2(f);
        }

        /// <summary>
        /// Multiplies functions.
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>a function that computes the product of the functions.</returns>
        public static UnivariateDifferentiableFunction multiply(params UnivariateDifferentiableFunction[] f)
        {
            return new UnivariateDifferentiableFunctionAnonymous3(f);
        }

        /// <summary>
        /// Multiplies functions. 
        /// </summary>
        /// <param name="f">List of functions.</param>
        /// <returns>a function that computes the product of the functions.</returns>
        [Obsolete("replaced by multiply(params UnivariateDifferentiableFunction)")]
        public static DifferentiableUnivariateFunction multiply(params DifferentiableUnivariateFunction[] f)
        {
            return new DifferentiableUnivariateFunctionAnonymous3(f);
        }

        /// <summary>
        /// Returns the univariate function <para/>
        /// <c>h(x) = combiner(f(x), g(x))</c>.
        /// </summary>
        /// <param name="combiner">Combiner function.</param>
        /// <param name="f">Function.</param>
        /// <param name="g">Function.</param>
        /// <returns>the composite function.</returns>
        public static UnivariateFunction combine(BivariateFunction combiner, UnivariateFunction f, UnivariateFunction g)
        {
            return new UnivariateFunctionAnonymous4(combiner, f, g);
        }

        /// <summary>
        /// Returns a MultivariateFunction h(x[]) defined by <code>
        /// h(x[]) = combiner(...combiner(combiner(initialValue,f(x[0])),f(x[1]))...),f(x[x.Length-1]))
        /// </code>
        /// </summary>
        /// <param name="combiner">Combiner function.</param>
        /// <param name="f">Function.</param>
        /// <param name="initialValue">Initial value.</param>
        /// <returns>a collector function.</returns>
        public static MultivariateFunction collector(BivariateFunction combiner, UnivariateFunction f, double initialValue)
        {
            return new MultivariateFunctionAnonymous(combiner, f, ref initialValue);
        }

        /// <summary>
        /// Returns a MultivariateFunction h(x[]) defined by <code>
        /// h(x[]) = combiner(...combiner(combiner(initialValue,x[0]),x[1])...),x[x.Length-1])
        /// </code>
        /// </summary>
        /// <param name="combiner">Combiner function.</param>
        /// <param name="initialValue">Initial value.</param>
        /// <returns>a collector function.</returns>
        public static MultivariateFunction collector(BivariateFunction combiner, double initialValue)
        {
            return collector(combiner, new Identity(), initialValue);
        }

        /// <summary>
        /// Creates a unary function by fixing the first argument of a binary function.
        /// </summary>
        /// <param name="f">Binary function.</param>
        /// <param name="dfixed">Value to which the first argument of <c>f</c> is set.</param>
        /// <returns>the unary function h(x) = f(dfixed, x)</returns>
        public static UnivariateFunction fix1stArgument(BivariateFunction f, double dfixed)
        {
            return new UnivariateFunctionAnonymous5(f, ref dfixed);
        }

        /// <summary>
        /// Creates a unary function by fixing the second argument of a binary function.
        /// </summary>
        /// <param name="f">Binary function.</param>
        /// <param name="dfixed">Value to which the second argument of <c>f</c> is set.</param>
        /// <returns>the unary function h(x) = f(x, dfixed)</returns>
        public static UnivariateFunction fix2ndArgument(BivariateFunction f, double dfixed)
        {
            return new UnivariateFunctionAnonymous6(f, ref dfixed);
        }

        /// <summary>
        /// Samples the specified univariate real function on the specified interval.
        /// <para/>
        /// The interval is divided equally into <c>n</c> sections and sample points
        /// are taken from <c>min</c> to <c>max - (max - min) / n</c>; therefore
        /// <c>f</c> is not sampled at the upper bound <c>max</c>.
        /// </summary>
        /// <param name="f">Function to be sampled</param>
        /// <param name="min">Lower bound of the interval (included).</param>
        /// <param name="max">Upper bound of the interval (excluded).</param>
        /// <param name="n">Number of sample points.</param>
        /// <returns>the array of samples.</returns>
        /// <exception cref="NumberIsTooLargeException"> if the lower bound <c>min</c> is
        /// greater than, or equal to the upper bound <c>max</c>.</exception>
        /// <exception cref="NotStrictlyPositiveException"> if the number of sample points
        /// <c>n</c> is negative.</exception>
        public static double[] sample(UnivariateFunction f, double min, double max, int n)
        {
            if (n <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("NOT_POSITIVE_NUMBER_OF_SAMPLES"), n);
            }
            if (min >= max)
            {
                throw new NumberIsTooLargeException<Double, Double>(min, max, false);
            }

            double[] s = new double[n];
            double h = (max - min) / n;
            for (int i = 0; i < n; i++)
            {
                s[i] = f.value(min + i * h);
            }
            return s;
        }

        /// <summary>
        /// Convert a <see cref="UnivariateDifferentiableFunction"/> into a 
        /// <see cref="DifferentiableUnivariateFunction"/>.
        /// </summary>
        /// <param name="f">function to convert</param>
        /// <returns>converted function</returns>
        [Obsolete("this conversion method is temporary in version 3.1, as the  DifferentiableUnivariateFunction interface itself is deprecated")]
        public static DifferentiableUnivariateFunction toDifferentiableUnivariateFunction(UnivariateDifferentiableFunction f)
        {
            return new DifferentiableUnivariateFunctionAnonymous4(f);
        }
        
        /// <summary>
        /// Convert a <see cref="DifferentiableUnivariateFunction"/> into a
        /// <see cref="UnivariateDifferentiableFunction"/>.
        /// <para>
        /// Note that the converted function is able to handle <see cref="DerivativeStructure"/>
        /// up to order one.
        /// If the function is called with higher order, a 
        /// <see cref="NumberIsTooLargeException"/> will be thrown.
        /// </para>
        /// </summary>
        /// <param name="f">function to convert</param>
        /// <returns>converted function</returns>
        [Obsolete("this conversion method is temporary in version 3.1, as the DifferentiableUnivariateFunction interface itself is deprecated")]
        public static UnivariateDifferentiableFunction toUnivariateDifferential(DifferentiableUnivariateFunction f)
        {
            return new UnivariateDifferentiableFunctionAnonymous4(f);
        }
        
        /// <summary>
        /// Convert a <see cref="MultivariateDifferentiableFunction"/> into a
        /// <see cref="DifferentiableMultivariateFunction"/>.
        /// </summary>
        /// <param name="f">function to convert</param>
        /// <returns>converted function</returns>
        [Obsolete("this conversion method is temporary in version 3.1, as the DifferentiableMultivariateFunction interface itself is deprecated")]
        public static DifferentiableMultivariateFunction toDifferentiableMultivariateFunction(MultivariateDifferentiableFunction f)
        {
            return new DifferentiableMultivariateFunctionAnonymous(f);
        }
        
        /// <summary>
        /// Convert a <see cref="DifferentiableMultivariateFunction"/> into a 
        /// <see cref="MultivariateDifferentiableFunction"/>.
        /// <para>
        /// Note that the converted function is able to handle <see cref="DerivativeStructure"/>
        /// elements  that all have the same number of free parameters and order, and with order 
        /// at most 1.
        /// If the function is called with inconsistent numbers of free parameters or higher order,
        /// a <see cref="DimensionMismatchException"/> or a 
        /// <see cref="NumberIsTooLargeException"/> will be thrown.
        /// </para>
        /// </summary>
        /// <param name="f">function to convert</param>
        /// <returns>converted function</returns>
        [Obsolete("this conversion method is temporary in version 3.1, as the DifferentiableMultivariateFunction interface itself is deprecated")]
        public static MultivariateDifferentiableFunction toMultivariateDifferentiableFunction(DifferentiableMultivariateFunction f)
        {
            return new MultivariateDifferentiableFunctionAnonymous(f);
        }
        
        /// <summary>
        /// Convert a <see cref="MultivariateDifferentiableVectorFunction"/> into a
        /// <see cref="DifferentiableMultivariateVectorFunction"/>.
        /// </summary>
        /// <param name="f">function to convert</param>
        /// <returns>converted function</returns>
        [Obsolete("this conversion method is temporary in version 3.1, as the   DifferentiableMultivariateVectorFunction interface itself is deprecated")]
        public static DifferentiableMultivariateVectorFunction toDifferentiableMultivariateVectorFunction(MultivariateDifferentiableVectorFunction f)
        {
            return new DifferentiableMultivariateVectorFunctionAnonymous(f);
        }
        
        /// <summary>
        /// Convert a <see cref="DifferentiableMultivariateVectorFunction"/> into a
        /// <see cref="MultivariateDifferentiableVectorFunction"/>.
        /// <para>
        /// Note that the converted function is able to handle <see cref="DerivativeStructure"/>
        /// elements that all have the same number of free parameters and order, and with order
        /// at most 1.
        /// If the function is called with inconsistent numbers of free parameters or higher order,
        /// a <see cref="DimensionMismatchException"/> or a 
        /// <see cref="NumberIsTooLargeException"/> will be thrown.
        /// </para>
        /// </summary>
        /// <param name="f">function to convert</param>
        /// <returns>converted function</returns>
        [Obsolete("this conversion method is temporary in version 3.1, as the DifferentiableMultivariateFunction interface itself is deprecated")]
        public static MultivariateDifferentiableVectorFunction toMultivariateDifferentiableVectorFunction(DifferentiableMultivariateVectorFunction f)
        {
            return new MultivariateDifferentiableVectorFunctionAnonymous(f);
        }

        private class UnivariateFunctionAnonymous1 : UnivariateFunction
        {
            private readonly UnivariateFunction[] f;
            public UnivariateFunctionAnonymous1(params UnivariateFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                double r = x;
                for (int i = f.Length - 1; i >= 0; i--)
                {
                    r = f[i].value(r);
                }
                return r;
            }
        }

        private class UnivariateFunctionAnonymous3 : UnivariateFunction
        {
            private readonly UnivariateFunction[] f;
            public UnivariateFunctionAnonymous3(params UnivariateFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                double r = f[0].value(x);
                for (int i = 1; i < f.Length; i++)
                {
                    r += f[i].value(x);
                }
                return r;
            }
        }

        private class UnivariateFunctionAnonymous2 : UnivariateFunction
        {
            private readonly UnivariateFunction[] f;
            public UnivariateFunctionAnonymous2(params UnivariateFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                double r = f[0].value(x);
                for (int i = 1; i < f.Length; i++)
                {
                    r *= f[i].value(x);
                }
                return r;
            }
        }

        private class UnivariateFunctionAnonymous4 : UnivariateFunction
        {
            private readonly BivariateFunction combiner;
            private readonly UnivariateFunction f;
            private readonly UnivariateFunction g;
            public UnivariateFunctionAnonymous4(BivariateFunction combiner, UnivariateFunction f, UnivariateFunction g)
            {
                this.combiner = combiner;
                this.f = f;
                this.g = g;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                return combiner.value(f.value(x), g.value(x));
            }
        }

        private class UnivariateFunctionAnonymous5 : UnivariateFunction
        {
            private readonly BivariateFunction f;
            private readonly double dfixed;

            public UnivariateFunctionAnonymous5(BivariateFunction f, ref double dfixed)
            {
                this.f = f;
                this.dfixed = dfixed;
            }

            /// <inheritdoc/>
            public double value(double x)
            {
                return f.value(dfixed, x);
            }
        }
        
        private class UnivariateFunctionAnonymous6 : UnivariateFunction
        {
            private readonly BivariateFunction f;
            private readonly double dfixed;

            public UnivariateFunctionAnonymous6(BivariateFunction f, ref double dfixed)
            {
                this.f = f;
                this.dfixed = dfixed;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                return f.value(x, dfixed);
            }
        }

        private class UnivariateDifferentiableFunctionAnonymous1 : UnivariateDifferentiableFunction
        {
            private readonly UnivariateDifferentiableFunction[] f;
            public UnivariateDifferentiableFunctionAnonymous1(params UnivariateDifferentiableFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double t)
            {
                double r = t;
                for (int i = f.Length - 1; i >= 0; i--)
                {
                    r = f[i].value(r);
                }
                return r;
            }

            /// <inheritdoc/>
            public DerivativeStructure value(DerivativeStructure t)
            {
                DerivativeStructure r = t;
                for (int i = f.Length - 1; i >= 0; i--)
                {
                    r = f[i].value(r);
                }
                return r;
            }
        }

        private class UnivariateDifferentiableFunctionAnonymous2 : UnivariateDifferentiableFunction
        {
            private readonly UnivariateDifferentiableFunction[] f;
            public UnivariateDifferentiableFunctionAnonymous2(params UnivariateDifferentiableFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double t)
            {
                double r = f[0].value(t);
                for (int i = 1; i < f.Length; i++)
                {
                    r += f[i].value(t);
                }
                return r;
            }

            /// <inheritdoc/>
            /// <exception cref="DimensionMismatchException"> if functions are not consistent with each other</exception>
            public DerivativeStructure value(DerivativeStructure t)
            {
                DerivativeStructure r = f[0].value(t);
                for (int i = 1; i < f.Length; i++)
                {
                    r = r.add(f[i].value(t));
                }
                return r;
            }
        }

        private class UnivariateDifferentiableFunctionAnonymous3 : UnivariateDifferentiableFunction
        {
            private readonly UnivariateDifferentiableFunction[] f;

            public UnivariateDifferentiableFunctionAnonymous3(params UnivariateDifferentiableFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double t)
            {
                double r = f[0].value(t);
                for (int i = 1; i < f.Length; i++)
                {
                    r *= f[i].value(t);
                }
                return r;
            }

            /// <inheritdoc/>
            public DerivativeStructure value(DerivativeStructure t)
            {
                DerivativeStructure r = f[0].value(t);
                for (int i = 1; i < f.Length; i++)
                {
                    r = r.multiply(f[i].value(t));
                }
                return r;
            }
        }

        private class UnivariateDifferentiableFunctionAnonymous4 : UnivariateDifferentiableFunction
        {
            private readonly DifferentiableUnivariateFunction f;
            public UnivariateDifferentiableFunctionAnonymous4(DifferentiableUnivariateFunction f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                return f.value(x);
            }

            /// <inheritdoc/>
            /// <exception cref="NumberIsTooLargeException"> if derivation order is greater
            /// than 1</exception>
            public DerivativeStructure value(DerivativeStructure t)
            {
                switch (t.getOrder())
                {
                    case 0:
                        return new DerivativeStructure(t.getFreeParameters(), 0, f.value(t.getValue()));
                    case 1:
                        {
                            int parameters = t.getFreeParameters();
                            double[] derivatives = new double[parameters + 1];
                            derivatives[0] = f.value(t.getValue());
                            double fPrime = f.derivative().value(t.getValue());
                            int[] orders = new int[parameters];
                            for (int i = 0; i < parameters; ++i)
                            {
                                orders[i] = 1;
                                derivatives[i + 1] = fPrime * t.getPartialDerivative(orders);
                                orders[i] = 0;
                            }
                            return new DerivativeStructure(parameters, 1, derivatives);
                        }
                    default:
                        throw new NumberIsTooLargeException<Int32, Int32>(t.getOrder(), 1, true);
                }
            }
        }

        private class DifferentiableUnivariateFunctionAnonymous1 : DifferentiableUnivariateFunction
        {
            private readonly DifferentiableUnivariateFunction[] f;
            public DifferentiableUnivariateFunctionAnonymous1(params DifferentiableUnivariateFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                double r = x;
                for (int i = f.Length - 1; i >= 0; i--)
                {
                    r = f[i].value(r);
                }
                return r;
            }

            /// <inheritdoc/>
            public UnivariateFunction derivative()
            {
                return new UnivariateFunctionAnonymous(f);
            }

            private class UnivariateFunctionAnonymous : UnivariateFunction
            {
                private readonly DifferentiableUnivariateFunction[] f;
                public UnivariateFunctionAnonymous(params DifferentiableUnivariateFunction[] f)
                {
                    this.f = f;
                }
                /// <inheritdoc/>
                public double value(double x)
                {
                    double p = 1;
                    double r = x;
                    for (int i = f.Length - 1; i >= 0; i--)
                    {
                        p *= f[i].derivative().value(r);
                        r = f[i].value(r);
                    }
                    return p;
                }
            }
        }

        private class DifferentiableUnivariateFunctionAnonymous2 : DifferentiableUnivariateFunction
        {
            private readonly DifferentiableUnivariateFunction[] f;

            public DifferentiableUnivariateFunctionAnonymous2(params DifferentiableUnivariateFunction[] f)
            {
                this.f = f;
            }

            /// <inheritdoc/>
            public double value(double x)
            {
                double r = f[0].value(x);
                for (int i = 1; i < f.Length; i++)
                {
                    r += f[i].value(x);
                }
                return r;
            }

            /// <inheritdoc/>
            public UnivariateFunction derivative()
            {
                return new UnivariateFunctionAnonymous(f);
            }

            private class UnivariateFunctionAnonymous : UnivariateFunction
            {
                private readonly DifferentiableUnivariateFunction[] f;
                public UnivariateFunctionAnonymous(params DifferentiableUnivariateFunction[] f)
                {
                    this.f = f;
                }
                /// <inheritdoc/>
                public double value(double x)
                {
                    double r = f[0].derivative().value(x);
                    for (int i = 1; i < f.Length; i++)
                    {
                        r += f[i].derivative().value(x);
                    }
                    return r;
                }
            }
        }

        private class DifferentiableUnivariateFunctionAnonymous3 : DifferentiableUnivariateFunction
        {
            private readonly DifferentiableUnivariateFunction[] f;
            public DifferentiableUnivariateFunctionAnonymous3(params DifferentiableUnivariateFunction[] f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                double r = f[0].value(x);
                for (int i = 1; i < f.Length; i++)
                {
                    r *= f[i].value(x);
                }
                return r;
            }

            /// <inheritdoc/>
            public UnivariateFunction derivative()
            {
                return new UnivariateFunctionAnonymous(f);
            }

            private class UnivariateFunctionAnonymous : UnivariateFunction
            {
                private readonly DifferentiableUnivariateFunction[] f;
                public UnivariateFunctionAnonymous(params DifferentiableUnivariateFunction[] f)
                {
                    this.f = f;
                }
                /// <inheritdoc/>
                public double value(double x)
                {
                    double sum = 0;
                    for (int i = 0; i < f.Length; i++)
                    {
                        double prod = f[i].derivative().value(x);
                        for (int j = 0; j < f.Length; j++)
                        {
                            if (i != j)
                            {
                                prod *= f[j].value(x);
                            }
                        }
                        sum += prod;
                    }
                    return sum;
                }
            }
        }


        private class DifferentiableUnivariateFunctionAnonymous4 : DifferentiableUnivariateFunction
        {
            private readonly UnivariateDifferentiableFunction f;
            public DifferentiableUnivariateFunctionAnonymous4(UnivariateDifferentiableFunction f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double value(double x)
            {
                return f.value(x);
            }

            /// <inheritdoc/>
            public UnivariateFunction derivative()
            {
                return new UnivariateFunctionAnonymous(f);
            }
            private class UnivariateFunctionAnonymous : UnivariateFunction
            {
                private readonly UnivariateDifferentiableFunction f;
                public UnivariateFunctionAnonymous(UnivariateDifferentiableFunction f)
                {
                    this.f = f;
                }
                public double value(double x)
                {
                    return f.value(new DerivativeStructure(1, 1, 0, x)).getPartialDerivative(1);
                }
            }
        }

        private class MultivariateFunctionAnonymous : MultivariateFunction
        {
            private readonly BivariateFunction combiner;
            private readonly UnivariateFunction f;
            private readonly double initialValue;

            public MultivariateFunctionAnonymous(BivariateFunction combiner, UnivariateFunction f, ref double initialValue)
            {
                this.combiner = combiner;
                this.f = f;
                this.initialValue = initialValue;
            }
            /// <inheritdoc/>
            public double value(double[] point)
            {
                double result = combiner.value(initialValue, f.value(point[0]));
                for (int i = 1; i < point.Length; i++)
                {
                    result = combiner.value(result, f.value(point[i]));
                }
                return result;
            }
        }

        private class MultivariateDifferentiableVectorFunctionAnonymous : MultivariateDifferentiableVectorFunction
        {
            private readonly DifferentiableMultivariateVectorFunction f;
            public MultivariateDifferentiableVectorFunctionAnonymous(DifferentiableMultivariateVectorFunction f)
            {
                this.f = f;
            }

            /// <inheritdoc/>
            public double[] value(double[] x)
            {
                return f.value(x);
            }

            /// <inheritdoc/>
            /// <exception cref="NumberIsTooLargeException"> if derivation order is higher 
            /// than 1</exception>
            /// <exception cref="DimensionMismatchException"> if numbers of free parameters 
            /// are inconsistent</exception>
            public DerivativeStructure[] value(DerivativeStructure[] t)
            {
                // check parameters and orders limits
                int parameters = t[0].getFreeParameters();
                int order = t[0].getOrder();
                int n = t.Length;
                if (order > 1)
                {
                    throw new NumberIsTooLargeException<Int32, Int32>(order, 1, true);
                }

                // check all elements in the array are consistent
                for (int i = 0; i < n; ++i)
                {
                    if (t[i].getFreeParameters() != parameters)
                    {
                        throw new DimensionMismatchException(t[i].getFreeParameters(), parameters);
                    }

                    if (t[i].getOrder() != order)
                    {
                        throw new DimensionMismatchException(t[i].getOrder(), order);
                    }
                }

                // delegate computation to underlying function
                double[] point = new double[n];
                for (int i = 0; i < n; ++i)
                {
                    point[i] = t[i].getValue();
                }
                double[] value = f.value(point);
                double[][] jacobian = f.jacobian().value(point);

                // merge value and Jacobian into a DerivativeStructure array
                DerivativeStructure[] merged = new DerivativeStructure[value.Length];
                for (int k = 0; k < merged.Length; ++k)
                {
                    double[] derivatives = new double[parameters + 1];
                    derivatives[0] = value[k];
                    int[] orders = new int[parameters];
                    for (int i = 0; i < parameters; ++i)
                    {
                        orders[i] = 1;
                        for (int j = 0; j < n; ++j)
                        {
                            derivatives[i + 1] += jacobian[k][j] * t[j].getPartialDerivative(orders);
                        }
                        orders[i] = 0;
                    }
                    merged[k] = new DerivativeStructure(parameters, order, derivatives);
                }

                return merged;
            }
        }

        private class DifferentiableMultivariateVectorFunctionAnonymous : DifferentiableMultivariateVectorFunction
        {
            private readonly MultivariateDifferentiableVectorFunction f;
            public DifferentiableMultivariateVectorFunctionAnonymous(MultivariateDifferentiableVectorFunction f)
            {
                this.f = f;
            }
            /// <inheritdoc/>
            public double[] value(double[] x)
            {
                return f.value(x);
            }

            public MultivariateMatrixFunction jacobian()
            {
                return new MultivariateMatrixFunctionAnonymous1(f);
            }

            private class MultivariateMatrixFunctionAnonymous1 : MultivariateMatrixFunction
            {
                private readonly MultivariateDifferentiableVectorFunction f;
                public MultivariateMatrixFunctionAnonymous1(MultivariateDifferentiableVectorFunction f)
                {
                    this.f = f;
                }
                /// <inheritdoc/>
                public double[][] value(double[] x)
                {

                    int n = x.Length;

                    // delegate computation to underlying function
                    DerivativeStructure[] dsX = new DerivativeStructure[n];
                    for (int i = 0; i < n; ++i)
                    {
                        dsX[i] = new DerivativeStructure(n, 1, i, x[i]);
                    }
                    DerivativeStructure[] y = f.value(dsX);

                    // extract Jacobian
                    double[][] jacobian = new double[y.Length][];
                    int[] orders = new int[n];
                    for (int i = 0; i < y.Length; ++i)
                    {
                        for (int j = 0; j < n; ++j)
                        {
                            orders[j] = 1;
                            jacobian[i][j] = y[i].getPartialDerivative(orders);
                            orders[j] = 0;
                        }
                    }

                    return jacobian;
                }
            }
        }

        private class MultivariateDifferentiableFunctionAnonymous : MultivariateDifferentiableFunction
        {
            private readonly DifferentiableMultivariateFunction f;
            public MultivariateDifferentiableFunctionAnonymous(DifferentiableMultivariateFunction f)
            {
                this.f = f;
            }

            /// <inheritdoc/>
            public double value(double[] x)
            {
                return f.value(x);
            }

            /// <inheritdoc/>
            /// <exception cref="NumberIsTooLargeException"> if derivation order is higher than 
            /// 1</exception>
            /// <exception cref="DimensionMismatchException"> if numbers of free parameters 
            /// are inconsistent</exception>
            public DerivativeStructure value(DerivativeStructure[] t)
            {

                // check parameters and orders limits
                int parameters = t[0].getFreeParameters();
                int order = t[0].getOrder();
                int n = t.Length;
                if (order > 1)
                {
                    throw new NumberIsTooLargeException<Int32, Int32>(order, 1, true);
                }

                // check all elements in the array are consistent
                for (int i = 0; i < n; ++i)
                {
                    if (t[i].getFreeParameters() != parameters)
                    {
                        throw new DimensionMismatchException(t[i].getFreeParameters(), parameters);
                    }

                    if (t[i].getOrder() != order)
                    {
                        throw new DimensionMismatchException(t[i].getOrder(), order);
                    }
                }

                // delegate computation to underlying function
                double[] point = new double[n];
                for (int i = 0; i < n; ++i)
                {
                    point[i] = t[i].getValue();
                }
                double value = f.value(point);
                double[] gradient = f.gradient().value(point);

                // merge value and gradient into one DerivativeStructure
                double[] derivatives = new double[parameters + 1];
                derivatives[0] = value;
                int[] orders = new int[parameters];
                for (int i = 0; i < parameters; ++i)
                {
                    orders[i] = 1;
                    for (int j = 0; j < n; ++j)
                    {
                        derivatives[i + 1] += gradient[j] * t[j].getPartialDerivative(orders);
                    }
                    orders[i] = 0;
                }

                return new DerivativeStructure(parameters, order, derivatives);

            }

        }

        private class DifferentiableMultivariateFunctionAnonymous : DifferentiableMultivariateFunction
        {
            private readonly MultivariateDifferentiableFunction f;
            public DifferentiableMultivariateFunctionAnonymous(MultivariateDifferentiableFunction f)
            {
                this.f = f;
            }

            /// <inheritdoc/>
            public double value(double[] x)
            {
                return f.value(x);
            }

            /// <inheritdoc/>
            public MultivariateFunction partialDerivative(int k)
            {
                return new MultivariateFunctionAnonymous2(f, ref k);
            }

            private class MultivariateFunctionAnonymous2 : MultivariateFunction
            {
                private readonly int k;
                private readonly MultivariateDifferentiableFunction f;
                public MultivariateFunctionAnonymous2(MultivariateDifferentiableFunction f, ref int k)
                {
                    this.k = k;
                    this.f = f;
                }
                /// <inheritdoc/>
                public double value(double[] x)
                {

                    int n = x.Length;

                    // delegate computation to underlying function
                    DerivativeStructure[] dsX = new DerivativeStructure[n];
                    for (int i = 0; i < n; ++i)
                    {
                        if (i == k)
                        {
                            dsX[i] = new DerivativeStructure(1, 1, 0, x[i]);
                        }
                        else
                        {
                            dsX[i] = new DerivativeStructure(1, 1, x[i]);
                        }
                    }
                    DerivativeStructure y = f.value(dsX);

                    // extract partial derivative
                    return y.getPartialDerivative(1);

                }
            }

            public MultivariateVectorFunction gradient()
            {
                return new MultivariateVectorFunctionAnonymous3(f);
            }

            private class MultivariateVectorFunctionAnonymous3 : MultivariateVectorFunction
            {
                private readonly MultivariateDifferentiableFunction f;
                public MultivariateVectorFunctionAnonymous3(MultivariateDifferentiableFunction f)
                {
                    this.f = f;
                }
                /// <inheritdoc/>
                public double[] value(double[] x)
                {

                    int n = x.Length;

                    // delegate computation to underlying function
                    DerivativeStructure[] dsX = new DerivativeStructure[n];
                    for (int i = 0; i < n; ++i)
                    {
                        dsX[i] = new DerivativeStructure(n, 1, i, x[i]);
                    }
                    DerivativeStructure y = f.value(dsX);

                    // extract gradient
                    double[] gradient = new double[n];
                    int[] orders = new int[n];
                    for (int i = 0; i < n; ++i)
                    {
                        orders[i] = 1;
                        gradient[i] = y.getPartialDerivative(orders);
                        orders[i] = 0;
                    }
                    return gradient;
                }
            }
        }
    }
}
