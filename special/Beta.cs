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
using Math3.exception;
using Math3.util;
using System;

namespace Math3.special
{
    /// <summary>
    /// <para>
    /// This is a utility class that provides computation methods related to the
    /// Beta family of functions.
    /// </para>
    /// <para>
    /// Implementation of <see cref="logBeta(double, double)"/> is based on the
    /// algorithms described in
    /// <list type="bullet">
    /// <item><a href="http://dx.doi.org/10.1145/22721.23109">Didonato and Morris
    /// (1986)</a>, Computation of the Incomplete Gamma Function Ratios
    /// and their Inverse, TOMS 12(4), 377-393,</item>
    /// <item><a href="http://dx.doi.org/10.1145/131766.131776">Didonato and Morris
    /// (1992)</a>, Algorithm 708: Significant Digit Computation of the
    /// Incomplete Beta Function Ratios, TOMS 18(3), 360-373,</item>
    /// </list>
    /// and implemented in the
    /// <a href="http://www.dtic.mil/docs/citations/ADA476840">NSWC Library of Mathematical 
    /// Functions</a>, available
    /// <a href="http://www.ualberta.ca/CNS/RESEARCH/Software/NumericalNSWC/site.html">here</a>.
    /// This library is "approved for public release", and the
    /// <a href="http://www.dtic.mil/dtic/pdf/announcements/CopyrightGuidance.pdf">Copyright 
    /// guidance</a> indicates that unless otherwise stated in the code, all FORTRAN functions
    /// in this library are license free. Since no such notice appears in the code these
    /// functions can safely be ported to Commons-Math.
    /// </para>
    /// </summary>
    public class Beta
    {
        /// <summary>
        /// Maximum allowed numerical error.
        /// </summary>
        private const double DEFAULT_EPSILON = 1E-14;

        /// <summary>
        /// The constant value of ½log 2π.
        /// </summary>
        private const double HALF_LOG_TWO_PI = .9189385332046727;

        /// <summary>
        /// <para>
        /// The coefficients of the series expansion of the Δ function. This function
        /// is defined as follows
        /// </para>
        /// <c>Δ(x) = log Γ(x) - (x - 0.5) log a + a - 0.5 log 2π</c>,
        /// <para>
        /// see equation (23) in Didonato and Morris (1992). The series expansion,
        /// which applies for x ≥ 10, reads
        /// </para>
        /// <code>
        ///                 14
        ///               ====
        ///            1  \                2 n
        ///    Δ(x) = ---  >    d  (10 / x)
        ///            x  /      n
        ///               ====
        ///               n = 0
        /// <code>
        /// </summary>
        private static readonly double[] DELTA =
        {
            .833333333333333333333333333333E-01,
            -.277777777777777777777777752282E-04,
            .793650793650793650791732130419E-07,
            -.595238095238095232389839236182E-09,
            .841750841750832853294451671990E-11,
            -.191752691751854612334149171243E-12,
            .641025640510325475730918472625E-14,
            -.295506514125338232839867823991E-15,
            .179643716359402238723287696452E-16,
            -.139228964661627791231203060395E-17,
            .133802855014020915603275339093E-18,
            -.154246009867966094273710216533E-19,
            .197701992980957427278370133333E-20,
            -.234065664793997056856992426667E-21,
            .171348014966398575409015466667E-22
        };

        /// <summary>
        /// Default constructor.  Prohibit instantiation.
        /// </summary>
        private Beta() { }

        /// <summary>
        /// Returns the
        /// <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
        /// regularized beta function</a> I(x, a, b).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Parameter <c>a</c>.</param>
        /// <param name="b">Parameter <c>b</c>.</param>
        /// <returns>the regularized beta function I(x, a, b).</returns>
        /// <exception cref="MaxCountExceededException">
        /// if the algorithm fails to converge.</exception>
        public static double regularizedBeta(double x, double a, double b)
        {
            return regularizedBeta(x, a, b, DEFAULT_EPSILON, Int32.MaxValue);
        }

        /// <summary>
        /// Returns the
        /// <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
        /// regularized beta function</a> I(x, a, b).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Parameter <c>a</c>.</param>
        /// <param name="b">Parameter <c>b</c>.</param>
        /// <param name="epsilon"> When the absolute value of the nth item in the
        /// series is less than epsilon the approximation ceases to calculate
        /// further elements in the series.</param>
        /// <returns>the regularized beta function I(x, a, b)</returns>
        /// <exception cref="MaxCountExceededException">
        /// if the algorithm fails to converge.</exception>
        public static double regularizedBeta(double x, double a, double b, double epsilon)
        {
            return regularizedBeta(x, a, b, epsilon, Int32.MaxValue);
        }

        /// <summary>
        /// Returns the regularized beta function I(x, a, b).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Parameter <c>a</c>.</param>
        /// <param name="b">Parameter <c>b</c>.</param>
        /// <param name="maxIterations">Maximum number of "iterations" to complete.</param>
        /// <returns>the regularized beta function I(x, a, b)</returns>
        /// <exception cref="MaxCountExceededException">
        /// if the algorithm fails to converge.</exception>
        public static double regularizedBeta(double x, double a, double b, int maxIterations)
        {
            return regularizedBeta(x, a, b, DEFAULT_EPSILON, maxIterations);
        }

        /// <summary>
        /// Returns the regularized beta function I(x, a, b).<para/>
        /// The implementation of this method is based on:
        /// <list type="bullet">
        /// <item>
        /// <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
        /// Regularized Beta Function</a>.</item>
        /// <item>
        /// <a href="http://functions.wolfram.com/06.21.10.0001.01">
        /// Regularized Beta Function</a>.</item>
        /// </list>
        /// </summary>
        /// <param name="x">the value.</param>
        /// <param name="a">Parameter <c>a</c>.</param>
        /// <param name="b">Parameter <c>b</c>.</param>
        /// <param name="epsilon">When the absolute value of the nth item in the
        /// series is less than epsilon the approximation ceases to calculate
        /// further elements in the series.</param>
        /// <param name="maxIterations">Maximum number of "iterations" to complete.</param>
        /// <returns>the regularized beta function I(x, a, b)</returns>
        /// <exception cref="MaxCountExceededException">
        /// if the algorithm fails to converge.</exception>
        public static double regularizedBeta(double x, double a, double b, double epsilon, int maxIterations)
        {
            double ret;

            if (Double.IsNaN(x) ||
                Double.IsNaN(a) ||
                Double.IsNaN(b) ||
                x < 0 ||
                x > 1 ||
                a <= 0 ||
                b <= 0)
            {
                ret = Double.NaN;
            }
            else if (x > (a + 1) / (2 + b + a) &&
                     1 - x <= (b + 1) / (2 + b + a))
            {
                ret = 1 - regularizedBeta(1 - x, b, a, epsilon, maxIterations);
            }
            else
            {
                ContinuedFraction fraction = new ContinuedFractionHelper(a, b);
                ret = FastMath.exp((a * FastMath.log(x)) + (b * FastMath.log1p(-x)) -
                    FastMath.log(a) - logBeta(a, b)) *
                    1.0 / fraction.evaluate(x, epsilon, maxIterations);
            }

            return ret;
        }

        private class ContinuedFractionHelper : ContinuedFraction
        {
            private double a;
            private double b;
            internal ContinuedFractionHelper(double a, double b)
            {
                this.a = a;
                this.b = b;
            }
            protected override double getB(int n, double x)
            {
                double ret;
                double m;
                if (n % 2 == 0)
                { // even
                    m = n / 2.0;
                    ret = (m * (b - m) * x) /
                        ((a + (2 * m) - 1) * (a + (2 * m)));
                }
                else
                {
                    m = (n - 1.0) / 2.0;
                    ret = -((a + m) * (a + b + m) * x) /
                            ((a + (2 * m)) * (a + (2 * m) + 1.0));
                }
                return ret;
            }

            protected override double getA(int n, double x)
            {
                return 1.0;
            }

        };

        /// <summary>
        /// Returns the natural logarithm of the beta function B(a, b).<para/>
        /// The implementation of this method is based on:
        /// <list type="bullet">
        /// <item><a href="http://mathworld.wolfram.com/BetaFunction.html">
        /// Beta Function</a>, equation (1).</item>
        /// </list>.
        /// </summary>
        /// <param name="a">Parameter <c>a</c>.</param>
        /// <param name="b">Parameter <c>b</c>.</param>
        /// <param name="epsilon">This parameter is ignored.</param>
        /// <param name="maxIterations">This parameter is ignored.</returns>
        /// <returns>log(B(a, b)).</returns>
        [Obsolete("this method is deprecated as the computation of the beta function is no longer iterative")]
        public static double logBeta(double a, double b, double epsilon, int maxIterations)
        {
            return logBeta(a, b);
        }


        /// <summary>
        /// Returns the value of log Γ(a + b) for 1 ≤ a, b ≤ 2. Based on the
        /// NSWC Library of Mathematics Subroutines double precision
        /// implementation, <c>DGSMLN</c>. In <c>BetaTest.testLogGammaSum()</c>,
        /// this private method is accessed through reflection.
        /// </summary>
        /// <param name="a">First argument.</param>
        /// <param name="b">Second argument.</param>
        /// <returns>the value of <c>log(Gamma(a + b))</c>.</returns>
        /// <exception cref="OutOfRangeException">if <c>a</c> or <c>b</c> is lower than
        /// <c>1.0</c> or greater than <c>2.0</c>.</exception>
        private static double logGammaSum(double a, double b)
        {
            if ((a < 1.0) || (a > 2.0))
            {
                throw new OutOfRangeException<Double>(a, 1.0, 2.0);
            }
            if ((b < 1.0) || (b > 2.0))
            {
                throw new OutOfRangeException<Double>(b, 1.0, 2.0);
            }

            double x = (a - 1.0) + (b - 1.0);
            if (x <= 0.5)
            {
                return Gamma.logGamma1p(1.0 + x);
            }
            else if (x <= 1.5)
            {
                return Gamma.logGamma1p(x) + FastMath.log1p(x);
            }
            else
            {
                return Gamma.logGamma1p(x - 1.0) + FastMath.log(x * (1.0 + x));
            }
        }

        /// <summary>
        /// Returns the value of log[Γ(b) / Γ(a + b)] for a ≥ 0 and b ≥ 10. Based on
        /// the NSWC Library of Mathematics Subroutines double precision
        /// implementation, <c>DLGDIV</c>. In
        /// <c>BetaTest.testLogGammaMinusLogGammaSum()</c>, this private method is
        /// accessed through reflection.
        /// </summary>
        /// <param name="a">First argument.</param>
        /// <param name="b">Second argument.</param>
        /// <returns>the value of <c>log(Gamma(b) / Gamma(a + b))</c>.</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>a < 0.0</c> or <c>b < 10.0</c>.</exception>
        private static double logGammaMinusLogGammaSum(double a, double b)
        {
            if (a < 0.0)
            {
                throw new NumberIsTooSmallException<Double, Double>(a, 0.0, true);
            }
            if (b < 10.0)
            {
                throw new NumberIsTooSmallException<Double, Double>(b, 10.0, true);
            }

            /*
             * d = a + b - 0.5
             */
            double d;
            double w;
            if (a <= b)
            {
                d = b + (a - 0.5);
                w = deltaMinusDeltaSum(a, b);
            }
            else
            {
                d = a + (b - 0.5);
                w = deltaMinusDeltaSum(b, a);
            }

            double u = d * FastMath.log1p(a / b);
            double v = a * (FastMath.log(b) - 1.0);

            return u <= v ? (w - u) - v : (w - v) - u;
        }

        /// <summary>
        /// Returns the value of Δ(b) - Δ(a + b), with 0 ≤ a ≤ b and b ≥ 10. Based
        /// on equations (26), (27) and (28) in Didonato and Morris (1992).
        /// </summary>
        /// <param name="a">First argument.</param>
        /// <param name="b">Second argument.</param>
        /// <returns>the value of <c>Delta(b) - Delta(a + b)</c></returns>
        /// <exception cref="OutOfRangeException"> if <c>a < 0</c> or <c>a > b</c></exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>b < 10</c></exception>
        private static double deltaMinusDeltaSum(double a, double b)
        {
            if ((a < 0) || (a > b))
            {
                throw new OutOfRangeException<Double>(a, 0, b);
            }
            if (b < 10)
            {
                throw new NumberIsTooSmallException<Double, Double>(b, 10, true);
            }

            double h = a / b;
            double p = h / (1.0 + h);
            double q = 1.0 / (1.0 + h);
            double q2 = q * q;
            /*
             * s[i] = 1 + q + ... - q**(2 * i)
             */
            double[] s = new double[DELTA.Length];
            s[0] = 1.0;
            for (int i = 1; i < s.Length; i++)
            {
                s[i] = 1.0 + (q + q2 * s[i - 1]);
            }
            /*
             * w = Delta(b) - Delta(a + b)
             */
            double sqrtT = 10.0 / b;
            double t = sqrtT * sqrtT;
            double w = DELTA[DELTA.Length - 1] * s[s.Length - 1];
            for (int i = DELTA.Length - 2; i >= 0; i--)
            {
                w = t * w + DELTA[i] * s[i];
            }
            return w * p / b;
        }

        /// <summary>
        /// Returns the value of Δ(p) + Δ(q) - Δ(p + q), with p, q ≥ 10. Based on
        /// the NSWC Library of Mathematics Subroutines double precision
        /// implementation, <c>DBCORR</c>. In
        /// <c>BetaTest.testSumDeltaMinusDeltaSum()</c>, this private method is
        /// accessed through reflection.
        /// </summary>
        /// <param name="p">First argument.</param>
        /// <param name="q">Second argument.</param>
        /// <returns>the value of <c>Delta(p) + Delta(q) - Delta(p + q)</c>.</returns>
        /// <exception cref="NumberIsTooSmallException"> if <c>p < 10.0</c> or <c>q < 10.0</c>.</exception>
        private static double sumDeltaMinusDeltaSum(double p, double q)
        {
            if (p < 10.0)
            {
                throw new NumberIsTooSmallException<Double, Double>(p, 10.0, true);
            }
            if (q < 10.0)
            {
                throw new NumberIsTooSmallException<Double, Double>(q, 10.0, true);
            }

            double a = FastMath.min(p, q);
            double b = FastMath.max(p, q);
            double sqrtT = 10.0 / a;
            double t = sqrtT * sqrtT;
            double z = DELTA[DELTA.Length - 1];
            for (int i = DELTA.Length - 2; i >= 0; i--)
            {
                z = t * z + DELTA[i];
            }
            return z / a + deltaMinusDeltaSum(a, b);
        }

        /// <summary>
        /// Returns the value of log B(p, q) for 0 ≤ x ≤ 1 and p, q > 0. Based on the
        /// NSWC Library of Mathematics Subroutines implementation, <c>DBETLN</c>.
        /// </summary>
        /// <param name="p">First argument.</param>
        /// <param name="q">Second argument.</param>
        /// <returns>the value of <c>log(Beta(p, q))</c>, <c>NaN</c> if
        /// <c>p <= 0</c> or <c>q <= 0</c>.</returns>
        public static double logBeta(double p, double q)
        {
            if (Double.IsNaN(p) || Double.IsNaN(q) || (p <= 0.0) || (q <= 0.0))
            {
                return Double.NaN;
            }

            double a = FastMath.min(p, q);
            double b = FastMath.max(p, q);
            if (a >= 10.0)
            {
                double w = sumDeltaMinusDeltaSum(a, b);
                double h = a / b;
                double c = h / (1.0 + h);
                double u = -(a - 0.5) * FastMath.log(c);
                double v = b * FastMath.log1p(h);
                if (u <= v)
                {
                    return (((-0.5 * FastMath.log(b) + HALF_LOG_TWO_PI) + w) - u) - v;
                }
                else
                {
                    return (((-0.5 * FastMath.log(b) + HALF_LOG_TWO_PI) + w) - v) - u;
                }
            }
            else if (a > 2.0)
            {
                if (b > 1000.0)
                {
                    int n = (int)FastMath.floor(a - 1.0);
                    double prod = 1.0;
                    double ared = a;
                    for (int i = 0; i < n; i++)
                    {
                        ared -= 1.0;
                        prod *= ared / (1.0 + ared / b);
                    }
                    return (FastMath.log(prod) - n * FastMath.log(b)) +
                            (Gamma.logGamma(ared) +
                             logGammaMinusLogGammaSum(ared, b));
                }
                else
                {
                    double prod1 = 1.0;
                    double ared = a;
                    while (ared > 2.0)
                    {
                        ared -= 1.0;
                        double h = ared / b;
                        prod1 *= h / (1.0 + h);
                    }
                    if (b < 10.0)
                    {
                        double prod2 = 1.0;
                        double bred = b;
                        while (bred > 2.0)
                        {
                            bred -= 1.0;
                            prod2 *= bred / (ared + bred);
                        }
                        return FastMath.log(prod1) +
                               FastMath.log(prod2) +
                               (Gamma.logGamma(ared) +
                               (Gamma.logGamma(bred) -
                                logGammaSum(ared, bred)));
                    }
                    else
                    {
                        return FastMath.log(prod1) +
                               Gamma.logGamma(ared) +
                               logGammaMinusLogGammaSum(ared, b);
                    }
                }
            }
            else if (a >= 1.0)
            {
                if (b > 2.0)
                {
                    if (b < 10.0)
                    {
                        double prod = 1.0;
                        double bred = b;
                        while (bred > 2.0)
                        {
                            bred -= 1.0;
                            prod *= bred / (a + bred);
                        }
                        return FastMath.log(prod) +
                               (Gamma.logGamma(a) +
                                (Gamma.logGamma(bred) -
                                 logGammaSum(a, bred)));
                    }
                    else
                    {
                        return Gamma.logGamma(a) +
                               logGammaMinusLogGammaSum(a, b);
                    }
                }
                else
                {
                    return Gamma.logGamma(a) +
                           Gamma.logGamma(b) -
                           logGammaSum(a, b);
                }
            }
            else
            {
                if (b >= 10.0)
                {
                    return Gamma.logGamma(a) +
                           logGammaMinusLogGammaSum(a, b);
                }
                else
                {
                    // The following command is the original NSWC implementation.
                    // return Gamma.logGamma(a) +
                    // (Gamma.logGamma(b) - Gamma.logGamma(a + b));
                    // The following command turns out to be more accurate.
                    return FastMath.log(Gamma.gamma(a) * Gamma.gamma(b) /
                                        Gamma.gamma(a + b));
                }
            }
        }
    }
}
