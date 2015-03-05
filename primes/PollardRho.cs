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

namespace Math3.primes
{
    /// <summary>
    /// Implementation of the Pollard's rho factorization algorithm.
    /// </summary>
    public class PollardRho
    {
        /// <summary>
        /// Hide utility class.
        /// </summary>
        private PollardRho() { }

        /// <summary>
        /// Factorization using Pollard's rho algorithm.
        /// </summary>
        /// <param name="n">number to factors, must be &gt; 0</param>
        /// <returns>the list of prime factors of n.</returns>
        public static List<Int32> primeFactors(int n)
        {
            List<Int32> factors = new List<Int32>();

            n = SmallPrimes.smallTrialDivision(n, factors);
            if (1 == n)
            {
                return factors;
            }

            if (SmallPrimes.millerRabinPrimeTest(n))
            {
                factors.Add(n);
                return factors;
            }

            int divisor = rhoBrent(n);
            factors.Add(divisor);
            factors.Add(n / divisor);
            return factors;
        }

        /// <summary>
        /// Implementation of the Pollard's rho factorization algorithm.
        /// <para>
        /// This implementation follows the paper "An improved Monte Carlo 
        /// factorization algorithm" by Richard P. Brent. This avoids the
        /// triple computation of f(x) typically found in Pollard's rho 
        /// implementations. It also batches several gcd computation into 1.
        /// </para><para>
        /// The backtracking is not implemented as we deal only with semi-primes.
        /// </para>
        /// </summary>
        /// <param name="n">number to factor, must be semi-prime.</param>
        /// <returns>a prime factor of n.</returns>
        internal static int rhoBrent(int n)
        {
            int x0 = 2;
            int m = 25;
            int cst = SmallPrimes.PRIMES_LAST;
            int y = x0;
            int r = 1;
            do
            {
                int x = y;
                for (int i = 0; i < r; i++)
                {
                    long y2 = ((long)y) * y;
                    y = (int)((y2 + cst) % n);
                }
                int k = 0;
                do
                {
                    int bound = FastMath.min(m, r - k);
                    int q = 1;
                    for (int i = -3; i < bound; i++)
                    { //start at -3 to ensure we enter this loop at least 3 times
                        long y2 = ((long)y) * y;
                        y = (int)((y2 + cst) % n);
                        long divisor = FastMath.abs(x - y);
                        if (0 == divisor)
                        {
                            cst += SmallPrimes.PRIMES_LAST;
                            k = -m;
                            y = x0;
                            r = 1;
                            break;
                        }
                        long prod = divisor * q;
                        q = (int)(prod % n);
                        if (0 == q)
                        {
                            return gcdPositive(FastMath.abs((int)divisor), n);
                        }
                    }
                    int outp = gcdPositive(FastMath.abs(q), n);
                    if (1 != outp)
                    {
                        return outp;
                    }
                    k += m;
                } while (k < r);
                r = 2 * r;
            } while (true);
        }

        /// <summary>
        /// Gcd between two positive numbers.
        /// <para>
        /// Gets the greatest common divisor of two numbers, using the "binary gcd" method,
        /// which avoids division and modulo operations. See Knuth 4.5.2 algorithm B.
        /// This algorithm is due to Josef Stein (1961).
        /// </para>
        /// Special cases:
        /// <list type="bullet">
        /// <item>The result of <c>gcd(x, x)</c>, <c>gcd(0, x)</c> and <c>gcd(x, 0)</c> is the value of <c>x</c>.</item>
        /// <item>The invocation <c>gcd(0, 0)</c> is the only one which returns <c>0</c>.</li>
        /// </list>
        /// </summary>
        /// <param name="a">first number, must be &ge; 0</param>
        /// <param name="b">second number, must be &ge; 0</param>
        /// <returns>gcd(a,b)</returns>
        internal static int gcdPositive(int a, int b)
        {
            // both a and b must be positive, it is not checked here
            // gdc(a,0) = a
            if (a == 0)
            {
                return b;
            }
            else if (b == 0)
            {
                return a;
            }

            // make a and b odd, keep in mind the common power of twos
            int mask = 1;
            int aTwos = 32;
            for (int i = 0; i < 32; ++i, mask <<= 1) //number of trailing zeros
            {
                if ((a & mask) != 0)
                {
                    aTwos = i;
                }
            }
            a >>= aTwos;
            mask = 1;
            int bTwos = 32;
            for (int i = 0; i < 32; ++i, mask <<= 1) //number of trailing zeros
            {
                if ((b & mask) != 0)
                {
                    bTwos = i;
                }
            }
            b >>= bTwos;
            int shift = FastMath.min(aTwos, bTwos);

            // a and b >0
            // if a > b then gdc(a,b) = gcd(a-b,b)
            // if a < b then gcd(a,b) = gcd(b-a,a)
            // so next a is the absolute difference and next b is the minimum of current values
            while (a != b)
            {
                int delta = a - b;
                b = FastMath.min(a, b);
                a = FastMath.abs(delta);
                // for speed optimization:
                // remove any power of two in a as b is guaranteed to be odd throughout all iterations
                mask = 1;
                int aTZ = 32;
                for (int i = 0; i < 32; ++i, mask <<= 1) //number of trailing zeros
                {
                    if ((a & mask) != 0)
                    {
                        aTZ = i;
                    }
                }
                a >>= aTZ;
            }

            // gcd(a,a) = a, just "add" the common power of twos
            return a << shift;
        }
    }
}
