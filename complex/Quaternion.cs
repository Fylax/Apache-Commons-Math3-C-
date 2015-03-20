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
using Math3.exception.util;
using Math3.util;
using System;
using System.Text;

namespace Math3.complex
{
    /// <summary>
    /// This class implements <a href="http://mathworld.wolfram.com/Quaternion.html">
    /// quaternions</a> (Hamilton's hypercomplex numbers).
    /// <para/>
    /// Instance of this class are guaranteed to be immutable.
    /// </summary>
    /// <remarks>I arbitrarily decided to add the operators +, -, *
    /// which use the corrispondent method. This maybe was wanted
    /// by the original developers too, but it's impossible in Java,
    /// as far as I know.</remarks>
    public class Quaternion
    {
        /// <summary>
        /// Identity quaternion.
        /// </summary>
        public static readonly Quaternion IDENTITY = new Quaternion(1, 0, 0, 0);

        /// <summary>
        /// Zero quaternion.
        /// </summary>
        public static readonly Quaternion ZERO = new Quaternion(0, 0, 0, 0);
        
        /// <summary>
        /// i
        /// </summary>
        public static readonly Quaternion I = new Quaternion(0, 1, 0, 0);
        
        /// <summary>
        /// j
        /// </summary>
        public static readonly Quaternion J = new Quaternion(0, 0, 1, 0);
        
        /// <summary>
        /// k
        /// </summary>
        public static readonly Quaternion K = new Quaternion(0, 0, 0, 1);

        /// <summary>
        /// First component (scalar part).
        /// </summary>
        private readonly double q0;

        /// <summary>
        /// Second component (first vector part).
        /// </summary>
        private readonly double q1;
        
        /// <summary>
        /// Third component (second vector part).
        /// </summary>
        private readonly double q2;
       
        /// <summary>
        /// Fourth component (third vector part).
        /// </summary>
        private readonly double q3;

        /// <summary>
        /// Builds a quaternion from its components.
        /// </summary>
        /// <param name="a">Scalar component.</param>
        /// <param name="b">First vector component.</param>
        /// <param name="c">Second vector component.</param>
        /// <param name="d">Third vector component.</param>
        public Quaternion(double a, double b, double c, double d)
        {
            this.q0 = a;
            this.q1 = b;
            this.q2 = c;
            this.q3 = d;
        }

        /// <summary>
        /// Builds a quaternion from scalar and vector parts.
        /// </summary>
        /// <param name="scalar">Scalar part of the quaternion.</param>
        /// <param name="v">Components of the vector part of the quaternion.</param>
        /// <exception cref="DimensionMismatchException"> if the array length is not 3.</exception>
        public Quaternion(double scalar, double[] v)
        {
            if (v.Length != 3)
            {
                throw new DimensionMismatchException(v.Length, 3);
            }
            this.q0 = scalar;
            this.q1 = v[0];
            this.q2 = v[1];
            this.q3 = v[2];
        }

        /// <summary>
        /// Builds a pure quaternion from a vector (assuming that the scalar
        /// part is zero).
        /// </summary>
        /// <param name="v">Components of the vector part of the pure quaternion.</param>
        public Quaternion(double[] v) : this(0, v) { }

        /// <summary>
        /// Returns the conjugate quaternion of the instance.
        /// </summary>
        /// <returns>the conjugate quaternion</returns>
        public Quaternion getConjugate()
        {
            return new Quaternion(q0, -q1, -q2, -q3);
        }

        /// <summary>
        /// Returns the Hamilton product of two quaternions.
        /// </summary>
        /// <param name="q1">First quaternion.</param>
        /// <param name="q2">Second quaternion.</param>
        /// <returns>the product <c>q1</c> and <c>q2</c>, in that order.</returns>
        public static Quaternion multiply(Quaternion q1, Quaternion q2)
        {
            // Components of the first quaternion.
            double q1a = q1.getQ0();
            double q1b = q1.getQ1();
            double q1c = q1.getQ2();
            double q1d = q1.getQ3();

            // Components of the second quaternion.
            double q2a = q2.getQ0();
            double q2b = q2.getQ1();
            double q2c = q2.getQ2();
            double q2d = q2.getQ3();

            // Components of the product.
            double w = q1a * q2a - q1b * q2b - q1c * q2c - q1d * q2d;
            double x = q1a * q2b + q1b * q2a + q1c * q2d - q1d * q2c;
            double y = q1a * q2c - q1b * q2d + q1c * q2a + q1d * q2b;
            double z = q1a * q2d + q1b * q2c - q1c * q2b + q1d * q2a;

            return new Quaternion(w, x, y, z);
        }

        /// <summary>
        /// Returns the Hamilton product of the instance by a quaternion.n 
        /// </summary>
        /// <param name="q">Quaternion.</param>
        /// <returns>the product of this instance with <c>q</c>, in that order.</returns>
        public Quaternion multiply(Quaternion q)
        {
            return multiply(this, q);
        }

        public static Quaternion operator *(Quaternion Quaternion1, Quaternion Quaternion2)
        {
            return multiply(Quaternion1, Quaternion2);
        }

        /// <summary>
        /// Computes the sum of two quaternions.
        /// </summary>
        /// <param name="q1">Quaternion.</param>
        /// <param name="q2">Quaternion.</param>
        /// <returns>the sum of <c>q1</c> and <c>q2</c>.</returns>
        public static Quaternion add(Quaternion q1, Quaternion q2)
        {
            return new Quaternion(q1.getQ0() + q2.getQ0(),
                                  q1.getQ1() + q2.getQ1(),
                                  q1.getQ2() + q2.getQ2(),
                                  q1.getQ3() + q2.getQ3());
        }

        /// <summary>
        /// Computes the sum of the instance and another quaternion. 
        /// </summary>
        /// <param name="q">Quaternion.</param>
        /// <returns>the sum of this instance and <c>q</c></returns>
        public Quaternion add(Quaternion q)
        {
            return add(this, q);
        }

        public static Quaternion operator +(Quaternion Quaternion1, Quaternion Quaternion2)
        {
            return add(Quaternion1, Quaternion2);
        }

        /// <summary>
        /// Subtracts two quaternions.
        /// </summary>
        /// <param name="q1">First Quaternion.</param>
        /// <param name="q2">Second quaternion.</param>
        /// <returns>the difference between <c>q1</c> and <c>q2</c>.</returns>
        public static Quaternion subtract(Quaternion q1, Quaternion q2)
        {
            return new Quaternion(q1.getQ0() - q2.getQ0(),
                                  q1.getQ1() - q2.getQ1(),
                                  q1.getQ2() - q2.getQ2(),
                                  q1.getQ3() - q2.getQ3());
        }

        /// <summary>
        /// Subtracts a quaternion from the instance.
        /// </summary>
        /// <param name="q">Quaternion.</param>
        /// <returns>the difference between this instance and <c>q</c>.</returns>
        public Quaternion subtract(Quaternion q)
        {
            return subtract(this, q);
        }

        public static Quaternion operator -(Quaternion Quaternion1, Quaternion Quaternion2)
        {
            return subtract(Quaternion1, Quaternion2);
        }

        /// <summary>
        /// Computes the dot-product of two quaternions.
        /// </summary>
        /// <param name="q1">Quaternion.</param>
        /// <param name="q2">Quaternion.</param>
        /// <returns>the dot product of <c>q1</c> and <c>q2</c>.</returns>
        public static double dotProduct(Quaternion q1, Quaternion q2)
        {
            return q1.getQ0() * q2.getQ0() +
                q1.getQ1() * q2.getQ1() +
                q1.getQ2() * q2.getQ2() +
                q1.getQ3() * q2.getQ3();
        }

        /// <summary>
        /// Computes the dot-product of the instance by a quaternion. 
        /// </summary>
        /// <param name="q">Quaternion.</param>
        /// <returns>the dot product of this instance and <c>q</c></returns>
        public double dotProduct(Quaternion q)
        {
            return dotProduct(this, q);
        }

        /// <summary>
        /// Computes the norm of the quaternion.
        /// </summary>
        /// <returns>the norm.</returns>
        public double getNorm()
        {
            return FastMath.sqrt(q0 * q0 +
                                 q1 * q1 +
                                 q2 * q2 +
                                 q3 * q3);
        }

        /// <summary>
        /// Computes the normalized quaternion (the versor of the instance).
        /// The norm of the quaternion must not be zero.
        /// </summary>
        /// <returns>a normalized quaternion.</returns>
        /// <exception cref="ZeroException"> if the norm of the quaternion is zero.</exception>
        public Quaternion normalize()
        {
            double norm = getNorm();

            if (norm < Precision.SAFE_MIN)
            {
                throw new ZeroException(new LocalizedFormats("NORM"), norm);
            }

            return new Quaternion(q0 / norm,
                                  q1 / norm,
                                  q2 / norm,
                                  q3 / norm);
        }

        /// <inheritdoc/>
        public override Boolean Equals(Object other)
        {
            if (this == other)
            {
                return true;
            }
            if (other is Quaternion)
            {
                Quaternion q = (Quaternion)other;
                return q0 == q.getQ0() &&
                    q1 == q.getQ1() &&
                    q2 == q.getQ2() &&
                    q3 == q.getQ3();
            }
            return false;
        }

        /// <inheritdoc/>
        public override Int32 GetHashCode()
        {
            // "Effective Java" (second edition, p. 47).
            int result = 17;
            foreach (double comp in new double[] { q0, q1, q2, q3 })
            {
                int c = MathUtils.hash(comp);
                result = 31 * result + c;
            }
            return result;
        }

        /// <summary>
        /// Checks whether this instance is equal to another quaternion
        /// within a given tolerance.
        /// </summary>
        /// <param name="q">Quaternion with which to compare the current quaternion.</param>
        /// <param name="eps">Tolerance.</param>
        /// <returns><c>true</c> if the each of the components are equal
        /// within the allowed absolute error.</returns>
        public Boolean Equals(Quaternion q, double eps)
        {
            return Precision.equals(q0, q.getQ0(), eps) &&
                Precision.equals(q1, q.getQ1(), eps) &&
                Precision.equals(q2, q.getQ2(), eps) &&
                Precision.equals(q3, q.getQ3(), eps);
        }

        /// <summary>
        /// Checks whether the instance is a unit quaternion within a given
        /// tolerance. 
        /// </summary>
        /// <param name="eps">Tolerance (absolute error).</param>
        /// <returns><c>true</c> if the norm is 1 within the given tolerance,
        /// <c>false</c> otherwise</returns>
        public Boolean isUnitQuaternion(double eps)
        {
            return Precision.equals(getNorm(), 1d, eps);
        }

        /// <summary>
        /// Checks whether the instance is a pure quaternion within a given
        /// tolerance.
        /// </summary>
        /// <param name="eps">Tolerance (absolute error).</param>
        /// <returns><c>true</c> if the scalar part of the quaternion is zero.</returns>
        public Boolean isPureQuaternion(double eps)
        {
            return FastMath.abs(getQ0()) <= eps;
        }

        /// <summary>
        /// Returns the polar form of the quaternion.
        /// </summary>
        /// <returns>the unit quaternion with positive scalar part.</returns>
        public Quaternion getPositivePolarForm()
        {
            if (getQ0() < 0)
            {
                Quaternion unitQ = normalize();
                // The quaternion of rotation (normalized quaternion) q and -q
                // are equivalent (i.e. represent the same rotation).
                return new Quaternion(-unitQ.getQ0(),
                                      -unitQ.getQ1(),
                                      -unitQ.getQ2(),
                                      -unitQ.getQ3());
            }
            else
            {
                return this.normalize();
            }
        }

        /// <summary>
        /// Returns the inverse of this instance.
        /// The norm of the quaternion must not be zero.
        /// </summary>
        /// <returns>the inverse.</returns>
        /// <exception cref="ZeroException"> if the norm (squared) of the quaternion is zero.</exception>
        public Quaternion getInverse()
        {
            double squareNorm = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
            if (squareNorm < Precision.SAFE_MIN)
            {
                throw new ZeroException(new LocalizedFormats("NORM"), squareNorm);
            }

            return new Quaternion(q0 / squareNorm,
                                  -q1 / squareNorm,
                                  -q2 / squareNorm,
                                  -q3 / squareNorm);
        }

        /// <summary>
        /// Gets the first component of the quaternion (scalar part).
        /// </summary>
        /// <returns>the scalar part.</returns>
        public double getQ0()
        {
            return q0;
        }

        /// <summary>
        /// Gets the second component of the quaternion (first component
        /// of the vector part).
        /// </summary>
        /// <returns>the first component of the vector part.</returns>
        public double getQ1()
        {
            return q1;
        }

        /// <summary>
        /// Gets the third component of the quaternion (second component
        /// the vector part).
        /// </summary>
        /// <returns>the second component of the vector part.</returns>
        public double getQ2()
        {
            return q2;
        }

        /// <summary>
        /// Gets the fourth component of the quaternion (third component
        /// of the vector part).
        /// </summary>
        /// <returns>the third component of the vector part.</returns>
        public double getQ3()
        {
            return q3;
        }

        /// <summary>
        /// Gets the scalar part of the quaternion.
        /// </summary>
        /// <returns>the scalar part.</returns>
        /// <remarks>
        /// See <see cref="getQ0()"/>
        /// </remarks>
        public double getScalarPart()
        {
            return getQ0();
        }

        /// <summary>
        /// Gets the three components of the vector part of the quaternion.
        /// </summary>
        /// <returns>the vector part.</returns>
        /// <remarks>
        /// See <see cref="getQ1()"/><para/>
        /// See <see cref="getQ2()"/><para/>
        /// See<see cref="getQ3()"/><para/>
        /// </remarks>
        public double[] getVectorPart()
        {
            return new double[] { getQ1(), getQ2(), getQ3() };
        }

        /// <summary>
        /// Multiplies the instance by a scalar.
        /// </summary>
        /// <param name="alpha">Scalar factor.</param>
        /// <returns>a scaled quaternion.</returns>
        public Quaternion multiply(double alpha)
        {
            return new Quaternion(alpha * q0,
                                  alpha * q1,
                                  alpha * q2,
                                  alpha * q3);
        }

        /// <inheritdoc/>
        public override String ToString()
        {
            String sp = " ";
            StringBuilder s = new StringBuilder();
            s.Append("[").Append(q0).Append(sp)
                .Append(q1).Append(sp)
                .Append(q2).Append(sp)
                .Append(q3).Append("]");
            return s.ToString();
        }
    }

}
