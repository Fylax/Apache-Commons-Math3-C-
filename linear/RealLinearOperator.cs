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
using System;

namespace Math3.linear
{
    /// <summary>
    /// This class defines a linear operator operating on real (<c>double</c>)
    /// vector spaces. No direct access to the coefficients of the underlying matrix
    /// is provided
    /// The motivation for such an interface is well stated by
    /// <a href="#BARR1994">Barrett et al. (1994)</a>:
    /// <blockquote>
    ///  We restrict ourselves to iterative methods, which work by repeatedly
    ///  improving an approximate solution until it is accurate enough. These
    ///  methods access the coefficient matrix A of the linear system only via the
    ///  matrix-vector product y = A &middot; x
    ///  (and perhaps z = A<sup>T</sup> &middot; x). Thus the user need only
    ///  supply a subroutine for computing y (and perhaps z) given x, which permits
    ///  full exploitation of the sparsity or other special structure of A.
    /// </blockquote>
    /// <para/>
    /// <list type="bullet">
    /// <item><a name="BARR1994">Barret et al. (1994)</a></item>
    ///  <list type="bullet">
    ///   <item>
    ///    R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. M. Donato, J. Dongarra,
    ///    V. Eijkhout, R. Pozo, C. Romine and H. Van der Vorst,
    ///    Templates for the Solution of Linear Systems: Building Blocks for
    ///    Iterative Methods, SIAM
    ///   </item>
    ///  </list>
    /// </list>
    /// </summary>
    public abstract class RealLinearOperator
    {
        /// <summary>
        /// Returns the dimension of the codomain of this operator.
        /// </summary>
        /// <returns>the number of rows of the underlying matrix</returns>
        public abstract int getRowDimension();

        /// <summary>
        /// Returns the dimension of the domain of this operator.
        /// </summary>
        /// <returns>the number of columns of the underlying matrix</returns>
        public abstract int getColumnDimension();

        /// <summary>
        /// Returns the result of multiplying <c>this</c> by the vector <c>x</c>.
        /// </summary>
        /// <param name="x">the vector to operate on</param>
        /// <returns>the product of <c>this</c> instance with <c>x</c></returns>
        /// <exception cref="DimensionMismatchException"> if the column dimension does not
        /// match the size of <c>x</c></exception>
        public abstract RealVector operate(RealVector x);

        /// <summary>
        /// Returns the result of multiplying the transpose of <c>this</c> operator
        /// by the vector <c>x</c> (optional operation). The default implementation
        /// throws an <see cref="UnsupportedOperationException"/>. Users overriding this
        /// method must also override <see cref="isTransposable()"/>.
        /// </summary>
        /// <param name="x">the vector to operate on</param>
        /// <returns>the product of the transpose of <c>this</c> instance with
        /// <c>x</c></returns>
        /// <exception cref="DimensionMismatchException">
        /// if the row dimension does not match the size of <c>x</c></exception>
        /// <exception cref="UnsupportedOperationException"> if this operation is not 
        /// supported by <c>this</c> operator</exception>
        public RealVector operateTranspose(RealVector x)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        /// Returns <c>true</c> if this operator supports
        /// <see cref="operateTranspose(RealVector)"/>. If <c>true</c> is returned,
        /// <see cref="operateTranspose(RealVector)"/> should not throw
        /// <c>UnsupportedOperationException</c>. The default implementation returns
        /// <c>false</c>.
        /// </summary>
        /// <returns><c>false</c></returns>
        public Boolean isTransposable()
        {
            return false;
        }
    }
}