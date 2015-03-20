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
    /// Interface handling decomposition algorithms that can solve A &times; X = B.
    /// <para>
    /// Decomposition algorithms decompose an A matrix has a product of several specific
    /// matrices from which they can solve A &times; X = B in least squares sense: they find X
    /// such that ||A &times; X - B|| is minimal.
    /// <para>
    /// Some solvers like <see cref="LUDecomposition"/> can only find the solution for
    /// square matrices and when the solution is an exact linear solution, i.e. when
    /// ||A &times; X - B|| is exactly 0. Other solvers can also find solutions
    /// with non-square matrix A and with non-null minimal norm. If an exact linear
    /// solution exists it is also the minimal norm solution.
    /// </summary>
    public interface DecompositionSolver
    {
        /// <summary>
        /// Solve the linear equation A &times; X = B for matrices A.
        /// <para>
        /// The A matrix is implicit, it is provided by the underlying
        /// decomposition algorithm.</para>
        /// </summary>
        /// <param name="b">right-hand side of the equation A &times; X = B</param>
        /// <returns>vector X that minimizes the two norm of A &times; X - B</returns>
        /// <exception cref="DimensionMismatchException">
        /// if the matrices dimensions do not match.</exception>
        /// <exception cref="SingularMatrixException"> if the decomposed matrix is singular.
        /// </exception>
        RealVector solve(RealVector b);

        /// <summary>
        /// Solve the linear equation A &times; X = B for matrices A.
        /// <para>
        /// The A matrix is implicit, it is provided by the underlying
        /// decomposition algorithm.</para>
        /// </summary>
        /// <param name="b">right-hand side of the equation A &times; X = B</param>
        /// <returns>a matrix X that minimizes the two norm of A &times; X - B</returns>
        /// <exception cref="DimensionMismatchException">
        /// if the matrices dimensions do not match.</exception>
        /// <exception cref="SingularMatrixException"> if the decomposed matrix is singular.
        /// </exception>
        RealMatrix solve(RealMatrix b);

        /// <summary>
        /// Check if the decomposed matrix is non-singular.
        /// </summary>
        /// <returns>true if the decomposed matrix is non-singular.</returns>
        Boolean isNonSingular();

        /// <summary>
        /// Get the <a href="http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse">pseudo-inverse</a>
        /// of the decomposed matrix.\
        /// <para>
        /// This is equal to the inverse  of the decomposed matrix, if such an inverse exists.
        /// </para>
        /// <para>
        /// If no such inverse exists, then the result has properties that resemble that of
        /// an inverse.
        /// </para>
        /// <para>
        /// In particular, in this case, if the decomposed matrix is A, then the system 
        /// of equations \( A x = b \) may have no solutions, or many. If it has no solutions, 
        /// then the pseudo-inverse \( A^+ \) gives the "closest" solution \( z = A^+ b \),
        /// meaning \( \left \| A z - b \right \|_2 \) is minimized. If there are many solutions,
        /// then \( z = A^+ b \) is the smallest solution, meaning \( \left \| z \right \|_2 \)
        /// is minimized.
        /// </para>
        /// <para>
        /// Note however that some decompositions cannot compute a pseudo-inverse for all matrices.
        /// For example, the <see cref="LUDecomposition"/> is not defined for non-square 
        /// matrices to begin with. The <see cref="QRDecomposition"/> can operate on non-square 
        /// matrices, but will throw <see cref="SingularMatrixException"/> if the decomposed
        /// matrix is singular. Refer to the doc of specific decomposition implementations 
        /// for more details.</para>
        /// </summary>
        /// <returns>pseudo-inverse matrix (which is the inverse, if it exists),
        /// if the decomposition can pseudo-invert the decomposed matrix</returns>
        /// <exception cref="SingularMatrixException"> if the decomposed matrix is singular
        /// and the decomposition can not compute a pseudo-inverse</exception>
        RealMatrix getInverse();
    }
}