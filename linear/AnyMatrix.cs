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
    /// Interface defining very basic matrix operations.
    /// </summary>
    public interface AnyMatrix
    {
        /// <summary>
        /// Is this a square matrix?
        /// </summary>
        /// <returns>true if the matrix is square (rowDimension = columnDimension)</returns>
        Boolean isSquare();

        /// <summary>
        /// Returns the number of rows in the matrix.
        /// </summary>
        /// <returns>rowDimension</returns>
        int getRowDimension();

        /// <summary>
        /// Returns the number of columns in the matrix. 
        /// </summary>
        /// <returns>columnDimension</returns>
        int getColumnDimension();
    }
}
