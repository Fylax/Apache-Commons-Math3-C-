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

namespace Math3
{
   /// <summary>
   /// Interface representing a <a href="http://mathworld.wolfram.com/Field.html">field</a>.
   /// <para>
   /// Classes implementing this interface will often be singletons.
   /// </para>
   /// </summary>
   /// <remarks>See <see cref="FieldElement"/></remarks>
   /// <typeparam name="T">the type of the field elements</typeparam>
    public interface Field<T>
    {
        /// <summary>
        /// Get the additive identity of the field.
        /// <para>
        /// The additive identity is the element e_0 of the field such that
        /// for all elements a of the field, the equalities a + e<sub>0</sub> =
        /// e_0 + a = a hold.
        /// </para>
        /// </summary>
        /// <returns>additive identity of the field</returns>
        T getZero();

        /// <summary>
        /// Get the multiplicative identity of the field.
        /// <para>
        /// The multiplicative identity is the element e<sub>1</sub> of the field such that
        /// for all elements a of the field, the equalities a &times; e_1 =
        /// e_1 &times; a = a hold.
        /// </para>
        /// </summary>
        /// <returns>multiplicative identity of the field</returns>
        T getOne();

        /// <summary>
        /// Returns the runtime type of the FieldElement.
        /// </summary>
        /// <returns>The <c>Type</c> object that represents the runtime
        /// Type of this object.</returns>
        Type getRuntimeClass();
    }
}
