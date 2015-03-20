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
using System;

namespace Math3.util
{
    /// <summary>
    /// <para>
    /// A variable length <see cref="DoubleArray"/> implementation that automatically
    /// handles expanding and contracting its internal storage array as elements
    /// are added and removed.
    /// </para>
    /// Important note: Usage should not assume that this class is thread-safe
    /// even though some of the methods are <c>synchronized</c>.
    /// This qualifier will be dropped in the next major release.
    /// <para>
    /// The internal storage array starts with capacity determined by the
    /// <c>initialCapacity</c> property, which can be set by the constructor.
    /// The default initial capacity is 16.  Adding elements using
    /// <see cref="addElement(double)"/> appends elements to the end of the array.
    /// When there are no open entries at the end of the internal storage array,
    /// the array is expanded.  The size of the expanded array depends on the
    /// <c>expansionMode</c> and <c>expansionFactor</c> properties.
    /// The <c>expansionMode</c> determines whether the size of the array is
    /// multiplied by the <c>expansionFactor</c>
    /// (<see cref="ExpansionMode.MULTIPLICATIVE"/>) or if the expansion is additive
    /// (<see cref="ExpansionMode.ADDITIVE"/> -- <c>expansionFactor</c> storage
    /// locations added).
    /// The default <c>expansionMode</c> is <c>MULTIPLICATIVE</c> and the default
    /// <c>expansionFactor</c> is 2.
    /// </para>
    /// <para>
    /// The <see cref="addElementRolling(double)"/> method adds a new element to the end
    /// of the internal storage array and adjusts the "usable window" of the
    /// internal array forward by one position (effectively making what was the
    /// second element the first, and so on).  Repeated activations of this method
    /// (or activation of <see cref="discardFrontElements(int)"/>) will effectively orphan
    /// the storage locations at the beginning of the internal storage array.  To
    /// reclaim this storage, each time one of these methods is activated, the size
    /// of the internal storage array is compared to the number of addressable
    /// elements (the <c>numElements</c> property) and if the difference
    /// is too large, the internal array is contracted to size
    /// <c>numElements + 1</c>.  The determination of when the internal
    /// storage array is "too large" depends on the <c>expansionMode</c> and
    /// <c>contractionFactor</c> properties.  If  the <c>expansionMode</c>
    /// is <c>MULTIPLICATIVE</c>, contraction is triggered when the
    /// ratio between storage array length and <c>numElements</c> exceeds
    /// <c>contractionFactor</c>.  If the <c>expansionMode</c>
    /// is <c>ADDITIVE</c>, the number of excess storage locations
    /// is compared to <c>contractionFactor</c>.
    /// </para>
    /// <para>
    /// To avoid cycles of expansions and contractions, the
    /// <c>expansionFactor</c> must not exceed the <c>contractionFactor</c>.
    /// Constructors and mutators for both of these properties enforce this
    /// requirement, throwing a <c>MathIllegalArgumentException</c> if it is
    /// violated.
    /// </para>
    /// </summary>
    public class ResizableDoubleArray : DoubleArray
    {
        /// <summary>
        /// Additive expansion mode.
        /// </summary>
        [Obsolete("Please use ExpansionMode.ADDITIVE instead")]
        public const int ADDITIVE_MODE = 1;

        /// <summary>
        /// Multiplicative expansion mode.
        /// </summary>
        [Obsolete("Please use ExpansionMode.MULTIPLICATIVE instead")]
        public const int MULTIPLICATIVE_MODE = 0;

        /// <summary>
        /// Default value for initial capacity.
        /// </summary>
        private const int DEFAULT_INITIAL_CAPACITY = 16;

        /// <summary>
        /// Default value for array size modifier.
        /// </summary>
        private const double DEFAULT_EXPANSION_FACTOR = 2.0;

        /// <summary>
        /// Default value for the difference between <see cref"contractionCriterion"/>
        /// and <see cref="expansionFactor"/>.
        /// </summary>
        private const double DEFAULT_CONTRACTION_DELTA = 0.5;

        /// <summary>
        /// The contraction criteria determines when the internal array will be
        /// contracted to fit the number of elements contained in the element
        /// array + 1.
        /// </summary>
        private double contractionCriterion = 2.5;

        /// <summary>
        /// The expansion factor of the array.  When the array needs to be expanded,
        /// <see cref="internalArray.length * expansionFactor"/>
        /// if <see cref="expansionMode"/> is set to MULTIPLICATIVE_MODE, or
        /// <see cref="internalArray.length + expansionFactor"/> if
        /// <see cref="expansionMode"/> is set to ADDITIVE_MODE.
        /// </summary>
        private double expansionFactor = 2.0;

        /// <summary>
        /// Determines whether array expansion by <c>expansionFactor</c>
        /// is additive or multiplicative.
        /// </summary>
        private ExpansionMode expansionMode = ExpansionMode.MULTIPLICATIVE;

        /// <summary>
        /// The internal storage array.
        /// </summary>
        private double[] internalArray;

        /// <summary>
        /// The number of addressable elements in the array.  Note that this
        /// has nothing to do with the length of the internal storage array.
        /// </summary>
        private int numElements = 0;

        /// <summary>
        /// The position of the first addressable element in the internal storage
        /// array.  The addressable elements in the array are
        /// <c>internalArray[startIndex],...,internalArray[startIndex + numElements - 1]</c>.
        /// </summary>
        private int startIndex = 0;

        /// <summary>
        /// Specification of expansion algorithm.
        /// </summary>
        public enum ExpansionMode
        {
            /// <summary>
            /// Multiplicative expansion mode.
            /// </summary>
            MULTIPLICATIVE,

            /// <summary>
            /// Additive expansion mode.
            /// </summary>
            ADDITIVE
        }

        /// <summary>
        /// Creates an instance with default properties.
        /// <list type="bullet">
        /// <item><c>initialCapacity = 16</c></item>
        /// <item><c>expansionMode = MULTIPLICATIVE</c></item>
        /// <item><c>expansionFactor = 2.0</c></item>
        /// <item><c>contractionCriterion = 2.5</c></item>
        /// </list>
        /// </summary>
        public ResizableDoubleArray() : this(DEFAULT_INITIAL_CAPACITY) { }

        /// <summary>
        /// Creates an instance with the specified initial capacity.
        /// Other properties take default values:
        /// <list type="bullet">
        /// <item><c>expansionMode = MULTIPLICATIVE</c></item>
        /// <item><c>expansionFactor = 2.0</c></item>
        /// <item><c>contractionCriterion = 2.5</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialCapacity">Initial size of the internal storage array.</param>
        /// <exception cref="MathIllegalArgumentException"> if <c>initialCapacity <= 0</c>
        /// .</exception>
        public ResizableDoubleArray(int initialCapacity) : this(initialCapacity, DEFAULT_EXPANSION_FACTOR) { }

        /// <summary>
        /// Creates an instance from an existing <c>double[]</c> with the
        /// initial capacity and numElements corresponding to the size of
        /// the supplied <c>double[]</c> array.
        /// If the supplied array is null, a new empty array with the default
        /// initial capacity will be created.
        /// The input array is copied, not referenced.
        /// Other properties take default values:
        /// <list type="bullet">
        /// <item><c>initialCapacity = 16</c></item>
        /// <item><c>expansionMode = MULTIPLICATIVE</c></item>
        /// <item><c>expansionFactor = 2.0</c></item>
        /// <item><c>contractionCriterion = 2.5</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialArray">initial array</param>
        public ResizableDoubleArray(double[] initialArray) : this(DEFAULT_INITIAL_CAPACITY, DEFAULT_EXPANSION_FACTOR, DEFAULT_CONTRACTION_DELTA + DEFAULT_EXPANSION_FACTOR, ExpansionMode.MULTIPLICATIVE, initialArray) { }

        /// <summary>
        /// Creates an instance with the specified initial capacity
        /// and expansion factor.
        /// The remaining properties take default values:
        /// <list type="bullet">
        /// <item><c>expansionMode = MULTIPLICATIVE</c></item>
        /// <item><c>contractionCriterion = 0.5 + expansionFactor</c></item>
        /// </list>
        /// <para/>
        /// Throws IllegalArgumentException if the following conditions are
        /// not met:
        /// <list type="bullet">
        /// <item><c>initialCapacity > 0</c></item>
        /// <item><c>expansionFactor > 1</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialCapacity">Initial size of the internal storage array.</param>
        /// <param name="expansionFactor">The array will be expanded based on this
        /// parameter.</param>
        /// <exception cref="MathIllegalArgumentException"> if parameters are not valid.
        /// </exception>
        [Obsolete("Please use ResizableDoubleArray(int,double) instead")]
        public ResizableDoubleArray(int initialCapacity, float expansionFactor) : this(initialCapacity, (double)expansionFactor) { }

        /// <summary>
        /// Creates an instance with the specified initial capacity
        /// and expansion factor.
        /// The remaining properties take default values:
        /// <list type="bullet">
        /// <item><c>expansionMode = MULTIPLICATIVE</c></item>
        /// <item><c>contractionCriterion = 0.5 + expansionFactor</c></item>
        /// </list>
        /// <para/>
        /// Throws IllegalArgumentException if the following conditions are
        /// not met:
        /// <list type="bullet">
        /// <item><c>initialCapacity > 0</c></item>
        /// <item><c>expansionFactor > 1</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialCapacity">Initial size of the internal storage array.</param>
        /// <param name="expansionFactor">The array will be expanded based on this
        /// parameter.</param>
        /// <exception cref="MathIllegalArgumentException"> if parameters are not valid.
        /// </exception>
        public ResizableDoubleArray(int initialCapacity, double expansionFactor) : this(initialCapacity, expansionFactor, DEFAULT_CONTRACTION_DELTA + expansionFactor) { }

        /// <summary>
        /// Creates an instance with the specified initialCapacity,
        /// expansionFactor, and contractionCriterion.
        /// The expansion mode will default to <c>MULTIPLICATIVE</c>.
        /// <para/>
        /// Throws IllegalArgumentException if the following conditions are
        /// not met:
        /// <list type="bullet">
        /// <item><c>initialCapacity > 0</c></item>
        /// <item><c>expansionFactor > 1</c></item>
        /// <item><c>contractionCriterion >= expansionFactor</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialCapacity">Initial size of the internal storage array.</param>
        /// <param name="expansionFactor">The array will be expanded based on this
        /// parameter.</param>
        /// <param name="contractionCriteria">Contraction criteria.</param>
        /// <exception cref="MathIllegalArgumentException"> if parameters are not valid.
        /// </exception>
        [Obsolete("Please use ResizableDoubleArray(int,double,double) instead")]
        public ResizableDoubleArray(int initialCapacity, float expansionFactor, float contractionCriteria) : this(initialCapacity, (double)expansionFactor, (double)contractionCriteria) { }

        /// <summary>
        /// Creates an instance with the specified initial capacity,
        /// expansion factor, and contraction criteria.
        /// The expansion mode will default to <c>MULTIPLICATIVE</c>.
        /// <para/>
        /// Throws IllegalArgumentException if the following conditions are
        /// not met:
        /// <list type="bullet">
        /// <item><c>initialCapacity > 0</c></item>
        /// <item><c>expansionFactor > 1</c></item>
        /// <item><c>contractionCriterion >= expansionFactor</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialCapacity">Initial size of the internal storage array.</param>
        /// <param name="expansionFactor">The array will be expanded based on this
        /// parameter.</param>
        /// <param name="contractionCriterion">Contraction criterion.</param>
        /// <exception cref="MathIllegalArgumentException"> if the parameters are not valid.
        /// </exception>
        public ResizableDoubleArray(int initialCapacity, double expansionFactor, double contractionCriterion) : this(initialCapacity, expansionFactor, contractionCriterion, ExpansionMode.MULTIPLICATIVE, null) { }

        /// <summary>
        /// <para>
        /// Create a ResizableArray with the specified properties.</para>
        /// <para>
        /// Throws IllegalArgumentException if the following conditions are
        /// not met:
        /// <list type="bullet">
        /// <item><c>initialCapacity > 0</c></item>
        /// <item><c>expansionFactor > 1</c></item>
        /// <item><c>contractionFactor >= expansionFactor</c></item>
        /// <item><c>expansionMode in {MULTIPLICATIVE_MODE, ADDITIVE_MODE}</c>
        /// </item>
        /// </list></para>
        /// </summary>
        /// <param name="initialCapacity">the initial size of the internal storage array</param>
        /// <param name="expansionFactor"> the array will be expanded based on this
        /// parameter</param>
        /// <param name="contractionCriteria">the contraction Criteria</param>
        /// <param name="expansionMode">the expansion mode</param>
        /// <exception cref="MathIllegalArgumentException"> if parameters are not valid.
        /// </exception>
        [Obsolete("Please use ResizableDoubleArray(int,double,double,ExpansionMode,double[]) instead")]
        public ResizableDoubleArray(int initialCapacity, float expansionFactor, float contractionCriteria, int expansionMode)
            : this(initialCapacity, expansionFactor, contractionCriteria, expansionMode == ADDITIVE_MODE ? ExpansionMode.ADDITIVE : ExpansionMode.MULTIPLICATIVE, null)
        {
            // XXX Just ot retain the expected failure in a unit test.
            // With the new "enum", that test will become obsolete.
            setExpansionMode(expansionMode);
        }

        /// <summary>
        /// Creates an instance with the specified properties.
        /// <para/>
        /// Throws MathIllegalArgumentException if the following conditions are
        /// not met:
        /// <list type="bullet">
        /// <item><c>initialCapacity > 0</c></item>
        /// <item><c>expansionFactor > 1</c></item>
        /// <item><c>contractionCriterion >= expansionFactor</c></item>
        /// </list>
        /// </summary>
        /// <param name="initialCapacity">Initial size of the internal storage array.</param>
        /// <param name="expansionFactor">The array will be expanded based on this
        /// parameter.</param>
        /// <param name="contractionCriterion">Contraction criteria.</param>
        /// <param name="expansionMode">Expansion mode.</param>
        /// <param name="data">Initial contents of the array.</param>
        /// <exception cref="MathIllegalArgumentException"> if the parameters are not valid.
        /// </exception>
        public ResizableDoubleArray(int initialCapacity, double expansionFactor, double contractionCriterion, ExpansionMode expansionMode, params double[] data)
        {
            if (initialCapacity <= 0)
            {
                throw new NotStrictlyPositiveException<Int32>(new LocalizedFormats("INITIAL_CAPACITY_NOT_POSITIVE"), initialCapacity);
            }
            checkContractExpand(contractionCriterion, expansionFactor);

            this.expansionFactor = expansionFactor;
            this.contractionCriterion = contractionCriterion;
            this.expansionMode = expansionMode;
            internalArray = new double[initialCapacity];
            numElements = 0;
            startIndex = 0;

            if (data != null && data.Length > 1)
            {
                addElements(data);
            }
        }

        /// <summary>
        /// Copy constructor.  Creates a new ResizableDoubleArray that is a deep,
        /// fresh copy of the original. Needs to acquire synchronization lock
        /// on original.  Original may not be null; otherwise a <see cref="NullArgumentException"/>
        /// is thrown.
        /// </summary>
        /// <param name="original">array to copy</param>
        /// <exception cref="NullArgumentException"> if original is null</exception>
        public ResizableDoubleArray(ResizableDoubleArray original)
        {
            MathUtils.checkNotNull(original);
            copy(original, this);
        }

        /// <summary>
        /// Adds an element to the end of this expandable array. 
        /// </summary>
        /// <param name="value">Value to be added to end of array.</param>
        public void addElement(double value)
        {
            lock (this)
            {
                if (internalArray.Length <= startIndex + numElements)
                {
                    expand();
                }
                internalArray[startIndex + numElements++] = value;
            }
        }

        /// <summary>
        /// Adds several element to the end of this expandable array.
        /// </summary>
        /// <param name="values">Values to be added to end of array.</param>
        public void addElements(double[] values)
        {
            lock (this)
            {
                double[] tempArray = new double[numElements + values.Length + 1];
                Array.Copy(internalArray, startIndex, tempArray, 0, numElements);
                Array.Copy(values, 0, tempArray, numElements, values.Length);
                internalArray = tempArray;
                startIndex = 0;
                numElements += values.Length;
            }
        }

        /// <summary>
        /// <para>
        /// Adds an element to the end of the array and removes the first
        /// element in the array.  Returns the discarded first element.
        /// The effect is similar to a push operation in a FIFO queue.
        /// </para>
        /// <para>
        /// Example: If the array contains the elements 1, 2, 3, 4 (in that order)
        /// and addElementRolling(5) is invoked, the result is an array containing
        /// the entries 2, 3, 4, 5 and the value returned is 1.
        /// </para>
        /// </summary>
        /// <param name="value">Value to be added to the array.</param>
        /// <returns>the value which has been discarded or "pushed" out of the array
        /// by this rolling insert.</returns>
        public double addElementRolling(double value)
        {
            lock (this)
            {
                double discarded = internalArray[startIndex];

                if ((startIndex + (numElements + 1)) > internalArray.Length)
                {
                    expand();
                }
                // Increment the start index
                startIndex += 1;

                // Add the new value
                internalArray[startIndex + (numElements - 1)] = value;

                // Check the contraction criterion.
                if (shouldContract())
                {
                    contract();
                }
                return discarded;
            }
        }

        /// <summary>
        /// Substitutes <c>value</c> for the most recently added value.
        /// Returns the value that has been replaced. If the array is empty (i.e.
        /// if <see cref="numElements"/> is zero), an IllegalStateException is thrown.
        /// </summary>
        /// <param name="value">New value to substitute for the most recently added value</param>
        /// <returns>the value that has been replaced in the array.</returns>
        /// <exception cref="MathIllegalStateException"> if the array is empty</exception>
        public double substituteMostRecentElement(double value)
        {
            if (numElements < 1)
            {
                throw new MathIllegalStateException(new LocalizedFormats("CANNOT_SUBSTITUTE_ELEMENT_FROM_EMPTY_ARRAY"));
            }
            lock (this)
            {
                int substIndex = startIndex + (numElements - 1);
                double discarded = internalArray[substIndex];

                internalArray[substIndex] = value;

                return discarded;
            }
        }

        /// <summary>
        /// Checks the expansion factor and the contraction criterion and throws an
        /// IllegalArgumentException if the contractionCriteria is less than the
        /// expansionCriteria
        /// </summary>
        /// <param name="contraction">expansion factor to be checked</param>
        /// <param name="expansion">contraction criteria to be checked</param>
        /// <exception cref="MathIllegalArgumentException"> if the contractionCriteria 
        /// is less than the expansionCriteria.</exception>
        [Obsolete("Please use checkContractExpand(double,double) instead")]
        protected void checkContractExpand(float contraction, float expansion)
        {
            checkContractExpand((double)contraction, (double)expansion);
        }

        /// <summary>
        /// Checks the expansion factor and the contraction criterion and raises
        /// expansion criterion.
        /// </summary>
        /// <param name="contraction">Criterion to be checked.</param>
        /// <param name="expansion">Factor to be checked.</param>
        /// <exception cref="NumberIsTooSmallException"> if <c>contraction < expansion</c>.
        /// </exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>contraction <= 1</c>.</exception>
        /// <exception cref="NumberIsTooSmallException"> if <c>expansion <= 1 </c>.</exception>
        protected void checkContractExpand(double contraction, double expansion)
        {
            if (contraction < expansion)
            {
                NumberIsTooSmallException<Double, Int32> e = new NumberIsTooSmallException<Double, Int32>(contraction, 1, true);
                e.getContext().addMessage(new LocalizedFormats("CONTRACTION_CRITERIA_SMALLER_THAN_EXPANSION_FACTOR"), contraction, expansion);
                throw e;
            }

            if (contraction <= 1)
            {
                NumberIsTooSmallException<Double, Int32> e = new NumberIsTooSmallException<Double, Int32>(contraction, 1, false);
                e.getContext().addMessage(new LocalizedFormats("CONTRACTION_CRITERIA_SMALLER_THAN_ONE"), contraction);
                throw e;
            }

            if (expansion <= 1)
            {
                NumberIsTooSmallException<Double, Int32> e = new NumberIsTooSmallException<Double, Int32>(contraction, 1, false);
                e.getContext().addMessage(new LocalizedFormats("EXPANSION_FACTOR_SMALLER_THAN_ONE"), expansion);
                throw e;
            }
        }

        /// <summary>
        /// Clear the array contents, resetting the number of elements to zero.
        /// </summary>
        public void clear()
        {
            numElements = 0;
            startIndex = 0;
        }

        /// <summary>
        /// Contracts the storage array to the (size of the element set) + 1 - to
        /// avoid a zero length array. This function also resets the startIndex to
        /// zero.
        /// </summary>
        public void contract()
        {
            lock (this)
            {
                double[] tempArray = new double[numElements + 1];

                // Copy and swap - copy only the element array from the src array.
                Array.Copy(internalArray, startIndex, tempArray, 0, numElements);
                internalArray = tempArray;

                // Reset the start index to zero
                startIndex = 0;
            }
        }

        /// <summary>
        /// Discards the <c>i</c> initial elements of the array.  For example,
        /// if the array contains the elements 1,2,3,4, invoking
        /// <c>discardFrontElements(2)</c> will cause the first two elements
        /// to be discarded, leaving 3,4 in the array.  Throws illegalArgumentException
        /// if i exceeds numElements.
        /// </summary>
        /// <param name="i">the number of elements to discard from the front of the array</param>
        /// <exception cref="MathIllegalArgumentException"> if i is greater than numElements.
        /// </exception>
        public void discardFrontElements(int i)
        {
            discardExtremeElements(i, true);
        }

        /// <summary>
        /// Discards the <c>i</c> last elements of the array.  For example,
        /// if the array contains the elements 1,2,3,4, invoking
        /// <c>discardMostRecentElements(2)</c> will cause the last two elements
        /// to be discarded, leaving 1,2 in the array.  Throws illegalArgumentException
        /// if i exceeds numElements.
        /// </summary>
        /// <param name="i">the number of elements to discard from the end of the array</param>
        /// <exception cref="MathIllegalArgumentException"> if i is greater than numElements.
        /// </exception>
        public void discardMostRecentElements(int i)
        {
            discardExtremeElements(i, false);
        }

        /// <summary>
        /// Discards the <c>i</c> first or last elements of the array,
        /// depending on the value of <c>front</c>.
        /// For example, if the array contains the elements 1,2,3,4, invoking
        /// <c>discardExtremeElements(2,false)</c> will cause the last two elements
        /// to be discarded, leaving 1,2 in the array.
        /// For example, if the array contains the elements 1,2,3,4, invoking
        /// <c>discardExtremeElements(2,true)</c> will cause the first two elements
        /// to be discarded, leaving 3,4 in the array.
        /// Throws illegalArgumentException
        /// if i exceeds numElements.
        /// </summary>
        /// <param name="i">the number of elements to discard from the front/end of the array
        /// </param>
        /// <param name="front">true if elements are to be discarded from the front
        /// of the array, false if elements are to be discarded from the end
        /// of the array</param>
        /// <exception cref="MathIllegalArgumentException"> if i is greater than numElements.
        /// </exception>
        private void discardExtremeElements(int i, Boolean front)
        {
            if (i > numElements)
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("TOO_MANY_ELEMENTS_TO_DISCARD_FROM_ARRAY"),
                        i, numElements);
            }
            else if (i < 0)
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("CANNOT_DISCARD_NEGATIVE_NUMBER_OF_ELEMENTS"),
                        i);
            }
            else
            {
                // "Subtract" this number of discarded from numElements
                numElements -= i;
                if (front)
                {
                    startIndex += i;
                }
            }
            if (shouldContract())
            {
                contract();
            }
        }

        /// <summary>
        /// Expands the internal storage array using the expansion factor.
        /// <para>
        /// if <c>expansionMode</c> is set to MULTIPLICATIVE_MODE,
        /// the new array size will be <c>internalArray.length * expansionFactor.</c>
        /// If <c>expansionMode</c> is set to ADDITIVE_MODE,  the length
        /// after expansion will be <c>internalArray.length + expansionFactor</c>
        /// </para>
        /// </summary>
        protected void expand()
        {
            lock (this)
            {
                // notice the use of FastMath.ceil(), this guarantees that we will always
                // have an array of at least currentSize + 1.   Assume that the
                // current initial capacity is 1 and the expansion factor
                // is 1.000000000000000001.  The newly calculated size will be
                // rounded up to 2 after the multiplication is performed.
                int newSize = 0;
                if (expansionMode == ExpansionMode.MULTIPLICATIVE)
                {
                    newSize = (int)FastMath.ceil(internalArray.Length * expansionFactor);
                }
                else
                {
                    newSize = (int)(internalArray.Length + FastMath.round(expansionFactor));
                }
                double[] tempArray = new double[newSize];

                // Copy and swap
                Array.Copy(internalArray, 0, tempArray, 0, internalArray.Length);
                internalArray = tempArray;
            }
        }

        /// <summary>
        /// Expands the internal storage array to the specified size.
        /// </summary>
        /// <param name="size">Size of the new internal storage array.</param>
        private void expandTo(int size)
        {
            lock (this)
            {
                double[] tempArray = new double[size];
                // Copy and swap
                Array.Copy(internalArray, 0, tempArray, 0, internalArray.Length);
                internalArray = tempArray;
            }
        }

        /// <summary>
        /// The contraction criteria defines when the internal array will contract
        /// to store only the number of elements in the element array.
        /// If  the <c>expansionMode</c> is <c>MULTIPLICATIVE_MODE</c>,
        /// contraction is triggered when the ratio between storage array length
        /// and <c>numElements</c> exceeds <c>contractionFactor</c>.
        /// If the <c>expansionMode</c> is <c>ADDITIVE_MODE</c>, the
        /// number of excess storage locations is compared to
        /// <c>contractionFactor.</c> 
        /// </summary>
        /// <returns>the contraction criteria used to reclaim memory.</returns>
        [Obsolete("Please use getContractionCriterion() instead")]
        public float getContractionCriteria()
        {
            return (float)getContractionCriterion();
        }

        /// <summary>
        /// The contraction criterion defines when the internal array will contract
        /// to store only the number of elements in the element array.
        /// If  the <c>expansionMode</c> is <c>MULTIPLICATIVE_MODE</c>,
        /// contraction is triggered when the ratio between storage array length
        /// and <c>numElements</c> exceeds <c>contractionFactor</c>.
        /// If the <c>expansionMode</c> is <c>ADDITIVE_MODE</c>, the
        /// number of excess storage locations is compared to
        /// <c>contractionFactor.</c>
        /// </summary>
        /// <returns>the contraction criterion used to reclaim memory.</returns>
        public double getContractionCriterion()
        {
            return contractionCriterion;
        }

        /// <summary>
        /// Returns the element at the specified index
        /// </summary>
        /// <param name="index">index to fetch a value from</param>
        /// <returns>stored at the specified index</returns>
        /// <exception cref="ArrayIndexOutOfBoundsException"> if <c>index</c> is less than
        /// zero or is greater than <c>getNumElements() - 1</c>.</exception>
        public double getElement(int index)
        {
            if (index >= numElements)
            {
                throw new IndexOutOfRangeException(String.Format("Array out of bounds at index {0}", index));
            }
            else if (index >= 0)
            {
                return internalArray[startIndex + index];
            }
            else
            {
                throw new IndexOutOfRangeException(String.Format("Array out of bounds at index {0}", index));
            }
        }

        /// <summary>
        /// Returns a double array containing the elements of this
        /// <c>ResizableArray</c>.  This method returns a copy, not a
        /// reference to the underlying array, so that changes made to the returned
        /// array have no effect on this <c>ResizableArray.</c>
        /// </summary>
        /// <returns>the double array.</returns>
        public double[] getElements()
        {
            double[] elementArray = new double[numElements];
            Array.Copy(internalArray, startIndex, elementArray, 0, numElements);
            return elementArray;
        }

        /// <summary>
        /// The expansion factor controls the size of a new array when an array
        /// needs to be expanded.  The <c>expansionMode</c>
        /// determines whether the size of the array is multiplied by the
        /// <c>expansionFactor</c> (MULTIPLICATIVE_MODE) or if
        /// the expansion is additive (ADDITIVE_MODE -- <c>expansionFactor</c>
        /// storage locations added).  The default <c>expansionMode</c> is
        /// MULTIPLICATIVE_MODE and the default <c>expansionFactor</c>
        /// is 2.0. 
        /// </summary>
        /// <returns>the expansion factor of this expandable double array</returns>
        [Obsolete]
        public float getExpansionFactor()
        {
            return (float)expansionFactor;
        }

        /// <summary>
        /// The expansion mode determines whether the internal storage
        /// array grows additively or multiplicatively when it is expanded.
        /// </summary>
        /// <returns>the expansion mode.</returns>
        [Obsolete]
        public int getExpansionMode()
        {
            switch (expansionMode)
            {
                case ExpansionMode.MULTIPLICATIVE:
                    return MULTIPLICATIVE_MODE;
                case ExpansionMode.ADDITIVE:
                    return ADDITIVE_MODE;
                default:
                    throw new MathInternalError(); // Should never happen.
            }
        }

        /// <summary>
        /// Notice the package scope on this method.   This method is simply here
        /// for the JUnit test, it allows us check if the expansion is working
        /// properly after a number of expansions.  This is not meant to be a part
        /// of the public interface of this class. 
        /// </summary>
        /// <returns>the length of the internal storage array.</returns>
        [Obsolete("Please use getCapacity() instead")]
        internal int getInternalLength()
        {
            lock (this)
            {
                return internalArray.Length;
            }
        }

        /// <summary>
        /// Gets the currently allocated size of the internal data structure used
        /// for storing elements.
        /// This is not to be confused with <see cref="getNumElements()"/> the number of
        /// elements actually stored}.
        /// </summary>
        /// <returns>the length of the internal array.</returns>
        public int getCapacity()
        {
            return internalArray.Length;
        }

        /// <summary>
        /// Returns the number of elements currently in the array.  Please note
        /// that this is different from the length of the internal storage array.
        /// </summary>
        /// <returns>the number of elements.</returns>
        public int getNumElements()
        {
            return numElements;
        }

        /// <summary>
        /// Returns the internal storage array.  Note that this method returns
        /// a reference to the internal storage array, not a copy, and to correctly
        /// address elements of the array, the <c>startIndex</c> is
        /// required (available via the <see cref="start"/> method).  This method should
        /// only be used in cases where copying the internal array is not practical.
        /// The <see cref="getElements"/> method should be used in all other cases.
        /// </summary>
        /// <returns>the internal storage array used by this object</returns>
        [Obsolete]
        public double[] getInternalValues()
        {
            lock (this)
            {
                return internalArray;
            }
        }

        /// <summary>
        /// Provides direct access to the internal storage array.
        /// Please note that this method returns a reference to this object's
        /// storage array, not a copy.
        /// <para/>
        /// To correctly address elements of the array, the "start index" is
        /// required (available via the <see cref="getStartIndex()">getStartIndex</see>
        /// method.
        /// <para/>
        /// This method should only be used to avoid copying the internal array.
        /// The returned value <em>must</em> be used for reading only; other
        /// uses could lead to this object becoming inconsistent.
        /// <para/>
        /// The <see cref="getElements"/> method has no such limitation since it
        /// returns a copy of this array's addressable elements.
        /// </summary>
        /// <returns>the internal storage array used by this object.</returns>
        protected double[] getArrayRef()
        {
            lock (this)
            {
                return internalArray;
            }
        }

        /// <summary>
        /// Returns the "start index" of the internal array.
        /// This index is the position of the first addressable element in the
        /// internal storage array.
        /// The addressable elements in the array are at indices contained in
        /// the interval [<see cref="getStartIndex()"/>,
        /// <see cref="getStartIndex()"/> + <see cref="getNumElements()"/> - 1].
        /// </summary>
        /// <returns>the start index.</returns>
        protected int getStartIndex()
        {
            return startIndex;
        }

        /// <summary>
        /// Sets the contraction criteria. 
        /// </summary>
        /// <param name="contractionCriteria">contraction criteria</param>
        /// <exception cref="MathIllegalArgumentException"> if the contractionCriteria is less
        /// than the expansionCriteria.</exception>
        [Obsolete]
        public void setContractionCriteria(float contractionCriteria)
        {
            checkContractExpand(contractionCriteria, getExpansionFactor());
            lock (this)
            {
                this.contractionCriterion = contractionCriteria;
            }
        }

        /// <summary>
        /// Performs an operation on the addressable elements of the array.
        /// </summary>
        /// <param name="f">Function to be applied on this array.</param>
        /// <returns>the result.</returns>
        public double compute(MathArrays.Function f)
        {
            double[] array;
            int start;
            int num;
            lock (this)
            {
                array = internalArray;
                start = startIndex;
                num = numElements;
            }
            return f.evaluate(array, start, num);
        }

        /// <summary>
        /// Sets the element at the specified index.  If the specified index is greater than
        /// <c>getNumElements() - 1</c>, the <c>numElements</c> property
        /// is increased to <c>index +1</c> and additional storage is allocated
        /// (if necessary) for the new element and all  (uninitialized) elements
        /// between the new element and the previous end of the array).
        /// </summary>
        /// <param name="index">index to store a value in</param>
        /// <param name="value">value to store at the specified index</param>
        /// <exception cref="ArrayIndexOutOfBoundsException"> if <c>index < 0</c>.</exception>
        public void setElement(int index, double value)
        {
            lock (this)
            {
                if (index < 0)
                {
                    throw new IndexOutOfRangeException(String.Format("Array out of bounds at index {0}", index));
                }
                if (index + 1 > numElements)
                {
                    numElements = index + 1;
                }
                if ((startIndex + index) >= internalArray.Length)
                {
                    expandTo(startIndex + (index + 1));
                }
                internalArray[startIndex + index] = value;
            }
        }

        /// <summary>
        /// Sets the expansionFactor. Throws IllegalArgumentException if the
        /// the following conditions are not met:
        /// <list type="bullet">
        /// <item><c>expansionFactor > 1</c></item>
        /// <item><c>contractionFactor >= expansionFactor</c></item>
        /// </list>
        /// </summary>
        /// <param name="expansionFactor">the new expansion factor value.</param>
        /// <exception cref="MathIllegalArgumentException"> if expansionFactor is <= 1
        /// or greater than contractionFactor</exception>
        [Obsolete]
        public void setExpansionFactor(float expansionFactor)
        {
            checkContractExpand(getContractionCriterion(), expansionFactor);
            // The check above verifies that the expansion factor is > 1.0;
            lock (this)
            {
                this.expansionFactor = expansionFactor;
            }
        }

        /// <summary>
        /// Sets the <c>expansionMode</c>. The specified value must be one of
        /// ADDITIVE_MODE, MULTIPLICATIVE_MODE.
        /// </summary>
        /// <param name="expansionMode">The expansionMode to set.</param>
        /// <exception cref="MathIllegalArgumentException"> if the specified mode value
        /// is not valid.</exception>
        [Obsolete]
        public void setExpansionMode(int expansionMode)
        {
            if (expansionMode != MULTIPLICATIVE_MODE &&
                expansionMode != ADDITIVE_MODE)
            {
                throw new MathIllegalArgumentException(new LocalizedFormats("UNSUPPORTED_EXPANSION_MODE"), expansionMode, MULTIPLICATIVE_MODE, "MULTIPLICATIVE_MODE", ADDITIVE_MODE, "ADDITIVE_MODE");
            }
            lock (this)
            {
                if (expansionMode == MULTIPLICATIVE_MODE)
                {
                    setExpansionMode(ExpansionMode.MULTIPLICATIVE);
                }
                else if (expansionMode == ADDITIVE_MODE)
                {
                    setExpansionMode(ExpansionMode.ADDITIVE);
                }
            }
        }

        /// <summary>
        /// Sets the <see cref="ExpansionMode">expansion mode</see>. 
        /// </summary>
        /// <param name="expansionMode">Expansion mode to use for resizing the array.</param>
        [Obsolete]
        public void setExpansionMode(ExpansionMode expansionMode)
        {
            this.expansionMode = expansionMode;
        }

        /// <summary>
        /// Sets the initial capacity.  Should only be invoked by constructors.
        /// </summary>
        /// <param name="initialCapacity">initial capacity of the array</param>
        /// <exception cref="MathIllegalArgumentException"> if <c>initialCapacity</c> is not
        /// positive.</exception>
        [Obsolete]
        protected void setInitialCapacity(int initialCapacity) { /* Body removed in 3.1. */ }

        /// <summary>
        /// This function allows you to control the number of elements contained
        /// in this array, and can be used to "throw out" the last n values in an
        /// array. This function will also expand the internal array as needed.
        /// </summary>
        /// <param name="i">a new number of elements</param>
        /// <exception cref="MathIllegalArgumentException"> if <c>i</c> is negative.</exception>
        public void setNumElements(int i)
        {
            lock (this)
            {
                // If index is negative thrown an error.
                if (i < 0)
                {
                    throw new MathIllegalArgumentException(new LocalizedFormats("INDEX_NOT_POSITIVE"),
                            i);
                }

                // Test the new num elements, check to see if the array needs to be
                // expanded to accommodate this new number of elements.
                int newSize = startIndex + i;
                if (newSize > internalArray.Length)
                {
                    expandTo(newSize);
                }

                // Set the new number of elements to new value.
                numElements = i;
            }
        }

        /// <summary>
        /// Returns true if the internal storage array has too many unused
        /// storage positions. 
        /// </summary>
        /// <returns>true if array satisfies the contraction criteria</returns>
        private Boolean shouldContract()
        {
            lock (this)
            {
                if (expansionMode == ExpansionMode.MULTIPLICATIVE)
                {
                    return (internalArray.Length / ((float)numElements)) > contractionCriterion;
                }
                else
                {
                    return (internalArray.Length - numElements) > contractionCriterion;
                }
            }
        }

        /// <summary>
        /// Returns the starting index of the internal array.  The starting index is
        /// the position of the first addressable element in the internal storage
        /// array.  The addressable elements in the array are <c>
        /// internalArray[startIndex],...,internalArray[startIndex + numElements -1]
        /// </c>
        /// </summary>
        /// <returns>the starting index.</returns>
        [Obsolete]
        public int start()
        {
            lock (this)
            {
                return startIndex;
            }
        }

        /// <summary>
        /// <para>Copies source to dest, copying the underlying data, so dest is
        /// a new, independent copy of source.  Does not contract before
        /// the copy.</para>
        /// <para>Obtains synchronization locks on both source and dest
        /// (in that order) before performing the copy.</para>
        /// <para>Neither source nor dest may be null; otherwise a <see 
        /// cref="NullArgumentException"/>
        /// is thrown</para>
        /// </summary>
        /// <param name="source">ResizableDoubleArray to copy</param>
        /// <param name="dest">ResizableArray to replace with a copy of the source array</param>
        /// <exception cref="NullArgumentException"> if either source or dest is null</exception>
        public static void copy(ResizableDoubleArray source, ResizableDoubleArray dest)
        {
            MathUtils.checkNotNull(source);
            MathUtils.checkNotNull(dest);
            lock (source)
            {
                lock (dest)
                {
                    dest.contractionCriterion = source.contractionCriterion;
                    dest.expansionFactor = source.expansionFactor;
                    dest.expansionMode = source.expansionMode;
                    dest.internalArray = new double[source.internalArray.Length];
                    Array.Copy(source.internalArray, 0, dest.internalArray, 0, dest.internalArray.Length);
                    dest.numElements = source.numElements;
                    dest.startIndex = source.startIndex;
                }
            }
        }

        /// <summary>
        /// Returns a copy of the ResizableDoubleArray.  Does not contract before
        /// the copy, so the returned object is an exact copy of this. 
        /// </summary>
        /// <returns>a new ResizableDoubleArray with the same data and configuration
        /// properties as this</returns>
        public ResizableDoubleArray copy()
        {
            lock (this)
            {
                ResizableDoubleArray result = new ResizableDoubleArray();
                copy(this, result);
                return result;
            }
        }

        /// <summary>
        /// Returns true iff object is a ResizableDoubleArray with the same properties
        /// as this and an identical internal storage array.
        /// </summary>
        /// <param name="obj">object to be compared for equality with this</param>
        /// <returns>true iff obj is a ResizableDoubleArray with the same data and
        /// properties as this</returns>
        public override Boolean Equals(Object obj)
        {
            if (obj == this)
            {
                return true;
            }
            if (obj is ResizableDoubleArray)
            {
                return false;
            }
            lock (this)
            {
                lock (obj)
                {
                    Boolean result = true;
                    ResizableDoubleArray other = (ResizableDoubleArray)obj;
                    result = result && (other.contractionCriterion == contractionCriterion);
                    result = result && (other.expansionFactor == expansionFactor);
                    result = result && (other.expansionMode == expansionMode);
                    result = result && (other.numElements == numElements);
                    result = result && (other.startIndex == startIndex);
                    if (!result)
                    {
                        return false;
                    }
                    else
                    {
                        return Array.Equals(internalArray, other.internalArray);
                    }
                }
            }
        }

        /// <summary>
        /// Returns a hash code consistent with equals.
        /// </summary>
        /// <returns>the hash code representing this <c>ResizableDoubleArray</c>.</returns>
        public override int GetHashCode()
        {
            lock (this)
            {
                int[] hashData = new int[6];
                hashData[0] = expansionFactor.GetHashCode();
                hashData[1] = contractionCriterion.GetHashCode();
                hashData[2] = expansionMode.GetHashCode();
                hashData[3] = internalArray.GetHashCode();
                hashData[4] = numElements;
                hashData[5] = startIndex;
                return hashData.GetHashCode();
            }
        }
    }
}
