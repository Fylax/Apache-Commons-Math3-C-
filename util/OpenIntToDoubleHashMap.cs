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
using System.Collections.Generic;

namespace Math3.util
{
    /// <summary>
    /// Open addressed map from int to double.
    /// <para>This class provides a dedicated map from integers to doubles with a
    /// much smaller memory overhead than standard <c>HashMap</c>.</para>
    /// <para>This class is not synchronized. The specialized iterators returned by
    /// <see cref="iterator()"/> are fail-fast: they throw a
    /// <c>ConcurrentModificationException</c> when they detect the map has been
    /// modified during iteration.</para>
    /// </summary>
    public class OpenIntToDoubleHashMap
    {
        /// <summary>
        /// Status indicator for free table entries.
        /// </summary>
        protected const byte FREE = 0;

        /// <summary>
        /// Status indicator for full table entries.
        /// </summary>
        protected const byte FULL = 1;

        /// <summary>
        /// Status indicator for removed table entries.
        /// </summary>
        protected const byte REMOVED = 2;

        /// <summary>
        /// Load factor for the map.
        /// </summary>
        private const float LOAD_FACTOR = 0.5f;

        /// <summary>
        /// Default starting size.
        /// <para>This must be a power of two for bit mask to work properly.</para>
        /// </summary>
        private const int DEFAULT_EXPECTED_SIZE = 16;

        /// <summary>
        /// Multiplier for size growth when map fills up.
        /// <para>This must be a power of two for bit mask to work properly. </para>
        /// </summary>
        private const int RESIZE_MULTIPLIER = 2;

        /// <summary>
        /// Number of bits to perturb the index when probing for collision resolution.
        /// </summary>
        private const int PERTURB_SHIFT = 5;

        /// <summary>
        /// Keys table.
        /// </summary>
        private int[] keys;

        /// <summary>
        /// Values table.
        /// </summary>
        private double[] values;

        /// <summary>
        /// States table.
        /// </summary>
        private byte[] states;

        /// <summary>
        /// Return value for missing entries.
        /// </summary>
        private readonly double missingEntries;

        /// <summary>
        /// Current size of the map.
        /// </summary>
        private int Size;

        /// <summary>
        /// Bit mask for hash values.
        /// </summary>
        private int mask;

        /// <summary>
        /// Modifications count.
        /// </summary>
        private int count;

        /// <summary>
        /// Build an empty map with default size and using NaN for missing entries.
        /// </summary>
        public OpenIntToDoubleHashMap() : this(DEFAULT_EXPECTED_SIZE, Double.NaN) { }

        /// <summary>
        /// Build an empty map with default size
        /// </summary>
        /// <param name="missingEntries">value to return when a missing entry is fetched</param>
        public OpenIntToDoubleHashMap(double missingEntries) : this(DEFAULT_EXPECTED_SIZE, missingEntries) { }

        /// <summary>
        /// Build an empty map with specified size and using NaN for missing entries.
        /// </summary>
        /// <param name="expectedSize">expected number of elements in the map</param>
        public OpenIntToDoubleHashMap(int expectedSize) : this(expectedSize, Double.NaN) { }

        /// <summary>
        /// Build an empty map with specified size.
        /// </summary>
        /// <param name="expectedSize">expected number of elements in the map</param>
        /// <param name="missingEntries">value to return when a missing entry is fetched</param>
        public OpenIntToDoubleHashMap(int expectedSize, double missingEntries)
        {
            int capacity = computeCapacity(expectedSize);
            keys = new int[capacity];
            values = new double[capacity];
            states = new byte[capacity];
            this.missingEntries = missingEntries;
            mask = capacity - 1;
        }

        /// <summary>
        /// Copy constructor.
        /// </summary>
        /// <param name="source">map to copy</param>
        public OpenIntToDoubleHashMap(OpenIntToDoubleHashMap source)
        {
            int length = source.keys.Length;
            keys = new int[length];
            Array.Copy(source.keys, 0, keys, 0, length);
            values = new double[length];
            Array.Copy(source.values, 0, values, 0, length);
            states = new byte[length];
            Array.Copy(source.states, 0, states, 0, length);
            missingEntries = source.missingEntries;
            Size = source.Size;
            mask = source.mask;
            count = source.count;
        }

        /// <summary>
        /// Compute the capacity needed for a given size.
        /// </summary>
        /// <param name="expectedSize">expected size of the map</param>
        /// <returns>to use for the specified size</returns>
        private static int computeCapacity(int expectedSize)
        {
            if (expectedSize == 0)
            {
                return 1;
            }
            int capacity = (int)FastMath.ceil(expectedSize / LOAD_FACTOR);
            int powerOfTwo = capacity.HighestOneBit();
            if (powerOfTwo == capacity)
            {
                return capacity;
            }
            return nextPowerOfTwo(capacity);
        }

        /// <summary>
        /// Find the smallest power of two greater than the input value
        /// </summary>
        /// <param name="i">input value</param>
        /// <returns>power of two greater than the input value</returns>
        private static int nextPowerOfTwo(int i)
        {
            return i.HighestOneBit() << 1;
        }

        /// <summary>
        /// Get the stored value associated with the given key
        /// </summary>
        /// <param name="key">key associated with the data</param>
        /// <returns>data associated with the key</returns>
        public double get(int key)
        {

            int hash = hashOf(key);
            int index = hash & mask;
            if (containsKey(key, index))
            {
                return values[index];
            }

            if (states[index] == FREE)
            {
                return missingEntries;
            }

            int j = index;
            for (int Perturb = perturb(hash); states[index] != FREE; Perturb >>= PERTURB_SHIFT)
            {
                j = probe(Perturb, j);
                index = j & mask;
                if (containsKey(key, index))
                {
                    return values[index];
                }
            }

            return missingEntries;

        }

        /// <summary>
        /// Check if a value is associated with a key.
        /// </summary>
        /// <param name="key">key to check</param>
        /// <returns>true if a value is associated with key</returns>
        public Boolean containsKey(int key)
        {

            int hash = hashOf(key);
            int index = hash & mask;
            if (containsKey(key, index))
            {
                return true;
            }

            if (states[index] == FREE)
            {
                return false;
            }

            int j = index;
            for (int Perturb = perturb(hash); states[index] != FREE; Perturb >>= PERTURB_SHIFT)
            {
                j = probe(Perturb, j);
                index = j & mask;
                if (containsKey(key, index))
                {
                    return true;
                }
            }

            return false;

        }

        /// <summary>
        /// Get an iterator over map elements.
        /// <para>The specialized iterators returned are fail-fast: they throw a
        /// <c>ConcurrentModificationException</c> when they detect the map
        /// has been modified during iteration.</para> 
        /// </summary>
        /// <returns>iterator over the map elements</returns>
        public Iterator iterator()
        {
            return new Iterator(this);
        }

        /// <summary>
        /// Perturb the hash for starting probing.
        /// </summary>
        /// <param name="hash">initial hash</param>
        /// <returns>perturbed hash</returns>
        private static int perturb(int hash)
        {
            return hash & 0x7fffffff;
        }

        /// <summary>
        /// Find the index at which a key should be inserted
        /// </summary>
        /// <param name="key">key to lookup</param>
        /// <returns>index at which key should be inserted</returns>
        private int findInsertionIndex(int key)
        {
            return findInsertionIndex(keys, states, key, mask);
        }

        /// <summary>
        /// Find the index at which a key should be inserted
        /// </summary>
        /// <param name="keys">keys table</param>
        /// <param name="states">states table</param>
        /// <param name="key">key to lookup</param>
        /// <param name="mask">bit mask for hash values</param>
        /// <returns>index at which key should be inserted</returns>
        private static int findInsertionIndex(int[] keys, byte[] states, int key, int mask)
        {
            int hash = hashOf(key);
            int index = hash & mask;
            if (states[index] == FREE)
            {
                return index;
            }
            else if (states[index] == FULL && keys[index] == key)
            {
                return changeIndexSign(index);
            }

            int Perturb = perturb(hash);
            int j = index;
            if (states[index] == FULL)
            {
                while (true)
                {
                    j = probe(Perturb, j);
                    index = j & mask;
                    Perturb >>= PERTURB_SHIFT;

                    if (states[index] != FULL || keys[index] == key)
                    {
                        break;
                    }
                }
            }

            if (states[index] == FREE)
            {
                return index;
            }
            else if (states[index] == FULL)
            {
                // due to the loop exit condition,
                // if (states[index] == FULL) then keys[index] == key
                return changeIndexSign(index);
            }

            int firstRemoved = index;
            while (true)
            {
                j = probe(Perturb, j);
                index = j & mask;

                if (states[index] == FREE)
                {
                    return firstRemoved;
                }
                else if (states[index] == FULL && keys[index] == key)
                {
                    return changeIndexSign(index);
                }

                Perturb >>= PERTURB_SHIFT;

            }

        }

        /// <summary>
        /// Compute next probe for collision resolution
        /// </summary>
        /// <param name="perturb">perturbed hash</param>
        /// <param name="j">previous probe</param>
        /// <returns>next probe</returns>
        private static int probe(int perturb, int j)
        {
            return (j << 2) + j + perturb + 1;
        }

        /// <summary>
        /// Change the index sign
        /// </summary>
        /// <param name="index">initial index</param>
        /// <returns>changed index</returns>
        private static int changeIndexSign(int index)
        {
            return -index - 1;
        }

        /// <summary>
        /// Get the number of elements stored in the map.
        /// </summary>
        /// <returns>number of elements stored in the map</returns>
        public int size()
        {
            return Size;
        }


        /// <summary>
        /// Remove the value associated with a key.
        /// </summary>
        /// <param name="key">key to which the value is associated</param>
        /// <returns>removed value</returns>
        public double remove(int key)
        {
            int hash = hashOf(key);
            int index = hash & mask;
            if (containsKey(key, index))
            {
                return doRemove(index);
            }

            if (states[index] == FREE)
            {
                return missingEntries;
            }

            int j = index;
            for (int Perturb = perturb(hash); states[index] != FREE; Perturb >>= PERTURB_SHIFT)
            {
                j = probe(Perturb, j);
                index = j & mask;
                if (containsKey(key, index))
                {
                    return doRemove(index);
                }
            }

            return missingEntries;

        }

        /// <summary>
        /// Check if the tables contain an element associated with specified key
        /// at specified index.
        /// </summary>
        /// <param name="key">key to check</param>
        /// <param name="index">index to check</param>
        /// <returns>true if an element is associated with key at index</returns>
        private Boolean containsKey(int key, int index)
        {
            return (key != 0 || states[index] == FULL) && keys[index] == key;
        }

        /// <summary>
        /// Remove an element at specified index.
        /// </summary>
        /// <param name="index">index of the element to remove</param>
        /// <returns>removed value</returns>
        private double doRemove(int index)
        {
            keys[index] = 0;
            states[index] = REMOVED;
            double previous = values[index];
            values[index] = missingEntries;
            --Size;
            ++count;
            return previous;
        }

        /// <summary>
        /// Put a value associated with a key in the map.
        /// </summary>
        /// <param name="key">key to which value is associated</param>
        /// <param name="value">value to put in the map</param>
        /// <returns>previous value associated with the key</returns>
        public double put(int key, double value)
        {
            int index = findInsertionIndex(key);
            double previous = missingEntries;
            Boolean newMapping = true;
            if (index < 0)
            {
                index = changeIndexSign(index);
                previous = values[index];
                newMapping = false;
            }
            keys[index] = key;
            states[index] = FULL;
            values[index] = value;
            if (newMapping)
            {
                ++Size;
                if (shouldGrowTable())
                {
                    growTable();
                }
                ++count;
            }
            return previous;

        }

        /// <summary>
        /// Grow the tables.
        /// </summary>
        private void growTable()
        {
            int oldLength = states.Length;
            int[] oldKeys = keys;
            double[] oldValues = values;
            byte[] oldStates = states;

            int newLength = RESIZE_MULTIPLIER * oldLength;
            int[] newKeys = new int[newLength];
            double[] newValues = new double[newLength];
            byte[] newStates = new byte[newLength];
            int newMask = newLength - 1;
            for (int i = 0; i < oldLength; ++i)
            {
                if (oldStates[i] == FULL)
                {
                    int key = oldKeys[i];
                    int index = findInsertionIndex(newKeys, newStates, key, newMask);
                    newKeys[index] = key;
                    newValues[index] = oldValues[i];
                    newStates[index] = FULL;
                }
            }

            mask = newMask;
            keys = newKeys;
            values = newValues;
            states = newStates;

        }

        /// <summary>
        /// Check if tables should grow due to increased size.
        /// </summary>
        /// <returns>true if  tables should grow</returns>
        private Boolean shouldGrowTable()
        {
            return Size > (mask + 1) * LOAD_FACTOR;
        }

        /// <summary>
        /// Compute the hash value of a key
        /// </summary>
        /// <param name="key">key to hash</param>
        /// <returns>value of the key</returns>
        private static int hashOf(int key)
        {
            int h = key ^ (unchecked((int)((uint)key >> 20)) ^ (unchecked((int)((uint)key >> 12))));
            return h ^ (unchecked((int)((uint)h >> 7))) ^ (unchecked((int)((uint)h >> 4)));
        }


        /// <summary>
        /// Iterator class for the map.
        /// </summary>
        public class Iterator : IEnumerator<KeyValuePair<int, double>>
        {
            /// <summary>
            /// Reference modification count.
            /// </summary>
            private int referenceCount;

            /// <summary>
            /// Index of current element.
            /// </summary>
            private int current;

            /// <summary>
            /// Index of next element.
            /// </summary>
            private int next;

            private Boolean disposed = false;
            private OpenIntToDoubleHashMap o;
            private KeyValuePair<int, double> currentElem;

            /// <summary>
            /// Simple constructor.
            /// </summary>
            /// <param name="oidhm"></param>
            internal Iterator(OpenIntToDoubleHashMap oidhm)
            {
                this.o = oidhm;
                // preserve the modification count of the map to detect concurrent modifications later
                referenceCount = o.count;

                // initialize current index
                next = -1;
                MoveNext();
            }


            public KeyValuePair<int, double> Current
            {
                get { return currentElem; }
            }

            object System.Collections.IEnumerator.Current
            {
                get { return currentElem; }
            }

            public bool MoveNext()
            {
                if (referenceCount != o.count)
                {
                    return false;
                }

                // advance on step
                current = next;

                // prepare next step
                try
                {
                    while (o.states[++next] != FULL)
                    { // NOPMD
                        // nothing to do
                    }
                }
                catch (IndexOutOfRangeException)
                {
                    next = -2;
                    if (current < 0)
                    {
                        return false;
                    }
                }
                currentElem = new KeyValuePair<int, double>(o.keys[current], o.values[current]);
                return true;
            }

            public void Reset()
            {
                return;
            }

            protected virtual void Dispose(Boolean disposing)
            {
                if(!this.disposed)
                {
                    if(disposing)
                    {
                        o = null;
                        disposed = true;
                    }
                }
            }

            public void Dispose()
            {
                Dispose(true);
                // This object will be cleaned up by the Dispose method.
                // Therefore, you should call GC.SupressFinalize to
                // take this object off the finalization queue
                // and prevent finalization code for this object
                // from executing a second time.
                GC.SuppressFinalize(this);
            }
        }
    }

    public static class Extensions
    {
        public static int HighestOneBit(this int number)
        {
            uint i = (uint)number;
            i |= (i >> 1);
            i |= (i >> 2);
            i |= (i >> 4);
            i |= (i >> 8);
            i |= (i >> 16);
            i -= (i >> 1);
            return (int)i;
        }
    }

}
