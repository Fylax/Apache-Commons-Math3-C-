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
using System;

namespace Math3.util
{
    /// <summary>
    /// Utility that increments a counter until a maximum is reached, at
    /// which point, the instance will by default throw a
    /// <see cref="MaxCountExceededException"/>.
    /// custom <see cref="MaxCountExceededCallback">callback</see>, in order to e.g.
    /// select which exception must be thrown.
    /// </summary>
    public class Incrementor
    {
        /// <summary>
        /// Upper limit for the counter.
        /// </summary>
        private int maximalCount;

        /// <summary>
        /// Current count.
        /// </summary>
        private int count = 0;

        /// <summary>
        /// Function called at counter exhaustion.
        /// </summary>
        private readonly MaxCountExceededCallback maxCountCallback;

        /// <summary>
        /// Default constructor.
        /// For the new instance to be useful, the maximal count must be set
        /// by calling <see cref="setMaximalCount(int)">setMaximalCount</see>.
        /// </summary>
        public Incrementor() : this(0) { }

        /// <summary>
        /// Defines a maximal count.
        /// </summary>
        /// <param name="max">Maximal count.</param>
        public Incrementor(int max) : this(max, new MaxCountExceededCallbackHelper()) { }

        private class MaxCountExceededCallbackHelper : MaxCountExceededCallback
        {
            /// <inheritdoc/>
            public void trigger(int max)
            {
                throw new MaxCountExceededException<Int32>(max);
            }
        }

        /// <summary>
        /// Defines a maximal count and a callback method to be triggered at
        /// counter exhaustion.
        /// </summary>
        /// <param name="max">Maximal count.</param>
        /// <param name="cb">Function to be called when the maximal count has been reached.</param>
        /// <exception cref="NullArgumentException"> if <c>cb</c> is <c>null</c></exception>
        public Incrementor(int max, MaxCountExceededCallback cb)
        {
            if (cb == null)
            {
                throw new NullArgumentException();
            }
            maximalCount = max;
            maxCountCallback = cb;
        }

        /// <summary>
        /// Sets the upper limit for the counter.
        /// This does not automatically reset the current count to zero (see
        /// <see cref="resetCount()"/>).
        /// </summary>
        /// <param name="max">Upper limit of the counter.</param>
        public void setMaximalCount(int max)
        {
            maximalCount = max;
        }

        /// <summary>
        /// Gets the upper limit of the counter.
        /// </summary>
        /// <returns>the counter upper limit.</returns>
        public int getMaximalCount()
        {
            return maximalCount;
        }

        /// <summary>
        /// Gets the current count.
        /// </summary>
        /// <returns>the current count.</returns>
        public int getCount()
        {
            return count;
        }

        /// <summary>
        /// Checks whether a single increment is allowed.
        /// </summary>
        /// <returns><c>false</c> if the next call to <see cref="incrementCount(int)">
        /// incrementCount</see> will trigger a <c>MaxCountExceededException</c>,
        /// <c>true</c> otherwise.</returns>
        public Boolean canIncrement()
        {
            return count < maximalCount;
        }

        /// <summary>
        /// Performs multiple increments.
        /// See the other <see cref="incrementCount()"/>incrementCount</see> method).
        /// </summary>
        /// <param name="value">Number of increments.</param>
        /// <exception cref="MaxCountExceededException"> at counter exhaustion.</exception>
        public void incrementCount(int value)
        {
            for (int i = 0; i < value; i++)
            {
                incrementCount();
            }
        }

        /// <summary>
        /// Adds one to the current iteration count.
        /// At counter exhaustion, this method will call the
        /// <see cref="MaxCountExceededCallback.trigger(int)">trigger</see> method of the
        /// callback object passed to the
        /// <see cref="Incrementor(int,MaxCountExceededCallback)">constructor</see>.
        /// If not explictly set, a default callback is used that will throw
        /// a <c>MaxCountExceededException</c>.
        /// </summary>
        /// <exception cref="MaxCountExceededException"> at counter exhaustion, unless a
        /// custom <see cref="MaxCountExceededCallback">callback</see> has been set at
        /// construction.</exception>
        public void incrementCount()
        {
            if (++count > maximalCount)
            {
                maxCountCallback.trigger(maximalCount);
            }
        }

        /// <summary>
        /// Resets the counter to 0.
        /// </summary>
        public void resetCount()
        {
            count = 0;
        }

        /// <summary>
        /// Defines a method to be called at counter exhaustion.
        /// The <see cref="trigger(int)">trigger</see> method should usually throw an exception.
        /// </summary>
        public interface MaxCountExceededCallback
        {
            /// <summary>
            /// Function called when the maximal count has been reached.
            /// </summary>
            /// <param name="maximalCount">Maximal count.</param>
            /// <exception cref="MaxCountExceededException"> at counter exhaustion</exception>
            void trigger(int maximalCount);
        }
    }
}
