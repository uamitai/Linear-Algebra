using System.Collections;
using System.Collections.Generic;


namespace Linear_Algebra
{
    class ColumnVector<F> : AbstractMatrix<F>, IEnumerable<F> where F : Field
    {
        public readonly int length;

        // @pre no null entries
        public ColumnVector(params F[] entries) : base(entries.Length, 1)
        {
            length = entries.Length;
            for (int i = 0; i < length; i++)
            {
                base.entries[i, 0] = entries[i];
            }
        }

        // @pre no null entries
        public ColumnVector(F[,] entries) : base(entries)
        {
            int cols = entries.GetLength(1);
            if(cols != 1) { throw new MatrixSizeException(string.Format(
                "Invalid width for column vector: {0}", cols)); }

            length = entries.Length;
        }

        // @pre 0 <= index < length
        public F this[int index]
        {
            get { return entries[index, 0]; }
            set { entries[index, 0] = value; }
        }

        public override ColumnVector<F> Clone()
        {
            return new ColumnVector<F>(entries);
        }

        public override ColumnVector<F> ToColumnVector() { return this; }

        public IEnumerator<F> GetEnumerator()
        {
            foreach(F item in entries)
            {
                yield return item;
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
