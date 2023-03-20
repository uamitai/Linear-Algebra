namespace Linear_Algebra
{
    class QuotientSpace<V, F> : Vector<F> where F : Field where V : Vector<F>
    {
        public readonly V vector;
        public readonly VectorSpace<V, F> subSpace;

        public QuotientSpace(V vector, VectorSpace<V, F> subSpace)
        {
            this.vector = vector;
            this.subSpace = subSpace;
        }

        public Vector<F> Add(Vector<F> vector)
        {
            QuotientSpace<V, F> other = vector as QuotientSpace<V, F>;
            return new QuotientSpace<V, F>((V)(this.vector + other.vector), subSpace);
        }

        public Vector<F> AddInverse()
        {
            return new QuotientSpace<V, F>((V)vector.AddInverse(), subSpace);
        }

        public bool Contains(V other)
        {
            return subSpace.Contains((V)(vector - other));
        }

        public int Length()
        {
            return vector.Length();
        }

        public Vector<F> Multiply(F scalar)
        {
            return new QuotientSpace<V, F>((V)(scalar * vector), subSpace);
        }

        public ColumnVector<F> ToColumnVector()
        {
            return vector.ToColumnVector();
        }

        public Vector<F> Zero()
        {
            return vector.Zero();
        }

        public override string ToString()
        {
            return string.Format("{0}\n+\n\n{1}\n", vector.ToString(), subSpace.ToString());
        }
    }
}
