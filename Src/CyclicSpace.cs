namespace Linear_Algebra
{
    class CyclicSpace<V, F> : VectorSpace<V, F> where F : Field where V : Vector<F>
    {
        private readonly V vector;
        private readonly Transform<V, F> transform;

        public CyclicSpace(V vector, Transform<V, F> transform) : base(vector.Length())
        {
            this.vector = vector;
            this.transform = transform;
            for (int i = 0; i < dim; i++)
            {
                Add(vector);
                vector = transform.ValueOf(vector);
            }
        }
    }
}
