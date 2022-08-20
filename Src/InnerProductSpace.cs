namespace Linear_Algebra
{
    class InnerProductSpace<V, F> : VectorSpace<V, F> where F : Complex where V : InnerProduct<F>
    {
        public InnerProductSpace(int dim) : base(dim) { }

        public InnerProductSpace(int dim, params V[] vectors) : base(dim, vectors) { }

        public override InnerProductSpace<V, F> Clone()
        {
            return new InnerProductSpace<V, F>(dim);
        }

        public InnerProductSpace<V, F> Normalized()
        {
            InnerProductSpace<V, F> norm = Clone();
            foreach(V vector in this)
            {
                norm.Add((V)vector.Normalize());
            }
            return norm;
        }

        public InnerProductSpace<V, F> GrahamSchmidt()
        {
            InnerProductSpace<V, F > GS = Clone();
            V vector;
            for (int i = 0; i < Dimension(); i++)
            {
                vector = basis[i];
                for (int j = 0; j < i; j++)
                {
                    vector = (V)(vector - basis[j].Multiply((F)(basis[i] * basis[j]).Multiply((basis[j] * basis[j]).MultInverse())));
                }
                GS.Add(vector);
            }
            return GS;
        }
    }
}
