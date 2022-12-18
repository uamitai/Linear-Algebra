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

        public InnerProductSpace<V, F> GramSchmidt()
        {
            InnerProductSpace<V, F> GS = Clone();
            int dimension = Dimension();
            V v;
            F scalar;
            foreach(V vector in this) { GS.Add(vector); }
            for (int i = 0; i < dimension; i++)
            {
                v = GS[i];
                GS[i] = (V)v.Normalize();
                scalar = (F)(v * v).MultInverse();
                for (int j = i + 1; j < dimension; j++)
                {
                    GS[j] = (V)(GS[j] - v.Multiply(GS[j] * v).Multiply(scalar));
                }
            }
            GS.basisMatrix = GS.BasisMatrix();
            return GS;
        }
    }
}
