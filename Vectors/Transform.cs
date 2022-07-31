namespace Linear_Algebra
{
    class Transform<V, F> : Vector<F> where F : Field where V : Vector<F>
    {
        public delegate V LinearTransform(V vector);
        public static readonly LinearTransform zero = (vector) => (V)vector.Zero();
        public static readonly LinearTransform identity = (vector) => vector;

        private readonly LinearTransform transform;

        public Transform(LinearTransform transform)
        {
            this.transform = transform;
        }

        public int Length()
        {
            throw new System.NotImplementedException();
        }

        public ColumnVector<F> ToColumnVector()
        {
            throw new System.NotImplementedException();
        }

        public V ValueOf(V vector)
        {
            return transform(vector);
        }

        public Vector<F> Add(Vector<F> vector)
        {
            Transform<V, F> other = vector as Transform<V, F>;
            return new Transform<V, F>(vec => (V)(transform(vec) + other.transform(vec)));
        }

        public Vector<F> Multiply(F scalar)
        {
            return new Transform<V, F>(vector => (V)(scalar * transform(vector)));
        }

        public Vector<F> Zero() { return new Transform<V, F>(zero); }

        public Vector<F> AddInverse()
        {
            return new Transform<V, F>(vector => (V)transform(vector).AddInverse());
        }

        public static Transform<V, F> operator * (Transform<V, F> T, Transform<V, F> S)
        {
            return new Transform<V, F>(vector => T.transform(S.transform(vector)));
        }

        public Transform<V, F> Power(int exp)
        {
            Transform<V, F> res = identity as Transform<V, F>;
            Transform<V, F> pow = this;
            while(exp > 0)
            {
                if(exp % 2 == 1)
                {
                    res *= pow;
                    exp--;
                }
                pow *= pow;
                exp /= 2;
            }
            return res;
        }

        public VectorSpace<V, F> Image(VectorSpace<V, F> vectorSpace)
        {
            VectorSpace<V, F> im = new VectorSpace<V, F>(vectorSpace.dim);
            foreach(V vector in vectorSpace)
            {
                im.Add(transform(vector));
            }
            return im;
        }

        // @pre vectorSpace.Dimension() == vectorSpace.dim
        public SquareMatrix<F> MatrixRep(VectorSpace<V, F> vectorSpace)
        {
            int dim = vectorSpace.dim;
            F[,] rep = new F[dim, dim];
            int col = 0;
            ColumnVector<F> colVec;
            foreach(V vector in vectorSpace)
            {
                colVec = vectorSpace.CoordinatesOf(transform(vector));
                for (int row = 0; row < dim; row++)
                {
                    rep[row, col] = colVec[row];
                }
                col++;
            }
            return new SquareMatrix<F>(rep);
        }

        public VectorSpace<V, F> Kernel(VectorSpace<V, F> vectorSpace)
        {
            VectorSpace<V, F> kernel = new VectorSpace<V, F>(vectorSpace.dim);
            VectorSpace<ColumnVector<F>, F> nullSpace = MatrixRep(vectorSpace).NullSpace();
            foreach(ColumnVector<F> vector in nullSpace)
            {
                kernel.Add(vectorSpace.LinearCombination(vector));
            }
            return kernel;
        }
    }
}
