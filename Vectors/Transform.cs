namespace Linear_Algebra
{
    class Transform<V, W, F> : Vector<F> where F : Field where V : Vector<F> where W : Vector<F>
    {
        public delegate W LinearTransform(V vector);

        protected readonly LinearTransform transform;
        private readonly VectorSpace<V, F> domain;
        private readonly VectorSpace<W, F> range;
        public readonly Matrix<F> matrixRep;

        public Transform(VectorSpace<V, F> domain, VectorSpace<W, F> range, LinearTransform transform)
        {
            this.domain = domain;
            this.range = range;
            this.transform = transform;
            this.matrixRep = GetMatrixRep();
        }

        public int Length()
        {
            return domain.dim * range.dim;
        }

        public ColumnVector<F> ToColumnVector()
        {
            throw new System.NotImplementedException();
        }

        // @pre domain.Contains(vector)
        public W ValueOf(V vector)
        {
            return transform(vector);
        }

        public virtual Vector<F> Add(Vector<F> vector)
        {
            Transform<V, W, F> other = vector as Transform<V, W, F>;
            return new Transform<V, W, F>(domain, range,
                vec => (W)(transform(vec) + other.transform(vec)));
        }

        public virtual Vector<F> Multiply(F scalar)
        {
            return new Transform<V, W, F>(domain, range,
                vector => (W)(scalar * transform(vector)));
        }

        public virtual Vector<F> Zero()
        {
            return new Transform<V, W, F>(domain, range,
                vector => (W)transform(vector).Zero());
        }

        public virtual Vector<F> AddInverse()
        {
            return new Transform<V, W, F>(domain, range,
                vector => (W)transform(vector).AddInverse());
        }

        public VectorSpace<W, F> Image()
        {
            return Image(domain);
        }

        public VectorSpace<W, F> Image(VectorSpace<V, F> vectorSpace)
        {
            VectorSpace<W, F> im = new VectorSpace<W, F>(range.dim);
            foreach (V vector in vectorSpace)
            {
                im.Add(transform(vector));
            }
            return im;
        }

        private Matrix<F> GetMatrixRep()
        {
            F[,] rep = new F[range.dim, domain.dim];
            int col = 0;
            ColumnVector<F> colVec;
            foreach (V vector in domain)
            {
                colVec = range.CoordinatesOf(transform(vector));
                for (int row = 0; row < range.dim; row++)
                {
                    rep[row, col] = colVec[row];
                }
                col++;
            }
            return new Matrix<F>(rep);
        }

        public VectorSpace<V, F> Kernel()
        {
            VectorSpace<V, F> kernel = new VectorSpace<V, F>(domain.dim);
            VectorSpace<ColumnVector<F>, F> nullSpace = matrixRep.NullSpace();
            foreach (ColumnVector<F> vector in nullSpace)
            {
                kernel.Add(domain.LinearCombination(vector));
            }
            return kernel;
        }
    }

    class Transform<V, F> : Transform<V, V ,F> where F : Field where V : Vector<F>
    {
        public static readonly LinearTransform zero = vector => (V)vector.Zero();
        public static readonly LinearTransform identity = vector => vector;
        private VectorSpace<V, F> vectorSpace;

        public Transform(VectorSpace<V, F> vectorSpace, LinearTransform transform)
            : base(vectorSpace, vectorSpace, transform) { this.vectorSpace = vectorSpace; }

        public override Vector<F> Add(Vector<F> vector)
        {
            Transform<V, F> other = vector as Transform<V, F>;
            return new Transform<V, F>(vectorSpace, vec => (V)(transform(vec) + other.transform(vec)));
        }

        public override Vector<F> Multiply(F scalar)
        {
            return new Transform<V, F>(vectorSpace, vector => (V)(scalar * vector));
        }

        public override Vector<F> Zero()
        {
            return new Transform<V, F>(vectorSpace, zero);
        }

        public override Vector<F> AddInverse()
        {
            return new Transform<V, F>(vectorSpace, vector => (V)vector.AddInverse());
        }

        public Transform<V, F> Identity()
        {
            return new Transform<V, F>(vectorSpace, identity);
        }

        public static Transform<V, F> operator * (Transform<V, F> T, Transform<V, F> S)
        {
            return new Transform<V, F>(T.vectorSpace, vector => T.transform(S.transform(vector)));
        }

        public Transform<V, F> Power(int exp)
        {
            Transform<V, F> res = identity as Transform<V, F>;
            Transform<V, F> pow = this;
            while (exp > 0)
            {
                if (exp % 2 == 1)
                {
                    res *= pow;
                    exp--;
                }
                pow *= pow;
                exp /= 2;
            }
            return res;
        }
    }
}
