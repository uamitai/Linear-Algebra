namespace Linear_Algebra
{
    class Transform<V, W, F> : Vector<F> where F : Field where V : Vector<F> where W : Vector<F>
    {
        public delegate W LinearTransform(V vector);

        protected readonly LinearTransform transform;
        protected readonly VectorSpace<V, F> domain;
        protected readonly VectorSpace<W, F> range;
        public readonly Matrix<F> matrixRep;

        // @pre for all vector v in domain, range.Contains(transform(v))
        public Transform(VectorSpace<V, F> domain, VectorSpace<W, F> range, LinearTransform transform)
        {
            this.domain = domain;
            this.range = range;
            this.transform = transform;
            this.matrixRep = GetMatrixRep();
        }

        public ColumnVector<F> ToColumnVector()
        {
            return matrixRep.ToColumnVector();
        }

        public virtual Transform<V, W, F> Clone(LinearTransform t)
        {
            return new Transform<V, W, F>(domain, range, t);
        }

        public Transform<V, W, F> Clone()
        {
            return Clone(transform);
        }

        // @pre domain.Contains(vector)
        public W ValueOf(V vector)
        {
            return transform(vector);
        }

        public virtual Vector<F> Add(Vector<F> vector)
        {
            Transform<V, W, F> other = vector as Transform<V, W, F>;
            return Clone(vec => (W)(transform(vec) + other.transform(vec)));
        }

        public virtual Vector<F> Multiply(F scalar)
        {
            return Clone(vector => (W)(scalar * transform(vector)));
        }

        public virtual Vector<F> Zero()
        {
            return Clone(vector => (W)transform(vector).Zero());
        }

        public virtual Vector<F> AddInverse()
        {
            return Clone(vector => (W)transform(vector).AddInverse());
        }

        public VectorSpace<W, F> Image()
        {
            return Image(domain);
        }

        // @pre for all vector v in vectorSpace, domain.Contains(v)
        public VectorSpace<W, F> Image(VectorSpace<V, F> vectorSpace)
        {
            VectorSpace<W, F> im = range.Clone();
            foreach (V vector in vectorSpace)
            {
                im.Add(transform(vector));
            }
            return im;
        }

        public virtual Matrix<F> GetMatrixRep()
        {
            return new Matrix<F>(MatRep(domain, range));
        }

        // @pre B == domain && C == range
        public virtual Matrix<F> GetMatrixRep(VectorSpace<V, F> B, VectorSpace<W, F> C)
        {
            return new Matrix<F>(MatRep(B, C));
        }

        protected F[,] MatRep(VectorSpace<V, F> domain, VectorSpace<W, F> range)
        {
            F[,] rep = new F[range.Dimension(), domain.Dimension()];
            int col = 0;
            ColumnVector<F> colVec;
            foreach (V vector in domain)
            {
                colVec = range.CoordinatesOf(transform(vector));
                for (int row = 0; row < range.Dimension(); row++)
                {
                    rep[row, col] = colVec[row];
                }
                col++;
            }
            return rep;
        }

        public VectorSpace<V, F> Kernel()
        {
            VectorSpace<V, F> kernel = domain.Clone();
            VectorSpace<ColumnVector<F>, F> nullSpace = matrixRep.NullSpace();
            foreach (ColumnVector<F> vector in nullSpace)
            {
                kernel.Add(domain.LinearCombination(vector));
            }
            return kernel;
        }

        // @pre domain.Dimension() == range.Dimension()
        public static Transform<V, W, F> FromBases(VectorSpace<V, F> domain, VectorSpace<W, F> range)
        {
            return new Transform<V, W, F>(domain, range, vector => range.LinearCombination(domain.CoordinatesOf(vector)));
        }

        public bool IsInvertible()
        {
            return Kernel().IsEmpty();
        }

        public bool IsIsomorphism()
        {
            return IsInvertible() && Image().Equals(range);
        }

        // @pre IsInvertible()
        public virtual Transform<W, V, F> Inverse()
        {
            return Transform<W, V, F>.FromBases(Image(), domain);
        }

        public override bool Equals(object obj)
        {
            Transform<V, W, F> T = obj as Transform<V, W, F>;
            if (!domain.Equals(T.domain)) { return false; }
            foreach(V vector in domain)
            {
                if (!transform(vector).Equals(T.transform(vector))){ return false; }
            }
            return true;
        }

        public override string ToString()
        {
            return matrixRep.ToString();
        }
    }

    class Transform<V, F> : Transform<V, V ,F> where F : Field where V : Vector<F>
    {
        public static readonly LinearTransform zero = vector => (V)vector.Zero();
        public static readonly LinearTransform identity = vector => vector;

        public new readonly SquareMatrix<F> matrixRep;

        // @pre vectorSpace is transform-invariant
        // meaning that for all vector v such that vectorSpace.Contains(v), vectorSpace.Contains(transform(v))
        public Transform(VectorSpace<V, F> vectorSpace, LinearTransform transform)
            : base(vectorSpace, vectorSpace, transform) { matrixRep = GetMatrixRep(); }

        public override Transform<V, V, F> Clone(LinearTransform t)
        {
            return new Transform<V, F>(domain, t);
        }

        public override Vector<F> Zero()
        {
            return new Transform<V, F>(domain, zero);
        }

        public Transform<V, F> Identity()
        {
            return new Transform<V, F>(domain, identity);
        }

        public static Transform<V, F> operator * (Transform<V, F> T, Transform<V, F> S)
        {
            return new Transform<V, F>(T.domain, vector => T.transform(S.transform(vector)));
        }

        // @pre IsInvertible()
        public override Transform<V, F> Inverse()
        {
            VectorSpace<V, F> im = Image();
            return new Transform<V, F>(im, vector => domain.LinearCombination(im.CoordinatesOf(vector)));
        }

        public override SquareMatrix<F> GetMatrixRep()
        {
            return new SquareMatrix<F>(MatRep(domain, range));
        }

        // @pre B == domain
        public SquareMatrix<F> GetMatrixRep(VectorSpace<V, F> B)
        {
            return new SquareMatrix<F>(MatRep(B, B));
        }

        // @pre B == C == domain
        public override SquareMatrix<F> GetMatrixRep(VectorSpace<V, F> B, VectorSpace<V, F> C)
        {
            return new SquareMatrix<F>(MatRep(B, C));
        }

        public Transform<V, F> Power(int exp)
        {
            Transform<V, F> res = Identity();
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

        public VectorSpace<V, F> EigenSpace(F eigenValue)
        {
            return ((this - Identity().Multiply(eigenValue)) as Transform<V, F>).Kernel();
        }
    }
}
