using System.Collections.Generic;


namespace Linear_Algebra
{
    class Transform<V, W, F> : Vector<F> where F : Field where V : Vector<F> where W : Vector<F>
    {
        public delegate W LinearTransform(V vector);

        protected readonly LinearTransform transform;
        public readonly VectorSpace<V, F> domain;
        protected readonly VectorSpace<W, F> range;

        // @pre for all vector v in domain, range.Contains(transform(v))
        // @pre transform is linear; for all vectors u, v and scalar a: transform(a*v + u) == a*transform(v) + transform(u)
        public Transform(VectorSpace<V, F> domain, VectorSpace<W, F> range, LinearTransform transform)
        {
            this.domain = domain;
            this.range = range;
            this.transform = transform;
        }

        public ColumnVector<F> ToColumnVector()
        {
            return MatrixRepresentation().ToColumnVector();
        }

        public virtual Transform<V, W, F> Clone(LinearTransform t)
        {
            return new Transform<V, W, F>(domain, range, t);
        }

        public virtual Transform<V, W, F> Clone()
        {
            return Clone(transform);
        }

        // @pre domain.Contains(vector)
        public W ValueOf(V vector)
        {
            return transform(vector);
        }

        // @pre transform.domain.Contains(vector)
        public static W operator * (Transform<V, W, F> transform, V vector)
        {
            return transform.ValueOf(vector);
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

        public virtual Matrix<F> MatrixRepresentation()
        {
            return MatrixRepresentation(domain, range);
        }

        // @pre B == domain && C == range
        public virtual Matrix<F> MatrixRepresentation(VectorSpace<V, F> B, VectorSpace<W, F> C)
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
            VectorSpace<ColumnVector<F>, F> nullSpace = MatrixRepresentation().NullSpace();
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
            return MatrixRepresentation().ToString();
        }

        public int Length()
        {
            return domain.dim * range.dim;
        }
    }

    class Transform<V, F> : Transform<V, V ,F> where F : Field where V : Vector<F>
    {
        public static readonly LinearTransform zero = vector => (V)vector.Zero();
        public static readonly LinearTransform identity = vector => vector;

        // @pre vectorSpace is transform-invariant
        // meaning that for all vector v such that vectorSpace.Contains(v), vectorSpace.Contains(transform(v))
        public Transform(VectorSpace<V, F> vectorSpace, LinearTransform transform)
            : base(vectorSpace, vectorSpace, transform) { }

        public override Transform<V, F> Clone(LinearTransform t)
        {
            return new Transform<V, F>(domain, t);
        }

        public override Transform<V, F> Clone()
        {
            return Clone(transform);
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

        public override SquareMatrix<F> MatrixRepresentation()
        {
            return MatrixRepresentation(domain);
        }

        // @pre B == domain
        public SquareMatrix<F> MatrixRepresentation(VectorSpace<V, F> B)
        {
            return new SquareMatrix<F>(MatRep(B, B));
        }

        // @pre B == C == domain
        public override SquareMatrix<F> MatrixRepresentation(VectorSpace<V, F> B, VectorSpace<V, F> C)
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

        public Transform<QuotientSpace<V, F>, F> QuotientTransform(VectorSpace<V, F> subSpace)
        {
            VectorSpace<QuotientSpace<V, F>, F> quotientSpace = new VectorSpace<QuotientSpace<V, F>, F>(domain.dim);
            foreach(V vector in domain)
            {
                quotientSpace.Add(vector + subSpace);
            }
            return new Transform<QuotientSpace<V, F>, F>(quotientSpace,
                q => new QuotientSpace<V, F>(transform(q.vector), subSpace));
        }

        // @pre this is a nilpotent transformation
        // meaning for some k, this.Power(k) == Zero()
        public Transform<V, F> NilpotentJordanForm()
        {
            Transform<V, F> T = Identity();
            VectorSpace<V, F> B = domain.Clone(), C, ker = domain.Clone();
            List<Transform<V, F>> powers = new List<Transform<V, F>>();
            List<VectorSpace<V, F>> kernels = new List<VectorSpace<V, F>>();

            // Make a list of powers of the transform, as well as the kernels of these powers
            powers.Add(Identity());
            kernels.Add(ker);
            int l = 0;
            while (!ker.IsSpanning())
            {
                T *= this;
                ker = T.Kernel();
                powers.Add(T);
                kernels.Add(ker);
                l++;
            }

            for (int i = l; i >= 1; i--)
            {
                // Step 1: Take C to be the intersection of B and ker and complete it to ker
                ker = kernels[i - 1];
                C = VectorSpace<V, F>.Intersection(B, ker);
                foreach(V vector in ker) { C.Add(vector); }

                // Step 2: Add the remaining vectors of the next kernel intersection
                ker = kernels[i];
                foreach(V vector in VectorSpace<V, F>.Intersection(B, ker)) { C.Add(vector); }

                // Step 3: each vector in ker which isn't in C starts a chain
                l = ker.Dimension() - C.Dimension();
                foreach(V vector in ker)
                {
                    if (l == 0) { break; }
                    if (C.Contains(vector)) { continue; }
                    for (int j = i-1; j >= 0; j--)
                    {
                        B.Add(powers[j] * vector);
                    }
                    l--;
                }
            }
            return new Transform<V, F>(B, transform);
        }
    }
}
