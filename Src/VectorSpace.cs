using System.Collections;
using System.Collections.Generic;
using System.Linq;


namespace Linear_Algebra
{
    // invariant: basis is a list of linearly independent vectors
    class VectorSpace<V, F> : IEnumerable<V> where F : Field where V : Vector<F>
    {
        protected List<V> basis;
        protected SquareMatrix<F> matrixRep;
        public readonly int dim;

        public VectorSpace(int dim)
        {
            this.dim = dim;
            basis = new List<V>();
        }

        // @pre no null vectors
        // @pre all vectors added to the space must have the size as dim
        public VectorSpace(int dim, params V[] vectors)
            : this(dim) { Add(vectors); }

        public virtual VectorSpace<V, F> Clone()
        {
            return new VectorSpace<V, F>(dim);
        }

        public int Dimension() { return basis.Count; }

        public bool IsEmpty() { return !basis.Any(); }

        public bool IsSpanning() { return Dimension() == dim; }

        // @pre no null vectors
        // @pre all vectors added to the space must have their size be the same set as the space dimension
        // for example "new VectorSpace<ColumnVector<Real>, Real>(3)" can only have column vectors with 3 entries
        public void Add(params V[] vectors)
        {
            foreach(V vec in vectors)
            {
                if(vec.Length() != dim) { continue; }
                if (Dimension() < dim && !Contains(vec))
                {
                    basis.Add(vec);
                    GetMatrixRep();
                }
            }
        }

        // @pre vector != null
        public bool Contains(V vector)
        {
            if (vector.Equals(vector.Zero())) { return true; }
            if (IsEmpty()) { return false; }
            return Matrix<F>.LinearSystemSolution(matrixRep, vector.ToColumnVector()) != null;
        }

        // @pre !IsEmpty() && vector.Length() == dim
        public ColumnVector<F> CoordinatesOf(V vector)
        {
            return Matrix<F>.LinearSystemSolution(matrixRep, vector.ToColumnVector());
        }

        private void GetMatrixRep()
        {
            if (IsEmpty()) { return; }
            F[,] matrix = new F[dim, dim];
            ColumnVector<F> colVec;
            int col = 0;
            F zero = (F)basis[0].ToColumnVector()[0].Zero();
            foreach (V vec in this)
            {
                colVec = vec.ToColumnVector();
                for (int row = 0; row < dim; row++)
                {
                    matrix[row, col] = colVec[row];
                }
                col++;
            }
            while(col < dim)
            {
                for (int row = 0; row < dim; row++)
                {
                    matrix[row, col] = zero;
                }
                col++;
            }
            matrixRep = new SquareMatrix<F>(matrix);
        }

        // @pre from.dim == to.dim && from.IsSpanning() && to.IsSpanning()
        public static SquareMatrix<F> TransitionMatrix(VectorSpace<V, F> from, VectorSpace<V, F> to)
        {
            return to.matrixRep.Inverse() * from.matrixRep;
        }

        // @pre scalars.length == Dimension() && !IsEmpty()
        public V LinearCombination(ColumnVector<F> scalars)
        {
            V linearComb = (V)basis[0].Zero();
            for (int i = 0; i < Dimension(); i++)
            {
                linearComb = (V)(linearComb + scalars[i] * basis[i]);
            }
            return linearComb;
        }

        // @pre U.dim == W.dim
        public static VectorSpace<V, F> operator + (VectorSpace<V, F> U, VectorSpace<V, F> W)
        {
            VectorSpace<V, F> V = U.Clone();
            V.Add(U.basis.ToArray());
            V.Add(W.basis.ToArray());
            return V;
        }

        // @pre U.dim == W.dim
        public static VectorSpace<V, F> Intersection(VectorSpace<V, F> U, VectorSpace<V, F> W)
        {
            F[,] matrix = new F[U.dim, U.Dimension() + W.Dimension()];
            ColumnVector<F> colVec;
            int col = 0;
            foreach (V vec in U)
            {
                colVec = vec.ToColumnVector();
                for (int row = 0; row < U.dim; row++)
                {
                    matrix[row, col] = colVec[row];
                }
                col++;
            }
            foreach (V vec in W)
            {
                colVec = vec.ToColumnVector();
                for (int row = 0; row < W.dim; row++)
                {
                    matrix[row, col] = (F)colVec[row].AddInverse();
                }
                col++;
            }

            VectorSpace<ColumnVector<F>, F> nullSpace = new Matrix<F>(matrix).NullSpace();
            VectorSpace<V, F> V = U.Clone();
            F[] scalars;
            foreach(ColumnVector<F> vec in nullSpace)
            {
                scalars = new F[U.Dimension()];
                for (int i = 0; i < U.Dimension(); i++)
                {
                    scalars[i] = vec[i];
                }
                V.Add(U.LinearCombination(new ColumnVector<F>(scalars)));
            }
            return V;
        }

        public static QuotientSpace<V, F> operator + (V vector, VectorSpace<V, F> subSpace)
        {
            return new QuotientSpace<V, F>(vector, subSpace);
        }

        public void Reverse() { basis.Reverse(); }

        public override string ToString() 
        {
            if (IsEmpty()) { return "{ }"; }
            return string.Join("\n", basis); 
        }

        public override bool Equals(object obj)
        {
            VectorSpace<V, F> U = obj as VectorSpace<V, F>;
            if(dim != U.dim || Dimension() != U.Dimension()) { return false; }
            foreach(V vec in U)
            {
                if (!Contains(vec)) { return false; }
            }
            return true;
        }

        public IEnumerator<V> GetEnumerator() { return basis.GetEnumerator(); }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }
}
