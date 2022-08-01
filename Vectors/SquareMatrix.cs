namespace Linear_Algebra
{
    class SquareMatrix<F> : Matrix<F> where F : Field
    {
        private readonly int size;

        private SquareMatrix(int n) : base(n, n) 
        {
            size = n;
        }

        // @pre no null items in entries
        public SquareMatrix(F[,] entries) : base(entries)
        {
            int rows = entries.GetLength(0), cols = entries.GetLength(1);
            if(rows != cols) { throw new MatrixSizeException(string.Format(
                "Invalid dimensions for square matrix: {0}, {1}", rows, cols)); }

            size = entries.GetLength(0);
        }

        public static SquareMatrix<F> operator * (SquareMatrix<F> mat1, SquareMatrix<F> mat2)
        {
            return new SquareMatrix<F>(mat1.MatMultiplication(mat2));
        }

        public override SquareMatrix<F> Clone()
        {
            return new SquareMatrix<F>(entries);
        }

        public SquareMatrix<F> Identity()
        {
            SquareMatrix<F> id = new SquareMatrix<F>(size);
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    id[i, j] = i == j ? FieldOne() : FieldZero();
                }
            }
            return id;
        }

        public F Trace()
        {
            F sum = FieldZero();

            for (int i = 0; i < size; i++)
            {
                sum = (F)(sum + this[i, i]);
            }
            return sum;
        }

        public Transform<ColumnVector<F>, F> ToTransform()
        {
            return new Transform<ColumnVector<F>, F>(ColumnSpace(), vector => this * vector);
        }

        // @pre no null entries
        public static SquareMatrix<F> Diag(params F[] entries)
        {
            SquareMatrix<F> mat = new SquareMatrix<F>(entries.Length);
            F zero = (F)entries[0].Zero();
            for (int i = 0; i < mat.size; i++)
            {
                for (int j = 0; j < mat.size; j++)
                {
                    mat[i, j] = i == j ? entries[i] : zero;
                }
            }
            return mat;
        }

        // @pre exp >= 0
        public SquareMatrix<F> Power(int exp)
        {
            SquareMatrix<F> res = Identity();
            SquareMatrix<F> pow = this;
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

        // @post $this will be upper triangular
        // @post $ret will be lower triangular
        // no swaps were performed $implies $prev == $ret * $this
        public SquareMatrix<F> LUDecomposition()
        {
            SquareMatrix<F> L = Identity();
            AddOp add = (src, dst, scalar) => L[dst, src] = (F)scalar.AddInverse();
            RowEchelonForm(this, AddRow + add, SwapRows);
            return L;
        }

        #region Determinants, Adjoints and Inverses

        public F Determinant()
        {
            F[,] M = entries.Clone() as F[,];
            Field inverse;
            int n = size;
            for (int k = 0; k < n - 1; k++)
            {
                inverse = k == 0 ? FieldOne() : M[k - 1, k - 1].MultInverse();
                for (int i = k + 1; i < n; i++)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        M[i, j] = (F)((M[i, j] * M[k, k] - M[i, k] * M[k, j]) * inverse);
                    }
                }
            }
            return M[n - 1, n - 1];
        }

        private F Minor(int row, int col)
        {
            SquareMatrix<F> minor = new SquareMatrix<F>(size - 1);
            int r = 0, c = 0;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if(i != row && j != col)
                    {
                        minor[r, c++] = this[i, j];
                        if(c == size - 1)
                        {
                            c = 0;
                            r++;
                        }
                    }
                }
            }
            return minor.Determinant();
        }

        public SquareMatrix<F> Adjoint()
        {
            SquareMatrix<F> adj = Clone();
            F sign;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    sign = (i + j) % 2 == 0 ? FieldOne() : (F)FieldOne().AddInverse();
                    adj[i, j] = (F)(sign * Minor(j, i));
                }
            }
            return adj;
        }

        public bool IsInvertible()
        {
            return !Determinant().Equals(FieldZero());
        }

        public SquareMatrix<F> Inverse()
        {
            SquareMatrix<F> clone = Clone(), id = Identity();
            AddOp add = clone.AddRow;
            MultiplyOp multiply = clone.MultiplyRow;
            SwapOp swap = clone.SwapRows;

            int rank = GaussianElimination(clone, 
                add + id.AddRow , multiply + id.MultiplyRow, swap + id.SwapRows);
            return rank == size ? id : null;
        }

        #endregion
    }
}
