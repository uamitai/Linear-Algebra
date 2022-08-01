using System;


namespace Linear_Algebra
{
    class MatrixSizeException : Exception
    {
        public MatrixSizeException() : base() { }
        public MatrixSizeException(string message) : base(message) { }
        public MatrixSizeException(string message, Exception inner) : base(message, inner) { }
    }

    abstract class AbstractMatrix<F> : Vector<F> where F : Field
    {
        protected readonly F[,] entries;
        protected readonly int rows, cols;

        protected AbstractMatrix(int m, int n)
        {
            if (m <= 0 || n <= 0) { throw new MatrixSizeException(string.Format(
                "Cannot create matrix with non-positive dimensions: {0}, {1}", m, n)); }

            rows = m;
            cols = n;
            entries = new F[m, n];
        }

        // @pre no null items in entries
        public AbstractMatrix(F[,] entries) : this(entries.GetLength(0), entries.GetLength(1))
        {
            this.entries = entries.Clone() as F[,];
        }

        abstract public AbstractMatrix<F> Clone();
        abstract public ColumnVector<F> ToColumnVector();
        public int Length() { return entries.Length; }

        public F this[int i, int j]
        {
            get { return entries[i, j]; }
            protected set { entries[i, j] = value; }
        }

        public Vector<F> Add(Vector<F> vec)
        {
            AbstractMatrix<F> matrix = vec as AbstractMatrix<F>;
            if (rows != matrix.rows || cols != matrix.cols) { throw new MatrixSizeException(string.Format(
                "Invalid sizes for matrix addition:\n{0}\n{1}", this, matrix)); }

            AbstractMatrix<F> res = Clone();
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    res[i, j] = (F)(this[i, j] + matrix[i, j]);
                }
            }
            return res;
        }

        public Vector<F> Multiply(F scalar)
        {
            AbstractMatrix<F> res = Clone();
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    res[i, j] = (F)(scalar * this[i, j]);
                }
            }
            return res;
        }

        public virtual Vector<F> Zero()
        {
            AbstractMatrix<F> zero = Clone();
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    zero[i, j] = FieldZero();
                }
            }
            return zero;
        }

        public Vector<F> AddInverse()
        {
            AbstractMatrix<F> res = Clone();
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    res[i, j] = (F)this[i, j].AddInverse();
                }
            }
            return res;
        }

        protected F[,] MatMultiplication(AbstractMatrix<F> other)
        {
            if(cols!= other.rows) { throw new MatrixSizeException(string.Format(
                "Invalid sizes for matrix multiplication:\n{0}\n{1}", this, other)); }

            F[,] res = new F[rows, other.cols];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < other.cols; j++)
                {
                    F sum = (F)(this[i, 0] * other[0, j]);
                    for (int k = 1; k < other.rows; k++)
                    {
                        sum = (F)(sum + this[i, k] * other[k, j]);
                    }
                    res[i, j] = sum;
                }
            }
            return res;
        }

        public static AbstractMatrix<F> operator * (AbstractMatrix<F> mat1, AbstractMatrix<F> mat2)
        {
            return new Matrix<F>(mat1.MatMultiplication(mat2));
        }

        public static ColumnVector<F> operator * (AbstractMatrix<F> matrix, ColumnVector<F> vector)
        {
            return new ColumnVector<F>(matrix.MatMultiplication(vector));
        }

        protected F FieldZero() { return (F)this[0, 0].Zero(); }

        protected F FieldOne() { return (F)this[0, 0].One(); }

        public override bool Equals(object obj)
        {
            AbstractMatrix<F> other = obj as AbstractMatrix<F>;
            if (rows != other.rows || cols != other.cols) { return false; }
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (!this[i, j].Equals(other[i, j])) { return false; }
                }
            }
            return true;
        }

        public override string ToString()
        {
            string res = "";
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    res += this[i, j] + " ";
                }
                res += "\n";
            }
            return res;
        }
    }
}
