namespace Linear_Algebra
{
    class Matrix<F> : AbstractMatrix<F> where F : Field
    {

        protected Matrix(int m, int n) : base(m, n) { }

        public Matrix(F[,] entries) : base(entries) { }

        public override Matrix<F> Clone()
        {
            return new Matrix<F>(entries);
        }

        public override ColumnVector<F> ToColumnVector()
        {
            F[] vec = new F[entries.Length];
            int k = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    vec[k++] = this[i, j];
                }
            }
            return new ColumnVector<F>(vec);
        }

        public Matrix<F> Transpose()
        {
            Matrix<F> tr = new Matrix<F>(cols, rows);
            for (int i = 0; i < cols; i++)
            {
                for (int j = 0; j < rows; j++)
                {
                    tr[i, j] = this[j, i];
                }
            }
            return tr;
        }

        public Transform<ColumnVector<F>, F> ToTransform()
        {
            return new Transform<ColumnVector<F>, F>(vec => this * vec);
        }

        #region Gaussian Elimination and matrix ranking

        protected static int RowEchelonForm(Matrix<F> matrix, AddOp add, SwapOp swap)
        {
            int pos, rank = 0;
            Field inverse;
            for (int i = 0; i < matrix.cols; i++)
            {
                // Find first non-zero entry in the column
                // If not found, continue to the next one
                pos = rank;
                while(pos < matrix.rows && matrix[pos, i].Equals(matrix.FieldZero())) { pos++; }
                if (pos == matrix.rows) { continue; }

                // if first non-zero entry is after the current row, swap them
                if (pos > rank)
                {
                    swap(pos, rank);
                }

                // Use the add operation to nullify every entry in the column below the current row 
                inverse = matrix[rank, i].MultInverse();
                for (int j = pos + 1; j < matrix.rows; j++)
                {
                    add(rank, j, (F)(matrix[j, i].AddInverse() * inverse));
                }
                rank++;
            }
            return rank;
        }

        public int Rank()
        {
            Matrix<F> clone = Clone();
            return RowEchelonForm(clone, clone.AddRow, clone.SwapRows);
        }

        public Matrix<F> EchelonForm()
        {
            Matrix<F> clone = Clone();
            RowEchelonForm(clone, clone.AddRow, clone.SwapRows);
            return clone;
        }

        protected static int GaussianElimination(Matrix<F> matrix, AddOp add, MultiplyOp multiply, SwapOp swap)
        {
            // First bring the matrix to a row echelon form
            int col = 0, rank = RowEchelonForm(matrix, add, swap);
            for (int row = 0; row < rank; row++)
            {
                // Find the leading entry in a row and multiply it by its inverse
                while(matrix[row, col].Equals(matrix.FieldZero())) { col++; }
                multiply(row, (F)matrix[row, col].MultInverse());

                // nullify each entry in the column above the row
                for (int i = 0; i < row; i++)
                {
                    add(row, i, (F)matrix[i, col].AddInverse());
                }
            }
            return rank;
        }

        public Matrix<F> CanonForm()
        {
            Matrix<F> clone = Clone();
            GaussianElimination(clone, clone.AddRow, clone.MultiplyRow, clone.SwapRows);
            return clone;
        }

        public VectorSpace<ColumnVector<F>, F> NullSpace()
        {
            // Use gaussian elimination to bring the matrix into a simpler form
            Matrix<F> matrix = Clone();
            int rank = GaussianElimination(matrix, matrix.AddRow, matrix.MultiplyRow, matrix.SwapRows);
            VectorSpace<ColumnVector<F>, F> nullSpace = new VectorSpace<ColumnVector<F>, F>(matrix.cols);
            F[] vector;

            // Keep the positions of the leading entries in each row
            int pos = 0;
            int[] leadingEntries = new int[rank];
            for (int i = 0; i < matrix.cols && pos < rank; i++)
            {
                if(!matrix[pos, i].Equals(matrix.FieldZero())) { leadingEntries[pos++] = i; }
            }

            for (int col = 0; col < matrix.cols; col++)
            {
                // Check whether the column has a non-zero entry
                // which isn't a leading one, otherwise continue
                pos = 0;
                while (pos < rank && matrix[pos, col].Equals(matrix.FieldZero())) { pos++; }
                if (pos == rank || leadingEntries[pos] == col) { continue; }

                // Initialize a vector with a length equal to the number of columns
                vector = new F[matrix.cols];
                for (int i = 0; i < matrix.cols; i++) { vector[i] = matrix.FieldZero(); }

                // For every row with a non-zero (a.k.a the first #rank rows) place an entry in the vector
                // also place a one in the entry according to the current column
                for (int i = 0; i < rank; i++)
                {
                    vector[leadingEntries[i]] = (F)matrix[i, col].AddInverse();
                }
                vector[col] = matrix.FieldOne();

                // Finally the vector is ready to be added to the null space
                nullSpace.Add(new ColumnVector<F>(vector));
            }
            return nullSpace;
        }

        #endregion

        public VectorSpace<ColumnVector<F>, F> ColumnSpace()
        {
            VectorSpace<ColumnVector<F>, F> colSpace = new VectorSpace<ColumnVector<F>, F>(rows);
            F[] vector;
            for (int i = 0; i < cols; i++)
            {
                vector = new F[rows];
                for (int j = 0; j < rows; j++)
                {
                    vector[j] = this[j, i];
                }
                colSpace.Add(new ColumnVector<F>(vector));
            }
            return colSpace;
        }

        public static ColumnVector<F> LinearSystemSolution(Matrix<F> matrix, ColumnVector<F> vector)
        {
            if(matrix.rows != vector.length) { throw new MatrixSizeException(string.Format(
                "Invalid configuration for a system of linear equations:\n{0}\n{1}", matrix, vector)); }

            matrix = matrix.Clone();
            vector = vector.Clone();
            F[] sol = new F[matrix.cols];
            for (int i = 0; i < matrix.cols; i++) { sol[i] = matrix.FieldZero(); }

            // Use Gaussian Elimination
            AddOp add = (src, dst, scalar) => vector[dst] = (F)(vector[dst] + scalar * vector[src]);
            MultiplyOp multiply = (row, scalar) => vector[row] = (F)(scalar * vector[row]);
            SwapOp swap = (row1, row2) => {
                F tmp = vector[row1];
                vector[row1] = vector[row2];
                vector[row2] = tmp;
            };
            int rank = GaussianElimination(matrix,
                matrix.AddRow + add, matrix.MultiplyRow + multiply, matrix.SwapRows + swap);

            // Check each entry according to a row of zeros to see whether a solution exists at all
            for (int i = rank; i < matrix.rows; i++)
            {
                if(!vector[i].Equals(matrix.FieldZero()))
                {
                    return null;
                }
            }

            // copy entries to the solution vector
            int pos = 0;
            for (int i = 0; i < matrix.cols && pos < matrix.rows; i++)
            {
                if (matrix[pos, i].Equals(matrix.FieldOne())) 
                { 
                    sol[i] = vector[pos++];
                }
            }
            return new ColumnVector<F>(sol);
        }
    }
}
