using System;
using System.Collections.Generic;


namespace Linear_Algebra
{
    class Program
    {
        static void Main(string[] args)
        {
            try
            {
#if false       // Matrix addition and multiplication
                Matrix<Real> A = new Matrix<Real>(new Real[,]
                {
                    {new Real(1), new Real(0), new Real(-3) },
                    {new Real(-2), new Real(4), new Real(1) }
                });
                Matrix<Real> B = new Matrix<Real>(new Real[,]
                {
                    {new Real(3), new Real(-2), new Real(0) },
                    {new Real(1), new Real(-1), new Real(-5) }
                });
                Console.WriteLine("A =\n" + A);
                Console.WriteLine("B =\n" + B);
                Console.WriteLine("A + B =");
                Console.WriteLine((Vector<Real>)A + B);

                Matrix<Real> B_T = B.Transpose();
                Console.WriteLine("B^T =");
                Console.WriteLine(B_T);
                Console.WriteLine("A * B^T =");
                Console.WriteLine(A * B_T);
                Console.WriteLine("B^T * A =");
                Console.WriteLine(B_T * A);
#endif

#if false       // Determinant and inverse matrix
                SquareMatrix<Complex> mat = new SquareMatrix<Complex>(new Complex[,]
                {
                    {new Real(7), Real.zero, new Complex(1, 1) },
                    {Real.zero, Real.one, new Complex(0, 9) },
                    {new Complex(1, -1), new Complex(0, -4), new Real(-10) }
                });
                Console.WriteLine("mat =\n" + mat);
                Console.WriteLine("Det(mat) = " + mat.Determinant() + "\n");
                SquareMatrix<Complex> inv = mat.Inverse();
                Console.WriteLine("mat^-1 =\n" + inv);
                Console.WriteLine("mat * inv =");
                Console.WriteLine(mat * inv);
#endif

#if false       // Matrix ranking algorithm
                Matrix<Real> A = new Matrix<Real>(new Real[,]
                {
                    {new Real(1), new Real(2), new Real(9) },
                    {new Real(2), new Real(1), new Real(12) },
                    {new Real(3), new Real(0), new Real(15) }
                });
                Console.WriteLine("A = ");
                Console.WriteLine(A);
                Console.WriteLine("Row echelon form:");
                Console.WriteLine(A.EchelonForm());
                Console.WriteLine("Canon form:");
                Console.WriteLine(A.CanonForm());
                Console.WriteLine("rank(A) = " + A.Rank());
#endif

#if false       // Linear equations
                Matrix<Real> A = new Matrix<Real>(new Real[,]
                {
                    {new Real(1), new Real(-2) },
                    {new Real(3), new Real(5) },
                    {new Real(4), new Real(3) }
                });
                ColumnVector<Real> b = new ColumnVector<Real>(new Real(-1), new Real(8), new Real(7));
                Console.WriteLine("A = \n" + A);
                Console.WriteLine("b = \n" + b);
                Console.WriteLine("Solution for Ax = b:");
                Console.WriteLine(Matrix<Real>.LinearSystemSolution(A, b));
#endif

#if false       // Nullspace
                Matrix<Real> C = new Matrix<Real>(new Real[,]
                {
                    {new Real(3), new Real(6), new Real(6), new Real(3), new Real(9) },
                    {new Real(6), new Real(12), new Real(13), new Real(0), new Real(3) }
                });
                Console.WriteLine("C = \n" + C);
                Console.WriteLine("Null(C) = ");
                Console.WriteLine(C.NullSpace());
#endif

#if false       // Row and column space
                Matrix<Real> A = new Matrix<Real>(new Real[,]
                {
                    {new Real(1), new Real(3), new Real(1), new Real(4) },
                    {new Real(2), new Real(7), new Real(3), new Real(9) },
                    {new Real(1), new Real(5), new Real(3), new Real(1) },
                    {new Real(1), new Real(2), new Real(0), new Real(8) }
                });
                Console.WriteLine("A = \n" + A);
                Console.WriteLine("Column space");
                VectorSpace<ColumnVector<Real>, Real> colspace = A.ColumnSpace();
                Console.WriteLine(colspace);
                Console.WriteLine("Row space:");
                Console.WriteLine(A.RowSpace());

                ColumnVector<Real> col = new ColumnVector<Real>(new Real(1), new Real(3), new Real(3), new Real(0));
                Console.WriteLine("Coordinates of third column:");
                ColumnVector<Real> coords = colspace.CoordinatesOf(col);
                Console.WriteLine(coords);
                Console.WriteLine("Third column:");
                Console.WriteLine(colspace.LinearCombination(coords));
#endif

#if false       // Kernel and image of a linear transformation
                VectorSpace<ColumnVector<Real>, Real> R2 = new VectorSpace<ColumnVector<Real>, Real>(2,
                    new ColumnVector<Real>(Real.one, Real.zero),
                    new ColumnVector<Real>(Real.zero, Real.one)
                    );
                VectorSpace<Polynomial<Real>, Real> Pol3 = new VectorSpace<Polynomial<Real>, Real>(4,
                    new Polynomial<Real>(3, Real.one),
                    new Polynomial<Real>(3, Real.zero, Real.one),
                    new Polynomial<Real>(3, Real.zero, Real.zero, Real.one),
                    new Polynomial<Real>(3, Real.zero, Real.zero, Real.zero, Real.one)
                    );
                Transform<Polynomial<Real>, ColumnVector<Real>, Real> T = new Transform<Polynomial<Real>, ColumnVector<Real>, Real>(
                    Pol3, R2,
                    p => new ColumnVector<Real>(p.ValueOf(Real.zero), p.ValueOf(Real.one))
                    );
                Console.WriteLine("T(p) = [p(0), p(1)]\n");
                Console.WriteLine("Kernel:\n" + T.Kernel() + "\n");
                Console.WriteLine("Image:\n" + T.Image());
#endif
            }
            catch(MatrixSizeException ex)
            {
                Console.WriteLine("Matrix size exception.");
                Console.WriteLine(ex.Message);
            }
        }
    }
}
