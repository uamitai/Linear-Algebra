using System;


namespace Linear_Algebra
{
    class Program
    {
        static void Main(string[] args)
        {
            try
            {
                VectorSpace<ColumnVector<Real>, Real> R2 = new VectorSpace<ColumnVector<Real>, Real>(2,
                    new ColumnVector<Real>(Real.one, Real.zero),
                    new ColumnVector<Real>(Real.zero, Real.one)
                    );
                VectorSpace<SquareMatrix<Real>, Real> M22 = new VectorSpace<SquareMatrix<Real>, Real>(4,
                    new SquareMatrix<Real>(new Real[2, 2] { 
                        { Real.one, Real.zero },
                        { Real.zero, Real.zero }
                    }),
                    new SquareMatrix<Real>(new Real[2, 2] {
                        { Real.zero, Real.one },
                        { Real.zero, Real.zero }
                    }),
                    new SquareMatrix<Real>(new Real[2, 2] {
                        { Real.zero, Real.zero },
                        { Real.one, Real.zero }
                    }),
                    new SquareMatrix<Real>(new Real[2, 2] {
                        { Real.zero, Real.zero },
                        { Real.zero, Real.one }
                    })
                    );
                VectorSpace<Polynomial<Real>, Real> P2 = new VectorSpace<Polynomial<Real>, Real>(3,
                    new Polynomial<Real>(3, Real.one),
                    new Polynomial<Real>(3, Real.zero, Real.one),
                    new Polynomial<Real>(3, Real.zero, Real.zero, Real.one)
                    );
                Transform<Polynomial<Real>, SquareMatrix<Real>, Real> T = new Transform<Polynomial<Real>, SquareMatrix<Real>, Real>(P2, M22,
                    p => new SquareMatrix<Real>(new Real[2, 2] {
                    { p[2].Add(p[1]) as Real, p[2].Add(p[0]) as Real },
                    { p[1].Add(p[0].AddInverse()) as Real, p[1].Add(p[0]) as Real },
                    })
                    );
            }
            catch(MatrixSizeException ex)
            {
                Console.WriteLine("Matrix size exception.");
                Console.WriteLine(ex.Message);
            }
        }
    }
}
