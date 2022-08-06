﻿using System;


namespace Linear_Algebra
{
    class Program
    {
        static void Main(string[] args)
        {
            try
            {
                SquareMatrix<Real> A = new SquareMatrix<Real>(new Real[3, 3] {
                { new Real(2.92f), new Real(0.86f), new Real(-1.15f) },
                { new Real(0.86f), new Real(6.51f), new Real(3.32f) },
                { new Real(-1.15f), new Real(3.32f), new Real(4.57f) }
                });

                SquareMatrix<Real> Q = new SquareMatrix<Real>(new Real[3, 3] {
                { new Real(0), new Real(-0.8f), new Real(-0.6f) },
                { new Real(0.8f), new Real(-0.36f), new Real(0.48f) },
                { new Real(0.6f), new Real(0.48f), new Real(-0.64f) }
                });

                SquareMatrix<Real> D = SquareMatrix<Real>.Diag(new Real(9), new Real(4), new Real(1));

                Console.WriteLine((Q * D * Q.Transpose() as Vector<Real>) - A);
            }
            catch(MatrixSizeException ex)
            {
                Console.WriteLine("Matrix size exception.");
                Console.WriteLine(ex.Message);
            }
        }
    }
}
