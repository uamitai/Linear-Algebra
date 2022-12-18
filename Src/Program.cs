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
                Bool one = Bool.one;
                Bool zero = Bool.zero;
                SquareMatrix<Bool> A = new SquareMatrix<Bool>(new Bool[5,5] {
                { zero, zero, zero, one, zero },
                { one, zero, one, zero, zero },
                { zero, zero, zero, zero, one },
                { zero, zero, one, zero, one },
                { one, one, zero, zero, zero }
                });
                Console.WriteLine(A.Power(4));
            }
            catch(MatrixSizeException ex)
            {
                Console.WriteLine("Matrix size exception.");
                Console.WriteLine(ex.Message);
            }
        }
    }
}
