using System;


namespace Linear_Algebra
{
    class Program
    {
        static void Main(string[] args)
        {
            try
            {
            }
            catch(MatrixSizeException ex)
            {
                Console.WriteLine("Matrix size exception.");
                Console.WriteLine(ex.Message);
            }
        }
    }
}
