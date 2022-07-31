namespace Linear_Algebra
{
    interface Vector<F> where F : Field
    {
        Vector<F> Add(Vector<F> vector);
        Vector<F> Multiply(F scalar);
        Vector<F> Zero();
        Vector<F> AddInverse();
        int Length();
        ColumnVector<F> ToColumnVector();

        public static Vector<F> operator + (Vector<F> vec1, Vector<F> vec2)
        {
            return vec1.Add(vec2);
        }

        public static Vector<F> operator - (Vector<F> vec1, Vector<F> vec2)
        {
            return vec1.Add(vec2.AddInverse());
        }

        public static Vector<F> operator * (F scalar, Vector<F> vec)
        {
            return vec.Multiply(scalar);
        }
    }
}
