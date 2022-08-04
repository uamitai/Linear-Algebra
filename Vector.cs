namespace Linear_Algebra
{
    interface Vector<F> where F : Field
    {
        // @post u.Add(v) == v.Add(u)
        // @post (u.Add(v)).Add(w) == u.Add(v.Add(w)
        Vector<F> Add(Vector<F> vector);

        // @post Multiply(a * b) == Multiply(a).Multiply(b) == Multiply(b).Multiply(a)
        // @post Multiply(a + b) == Multiply(a).Add(Multiply(b))
        // @post (u.Add(v)).Multiply(a) == u.Multiply(a).Add(v.Multiply(a))
        Vector<F> Multiply(F scalar);

        // @post this.Add(Zero()) == this
        Vector<F> Zero();

        // @post this.Add(AddInverse()) == Zero()
        Vector<F> AddInverse();

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
