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

    interface InnerProduct<F> : Vector<F> where F : Complex
    {
        // @post u.DotProduct(v) == v.DotProduct(u).Complement()
        // @post vec1 == vec2 $implies $ret >= 0 (in particular $ret is Real
        // @post vec1 == vec2 $implies $ret == 0 iff vec1, vec2 == Zero()
        // @post (a * u + b * v).DotProduct(w) == (a * u).DotProduct(w) + (b * v).DotProduct(w)
        F DotProduct(InnerProduct<F> vec);

        public static F operator * (InnerProduct<F> vec1, InnerProduct<F> vec2)
        {
            return vec1.DotProduct(vec2);
        }

        Vector<F> Normalize()
        {
            return Multiply((F)DotProduct(this).Magnitude().MultInverse());
        }
    }
}
