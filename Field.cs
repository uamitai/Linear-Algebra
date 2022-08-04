namespace Linear_Algebra
{
    interface Field
    {
        // @post a.Add(b) == b.Add(a)
        // @post (a.Add(b)).Add(c) == a.Add(b.Add(c)
        Field Add(Field other);

        // @post a.Multiply(b) == b.Multiply(a)
        // @post (a.Multiply(b)).Multiply(c) == a.Multiply(b.Multiply(c)
        // @post a.Multiply(b.Add(c)) == a.Multiply(b).Add(a.Multiply(c))
        Field Multiply(Field other);

        // @post this.Add(Zero()) == this
        Field Zero();

        // @post this.Multiply(One()) == this
        Field One();

        // @post this.Add(AddInverse()) == Zero()
        Field AddInverse();

        // @post this.Multiply(MultInverse()) == One()
        Field MultInverse();

        bool IsZero() { return Equals(Zero()); }

        public static Field operator + (Field a, Field b)
        {
            return a.Add(b);
        }

        public static Field operator - (Field a, Field b)
        {
            return a.Add(b.AddInverse());
        }

        public static Field operator * (Field a, Field b)
        {
            return a.Multiply(b);
        }
    }
}
