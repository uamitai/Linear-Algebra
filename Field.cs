namespace Linear_Algebra
{
    interface Field
    {
        Field Add(Field other);
        Field Multiply(Field other);
        Field Zero();
        Field One();
        Field AddInverse();
        Field MultInverse();

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
