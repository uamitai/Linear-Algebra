namespace Linear_Algebra
{
    class Real : Complex
    {
        public const float eps = 0.0003f;
        public static readonly Real zero = new Real(0f);
        public static readonly Real one = new Real(1f);

        public Real(float value) : base(value, 0) { }

        public override Field Add(Field other)
        {
            Real res = new Real(real + (other as Real).real);
            return res.Equals(zero) ? zero : res.Equals(one) ? one : res;
        }

        public override Field Multiply(Field other)
        {
            Real res = new Real(real * (other as Real).real);
            return res.Equals(one) ? one : res;
        }

        public override Field AddInverse()
        {
            return new Real(-real);
        }

        public override Field MultInverse()
        {
            return new Real(1 / real);
        }

        public override Real Conjugate()
        {
            return this;
        }
    }
}
