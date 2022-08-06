using System;


namespace Linear_Algebra
{
    class Real : Complex
    {
        public const float eps = 0.005f;
        public static readonly new Real zero = new Real(0f);
        public static readonly new Real one = new Real(1f);

        public Real(float value) : base(value, 0) { }

        public override Field Add(Field other)
        {
            float res = real + (other as Real).real;
            if(Math.Abs(res) < eps) { return zero; }
            return new Real(res);
        }

        public override Field Multiply(Field other)
        {
            return new Real(real * (other as Real).real);
        }

        public override Field Zero() { return zero; }

        public override Field One() { return one; }

        public override Field AddInverse()
        {
            return new Real(-real);
        }

        public override Field MultInverse()
        {
            return new Real(1 / real);
        }

        public override Complex Complement()
        {
            return this;
        }
    }
}
