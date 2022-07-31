using System;


namespace Linear_Algebra
{
    class Real : Field
    {
        public static readonly Real zero = new Real(0f);
        public static readonly Real one = new Real(1f);

        private readonly float value;

        public Real(float value)
        {
            this.value = value;
        }

        public Field Add(Field other)
        {
            return new Real(value + (other as Real).value);
        }

        public Field Multiply(Field other)
        {
            return new Real(value * (other as Real).value);
        }

        public Field Zero() { return zero; }

        public Field One() { return one; }

        public Field AddInverse()
        {
            return new Real(-value);
        }

        public Field MultInverse()
        {
            return new Real(1 / value);
        }

        public override bool Equals(object obj)
        {
            return Math.Abs(value - (obj as Real).value) < 0.0001f;
        }

        public override string ToString()
        {
            if (value.ToString().Equals("-0")) { return "0"; }
            return value.ToString();
        }
    }
}
