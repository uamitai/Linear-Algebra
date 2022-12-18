namespace Linear_Algebra
{
    class Bool : Field
    {
        public static readonly Bool one = new Bool(1);
        public static readonly Bool zero = new Bool(0);

        private readonly bool value;

        public Bool(bool value)
        {
            this.value = value;
        }

        public Bool(int value)
        {
            this.value = value != 0;
        }

        public Field Add(Field other)
        {
            return new Bool(value || (other as Bool).value);
        }

        public Field Multiply(Field other)
        {
            return new Bool(value && (other as Bool).value);
        }

        public Field One()
        {
            return one;
        }

        public Field Zero()
        {
            return zero;
        }

        public Field AddInverse()
        {
            return Not();
        }

        public Field MultInverse()
        {
            return Not();
        }

        public Bool Not()
        {
            return Equals(one) ? zero : one;
        }

        public override bool Equals(object obj)
        {
            return value == (obj as Bool).value;
        }

        public override string ToString()
        {
            return value ? "1" : "0";
        }
    }
}
