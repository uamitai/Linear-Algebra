namespace Linear_Algebra
{
    class ModuloP
    {
        private readonly int p;

        // @pre p must be positive and prime
        public ModuloP(int p)
        {
            this.p = p;
        }

        public Integer CreateInt(int value)
        {
            return new Integer(this, value);
        }

        public Integer Zero() { return CreateInt(0); }
        public Integer One() { return CreateInt(1); }

        public bool Equals(ModuloP obj)
        {
            return p == obj.p;
        }

        public class Integer : Field
        {
            private readonly ModuloP modulus;
            private readonly int value;

            public Integer(ModuloP modulus, int value)
            {
                this.modulus = modulus;
                this.value = value >= 0 ? value % modulus.p : modulus.p - (-value % modulus.p);
            }

            // @pre this and other are modulo the same prime
            public Field Add(Field other)
            {
                return modulus.CreateInt((value + (other as Integer).value) % modulus.p);
            }

            // @pre this and other are modulo the same prime
            public Field Multiply(Field other)
            {
                return modulus.CreateInt((value * (other as Integer).value) % modulus.p);
            }

            public Field Zero() { return modulus.Zero(); }

            public Field One() { return modulus.One(); }

            public Field AddInverse()
            {
                return modulus.CreateInt((modulus.p - value) % modulus.p);
            }

            public Field MultInverse()
            {
                return modulus.CreateInt(ModPower(value, modulus.p - 2, modulus.p));
            }

            public override string ToString()
            {
                return value.ToString();
            }

            public override bool Equals(object obj)
            {
                Integer n = obj as Integer;
                return value % n.modulus.p == n.value || n.value % modulus.p == value;
            }

            private static int ModPower(int a, int b, int c)
            {
                int res = 1;
                while (b > 0)
                {
                    if (b % 2 == 1)
                    {
                        res = (res * a) % c;
                        b--;
                    }
                    a = (a * a) % c;
                    b /= 2;
                }
                return res;
            }
        }
    }
}
