namespace Linear_Algebra
{
    class Polynomial<F> : Vector<F> where F : Field
    {
        private readonly ColumnVector<F> coefficients;
        public readonly int maxDeg;

        // @pre no null coefficients
        public Polynomial(params F[] coeffs) : this(new ColumnVector<F>(coeffs)) { }

        // @pre no null coefficients
        // @pre length >= coeffs.length >= 1
        public Polynomial(int maxDeg, params F[] coeffs)
        {
            this.maxDeg = maxDeg;
            F[] coeffsCopy = new F[maxDeg + 1];
            F zero = (F)coeffs[0].Zero();
            for (int i = 0; i < coeffs.Length && i <= maxDeg; i++)
            {
                coeffsCopy[i] = coeffs[i];
            }
            for (int i = coeffs.Length; i <= maxDeg; i++)
            {
                coeffsCopy[i] = zero;
            }
            coefficients = new ColumnVector<F>(coeffsCopy);
        }

        // @pre no null coefficients
        public Polynomial(ColumnVector<F> coeffs)
        {
            coefficients = coeffs.Clone();
            maxDeg = coeffs.length - 1;
        }

        // @pre 0 <= index < length
        public F this[int index]
        {
            get { return coefficients[index]; }
        }

        // @pre scalar != null
        public F ValueOf(F scalar)
        {
            F res = (F)scalar.Zero();
            F pow = (F)scalar.One();
            foreach(F coeff in coefficients)
            {
                res = (F)(res + coeff * pow);
                pow = (F)(pow * scalar);
            }
            return res;
        }

        // @pre matrix != null
        public SquareMatrix<F> ValueOf(SquareMatrix<F> matrix)
        {
            SquareMatrix<F> res = matrix.Zero() as SquareMatrix<F>;
            SquareMatrix<F> pow = matrix.Identity();
            foreach(F coeff in coefficients)
            {
                res = (res + pow.Multiply(coeff)) as SquareMatrix<F>;
                pow *= matrix;
            }
            return res;
        }

        // @pre T != null
        public Transform<V, F> ValueOf<V>(Transform<V, F> transform) where V : Vector<F>
        {
            Transform<V, F> res = transform.Zero() as Transform<V, F>;
            Transform<V, F> pow = transform.Identity();
            foreach(F coeff in coefficients)
            {
                res = (res + pow.Multiply(coeff)) as Transform<V, F>;
                pow *= transform;
            }
            return res;
        }

        public Vector<F> Add(Vector<F> vector)
        {
            Polynomial<F> other = vector as Polynomial<F>;
            Vector<F> coeffs1 = coefficients;
            Vector<F> coeffs2 = other.coefficients;
            return new Polynomial<F>((coeffs1 + coeffs2) as ColumnVector<F>);
        }

        public Vector<F> Multiply(F scalar)
        {
            return new Polynomial<F>(scalar * (coefficients as Vector<F>) as ColumnVector<F>);
        }

        public Vector<F> AddInverse()
        {
            return new Polynomial<F>((coefficients as Vector<F>).AddInverse() as ColumnVector<F>);
        }

        public Vector<F> Zero()
        {
            return new Polynomial<F>((coefficients as Vector<F>).Zero() as ColumnVector<F>);
        }

        public Polynomial<F> Clone()
        {
            return new Polynomial<F>(coefficients);
        }

        public ColumnVector<F> ToColumnVector() { return coefficients; }

        public override string ToString()
        {
            if(maxDeg == 0) { return "0"; }
            F coeff = this[0];
            string res = coeff.Equals(coeff.Zero()) ? "" : coeff.ToString();
            for (int i = 1; i <= maxDeg; i++)
            {
                coeff = this[i];
                if (coeff.Equals(coeff.Zero())) { continue; }
                res += res == "" ? "" : " + ";
                res += term(coeff, i);
            }
            return res == "" ? "0" : res;
        }

        private string term(F coeff, int pow)
        {
            if (coeff.Equals(coeff.One()))
            {
                if(pow == 1)
                {
                    return "x";
                }
                return string.Format("x^{0}", pow);
            }
            if (coeff.ToString().Length == 1)
            {
                if (pow == 1)
                {
                    return string.Format("{0}x", coeff);
                }
                return string.Format("{0}x^{1}", coeff, pow);
            }
            if (pow == 1)
            {
                return string.Format("({0})x", coeff);
            }
            return string.Format("({0})x^{1}", coeff, pow);
        }

        public int Length()
        {
            return coefficients.length;
        }
    }
}
