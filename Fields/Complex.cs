using System;


namespace Linear_Algebra
{
    class Complex : Field
    {
        public static readonly Complex zero = new Complex(0f, 0f);
        public static readonly Complex one = new Complex(1f, 0f);

        private readonly float re;
        private readonly float im;

        public Complex(float a, float b)
        {
            re = a;
            im = b;
        }

        public Complex(float a) : this(a, 0) { }

        public Field Add(Field other)
        {
            Complex complex = other as Complex;
            return new Complex(re + complex.re, im + complex.im);
        }

        public Field Multiply(Field other)
        {
            Complex complex = other as Complex;
            return new Complex(re * complex.re - im * complex.im, re * complex.im + im * complex.re);
        }

        public Field Zero() { return zero; }

        public Field One() { return one; }

        public Field AddInverse()
        {
            return new Complex(-re, -im);
        }

        public Field MultInverse()
        {
            return new Complex(re / sqrMagnitude(), -im / sqrMagnitude());
        }

        private float sqrMagnitude()
        {
            return re * re + im * im;
        }

        public override bool Equals(object obj)
        {
            Complex complex = obj as Complex;
            return Math.Abs(re - complex.re) < 0.0001f && Math.Abs(im - complex.im) < 0.0001f;
        }

        public override string ToString()
        {
            string imaginary;
            switch (im)
            {
                case 0: 
                    imaginary = "";
                    break;
                case 1:
                    imaginary = "i";
                    break;
                case -1:
                    imaginary = "-i";
                    break;
                default:
                    imaginary = im + "i";
                    break;
            }
            if(re == 0)
            {
                return im == 0 ? "0" : imaginary;
            }
            return im > 0 ? re + "+" + imaginary : re + imaginary;
        }
    }
}
