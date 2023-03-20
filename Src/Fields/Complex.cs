using System;


namespace Linear_Algebra
{
    class Complex : Field
    {
        protected readonly float real;
        private readonly float imaginary;

        public Complex(float a, float b)
        {
            real = a;
            imaginary = b;
        }

        public virtual Field Add(Field other)
        {
            Complex complex = other as Complex;
            float re = real + complex.real;
            float im = imaginary + complex.imaginary;
            if(Math.Abs(re) < Real.eps) { re = 0; }
            if(Math.Abs(im) < Real.eps) { im = 0; }
            if(Math.Abs(re - 1) < Real.eps) { re = 1; }
            return new Complex(re, im);
        }

        public virtual Field Multiply(Field other)
        {
            Complex complex = other as Complex;
            return new Complex(real * complex.real - imaginary * complex.imaginary, real * complex.imaginary + imaginary * complex.real);
        }

        public Field Zero() { return Real.zero; }

        public Field One() { return Real.one; }

        public virtual Field AddInverse()
        {
            return new Complex(-real, -imaginary);
        }

        public virtual Field MultInverse()
        {
            return new Complex(real / SqrMagnitude(), -imaginary / SqrMagnitude());
        }

        private float SqrMagnitude()
        {
            return real * real + imaginary * imaginary;
        }

        public float Magnitude()
        {
            return (float)Math.Sqrt(SqrMagnitude());
        }

        public virtual Complex Conjugate()
        {
            return new Complex(real, -imaginary);
        }

        public override bool Equals(object obj)
        {
            Complex complex = obj as Complex;
            return Math.Abs(real - complex.real) < Real.eps && Math.Abs(imaginary - complex.imaginary) < Real.eps;
        }

        public override string ToString()
        {
            string im;
            switch (imaginary)
            {
                case 0: 
                    im = "";
                    break;
                case 1:
                    im = "i";
                    break;
                case -1:
                    im = "-i";
                    break;
                default:
                    im = imaginary + "i";
                    break;
            }
            if(real == 0)
            {
                return imaginary == 0 ? "0" : im;
            }
            return imaginary > 0 ? real + "+" + im : real + im;
        }
    }
}
