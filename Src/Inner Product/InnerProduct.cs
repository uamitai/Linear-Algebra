using System;


namespace Linear_Algebra
{
    interface InnerProduct<F> : Vector<F> where F : Complex
    {
        // @post u.DotProduct(v) == v.DotProduct(u).Complement()
        // @post vec1 == vec2 $implies $ret >= 0 (in particular $ret is Real
        // @post vec1 == vec2 $implies $ret == 0 iff vec1, vec2 == Zero()
        // @post (a * u + b * v).DotProduct(w) == (a * u).DotProduct(w) + (b * v).DotProduct(w)
        F InnerProduct(InnerProduct<F> vector);

        public static F operator * (InnerProduct<F> vec1, InnerProduct<F> vec2)
        {
            return vec1.InnerProduct(vec2);
        }

        // @pre this != Zero()
        Vector<F> Normalize()
        {
            float mag = InnerProduct(this).Magnitude();
            return Multiply((F)new Real((float)Math.Sqrt(mag)).MultInverse());
        }
    }
}
