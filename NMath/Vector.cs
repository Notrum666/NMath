using System;

namespace NMath
{
    public class Vector
    {
        public double[] values { get; private set; }
        public int size { get { return values.Length; } }
        public double magnitude { get { double val = 0; foreach (double v in values) val += v * v; return Math.Sqrt(val); } }
        public double sqrMagnitude { get { double val = 0; foreach (double v in values) val += v * v; return val; } }
        public double this[int i] { get { return values[i]; } set { values[i] = value; } }
        public Vector(int size)
        {
            values = new double[size];
        }
        public Vector(params double[] values)
        {
            this.values = values;
        }
        public Vector(Vector v)
        {
            values = new double[v.size];
            for (int i = 0; i < values.Length; i++)
                values[i] = v.values[i];
        }
        public static double operator*(Vector lhs, Vector rhs)
        {
            if (lhs.size != rhs.size)
                throw new ArgumentException("Vectors must be same size.");
            double result = 0;
            for (int i = 0; i < lhs.size; i++)
                result += lhs[i] * rhs[i];
            return result;
        }
        public static Vector operator+(Vector lhs, Vector rhs)
        {
            if (lhs.size != rhs.size)
                throw new ArgumentException("Vectors must be same size.");
            double[] values = new double[lhs.size];
            for (int i = 0; i < lhs.size; i++)
                values[i] = lhs[i] + rhs[i];
            return new Vector(values);
        }
        public static Vector operator-(Vector lhs, Vector rhs)
        {
            if (lhs.size != rhs.size)
                throw new ArgumentException("Vectors must be same size.");
            double[] values = new double[lhs.size];
            for (int i = 0; i < lhs.size; i++)
                values[i] = lhs[i] - rhs[i];
            return new Vector(values);
        }
        public override string ToString()
        {
            return "(" + string.Join(", ", values) + ")";
        }
        public string ToString(string format)
        {
            string res = "(";
            for (int i = 0; i < size; i++)
                res += values[i].ToString(format) + (i == size - 1 ? "" : ", ");
            return res + ")";
        }
        public override bool Equals(object obj)
        {
            if (obj is Vector)
            {
                if ((obj as Vector).size != size)
                    return false;
                double[] otherValues = (obj as Vector).values;
                for (int i = 0; i < size; i++)
                    if (Math.Abs(otherValues[i] - values[i]) > SLAE.epsilon)
                        return false;
                return true;
            }
            return base.Equals(obj);
        }
    }
    public class Vector2 : Vector
    {
        public double X { get { return values[0]; } set { values[0] = value; } }
        public double Y { get { return values[1]; } set { values[1] = value; } }
        public new readonly int size = 2;
        public Vector2() : base(2)
        {

        }
        public Vector2(double X, double Y) : base(X, Y)
        {

        }
        public double vecMul(Vector2 vec)
        {
            return this % vec;
        }
        public static double operator%(Vector2 lhs, Vector2 rhs)
        {
            return lhs.X * rhs.Y - lhs.Y * rhs.X;
        }
    }
    public class Vector3 : Vector
    {
        public double X { get { return values[0]; } set { values[0] = value; } }
        public double Y { get { return values[1]; } set { values[1] = value; } }
        public double Z { get { return values[2]; } set { values[2] = value; } }
        public new readonly int size = 3;
        public Vector3() : base(3)
        {

        }
        public Vector3(double X, double Y, double Z) : base(X, Y, Z)
        {

        }
        public Vector3 vecMul(Vector3 vec)
        {
            return this % vec;
        }
        public static Vector3 operator%(Vector3 lhs, Vector3 rhs)
        {
            return new Vector3(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X);
        }
    }
}