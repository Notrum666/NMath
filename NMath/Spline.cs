using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace NMath
{
    public class Spline
    {
        public double LeftBorder { get { return points[0]; } }
        public double RightBorder { get { return points[intervalsCount]; } }
        private double[] points;
        private double[] coeffs;
        private int intervalsCount;
        public double Get(double arg)
        {
            for (int i = 0; i < intervalsCount; i++)
                if (arg >= points[i] && arg <= points[i + 1])
                {
                    return coeffs[i * 4] * arg * arg * arg + coeffs[i * 4 + 1] * arg * arg + coeffs[i * 4 + 2] * arg + coeffs[i * 4 + 3];
                }
            throw new ArgumentOutOfRangeException();
        }
        public static Spline Interpolate(Dictionary<double, double> pointsAndValues)
        {
            List<KeyValuePair<double, double>> points = pointsAndValues.ToList();
            int size = points.Count - 1;
            if (size < 1)
                throw new ArgumentException("There must be at least 2 points.");
            points.OrderBy(point => point.Key);
            for (int i = 0; i < size; i++)
                for (int j = i + 1; j <= size; j++)
                    if (points[i].Key == points[j].Key)
                        throw new ArgumentException("There can't be two values at the same points.");
            Matrix_Full mat = new Matrix_Full(size);
            Vector vec = new Vector(size);
            mat[0, 0] = 6 * points[0].Key;
            mat[0, 1] = 2;
            int offset = 0;
            for (int i = 1; i < size; i++, offset += 4)
            {
                mat[offset + 1, offset] = points[i - 1].Key * points[i - 1].Key * points[i - 1].Key;
                mat[offset + 1, offset + 1] = points[i - 1].Key * points[i - 1].Key;
                mat[offset + 1, offset + 2] = points[i - 1].Key;
                mat[offset + 1, offset + 3] = 1;

                mat[offset + 2, offset] = points[i].Key * points[i].Key * points[i].Key;
                mat[offset + 2, offset + 1] = points[i].Key * points[i].Key;
                mat[offset + 2, offset + 2] = points[i].Key;
                mat[offset + 2, offset + 3] = 1;

                mat[offset + 3, offset] = 6 * points[i].Key;
                mat[offset + 3, offset + 1] = 2;
                mat[offset + 3, offset + 4] = -6 * points[i].Key;
                mat[offset + 3, offset + 5] = -2;

                mat[offset + 4, offset] = 3 * points[i].Key * points[i].Key;
                mat[offset + 4, offset + 1] = 2 * points[i].Key;
                mat[offset + 4, offset + 2] = 1;
                mat[offset + 4, offset + 4] = -3 * points[i].Key * points[i].Key;
                mat[offset + 4, offset + 5] = -2 * points[i].Key;
                mat[offset + 4, offset + 6] = -1;

                vec[offset + 1] = points[i - 1].Value;
                vec[offset + 2] = points[i].Value;
            }

            mat[offset + 1, offset] = points[size - 1].Key * points[size - 1].Key * points[size - 1].Key;
            mat[offset + 1, offset + 1] = points[size - 1].Key * points[size - 1].Key;
            mat[offset + 1, offset + 2] = points[size - 1].Key;
            mat[offset + 1, offset + 3] = 1;
            vec[offset + 1] = points[size - 1].Value;

            mat[offset + 2, offset] = points[size].Key * points[size].Key * points[size].Key;
            mat[offset + 2, offset + 1] = points[size].Key * points[size].Key;
            mat[offset + 2, offset + 2] = points[size].Key;
            mat[offset + 2, offset + 3] = 1;
            vec[offset + 2] = points[size].Value;

            mat[offset + 3, offset] = 6 * points[size].Key;
            mat[offset + 3, offset + 1] = 2;

            //Console.WriteLine(mat.ToString("F4"));
            //Console.WriteLine(vec);
            SLAE.Solve_Gauss(mat, vec);
            //Console.WriteLine(vec.ToString());

            Spline spline = new Spline();
            spline.points = points.Select(point => point.Key).ToArray();
            spline.coeffs = vec.values;
            spline.intervalsCount = size;

            return spline;
        }
    }
}
