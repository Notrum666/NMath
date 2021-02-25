using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace NMath
{
    public class Matrix_Sparse : Matrix
    {
        public int[] ig { get; private set; }
        public int[] jg { get; private set; }
        public double[] di { get; private set; }
        public double[] ggl { get; private set; }
        public double[] ggu { get; private set; }
        public override int size { get; protected set; }
        public Matrix_Sparse(int size, int[] ig, int[] jg, double[] di, double[] ggu, double[] ggl)
        {
            this.size = size;
            this.ig = ig;
            this.jg = jg;
            this.di = di;
            this.ggu = ggu;
            this.ggl = ggl;
        }
        public static Matrix_Sparse Parse(string size, string ig, string jg, string di, string ggu, string ggl)
        {
            return new Matrix_Sparse(int.Parse(size),
                ig.Split(' ').Select(value => int.Parse(value)).ToArray(),
                jg.Split(' ').Select(value => int.Parse(value)).ToArray(),
                di.Split(' ').Select(value => double.Parse(value)).ToArray(),
                ggu.Split(' ').Select(value => double.Parse(value)).ToArray(),
                ggl.Split(' ').Select(value => double.Parse(value)).ToArray());
        }
        public static Vector operator *(Matrix_Sparse mat, Vector vec)
        {
            if (mat.size != vec.size)
                throw new ArgumentException("Matrix and vector must be same size.");
            int size = mat.size;
            int[] ig = mat.ig;
            int[] jg = mat.jg;

            double[] di = mat.di;
            double[] ggu = mat.ggu;
            double[] ggl = mat.ggl;

            Vector result = new Vector(size);

            for (int i = 0; i < size; i++)
            {
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    result[jg[j]] += vec[i] * ggu[j];
                    result[i] += vec[jg[j]] * ggl[j];
                }
                result[i] += di[i] * vec[i];
            }
            return result;
        }
    }
    public class MatrixSymm_Sparse : MatrixSymm
    {
        public int[] ig { get; private set; }
        public int[] jg { get; private set; }
        public double[] di { get; private set; }
        public double[] gg { get; private set; }
        public override int size { get; protected set; }
        public MatrixSymm_Sparse(int size, int[] ig, int[] jg, double[] di, double[] gg)
        {
            this.size = size;
            this.ig = ig;
            this.jg = jg;
            this.di = di;
            this.gg = gg;
        }
        public static MatrixSymm_Sparse GetGilbert(int size)
        {
            double[] di = new double[size];
            int[] ig = new int[size + 1];

            int ggSize = size * (size - 1) / 2;
            int[] jg = new int[ggSize];
            double[] gg = new double[ggSize];

            ig[0] = 0;

            int i, j;
            for (i = 0; i < size; i++)
            {
                di[i] = 1.0 / (i + i + 1);

                ig[i + 1] = ig[i] + i;
                for (j = 0; j < i; j++)
                {
                    jg[ig[i] + j] = j;
                    gg[ig[i] + j] = 1.0 / (i + j + 1);
                }
            }

            return new MatrixSymm_Sparse(size, ig, jg, di, gg);
        }
        public static MatrixSymm_Sparse Parse(string size, string ig, string jg, string di, string gg)
        {
            return new MatrixSymm_Sparse(int.Parse(size),
                ig.Split(' ').Select(value => int.Parse(value)).ToArray(),
                jg.Split(' ').Select(value => int.Parse(value)).ToArray(),
                di.Split(' ').Select(value => double.Parse(value)).ToArray(),
                gg.Split(' ').Select(value => double.Parse(value)).ToArray());
        }
        public static Vector operator *(MatrixSymm_Sparse mat, Vector vec)
        {
            int size = mat.size;
            int[] ig = mat.ig;
            int[] jg = mat.jg;

            double[] di = mat.di;
            double[] gg = mat.gg;

            Vector result = new Vector(size);

            for (int i = 0; i < size; i++)
            {
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    result[jg[j]] += vec[i] * gg[j];
                    result[i] += vec[jg[j]] * gg[j];
                }
                result[i] += di[i] * vec[i];
            }
            return result;
        }
        public override string ToString()
        {
            string res = "";
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (i == j)
                        res += di[i].ToString() + " ";
                    else
                    {
                        int k;
                        if (i < j)
                        {
                            for (k = ig[j]; k < ig[j + 1]; k++)
                                if (jg[k] == i)
                                {
                                    res += gg[k].ToString() + " ";
                                    break;
                                }
                            if (k == ig[j + 1])
                                res += "0 ";
                        }
                        else
                        {
                            for (k = ig[i]; k < ig[i + 1]; k++)
                                if (jg[k] == j)
                                {
                                    res += gg[k].ToString() + " ";
                                    break;
                                }
                            if (k == ig[i + 1])
                                res += "0 ";
                        }
                    }
                }
                res += "\n";
            }
            return res;
        }
        public string ToString(string format)
        {
            string res = "";
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (i == j)
                        res += di[i].ToString(format) + " ";
                    else
                    {
                        int k;
                        if (i < j)
                        {
                            for (k = ig[j]; k < ig[j + 1]; k++)
                                if (jg[k] == i)
                                {
                                    res += gg[k].ToString(format) + " ";
                                    break;
                                }
                            if (k == ig[j + 1])
                                res += 0.ToString(format) + " ";
                        }
                        else
                        {
                            for (k = ig[i]; k < ig[i + 1]; k++)
                                if (jg[k] == j)
                                {
                                    res += gg[k].ToString(format) + " ";
                                    break;
                                }
                            if (k == ig[i + 1])
                                res += 0.ToString(format) + " ";
                        }
                    }
                }
                res += "\n";
            }
            return res;
        }
    }
}
