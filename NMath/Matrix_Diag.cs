using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace NMath
{
    public class Matrix_Diag : Matrix
    {
        public double[,] values { get; private set; }
        public int[] offsets { get; private set; }
        public int diagCount { get { return offsets.Length; } }
        public override int size { get { return values.GetLength(1); } protected set { } }
        public Matrix_Diag(int size, int[] offsets)
        {
            this.offsets = offsets;
            values = new double[offsets.Length, size];
        }
        public Matrix_Diag(int[] offsets, double[,] values)
        {
            this.offsets = offsets;
            this.values = values;
        }
        public Matrix_Diag(Matrix_Diag mat)
        {
            int size = mat.size;
            int diagCount = mat.diagCount;
            values = new double[diagCount, size];
            offsets = new int[diagCount];
            double[,] origValues = mat.values;
            int[] origOffsets = mat.offsets;
            for (int i = 0; i < diagCount; i++)
            {
                offsets[i] = origOffsets[i];
                for (int j = 0; j < size; j++)
                    values[i, j] = origValues[i, j];
            }
        }
        public static Matrix_Diag Parse(Matrix mat)
        {
            if (mat is Matrix_Full)
                return Parse(mat as Matrix_Full);
            if (mat is Matrix_Diag)
                return new Matrix_Diag(mat as Matrix_Diag);
            if (mat is Matrix_Sparse)
                return Parse(mat as Matrix_Sparse);
            throw new NotImplementedException("This parsing method is not implemented yet.");
        }
        public static Matrix_Diag Parse(MatrixSymm mat)
        {
            if (mat is MatrixSymm_Sparse)
                return Parse(mat as MatrixSymm_Sparse);
            if (mat is MatrixSymm_Diag)
                return Parse(mat as MatrixSymm_Diag);
            if (mat is MatrixSymm_Full)
                return Parse(mat as MatrixSymm_Full);
            throw new NotImplementedException("This parsing method is not implemented yet.");
        }
        public static Matrix_Diag Parse(Matrix_Full mat)
        {
            int size = mat.size;
            double[,] values = mat.values;
            List<int> offsets = new List<int>();
            List<double[]> res = new List<double[]>();
            int diagIndex;
            int offset;
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (values[i, j] != 0)
                    {
                        offset = j - i;
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset);
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][i] = values[i, j];
                    }
            int diagCount = offsets.Count;
            double[,] result = new double[diagCount, size];
            for (int i = 0; i < diagCount; i++)
                for (int j = 0; j < size; j++)
                    result[i, j] = res[i][j];
            return new Matrix_Diag(offsets.ToArray(), result);
        }
        public static Matrix_Diag Parse(Matrix_Sparse mat)
        {
            int size = mat.size;
            int[] ig = mat.ig;
            int[] jg = mat.jg;
            double[] ggu = mat.ggu;
            double[] ggl = mat.ggl;
            List<int> offsets = new List<int>();
            List<double[]> res = new List<double[]>();
            int diagIndex;
            int offset;
            for (int i = 0; i < size; i++)
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    offset = jg[i] - i;
                    if (ggl[j] != 0)
                    {
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset);
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][i] = ggl[j];
                    }
                    offset = -offset;
                    if (ggu[j] != 0)
                    {
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset) ;
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][jg[i]] = ggu[j];
                    }
                }
            int diagCount = offsets.Count;
            double[,] result = new double[diagCount, size];
            for (int i = 0; i < diagCount; i++)
                for (int j = 0; j < size; j++)
                    result[i, j] = res[i][j];
            return new Matrix_Diag(offsets.ToArray(), result);
        }
        public static Matrix_Diag Parse(MatrixSymm_Full mat)
        {
            int size = mat.size;
            double[][] values = mat.values;
            List<int> offsets = new List<int>();
            List<double[]> res = new List<double[]>();
            int diagIndex;
            int offset;
            for (int i = 0; i < size; i++)
                for (int j = 0; j <= i; j++)
                    if (values[i][j] != 0)
                    {
                        offset = j - i;
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset) ;
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][i] = values[i][j];
                        if (offset == 0)
                            continue;
                        offset = -offset;
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset) ;
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][j - offset] = values[i][j];
                    }
            int diagCount = offsets.Count;
            double[,] result = new double[diagCount, size];
            for (int i = 0; i < diagCount; i++)
                for (int j = 0; j < size; j++)
                    result[i, j] = res[i][j];
            return new Matrix_Diag(offsets.ToArray(), result);
        }
        public static Matrix_Diag Parse(MatrixSymm_Diag mat)
        {
            int size = mat.size;
            int diagCount = mat.diagCount;
            double[,] values = mat.values;
            int[] offsets = mat.offsets;
            int resDiagCount = diagCount * 2 - (offsets.Contains(0) ? 1 : 0);
            int[] resOffsets = new int[resDiagCount];
            double[,] res = new double[resDiagCount, size];
            int inverseDiagIndex;
            for (int i = 0; i < diagCount; i++)
            {
                inverseDiagIndex = resDiagCount - 1 - i;
                resOffsets[i] = offsets[i];
                resOffsets[inverseDiagIndex] = -offsets[i];
                for (int j = 0; j < size; j++)
                {
                    res[i, j] = values[i, j];
                    res[inverseDiagIndex, j + offsets[i]] = values[i, j];
                }
            }
            return new Matrix_Diag(resOffsets, res);
        }
        public static Matrix_Diag Parse(MatrixSymm_Sparse mat)
        {
            int size = mat.size;
            int[] ig = mat.ig;
            int[] jg = mat.jg;
            double[] gg = mat.gg;
            List<int> offsets = new List<int>();
            List<double[]> res = new List<double[]>();
            int diagIndex;
            int offset;
            for (int i = 0; i < size; i++)
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    if (gg[j] != 0)
                    {
                        offset = jg[i] - i;
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset) ;
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][i] = gg[j];
                        offset = -offset;
                        diagIndex = offsets.IndexOf(offset);
                        if (diagIndex == -1)
                        {
                            while ((++diagIndex) < offsets.Count && offsets[diagIndex] < offset) ;
                            offsets.Insert(diagIndex, offset);
                            res.Insert(diagIndex, new double[size]);
                        }
                        res[diagIndex][jg[i]] = gg[j];
                    }
                }
            int diagCount = offsets.Count;
            double[,] result = new double[diagCount, size];
            for (int i = 0; i < diagCount; i++)
                for (int j = 0; j < size; j++)
                    result[i, j] = res[i][j];
            return new Matrix_Diag(offsets.ToArray(), result);
        }
        public static Vector operator *(Matrix_Diag mat, Vector vec)
        {
            int size = mat.size;
            int i, j, k;
            int diagCount = mat.diagCount;
            double[,] values = mat.values;
            int[] offsets = mat.offsets;
            double[] x = vec.values;

            Vector res = new Vector(size);

            for (i = 0; i < diagCount; i++)
                for (j = Math.Max(-offsets[i], 0), k = Math.Max(offsets[i], 0); j < size && k < size; j++, k++)
                    res[j] += x[k] * values[i, j];
            return res;
        }
    }
    public class MatrixSymm_Diag : MatrixSymm
    {
        public double[,] values { get; private set; }
        public int[] offsets { get; private set; }
        public int diagCount { get { return offsets.Length * 2 - (offsets.Contains(0) ? 1 : 0); } }
        public override int size { get { return values.GetLength(1); } protected set { } }
        public MatrixSymm_Diag(int size, int[] offsets)
        {
            this.offsets = offsets;
            values = new double[offsets.Length, size];
        }
        public MatrixSymm_Diag(int[] offsets, double[,] values)
        {
            this.offsets = offsets;
            this.values = values;
        }
        public MatrixSymm_Diag(MatrixSymm_Diag mat)
        {
            int size = mat.size;
            int diagCount = mat.diagCount;
            values = new double[diagCount, size];
            offsets = new int[diagCount];
            double[,] origValues = mat.values;
            int[] origOffsets = mat.offsets;
            for (int i = 0; i < diagCount; i++)
            {
                offsets[i] = origOffsets[i];
                for (int j = 0; j < size; j++)
                    values[i, j] = origValues[i, j];
            }
        }
    }
}
