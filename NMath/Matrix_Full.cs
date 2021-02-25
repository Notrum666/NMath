using System;
using System.Collections.Generic;
using System.Text;

namespace NMath
{
    public class Matrix_Full : Matrix
    {
        public double[,] values { get; private set; }
        public double this[int i, int j] { get => values[i, j]; set { values[i, j] = value; } }
        public override int size { get; protected set; }
        public Matrix_Full(int size)
        {
            values = new double[size, size];
            this.size = size;
        }
        public Matrix_Full(Matrix_Full mat)
        {
            size = mat.size;
            values = new double[size, size];
            double[,] val = mat.values;
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    values[i, j] = val[i, j];
        }
        public static Vector operator *(Matrix_Full mat, Vector vec)
        {
            if (mat.size != vec.size)
                throw new ArgumentException("Matrix and vector must be same size.");
            int size = vec.size;
            Vector result = new Vector(size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    result[i] += mat[i, j] * vec[j];
            return result;
        }
        public static Matrix_Full Parse(Matrix mat)
        {
            if (mat is Matrix_Full)
                return new Matrix_Full(mat as Matrix_Full);
            if (mat is Matrix_Diag)
                return Parse(mat as Matrix_Diag);
            if (mat is Matrix_Sparse)
                return Parse(mat as Matrix_Sparse);
            throw new NotImplementedException("This parsing method is not implemented yet.");
        }
        public static Matrix_Full Parse(MatrixSymm mat)
        {
            if (mat is MatrixSymm_Sparse)
                return Parse(mat as MatrixSymm_Sparse);
            if (mat is MatrixSymm_Diag)
                return Parse(mat as MatrixSymm_Diag);
            if (mat is MatrixSymm_Full)
                return Parse(mat as MatrixSymm_Full);
            throw new NotImplementedException("This parsing method is not implemented yet.");
        }
        public static Matrix_Full Parse(MatrixSymm_Full mat)
        {
            int size = mat.size;
            Matrix_Full result = new Matrix_Full(size);
            double[,] res = result.values;
            double[][] values = mat.values;
            for (int i = 0; i < size; i++)
            {
                res[i, i] = values[i][i];
                for (int j = 0; j < i; j++)
                {
                    res[i, j] = values[i][j];
                    res[j, i] = values[i][j];
                }
            }
            return result;
        }
        public static Matrix_Full Parse(Matrix_Diag mat)
        {
            int size = mat.size;
            Matrix_Full result = new Matrix_Full(size);
            double[,] res = result.values;
            double[,] origValues = mat.values;

            for (int i = 0; i < mat.diagCount; i++)
                for (int j = Math.Max(0, -mat.offsets[i]); j < size - Math.Max(0, mat.offsets[i]); j++)
                    res[j, mat.offsets[i] + j] = origValues[i, j];

            return result;
        }
        public static Matrix_Full Parse(Matrix_Sparse mat)
        {
            int size = mat.size;
            Matrix_Full result = new Matrix_Full(size);
            double[,] res = result.values;

            int[] ig = mat.ig;
            int[] jg = mat.jg;
            double[] ggu = mat.ggu;
            double[] ggl = mat.ggl;
            double[] di = mat.di;

            for (int i = 0; i < size; i++)
            {
                res[i, i] = di[i];
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    res[i, jg[j]] = ggl[j];
                    res[jg[j], i] = ggu[j];
                }
            }

            return result;
        }
        public static Matrix_Full Parse(MatrixSymm_Sparse mat)
        {
            int size = mat.size;
            Matrix_Full result = new Matrix_Full(size);
            double[,] res = result.values;

            int[] ig = mat.ig;
            int[] jg = mat.jg;
            double[] gg = mat.gg;
            double[] di = mat.di;

            for (int i = 0; i < size; i++)
            {
                res[i, i] = di[i];
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    res[i, jg[j]] = gg[j];
                    res[jg[j], i] = gg[j];
                }
            }

            return result;
        }
        public static Matrix_Full Parse(MatrixSymm_Diag mat)
        {
            int size = mat.size;
            Matrix_Full result = new Matrix_Full(size);
            int diagCount = mat.diagCount;
            double[,] res = result.values;
            double[,] values = mat.values;
            int[] offsets = mat.offsets;

            for (int i = 0; i < diagCount; i++)
                for (int j = -offsets[i]; j < size; j++)
                {
                    res[j, j + mat.offsets[i]] = values[i, j];
                    res[j + offsets[i], j] = values[i, j];
                }
            return result;
        }
        public override string ToString()
        {
            string[] result = new string[size];
            for (int i = 0; i < size; i++)
            {
                result[i] += "|";
                for (int j = 0; j < size; j++)
                {
                    result[i] += values[i, j].ToString();
                    if (j != size - 1)
                        result[i] += " ";
                }
                result[i] += "|";
            }
            return string.Join("\n", result);
        }
        public string ToString(string format)
        {
            string[] result = new string[size];
            for (int i = 0; i < size; i++)
            {
                result[i] += "|";
                for (int j = 0; j < size; j++)
                {
                    result[i] += values[i, j].ToString(format);
                    if (j != size - 1)
                        result[i] += " ";
                }
                result[i] += "|";
            }
            return string.Join("\n", result);
        }
    }
    public class MatrixSymm_Full : MatrixSymm
    {
        public double[][] values { get; private set; }
        public double this[int i, int j] { get { return i > j ? values[i][j] : values[j][i]; } set { if (i > j) values[i][j] = value; else values[j][i] = value; } }
        public override int size { get; protected set; }
        public MatrixSymm_Full(int size)
        {
            values = new double[size][];
            for (int i = 0; i < size; i++)
                values[i] = new double[i + 1];
            this.size = size;
        }
        public MatrixSymm_Full(MatrixSymm_Full mat)
        {
            size = mat.size;
            values = new double[size][];
            double[][] val = mat.values;
            for (int i = 0; i < size; i++)
            {
                values[i] = new double[i + 1];
                for (int j = 0; j <= i; j++)
                    values[i][j] = val[i][j];
            }
        }
        public static Vector operator *(MatrixSymm_Full mat, Vector vec)
        {
            if (mat.size != vec.size)
                throw new ArgumentException("Matrix and vector must be same size.");
            int size = vec.size;
            double[][] values = mat.values;
            double[] v = vec.values;
            Vector result = new Vector(size);
            double[] res = result.values;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    res[i] += values[i][j] * v[j];
                    res[j] += values[i][j] * v[i];
                }
            }
            return result;
        }
        public static bool TryParse(Matrix mat, out MatrixSymm_Full result)
        {
            if (mat is Matrix_Full)
                return TryParse(mat as Matrix_Full, out result);
            if (mat is Matrix_Diag)
                return TryParse(mat as Matrix_Diag, out result);
            if (mat is Matrix_Sparse)
                return TryParse(mat as Matrix_Sparse, out result);
            throw new NotImplementedException("This parsing method is not implemented yet.");
        }
        public static MatrixSymm_Full Parse(MatrixSymm mat)
        {
            if (mat is MatrixSymm_Full)
                return new MatrixSymm_Full(mat as MatrixSymm_Full);
            if (mat is MatrixSymm_Diag)
                return Parse(mat as MatrixSymm_Diag);
            if (mat is MatrixSymm_Sparse)
                return Parse(mat as MatrixSymm_Sparse);
            throw new NotImplementedException("This parsing method is not implemented yet.");
        }
        public static bool TryParse(Matrix_Full mat, out MatrixSymm_Full result)
        {
            int size = mat.size;
            result = new MatrixSymm_Full(size);
            double[][] res = result.values;
            double[,] val = mat.values;
            for (int i = 0; i < size; i++)
                for (int j = 0; j <= i; j++)
                {
                    if (val[i, j] != val[j, i])
                    {
                        result = null;
                        return false;
                    }
                    res[i][j] = val[i, j];
                }
            return true;
        }
        public static bool TryParse(Matrix_Diag mat, out MatrixSymm_Full result)
        {
            int size = mat.size;
            result = new MatrixSymm_Full(size);
            double[][] res = result.values;
            double[,] val = mat.values;
            int[] offsets = mat.offsets;
            int symmDiagIndex;
            int i;
            for (i = 0; offsets[i] <= 0; i++)
            {
                symmDiagIndex = Array.IndexOf(offsets, -offsets[i]);
                if (symmDiagIndex == -1)
                {
                    result = null;
                    return false;
                }
                for (int j = -offsets[i]; j < size; j++)
                {
                    if (val[i, j] != val[symmDiagIndex, j + offsets[i]])
                    {
                        result = null;
                        return false;
                    }
                    res[j][offsets[i] + j] = val[i, j];
                }
            }
            if (i * 2 - 1 != offsets.Length)
            {
                result = null;
                return false;
            }
            return true;
        }
        public static bool TryParse(Matrix_Sparse mat, out MatrixSymm_Full result)
        {
            int size = mat.size;
            result = new MatrixSymm_Full(size);
            double[][] res = result.values;
            int[] ig = mat.ig;
            int[] jg = mat.jg;
            double[] ggu = mat.ggu;
            double[] ggl = mat.ggl;
            double[] di = mat.di;
            for (int i = 0; i < size; i++)
            {
                res[i][i] = di[i];
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    if (ggu[j] != ggl[j])
                    {
                        result = null;
                        return false;
                    }
                    res[i][jg[j]] = ggl[j];
                }
            }
            return true;
        }
        public static MatrixSymm_Full Parse(MatrixSymm_Diag mat)
        {
            int size = mat.size;
            MatrixSymm_Full result = new MatrixSymm_Full(size);
            double[][] res = result.values;
            int diagCount = mat.diagCount;
            int[] offsets = mat.offsets;
            double[,] values = mat.values;
            for (int i = 0; i < diagCount; i++)
                for (int j = -offsets[i]; j < size; j++)
                    res[j][j + offsets[i]] = values[i, j];
            return result;
        }
        public static MatrixSymm_Full Parse(MatrixSymm_Sparse mat)
        {
            int size = mat.size;
            MatrixSymm_Full result = new MatrixSymm_Full(size);
            double[][] res = result.values;
            int[] ig = mat.ig;
            int[] jg = mat.jg;
            double[] gg = mat.gg;
            double[] di = mat.di;
            for (int i = 0; i < size; i++)
            {
                res[i][i] = di[i];
                for (int j = ig[i]; j < ig[i + 1]; j++)
                    res[i][jg[j]] = gg[j];
            }
            return result;
        }
        public override string ToString()
        {
            string[] result = new string[size];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    result[i] += values[i][j].ToString();
                    result[j] += values[i][j].ToString();
                    if (i != size - 1)
                        result[j] += " ";
                }
                result[i] += values[i][i];
                if (i != size - 1)
                    result[i] += " ";
            }
            return "|" + string.Join("|\n|", result) + "|";
        }
        public string ToString(string format)
        {
            string[] result = new string[size];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    result[i] += values[i][j].ToString(format);
                    result[j] += values[i][j].ToString(format);
                    if (i != size - 1)
                        result[j] += " ";
                }
                result[i] += values[i][i];
                if (i != size - 1)
                    result[i] += " ";
            }
            return "|" + string.Join("|\n|", result) + "|";
        }
    }
}
