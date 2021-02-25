using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace NMath
{
    public static partial class SLAE
    {
        public static void Solve_Gauss(Matrix matrix, Vector vec)
        {
            if (matrix is Matrix_Full)
                Solve_Gauss(matrix as Matrix_Full, vec);
            else
                throw new NotImplementedException("This solving method does not support this matrix type.");
        }
        public static void Solve_Gauss(Matrix_Full matrix, Vector vec)
        {
            double[,] mat = matrix.values;
            int size = matrix.size;

            int i, j, k;
            double max;
            int maxIndex;
            double tmp;

            File.WriteAllText("tmp.txt", matrix.ToString("f4") + "\n\n");
            for (i = 0; i < size; i++)
            {
                max = Math.Abs(mat[i, i]);
                maxIndex = i;
                for (j = i + 1; j < size; j++)
                {
                    tmp = Math.Abs(mat[j, i]);
                    if (tmp > max)
                    {
                        max = tmp;
                        maxIndex = j;
                    }
                }
                if (maxIndex != i)
                {
                    for (j = i; j < size; j++)
                    {
                        tmp = mat[i, j];
                        mat[i, j] = mat[maxIndex, j];
                        mat[maxIndex, j] = tmp;
                    }
                    tmp = vec[i];
                    vec[i] = vec[maxIndex];
                    vec[maxIndex] = tmp;
                }
                File.AppendAllText("tmp.txt", matrix.ToString("f4") + "\n\n");

                for (j = i + 1; j < size; j++)
                {
                    if (mat[j, i] == 0)
                        continue;
                    tmp = mat[j, i] / mat[i, i];
                    for (k = i; k < size; k++)
                        mat[j, k] -= tmp * mat[i, k];
                    vec[j] -= tmp * vec[i];
                }
                File.AppendAllText("tmp.txt", matrix.ToString("f4") + "\n\n");
            }
            for (i = size - 1; i >= 0; i--)
            {
                for (j = i + 1; j < size; j++)
                    vec[i] -= mat[i, j] * vec[j];
                vec[i] /= mat[i, i];
            }
            File.AppendAllText("tmp.txt", matrix.ToString("f4") + "\n\n");
        }
    }
}
