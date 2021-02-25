using System;
using System.Collections.Generic;
using System.Text;

namespace NMath
{
    public static partial class SLAE
    {
        public static void Solve_Seidel(Matrix matrix, Vector vec, Vector start, double weight = 0.75)
        {
            if (matrix is Matrix_Diag)
                Solve_Seidel(matrix as Matrix_Diag, vec, start, weight);
            else
                throw new NotImplementedException("This solving method does not support this matrix type.");
        }
        public static void Solve_Seidel(Matrix_Diag matrix, Vector vec, Vector start, double weight = 0.75)
        {
            int size = matrix.size;
            int diagCount = matrix.diagCount;
            int mainDiagIndex = Array.IndexOf(matrix.offsets, 0);
            double[,] values = matrix.values;
            int[] offsets = matrix.offsets;
            double[] v = vec.values;
            double[] x = start.values;
            int i, j;
            double collector;
            int step;
            double norm = 0;
            double curNorm;
            int offset;
            for (i = 0; i < size; i++)
                norm += v[i] * v[i];

            for (step = 0; step < maxSteps; step++)
            {
                curNorm = 0;
                for (i = 0; i < size; i++)
                {
                    collector = 0;
                    for (j = 0; j < diagCount; j++)
                    {
                        offset = i + offsets[j];
                        if (offset < 0)
                            continue;
                        if (offset >= size)
                            break;
                        collector += values[j, i] * x[offset];
                    }
                    x[i] += weight * (v[i] - collector) / values[mainDiagIndex, i];
                    curNorm += (v[i] - collector) * (v[i] - collector);
                }
                if (curNorm / norm < epsilon * epsilon)
                    break;
            }
        }
    }
}
