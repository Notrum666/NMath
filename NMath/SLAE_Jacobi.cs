using System;
using System.Collections.Generic;
using System.Text;

namespace NMath
{
    public static partial class SLAE
    {
        private class MemBlock_Jacobi
        {
            public double[] collector;
            public int size;
            public MemBlock_Jacobi(int size)
            {
                collector = new double[size];
                this.size = size;
            }
        }
        private static MemBlock_Jacobi Mem_Jacobi = null;
        public static void ClearMem_Jacobi()
        {
            Mem_Jacobi = null;
        }
        public static void Solve_Jacobi(Matrix matrix, Vector vec, Vector start, double weight = 0.75)
        {
            if (matrix is Matrix_Diag)
                Solve_Jacobi(matrix as Matrix_Diag, vec, start, weight);
            else
                throw new NotImplementedException("This solving method does not support this matrix type.");
        }
        public static void Solve_Jacobi(Matrix_Diag matrix, Vector vec, Vector start, double weight = 0.75)
        {
            int size = matrix.size;
            if (Mem_Jacobi == null || Mem_Jacobi.size != size)
                Mem_Jacobi = new MemBlock_Jacobi(size);

            int diagCount = matrix.diagCount;
            int mainDiagIndex = Array.IndexOf(matrix.offsets, 0);
            double[,] values = matrix.values;
            int[] offsets = matrix.offsets;
            double[] v = vec.values;
            double[] x = start.values;
            int i, j, k;
            double[] collector = Mem_Jacobi.collector;
            int step;
            double norm = 0;
            double curNorm;

            for (i = 0; i < size; i++)
            {
                norm += v[i] * v[i];
                collector[i] = 0;
            }

            for (step = 0; step < maxSteps; step++)
            {
                curNorm = 0;
                for (i = 0; i < diagCount; i++)
                    for (j = Math.Max(-offsets[i], 0), k = Math.Max(offsets[i], 0); j < size && k < size; j++, k++)
                        collector[j] += x[k] * values[i, j];
                for (i = 0; i < size; i++)
                {
                    curNorm += (v[i] - collector[i]) * (v[i] - collector[i]);
                    x[i] += weight * (v[i] - collector[i]) / values[mainDiagIndex, i];
                    collector[i] = 0;
                }
                if (curNorm / norm < epsilon * epsilon)
                    break;
            }
        }
    }
}
