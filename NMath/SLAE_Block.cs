using System;
using System.Collections.Generic;
using System.Text;

namespace NMath
{
    public static partial class SLAE
    {
        private class MemBlock_Block
        {
            public double[][,] blocks;
            public double[] r;
            public int blocksCount, blockSize;
            public MemBlock_Block(int blocksCount, int blockSize)
            {
                blocks = new double[blocksCount][,];
                for (int i = 0; i < blocksCount; i++)
                    blocks[i] = new double[blockSize, blockSize];
                r = new double[blockSize];
                this.blocksCount = blocksCount;
                this.blockSize = blockSize;
            }
        }
        private static MemBlock_Block Mem_Block = null;
        public static void ClearMem_Block()
        {
            Mem_Block = null;
        }
        public static void Solve_Block(Matrix matrix, Vector vec, Vector start, double weight = 0.75, int blockSize = 2)
        {
            if (matrix is Matrix_Diag)
                Solve_Block(matrix as Matrix_Diag, vec, start, weight, blockSize);
            else
                throw new NotImplementedException("This solving method does not support this matrix type.");
        }
        public static void Solve_Block(Matrix_Diag matrix, Vector vec, Vector start, double weight = 0.75, int blockSize = 2)
        {
            if (matrix.size % blockSize != 0)
                throw new ArgumentException("Matrix should be equally divided by block size.");

            int i, j, k, b, blockOffset;

            int blocksCount = matrix.size / blockSize;

            if (Mem_Block == null || Mem_Block.blockSize != blockSize || Mem_Block.blocksCount != blocksCount)
                Mem_Block = new MemBlock_Block(blocksCount, blockSize);

            double[][,] diagBlocks = Mem_Block.blocks;
            for (b = 0; b < blocksCount; b++)
            {
                blockOffset = b * blockSize;
                for (i = 0; i < blockSize; i++)
                    for (j = 0; j < blockSize; j++)
                        if ((k = Array.IndexOf(matrix.offsets, j - i)) != -1)
                        diagBlocks[b][i, j] = matrix.values[k, blockOffset + i];
                for (i = 1; i < blockSize; i++)
                {
                    diagBlocks[b][i, 0] = diagBlocks[b][i, 0] / diagBlocks[b][0, 0];
                    for (k = 0; k < i; k++)
                        diagBlocks[b][i, i] -= diagBlocks[b][i, k] * diagBlocks[b][k, i];
                    for (j = i + 1; j < blockSize; j++)
                    {
                        for (k = 0; k < i; k++)
                        {
                            diagBlocks[b][i, j] -= diagBlocks[b][i, k] * diagBlocks[b][k, j];
                            diagBlocks[b][j, i] -= diagBlocks[b][j, k] * diagBlocks[b][k, i];
                        }
                        diagBlocks[b][j, i] /= diagBlocks[b][i, i];
                    }
                }
            }

            double norm = 0;
            for (i = 0; i < matrix.size; i++)
                norm += vec[i] * vec[i];

            double[] R = Mem_Block.r;
            int step;
            double curNorm;
            int verticalOffset;
            int diagIndex;
            for (step = 0; step < maxSteps; step++)
            {
                curNorm = 0;
                for (i = 0; i < blocksCount; i++)
                {
                    verticalOffset = i * blockSize;
                    for (j = 0; j < blockSize; j++)
                        R[j] = vec[verticalOffset + j];
                    for (diagIndex = 0; diagIndex < matrix.diagCount; diagIndex++)
                        for (j = verticalOffset; j < verticalOffset + blockSize; j++)
                        {
                            if (matrix.values[diagIndex, j] == 0)
                                continue;
                            R[j - verticalOffset] -= matrix.values[diagIndex, j] * start[j + matrix.offsets[diagIndex]];
                        }
                    for (j = 0; j < blockSize; j++)
                    {
                        curNorm += R[j] * R[j];
                        R[j] *= weight;
                    }

                    // straight
                    for (j = 1; j < blockSize; j++)
                        for (k = 0; k < j; k++)
                            R[j] -= diagBlocks[i][j, k] * R[k];

                    // reverse
                    for (j = blockSize - 1; j >= 0; j--)
                    {
                        for (k = j + 1; k < blockSize; k++)
                            R[j] -= diagBlocks[i][j, k] * R[k];
                        R[j] /= diagBlocks[i][j, j];
                    }

                    // X(k+1) = X(k) + Y
                    for (j = 0; j < blockSize; j++)
                        start[verticalOffset + j] += R[j];
                }
                if (curNorm / norm < epsilon * epsilon)
                    break;
            }
        }
    }
}
