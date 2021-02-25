using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using NMath;

namespace NMath_Tests
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void SLAE_Gauss_Full_9x9()
        {
            Matrix_Full mat = new Matrix_Full(9);
            for (int i = 0; i < 9; i++)
                mat[i, i] = 1;
            mat[4, 1] = -4.0;
            mat[4, 3] = -6.0;
            mat[4, 4] = 18.0;
            mat[4, 5] = -3.0;
            mat[4, 7] = -4.0;
            Vector vec = new Vector(1, 1, 1, 1, 1, 1, 1, 1, 1);
            SLAE.Solve_Gauss(mat, vec);

            Assert.AreEqual(new Vector(1, 1, 1, 1, 1, 1, 1, 1, 1), vec);
        }
        //[TestMethod]
        //public void Matrix_Full_Parse_ToDiag_Empty()
        //{
        //    Matrix_Full mat = new Matrix_Full(10);
        //
        //}
    }
}
