using System;
using System.Collections.Generic;
using System.Text;

namespace NMath
{
    public static partial class SLAE
    {
        public static double epsilon { get; set; } = 0.0000001;
        public static int maxSteps { get; set; } = 100000;
        public static void Clear_Mem()
        {
            ClearMem_Jacobi();
            ClearMem_Block();
            ClearMem_CGM();
            ClearMem_CGM_LLt();
        }
    }
}
