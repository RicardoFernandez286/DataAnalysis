// anyEq.c
// anyEq - Fast check if two arrays have any common element
// R = anyEq(X, Y)
// INPUT:
//   X, Y: Arrays of any size. Complex or sparse array are rejected.
//         Types: DOUBLE, SINGLE, (U)INT8/16/32/64, CHAR, LOGICAL.
// OUTPUT:
//   R:    TRUE is replied if any element of X occurs in Y, FALSE otherwise.
//         NaNs are treated as equal.
//
// NOTES:
// - This is equivalent to:
//     R = any(X(:) == Y(1)) || any(X(:) == Y(2)) || ...
// - This MEX version is faster than the Matlab method, because the creation of
//   the intermediate logical array is avoided:
//     Worst case (no matching element):  25% to 60% faster
//     Best case (first element matches): 99.99% faster for 1e6 elements
//     (Matlab 2011b/64, MSVC 2008).
// - For small LOGICAL arrays (< 5'000 elements) anyEq is up to 65% slower than
//   any(X) or ~all(X). For larger arrays the speedup depends on the position of
//   the first match and can reach a factor of 7.
// - X and Y are ordered automatically such, that the elements of the smaller
//   array are searched in the larger one.
//
// EXAMPLES:
//   anyEq(0:0.1:1, 0.3)   % FALSE: Effect of limited precision
//   anyEq(1:4, [5,2])     % TRUE:  2 is found in 1:4
//
// COMPILATION:
// MEX files m,ust be compiled before using:
//   mex -O anyEq.c
// Consider C99 comments on Linux:
//   mex -O CFLAGS="\$CFLAGS -std=c99" anyEq.c
// Pre-compiled Mex: http://www.n-simon.de/mex
//
// TEST: Run uTest_anyEq to check validity and speed of the Mex function.
//
// NOTE: This function can be used as template for other operators easily, e.g.
//   ~=, >, <, >=, <=. Feel free to mail me if you need support.
//
// Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
// Assumed Compatibility: higher Matlab versions, Mac, Linux
// Author: Jan Simon, Heidelberg, (C) 2006-2013 matlab.THISYEAR(a)nMINUSsimon.de
//
// See also ANY, ALL, anyExceed.

/*
% $JRev: R-E V:042 Sum:7O5ftPaVsFIe Date:11-Sep-2013 00:48:24 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_anyEq $
% $File: Tools\Mex\Source\anyEq.c $
% History:
% 001: 04-Mar-2010 01:45, First version.
% 041: 08-Sep-2013 23:31, LOGICAL input accepted.
*/

#include "mex.h"

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_HEAD "*** anyEq[mex]: "
#define ERR_ID   "JSimon:anyEq:"
#define ERROR(id,msg) mexErrMsgIdAndTxt(ERR_ID id, ERR_HEAD msg);

// LCC v2.4 and v3.8 have problems with int64_T, which I do not understand.
// Using the equivalent LONG LONG instead works well, but slower than DOUBLEs
// inspite of the additional tests for NaN!
#ifdef __LCC__
#define int64_T long long
#endif

// 64 bit constants to compare logicals in blocks of 8 elements:
#if defined(__BORLANDC__)     // No "LL" for BCC 5.5
#  define ALLZERO_int64 (int64_T) 0
#  define ALLONE_int64  (int64_T) 0x0101010101010101L
#else
#  define ALLZERO_int64 0LL
#  define ALLONE_int64  0x0101010101010101LL
#endif

// No mxIsNaN for float (_isnanf for MSVC, isnanf for LCC, not defined in
// Open Watcom 1.8 and BCC 5.5):
#define ISNANF(x) (((*(int32_T *)&(x) & 0x7f800000L) == 0x7f800000L) && \
                   ((*(int32_T *)&(x) & 0x007fffffL) != 0x00000000L))

// Prototypes:
mxLogical Eq_Double(double *X, mwSize nX, double  *Y, mwSize nY);
mxLogical Eq_Single(float  *X, mwSize nX, float   *Y, mwSize nY);
mxLogical Eq_8Byte(int64_T *X, mwSize nX, int64_T *Y, mwSize nY);
mxLogical Eq_4Byte(int32_T *X, mwSize nX, int32_T *Y, mwSize nY);
mxLogical Eq_2Byte(int16_T *X, mwSize nX, int16_T *Y, mwSize nY);
mxLogical Eq_1Byte(int8_T  *X, mwSize nX, int8_T  *Y, mwSize nY);
mxLogical Eq_Logical(mxLogical *X, mwSize nX, mxLogical *Y, mwSize nY);
mxLogical Eq_Logical_true(mxLogical  *X, mwSize nX);
mxLogical Eq_Logical_false(mxLogical *X, mwSize nX);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxLogical Eq;
  mwSize    nX, nY, nS;
  void      *X, *Y;
  mxClassID XClass;
  
  // Proper number of arguments:
  if (nrhs != 2) {
     ERROR("BadNInput", "2 inputs required.");
  }
  if (nlhs > 1) {
     ERROR("BadNOutput", "1 output allowed.");
  }
  
  // Both arrays must have the same type:
  XClass = mxGetClassID(prhs[0]);
  if (XClass != mxGetClassID(prhs[1])) {
     ERROR("NotMatchingClass", "Both inputs must have the same class.");
  }
  if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])) {
     // If anybody needs this, please mail me!
     ERROR("Complex", "Complex array are not handled.");
  }
  
  // Let X be the longer array:
  nX = mxGetNumberOfElements(prhs[0]);
  nY = mxGetNumberOfElements(prhs[1]);
  if (nX >= nY) {
     X = mxGetData(prhs[0]);
     Y = mxGetData(prhs[1]);
  } else {
     nS = nY;      // Swap nX and nY
     nY = nX;
     nX = nS;
     X  = mxGetData(prhs[1]);
     Y  = mxGetData(prhs[0]);
  }
  
  // Call subroutines for the different classes:
  // (Multi-threading: Split X and nX to parts here!)
  switch (XClass) {
     case mxDOUBLE_CLASS:   // Consider NaNs:
        Eq = Eq_Double((double *) X, nX, (double *) Y, nY);
        break;
     case mxSINGLE_CLASS:   // Consider NaNs:
        Eq = Eq_Single((float *) X, nX, (float *) Y, nY);
        break;
     case mxLOGICAL_CLASS:
        Eq = Eq_Logical((mxLogical *) X, nX, (mxLogical *) Y, nY);
        break;
     case mxINT64_CLASS:
     case mxUINT64_CLASS:
        Eq = Eq_8Byte((int64_T *) X, nX, (int64_T *) Y, nY);
        break;
     case mxINT32_CLASS:
     case mxUINT32_CLASS:
        Eq = Eq_4Byte((int32_T *) X, nX, (int32_T *) Y, nY);
        break;
     case mxINT16_CLASS:
     case mxUINT16_CLASS:
     case mxCHAR_CLASS:
        Eq = Eq_2Byte((int16_T *) X, nX, (int16_T *) Y, nY);
        break;
     case mxINT8_CLASS:
     case mxUINT8_CLASS:
        Eq = Eq_1Byte((int8_T *) X, nX, (int8_T *) Y, nY);
        break;
               
     default: ERROR("BadClass", "Inputs must be numeric, CHAR or LOGICAL.");
  }
  
  // Create the output:
  plhs[0] = mxCreateLogicalScalar(Eq);
  
  return;
}

// Subroutines =================================================================
mxLogical Eq_Double(double *X, mwSize nX, double *Y, mwSize nY)
{
  // Works for DOUBLE considering NaNs are equal to each other.
  // The comparison with NaNs should reply FALSE ever, but for some compilers
  // it doesn't.
  
  double *Xp, *Xf, *Yf;
  
  Xf = X + nX;
  for (Yf = Y + nY; Y < Yf; Y++) {  // Loop through shorter input array
     Xp = X - 1;                    // Start before X due to ++X increment
     if (mxIsNaN(*Y)) {             // If Y is NaN, the == operator fails
        while (++Xp < Xf) {         // Loop through longer input array
           if (mxIsNaN(*Xp)) {      // Check if X element is NaN also
              return (true);
           }
        }
        
     } else {                       // Searched value is not NaN:
        while (++Xp < Xf) {         // Loop through longer input array
           if (*Xp == *Y) {         // Check for equality
              if (!mxIsNaN(*Xp)) {  // Success only if X is not NaN
                 return (true);
              }
           }
        }
     }
  }

  return (false);
}

// =============================================================================
mxLogical Eq_Single(float *X, mwSize nX, float *Y, mwSize nY)
{
  // Works for SINGLE considering NaNs are equal to each other.
  // This equals Eq_Double except for the pointer types.
  
  float *Xp, *Xf, *Yf;
  
  Xf = X + nX;
  for (Yf = Y + nY; Y < Yf; Y++) {
     Xp = X - 1;
     if (ISNANF(*Y)) {  // Search for a NaN:
        while (++Xp < Xf) {
           if (ISNANF(*Xp)) {
              return (true);
           }
        }
        
     } else {           // Searched value is not NaN:
        while (++Xp < Xf) {
           if (*Xp == *Y) {
              if (!ISNANF(*Xp)) {  // Needed for some compilers!
                 return (true);
              }
           }
        }
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_8Byte(int64_T *X, mwSize nX, int64_T *Y, mwSize nY)
{
  // Works for INT64, UINT64.
  // Eq_Double replies the same result except for NaNs.
  
  int64_T *Xp, *Xf, *Yf;
  
  Xf = X + nX;
  for (Yf = Y + nY; Y < Yf; Y++) {  // Loop through shorter input array
     Xp = X;                        // Start loop through longer input array
     while (Xp < Xf) {              // Loop until end of X
        if (*Xp++ == *Y) {          // Check for equality, no NaN in INT*
           return (true);
        }
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_4Byte(int32_T *X, mwSize nX, int32_T *Y, mwSize nY)
{
  // Works for INT32, UINT32.
  // Eq_Single replies the same result except for NaNs.
  
  int32_T *Xp, *Xf, *Yf;
  
  Xf = X + nX;
  for (Yf = Y + nY; Y < Yf; Y++) {
     Xp = X;
     while (Xp < Xf) {
        if (*Xp++ == *Y) {
           return (true);
        }
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_2Byte(int16_T *X, mwSize nX, int16_T *Y, mwSize nY)
{
  // Works for INT16 and UINT16 and CHAR.
  
  int16_T *Xp, *Xf, *Yf;
  
  Xf = X + nX;
  for (Yf = Y + nY; Y < Yf; Y++) {
     Xp = X;
     while (Xp < Xf) {
        if (*Xp++ == *Y) {
           return (true);
        }
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_1Byte(int8_T *X, mwSize nX, int8_T *Y, mwSize nY)
{
  // Works for INT8 and UINT8.
  int8_T *Xp, *Xf, *Yf;
  
  Xf = X + nX;
  for (Yf = Y + nY; Y < Yf; Y++) {
     Xp = X;
     while (Xp < Xf) {
        if (*Xp++ == *Y) {
           return (true);
        }
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_Logical(mxLogical *X, mwSize nX, mxLogical *Y, mwSize nY)
{
  // Special treatment of logicals.
    
  while (nY-- != 0) {
     if (*Y++ != 0) {
        if (Eq_Logical_true(X, nX)) {
           return (true);
        }
     } else {
        if (Eq_Logical_false(X, nX)) {
           return (true);
        }
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_Logical_true(mxLogical *X, mwSize nX)
{
  // Check if any element is TRUE.
   
  mxLogical *Xf = X + nX;
  int64_T   *X64, *X64f;

  // Skip initial bytes if X is not aligned to 8 bytes - not expected:
  while (X < Xf && (((size_t) X & 7) != 0)) {
     if (*X++ != 0) {
        return (true);
     }
  }
  
  // Process aligned data block as uint64 elements:
  X64  = (int64_T *) X;
  X64f = (int64_T *) ((size_t) Xf - ((size_t) Xf & 7));
  // Perhaps this is safer:
  // #define PTR_OFF (uintptr_t)(const void *)
  // X64f = (int64_T *) (PTR_OFF Xf - (PTR_OFF Xf & 7));
  while (X64 < X64f) {
     if (*X64++ != ALLZERO_int64) {
        return (true);
     }
  }
  
  // Trailing bytes:
  X = (mxLogical *) X64;
  while (X < Xf) {
     if (*X++ != 0) {
        return (true);
     }
  }
  
  return (false);
}

// =============================================================================
mxLogical Eq_Logical_false(mxLogical *X, mwSize nX)
{
  // Check if any element is TRUE.
   
  mxLogical *Xf = X + nX;
  int64_T   *X64, *X64f;

  // Leading bytes if X is not aligned to 8 bytes - not expected:
  while (X < Xf && (((size_t) X & 7) != 0)) {
     if (*X++ == 0) {
        return (true);
     }
  }
  
  // Process aligned data block as uint64 elements:
  X64  = (int64_T *) X;
  X64f = (int64_T *) ((size_t) Xf - ((size_t) Xf & 7));
  // Perhaps this is safer:
  // #define PTR_OFF (uintptr_t)(const void *)
  // X64f = (int64_T *) (PTR_OFF Xf - (PTR_OFF Xf & 7));
  while (X64 < X64f) {
     if (*X64++ != ALLONE_int64) {
        return (true);
     }
  }
  
  // Trailing bytes:
  X = (mxLogical *) X64;
  while (X < Xf) {
     if (*X++ == 0) {
        return (true);
     }
  }
  
  return (false);
}
