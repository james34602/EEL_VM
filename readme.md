EEL2 is a language that powers [REAPER](https://www.reaper.fm/)'s [JSFX](https://www.reaper.fm/sdk/js), LiveProgVST, and is also used by [JamesDSP](https://github.com/james34602/JamesDSPManager). This document is purely for the core EEL2 functionality.

* Basic Language Attributes
* Operator reference
* Simple math functions
* Loops
* User defined functions and namespace pseudo-objects
* Advanced Mathematics, Linear algebra and Signal Processing Functions
* Strings
* String functions
* Misc utility

***

**Basic Language Attributes**

The core of EEL2 has many similarities to C but is distictly different. Some notable qualities of this language are:

- Variables do not need to be declared, are by default global, and are all single-precision floating point.

- Parentheses "(" and ")" can be used to clarify precidence, contain parameters for functions, and collect multiple statements into a single statement.

- A semicolon ";" is used to separate statements from eachother (including within parentheses).

- A virtual local address space of about 8 million words (queryable at runtime via [__memtop()](#func---memtop)) can be accessed via brackets "[" and "]".

- [User definable functions](#function), which can define private variables, parameters, and also can optionally access namespaced instance variables. Recursion is NOT supported.

- Numbers are in normal decimal, however if you prefix a '$x' or '0x' to them, they will be hexadecimal (e.g. $x90, 0xDEADBEEF, etc).

- You may specify the ASCII value of a character using $'c' or 'c' (where c is the character). Multibyte characters are also supported using 'abc'.

- If you wish to generate a mask of 1 bits in integer, you can use $~X, for example $~7 is 127, $~8 is 255, $~16 is 65535, etc.

- Comments can be specified using:
  
  + // comments to end of line
  + /* comments block of code that span lines or be part of a line */

***

**Operator reference**

Listed from highest precedence to lowest (but one should use parentheses whenever there is doubt!):

- **[ ]**

```
z=x[y];
x[y]=z;
```

- **!value** -- returns the logical NOT of the parameter (if the parameter is 0.0, returns 1.0, otherwise returns 0.0).

- **-value** -- returns value with a reversed sign (-1 * value).

- **+value** -- returns value unmodified.

- **base ^ exponent** -- returns the first parameter raised to the power of the second parameter. This is also available the function pow(x,y)

- **numerator % denominator** -- divides two values as integers and returns the remainder.

- **value << shift_amt** -- converts both values to 32 bit integers, bitwise left shifts the first value by the second. Note that shifts by more than 32 or less than 0 produce undefined results.

- **value >> shift_amt** -- converts both values to 32 bit integers, bitwise right shifts the first value by the second, with sign-extension (negative values of y produce non-positive results). Note that shifts by more than 32 or less than 0 produce undefined results.

- **value / divisor** -- divides two values and returns the quotient.

- **value * another_value** -- multiplies two values and returns the product.

- **value - another_value** -- subtracts two values and returns the difference.

- **value + another_value** -- adds two values and returns the sum.

- **a | b** -- converts both values to integer, and returns bitwise OR of values.

- **a & b** -- converts both values to integer, and returns bitwise AND of values.

- **a ~ b** -- converts both values to 32 bit integers, bitwise XOR the values.

- **value1 == value2** -- compares two values, returns 1 if difference is less than 0.00001, 0 if not.

- **value1 === value2** -- compares two values, returns 1 if exactly equal, 0 if not.

- **value1 != value2** -- compares two values, returns 0 if difference is less than 0.00001, 1 if not.

- **value1 !== value2** -- compares two values, returns 0 if exactly equal, 1 if not.

- **value1 < value2** -- compares two values, returns 1 if first parameter is less than second.

- **value1 > value2** -- compares two values, returns 1 if first parameter is greater than second.

- **value1 <= value2** -- compares two values, returns 1 if first is less than or equal to second.

- **value1 >= value2** -- compares two values, returns 1 if first is greater than or equal to second.

- **y || z** -- returns logical OR of values. If 'y' is nonzero, 'z' is not evaluated.

- **y && z** -- returns logical AND of values. If 'y' is zero, 'z' is not evaluated.

- **y ? z** _-- how conditional branching is done -- similar to C's if/else_
  **y ? z : x**

If y is non-zero, executes and returns z, otherwise executes and returns x (or 0.0 if _': x'_ is not specified).

Note that the expressions used can contain multiple statements within parentheses, such as:

```
x % 5 ? (
f += 1;
x *= 1.5;
) : (
f = max(3,f);
x = 0;
);
```

- **y = z** -- assigns the value of 'z' to 'y'. 'z' can be a variable or an expression.
- **y \*= z** -- multiplies two values and stores the product back into 'y'.
- **y /= divisor** -- divides two values and stores the quotient back into 'y'.
- **y %= divisor** -- divides two values as integers and stores the remainder back into 'y'.
- **base ^= exponent** -- raises first parameter to the second parameter-th power, saves back to 'base'
- **y += z** -- adds two values and stores the sum back into 'y'.
- **y -= z** -- subtracts 'z' from 'y' and stores the difference into 'y'.
- **y |= z** -- converts both values to integer, and stores the bitwise OR into 'y'
- **y &= z** -- converts both values to integer, and stores the bitwise AND into 'y'
- **y ~= z** -- converts both values to integer, and stores the bitwise XOR into 'y'

Some key notes about the above, especially for C programmers:

- ( and ) (vs { } ) -- enclose multiple statements, and the value of that expression is the last statement within the block:

```
z = (a = 5; b = 3; a+b;); // z will be set to 8, for example
```

- Conditional branching is done using the ? or ? : operator, rather than if()/else.

```
a < 5 ? b = 6; // if a is less than 5, set b to 6
a < 5 ? b = 6 : c = 7; // if a is less than 5, set b to 6, otherwise set c to 7
a < 5 ? ( // if a is less than 5, set b to 6 and c to 7
b = 6;
c = 7; );
```

- The ? and ?: operators can also be used as the lvalue of expressions:

```
(a < 5 ? b : c) = 8; // if a is less than 5, set b to 8, otherwise set c to 8
```

***

**Simple math functions**

- **sin(angle)** -- returns the Sine of the angle specified (specified in radians -- to convert from degrees to radians, multiply by $pi/180, or 0.017453)
- **cos(angle)** -- returns the Cosine of the angle specified (specified in radians).
- **tan(angle)** -- returns the Tangent of the angle specified (specified in radians).
- **asin(x)** -- returns the Arc Sine of the value specified (return value is in radians).
- **acos(x)** -- returns the Arc Cosine of the value specified (return value is in radians).
- **atan(x)** -- returns the Arc Tangent of the value specified (return value is in radians).
- **atan2(x,y)** -- returns the Arc Tangent of x divided by y (return value is in radians).
- **hypot(x,y)** -- returns the hypotenuse of a right-angled triangle whose legs are x and y.
- **hypotFast(x,y)** -- returns the hypotenuse of a right-angled triangle whose legs are x and y using sqrt(x * x + y * y).
- **sqr(x)** -- returns the square of the parameter (similar to x*x, though only evaluating x once).
- **sqrt(x)** -- returns the square root of the parameter.
- **pow(x,y)** -- returns the first parameter raised to the second parameter-th power. Identical in behaviour and performance to the ^ operator.
- **exp(x)** -- returns the number e (approx 2.718) raised to the parameter-th power. This function is significantly faster than pow() or the ^ operator
- **log(x)** -- returns the natural logarithm (base e) of the parameter.
- **log10(x)** -- returns the logarithm (base 10) of the parameter.
- **abs(x)** -- returns the absolute value of the parameter.
- **min(x,y)** -- returns the minimum value of the two parameters.
- **max(x,y)** -- returns the maximum value of the two parameters.
- **sign(x)** -- returns the sign of the parameter (-1, 0, or 1).
- **rand(x)** -- returns a psuedorandom number between 0 and the parameter.
- **round(x)** -- rounds the value to the nearest integer (round(-5.6) == 3, round(2.1) == 2).
- **floor(x)** -- rounds the value to the lowest integer possible (floor(3.9) == 3, floor(-3.1) == -4).
- **ceil(x)** -- rounds the value to the highest integer possible (ceil(3.1) == 4, ceil(-3.9) == -3).
- **expint(x)** -- returns the exponential integral of the parameter.
- **expintFast(x)** -- returns exponential integral approximation of the parameter using lookup table.
- **invsqrt(x)** -- returns inverse square root (1/sqrt(x)) of the parameter.
- **invsqrtFast(x)** -- returns a fast inverse square root (1/sqrt(x)) approximation of the parameter.

***

**Loops**

Looping is supported in EEL2 via the following functions:

- **loop(count,code)**

```
loop(32,
    r += b;
    b = var * 1.5;
    );
```

Evaluates the first parameter once in order to determine a loop count. If the loop count is less than 1, the second parameter is not evaluated.
Implementations may choose to limit the number of iterations a loop is permitted to execute (usually such limits are in the millions and should rarely be encountered).

The first parameter is only evaluated once (so modifying it within the code will have no effect on the number of loops). For a loop of indeterminate length, see [while()](#while) below.

- **while(code)**

```
while(
a += b;
b *= 1.5;
a < 1000; // as long as a is below 1000, we go again.
);
```

Evaluates the first parameter until the last statement in the code block evaluates to zero.

Implementations may choose to limit the number of iterations a loop is permitted to execute (usually such limits are in the millions and should rarely be encountered).

- **while(condition) ( code )**

```
while ( a < 1000 ) (
a += b;
b *= 1.5;
);
```

Evaluates the parameter, and if nonzero, evaluates the following code block, and repeats. This is similar to a C style while() construct.

Implementations may choose to limit the number of iterations a loop is permitted to execute (usually such limits are in the millions and should rarely be encountered).

***

**User defined functions and namespace pseudo-objects**

EEL2 supports user defined functions, as well as some basic object style data access.

Functions can be defined anywhere in top level code (i.e. not within an existing () block, but before or after existing code). Functions are not able to be called recursively -- this is enforced by functions only being able to call functions that are declared before the current function, and functions not being able to call themselves. Functions may have 0 to 40 parameters. To define a function, use the following syntax:

```
function getSampleRate()
(
  srate; // return srate
);

function mySine(x)
(
  // taylor approximation
  x - (x^3)/(3*2) + (x^5)/(5*4*3*2) - (x^7)/(7*6*5*4*3*2) + (x^9)/(9*8*7*6*5*4*3*2);
);

function calculateSomething(x y)
(
  x += mySine(y);
  x/y;
);
```

Which would then be callable from other code, such as:

```
y = mySine($pi * 18000 / getSampleRate());
z = calculateSomething(1,2);
```

Note that the parameters for functions are private to the function, and will not affect global variables. If you need more private variables for a function, you can declare additional variables using a local() statement between the function declaration and the body of the function. Variables declared in the local() statement will be local to that function, and persist across calls of the function. Example:

```
function mySine(x) local(lastreq lastvalue)
(
  lastreq != x ? (
  lastreq = x; // save last input
  // taylor approximation
  lastvalue = x - (x^3)/(3*2) + (x^5)/(5*4*3*2) - (x^7)/(7*6*5*4*3*2) + (x^9)/(9*8*7*6*5*4*3*2);
  );
  lastvalue; // result of function is cached value
);
```

In the above example, mySine() will cache the last value used and not perform the calculation if the cached value is available. Note that the local variables are initialized to 0, which happens to work for this demonstration but if it was myCosine(), additional logic would be needed.

EEL2 also supports relative namespaces on global variables, allowing for pseudo object style programming. Accessing the relative namespace is accomplished either by using a "this." prefix for variable/function names, or by using the instance() declaration in the function definition for variable names:

```
function set_foo(x) instance(foo)
(
  foo = x;
);
// or
function set_foo(x)
(
  this.foo = x;
);

whatever.set_foo(32); // whatever.foo = 32;
set_foo(32); // set_foo.foo = 32;

function test2()
(
  this.set_foo(32);
);
whatever.test2(); // whatever.foo = 32
```

Additionally functions can use the "this.." prefix for navigating up the namespace hierarchy, such as:

```
function set_par_foo(x)
(
  this..foo = x;
);
a.set_par_foo(1); // sets foo (global) to 1
a.b.set_par_foo(1); // sets  a.foo to 1
```

***

**Advanced Mathematics, Optimization, Linear algebra and Signal Processing Functions**

**Short time Fourier transform**

_Compute STFT on equispace sampled signal, following function's state variables occupy virtual local address space_

- **stftCheckMemoryRequirement(start_index_idx, fft_length, overlap, wndPow)**
  Return memory requirement from given spectral transform parameters.
- **stftInit(start_index_idx, start_indexTransform)**
  Return latency of the STFT, the latency value provide the hint that how big the input frame is.
- **stftForward(start_index_buffer, start_index_idx, start_indexTransform, retPolarCart)**
  Perform forward STFT on a framed signal at the offset specified by the first parameter(**start_index_buffer**), function accept return either Polar grid or Cartesian grid.
- **stftBackward(start_index_buffer, start_index_idx, start_indexTransform, retPolarCart)**
  Perform inverse STFT on the spectrum at the offset specified by the first parameter( )**start_index_buffer**), **retPolarCart** must be the same as the one used at **stftForward()**.

Example:

```
fftsize = 1024;
stftIndexLeft = 0;
stftStructLeft = 50;
requiredSamples = stftInit(stftIndexLeft, stftStructLeft);
memreq = stftCheckMemoryRequirement(stftIndexLeft, fftsize, 4, 1.5);
// Fill up buffer
...
// Perform forward transform
spectralLen = stftForward(inBufLeft, stftIndexLeft, stftStructLeft, 1);
// Modify spectrum
...
// Perform inverse transform
error = stftBackward(inBufLeft, stftIndexLeft, stftStructLeft, 1);
```

The STFT provided also provide windowing, so your code is not required to window the overlapped results, but simply fill up the buffer. See the example codes for more information.

**Fast Fourier transform**

- **fft(start_index, size), ifft(start_index, size)**
  **fft_real(start_index, size), ifft_real(start_index, size)**
  **fft_permute(start_index, size), fft_ipermute(start_index, size)**
  Example:

```
buffer=0;
fft(buffer, 512);
fft_permute(buffer, 512);
buffer[32]=0;
fft_ipermute(buffer, 512);
ifft(buffer, 512);
// need to scale output by 1/512.0, too.
```

Performs a FFT (or inverse in the case of ifft()) on the data in the local memory buffer at the offset specified by the first parameter. The size of the FFT is specified by the second parameter, which must be 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, or 32768. The outputs are permuted, so if you plan to use them in-order, call **fft_permute(buffer, size)** before and **fft_ipermute(buffer,size)** after your in-order use. Your inputs or outputs will need to be scaled down by 1/size, if used.

Note that the FFT/IFFT require real/imaginary input pairs (so a 256 point FFT actually works with 512 items).

Note that the FFT/IFFT must NOT cross a 65,536 item boundary, so be sure to specify the offset accordingly.

The fft_real()/ifft_real() variants operate on a set of size real inputs, and produce size/2 complex outputs. The first output pair is DC,nyquist. Normally this is used with fft_permute(buffer,size/2).

- **convolve_c(dest,src,size)**
  Used to convolve two buffers, typically after FFTing them. convolve_c works with complex numbers. The sizes specify number of items (the number of complex number pairs).

Note that the convolution must NOT cross a 65,536 item boundary, so be sure to specify the offset accordingly.

**Fractional delay line**

_Delay input signal by arbitrary amount_

- **fractionalDelayLineInit(start_index, maxDelay)**
  Return memory requirement from given max delay parameters.
- **fractionalDelayLineSetDelay(start_index, delaySmps)**
  Set fractional delay in samples.
- **fractionalDelayLineClear(start_index)**
  Clear all internal delay.
- **fractionalDelayLineProcess(start_index, xn)**
  Save **xn** sample value into internal buffer, and return delayed value

Example:

```
ptr1 = 0;
// Initialize
req = fractionalDelayLineInit(ptr1, 1024);
// Set delay
fractionalDelayLineSetDelay(ptr1, 15.1);
// Performing some sample-by-sample processing
delayedOutput = fractionalDelayLineProcess(ptr1, input);
```

The fractional delay can be change anytime during processing. See the example codes for more information.

**Direct form FIR Filter**

_Perform FIR filtering on purely time domain_

- **FIRInit(start_index, hLen)**
  Return memory requirement for the FIR filter from given impulse response length.
- **FIRProcess(start_index, xn, coefficients)**
  Perform FIR filtering on **xn** using provided coefficients.

Example:

```
coefficients = 0;
// Import coefficients from string
hLen = importFLTFromStr("0.0446871259733802,0.207161817163176,0.333544729155861,0.333544729155861,0.207161817163176,0.0446871259733802,-0.0538264189367632,-0.05743465917334062", coefficients);
ptr1 = coefficients + hLen;
req = FIRInit(ptr1, hLen);
// Processing
output = FIRProcess(ptr1, input, coefficients);
```

**Real time FFT convolution**

*Perform FIR filtering on frequency domain with non-uniform partitioned convolution algorithm*

- **Conv1DInit(buffer_len, hLen, h1, h2, h3, h4)**
  
  Initialization of FFT convolver.
  
  User can specify one or two or four impulse responses, if user specify one, the function assume there is only one input, if specify two, the function assume the input is two channel, if specify four, the function assume the input is two channel and a cross-channel summation will be perform after convolving each input with two impulse responses
  
  Specify the buffer length with buffer_len for buffer that stored for FFT convolver to read and write, user will have to allocate a buffer with buffer_len during processing
  
  hLen is impulse response length
  
  Memory of FFT convolver are being stored with limited slot, user can allocate 1024 FFT convolvers
  Return slot ID

- **Conv1DProcess(start_index, x1, x2)**
  Perform FIR filtering using FFT convolution.
  
  If user specify one impulse response during initialization, then only one input needs to be supplied.
  
  If user specify two or four impulse responses during initialization, then two input needs to be supplied.

- **Conv1DFree(slot_id)**
  Free up occupied memory
  
  The specify slot_id must match the one return by previous initialization

Example:

```
// Initialization
requiredSamples = 1024;
h1 = 0;
impulseResponseLength = importFLTFromStr("IMPULSE1", h1);
h2 = h1 + impulseResponseLength;
impulseResponseLength = importFLTFromStr("IMPULSE2", h2);
h3 = h2 + impulseResponseLength;
impulseResponseLength = importFLTFromStr("IMPULSE3", h3);
h4 = h3 + impulseResponseLength;
impulseResponseLength = importFLTFromStr("IMPULSE4", h4);
convId = Conv1DInit(requiredSamples, impulseResponseLength, h1, h2, h3, h4);
inBufLeft = h4 + impulseResponseLength;
outBufLeft = inBufLeft + requiredSamples;
inBufRight = outBufLeft + requiredSamples;
outBufRight = inBufRight + requiredSamples;
// Processing
inBufLeft[bufpos] = spl0;
spl0 = outBufLeft[bufpos];
inBufRight[bufpos] = spl1;
spl1 = outBufRight[bufpos];
bufpos += 1;
bufpos >= requiredSamples ?
(
  Conv1DProcess(convId, inBufLeft, inBufRight);
  idx = 0;
  loop(requiredSamples,
  outBufLeft[idx] = inBufLeft[idx];
  outBufRight[idx] = inBufRight[idx];
  idx+=1);
  bufpos = 0;
);
```

**Butterworth subband transform**

_Decompose signal to N channels time domain subbands_

- **IIRBandSplitterInit(start_index, fs, band1, band2, ...)**
  Algorithm support [1...7] bands decomposition, user pass are allowed passing maximum 7 frequency cut-off to function arguments, return memory requirement for the internal band splitting filters.
- **IIRBandSplitterClearState(start_index)**
  Reset the band splitting filters states.
- **IIRBandSplitterProcess(start_index, xn, band1, band2, ...)**
  Perform band splitting on xn, algorithm support [1...7] bands decomposition, user pass are allowed passing maximum 8 variable for transformed value returning.

Example:

```
iirBPS = 0;
// Initialize 3 bands
reqSize = IIRBandSplitterInit(iirBPS, 48000, 4000, 12000); // Specify 2 cut-off frequencies
kDelta = iirBPS + reqSize;
kDelta[0] = 1; // Kronecker delta
sigLen = 1024;
low = 0;
mid = 0;
high = 0;
// Processing
// ...
IIRBandSplitterProcess(iirBPS, kDelta[idx], low, mid, high);
// ...
// Reset states
IIRBandSplitterClearState(iirBPS);
```

**Constant-Q polyphase filterbank**

_Decompose signal to N channels decimated time domain subbands with logarithmic centre frequencies_

- **InitPolyphaseFilterbank(fs, N, m, start_index1, start_index2)**
  
  Initialization of polyphase filterbank
  User supply fs as sample rate, N as channels, m as number of polyphase components, function return memory requirement of all the internals buffer.

- **PolyphaseFilterbankChangeWarpingFactor(start_index, sample_rate, warpingFactor)**
  Change dispersive delay line warping factor
  
  start_index is filterbank memory block, sample_rate is sample rate of signal user like to break down into subbands, warpingFactor is dispersive delay line warping factor, range is [-1, 1]

- **PolyphaseFilterbankGetPhaseCorrector(start_index1, percentage, start_index2)**
  
  Function to obtain impulse response that correct dispersive delay line
  start_index1 is filterbank memory block, percentage corresponding to how much the nonlinear phase user require to correct, start_index2 is memory area for phase correction impulse response.

- **PolyphaseFilterbankGetDecimationFactor(start_index1, start_index2)**
  
  Nonuniform polyphase subbands has nonuniform sampling periods, obtaining sampling periods is required to perform accurate analysis of signal and re-synthesis
  
  start_index1 is filterbank memory block, where start_index2 is decimation factor.

- **PolyphaseFilterbankAnalysisMono(start_index1, start_index2, start_index3, x1)**
  
  Function to perform analysis of a signal
  
  start_index1 is filterbank memory block, start_index2 is filterbank output of x1, start_index3 is current indices of decimation, x1 is input sample

- **PolyphaseFilterbankAnalysisStereo(start_index1, start_index2, start_index3, start_index4, start_index5, x1, x2)**
  
  Function to perform analysis of two input signals
  
  start_index1 is filterbank memory block 1, start_index2 is filterbank memory block2, start_index3 is filterbank output of x1, start_index4 is filterbank output of x2, start_index5 is current indices of decimation, x1 and x2 are input sample

- **PolyphaseFilterbankSynthesisMono(start_index1, start_index2, y1)**
  
  Function to perform synthesis of subbands of a signal
  
  start_index1 is filterbank memory block 1, start_index2 is filterbank output of x1, y1 is synthesis output of filterbank

- **PolyphaseFilterbankSynthesisStereo(start_index1, start_index2, start_index3, start_index4, y1, y2)**
  
  Function to perform synthesis of subbands of two input signal
  
  start_index1 is filterbank memory block 1, start_index2 is filterbank memory block2, start_index3 is filterbank output of x1, start_index4 is filterbank output of x2, y1 and y2 are synthesis output of filterbank

Example:

```
// Initialization
N = 8;
m = 3;
pfb1 = 0;
memSize = InitPolyphaseFilterbank(srate, N, m, pfb1, pfb2);
PolyphaseFilterbankChangeWarpingFactor(pfb1, srate, 0.99);
phaseCorrFilt = pfb2 + memSize;
corrFiltLen = PolyphaseFilterbankGetPhaseCorrector(pfb1, 0.5, phaseCorrFilt);
i = 0;
loop(corrFiltLen, 
//printf("%1.7f,", phaseCorrFilt[i]);
i += 1);
Sk = phaseCorrFilt + corrFiltLen;
PolyphaseFilterbankGetDecimationFactor(pfb1, Sk);
subbandOutputLeft = Sk + N;
subbandOutputRight = subbandOutputLeft + N;
decimationCnt = subbandOutputRight + N;
eqGain = decimationCnt + N;
eqGain[0] = 1;
eqGain[1] = 0;
eqGain[2] = 1;
eqGain[3] = 0;
eqGain[4] = 1;
eqGain[5] = 0;
eqGain[6] = 1;

// Processing
PolyphaseFilterbankAnalysisStereo(pfb1, pfb2, subbandOutputLeft, subbandOutputRight, decimationCnt, spl0, spl1);
// Modify subbands
i = 0;
loop(N, 
decimationCnt[i] == Sk[i] ? (
subbandOutputLeft[i] = subbandOutputLeft[i] * eqGain[i];
subbandOutputRight[i] = subbandOutputRight[i] * eqGain[i];
);
i += 1);
PolyphaseFilterbankSynthesisStereo(pfb1, pfb2, subbandOutputLeft, subbandOutputRight, y1, y2);
spl0 = y1;
spl1 = y2;

// Reset states
PolyphaseFilterbankClearState(iirBPS);
```

The subband transform supposed work on time series. See the example codes for more information.

**Autoregressive model**

_Decompose signal to N channels time domain subbands, following function's state variables occupy virtual local address space_

- **arburgCheckMemoryRequirement(model_order, require_predictionState)**
  Check memory requirement for autoregressive calculation
  
  model_order is autoregressive model order, require_predictionState is boolean option that specific user require to get forward and backward state

- **arburgTrainModel(start_index1, model_order, require_predictionState, start_index2, input_len)**
  
  Fitting Burg model with input
  
  start_index1 is memory block for ARBurg, model_order is autoregressive model order, require_predictionState is boolean option that specific user require to get forward and backward state, start_index2 is input signal vector, input_len is signal length

- **arburgGetPredictionReflectionCoeff(start_index1, start_index2, start_index3)**
  
  Obtain model coefficient from ARBurg memory
  start_index1 is memory block for ARBurg, start_index2 is vector of prediction coefficients, start_index3 is vector of reflection coefficients.

- arburgPredictForward(start_index)
  
  Perform forward prediction of signal just being fitted
  
  start_index is memory block for ARBurg

- arburgPredictBackward(start_index)
  
  Perform backward predictionof signal just being fitted
  
  start_index is memory block for ARBurg

Example:

```
r_ = 0;
// Define time domain signal r_
inLen = importFLTFromStr("ARRAY OF SAMPLE", r_);
order = 31;
getPredictionState = 1;
memReq = arburgCheckMemoryRequirement(order, getPredictionState);
burg = r_ + inLen;
predictionCoefficients = burg + memReq;
reflectionCoefficient = predictionCoefficients + (order + 1);
arburgTrainModel(burg, order, getPredictionState, r_, inLen);
arburgGetPredictionReflectionCoeff(burg, predictionCoefficients, reflectionCoefficient);
printf("predictionCoefficients: ");
idx = 0;
loop((order + 1),
printf("%1.8f, ", predictionCoefficients[idx]);
idx += 1;
);
printf("\nreflectionCoefficient: \n");
idx = 0;
loop((order + 1),
printf("%1.8f, ", reflectionCoefficient[idx]);
idx += 1;
);
printf("\nBackward regression: \n");
idx = 0;
loop(1000,
printf("%1.8f, ", arburgPredictBackward(burg));
idx += 1;
);
printf("\nForward regression: \n");
idx = 0;
loop(1000,
printf("%1.8f, ", arburgPredictForward(burg));
idx += 1;
);
```

** Miscellaneous

- unwrap

- cplxpair

- zp2sos

- tf2sos

- roots

- eqnerror

- firls

** Linear algebra

- **det(start_index, m, n)** -- Compute the determinant of the input matrix **start_index**, m and n must be identical, otherwise return -1

- **rank(start_index, m, n)** -- returns the rank of the input matrix **start_index**.

- **transpose(start_indexIn, start_indexOut, m, n)** -- Out-of-place transpose of input matrix **start_indexIn**.

- **inv(start_indexIn, m, n, start_indexOut)** -- Compute matrix inversion of input matrix **start_indexIn**, m and n must be identical, otherwise return -1

- **pinv(start_indexIn, m, n, start_indexOut, start_indexMatSize)** -- Compute Moore–Penrose inverse of input matrix **start_indexIn**

- **mldivide(A, m1, n1, B, m, n, Out, MatSize)** -- Solve linear system or solves the system of linear equations Ax = B for x

- **mrdivide(A, m1, n1, B, m, n, Out, MatSize)** -- Solve linear system or solves the system of linear equations xA = B for x

- cholesky

- inv_chol

- pinv_svd

- pinv_fast

** Optimization

Following are constrained optimization, altough ordinary least squares(OLS) is a type of optimization, the built-in routines are included in Linear algebra section.
The actual usage is somehow similar to the Matlab counterparts.

- **lsqlin(C, d, A, b, problemLen, inequalityLen, solution)** -- Solves the linear system C*x = d in the least-squares sense, subject to A*x ≤ b, we don't support equalities input, user need to convert linear equalities to inequalities and deal with input arguments manually, routine will return objective function value at the solution, the **fval** to stack.
- **quadprog(problemLen, H, f, inequalityLen, A, b, equalityLen, Aeq, beq, lb, ub, solution)** -- minimizes 1/2*x'*H*x + f'*x subject to the restrictions A*x ≤ b. Routine will return objective function value at the solution, the **fval** to stack.

The function description is not intuitive, for actual use, please check out the examples.
**Important**
The solution computed by the optimization algorithm is not identical to Matlab counterpart, but the solution is quite close, the optimization routine author had use the routine on optimization-based digital filter design, the design result is almost identical to the filter design Matlab routines.

***

**Strings**

Strings can be specified as literals using quotes, such as "This is a test string". Much of the syntax mirrors that of C: you must escape quotes with backslashes to put them in strings ("He said \"hello, world\" to me"), multiple literal strings will be automatically concatenated by the compiler. Unlike C, quotes can span multiple lines. Internal strings lookup was implemented in a growable container.

Strings are always referred to by a number, so one can reference a string using a normal EEL2 variable:

```
x = "hello world";
printf("%s\n", x);
printf("Val = %f", 3.4488);
```

**String functions**

- **resetStringContainers()** -- remove all string from VM in case of dangling string bother you

- **importFLTFromStr(str, start_index)** -- import string characters to double precision floating pointer array

- **strlen(str)** -- returns length of string

- **strcmp(str, str2)** -- compares str to str2, case sensitive, returns -1, 0, or 1

- **stricmp(str, str2)** -- compares str to str2, ignoring case, returns -1, 0, or 1

- **strncmp(str, str2, maxlen)** -- compares str to str2 up to maxlen bytes, case sensitive, returns -1, 0, or 1

- **strnicmp(str, str2, maxlen)** -- compares str to str2 up to maxlen bytes, ignoring case, returns -1, 0, or 1

- **printf(str, format, ...)** -- print formated str to CLI, converting format strings:

- **sprintf(str, format, ...)** -- copies format to str, converting format strings:
  
  + %% = %
  + %s = string from parameter
  + %d = parameter as integer
  + %i = parameter as integer
  + %u = parameter as unsigned integer
  + %x = parameter as hex (lowercase) integer
  + %X = parameter as hex (uppercase) integer
  + %c = parameter as character
  + %f = parameter as floating point
  + %e = parameter as floating point (scientific notation, lowercase)
  + %E = parameter as floating point (scientific notation, uppercase)
  + %g = parameter as floating point (shortest representation, lowercase)
  + %G = parameter as floating point (shortest representation, uppercase)

Many standard [C printf()](http://www.cplusplus.com/reference/cstdio/printf/) modifiers can be used, including:

    + %.10s = string, but only print up to 10 characters
    + %-10s = string, left justified to 10 characters
    + %10s = string, right justified to 10 characters
    + %+f = floating point, always show sign
    + %.4f = floating point, minimum of 4 digits after decimal point
    + %10d = integer, minimum of 10 digits (space padded)
    + %010f = integer, minimum of 10 digits (zero padded)

The **printf** will print the value to **stdout** if the Host application is a CLI application, the string will be printed to logcat if the Host application have no GUI and not necessary limit to CLI application, however, if Host application have a GUI, for example, a executable or VST plugin, author had implement a string circular buffer to simulate the way CLI behave.

- **match(needle, haystack, ...)** -- search for needle in haystack
  **matchi(needle, haystack, ...)** -- search for needle in haystack (case insensitive)

For these you can use simplified regex-style wildcards:

    * = match 0 or more characters
    *? = match 0 or more characters, lazy
    + = match 1 or more characters
    +? = match 1 or more characters, lazy
    ? = match one character

Examples:

    match("*blah*", "this string has the word blah in it") == 1
    match("*blah", "this string ends with the word blah") == 1

You can also use format specifiers to match certain types of data, and optionally put that into a variable:

    %s means 1 or more chars
    %0s means 0 or more chars
    %5s means exactly 5 chars
    %5-s means 5 or more chars
    %-10s means 1-10 chars
    %3-5s means 3-5 chars.
    %0-5s means 0-5 chars.
    %x, %d, %u, and %f are available for use similarly
    %c can be used, but can't take any length modifiers
    Use uppercase (%S, %D, etc) for lazy matching 

The variables can be specified as additional parameters to match(), or directly within {} inside the format tag (in this case the variable will always be a global variable):

    match("*%4d*","some four digit value is 8000, I say",blah)==1 && blah == 8000
    match("*%4{blah}d*","some four digit value is 8000, I say")==1 && blah == 8000

***

**Misc Utility**
**Memory Utility**

- **freembuf(top)**
  The freembuf() function provides a facility for you to notify the memory manager that you are no longer using a portion of the local memory buffer.

For example, if the user changed a parameter on your effect halving your memory requirements, you should use the lowest indices possible, and call this function with the highest index you are using plus 1, i.e. if you are using 128,000 items, you should call freembuf(128001); If you are no longer using any memory, you should call freembuf(0);

Note that calling this does not guarantee that the memory is freed or cleared, it just provides a hint that it is OK to free it.

- **memcpy(dest,source,length)**
  The memcpy() function provides the ability to quickly copy regions of the local memory buffer. If the buffers overlap and either buffer crosses a 65,536 item boundary, the results may be undefined.

- **memset(dest,value,length)**
  The memset() function provides the ability to quickly set a region of the local memory buffer to a particular value.

- **__memtop()** -- returns the maximum memory words available to the script (read/writes to __memtop()[-1] will succeed, but __memtop()[0] will not)

**Stack**
A small (approximately 32768 item) user stack is available for use in code:

- **stack_push(value)**
  Pushes value onto the user stack, returns a reference to the value.

- **stack_pop(value)**
  Pops a value from the user stack into value, or into a temporary buffer if value is not specified, and returns a reference to where the stack was popped. Note that no checking is done to determine if the stack is empty, and as such stack_pop() will never fail.

- **stack_peek(index)**
  Returns a reference to the item on the top of the stack (if index is 0), or to the Nth item on the stack if index is greater than 0.

- **stack_exch(value)**
  Exchanges a value with the top of the stack, and returns a reference to the parameter (with the new value).
