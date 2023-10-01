//
//  main.h
//
//  Authors: Nicolae Mihalache & François Vigneron
//
//  This software is released under BSD licence, with an attribution clause (see Licence file).
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] R. Anton, N. Mihalache & F. Vigneron. Fast evaluation of real and complex polynomials. 2022.
//        https://hal.archives-ouvertes.fr/hal-03820369
//
//  Copyright 2022 Univ. Paris-Est Créteil & Univ. de Reims Champagne-Ardenne.
//

/**
 \mainpage Fast Evaluation of Real and Complex Polynomials
 
 \tableofcontents
 
 This documentation is availlable online at https://fvigneron.github.io/FastPolyEval.
 It also exists in PDF form: [FastPolyEval_doc.pdf](https://fvigneron.github.io/FastPolyEval/FastPolyEval_doc.pdf).
 The source code can be found on [GitHub](https://github.com/fvigneron/FastPolyEval/).

 \section sectIntro Abstract
 <b> @c FastPolyEval is a library, written in C, that aims at evaluating polynomials very efficiently, without compromising the accuracy of the result.
 It is based on the @c FPE algorithm (Fast Polynomial Evaluator) introduced in [reference [1]](https://hal.archives-ouvertes.fr/hal-03820369)</b> (see Section \ref sectRefCopy).
 
 In  @c FastPolyEval, the computations are done for real or complex numbers, in floating point arithmetic, with a fixed precision, which can be one
 of the machine types [FP32] (https://en.wikipedia.org/wiki/Single-precision_floating-point_format),
 [FP64] (https://en.wikipedia.org/wiki/Double-precision_floating-point_format),
 [FP80] (https://en.wikipedia.org/wiki/Extended_precision)
 or an arbitrary precision based on [MPFR] (https://www.mpfr.org).
 
 Evaluations are performed on arbitrary (finite...) sets of points in the complex plane or along the real line, without geometrical constraints.
 
 The average speed-up achieved by @c FastPolyEval over Hörner's scheme is illustrated on Figure 1.

 \image html HornerFPE.png Figure 1: Speed gain of @c FastPolyEval versus Hörner for O(d) evaluations in precision p=53 MPFR.
 \image latex HornerFPE.png Speed gain of @c FastPolyEval versus Hörner for O(d) evaluations  in precision p=53 MPFR.
  
 \htmlonly <hr> \endhtmlonly

 \section sectTheory Complexity & accuracy

 @c FastPolyEval splits the evaluation process in two phases:
 - A first phase, called \b pre-processing, analyses the coefficients of the polynomial (actually only the exponents) and determines an evaluator, based on a parcimonious representation of the polynomial.
 - Subsequently, the \b evaluator is applied to each of the requested evaluation points. A second reduction is performed. Then the final result is computed.

 The complexity of the pre-processing phase is bounded by \f$ O(d \log d)\f$ where \f$ d\f$ is the degree of the polynomial \f$ P(z)\f$.
 It is independent of the precision used to express the coefficients or requested for the rest of the computations.
 The @c FastPolyEval library is backed by theoretical results that guaranty that the average arithmetic complexity of the final evaluator is of order
 \f[
 O\left( \sqrt{d(p+\log d)} \right)
 \f]
 where \f$ p\f$ is the precision of the computation, in bits.
 The averaging process corresponds to evaluation points \f$ z\f$ uniformly distributed either on the Riemann sphere or
 on the unit disk of the complex plane.
 The worst complexity of the evaluator does not exceed that of Hörner, i.e. \f$ O(d)\f$, and can drop as low as \f$ O(\log^2 d)\f$
 in some favorable cases (see Figure 1 in \ref sectIntro).

 The memory requirement of @c FastPolyEval is minimal.
 To handle a polynomial of degre \f$ d\f$ with \f$ p\f$ bits of precision, the memory requirement is \f$ O(dp)\f$.
 An additional memory \f$ O(kp)\f$ is required to perform \f$ k\f$ evaluations in one call.

 Regarding accuracy, the theory guaranties that the relative error between the exact value of \f$ P(z)\f$ and
 the computed value does not exced \f$ 2^{-p-c-1}\f$ where \f$ c\f$ is the number of cancelled leading bits, i.e.
 \f[
 c = \begin{cases}
 0 & \text{if }  |P(z)| \geq  M(z),\\
 \left\lfloor \log_2\left| M(z) \right| \right\rfloor - \left\lfloor \log_2\left| P(z) \right| \right\rfloor & \text{else,}
 \end{cases}
 \f]
 where \f$ M(z)=\max |a_j z^j|\f$ is the maximum modulus of the monomials of \f$ P(z)\f$ and \f$ \lfloor\cdot\rfloor\f$ is the floor function.
  
 \note The geometric preprocessing uses only the exponents of the coefficients, which is a significantly smaller amount of data to process
 than reading all the bits of the coefficients.
 For a high-precision high-degree evaluation at a \b single point, @c FastPolyEval can often outperform Hörner.
 The preprocessing of the exponents followed by the evaluation with a high precision of a parcimonious representation of the polynomial
 is more efficient than handling indiscriminately all the coefficients.
 This fact does not contradict that Hörner is a theoretical best for one single evaluation. Hörner has the best estimate on complexity that holds for \b any evaluation point.
 The advantage of @c FastPolyEval can be substantial, but it holds \b on \b average. Note that the worst case for the evaluator (no coefficient dropped) has the same
 complexity as Hörner. However, in that case, one can prove that the average complexity is much better, typicaly \f$ O(\log^2 d) \f$.
 
 For further details, precise statements and proofs, please see [reference [1]] (https://hal.archives-ouvertes.fr/hal-03820369) in section \ref sectRefCopy.
 
 \htmlonly <hr> \endhtmlonly

 \section sectAlgo Description of the Fast Polynomial Evaluator (FPE) algorithm
 
 In this section, we describe briefly the mathematical principles at the foundation of @c FastPolyEval.
 
 \subsection sectIdeas General principles
 
 @c FastPolyEval relies on two general principles.
 -# <b> Lazy evaluation </b> in finite precision. When adding two floating points numbers with p bits, one can simply discard the smaller one if the orders of magnitude
 are so far apart that the leading bit of the smaller number cannot affect the bit in the last position of the larger number. When evaluating a polynomial of large
 degree, the orders of magnitude of the monomials tends to be extremely varied, which is a strong incentive to use lazy evaluations.
 -# <b> Geometric selection principle</b>. For a given evaluation point, evaluating each monomial and sorting them in decreasing order of magnitude would
 be a very inefficient way of implementing lazy additions. Instead, @c FastPolyEval can identify the leading monomials using a simple geometric method
 that allows us to factor most of the work in a pre-processing phase, which is only needed once. At each evaluation point, a subsequent reduction
 produces a very parcimonious representation of the polynomial (see Figure 2).
 
 \image html parsimonious.png Figure 2: Parsimony of the  representation used by @c FastPolyEval to evaluate a polynomial of degree 1000 (half-circle family).
 \image latex parsimonious.png  Parsimony of the  representation used by @c FastPolyEval to evaluate a polynomial of degree 1000 (half-circle family).
 
 \subsection sectGeometry Geometric selection principle

 For a given polynomial \f$ P\f$ of degree \f$ d\f$, one represents the modulus of the coefficients \f$ a_j\f$ in logarithmic coordinates,
 that is the scatter plot \f$ E_P\f$ of \f$ \log_2 |a_j|\f$ in function of \f$ j\in \{0,1,\ldots,d\}\f$.
 One then computes the concave envelope \f$ \hat{E}_P\f$ of \f$ E_P\f$ (obviously piecewise linear)
 and one identifies the strip \f$ S_\delta(\hat{E}_P)\f$ situated below \f$ \hat{E}_P\f$ and of vertical thickness
 \f[
 \delta = p + \lfloor \log_2 d \rfloor + 4.
 \f]
 Note that the thickness of this strip is mostly driven by the precision \f$ p\f$ (in bits) that will be used to evaluate \f$ P\f$.
 On Figure 3, the blue points represent \f$ E_P\f$, the cyan line is \f$ \hat{E}_P\f$ and the strip \f$ S_\delta\f$ is colored in light pink.
 
 \image html concave.png Figure 3: Concave cover of the oefficients of a polynomial (log-scale) and selection principle for a given precision.
 \image latex concave.png Concave cover of the oefficients of a polynomial (log-scale) and selection principle for a given precision.

 The coefficients of \f$ P\f$ that lie in the strip \f$ S_\delta\f$ are pertinent for an evaluation with precision \f$ p\f$.
 The others can safely be discarded because they will never be used at this precision.

 \subsection sectEval Parcimonious representation at an arbitrary point
 
 For a given \f$ z\in\mathbb{C}\f$, a parcimonious representation of \f$ P(z)=\sum\limits_{j=0}^d a_j z^j\f$ is obtained by selecting a proper subset of the indices \f$ J_p(z) \subset\{0,1,\ldots,d\}\f$
 and discarding the other monomials:
 \f[
 Q_p(z) = \sum_{j\in J_p(z)} a_j z^j.
 \f]
 For computations with a precision of \f$ p\f$ bits, the values of \f$ Q_p(z)\f$ and \f$ P(z)\f$ are equivalent because the relative error does not exceed \f$ 2^{-p-c-1}\f$
 where c is the number of cancelled leading bits (defined in Section \ref sectTheory).

 The subset \f$ J_p(z)\f$ is obtained as a subset of the (indices of the) coefficients that respect the two following criteria in the logarithmic representation (see Figure 3):
 -# the coefficient belongs to the strip \f$ S_\delta\f$,
 -# the coefficient is above the longest segment of slope \f$ \tan\theta = -\log_2|z| \f$ contained in \f$ S_\delta\f$.

 The selecting segment is always tangent to the lower edge of \f$ S_\delta\f$ or, in case of ambiguity (e.g. if \f$ S_\delta\f$ is a parallelogram), it is required to be.
 For \f$ |z|=1\f$, the selecting segment is an horizontal line illustrated in solid red in Figure 3 above. The red dashed line corresponds to some value \f$ |z|<1\f$.
 Figure 4 illustrates how the selection principle at an arbitrary point \f$ z\in\mathbb{C}\f$ can be transposed from the analysis performed
 on the raw coefficients, i.e. for \f$ |z|=1\f$. In log-coordinates, one has indeed:
 \f$ \log_2\left( |a_j z^j | \right) = \log_2 |a_j| - j \tan\theta.\f$
 Selecting the coefficients such that \f$ \log_2\left( |a_j z^j | \right)\f$ differs from its maximum value by less than \f$ \delta\f$ (solid red in Figure 4) is equivalent to selecting
 those above the segment of slope \f$ \tan\theta\f$ in the original configuration (dashed line in Figure 3).
 
 \image html shear.png Figure 4: Transposition of the selection principle at an arbitrary evaluation point
 \image latex shear.png Transposition of the selection principle at an arbitrary evaluation point

 \subsection sectImplement Implementation of the algorithm
 
 @c FastPolyEval computes Figure 3 and the strip \f$ S_\delta\f$ in the pre-processing phase and keeps only the subset of the coefficients that belong to that strip
 (the so-called "good" set \f$ G_p\f$ in [reference [1]] (https://hal.archives-ouvertes.fr/hal-03820369)).
 In the evaluation phase, the subset of coefficients is thinned to \f$ J_p(z)\f$ for each evaluation point \f$ z\f$ using the 2nd criterion, thus providing
 the parcimonious representation \f$ Q_p(z)\f$ of \f$ P(z)\f$ at that point. The value of the polynomial is then computed using Hörner's method for \f$ Q_p\f$.
 Theoretical results guaranty that, on average, the number of coefficients kept in \f$ J_p(z)\f$ is only of order \f$ \sqrt{d \delta}\f$, which explains why the complexity
 of @c FastPolyEval is so advantageous.
 
 For further details, see [reference [1]] (https://hal.archives-ouvertes.fr/hal-03820369) in Section \ref sectRefCopy.
 
 \htmlonly <hr> \endhtmlonly

 \section sectInstall Installation of FastPolyEval

 @c FastPolyEval is distributed as source code, written in the @c C language.
 The source code can be downloaded from our public [github] (https://github.com/fvigneron/FastPolyEval) repository.
 
 \subsection sectPrereq Prerequisites
 
 First of all, don't panic, it's easy...
 
 Obvioulsy, you will need a computer and a @c C compiler, for example [GCC] (https://gcc.gnu.org), though any other @c C compiler will be fine.
 
 @c FastPolyEval relies on [MPFR] (https://www.mpfr.org) to handle floating point numbers with arbitrary precision.
 You need to install this library on your system if it is not already installed. Please follow the MPFR installation [instructions] (https://www.mpfr.org/mpfr-current/mpfr.html#Installing-MPFR) to do so.
 Note that MPFR has one main depency, [GMP] (https://gmplib.org), which will also need to be installed (see the GMP installation [instructions] (https://gmplib.org/manual/Installing-GMP)).
 If you do not have admin privileges on your machine, you can still build those libraries in a local directory and instruct your compiler to look for the appropriate library path.
 For example, if you use GCC, see the [documentation] (https://gcc.gnu.org/onlinedocs/gcc/Directory-Options.html) of option @c -L.
 
 Note that on most OS, one can also install those libraries through a package manager (see \ref sectOS).

 \note The libraries GMP and MPFR are \b free to download and use and are distributed under the Lesser General Public License, which enables developers of non-free programs to use GMP/MPFR in their programs.
 @c FastPolyEval is released under the BSD licence, with an attribution clause (see \ref sectRefCopy), which means that you are <b> free to download, use, modify, distribute or sell, provided you give us proper credit</b>.

 \subsection sectCompile Compile FastPolyEval
 
 Download the [source code] (https://github.com/fvigneron/FastPolyEval). In the command line, go to the @c code subfolder then type:
 \code
 make
 \endcode
 to get a basic help message. The default build command is:
 \code
 make fpe
 \endcode
 This will create and populate a @c bin subdirectory of the downloaded folder with an app, called @c FastPolyEval.
 You may move the executable file to a convenient location. Alternatively, consider adding this folder to your @c PATH.
 The @c Makefile does not attempt this step for you.
 See \ref sectUsage for additional help on how to use this app.
 
 If you prefer to compile by hand, the structure of your command line should be:
 \code
 gcc all_source_files.c -o bin/FastPolyEval -I each_source_sub_folder_containing_header_files -Wall -lm -lgmp -lmpfr -O3
 \endcode
 Note that the @c file.c and @c -I @c folder parameters must be repeated as necessary.
 The @c Makefile is simple enough and understanding it can provide additional guidance.
 
 @c FastPolyEval switches automatically from hardware machine precision to emulated high-precision using MPFR (see \ref sectNumberFormats). The format used for hardware numbers
 is chosen at the time of the compilation. You may edit the source file @c code/numbers/ntypes.h before the default build. Alternatively,
 you can execute
 \code
 make hardware
 \endcode
 This should generate an automatically edited copy of this header file in a temporary folder and launch successive compilations.
 This build will populate the @c bin subdirectory with three apps, called @c FastPolyEval_FP32, @c FastPolyEval_FP64
 and @c FastPolyEval_FP80 corresponding to each possible choice.
 
 \subsection sectOS Additional help for Linux, MacOS, Windows
 
 Github offers access to virtual machines that run on their servers. The corresponding tasks are called [workflows](https://docs.github.com/en/actions/using-workflows). Each workflow serves as a guidline (read the subsection @c job:steps in the @c yml files) to perform similar actions on your computer.
 The directory [.github/workflows](https://github.com/fvigneron/FastPolyEval/tree/main/.github/workflows) on our [GitHub repository](https://github.com/fvigneron/FastPolyEval) contains a @c build workflow for each major OS : Linux, MacOS and Windows.

 \subsection sectUninstall Uninstall FastPolyEval

 As the @c Makefile does not move the binaries and does not alter your @c PATH, cleanup is minimal.
 To remove the binary files generated in the build process, execute
 \code
 make clean
 \endcode
 To uninstall @c FastPolyEval from your system, simply delete the executable file and remove the downloaded sources.

 \htmlonly <hr> \endhtmlonly

 \section sectUsage Basic usage of FastPolyEval
 
@c FastPolyEval is used at the command line, either directly or through a shell script.
 A typical command line call has the following structure:
 \code
 FastPolyEval -task prec parameter [optional_parameter]
 \endcode
 where @c prec is the requested precision for the computations.
 
 \subsection sectHelp Onboard help system
 
 To call the onboard help with a list and brief description of all possible tasks, type
 \code
 FastPolyEval -help
 \endcode
 
 A detailed help exists for each task. For example, try @c FastPolyEval @c -eval @c -help.
 
 The different tasks belong to 3 groups:
 
 - \ref sectTasksPoly,
 - \ref sectTasksComplex,
 - \ref sectTasksFPE
 
 and will be detailed below.
 
 \note In absence of argument, the default behavior is a call for help. The previous messaged can also be displayed respectively with @c FastPolyEval  and  @c FastPolyEval @c -eval .
  
 \subsection sectNumberFormats Number formats
 
 The flags @c MACHINE_LOW_PREC and @c MACHINE_EXTRA_PREC in ntypes.h
 define what kind of machine numbers @c FastPolyEval defaults to for low-precision computations, and the corresponding precision threshold that defines how low is "low". Note that the main limitation of
 machine numbers is the limited range of exponents. The numerical limits can easily be reached when evaluating a polynomial of high degree. On the contrary, MPFR number
 are virtually unlimited, with a default ceiling set to \f$ \simeq 2^{\pm 4\times 10^{18}}\f$.
 
 <table>
 <caption id="multi_row">Number formats</caption>
 <tr><th> @c numbers/ntypes.h <th> Precision requested       <th> Format  <th> Numerical limit
 <tr><td rowspan="2"> MACHINE_LOW_PREC <td> up to 24 <td> FP32, @c float <td> about \f$ 10^{-38}\f$ to \f$ 10^{+ 38}\f$
 <tr> <td> 25 and up <td> MPFR <td>  practically unlimited
 <tr><td rowspan="2"> no flag <td> up to 53 <td> FP64, @c double <td> about \f$ 10^{- 308}\f$ to \f$ 10^{+ 308}\f$
 <tr> <td> 54 and up <td> MPFR <td> practically unlimited
 <tr><td rowspan="2"> MACHINE_EXTRA_PREC <td> up to 64 <td> FP80, @c long @c double <td> about \f$ 10^{- 4930}\f$ to \f$ 10^{+ 4930}\f$
 <tr> <td> 65 and up <td> MPFR <td> practically unlimited
 </table>
 

 \warning
 If a number either read as input from a file or computed and destined to the output cannot be represented with the selected precision (@c inf or @c NaN exception),
 a warning is issued and the operation fails. Usually, the solution consists in using a higher precision (i.e. MPFR numbers).
 Anoter possible issue (if you are running Newton steps) is that you have entered the immediate neigborhood of a critical point,
 which induces a division by zero or almost zero. Try neutralizing the offending point.
 
 Other numerical settings can be adjusted in ntypes.h, in particular :
 - @c HUGE_DEGREES sets the highest polynomial degree that can be handled (default \f$ 2^{64}-2\f$),
 - @c HUGE_MP_EXP sets the highest exponent that MPFR numbers can handle (default \f$ \pm 4 \times 10^{18} \f$).
 
 \subsection sectInputOutput Input and outputs
 
 @c FastPolyEval reads its input data from files and produces its output as file. The file format is @c CSV.
 To be valid, the rest of the @c CSV must contains a listing of complex numbers, written
 \code
 a, b
 \endcode
 where @c a is the real part and @c b is the imaginary part, in decimal form.
 
 Depending on the context, the file will be interpreted as a set of evaluation points, a set of values, or the coefficients of a polynomial.
 For example, here is a valid file:
 \code
 // An approximation of pi
 3.14, 0
 // A complex number close to the real axis
 1.89755593218329076e+99, 2.35849881747969046e-427
 \endcode
 Polynomials are written with the lowest degree first, so \f$ P(z) = 1-2 i z\f$ is represented by the file
 \code
 1, 0
 0, -2
 \endcode
 Comments line are allowed: they start either with @c #, @c // or @c ;. Comment lines will be ignored. Be mindful of the fact that comments can induce an \b offset between
 the line count of the file and the internal line count of @c FastPolyEval.  Blank lines are not allowed.
 Lines are limited to 10 000 symbols (see \ref array_read), which puts a practical limit to the precision at around 15 000 bits. You can easily edit the code to push this
 limitation if you have the need for it.
 
 \warning As a general rule for the tasks of @c FastPolyEval, if an output file name matches the one of an input file, the operation is \b safe. However, the desired output will overwrite the previous version of the file.

 \note
 For users that want a turnkey solution to high-performance polynomial evaluations where the
 outputs of some computations are fed back as input to others,
 we recommend the use of @c FastPolyEval with data residing in a virtual RAM disk,
 on local solid-state drives or on similar low latency/high bandwidth storage solutions.
 If you have doubts on how to do this, please contact your system administrator.
 
 \subsection sectTasksPoly Tasks for generating and handling polynomials
 
 @c FastPolyEval contains a comprehensive set of tools to generate new polynomials, either from scratch or by performing simple operations on other ones.
  
 - Classical polynomials are generated with the tasks @c -Chebyshev, @c -Legendre, @c -Hermite, @c -Laguerre.
 The family @c -hyperbolic is defined recursively by
 \f[
 p_1(z)=z, \qquad p_{n+1}(z) = p_n(z)^2 + z
 \f]
 and plays a central role in the study of the Mandelbrot set (see e.g. reference [2] in \ref sectRefCopy).
 
 - One can build a polynomial from a predetermined set of roots (listed with multiplicity) by invoquing the task @c -roots.
 
 - One can compute the sum (@c -sum), difference (@c -diff), product (@c -prod) or the derivative (@c -der) of polynomials.
 
 \subsection sectTasksComplex Tasks for generating and handling sets of complex numbers
 
 @c FastPolyEval contains a comprehensive set of tools to generate new sets of real and complex numbers, either from scratch or by performing simple operations on other ones.

 - The real (@c -re) and imaginary (@c -im) parts of a list of complex numbers can be extracted. The results are real and therefore of the form @c x, @c 0 .
 The complex conjugate (sign change of the imaginary part) is computed with @c -conj and the complex exponential with @c -exp.

 - @c FastPolyEval can perform various interactions between valid @c CSV files.
 
     1. The concatenation (@c -cat) of two files writes a copy of the second one after a copy of the first one.
     For example, to rewrite a high-precision output @c file_high_prec.csv with a lower precision @c prec, use
     \code
     FastPolyEval prec file_high_prec.csv /dev/null file_low_prec.csv
     \endcode
 
     2. The @c -join task takes the real parts of two sequences and builds a complex number out of it. The imaginary parts are lost.
     \verbatim
     a, *    (entry k from file_1)
     b, *    (entry k from file_2)
     ----
     a, b    (output k)
     \endverbatim
     If the files have unequal lengh, the operation succeeds but stops once the shortest file runs out of entries. Be mindfull of offsets induced by commented lines (see \ref sectInputOutput).
 
     3. The @c -tensor task computes a tensor product (i.e. complex multiplication line by line).
     \verbatim
     a, b    (entry k from file_1)
     c, d    (entry k from file_2)
     ----
     ac-bd, ad+bc  (output k)
     \endverbatim
     If the files have unequal lengh, the operation succeeds but stops once the shortest file runs out of entries. Be mindfull of offsets induced by commented lines (see \ref sectInputOutput).

     4. The @c -grid task computes the grid product of the real part of the entries. The output file contains nb_lines_file_1 \f$ \times\f$ nb_lines_file_2 lines and lists all the products possibles between
     the real part of entries in the first file with the real part of entries in the second file. The imaginary parts are ignored.
 
 - Rotations (@c -rot) map the complex number \f$ a + i b\f$ to \f$ a \times e^{i b}\f$.
 For example, to generate complex values from a real profile @c profile.csv and a set of phases @c phases.csv,  you can use
 \code
 FastPolyEval -join prec profile.csv phases.csv rot_tmp.csv
 FastPolyEval -rot prec rot_tmp.csv complex_output.csv
 rm -f rot_tmp.csv
 \endcode
 
 - The @c -unif task writes real numbers in an arithmetic progression whose characteristics are specified by the parameters.

 - The @c -sphere task writes polar coordinates approximating an uniform distribution on the [Riemann sphere](https://en.wikipedia.org/wiki/Riemann_sphere).
 The output value \f$ (a, b)\f$ represents Euler angles in the (vertical, horizontal) convention. It represents the point
 \f[
 (\cos a \cos b, \cos a \sin b, \sin a) \in \mathbb{S}^2 \subset \mathbb{R}^3.
 \f]
 The @c -polar task maps a pair of (vertical, horizontal) angles onto the complex plane by [stereographical projection](https://en.wikipedia.org/wiki/Stereographic_projection).
 To generate complex points that are uniformly distributed on the Riemann sphere, starting with @c SEED points at the equator, you can use the following code that combines those tasks.
 Note that there will be about \f$ SEED^2/\pi\f$ points on the sphere and that the points are computed in a deterministic way (they will be identical from one call to the next).
 A SEED of 178 produces about 10 000 complex points.
 \code
 SEED=178
 FastPolyEval -sphere prec $SEED tmp_file.csv
 FastPolyEval -polar prec tmp_file.csv sphere.csv
 rm -f tmp_file.csv
 \endcode
  
 - @c FastPolyEval  also contains two random generators (actually based on MPFR).
 The @c -rand task produces real random numbers that are uniformly distributed in an interval.
 The @c -normal task produces real random numbers with a [Gaussian distribution](https://en.wikipedia.org/wiki/Normal_distribution) whose characteristics are specified by the parameters.
 Use the @c -join task to merge the outputs of two @c -normal tasks into a [complex valued normal distribution](https://en.wikipedia.org/wiki/Complex_normal_distribution).
 
 
 \subsection sectTasksFPE Tasks using the FPE algorithm for production use and benchmarking

 The main tasks of the @c FastPolyEval library are the following:
     - @c -eval evaluates a polynomial on a set of points using the FPE algorithm.
     - @c -evalD  evaluates the derivative of a polynomial on a set of points using the FPE algorithm.
 
 The precision of the computation and the name of the inputs (polynomial and points) & output (values) files are given as parameters.
 The first optional parameter specifies a second output file that contains  a complementary report on the estimated quality of the evaluation at the precision chosen.
 For each evaluation point, this report contains an upper bound for the evaluation errors (in bits), a conservative estimate on the number of correct bits of the result,
 and the number of terms that where kept by the FPE algorithm.
 
 For benchmark purposes (see [reference [1]] (https://hal.archives-ouvertes.fr/hal-03820369) in Section \ref sectRefCopy), we propose a set of optional parameters
 that specify how many times the preprocessing and the FPE algorithm should be run at each point (to improve the accuracy of the time measurement), whether the Hörner algorithm
 should also be run and, if so, how many times it should run.
 
 The task @c -evalN  evaluates one @b Newton @b step of a polynomial on a set of points, i.e.
 \f[
 \textrm{NS}(P,z) = \frac{P(z)}{P'(z)}.
 \f]
 Recall that [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method) for finding roots of \f$ P\f$ consists in computing
 the sequence defined recursively by
 \f[
 z_{n+1} = z_n -\textrm{NS}(P,z_n).
 \f]
 The @c -iterN task iterates the Newton method  a certain amount of times and with a given precision. It thus realizes a partial search of the roots of the polynomial.
 For best results, we recommend successive calls where the precision is gradually increased and the maximum number of iterations is reduced.
 For guidance in the choice of the starting points, please refer to reference [3] in Section \ref sectRefCopy.
 
 -analyse computes the concave cover and the intervals of |z| for which the evaluation  strategy changes
 
 
 The @c -analyse task computes the concave cover \f$ \hat{E}_P\f$ , the strip \f$ S_\delta\f$ (see Section \ref sectAlgo)
 and the intervals of |z| for which the evaluation strategy of FPE (i.e. the reduced polynomial \f$ Q_p(z)\f$) changes.
 It is intended mostly for an illustrative purpose on low degrees, when the internals of the FPE algorithm can still be checked by hand.
 However, the intervals where a parsimonious representation is valid may also be of practical use
 (see [reference [1]] (https://hal.archives-ouvertes.fr/hal-03820369) in Section \ref sectRefCopy).
 
 \htmlonly <hr> \endhtmlonly

 \section sectCode Notes about the implementation

 For further details, please read the documentation associated to each [<b> source file</b>] (files.html).
  
 Do not hesitate to contact us (see \ref sectRefCopy) to signal bugs and to request new features.
 We are happy to help the development of the community of the users of @c FastPolyEval.
  
 \htmlonly <hr> \endhtmlonly

 \section sectRefCopy References, Contacts and Copyright

 Please cite the reference [1] below if you use or distribute this software. The other citations are in order of appearance in the text:

   - <b>[1] R. Anton, N. Mihalache & F. Vigneron. Fast evaluation of real and complex polynomials. 2022.
        [[hal-03820369] (https://hal.archives-ouvertes.fr/hal-03820369)].</b>
 
   - [2] N. Mihalache & F. Vigneron. How to split a tera-polynomial. In preparation.
 
   - [3] J.H. Hubbard, D. Schleicher, S. Sutherland.  How to find all roots of complex polynomials by Newton's method. Invent. math., 146:1-33, 2001.
      [[Llink](http://pi.math.cornell.edu/~hubbard/NewtonInventiones.pdf)] on author's page.

  @authors Nicolae Mihalache -- Université Paris-Est Créteil (nicolae.mihalache@u-pec.fr)\n
 François Vigneron -- Université de Reims Champagne-Ardenne (francois.vigneron@univ-reims.fr)
  
  \copyright This software is released under the BSD licence, with an attribution clause (see \ref sectLicence).\n
            Copyright 2022 -- Université Paris-Est Créteil & Université de Reims Champagne-Ardenne.\n
 
 \section sectThanks Thanks
 It is our pleasure to thank our families who supported us in the long and, at times stressful, process that gave birth to the @c FastPolyEval library,
 in particular Ramona who co-authored reference [1], Sarah for a thorough proof-reading of [1] and Anna who suggested the lightning pattern for our logo.
 
 We also thank [Romeo](https://romeo.univ-reims.fr), the HPC center of the University of Reims Champagne-Ardenne, on which
 we where able to benchmark our code graciously.
 
 \image html romeo.jpg

 \section sectLicence License

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.

 3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

 4. Redistributions of any form whatsoever must retain the following
    acknowledgment:
    
   'This product includes software developed by Nicolae Mihalache
    & François Vigneron at Univ. Paris-Est Créteil & Univ. de Reims
    Champagne-Ardenne. Please cite the following reference if you
    use or distribute this software:
    R. Anton, N. Mihalache & F. Vigneron.
    Fast evaluation of real and complex polynomials. 2022.
    [[arxiv-2211.06320](https://arxiv.org/abs/2211.06320)] or
    [[hal-03820369](https://hal.archives-ouvertes.fr/hal-03820369)]'

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 */

/**
 \dir numbers
 \brief This folder contains type definitions, data structures and basic functions to compute with complex numbers.
 
 FastPolyEval can handle two types of numbers :
  - machine floating point numbers:  [FP32] (https://en.wikipedia.org/wiki/Single-precision_floating-point_format),
 [FP64] (https://en.wikipedia.org/wiki/Double-precision_floating-point_format) and
 [FP80] (https://en.wikipedia.org/wiki/Extended_precision),
  - arbitary precision numbers that are  based on [MPFR] (https://www.mpfr.org).
 */

/**
 \dir poly
 \brief This folder contains type definitions, data structures and basic functions to compute with complex polynomials.
 
 There are four types of polynomials : real or complex, and depending of the type of numbers used, machine floating point numbers,
 or arbitary precision numbers that are based on [mpfr] (https://www.mpfr.org).
*/

/**
 \dir eval
 \brief This folder contains type definitions, data structures and basic functions for fast evaluation of complex polynomials.
 
 There are two types of evaluators, depending of the type of coefficients of polynomials, machine floating point numbers,
 or arbitary precision numbers that are based on [mpfr] (https://www.mpfr.org).
*/

/**
 \dir tools
 \brief This folder contains some tools, like arrays to do IO and conversions of polynomials, time measuring tools, etc.
 
 There are two types of arrays, depending of the type of coefficients of polynomials, machine floating point numbers,
 or arbitary precision numbers that are based on [mpfr] (https://www.mpfr.org).
*/

/**
 \dir apps
 \brief This folder contains the production apps.
 
 The role of the main app is to parse the parameters and call the appropriate specialized app.
 
 The specialized apps are classified in two families ...
*/

/**
 \file main.h
 \brief This is the entry point in the main app.
 
 Depending on the first command line parameter,
 it launches an internal app or presents a help screen to guide the user.
*/


#ifndef main_h
#define main_h

#define APP_OK         0
#define APP_ERROR      1
#define APP_PARAMS     2

/// The signature of a main function of an app in this project with machine numbers.
typedef int (* mainFuncf)(int argv, const char **args);

/// The signature of a main function of an app in this project with arbitrary precision.
typedef int (* mainFunc)(long prec, int argv, const char **args);

#endif /* main_h */
