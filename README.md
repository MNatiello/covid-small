# covid-small
Markov Jump Process simulating an epidemic outbreak
Programme covid-small -- version 11 -- 2020-11-3
Language: C
Files: Makefile, corona.c, corona.h, README.txt

Uses: dSFMT
The code uses random numbers generated by SIMD-oriented Fast Mersenne Twister
Mutsuo Saito (Hiroshima University) and Makoto Matsumoto (Hiroshima University).
The code is available from
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/dSFMT-src-2.2.3.tar.gz

Implements Kendall's algorithm for a Markov Jump Process inspired in Feller's se
tup. See e.g.,
* Solari, Hernán; Natiello, Mario (2014). Linear Processes in Stochastic Populat
ion Dynamics: Theory and Application to Insect Development.
   Scientific World Journal.
* Natiello, Mario A.; Solari, Hernán G. (2020). Modelling population dynamics ba
sed on experimental trials with genetically modified (RIDL)
    mosquitoes. Ecological Modelling, 424,
* Solari, Hernán; Natiello, Mario (2020) Stochastic model for COVID-19 in slums:
 interaction between biology and public policies. Preprint.

Time-unit is days
Input variables:

NAME      Name of run                               char
Stat      Number of realisations                    positive integer
Duration  Simulation length in days                 positive integer
ext       Contagion rate external to the system     positive double
betaT     Internal contagion rate, class T          positive double
betaU     Internal contagion rate, class U          positive double
detT      Removal intensity, class T                positive double
detU      Removal intensity, class U                positive double
delayT    First day of removal effort, class T      positive integer
delayU    First day of removal effort, class U      positive integer
S         Initial number of susceptibles            positive integer
T         Initial number of T-individuals           nonnegative integer
U         Initial number of U-individuals           nonnegative integer
seed      For random number generator               positive integer
wipe      If non-zero, delete short-lived epidemies non-negative integer
Tmin      upper limit for short-lived epidemies     positive integer

Usage (command line):

./corona NAME Stat Duration ext betaT betaU detT detU delayT delayU S T U seed wipe Tmin
e.g.,
./corona Test 100 300 0.05 2.5 2.5 1.7 1.0 1 3 998 2 0 1977 19 25 

Optional:
If file IConditions exists  reads date and number of individuals in each compartment for each realisation

Output:
Daily record of population in each compartment in file NAMEnn
Final date and population distribution in each compartment for each realisation in file FConditions

Example of continuation run:

../corona Test 100   50 0.05 2.5 2.5 1.7 1.0 1 3 998 2 0 1977 10 19 
echo  "continuing the first 50-days run up to 300 days..."
cp FConditions IConditions
../corona Test 100 300 0.05 2.0 2.0 1.9 1.5 1 1 998 2 0 1997 10 19


