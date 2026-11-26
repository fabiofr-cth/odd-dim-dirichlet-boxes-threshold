# Certified computation of the odd dimension threshold  

Code accompanying Section 5 of
Fabio Francesconi & Julie Rowlett, *<title>*, arXiv:<number>.

The program computes certified enclosures of M(23), M(24) and M(25), defined in
eq. (66) of the paper, using the ball arithmetic of Arb (part of FLINT):
interval Riemann sums on [1, 20] with step 2^-10, theta series truncated at
K = 6 with explicit remainder bounds, and the tail bound of eq. (69).

## Requirements

FLINT >= 3.0 and a C++ compiler. On Ubuntu, Debian or Linux Mint:

    sudo apt install build-essential libflint-dev

## Usage

    g++ -O2 threshold.cpp -o threshold -lflint
    ./threshold        # step h = 2^-10;  ./threshold 12  uses h = 2^-12

## Output

Obtained with FLINT 3.0.1, g++ 13.3, Linux Mint 22.2 (running time < 1 s):

    FLINT 3.0.1, 128 bits, K = 6, a = 20, h = 2^-10

    m( 3) in 0.1349603193 +/- 0.00048956   < 1
    m(23) in 0.9113976983 +/- 0.0027996   < 1
    m(24) in 1.485274531 +/- 0.0045625   > 1
    m(25) in 2.476626428 +/- 0.0076076   > 1

Each line gives a ball [midpoint +/- radius] that is guaranteed to contain the
exact value. With other versions of FLINT the last digits of the midpoints may
change; the enclosures and the conclusions do not.

## License

MIT
