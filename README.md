# odd-dim-dirichlet-boxes-threshold

Code accompanying Section 5 of
Fabio Francesconi and Julie Rowlett, ..., arXiv: ...

 The program computes certified enclosures of m(23), m(24) and m(25), defined in
eq. <(5.2)> of the paper, using the ball arithmetic of Arb (part of FLINT):
interval Riemann sums on [1, 20] with step 2^-10, theta series truncated at
K = 6 with explicit remainder bounds, and the tail bound of eq. <(5.5)>.

## Requirements

FLINT >= 3.0 and a C++ compiler. On Ubuntu, Debian or Linux Mint:

    sudo apt install build-essential libflint-dev

## Usage

    g++ -O2 threshold.cpp -o threshold -lflint
    ./threshold        # step h = 2^-10;  ./threshold 12  uses h = 2^-12

## Output

Obtained with FLINT 3.0.1, g++ 13.3, Linux Mint 22.2 (running time < 1 s):

    ...
    
Each line gives a ball [midpoint +/- radius] that is guaranteed to contain the
exact value. With other versions of FLINT the last digits of the midpoints may
change; the enclosures and the conclusions do not.

## License

MIT
