/* fpaqc - Stationary order 0 file compressor.
(C) 2007, Matt Mahoney under GPL ver. 3.
Released Dec. 24, 2007.

    LICENSE

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 3 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details at
    Visit <http://www.gnu.org/copyleft/gpl.html>.

This program is a simple order-0 model implementing Jarek Duda's asymmetric
binary coder.  The coder uses the same interface as the arithmetic coder
in better compresssors such as PAQ, LPAQ, BBB, or SR2, so it could be
used in those programs as well.

To compress:   fpaqc c input output
To decompress: fpaqc d input output

To compile:

  g++ -O3 -march=pentiumpro -fomit-frame-pointer -s fpaqc.cpp -o fpaqc

fpaqc is an optimization of fpaqb as suggested by Piotr Tarsa.  The format
differs from fpaqb in that the block size (n) is not stored.

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include "my_predictor.h"

// Create an array p of n elements of type T
template <class T> void alloc(T*&p, int n) {
  p=(T*)calloc(n, sizeof(T));
  if (!p) printf("Out of memory"), exit(1);
}

/* This version of the encoder (if compiled without -DARITH) implements
Jarek Duda's asymmetric binary system, as described in

Asymmetric binary encoding is described in chapter 3 of Jarek Duda,
"Optimal encoding on discrete lattice with translational invariant constrains
using statistical algorithms", 2007.
http://uk.arxiv.org/PS_cache/arxiv/pdf/0710/0710.3861v1.pdf

The interface is the same as the arithmetic version (see above).

The coding method has an advantage over arithmetic coding in that it
uses only one state variable instead of two (range bounds), so the
coder can be implemented with lookup tables (although it is not).  
However, the coding method as described is unsuitable for predictive 
methods (PPM, PAQ, LZMA, etc) because the bit stream must be decoded in 
the reverse order that it was encoded, which means that the model would 
see different contexts during compression and decompression.  This was not 
a problem in the paper because the application described, the bit predictions 
were fixed in advance.

This implementation solves the general problem by saving predictions in
a stack and then coding them in reverse order to a second stack, which
is then popped and written to disk in blocks.  Decoding is straightforward
and fast.  The asymmetric coder can be used as a drop-in replacement for
an arithmetic coder.

Asymmetric binary encoding encodes bit d in {0,1}, having probability q
in (0,1), as an integer state x.  If w is the previous state, then d is
encoded:

  if d is 0 then x = ceil((w + 1)/(1 - q)) - 1
  if d is 1 then x = floor(w/q)

To decode, given state x and prediction q:

  d = ceil((x + 1)*q) - ceil(x*q)

which is 1 if the fractional part of x*q >= 1 - q, else 0.  Then the prior
state w is obtained from x:

  if d is 0 then w = x - ceil(x*q)
  if d is 1 then w = ceil(x*q)

The operations may be carried out in N bit arithmetic.  Let Q = pow(2,N).
Then we restrict q to the range [1..Q-1] and x to the range [Q..256*Q-1].

Define x = encode(w,q,d) as:

  if d is 0 then x = ((w + 1)*Q - 1)/(Q - q)
  if d is 1 then x = (w*Q)/q

To keep x in the desired range, we write the low byte (or bytes)
of w first as needed, then update w := x.

Define w,d = decode(x,q) as follows:

  d = floor(((x+1)*q - 1)/Q) - floor((x*q - 1)/Q)
  xq = floor((x*q - 1)/Q) + 1
  if d = 0 then w := x - xq
  if d = 1 then w := xq

If w is not in the desired range, then we read into the low bytes
of w until it is, then update x := w.

Both reading and writing during compression must be in the reverse
order as decompression.  During compression, the input is divided into
blocks of size B bits.  When d is encoded with probability q, the
encoder just saves d and q on a stack of size B until full.  When the
stack is ready to be written, either because it is full or the end
of input is reached, then x is reset to state x = Q, the d and q
are popped and encoded to a second stack of bytes, then the final state
x is written to the output file, and then the bytes are popped and written.

The decoder does not use any stacks.  To decode, the decoder reads x,
then decodes B bits of output.  As it counts output bits down
to 0, it expects to find x = Q, then reads a new value of x as a 3 byte
big-endian integer.

Compression algorithm:

  read d
  push q,d on input stack
  if input stack is full or end of input then
    n := input stack size
    w := Q
    while input stack is not empty
      pop q,d
      while (x = encode(q,d,w)) >= 256*Q
        push and shift out low byte of w to output stack
        x := floor(x/256)
      w := x
    write w
    pop and write output stack

  q is input as a 12 bit number in the range [1..4095].  In the input
  stack, d is appended to the low bit, i.e. push(q*2+d).  Writing the
  output stack is a simple block write.

Decompression:

  n := 0
  while not end of output
    if n = 0 then read x
    x,d := decode(x,q)
    write d
    while x < Q
      read and shift into low byte of x
    n:= n - 1

As an optimization, the encoding function replaces division by q with
a 64-bit scaled multiplication by 1/q from a lookup table.
*/

typedef enum {COMPRESS, DECOMPRESS} Mode;
class Encoder {
private:
  Predictor predictor;
  const Mode mode;     // Compress or decompress?
  enum {B=0x80000};    // Input stack size, max 0xffffff
  enum {BO=B/4+10};    // Output stack size in bytes, >B/8 to allow expansion
  enum {N=12};         // Number of bits in q and min number in x
  enum {FRAC_BITS=32}; // Number of fractional bits
  FILE* archive;       // Compressed data file
  unsigned int x;      // Encoder state
  int n;               // Number of bits in input stack/block
  unsigned short *ins; // Input stack, q*2+d (compression)
  unsigned char *outs; // Output stack (compression)
  unsigned long long qinv[2<<N]; // lookup table of inverses of q
public:
  Encoder(Mode m, FILE* f);
  int decode();        // Uncompress and return bit y
  void encode(int y);  // Compress bit y
  void flush();        // Call when done compressing
  void disp_stack(unsigned char *p, unsigned int w, int q);

  template<typename T>
  bool pop_instack(T *x) {
    if (!n) {
      return false;
    }
    *x = ins[--n];
    return true;
  }

  void push_instack(unsigned short x) {
    ins[n++] = x;
  }
};

#define L (1<<N)
#define LOWBIT(x) ((x) & 1)
#define b 256
#define PUSH_OUTSTACK(stack) *(--stack)

// Initialize
Encoder::Encoder(Mode m, FILE* f): mode(m), archive(f), x(1<<N), 
    n(mode==DECOMPRESS), ins(0), outs(0) {
  if (mode==COMPRESS) {
    alloc(ins, B);
    alloc(outs, BO);
#define DEST_STATE(currstate, q) ((currstate + LOWBIT(~q)) * qinv[q] - 1 >> FRAC_BITS)
    for (int i=1; i<1<<N; ++i) {
      qinv[i*2+1]=qinv[(2<<N)-2*i]=(1ull<<FRAC_BITS+N)/i;
      ++qinv[i*2+1];
    }

    // look for errors
    for (int i=1; i<1<<N; ++i) {
      if(i<=20) {
        int q = i*2+1;
        int qs;
        for(qs=q; qs<=L; qs<<=1);
        int dest1 = DEST_STATE(qs, q);
        int dest2 = DEST_STATE(qs+2, q);
        printf("dest1 = %6x\n", dest1);
        printf("dest2 = %6x\n", dest2);
        putchar('\n');
      }
    }
  }
}

// Return an uncompressed bit
inline int Encoder::decode() {

  // Read x and n from block header
  if (!--n) {
    if (x!=1<<N) {
      printf("Archive error: x=%X at %ld\n", x, ftell(archive));
      exit(1);
    }
    x=getc(archive);
    x=x*256+getc(archive);
    x=x*256+getc(archive);
    n=B;
    if (x<1<<N || x>=256<<N) {
      printf("Archive error: x=%X at %ld\n", x, ftell(archive));
      exit(1);
    }
  }

  // Decode
  unsigned int q=predictor.p();
  q+=q<2048;
  assert(q>=1 && q<1<<N);
  unsigned int xq=x*q-1;
  int d=(xq&((1<<N)-1))+q>>N;
  assert(d==0 || d==1);
  predictor.update(d);
  xq=(xq>>N)+1;
  assert(xq>0 && xq<x);
  if (d) x=xq;
  else x-=xq;
  while (x<1<<N)
    x=(x<<8)+getc(archive);
  assert(x>=1<<N && x<256<<N);
  return d;
}

// Compress bit d
inline void Encoder::encode(int d) {
  if (n>=B) flush();  // sets n=0;
  int q=predictor.p();
  q+=q<2048;
  ins[n++]=q*2+d;
  predictor.update(d);
}

#define blue "\e[1;34m"
#define green "\e[1;32m"
#define red "\e[1;31m"
#define bold "\e[1m"
#define normal "\e[0m"
#define bkgrd_magenta "\e[1;45m"
#define bkgrd_red "\e[1;41m"
#define bkgrd_blue "\e[1;44m"
#define bkgrd_cyan "\e[1;46m"
#define underscore "\e[4m"

#define DISP_BYTES 50
void Encoder::disp_stack(unsigned char *p, unsigned int w, int q) {
  unsigned char *pend = &outs[BO];
  if(pend-p >= DISP_BYTES) {
    return;
  }
  putchar('\r');
  printf("%*s" bkgrd_blue underscore "%2x" normal bkgrd_blue "%03x" normal, 2*(DISP_BYTES-(int)(pend-p)), "", w >> N, w & (L - 1));
  while(p < pend) {
    printf("%02x", *(p++));
  }
  printf(" * %*llx.%0*llx", N/4, qinv[q]>>FRAC_BITS, FRAC_BITS/4, qinv[q]&((1ull<<FRAC_BITS)-1));
  usleep(0.3 * 1000000);
}

void Encoder::flush() {
  unsigned int w=L;  // state
  unsigned char *p=&outs[BO];  // output stack pointer
  int q;

  while (pop_instack(&q)) {
    assert(q>=2 && q<2<<N);

    disp_stack(p, w, q);

//#define DEST_STATE(currstate, q) ((currstate + LOWBIT(~q)) * qinv[q] - 1 >> FRAC_BITS)
    while ((x=DEST_STATE(w, q)) >= b*L) {
      assert(p>outs+3 && p<=outs+BO);
      PUSH_OUTSTACK(p) = (unsigned char)w;
      w /= b;

      //disp_stack(p, w, q);
    }
    assert(x >= L && x < b*L);
    w = x;
  }

  // Append w to output stack and write
  assert(p>outs+3 && p<=outs+BO);
  PUSH_OUTSTACK(p)=w;
  PUSH_OUTSTACK(p)=w>>8;
  PUSH_OUTSTACK(p)=w>>16;
  fwrite(p, 1, outs+BO-p, archive);
}
#undef b
#undef L

//////////////////////////// main ////////////////////////////

int main(int argc, char** argv) {

  // Check arguments: fpaqc c/d input output
  if (argc!=4 || (argv[1][0]!='c' && argv[1][0]!='d')) {
    printf("fpaqc (C) 2007, Matt Mahoney\n"
      "This is free software under GPL, http://www.gnu.org/copyleft/gpl.html\n\n"
      "To compress:   fpaqc c input output\n"
      "To decompress: fpaqc d input output\n");
    exit(1);
  }

  // Start timer
  clock_t start = clock();

  setvbuf(stdout, NULL, _IONBF, 0);

  // Open files
  FILE *in=fopen(argv[2], "rb");
  if (!in) perror(argv[2]), exit(1);
  FILE *out=fopen(argv[3], "wb");
  if (!out) perror(argv[3]), exit(1);
  int c;

  // Compress
  if (argv[1][0]=='c') {
    Encoder e(COMPRESS, out);
    while ((c=getc(in))!=EOF) {
      e.encode(0);
      for (int i=7; i>=0; --i)
        e.encode((c>>i)&1);
    }
    e.encode(1);  // EOF code
    e.flush();
  }

  // Decompress
  else {
    Encoder e(DECOMPRESS, in);
    while (!e.decode()) {
      int c=1;
      while (c<256)
        c+=c+e.decode();
      putc(c-256, out);
    }
  }
  putchar('\n');

  // Print results
  printf("%s (%ld bytes) -> %s (%ld bytes) in %1.2f s.\n",
    argv[2], ftell(in), argv[3], ftell(out),
    ((double)clock()-start)/CLOCKS_PER_SEC);
  return 0;
}

