#ifndef _PREDICTOR_H
#define _PREDICTOR_H

//////////////////////////// Predictor /////////////////////////

// Adaptive order-0 model, derived from Ilia Muraviev's fpaq0p

class Predictor {
private:
  int cxt; // Context: last 0-8 bits with a leading 1
  unsigned int t[512]; // Probability of 1

public:
  Predictor(): cxt(1) {
    for (int i=0; i<512; i++) t[i]=32768;
  }

  // Assume a stationary order 0 stream of 9-bit symbols
  int p() const {
    return t[cxt]>>4;
  }

  void update(int y) {
    if (y) t[cxt]+=65536-t[cxt]>>5;
    else t[cxt]-=t[cxt]>>5;
    if ((cxt+=cxt+y)>=512) cxt=1;
  }
};

#endif
