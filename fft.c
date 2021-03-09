/*****************************************************************************
 * FFT for HP Agilent 16700 series logic analyzers with 16533A, 16534A
 * oscilloscope module.
 * Andrei A. Semenov <rtr989@gmail.com>
 *
 * See the LICENSE file accompying this source file for copyright and
 * redistribution information. Copyright (c) 2020 by Andrei A. Semenov.
 *
 * For use with the Tool Development Kit.
 *
 *****************************************************************************/

const int size = 32768;
double lagrange(double *x_values, double *y_values, int size, double x);

void execute(TDKDataGroup &dg, TDKBaseIO &io) {
  TDKDataSet ds;
  TDKDataSet fftDS;
  TDKLabelEntry fft;
  TDKLabelEntry le;
  TDKLabelEntry freqTXT;
  TDKLabelEntry freqINT;

  const double TwoPi = 6.283185307179586;
  double fftA1[size];
  float frequency[100];
  float db[100];
  double max, min_max, value, min, min_min;
  long long max_place[100];
  long long state, first_pos;
  char freqVal[100];
  char dbVal[100];

  int c, a, i, k, err, j, n, m, Mmax, Istp, threshold;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  double *Tmvl;
  String f;
  int range = size;
  float baseFreq;
  double xVal[30];
  double yVal[30];
  float refV;
  

  for (i = 0; i < 100; i++) {
    frequency[i] = 0;
    max_place[i] = 0;
    db[i] = 0;
  }

  String channel = io.getArg(0);

  err = sscanf(io.getArg(1), "%i", &threshold);
  if (err != 1) {
    io.print("Unable to convert Threshold parameter");
    return;
  }

  err = sscanf(io.getArg(2), "%f", &refV);
  if (err != 1) {
    io.print("Unable to convert Threshold parameter");
    return;
  }
  if(refV==0) refV=1;
  
  err = ds.attach(dg);
  if (err) {
    io.printError(err);
    return;
  }

  first_pos = ds.getPosition();

  err = ds.setTimeBias();
  if (err) {
    io.printError(err);
    return;
  }

  long long Tsample = (ds.lastPosition() - ds.firstPosition()) / (range - 1);

  fftDS.createTimePeriodic(dg, "fftDataSet", 16384, 0, ds.getCorrelationTime(),
                           Tsample);

  err = fftDS.setStateBias();

  baseFreq = 1000000000000 / (size * 1.0 * Tsample);

  err = le.attach(ds, channel);
  if (err) {
    io.print("Can not connect to Oscilloscope channel");
    io.printError(err);
    return;
  }

  err = freqTXT.createTextData(fftDS, "Frequency", 20);
  if (err) {
    io.printError(err);
    return;
  }

  ////////Get FFT data////////////////////////////////

  n = range * 2;
  Tmvl = new double[n];

  le.setPosition(le.firstPosition());
  for (i = 0; i < n; i += 2) {
    err = le.next(value);
    Tmvl[i] = 0;
    Tmvl[i + 1] = value;
  }

  i = 1;
  j = 1;
  while (i < n) {
    if (j > i) {
      Tmpr = Tmvl[i];
      Tmvl[i] = Tmvl[j];
      Tmvl[j] = Tmpr;
      Tmpr = Tmvl[i + 1];
      Tmvl[i + 1] = Tmvl[j + 1];
      Tmvl[j + 1] = Tmpr;
    }
    i = i + 2;
    m = range;
    while ((m >= 2) && (j > m)) {
      j = j - m;
      m = m >> 1;
    }
    j = j + m;
  }

  Mmax = 2;
  while (n > Mmax) {
    Theta = -TwoPi / Mmax;
    Wpi = sin(Theta);
    Wtmp = sin(Theta / 2);
    Wpr = Wtmp * Wtmp * 2;
    Istp = Mmax * 2;
    Wr = 1;
    Wi = 0;
    m = 1;

    while (m < Mmax) {
      i = m;
      m = m + 2;
      Tmpr = Wr;
      Tmpi = Wi;
      Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
      Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

      while (i < n) {
        j = i + Mmax;
        Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j - 1];
        Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j - 1];

        Tmvl[j] = Tmvl[i] - Tmpr;
        Tmvl[j - 1] = Tmvl[i - 1] - Tmpi;
        Tmvl[i] = Tmvl[i] + Tmpr;
        Tmvl[i - 1] = Tmvl[i - 1] + Tmpi;
        i = i + Istp;
      }
    }

    Mmax = Istp;
  }

  for (i = 0; i < range; i++) {
    j = i * 2;
    fftA1[i] = 2 * sqrt(pow(Tmvl[j], 2) + pow(Tmvl[j + 1], 2)) / range;
  }

  delete[] Tmvl;

  ////////Signal threshold calculation/////////////////
  int count = 0;
  max = 0;
  min = 1000;
  for (i = 0; i < 16384; i++) {

    if (fftA1[i] > max) {
      max = fftA1[i];
    }
    /*
    min_min += fftA1[i];
    count++;

    if (count == 9) {
      min_min = min_min / 10;
      count = 0;
      if (min_min < min){
        min = min_min;
      }
      min_min=0;
    }
*/

  }
  min_max = max / 100 * threshold;
  //io.printf("min: %f", min);



  ////////Put FFT on Waveform///////////////////////////

  err = fft.createAnalogData(fftDS, "FFT", max, max * 2.5);

  for (i = 0; i < range / 2; i++) {
    fft.replaceNext(fftA1[i] * 2);
  }

  err = freqINT.createIntegralData(fftDS, "Frequency_#", 8);

  fft.setPosition(fft.firstPosition());
  while (fft.next(value)) {
    freqINT.replaceNext((unsigned int)(0));
  }

  ////////Highlight all signals/////////////////////////
  k = 0;
  for (i = 1; i < 16383; i++) {

    if (k < 100 && (fftA1[i] > min_max) && (fftA1[i - 1] < fftA1[i]) &&
        (fftA1[i + 1] < fftA1[i])) {

      max_place[k] = i + 1;
      frequency[k] = (i + 1) * baseFreq;

      c = i - 15;
      for (a = 0; a < 30; a++) {
        xVal[a] = a;
        yVal[a] = fftA1[c];
        c++;
      }

      float o = 15, p;
      double maxLagr = 0;

      if (fftA1[i - 1] > fftA1[i + 1]) {

        for (a = 0; a < 7; a++) {

          value = lagrange(xVal, yVal, 30, o);
          if (value > maxLagr) {
            maxLagr = value;
            p = a;
          }
          o -= 0.1;
        }
        if (p != 0)
          p = p + 2;
        frequency[k] = ((i + 1) * baseFreq) - (baseFreq / 10 * p);
        
         
      }

      if (fftA1[i + 1] > fftA1[i - 1]) {

        for (a = 0; a < 7; a++) {

          value = lagrange(xVal, yVal, 30, o);
          if (value > maxLagr) {
            maxLagr = value;
            p = a;
          }
          o += 0.1;
        }
        if (p != 0)
          p = p + 2;
        frequency[k] = ((i + 1) * baseFreq) + (baseFreq / 10 * p);
      }
      db[k]= 20 * log10(fftA1[i]/refV);
      k++;
    }
  }


  ////////Frequency calculation/////////////////////////

  k = 1;

  if (baseFreq / 10 < 1000) {
    io.printf("Measurement accuracy: %.3f Hz", baseFreq / 10);
    if (baseFreq * 16384 < 1000)
      io.printf("Minimum Frequency: %.3f Hz, Max Frequency: %.3f Hz", baseFreq,
                baseFreq * 16384);
    else
      io.printf("Minimum Frequency: %.3f Hz, Max Frequency: %.3f kHz", baseFreq,
                baseFreq * 16384 / 1000);

  } else {
    io.printf("Measurement accuracy: %.3f kHz", baseFreq / 10000);
    io.printf("Minimum Frequency: %.3f kHz, Max Frequency: %.3f MHz",
              baseFreq / 1000, baseFreq * 16384 / 1000000);
  }
  io.printf("--------------------------------------------");

  for (i = 0; i < 100; i++) {
    if (max_place[i] != 0) {

      if (frequency[i] < 100) {
        sprintf(freqVal, "#%i: %.4f Hz", k, frequency[i]);
      }

      if (frequency[i] >= 100 && frequency[i] < 1000) {
        sprintf(freqVal, "#%i: %.2f Hz", k, frequency[i]);
      }

      if (frequency[i] >= 1000 && frequency[i] < 10000) {
        sprintf(freqVal, "#%i: %.4f kHz", k, frequency[i] / 1000);
      }

      if (frequency[i] >= 10000 && frequency[i] < 100000) {
        sprintf(freqVal, "#%i: %.3f kHz", k, frequency[i] / 1000);
      }

      if (frequency[i] >= 100000 && frequency[i] < 1000000) {
        ;
        sprintf(freqVal, "#%i: %.2f kHz", k, frequency[i] / 1000);
      }

      if (frequency[i] >= 1000000 && frequency[i] < 10000000) {
        ;
        sprintf(freqVal, "#%i: %.4f MHz", k, frequency[i] / 1000000);
      }

      if (frequency[i] >= 10000000 && frequency[i] < 100000000) {
        ;
        sprintf(freqVal, "#%i: %.3f MHz", k, frequency[i] / 1000000);
      }

      if (frequency[i] >= 100000000) {
        ;
        sprintf(freqVal, "#%i: %.4f MHz", k, frequency[i] / 1000000);
      }
      sprintf(dbVal, "  %.2f dB", db[i]);
      f="";
      f += freqVal;
      f += dbVal;

      state = max_place[i] - 1;
      freqTXT.setPosition(state);
      freqTXT.replaceNext(f);
      freqTXT.setHighlight(state);

      freqINT.setPosition(state);
      freqINT.replaceNext((unsigned int)(k));
      freqINT.setHighlight(state);

      io.printf(f);
      k++;
    }
  }

  /////////////////////////////////////////////////////
  ds.removeLabelEntry(le);
  dg.removeDataSet(ds);
  delete[] frequency;
  delete[] max_place;
}

StringList getLabelNames() {
  StringList labels;
  labels.put("Input channel: ");
  labels.put("Signal threshold % of the highest: ");
  labels.put("dB Reference voltage Volt: ");
  return labels;
}

// Assign default values to runtime arguments
StringList getDefaultArgs() {
  StringList defs;
  defs.put("Channel C1");
  defs.put("3");
  defs.put("1");
  return defs;
}

double lagrange(double *x_values, double *y_values, int size, double x) {
  double lagrange_pol = 0;
  double basics_pol;

  for (int i = 0; i < size; i++) {
    basics_pol = 1;
    for (int j = 0; j < size; j++) {
      if (j == i)
        continue;
      basics_pol *= (x - x_values[j]) / (x_values[i] - x_values[j]);
    }
    lagrange_pol += basics_pol * y_values[i];
  }
  return lagrange_pol;
}
