
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

const double TwoPi = 6.283185307179586;
void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft);
double A1[32768];
double fftA1[32768];



void execute(TDKDataGroup& dg, TDKBaseIO& io)
{
  TDKDataSet ds;
  TDKLabelEntry fft;
  TDKLabelEntry le;
  TDKLabelEntry freqTXT;
  TDKLabelEntry freqINT;

  double value;
  int i=0;
  int k=0;
  int err;
  int gain, bias;
  int threshold=10;
  int rate = 2000000000;
  float frequency[100];
  double max;
  double min_max;
  int max_place[100];
  char freqVal[100];
  String f;
  int state;
  
  enum { white, white2, scarlet, pumpkin, yellow, lime, turquoise, lavender };
  
  
  for (i = 0; i < 100; i++) {
   frequency[i]=0;
   max_place[i]=0;
  }
  
  String channel = io.getArg(0);
  err = sscanf( io.getArg( 1 ), "%i", &rate );
    if( err != 1 ){
        io.print( "Unable to convert Sample Rate parameter" );
        return;
    }
  err = sscanf( io.getArg( 2 ), "%i", &threshold );
    if( err != 1 ){
        io.print( "Unable to convert Threshold parameter" );
        return;
    }  
  err = sscanf( io.getArg( 3 ), "%i", &gain );
    if( err != 1 ){
        io.print( "Unable to convert Gain parameter" );
        return;
    }  
  err = sscanf( io.getArg( 4 ), "%i", &bias );
    if( err != 1 ){
        io.print( "Unable to convert Bias parameter" );
        return;
    }  
  
  err = ds.attach(dg);
  if( err )
  {
    io.printError( err );
    return;
  }

  err = le.attach(ds, channel );
  if( err )
  {
    io.print( "Can not connect to Oscilloscope channel" );
    io.printError( err );
    return;
  }
  
  i=0;
  while (le.next( value ))
  {
    A1[i] = value;
    i++;
  }
  
  
  err = fft.create( ds, "FFT", le );
  if( err )
  {
    io.printError( err );
    return;
  }
  
  
  err = freqTXT.createTextData( ds, "Frequency", 20 );
  if( err )
  {
    io.printError( err );
    return;
  }
  
  err = freqINT.createIntegralData( ds, "Frequency_#", 8 );
  if( err )
  {
    io.printError( err );
    return;
  }
  
  le.setPosition(le.firstPosition());
  while (le.next( value ))
  {
    freqINT.replaceNext( (unsigned int)(0) );
  }
  
////////Get FFT array////////////////////////////////
  
  FFTAnalysis(A1, fftA1, 32768, 32768);
  
////////Put FFT to Waveform//////////////////////////
delete []A1;


  le.setPosition(le.firstPosition());
  i=0;
  while( le.next( value ) )
  {
    //value = fftA1[i];
    fft.replaceNext( (fftA1[i]*gain)+bias );
    i++;
  }

////////Signal threshold calculation/////////////////

  max=0;
  for (i = 0; i < 16384; i++) {
  
    if(fftA1[i]>max){
    max=fftA1[i];
    }
  }
  min_max = max/100*threshold;
    
////////Highlight all signals/////////////////////////
  for (i = 1; i < 16383; i++) {
  
    if(k<100&&(fftA1[i]>min_max)&&(fftA1[i-1]<fftA1[i])&&(fftA1[i+1]<fftA1[i])){
      max_place[k]=i+1;
      k++;
    }
  }
////////Frequency calculation/////////////////////////
  k=1;
  for (i = 0; i < 100; i++) {
    if(max_place[i]>0){
    
      frequency[i] = max_place[i] * (rate / 32768.0);
      
      if(frequency[i]<100) {
      sprintf ( freqVal, "#%i %.4f Hz", k, frequency[i] );
      }
      
      if(frequency[i]>=100 && frequency[i]<1000) {
      sprintf ( freqVal, "#%i %.2f Hz", k, frequency[i] );
      }
      
      if(frequency[i]>=1000 && frequency[i]<10000) {
      sprintf ( freqVal, "#%i %.4f kHz", k, frequency[i]/1000 );
      }
      
      if(frequency[i]>=10000 && frequency[i]<100000) {
      sprintf ( freqVal, "#%i %.3f kHz", k, frequency[i]/1000 );
      }
      
      if(frequency[i]>=100000 && frequency[i]<1000000) {;
      sprintf ( freqVal, "#%i %.2f kHz", k, frequency[i]/1000 );
      }
      
      if(frequency[i]>=1000000 && frequency[i]<10000000) {;
      sprintf ( freqVal, "#%i %.4f MHz", k, frequency[i]/1000000 );
      }
      
      if(frequency[i]>=10000000 && frequency[i]<100000000) {;
      sprintf ( freqVal, "#%i %.3f MHz", k, frequency[i]/1000000 );
      }
      
      if(frequency[i]>=100000000) {;
      sprintf ( freqVal, "#%i %.2f MHz", k, frequency[i]/1000000 );
      }
      
      
      
      f=freqVal;
      
      state = max_place[i] - 16384;
      freqTXT.setPosition( state );
      freqTXT.setColor( state, yellow );      
      freqTXT.replaceNext( f );
    
      freqINT.setPosition( state );  
      freqINT.setColor( state, yellow );
      freqINT.replaceNext( (unsigned int)(k) );
      
      io.printf( f );
      k++;
    }
  }

/////////////////////////////////////////////////////

delete []frequency;
delete []max_place;
delete []fftA1;

}





void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft) {
  int i, j, n, m, Mmax, Istp;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  double *Tmvl;

  n = Nvl * 2; 
  Tmvl = new double[n];

  for (i = 0; i < n; i+=2) {
   Tmvl[i] = 0;
   Tmvl[i+1] = AVal[i/2];
 }

 i = 1; j = 1;
 while (i < n) {
  if (j > i) {
    Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
    Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
  }
  i = i + 2; m = Nvl;
  while ((m >= 2) && (j > m)) {
    j = j - m; m = m >> 1;
  }
  j = j + m;
}

Mmax = 2;
while (n > Mmax) {
  Theta = -TwoPi / Mmax; Wpi = sin(Theta);
  Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
  Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

  while (m < Mmax) {
    i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
    Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
    Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

    while (i < n) {
      j = i + Mmax;
      Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
      Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

      Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
      Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
      i = i + Istp;
    }
  }

  Mmax = Istp;
}

for (i = 0; i < Nft; i++) {
  j = i * 2; 
  FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;

}

delete []Tmvl;

}


StringList getLabelNames()
{
  StringList labels;
  labels.put("Input channel: ");
  labels.put("Sample Rate Sa/s: ");
  labels.put("Signal threshold % of the highest: ");
  labels.put("Waveform signal gain: ");
  labels.put("Waveform signal bias (0-100): ");
  return labels;
}


// Assign default values to runtime arguments
StringList getDefaultArgs()
{
  StringList defs;
  defs.put("Channel A1");
  defs.put("2000000000");
  defs.put("10");
  defs.put("1");
  defs.put("100");
  return defs;
}
