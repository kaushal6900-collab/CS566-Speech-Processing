


#include "stdafx.h"
#include <iostream>
#include <Windows.h>
#include <stdio.h>    
#include <process.h>  

using namespace std;

#pragma comment(lib, "winmm.lib")


// A global array to store the recorded audio samples.
// It's sized for 10 seconds of audio at a 16025 Hz sample rate.
// short int is used because samples are 16-bit.
short int waveIn[16025 * 10]; 

const long long tnrg = 10000; //threshold enrgy which will help to differentiate between noise and speech 


void StartRecord();

int _tmain(int argc, _TCHAR* argv[])
{
    
    StartRecord();

    int x, i = 0; 

    int a[60001] = {0}; //storage for first word
    int b[60001] = {0}; //storage for second word
   
    int p = 0, q = 0; 

    
    int c = 0; /// to track which word we are recording 

    int silenceCounter = 0; 

    const int sampleRate = 16025;
    const int silenceTh = sampleRate * 0.2;/// threshold for scielence counter after which second word will start recording

    
    const int totalS = 16025 * 10; 

    printf("\nProcessing recorded audio...\n");

    
    for (i = 0; i < totalS; i++) {
        x = waveIn[i]; // Get the current audio sample.

        
        long long energy = (long long)x * x;/// energy calculation of every sample

        if (c == 0) { // first word.
            if (energy > tnrg) { // checking for valid speech
                a[p] = x; 
                p++; 
                silenceCounter = 0; // it stay 0 until we are speeking
            } else {
                if (p > 0) { 
                    silenceCounter++; // low energy means silence or noise
                    if (silenceCounter > silenceTh) { 
                        c = 1; //end of 1st word
                        silenceCounter = 0; // reset to zero
                    }
                }
            }
        } else if (c == 1) { //  second word.
            if (energy > tnrg) { // checking for valid speech.
                b[q] = x; 
                q++; 
                silenceCounter = 0; // it stay 0 until we are speeking.
            } else { 
                if (q > 0) { 
                    silenceCounter++; // low energy means silence or noise.
                    if (silenceCounter > silenceTh) { //end of speech
                        c = 2; 
                    }
                }
            }
        } else if (c == 2) { // both words recorded.
            break; 
        }
    }

    
    int end1 = p; 
    int end2 = q; 

    
    int zcr1 = 0, zcr2 = 0;

   
    for (int l = 0; l < end1 - 1; l++) {
       
        if ((a[l] >= 0 && a[l + 1] < 0)) zcr1++;
        else if((a[l + 1] >= 0 && a[l] < 0)) {
            zcr1++;
        }
    }

   
    for (int l = 0; l < end2 - 1; l++) {
        if ((b[l] >= 0 && b[l + 1] < 0)) zcr2++;
        else if((b[l + 1] >= 0 && b[l] < 0)) {
            zcr2++;
        }
    }

   
    printf("ZCR1- %d\n", zcr1);
    printf("ZCR2- %d\n", zcr2);

    
    if (zcr1 > zcr2) {
        printf("First word is 'yes'\n");
    } else if(zcr1 < zcr2) {
        printf("Second word is 'yes'\n");
    }
	else printf("no word spoken \n");

    return 0; 
}



void StartRecord()
{
    const int NUMPTS = 16025 * 10; 
    int sampleRate = 16025; 
    
    HWAVEIN hWaveIn; 
    WAVEFORMATEX pFormat; 
    pFormat.wFormatTag = WAVE_FORMAT_PCM;     
    pFormat.nChannels = 1;                     
    pFormat.nSamplesPerSec = sampleRate;       
    pFormat.nAvgBytesPerSec = sampleRate * 2;  
    pFormat.nBlockAlign = 2;                   
    pFormat.wBitsPerSample = 16;               
    pFormat.cbSize = 0;                        

   
    waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);

    WAVEHDR WaveInHdr; 
    
   
    WaveInHdr.lpData = (LPSTR)waveIn;

    WaveInHdr.dwBufferLength = NUMPTS * 2;

    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
  
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    
    waveInStart(hWaveIn);

    cout << "start speaking" << endl;
   
    Sleep(10 * 1000); 

  
    waveInClose(hWaveIn);
    cout << "Recording finished." << endl;
}
