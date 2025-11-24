#include "stdafx.h"
#include "FeatureExtraction.h"
#include "Codebook.h"
#include "HMM.h"
#include <Windows.h>
#include <mmsystem.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

#pragma comment(lib, "winmm.lib")

using namespace std;

// Global buffer for recording
short int waveIn[16025 * 3];

// Function Prototypes
void StartRecord();
string GenerateFilename(int digit, int utterance);
vector<double> ReadSignal(string filename);

int _tmain(int argc, _TCHAR *argv[])
{
	int option;
	Codebook cb;
	HMM hmms[10];
	bool modelsLoaded = false;

	while (true)
	{
		cout << "\n--- Digit Recognition System ---\n";
		cout << "1. Train\n";
		cout << "2. Test (Offline)\n";
		cout << "3. Live Test\n";
		cout << "4. Exit\n";
		cout << "Enter option: ";
		cin >> option;

		if (option == 1)
		{
			// Train
			vector<vector<double>> allFeatures;
			cout << "Collecting features for Codebook training...\n";

			// 1. Collect all features for Codebook
			for (int d = 0; d <= 9; d++)
			{
				for (int u = 1; u <= 30; u++)
				{ // Utterances 1-30 for training
					string filename = GenerateFilename(d, u);
					vector<double> signal = ReadSignal(filename);
					if (signal.empty())
						continue;

					vector<vector<double>> features;
					ComputeFeatures(signal, features);
					allFeatures.insert(allFeatures.end(), features.begin(), features.end());
				}
			}

			if (allFeatures.empty())
			{
				cout << "No features collected. Check file paths.\n";
				continue;
			}

			cout << "Training Codebook with " << allFeatures.size() << " vectors...\n";
			cb.Train(allFeatures);
			cb.Save("codebook.txt");
			cout << "Codebook saved.\n";

			// 2. Train HMMs
			cout << "Training HMMs...\n";
			for (int d = 0; d <= 9; d++)
			{
				vector<vector<int>> observations;
				for (int u = 1; u <= 30; u++)
				{
					string filename = GenerateFilename(d, u);
					vector<double> signal = ReadSignal(filename);
					if (signal.empty())
						continue;

					vector<vector<double>> features;
					ComputeFeatures(signal, features);
					vector<int> obsSeq;
					for (size_t i = 0; i < features.size(); i++)
					{
						obsSeq.push_back(cb.Quantize(features[i]));
					}
					observations.push_back(obsSeq);
				}

				if (observations.empty())
					continue;

				hmms[d].BaumWelch(observations);
				stringstream ss;
				ss << "model_digit_" << d << ".txt";
				hmms[d].Save(ss.str());
				cout << "HMM for digit " << d << " trained and saved.\n";
			}
			modelsLoaded = true;
		}
		else if (option == 2)
		{
			// Test Offline
			if (!modelsLoaded)
			{
				cb.Load("codebook.txt");
				for (int d = 0; d <= 9; d++)
				{
					stringstream ss;
					ss << "model_digit_" << d << ".txt";
					hmms[d].Load(ss.str());
				}
				modelsLoaded = true;
			}

			int totalTests = 0;
			int correct = 0;

			for (int d = 0; d <= 9; d++)
			{
				for (int u = 31; u <= 40; u++)
				{ // Utterances 31-40 for testing
					string filename = GenerateFilename(d, u);
					vector<double> signal = ReadSignal(filename);
					if (signal.empty())
						continue;

					vector<vector<double>> features;
					ComputeFeatures(signal, features);
					vector<int> obsSeq;
					for (size_t i = 0; i < features.size(); i++)
					{
						obsSeq.push_back(cb.Quantize(features[i]));
					}

					double maxProb = -1e308;
					int recognizedDigit = -1;

					for (int m = 0; m <= 9; m++)
					{
						double prob = hmms[m].GetProbability(obsSeq);
						if (prob > maxProb)
						{
							maxProb = prob;
							recognizedDigit = m;
						}
					}

					cout << "File: " << filename << " -> Recognized: " << recognizedDigit << " (Actual: " << d << ")\n";
					if (recognizedDigit == d)
						correct++;
					totalTests++;
				}
			}
			if (totalTests > 0)
				cout << "Accuracy: " << (double)correct / totalTests * 100.0 << "%\n";
			else
				cout << "No test files found.\n";
		}
		else if (option == 3)
		{
			// Live Test
			if (!modelsLoaded)
			{
				cb.Load("codebook.txt");
				for (int d = 0; d <= 9; d++)
				{
					stringstream ss;
					ss << "model_digit_" << d << ".txt";
					hmms[d].Load(ss.str());
				}
				modelsLoaded = true;
			}

			StartRecord();

			// Silence Removal (Energy based)
			vector<double> speechSamples;
			int frameSize = 100;
			long long energyThreshold = 100000;

			int totalSamples = 16025 * 3;
			int numFrames = totalSamples / frameSize;

			for (int i = 0; i < numFrames; i++)
			{
				long long sum = 0;
				for (int j = 0; j < frameSize; j++)
				{
					int idx = i * frameSize + j;
					sum += (long long)waveIn[idx] * waveIn[idx];
				}
				long long avgEnergy = sum / frameSize;

				if (avgEnergy > energyThreshold)
				{
					for (int j = 0; j < frameSize; j++)
					{
						speechSamples.push_back((double)waveIn[i * frameSize + j]);
					}
				}
			}

			if (speechSamples.size() < 2500)
			{
				cout << "No speech detected or too short.\n";
				continue;
			}

			// Extract features from speechSamples
			vector<vector<double>> features;
			ComputeFeatures(speechSamples, features);
			vector<int> obsSeq;
			for (size_t i = 0; i < features.size(); i++)
			{
				obsSeq.push_back(cb.Quantize(features[i]));
			}

			double maxProb = -1e308;
			int recognizedDigit = -1;

			for (int m = 0; m <= 9; m++)
			{
				double prob = hmms[m].GetProbability(obsSeq);
				if (prob > maxProb)
				{
					maxProb = prob;
					recognizedDigit = m;
				}
			}

			cout << "Recognized Digit: " << recognizedDigit << "\n";
		}
		else if (option == 4)
		{
			break;
		}
	}
	return 0;
}

void StartRecord()
{
	const int NUMPTS = 16025 * 3; // 3 seconds
	int sampleRate = 16025;
	HWAVEIN hWaveIn;
	MMRESULT result;
	WAVEFORMATEX pFormat;
	pFormat.wFormatTag = WAVE_FORMAT_PCM;
	pFormat.nChannels = 1;
	pFormat.nSamplesPerSec = sampleRate;
	pFormat.nAvgBytesPerSec = sampleRate * 2;
	pFormat.nBlockAlign = 2;
	pFormat.wBitsPerSample = 16;
	pFormat.cbSize = 0;

	result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
	WAVEHDR WaveInHdr;
	WaveInHdr.lpData = (LPSTR)waveIn;
	WaveInHdr.dwBufferLength = NUMPTS * 2;
	WaveInHdr.dwBytesRecorded = 0;
	WaveInHdr.dwUser = 0L;
	WaveInHdr.dwFlags = 0L;
	WaveInHdr.dwLoops = 0L;
	waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	result = waveInStart(hWaveIn);
	printf("recording for 3 seconds...\n");
	Sleep(3 * 1000);
	waveInReset(hWaveIn);
	waveInUnprepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	waveInClose(hWaveIn);
}

string GenerateFilename(int digit, int utterance)
{
	stringstream ss;
	ss << "Digits/txt/254101020_E_" << digit << "_" << setfill('0') << setw(2) << utterance << ".txt";
	return ss.str();
}

vector<double> ReadSignal(string filename)
{
	vector<double> signal;
	ifstream f(filename.c_str());
	if (!f.is_open())
		return signal;
	double val;
	while (f >> val)
	{
		signal.push_back(val);
	}
	f.close();
	return signal;
}
