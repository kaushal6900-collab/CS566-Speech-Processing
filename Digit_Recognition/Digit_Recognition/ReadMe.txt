HMM-Based Isolated Digit Recognition System

This project is a C++ implementation of an isolated digit recognition system (0-9). It utilizes Vector Quantization (Codebook) for feature mapping and Hidden Markov Models (HMM) for temporal pattern recognition. The system supports offline training/testing via file input and real-time testing via microphone recording.
📋 Features

    Feature Extraction: Converts raw audio signals into feature vectors (implied LPC/Cepstral analysis via FeatureExtraction.h).

    Vector Quantization (Codebook): clustering of feature vectors to generate observation sequences.

    HMM Training: Uses the Baum-Welch algorithm to train a specific model for each digit (0-9).

    Offline Testing: Evaluates accuracy using pre-recorded datasets.

    Live Testing: Records audio in real-time (3 seconds) using the Windows Multimedia API (winmm.lib), performs energy-based silence removal, and predicts the digit.

    Model Persistence: Saves and loads trained Codebooks and HMM parameters to text files.

🛠️ Prerequisites

    Operating System: Windows (Required for Windows.h and mmsystem.h).

    IDE/Compiler: Visual Studio (Recommended due to stdafx.h and _tmain usage).

    Libraries: winmm.lib (Linked automatically via #pragma comment).

📂 Project Structure & Data Organization

To run this system effectively, you must organize your source files and data directories as follows:
Plaintext

ProjectRoot/
├── main.cpp                // The core logic (provided code)
├── FeatureExtraction.h     // Header for feature computation
├── Codebook.h              // Header for VQ/Codebook logic
├── HMM.h                   // Header for HMM logic
├── stdafx.h                // Visual Studio precompiled header
├── codebook.txt            // Generated after training
├── model_digit_0.txt       // Generated after training
│   ...
└── Digits/                 // Data Directory
    └── txt/                // Input signal files
        ├── 254101020_E_0_01.txt
        ├── 254101020_E_0_02.txt
        └── ...

⚠️ Important: Data Naming Convention

The system is hardcoded to look for training/testing files with a specific naming pattern based on the Roll Number 254101020.

Format: Digits/txt/254101020_E_[Digit]_[Utterance].txt

    Digit: 0-9

    Utterance: 01-30 (Used for Training), 31-40 (Used for Offline Testing)

    File Content: The text files should contain raw signal amplitude values (floating point or integer), separated by whitespace or newlines.

🚀 How to Build and Run

    Open Visual Studio: Create a new C++ Console Application.

    Add Files: Add main.cpp, FeatureExtraction.h, Codebook.h, and HMM.h to your project.

    Disable Precompiled Headers (Optional): If you are not using stdafx.h in your other files, remove the #include "stdafx.h" line or configure your project settings to "Not Using Precompiled Headers".

    Build: Compile the solution (Ctrl+Shift+B).

    Run: Start the application (Ctrl+F5).

🎮 Usage Guide

Upon running the application, you will see the following menu:
Plaintext

--- Digit Recognition System ---
1. Train
2. Test (Offline)
3. Live Test
4. Exit

1. Training

Select Option 1.

    The system reads utterances 1 through 30 for all digits (0-9).

    It generates a global Codebook.

    It trains 10 separate HMMs (one for each digit).

    Output: codebook.txt and model_digit_X.txt files are saved to the disk.

2. Offline Testing

Select Option 2.

    Requirement: Models must be trained or loaded first.

    The system reads utterances 31 through 40.

    It compares the predicted digit against the actual filename digit.

    Output: Displays individual file recognition results and the overall accuracy percentage.

3. Live Testing

Select Option 3.

    Requirement: Microphone connected.

    The system records audio for 3 seconds.

    It applies Energy-Based Silence Removal to isolate the speech.

    Output: Prints the recognized digit to the console.

⚙️ Configuration Details

If you need to tweak the system performance, look for these variables in main.cpp:

    Sampling Rate: 16025 Hz

    Recording Duration: 3 seconds (16025 * 3 samples)

    Silence Removal Frame Size: 100 samples

    Energy Threshold: 100000 (Adjust this if live recording is not detecting your voice or picking up too much noise).

📝 Authors

    Student ID/Roll No: 254101020 (As inferred from the file naming convention).
