Steps for vowel recognition:
Generating the reference file for a vowel:
1. Take the first vowel recording, and perform a DC shift and normalization.
2. Select 5 frames from the steady part.
3. Compute the Ri's, Ai's, and Ci's for these frames, and apply a raised Sine window to the Ci's.
4. Repeat these steps for 19 more recordings of the same vowel.
5. This will give you 100 rows of Ci values for the vowel:
(20 recordings × 5 frames per recording × 12 Ci values per frame) = 100 rows.
6. Average the Ci values for corresponding frames across all 20 recordings (e.g., average all first frames, all second frames, etc.).
7. This results in 5 rows of Ci values (5 rows × 12 columns).
8. Save these values in a text file.
Repeat the process for the remaining vowels, resulting in 5 text files—one for each vowel.
Now for testing:
1. Take 10 test files per vowel, and process them in a loop to check how many are recognized correctly.
2. For each test file take 5 frames from the stable part.
3. Calculate Tokhura's distance for each test file using the reference files.
4. Since each test file has 5 frames and each reference file has 5 rows, calculate Tokhura’s distance for corresponding frames. Average the 5 distances to get the final distance.
5. Repeat this for each reference file.
6. The vowel with the smallest distance will be recognized as the correct vowel.
Tokhura weights:
1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0