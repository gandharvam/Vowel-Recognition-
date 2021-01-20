	204101034_vowelRecognition.cpp

This project is about Vowel Recognition.
 
I used several method to implement it.Those instruction are written in comment in the project code file.
For this , I recorded several utterances of each vowel.
Procedure :
Recorded each vowel at least 20 times. 10 will be used for training and 10 for testing.
Trimmed the vowel part manually through Cooledit, kept little bit of silence before and after the vowel.
Naming is like.
Eg: RollNo_vowel_utteranceNo.txt

....................................................................
Steps for vowel recognition:
Steps to generate reference file for 1 vowel
1) Took first vowel recording, do DC shift and normalization.
2) Selected 5 frames from the steady part and applied hamming window.
3) Computed Ri's, Ai's and Ci's of these frames and apply raised Sine window on Ci's
4) Repeated the above steps for 9 more recording of the same vowel (so now we have 50 rows of Ci values) and took the average of these recordings with respect to frames, i.e. frame 1 of all the 10 recordings should be considered for average and hence 10 rows will produce 1 row of Ci's. Similarly for all the frames, and so finally we will get 5 rows of Ci values (5 rows and 12 columns). Dumped these values in a text file.
5)Similarly we will get text file for remaining 9 vowels and these 10 text files will be used as reference file for vowels.   
6)Finally we have 5 row of ci for each vowel (ie.25 row and 12 column).

Steps for testing:
1)Took input files for testing (10 test files per vowel) and pass these test files in a loop so that we can check out of 10 files how many are recognized correctly
2) for each test file took 5 frames from stable part
3) Calculated the Tokhura's distance from each reference file.
4) the one with minimum distance will be recognized as a vowel.   

Please find the Tokhura weights:
1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0

I have implemented the live recornding to but it doesnot work properly.
So I didn't get the issue.