HybridNET
Provided files:
HybridNET.py:The similarity between the parent spectrum and the child spectra is calculated using a hybrid algorithm. After the first round of matching yields the highest similarity data for each pair of spectra, the matched child spectra are used as the parent spectra for the second round, and the hybrid algorithm is continuously applied to calculate their similarity with the unmatched child spectra. Through this iterative process, more compounds can be identified.

Requirements
The program requires the following Python libraries:
numpy
pandas
matchms
spec2vec
math
logging
pathlib
datetime

Installation command:
```
pip install numpy pandas matchms spec2vec math logging pathlib datetime
```

Important Notes
Input Data Requirements:
The script requires two input files in the working directory:
`parent.xlsx` - Parent spectra data
`product.xlsx` - Product spectra data
Parameter Configuration:
Similarity threshold can be adjusted by modifying the condition:
```
  THRESHOLD = 0.7  # Default threshold set to 0.7 for identifying similar spectra
```

Output Files
Primary Output:
`result` - Includes data processing result files and the calculation results of the similarity algorithms.
