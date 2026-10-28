Code accompanying the paper "Bilayer Molecular Network Modeling Driven by Hybrid Spectral Similarity for Nontarget Annotation of Pharmaceuticals and Pesticides".

Files provided:

A. HybridNET.py: The main script, used to calculate the hybrid similarity between parent spectra and product spectra, and by using the product spectra identified in the first round as the parent spectra for the second round of hybrid calculation, to continue calculating similarity with the product spectra not identified in the first round.

B. test_parent.xlsx: A sample input file containing parent spectral data. The columns include:
    ID: A unique identifier for each parent spectrum.
    mz: Precursor mass-to-charge ratio (m/z).
    MS2: Fragmentation spectrum data, in the format peak intensity:mass-to-charge ratio.

C. test_product.xlsx: A sample input file containing product spectral data, in the same format as parent.xlsx.

1. Requirements:

This program requires the following Python libraries: numpy pandas matchms spec2vec math logging pathlib datetime.
You can install them using the following command:
pip install numpy pandas matchms spec2vec math logging pathlib datetime

2. Run the script:

To run the program, execute the following command in the terminal:
python HybridNET.py

3. Notes:

Input data: The script requires the files test_parent.xlsx and test_product.xlsx to exist in the working directory.
Parameters: You can adjust the similarity threshold by modifying THRESHOLD=0.7 in the code. The default threshold used to identify similar spectra is set to 0.7.

4. Output:

result_test.csv: This file contains the calculated similarity scores for each pair of parent and product spectra.

Citation: If you find the code useful, please cite the following paper:
"Bilayer Molecular Network Modeling Driven by Hybrid Spectral Similarity for Nontarget Annotation of Pharmaceuticals and Pesticides".
