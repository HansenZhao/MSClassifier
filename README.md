### MS Classifier
### prepare your data
Please refer to example\PCA.csv
Note that the file format should be **.csv**
### begin with analysis
To begin with your classification analysis, your should create a MSPCA instance, which will read and standardize your original data. In MATLAB consol, type
```matlab
s = MSPCA(0.01,0.05);
```
to initialize your MSPCA instance. The first parameter stands for analysis accuracy, which is set to 0.01 m/z in this case; the second parameter represents minimum tolerance, which mean peaks intensity smaller than a specific fraction of the max intensity will be ignored.**s** is your instance name, you can assign any name you like such as **apple** or **banana**. In this case, we just use **s** for convenience.

You can add MS data directly from your CSV file by typing
```matlab
s.addMSFile();
```
then you can select your prepared CSV file in an UI window.

To plot PCA result, just type
```matlab
[coeff,score,lat] = s.plotPCA(0,0.95);
```
The first input parameter is a normalization tag. Setting the tag to 1 means perform the normalization before PCA analysis, which will get rid of the intensity difference between samples by normalize the maximum intensity in each sample to **1**. Setting the parameter to 0 will cancel the normalization. The second input parameter is minimum information lost, after the analysis, MSPCA will tell you how many PCs you need to reach the minimum information lost.

To visualize the contribute of each peak in the first two principle components, just type
```matlab
s.coefBar(coeff, 3)
```
coeff is the analysis result from *plotPCA* function in the last step. the second parameter stands for *mark tolerance*. Mark tolerance can be an integer larger than 1 to tell the program mark the m/z of the mass with top *n* contribution(when mark tolerance is set to *n*). Mark torlerance also can be a fraction ranging from 0 to 1 to tell the ratio to the maximum contribution need to be marked.

To check your original data, use
```matlab
s.plotMSByName('Hela-5')
```
the only parameter is the sample name, which is defined by your CSV file.

Have fun with it!
