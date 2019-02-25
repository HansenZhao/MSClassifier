## MS Classifier
### Fields
#### accuracy
The accuracy of peak alignment. Example: 0.01
#### minTor
The minimum relative peak intensity threshold. Example: 0.05
#### nMS
The number of sample.
#### nPeaks
The number of peaks.
#### data
The nMS x nPeaks matrix of raw data.
#### pks
The mz list of peaks.
### Methods
#### MSPCA
The construction method for MSPCA, two parameters to specify the alignment accuracy and relative intensity tolerence.
```matlab
instance = MSPCA(accuracy,minTor)
```
Example:
```matlab
s = MSPCA(0.01,0.05)
```
#### addMS
Add a new MS sample by **MS profile data**
```matlab
addMS(MSLoc,MSInts,name)
```
Example:
```matlab
s.addMS(mz,intensity,'HeLa-1');
```
#### subMSLoc
Delete a specific MS peak
```matlab
subMSLoc(MSLoc)
```
Example:
```matlab
s.subMSLoc(108.1);
```
#### addMSFile
Add MS samples from a csv file.
```matlab
addMSFile(isView=false)
```
Example:
```matlab
s.addMSFile(); %disable view by default
s.addMSFile(1); % view all load MS sample profile
```
#### getMSById
Get MS profile by ID.
+ raw: raw data input
+ pksLoc: peak location after alignment
+ pksInts: peak intensity
```matlab
[raw,pksLoc,pksInts] = getMSById(index)
```
Example:
```matlab
[a,b,c] = s.getMSById(1);
```
#### getMSByName
Get MS profile by sample name
+ raw: raw data input
+ pksLoc: peak location after alignment
+ pksInts: peak intensity
```matlab
[raw,pksLoc,pksInts] = getMSByName(name)
```
Example:
```matlab
[a,b,c] = s.getMSByName('Hela-1');
```
#### sortMS
Sort MS data by peak locations in ascend way.
```matlab
sortMS()
```
Example:
```matlab
s.sortMS()
```
#### getTag
Get sample group tag
+ name: sample names, string
+ tag: a nMS x 1 vector for sample group tag
+ dic: a dictionary map sample name to tag, Key: sample name identity Value: tag index
```matlab
[name,tag,dic] = getTag()
```
Example:
```matlab
[sampleNames,sampleTag,dictionary] = s.getTag();
```
#### pPCA
Process PCA algorithm.
+ isNormal: is normalize the raw data before processing.
+ res: a structure containing processing results
+ res.method: string, processing algorithm
+ res.score: nMS x p matrix for new p-dimension representation of the sample
+ res.coeff: loading for PCA
+ res.latent: p x 1 vector for weights of each Principle Component
```matlab
res = pPCA(isNormal)
```
Example:
```matlab
r = s.pPCA(1);
```
#### pTSNE
Process tSNE algorithm.
+ isNormal: is normalize the raw data before processing.
+ perp: Perplexity of tSNE (5~50)
+ dist: distance used in tSNE, refer to the MATLAB doc of tsne function.
+ res: a structure containing processing results
+ res.method: string, processing algorithm
+ res.score: nMS x p matrix for new p-dimension representation of the sample
```matlab
res = pTSNE(isNormal,perp=10,dist='euclidean')
```
Example:
```matlab
res = s.pTSNE(1); %Perplexity = 10 and dist = 'euclidean' by default
res = s.pTSNE(1,5); %Perplexity = 5 and dist = 'euclidean'
res = s.pTSNE(1,5,'cosine'); %Perplexity = 5 and dist = 'cosine'
```
#### coefBar
Plot loading bar
+ res: res structure output from pPCA method
+ markTor: if markTor > 1, mark the largest first markTor peak
+ if 0 < markTor < 1, mark the peaks whose loading weights > maxWeight x markTor
```matlab
coefBar(res,markTor)
```
Example:
```matlab
s.coefBar(res,5);
s.coefBar(res,0.5);
```
#### coefScatter
Plot loading scatter
+ coeff: coeff from pPCA method
+ markTor: if markTor > 1, mark the largest first markTor peak
+ if 0 < markTor < 1, mark the peaks whose loading weights > maxWeight x markTor
```matlab
coefScatter(coeff,markTor)
```
Example:
```matlab
s.coefScatter(res.coeff,5);
s.coefScatter(res.coeff,0.5);
```
#### plotMSByName
Plot MS profile by Name
+ name: string, sample name to be plotted
+ hAxes: the axes handle to plot
```matlab
plotMSByName(name,hAxes = gca)
```
Example:
```matlab
s.plotMSByName('HeLa-1') %create new figure
s.plotMSByName('HeLa-1',subplot(121)) %plot in specified axes
```
#### scatterRes
Plot MS profile by Name
+ res: res structure from pPCA or pTSNE
+ dim: plot dimension, default [1,2]
+ CRegOpt: Confidential Region option, default disable
+ CRegOpt.alpha: the confidential threshold
+ CRegOpt.norm: the assuming distribution of the data: defalut: 'norm' for normal distribution. Another option is 'exp'
```matlab
scatterRes(res,dim = [1,2],CRegOpt = [])
```
Example:
```matlab
s.scatterRes(res) %plot scatter figure with dim 1, dim 2, no confidential region
s.scatterRes(res,[1,3]) %plot scatter figure with dim 1, dim 3, no confidential region
s.scatterRes(res,[1,3],struct('alpha',0.05,'dist','norm')) %plot scatter figure with dim 1, dim 3, 95% confidence, assuming normal distribution
```
### Tutorial
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

process PCA with normalized data
```matlab
res = s.pPCA(1);
```

### visualize the result
loading plot
```matlab
s.coefScatter(res.coeff, 3)
```

check your original data, use
```matlab
s.plotMSByName('Hela-5')
```

scatter plot
```matlab
s.scatterRes(res);
```
or
```matlab
s.scatterRes(res,[1,2],struct('alpha',0.05,'dist','norm'));
```
Have fun with it!
