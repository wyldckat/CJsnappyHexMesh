//reading the dictionary for the layer properties.
scalar nLayers=readScalar(meshDict.lookup("nLayers"));
scalar stretching=readScalar(meshDict.lookup("stretching"));
nLayers++;
// label patchName=readLabel(meshDict.lookup("patch"));
//Calculating the length of the cell that's going to be refined (in weighted terms).
scalar L[100];
scalar Li[100];
L[0]=0;
L[1]=1;
for(int j=1;j<nLayers-1;j++)
{
       L[j+1]=stretching*(L[j]-L[j-1])+L[j];
}
int k=0;
for(int i=nLayers-2;i>=1;i--)
{
       Li[k]=1-(L[i+1]-L[i])/L[i+1];
       Info << "weightedthings =" << Li[k] << endl;
       k=k+1;
}

