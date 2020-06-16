#!/opt/local/bin/python2.7
# Written by Gregg Rice 09/12/11 using functions from QuSHAPE
#



def scaleShapeData(data0,data1,rate=0.25):
    """   Scale Shape Data
    """
    # Data1 is scaled to Data1. 
    # Data0 is sorted and then the lower ones are used using the ratio. Data0[:N*rate] 
    N=len(data0) #100
    if rate>=1:
        A=data0.copy()
        B=data1.copy()
    else:
        A,B=selectDataForScale1(data0,data1,rate)
    
  #  A,B=removeOutlier(A,B)
    #newFactor= optimizeScaleFactor(A,B)
    
    aver=np.average(A)/np.average(B)
    k=0
    while k<3: 
        s=aver*0.8
        e=aver*1.2
        NScore=40
        testFactors=np.linspace(s,e,NScore)
        score=np.zeros(NScore)
        for i in np.arange(NScore):
            score[i]=scaleFactorFunc(testFactors[i],A,B)
        aver=testFactors[np.argmin(score)]
        k+=1
    
    newFactor=aver
      
    return newFactor

def selectDataForScale1(data0,data1,rate=0.25):
    """ Select the lowest RX  area with corresponding BG area
    """
    NData=len(data0)
    argSorted0=np.argsort(data0)
    NSelect=int(NData*rate)
    #s=int(NData*0.5)
    #e=int(NData*rate)
    selectedArgSortAreaRX=argSorted0[:NSelect]
    A=np.zeros(NSelect)
    B=np.zeros(NSelect)
    for i in range(len(selectedArgSortAreaRX)):
        ind=selectedArgSortAreaRX[i]
        A[i]=data0[ind]
        B[i]=data1[ind]
    return A,B

reportKeys=['seqNum','normDiff']
def DReport():
    dReport={}
    dReport['seqNum']=np.array([],dtype='i4')
    dReport['normDiff']=np.array([],dtype='i4')
    
    return dReport

def writeReportFile(dReport,fName):
    myfile=open(fName,'w')    
##    for key in reportKeys:
##        myfile.write(str(key)+'\t')
##    myfile.write('\n')
    for i in range(len(dReport['seqNum'])):
        myfile.write(str(dReport['seqNum'][i])+'\t'+str(round(dReport['normDiff'][i],4)))
        myfile.write('\n')
 
def getReportFromTxt(fName):
    fl=open(fName, "r")
    a,data=[],[]
    lines=fl.readlines()
    dReport=DReport()
    #for i in range(len(lines)-1,0,-1): #fixed error where the files were read in backwards.
		#                                   #I'll have to ask F why he did this...
		# first line doesn't have a header in these files
    for i in range(0,len(lines)):
        a= lines[i].split('\t')
        if minDataCut <= float(a[1]) < 0:
            a[1] = str(0.00) #NEGATIVE SHAPE NON-PHYSICAL
        if float(a[1]) < minDataCut:
            a[1] = str(0) # NO DATA!!! 
            noData.append(int(a[0]))
        dReport['seqNum']=np.append(dReport['seqNum'],int(a[0]))
        dReport['normDiff']=np.append(dReport['normDiff'],float(a[1].replace('\n','')))
        
    fl.close()
    return dReport
       
def normSimple(dataIn,POutlier=2.0, PAver=10.0):
    NData=len(dataIn)
    NOutlier=int(float(NData)*float(POutlier)/100.0)
    if NOutlier<1:
        NOutlier=1
    NAver=int(float(NData)*float(PAver)/100.0)
    dataSorted=np.sort(dataIn)
    aver=np.average(dataSorted[-NAver:-NOutlier])
    dataNormed=dataIn/aver
    return dataNormed

def findPOutlierStat(dataIn):
    # Methods : stats , box
    NData=len(dataIn)
    dataNormed=normStat(dataIn)
    outlierA=np.array([])
    averA=np.array([])
    for i in range(NData):
        if dataNormed[i]>3:
            outlierA=np.append(outlierA, dataNormed[i])
        elif dataNormed[i]>1:
            averA=np.append(averA, dataNormed[i])
        else:
            pass
    NOutlier=len(outlierA)
    NAver=len(averA)
    POutlier=float(NOutlier)/float(NData)*100.0
    NOutlier=float(NAver)/float(NData)*100.0+POutlier
    
    return POutlier, NOutlier

def normStat(data):
    normalized=np.zeros(len(data))
    mean=np.mean(data)
    std=np.std(data)
    normalized=(data-(mean))/std
    #normalized=normalized+1
    return normalized



def smoothRect(dataIn,degree=1):
    NData=len(dataIn)  
    dataOut=np.zeros(NData)
    window=degree*2+1
    for i in range(degree):
        dataOut[i]=np.average(dataIn[:i+degree+1])
    for i in range(1,degree+1):
        dataOut[-i]=np.average(dataIn[-(i+degree):])
    for i in range(degree,NData-degree):
        dataOut[i]=np.average(dataIn[i-degree:i+degree+1])
    return dataOut  

def fitLinear(x,y,NData):
    fittedData=np.zeros(NData)
    fittedData[0:x[0]]=y[0]
    fittedData[x[-1]:]=y[-1]
    NPoint=len(x)
    
    for i in range(NPoint-1):
        x1=np.array([x[i],x[i+1]])
        y1=np.array([y[i],y[i+1]])
        coeff=np.polyfit(x1,y1,1)
        poly=np.poly1d(coeff)
        xNew=np.arange(x[i],x[i+1])
        xNew=np.array(xNew,int)
        yNew=np.polyval(poly, xNew)
        fittedData[xNew]=yNew
    
    return fittedData

def scaleShapeDataWindow(data0,data1,deg=25,rate=1,step=10,fit=None,ref=None):
    win=2*deg+1
    N=len(data0)
    if N<win+step:
        scaleFactor=scaleShapeData(data0,data1,rate)
        #data11=data1*scaleFactor
        return scaleFactor
    
    aScaleFactor=np.array([])
    aX=np.array([])
    for i in range(0,N,step):
        if i<deg:
            s=0
            e=win
        elif i>N-deg:
            e=N
            s=N-win
        else:
            s=i-deg
            e=i+deg+1      
        partData0=data0[s:e]
        partData1=data1[s:e]
        scaleFactor=scaleShapeData(partData0,partData1,rate)
        aScaleFactor=np.append(aScaleFactor,scaleFactor)
        aX=np.append(aX,i)
  
    #aY=scipy.signal.medfilt(aScaleFactor,5) 
    aY=smoothRect(aScaleFactor,degree=2)
    aX=aX[1:-1]
    aY=aY[1:-1]
    fittedSig=fitLinear(aX,aY,len(data1))
    # data11=data1*fittedSig
    if fit=='linear':
        newX=np.arange(len(fittedSig))
        coeff=np.polyfit(newX,fittedSig,1)
        poly=np.poly1d(coeff)
        fittedSig=np.polyval(poly, newX)
    if fit=='exp':
        newX=np.arange(len(fittedSig))
    if ref==0:
        data11=data1*fittedSig
        return data11
    if ref==1:
        data00=data0/fittedSig
        return data00
    
    return fittedSig

 
def findRoiReports(seqNum0,seqNum1):
    N0=len(seqNum0)
    N1=len(seqNum1)
    s0,e0=0,N0  
    s1,e1=0,N1
    
    ok=True
    i=0
    while ok and i<N0:
        j=0
        while ok and j<N1:
            if seqNum0[i]==seqNum1[j]:
                s0,s1=i,j
                ok=False
            j+=1
        i+=1
    
    ok=True
    i=N0-1
    while ok  and i>=0:
        j=N1-1
        while ok  and j>=0:
            if seqNum0[i]==seqNum1[j]:
                e0,e1=i+1,j+1
                ok=False
            j-=1
        i-=1
    
    return s0,e0,s1,e1

def removeOutlier(A,B):
    A0=np.array([])
    B0=np.array([])
    fark=np.subtract(A,B)
    # fark=np.argsort(fark)
    fark=normStat(fark)
    for i in range(len(fark)):
        if fark[i]<2 and fark[i]>-2:
            A0=np.append(A0,A[i])
            B0=np.append(B0,B[i])
    return A0,B0

def addNoData(diffReport,noDataArray):
    # make a unique array for noData positions
    nd, diffReportOut = [], diffReport.copy()
    for i in noDataArray:
        ndata = set(nd)
        if not i in ndata:
            nd.append(i)

    # go through the diff report and set those positions
    # to -999 for the output file

    for i in nd:
        index = np.nonzero( diffReportOut['seqNum'] == i )[0][0]
        diffReportOut['normDiff'][index] = -999

    return diffReportOut
    

def optimizeScaleFactor(A,B):
    factor=1.0
    resultList= fmin(scaleFactorFunc, factor, args=(A,B),full_output=1,disp=0)
    if resultList[4]==0:
        scaleFactor=resultList[0]
    else:
        scaleFactor=1
    return float(scaleFactor)

def scaleFactorFunc(factor,A,B):
    err=np.sum(np.abs(A-factor*B))
    return err

def scaleSampleReactReport(dReport0,dReport1,isScale=True,window=25):
    dReport00=dReport0.copy()
    dReport11=dReport1.copy()
    
    s0,e0,s1,e1= findRoiReports(dReport0['seqNum'],dReport1['seqNum'])
    #print s0,e0,s1,e1  # above is region of intrest from both traces
    #   is this necessary if the data is avail and already curaited?
    for key in dReport00.keys():
        dReport00[key]=dReport00[key][s0:e0]
        dReport11[key]=dReport11[key][s1:e1] 
    if isScale:
        aScale=scaleShapeDataWindow(dReport0['normDiff'],dReport1['normDiff'],deg=window)

        dReport11['normDiff']=dReport11['normDiff']*aScale
        
        #aScale=scaleShapeDataWindow(dReport0['areaBG'],dReport1['areaBG'])
        #dReport11['areaBG']=dReport11['areaBG']*aScale
        
        #dReport11['areaDiff']=dReport11['areaRX']-dReport11['areaBG']
        POutlier,PAver=findPOutlierStat(dReport11['normDiff'])
        dReport11['normDiff']=normSimple(dReport11['normDiff'],POutlier,PAver)
       
        partData00,partData11=removeOutlier(dReport00['normDiff'],dReport11['normDiff'])
        scaleFactor=optimizeScaleFactor(partData00,partData11)
        dReport11['normDiff']=dReport11['normDiff']*scaleFactor
         
    return dReport00,dReport11

         
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 5:
        print 'Usage: <nmia.txt> <1m6.txt> <difference.dif> <i>'
        print 'window = 2*i+1 ... good place to start is 25'
        quit()
    import numpy as np 
    from pylab import figure,show,savefig,title
    from matplotlib.pyplot import setp
    from scipy.optimize import fmin
  #  from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    minDataCut = -0.4
    noData = []
    
    fig0 = plt.figure()
    ax00 = fig0.add_subplot(211)
    ax01 = fig0.add_subplot(212)
    
    fig1 = plt.figure()
    ax10 = fig1.add_subplot(111)
#   ax11 = fig1.add_subplot(212)

### SPECIFY THE FILE NAMES    
    fName0=sys.argv[1]
    fName1=sys.argv[2]
### GET THE DATA FROM THE REPORT FILES    
###    ##reads in the data from the output files from the qushape output files
    dReport0=getReportFromTxt(fName0)
    dReport1=getReportFromTxt(fName1)
    #print dReport0
### SCALE AND NORMALIZE SAMPLE DATA (dReport1) to REFERENCE DATA(dReport0)    
###     ## the 'meat and potatoes' I guess of this script
    dReport00,dReport11=scaleSampleReactReport(dReport0,dReport1,isScale=True,window=int(sys.argv[4]))

   
    
### SPECIFY THE FILE NAMES TO WRITE THE DATA TO TXT FILES
    diffReportfName=sys.argv[3]


    #make the diffReport file from subtraction
    diffReport = DReport()
    diffReport['seqNum'] = dReport00['seqNum']
    diffReport['normDiff'] = dReport00['normDiff']-dReport11['normDiff']

    #Go back and put in NODATA points into the output file and write it
    diffReportOut = addNoData(diffReport,noData)

    writeReportFile(diffReportOut,diffReportfName)
    
### PLOTTING FUNCTIONS   
    ax00.plot(dReport0['normDiff'],'r',linestyle='steps')
    ax00.plot(dReport1['normDiff'],'b',linestyle='steps')
    ax00.set_title('Reactivity before scaling')
    ax00.legend(['Reference', 'Sample'])
    
    ax01.plot(dReport00['normDiff'],'r',linestyle='steps')
    ax01.plot(dReport11['normDiff'],'b',linestyle='steps')
    ax01.legend(['Reference', 'Sample'])
    ax01.set_title('Reactivity after scaling')
    
    diff0=dReport0['normDiff']-dReport1['normDiff']
    diff1=dReport00['normDiff']-dReport11['normDiff']
    #print np.sum(np.abs(diff0)),np.sum(np.abs(diff1))
    zeros = np.zeros(120)   
    
    ax10.plot(diff0,'r',linestyle='steps')
    ax10.plot(diff1,'b',linestyle='steps')
    ax10.plot(zeros,'g')
  
    ax10.legend(['Before', 'After'])
    ax10.set_title('Difference')
    
    show()

    
    
