import numpy as np
import pandas as pd
from scipy import stats
from statistics import variance
# np.set_printoptions(threshold=np.inf)


data1 = pd.read_csv(r'E:\科研\Homo sapiens\BioGRID\分析\geneExpressionBladder.csv')
data2 = pd.read_csv(r'E:\科研\Homo sapiens\BioGRID\dynamic\HumoBioNet3.csv')

timeData = data1.iloc[:, 1:12].values
row = data2["row"].values
col = data2["col"].values

threshold=[]
for i in range(0, len(timeData)):
    nProtein=timeData[i]
    sum = 0
    avg = np.mean(nProtein)
    sigma = np.std(nProtein, ddof=1)
    var = variance(nProtein)
    # f = 1/(1+sigma*sigma)
    # threshold.append(avg+2*sigma*(1-f))
    f = var / (1 + var)
    threshold.append(avg+2*sigma*f)


threTimeDate = []
for i in range(0, len(timeData)):
    nProtein=timeData[i]
    thres = threshold[i]
    for j in range(0, len(timeData[0])):
        if(nProtein[j]<=thres):
            nProtein[j] = 0
        # else:
        #     nProtein[j] = 1
    threTimeDate.append(nProtein)
print("动态阈值计算完成")
thresholdM = np.array(threTimeDate)


def subNetBuild(subperiod, row, col):
    res=[]
    for i in range(0, len(row)):
        curRow = row[i]
        curCol = col[i]
        rowData = subperiod[curRow-1]
        colData = subperiod[curCol-1]
        flag = 0
        for j in range(0, len(subperiod[0])):
            if(flag==0):
                if(rowData[j]!=0 and colData[j]!=0):
                    res.append(1)
                    flag=1
            if(flag==1):
                break
        if(flag==0):
            res.append(0)

    row1=[]
    col1=[]
    for i in range(0, len(res)):
        if(res[i]==1):
            row1.append(row[i])
            col1.append(col[i])
    print("构建子网完成")
    print("对应关系：", len(row1))

    nodeNum=[]
    for i in range(1, len(threTimeDate)+1):
        flag=0
        for j in range(0, len(row1)):
            if(i == row1[j]):
                nodeNum.append(1)
                flag=1
                break
            if (i == col1[j]):
                nodeNum.append(1)
                flag=1
                break
        if(flag==0):
            nodeNum.append(0)

    return row1, col1, nodeNum

def ecc(row, col):
    eccScore = []
    for i in range(0, len(row)):
        targetRow = row[i]
        targetCol = col[i]
        connectRow = []
        connectCol = []
        commonNodeNum = 0
        # print(i)
        for j in range(0, len(row)):
            if (targetRow == row[j]):
                connectRow.append(col[j])
            elif (targetRow == col[j]):
                connectRow.append(row[j])
            if (targetCol == row[j]):
                connectCol.append(col[j])
            elif (targetCol == col[j]):
                connectCol.append(row[j])
        connectRowNum = len(connectRow)
        connectColNum = len(connectCol)
        minNum = connectRowNum if connectRowNum < connectColNum else connectColNum
        if (minNum == 1):
            eccScore.append(0)
            continue
        commonNode = [x for x in connectRow if x in connectCol]
        commonNodeNum = len(commonNode)
        ecc = commonNodeNum / (minNum - 1)
        eccScore.append(ecc)
    print("ecc")

    return eccScore


def n5(subperiod, eccScore, row, col):
    #person
    perScore = []
    for i in range(0,len(row)):
        curRow = row[i]
        curCol = col[i]
        c1 = subperiod[curRow-1]
        c2 = subperiod[curCol-1]
        per = stats.pearsonr(c1, c2)[0]
        perScore.append(per)
    print("person")

    multip=[]
    sum=[]
    for i in range(0, len(row)):
        m = float(eccScore[i])*float(perScore[i])
        multip.append(m)
        s = 0.5*(float(eccScore[i])+float(perScore[i]))
        sum.append(s)
    print("mutip,sum")

    pSocre = []
    wScore = []
    for i in range(1, len(subperiod)+1):
        ps=0
        ws=0
        # print(i)
        for j in range(0, len(multip)):
            if(i==row[j]):
                pNum = float(multip[j])
                ps += pNum
                wNum = float(sum[j])
                ws += wNum
            if(i==col[j]):
                pNum = float(multip[j])
                ps += pNum
                wNum = float(sum[j])
                ws += wNum
        pSocre.append(ps)
        wScore.append(ws)
    return pSocre, wScore


##四舍五入
period = len(timeData[0])
interval = int(np.rint(period/10))
print(interval)
##切分子周期
subPeriod1 = thresholdM[:, 0: period-4*interval]
subPeriod2 = thresholdM[:, interval: period-3*interval]
subPeriod3 = thresholdM[:, 2*interval: period-2*interval]
subPeriod4 = thresholdM[:, 3*interval: period-interval]
subPeriod5 = thresholdM[:, 4*interval: period]

#子网1
subRow1, subCol1, nodeNum1 = subNetBuild(subPeriod1, row, col)
eccScore1 = ecc(subRow1, subCol1)
p1, w1 = n5(subPeriod1, eccScore1, subRow1, subCol1)
print('子网1计算完成')
#子网2
subRow2, subCol2, nodeNum2 = subNetBuild(subPeriod2, row, col)
eccScore2 = ecc(subRow2, subCol2)
p2, w2 = n5(subPeriod2, eccScore2, subRow2, subCol2)
print('子网2计算完成')
#子网3
subRow3, subCol3, nodeNum3 = subNetBuild(subPeriod3, row, col)
eccScore3 = ecc(subRow3, subCol3)
p3, w3 = n5(subPeriod3, eccScore3, subRow3, subCol3)
print('子网3计算完成')
#子网4
subRow4, subCol4, nodeNum4 = subNetBuild(subPeriod4, row, col)
eccScore4 = ecc(subRow4, subCol4)
p4, w4 = n5(subPeriod4, eccScore4, subRow4, subCol4)
print('子网4计算完成')
subRow5, subCol5, nodeNum5 = subNetBuild(subPeriod5, row, col)
eccScore5 = ecc(subRow5, subCol5)
p5, w5 = n5(subPeriod5, eccScore5, subRow5, subCol5)
print('子网5计算完成')

pScore=[]
wScore=[]
nodeN=[]
for i in range(0, len(timeData)):
    p = p1[i]+p2[i]+p3[i]+p4[i]+p5[i]
    w = w1[i]+w2[i]+w3[i]+w4[i]+w5[i]
    nn = nodeNum1[i]+nodeNum2[i]+nodeNum3[i]+nodeNum4[i]+nodeNum5[i]
    pScore.append(p)
    wScore.append(w)
    nodeN.append(nn)

p5n=[]
w5n=[]
for i in range(0 , len(timeData)):
    if(nodeN[i] == 0):
        p5n.append(0)
        w5n.append(0)
        continue
    p = pScore[i]/nodeN[i]
    w = wScore[i]/nodeN[i]
    p5n.append(p)
    w5n.append(w)

data = {"p5n": p5n, "w5n": w5n}
data1 = pd.DataFrame(data)
data1.to_csv(r"E:\科研\Homo sapiens\BioGRID\dynamic\k=2\bladder\score.csv", index=False)
