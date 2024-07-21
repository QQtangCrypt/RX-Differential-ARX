"""
In this code we limit that the difference in key schdule is not over 2^(-64)
"""
import os
import time
import random
import itertools

FullRound = 4

BlockSize = 48
HalfBlockSize = 24

SearchRoundStart = 1
SearchRoundEnd = 11
InitialLowerBound = 0

RotateR = 8
RotateL = 3

# the round with the minuium
Fix = 0
Fixvalue = 2

GroupConstraintChoice = 1

# Parameters for choice 1
GroupNumForChoice1 = 1

DifferentialProbabilityBound = list([])
for i in range(20):
    DifferentialProbabilityBound += [0]


def GenSequentialEncoding(x, u, main_var_num, cardinalitycons, fout):
    n = main_var_num  # the number of wi
    k = cardinalitycons
    if (k > 0):
        clauseseq = "-" + str(x[0] + 1) + " " + str(u[0][0] + 1) + " 0" + "\n"
        fout.write(clauseseq)
        for j in range(1, k):
            clauseseq = "-" + str(u[0][j] + 1) + " 0" + "\n"
            fout.write(clauseseq)

        for i in range(1, n - 1):
            clauseseq = "-" + str(x[i] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"
            fout.write(clauseseq)
            clauseseq = "-" + str(u[i - 1][0] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"
            fout.write(clauseseq)
            clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][k - 1] + 1) + " 0" + "\n"
            fout.write(clauseseq)

        for j in range(1, k):
            for i in range(1, n - 1):
                clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][j - 1] + 1) + " " + str(
                    u[i][j] + 1) + " 0" + "\n"
                fout.write(clauseseq)
                clauseseq = "-" + str(u[i - 1][j] + 1) + " " + str(u[i][j] + 1) + " 0" + "\n"
                fout.write(clauseseq)

        clauseseq = "-" + str(x[n - 1] + 1) + " " + "-" + str(u[n - 2][k - 1] + 1) + " 0" + "\n"
        fout.write(clauseseq)

    if (k == 0):
        for i in range(n):
            clauseseq = "-" + str(x[i] + 1) + " 0" + "\n"
            fout.write(clauseseq)


def CountClausesInSequentialEncoding(main_var_num, cardinalitycons, clause_num):
    count = clause_num
    n = main_var_num
    k = cardinalitycons
    # print(count,n,k)
    if (k > 0):
        count += 1
        # print('1',count,k)
        for j in range(1, k):
            count += 1
        # print('2',count)
        for i in range(1, n - 1):
            count += 3
        # print('3',count)
        for j in range(1, k):
            for i in range(1, n - 1):
                count += 2
        # print('1',count)
        count += 1
    if (k == 0):
        for i in range(n):
            count += 1
    # print(count)
    return count


def GenMatsuiConstraint(x, u, n, k, left, right, m, g, fout):
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            if m < k:
                for i in range(right, right + 1):
                    clauseseq = "-" + str(u[i][m] + 1) + " 0" + "\n"
                    fout.write(clauseseq)
            if g > 0:
                for i in range(right, right + 1):
                    clauseseq = str(u[i][g - 1] + 1) + " 0" + "\n"
                    fout.write(clauseseq)

        if ((left > 0) and (right == n - 1)):
            for i in range(0, k - m):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(u[right - 1][i + m] + 1) + " 0" + "\n"
                fout.write(clauseseq)
            for i in range(0, k - m + 1):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(x[right] + 1) + " " + "-" + str(
                    u[right - 1][i + m - 1] + 1) + " 0" + "\n"
                fout.write(clauseseq)
        if ((left > 0) and (right < n - 1)):
            for i in range(0, k - m):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(u[right][i + m] + 1) + " 0" + "\n"
                fout.write(clauseseq)
    if (m == 0):
        for i in range(left, right + 1):
            clauseseq = "-" + str(x[i] + 1) + " 0" + "\n"
            fout.write(clauseseq)


def CountClausesForMatsuiStrategy(n, k, left, right, m, g, clausenum):
    count = clausenum
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            if m < k:
                for i in range(right, right + 1):
                    count += 1
            if g > 0:
                for i in range(right, right + 1):
                    count += 1

        if ((left > 0) and (right == n - 1)):
            for i in range(0, k - m):
                count += 1
            for i in range(0, k - m + 1):
                count += 1
        if ((left > 0) and (right < n - 1)):
            for i in range(0, k - m):
                count += 1
    if (m == 0):
        for i in range(left, right + 1):
            count += 1
    return count


def GenBitModularAddition(x, y, z, xx, yy, zz, fout):
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(xx + 1) + " " + "-" + str(
        yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(
        xx + 1) + " " + "-" + str(yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(
        xx + 1) + " " + "-" + str(yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + "-" + str(
        xx + 1) + " " + "-" + str(yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)


def GenBitXOR(x, y, z, fout):
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)

def GenEqual(x, y,fout):
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)

def GenEqualNot(x, y,fout):
    clauseseq = "-" +str(x + 1) + " " + "-" + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)

def GenBitwToAdditionPr(x,y,z,p0,p1,p2,e,fout):
    clauseseq = str(p0 + 1) + " " + "-" + str(p1 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(p2 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(p1 + 1) + " " + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)


def GenWToPr(x,y,z,p0,p1,p2,e,n,fout):
    for i in range(0, 1):
        GenBitwToAdditionPr(x[n-2], y[n-2], z[n-2], p0[i],p1[i],p2[i],e[i], fout)


def ToBinary(x,n):
    out=[0 for i in range(n)]
    t=x
    for i in range(n):
        out[n-1-i]=t%2
        t=int(t/2)
    return out

def GenXORC(x,y,vc,n,ro,fout):
    c=ToBinary(vc,n)
    rc=[0 for i in range(n)]
    ROL(c,rc,ro,n)#ro is the gama
    oc=[0 for i in range(n)]
    for i in range(n):
        oc[i]=c[i]^rc[i]
    for i in range(n):
        if oc[i]:
            GenEqualNot(x[i],y[i],fout)
        else:
            GenEqual(x[i],y[i],fout)



def GenAddition(x, y, z, n, fout):
    for i in range(0, n - 2):
        GenBitModularAddition(x[i], y[i], z[i], x[i + 1], y[i + 1], z[i + 1], fout)
    #GenBitXOR(x[n - 1], y[n - 1], z[n - 1], fout)


def GenBitAddPr(x, y, z, w, fout):
    clauseseq = "-" + str(x + 1) + " " + str(z + 1) + " " + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(y + 1) + " " + "-" + str(z + 1) + " " + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(
        w + 1) + " 0" + "\n"
    fout.write(clauseseq)


def GenAddPr(x, y, z, w, n, fout):
    for i in range(0, n - 2):
        GenBitAddPr(x[i + 1], y[i + 1], z[i + 1], w[i], fout)



def ROR(x, rx, t, n):
    for i in range(0, n):
        rx[i] = x[(i - t) % n]


def CountCluseRoundFunction(Round, clausesum):
    count = clausesum
    # nonzero input
    count = count + 1
    for r in range(0, Round):
        for i in range(0, HalfBlockSize - 2):
            count = count + 8
        #count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 5
        for i in range(0, 1):
            count = count + 11
        for i in range(0, HalfBlockSize):
            count = count + 4
        for i in range(0, HalfBlockSize):
            count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 8
        #count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 5
        for i in range(0, 1):
            count = count + 11
        for i in range(0, HalfBlockSize):
            count = count + 4
        for i in range(0, HalfBlockSize):
            count = count + 2
    return count


"""
Rotation to left Model
rx=x<<<t :rx and x are n bits(note the (n-1)-th bit is the least significant bit)
"""


def ROL(x, rx, t, n):
    for i in range(0, n):
        rx[i] = x[(i + t) % n]


def Per(x, px, P, n):
    for i in range(0, n):
        px[i] = x[P[i]]


def Decision(Round, n, ki,MatsuiRoundIndex, MatsuiCount, flag):
    num = 0
    count = 0
    NI = (n+1) * Round+Round
    MI = ki
    kNI=(n+1) * Round
    #ND= Round
    #MD=kd
    km=4
    # sl=8
    sr = 3
    ###########################################################################################################
    # the defination of variable in the round function
    x = []  # the left input
    y = []  # the right input
    w = []  # the auxiliary variable in computing the differential probability for modular addition
    p0=[]
    p1=[]
    p2=[]
    e=[]
    ui = []  # the auxiliary variable in integer
    #ud=[]#the decimal auxiliary variable
    kui = []  # the auxiliary variable in integer
    l=[]# the input in key schedule note the length is r+m-1
    k=[]#the subkey
    z=[]#the output of modular addition in round function
    kz=[]#the output of modualr addition in key schedule
    kw=[]#the variable of computing the probability in modular addition
    kp0=[]
    kp1=[]
    kp2=[]
    ke=[]
    #ki=11
    #kd=11
    # declar the variable in 1~Round round function
    for r in range(km-1):
        ltemp = []
        for i in range(0, n):
            ltemp.append(count)
            count = count + 1
        l.append(ltemp)
    for r in range(0, Round):
        xtemp = []
        for i in range(0, n):
            xtemp.append(count)
            count = count + 1
        x.append(xtemp)
        ytemp = []
        for i in range(0, n):
            ytemp.append(count)
            count = count + 1
        y.append(ytemp)
        wtemp = []
        for i in range(0, n - 1):
            wtemp.append(count)
            count = count + 1
        w.append(wtemp)
        p0temp = []
        for i in range(0, 1):
            p0temp.append(count)
            count = count + 1
        p0.append(p0temp)
        p1temp = []
        for i in range(0, 1):
            p1temp.append(count)
            count = count + 1
        p1.append(p1temp)
        p2temp = []
        for i in range(0, 1):
            p2temp.append(count)
            count = count + 1
        p2.append(p2temp)
        etemp = []
        for i in range(0, 1):
            etemp.append(count)
            count = count + 1
        e.append(etemp)
        ktemp = []
        for i in range(0, n):
            ktemp.append(count)
            count = count + 1
        k.append(ktemp)
        ltemp = []
        for i in range(0, n):
            ltemp.append(count)
            count = count + 1
        l.append(ltemp)
        kztemp = []
        for i in range(0, n):
            kztemp.append(count)
            count = count + 1
        kz.append(kztemp)
        ztemp = []
        for i in range(0, n):
            ztemp.append(count)
            count = count + 1
        z.append(ztemp)
        kwtemp = []
        for i in range(0, n - 1):
            kwtemp.append(count)
            count = count + 1
        kw.append(kwtemp)
        kp0temp = []
        for i in range(0, 1):
            kp0temp.append(count)
            count = count + 1
        kp0.append(kp0temp)
        kp1temp = []
        for i in range(0, 1):
            kp1temp.append(count)
            count = count + 1
        kp1.append(kp1temp)
        kp2temp = []
        for i in range(0, 1):
            kp2temp.append(count)
            count = count + 1
        kp2.append(kp2temp)
        ketemp = []
        for i in range(0, 1):
            ketemp.append(count)
            count = count + 1
        ke.append(ketemp)

    # declar the output variable in the Round-round function
    xtemp = []
    for i in range(0, n):
        xtemp.append(count)
        count = count + 1
    x.append(xtemp)
    ytemp = []
    for i in range(0, n):
        ytemp.append(count)
        count = count + 1
    y.append(ytemp)
    ltemp = []
    for i in range(0, n):
        ltemp.append(count)
        count = count + 1
    l.append(ltemp)
    ktemp = []
    for i in range(0, n):
        ktemp.append(count)
        count = count + 1
    k.append(ktemp)

    # declar the auxility variable in the boolean cardinality constraint
    for i in range(0, NI - 1):
        uitemp = []
        for j in range(0, MI):
            uitemp.append(count)
            count = count + 1
        ui.append(uitemp)
    """
    for i in range(0, ND - 1):
        udtemp = []
        for j in range(0, MD):
            udtemp.append(count)
            count = count + 1
        ud.append(udtemp)
    """
    # declar the auxility variable in the boolean cardinality constraint
    for i in range(0, NI - 1):
        kuitemp = []
        for j in range(0, 96):
            kuitemp.append(count)
            count = count + 1
        kui.append(kuitemp)

    countclause = 0
    # print(Round)
    # print(HalfBlockSize)
    # print(k)
    # print(N)
    countclause = CountCluseRoundFunction(Round, countclause)
    countclause = CountClausesInSequentialEncoding(NI, ki, countclause)
    #countclause = CountClausesInSequentialEncoding(ND, kd, countclause)

    countclause = CountClausesInSequentialEncoding(kNI, 96, countclause)
    # print(countclause)


    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = (HalfBlockSize + 2) * StartingRound
        RightNode = (HalfBlockSize + 2) * EndingRound - 1
        g = DifferentialProbabilityBound[EndingRound]
        PartialCardinalityCons = ki - DifferentialProbabilityBound[StartingRound] - DifferentialProbabilityBound[
            Round - EndingRound]
        countclause = CountClausesForMatsuiStrategy(NI, ki, LeftNode, RightNode, PartialCardinalityCons, g, countclause)
    # Open file


    # countclause = CountClausesForMatsuiStrategy(N, k, LeftNode, RightNode, fixvalue, countclause)

    kik=11
    kdk=11
    time_start = time.time()
    file = open("Problem-Round" + str(Round) + "-Probability" + str(ki)+ ".cnf", "w")
    file.write("p cnf " + str(count) + " " + str(countclause) + "\n")
    # print('===============================')

    # Add constraints for the round function
    for r in range(0, Round):
        rx = [0 for i in range(0, HalfBlockSize)]
        ry = [0 for i in range(0, HalfBlockSize)]
        rl = [0 for i in range(0, HalfBlockSize)]
        rk = [0 for i in range(0, HalfBlockSize)]
        ROR(x[r], rx, RotateR, HalfBlockSize)
        GenAddition(rx, y[r], z[r], HalfBlockSize, file)
        GenAddPr(rx, y[r], z[r], w[r], HalfBlockSize, file)
        GenWToPr(rx,y[r],z[r],p0[r],p1[r],p2[r],e[r],HalfBlockSize,file)
        for i in range(0,HalfBlockSize):
            GenBitXOR(k[r][i],z[r][i],x[r+1][i],file)
        ROL(y[r], ry, RotateL, HalfBlockSize)
        for i in range(0, HalfBlockSize):
            num = GenBitXOR(x[r + 1][i], ry[i], y[r + 1][i], file)

        #the propagation in keyschedule
        ROR(l[r], rl, RotateR, HalfBlockSize)
        GenAddition(rl, k[r], kz[r], HalfBlockSize, file)
        GenAddPr(rl, k[r], kz[r], kw[r], HalfBlockSize, file)
        GenWToPr(rl,k[r],kz[r],kp0[r],kp1[r],kp2[r],ke[r],HalfBlockSize,file)
        #there the l[i+m-1] m=4
        GenXORC(kz[r],l[r+3],r,HalfBlockSize,1,file)
        ROL(k[r], rk, RotateL, HalfBlockSize)
        for i in range(0, HalfBlockSize):
            num = GenBitXOR(l[r + 3][i], rk[i], k[r + 1][i], file)



    C = []
    for r in range(0, Round):
        for i in range(0, HalfBlockSize - 2):
            C.append(w[r][i])
        C.append(p0[r][0])
        C.append(p1[r][0])
        C.append(p2[r][0])
        C.append(e[r][0])
    GenSequentialEncoding(C, ui, NI, ki, file)
    kC = []
    for r in range(0, Round):
        for i in range(0, HalfBlockSize - 2):
            kC.append(kw[r][i])
        kC.append(kp0[r][0])
        kC.append(kp1[r][0])
        kC.append(kp2[r][0])
    GenSequentialEncoding(kC, kui, kNI, 96 , file)


    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = (HalfBlockSize +2) * StartingRound
        RightNode = (HalfBlockSize +2) * EndingRound - 1
        g = DifferentialProbabilityBound[EndingRound]
        PartialCardinalityCons = ki - DifferentialProbabilityBound[StartingRound] - DifferentialProbabilityBound[
            Round - EndingRound]
        GenMatsuiConstraint(C, u, NI, ki, LeftNode, RightNode, PartialCardinalityCons, g, file)


    # GenMatsuiConstraint(C, u, N, k, LeftNode, RightNode, fixvalue, file)

    # Add constraints to claim nonzero input difference
    clauseseq = ""
    for i in range(0, HalfBlockSize):
        clauseseq += str(x[0][i] + 1) + " " + str(y[0][i] + 1) + " "+ str(k[0][i] + 1) + " "+ str(l[0][i] + 1) + " "+ str(l[1][i] + 1) + " "+ str(l[2][i] + 1) + " "
    clauseseq += "0" + "\n"
    file.write(clauseseq)
    file.close()
    # Call solver cadical
    order = "/home/zhangyh/Software/cadical/build/cadical " + "Problem-Round" + str(Round) + "-Probability" + str(ki)+"_"+str(kd) + ".cnf > Round" + str(Round) + "-Probability" + str(ki)+"_"+str(kd)   + "-solution.out"
    os.system(order)
    # Extracting results
    order = "sed -n '/s SATISFIABLE/p' Round" + str(Round) + "-Probability" +  str(ki)  + "-solution.out > SatSolution.out"
    os.system(order)
    order = "sed -n '/s UNSATISFIABLE/p' Round" + str(Round) + "-Probability" +  str(ki)+ "-solution.out > UnsatSolution.out"
    os.system(order)
    satsol = open("SatSolution.out")
    unsatsol = open("UnsatSolution.out")
    satresult = satsol.readlines()
    unsatresult = unsatsol.readlines()
    satsol.close()
    unsatsol.close()
    if ((len(satresult) == 0) and (len(unsatresult) > 0)):
        flag = False
    if ((len(satresult) > 0) and (len(unsatresult) == 0)):
        flag = True
    order = "rm SatSolution.out"
    os.system(order)
    order = "rm UnsatSolution.out"
    os.system(order)
    # Removing cnf file
    order = "rm Problem-Round" + str(Round) + "-Probability" +  str(ki)+ ".cnf"
    os.system(order)
    time_end = time.time()
    # Printing solutions

    if (flag == True):
        print(
            "Round:" + str(Round) + "; Active: " +  str(ki)+ "; Sat; TotalCost: " + str(time_end - time_start))
    else:
        print("Round:" + str(Round) + "; Active: " +  str(ki)+ "; Unsat; TotalCost: " + str(
            time_end - time_start))

    return flag


# main function

TotalTimeStart = time.time()
#Bou=[[0,0],[1,1],[2,2],[3,0],[3,3],[4,1],[4,4],[5,2],[5,5],[6,3],[6,0]]
#Bou=[[0, 0], [1, 0], [1, 1], [2, 0], [2, 1], [2, 2], [3, 0], [3, 1], [3, 2], [4, 0], [3, 3], [4, 1], [4, 2], [5, 0], [4, 3], [5, 1], [4, 4], [5, 2], [6, 0], [5, 3], [6, 1], [5, 4], [6, 2], [7, 0], [5, 5], [6, 3], [7, 1], [6, 4], [7, 2], [8, 0], [6, 5], [7, 3], [8, 1], [6, 6], [7, 4], [8, 2], [9, 0], [7, 5], [8, 3], [9, 1], [7, 6], [8, 4], [9, 2], [7, 7], [10, 0], [8, 5], [9, 3], [10, 1], [8, 6], [9, 4], [10, 2], [8, 7], [11, 0], [9, 5], [10, 3], [8, 8], [11, 1], [9, 6], [10, 4], [11, 2], [9, 7], [12, 0], [10, 5], [11, 3], [9, 8], [12, 1], [10, 6], [11, 4], [9, 9], [12, 2], [10, 7], [13, 0], [11, 5], [12, 3], [10, 8], [13, 1], [11, 6], [12, 4], [10, 9], [13, 2], [11, 7], [14, 0], [12, 5], [10, 10], [13, 3], [11, 8], [14, 1], [12, 6], [13, 4], [11, 9], [14, 2], [12, 7], [15, 0], [13, 5], [11, 10], [14, 3], [12, 8], [15, 1], [13, 6], [11, 11], [14, 4], [12, 9], [15, 2], [13, 7], [16, 0], [14, 5], [12, 10], [15, 3], [13, 8], [16, 1], [14, 6], [12, 11], [15, 4], [13, 9], [16, 2], [14, 7], [12, 12], [17, 0], [15, 5], [13, 10], [16, 3], [14, 8], [17, 1], [15, 6], [13, 11], [16, 4], [14, 9], [17, 2], [15, 7], [13, 12], [18, 0], [16, 5], [14, 10], [17, 3], [15, 8], [13, 13], [18, 1], [16, 6], [14, 11], [17, 4], [15, 9], [18, 2], [16, 7], [14, 12], [19, 0], [17, 5], [15, 10], [18, 3], [16, 8], [14, 13], [19, 1], [17, 6], [15, 11], [18, 4], [16, 9], [14, 14], [19, 2], [17, 7], [15, 12], [20, 0], [18, 5], [16, 10], [19, 3], [17, 8], [15, 13], [20, 1], [18, 6], [16, 11], [19, 4], [17, 9], [15, 14], [20, 2], [18, 7], [16, 12], [21, 0], [19, 5], [17, 10], [15, 15], [20, 3], [18, 8], [16, 13], [21, 1], [19, 6], [17, 11], [20, 4], [18, 9], [16, 14], [21, 2], [19, 7], [17, 12], [22, 0], [20, 5], [18, 10], [16, 15], [21, 3], [19, 8], [17, 13], [22, 1], [20, 6], [18, 11], [16, 16], [21, 4], [19, 9], [17, 14], [22, 2], [20, 7], [18, 12], [23, 0], [21, 5], [19, 10], [17, 15], [22, 3], [20, 8], [18, 13], [23, 1], [21, 6], [19, 11], [17, 16], [22, 4], [20, 9], [18, 14], [23, 2], [21, 7], [19, 12], [24, 0], [17, 17], [22, 5], [20, 10], [18, 15], [23, 3], [21, 8], [19, 13], [24, 1], [22, 6], [20, 11], [18, 16], [23, 4], [21, 9], [19, 14], [24, 2], [22, 7], [20, 12], [25, 0], [18, 17], [23, 5], [21, 10], [19, 15], [24, 3], [22, 8], [20, 13], [25, 1], [18, 18], [23, 6], [21, 11], [19, 16], [24, 4], [22, 9], [20, 14], [25, 2], [23, 7], [21, 12], [26, 0], [19, 17], [24, 5], [22, 10], [20, 15], [25, 3], [23, 8], [21, 13], [26, 1], [19, 18], [24, 6], [22, 11], [20, 16], [25, 4], [23, 9], [21, 14], [26, 2], [19, 19], [24, 7], [22, 12], [27, 0], [20, 17], [25, 5], [23, 10], [21, 15], [26, 3], [24, 8], [22, 13], [27, 1], [20, 18], [25, 6], [23, 11], [21, 16], [26, 4], [24, 9], [22, 14], [27, 2], [20, 19], [25, 7], [23, 12], [28, 0], [21, 17], [26, 5], [24, 10], [22, 15], [27, 3], [20, 20], [25, 8], [23, 13], [28, 1], [21, 18], [26, 6], [24, 11], [22, 16], [27, 4], [25, 9], [23, 14], [28, 2], [21, 19], [26, 7], [24, 12], [29, 0], [22, 17], [27, 5], [25, 10], [23, 15], [28, 3], [21, 20], [26, 8], [24, 13], [29, 1], [22, 18], [27, 6], [25, 11], [23, 16], [28, 4], [21, 21], [26, 9], [24, 14], [29, 2], [22, 19], [27, 7], [25, 12], [23, 17], [28, 5], [26, 10], [24, 15], [29, 3], [22, 20], [27, 8], [25, 13], [23, 18], [28, 6], [26, 11], [24, 16], [29, 4], [22, 21], [27, 9], [25, 14], [23, 19], [28, 7], [26, 12], [24, 17], [29, 5], [22, 22], [27, 10], [25, 15], [23, 20], [28, 8], [26, 13], [24, 18], [29, 6], [27, 11], [25, 16], [23, 21], [28, 9], [26, 14], [24, 19], [29, 7], [27, 12], [25, 17], [23, 22], [28, 10], [26, 15], [24, 20], [29, 8], [27, 13], [25, 18], [23, 23], [28, 11], [26, 16], [24, 21], [29, 9], [27, 14], [25, 19], [28, 12], [26, 17], [24, 22], [29, 10], [27, 15], [25, 20], [28, 13], [26, 18], [24, 23], [29, 11], [27, 16], [25, 21], [28, 14], [26, 19], [24, 24], [29, 12], [27, 17], [25, 22], [28, 15], [26, 20], [29, 13], [27, 18], [25, 23], [28, 16], [26, 21], [29, 14], [27, 19], [25, 24], [28, 17], [26, 22], [29, 15], [27, 20], [25, 25], [28, 18], [26, 23], [29, 16], [27, 21], [28, 19], [26, 24], [29, 17], [27, 22], [28, 20], [26, 25], [29, 18], [27, 23], [28, 21], [26, 26], [29, 19], [27, 24], [28, 22], [29, 20], [27, 25], [28, 23], [29, 21], [27, 26], [28, 24], [29, 22], [27, 27], [28, 25], [29, 23], [28, 26], [29, 24], [28, 27], [29, 25], [28, 28], [29, 26], [29, 27], [29, 28], [29, 29]]
for u in range(1):
    CountProbability = InitialLowerBound
    #CountProbabilityd = InitialLowerBound
    Countindex=0
    for totalround in range(2,13):
        flag = False
        time_start = time.time()
        MatsuiRoundIndex = []
        MatsuiCount = 0

        if (GroupConstraintChoice == 1):
            for group in range(0, GroupNumForChoice1):
                for round in range(1, totalround - group + 1):
                    MatsuiRoundIndex.append([])
                    MatsuiRoundIndex[MatsuiCount].append(group)
                    MatsuiRoundIndex[MatsuiCount].append(group + round)
                    MatsuiCount += 1

            for group in range(totalround, totalround):
                for round in range(1, totalround):
                    MatsuiRoundIndex.append([])
                    MatsuiRoundIndex[MatsuiCount].append(group-round)
                    MatsuiRoundIndex[MatsuiCount].append(group)
                    MatsuiCount += 1

        # Printing Matsui conditions
        file = open("MatsuiCondition.out", "a")
        resultseq = "Round: " + str(totalround) + "; Partial Constraint Num: " + str(MatsuiCount) + "\n"
        file.write(resultseq)
        file.write(str(MatsuiRoundIndex) + "\n")
        file.close()
        while (flag == False):
            flag = Decision(totalround, HalfBlockSize, CountProbability, MatsuiRoundIndex, MatsuiCount, flag)
            CountProbability=CountProbability+1

        DifferentialProbabilityBound[totalround] = CountProbability - 1
        CountProbability = CountProbability-1
        time_end = time.time()
        file = open("RunTimeSummarise.out", "a")
        resultseq = "Round: " + str(totalround) + "; Differential Probability: " + str(
            DifferentialProbabilityBound[totalround]) + "; Runtime: " + str(time_end - time_start) + "\n"
        file.write(resultseq)
        file.close()

    # print(str(DifferentialProbabilityBound))
    # ooot.write(str(DifferentialProbabilityBound))
    # ooot.write("\n")
    TotalTimeEnd = time.time()
    # print("Total Runtime: " + str(TotalTimeEnd - TotalTimeStart))
    file = open("RunTimeSummarise.out", "a")
    resultseq = "Total Runtime: " + str(TotalTimeEnd - TotalTimeStart)
    file.write(resultseq)

